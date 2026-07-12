/*                  H O S T _ F A C T O R Y . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/host_factory.h"

#include <algorithm>
#include <cstring>
#include <mutex>
#include <new>
#include <string>
#include <vector>

struct brlobol_host_factory_token {
    brlobol_host_factory_token(void) : references(0), registered(false)
    {
	std::memset(&factory, 0, sizeof(factory));
    }

    struct brlobol_host_factory factory;
    std::string name;
    size_t references;
    bool registered;
};

static std::mutex registry_mutex;
static std::vector<brlobol_host_factory_token_t *> registry;

static size_t
factory_minimum_size(void)
{
    return offsetof(struct brlobol_host_factory, resize) +
	sizeof(((struct brlobol_host_factory *)0)->resize);
}

static struct brlobol_host_desc
sanitize_desc(const struct brlobol_host_desc *desc)
{
    struct brlobol_host_desc actual;
    std::memset(&actual, 0, sizeof(actual));
    actual.struct_size = sizeof(actual);
    actual.mode = BRLOBOL_HOST_MODE_HEADLESS;
    actual.width = 1;
    actual.height = 1;
    actual.device_pixel_ratio = 1.0;
    if (desc) {
	const size_t copy_size = std::min<size_t>(desc->struct_size,
	    sizeof(actual));
	if (copy_size)
	    std::memcpy(&actual, desc, copy_size);
	actual.struct_size = sizeof(actual);
    }
    if (!actual.width)
	actual.width = 1;
    if (!actual.height)
	actual.height = 1;
    if (actual.device_pixel_ratio <= 0.0)
	actual.device_pixel_ratio = 1.0;
    return actual;
}

static bool
factory_compatible(const brlobol_host_factory_token_t *token,
	const struct brlobol_host_desc &desc)
{
    return token &&
	(token->factory.capabilities & desc.required_capabilities) ==
	    desc.required_capabilities;
}

extern "C" brlobol_host_factory_token_t *
brlobol_host_factory_register(const struct brlobol_host_factory *factory)
{
    if (!factory || factory->abi_version != BRLOBOL_HOST_FACTORY_ABI_VERSION ||
	factory->struct_size < factory_minimum_size() || !factory->name ||
	!factory->name[0] || !factory->create || !factory->destroy ||
	!factory->bind_controller || !factory->request_frame || !factory->resize)
	return NULL;

    if ((factory->capabilities & BRLOBOL_HOST_CAP_READBACK) &&
	(!factory->capture || factory->struct_size < sizeof(*factory)))
	return NULL;

    std::lock_guard<std::mutex> lock(registry_mutex);
    for (const brlobol_host_factory_token_t *entry : registry) {
	if (entry->name == factory->name)
	    return NULL;
    }

    brlobol_host_factory_token_t *token =
	new (std::nothrow) brlobol_host_factory_token_t;
    if (!token)
	return NULL;

    const size_t copy_size = std::min<size_t>(factory->struct_size,
	sizeof(token->factory));
    std::memcpy(&token->factory, factory, copy_size);
    token->factory.struct_size = sizeof(token->factory);
    token->name = factory->name;
    token->factory.name = token->name.c_str();
    token->registered = true;
    registry.push_back(token);
    return token;
}

extern "C" int
brlobol_host_factory_unregister(brlobol_host_factory_token_t *token)
{
    if (!token)
	return 0;

    std::lock_guard<std::mutex> lock(registry_mutex);
    std::vector<brlobol_host_factory_token_t *>::iterator it =
	std::find(registry.begin(), registry.end(), token);
    if (it == registry.end() || !token->registered || token->references)
	return 0;

    registry.erase(it);
    token->registered = false;
    delete token;
    return 1;
}

extern "C" size_t
brlobol_host_factory_registry_count(void)
{
    std::lock_guard<std::mutex> lock(registry_mutex);
    return registry.size();
}

extern "C" size_t
brlobol_host_factory_registry_name(size_t index, char *name, size_t name_size)
{
    std::lock_guard<std::mutex> lock(registry_mutex);
    if (index >= registry.size())
	return 0;
    const std::string &factory_name = registry[index]->name;
    const size_t required = factory_name.size() + 1;
    if (name && name_size) {
	const size_t copied = std::min(factory_name.size(), name_size - 1);
	std::memcpy(name, factory_name.data(), copied);
	name[copied] = '\0';
    }
    return required;
}

extern "C" uint64_t
brlobol_host_factory_registry_capabilities(size_t index)
{
    std::lock_guard<std::mutex> lock(registry_mutex);
    return index < registry.size() ? registry[index]->factory.capabilities : 0;
}

extern "C" const char *
brlobol_host_factory_name(const brlobol_host_factory_token_t *token)
{
    return token ? token->name.c_str() : NULL;
}

extern "C" uint64_t
brlobol_host_factory_capabilities(const brlobol_host_factory_token_t *token)
{
    return token ? token->factory.capabilities : 0;
}

extern "C" brlobol_host_factory_token_t *
brlobol_host_factory_acquire(const char *name,
	const struct brlobol_host_desc *desc)
{
    const struct brlobol_host_desc actual = sanitize_desc(desc);
    std::vector<brlobol_host_factory_token_t *> candidates;
    {
	std::lock_guard<std::mutex> lock(registry_mutex);
	for (brlobol_host_factory_token_t *entry : registry) {
	    if (name && name[0] && entry->name != name)
		continue;
	    if (!factory_compatible(entry, actual))
		continue;
	    entry->references++;
	    candidates.push_back(entry);
	}
    }

    std::sort(candidates.begin(), candidates.end(),
	[](const brlobol_host_factory_token_t *a,
	   const brlobol_host_factory_token_t *b) {
	    if (a->factory.priority != b->factory.priority)
		return a->factory.priority > b->factory.priority;
	    return a->name < b->name;
	});

    brlobol_host_factory_token_t *selected = NULL;
    for (brlobol_host_factory_token_t *candidate : candidates) {
	if (!candidate->factory.probe ||
	    candidate->factory.probe(&actual,
		candidate->factory.user_data) > 0) {
	    selected = candidate;
	    break;
	}
    }

    for (brlobol_host_factory_token_t *candidate : candidates) {
	if (candidate != selected)
	    brlobol_host_factory_release(candidate);
    }
    return selected;
}

extern "C" void
brlobol_host_factory_release(brlobol_host_factory_token_t *token)
{
    if (!token)
	return;
    std::lock_guard<std::mutex> lock(registry_mutex);
    if (token->references)
	token->references--;
}

extern "C" int
brlobol_host_factory_instance_create(brlobol_host_factory_token_t *token,
	const struct brlobol_host_desc *desc, void *controller, void **instance)
{
    if (instance)
	*instance = NULL;
    if (!token || !controller || !instance)
	return 0;

    const struct brlobol_host_desc actual = sanitize_desc(desc);
    void *created = token->factory.create(&actual, token->factory.user_data);
    if (!created)
	return 0;
    if (token->factory.bind_controller(created, controller,
	    token->factory.user_data) <= 0) {
	token->factory.destroy(created, token->factory.user_data);
	return 0;
    }
    if (token->factory.open &&
	token->factory.open(created, &actual, token->factory.user_data) <= 0) {
	(void)token->factory.bind_controller(created, NULL,
	    token->factory.user_data);
	token->factory.destroy(created, token->factory.user_data);
	return 0;
    }

    *instance = created;
    return 1;
}

extern "C" void
brlobol_host_factory_instance_destroy(brlobol_host_factory_token_t *token,
	void *instance)
{
    if (!token || !instance)
	return;
    if (token->factory.close)
	token->factory.close(instance, token->factory.user_data);
    (void)token->factory.bind_controller(instance, NULL,
	token->factory.user_data);
    token->factory.destroy(instance, token->factory.user_data);
}

extern "C" int
brlobol_host_factory_instance_request_frame(
	brlobol_host_factory_token_t *token, void *instance, const char *reason)
{
    if (!token || !instance || !token->factory.request_frame)
	return 0;
    return token->factory.request_frame(instance, reason,
	token->factory.user_data);
}

extern "C" int
brlobol_host_factory_instance_resize(brlobol_host_factory_token_t *token,
	void *instance, unsigned int width, unsigned int height,
	double device_pixel_ratio)
{
    if (!token || !instance || !token->factory.resize || !width || !height ||
	device_pixel_ratio <= 0.0)
	return 0;
    return token->factory.resize(instance, width, height,
	device_pixel_ratio, token->factory.user_data);
}

extern "C" int
brlobol_host_factory_instance_capture(brlobol_host_factory_token_t *token,
	void *instance, unsigned char **pixels, size_t *size,
	unsigned int *width, unsigned int *height, unsigned int *components)
{
    if (!token || !instance || !token->factory.capture || !pixels || !size ||
	!width || !height || !components)
	return 0;
    return token->factory.capture(instance, pixels, size, width, height,
	components, token->factory.user_data);
}

extern "C" int
brlobol_host_factory_instance_dimensions(brlobol_host_factory_token_t *token,
	void *instance, unsigned int *width, unsigned int *height,
	double *device_pixel_ratio)
{
    if (!token || !instance || !token->factory.dimensions || !width ||
	!height || !device_pixel_ratio)
	return 0;
    return token->factory.dimensions(instance, width, height,
	device_pixel_ratio, token->factory.user_data);
}
