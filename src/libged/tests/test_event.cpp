/*                    T E S T _ E V E N T . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 */

#include "common.h"

#include <atomic>
#include <string>
#include <thread>
#include <vector>

#include "bu/app.h"
#include "ged.h"

static int failures = 0;

#define CHECK(expression, message) do { \
    if (!(expression)) { \
	bu_log("FAIL: %s\n", (message)); \
	failures++; \
    } \
} while (0)

struct observation {
    std::vector<std::string> calls;
    std::vector<std::string> names;
    std::vector<enum ged_event_kind> kinds;
    ged_event_observer_token self = 0;
    ged_event_observer_token other = 0;
    bool publish_followup = false;
};

static void
record_events(struct ged *UNUSED(gedp), const struct ged_event *events,
	      size_t event_count, const struct ged_event_result *result,
	      void *client_data)
{
    observation *state = static_cast<observation *>(client_data);
    state->calls.push_back("record");
    CHECK(ged_event_result_coalesced_count(result) == event_count,
	"observer result describes delivered events");
    for (size_t i = 0; i < event_count; i++) {
	state->names.push_back(events[i].name ? events[i].name : "");
	state->kinds.push_back(events[i].kind);
    }
}

static void
record_internal(struct ged *UNUSED(gedp),
		const struct ged_event *UNUSED(events),
		size_t UNUSED(event_count),
		const struct ged_event_result *UNUSED(result),
		void *client_data)
{
    observation *state = static_cast<observation *>(client_data);
    state->calls.push_back("internal");
}

static void
record_post(struct ged *UNUSED(gedp),
	    const struct ged_event *UNUSED(events),
	    size_t UNUSED(event_count),
	    const struct ged_event_result *UNUSED(result),
	    void *client_data)
{
    observation *state = static_cast<observation *>(client_data);
    state->calls.push_back("post");
}

static void
remove_observers(struct ged *gedp,
		 const struct ged_event *UNUSED(events),
		 size_t UNUSED(event_count),
		 const struct ged_event_result *UNUSED(result),
		 void *client_data)
{
    observation *state = static_cast<observation *>(client_data);
    state->calls.push_back("remove");
    CHECK(ged_event_observer_remove(gedp, state->self) == 1,
	"observer removes itself during dispatch");
    CHECK(ged_event_observer_remove(gedp, state->other) == 1,
	"observer removes a later observer during dispatch");
}

static void
publish_followup(struct ged *gedp, const struct ged_event *events,
		 size_t event_count,
		 const struct ged_event_result *UNUSED(result),
		 void *client_data)
{
    observation *state = static_cast<observation *>(client_data);
    state->calls.push_back(events[0].name ? events[0].name : "");
    if (!state->publish_followup || event_count != 1 ||
	!events[0].name || !BU_STR_EQUAL(events[0].name, "first"))
	return;
    state->publish_followup = false;
    struct ged_event event;
    ged_event_init(&event);
    event.kind = GED_EVENT_OBJECT_MODIFIED;
    event.name = "second";
    CHECK(ged_event_publish(gedp, &event, 1, nullptr) == GED_EVENT_OK,
	"reentrant publication is accepted");
}

static void
test_batching_and_payloads(struct ged *gedp)
{
    observation state;
    ged_event_observer_token observer = ged_event_observer_add(gedp,
	GED_EVENT_OBSERVER_POST_RECONCILE, record_events, &state);
    CHECK(observer != 0, "post-reconcile observer is registered");

    struct ged_event_result *result = ged_event_result_create();
    char copied_name[] = "alpha";
    struct ged_event events[3];
    for (size_t i = 0; i < 3; i++)
	ged_event_init(&events[i]);
    events[0].kind = GED_EVENT_OBJECT_MODIFIED;
    events[0].name = copied_name;
    events[1] = events[0];
    events[2].kind = GED_EVENT_COMB_INSTANCE_REMOVED;
    events[2].parent_name = "comb";
    events[2].child_name = "alpha";
    events[2].path = "comb/alpha";

    CHECK(ged_event_batch_begin(gedp) == GED_EVENT_OK,
	"outer batch begins");
    CHECK(ged_event_batch_begin(gedp) == GED_EVENT_OK,
	"nested batch begins");
    CHECK(ged_event_publish(gedp, events, 3, result) == GED_EVENT_OK,
	"batched events are accepted");
    bu_strlcpy(copied_name, "other", sizeof(copied_name));
    CHECK(state.calls.empty(), "nested batch does not dispatch early");
    CHECK(ged_event_batch_end(gedp, result) == GED_EVENT_OK,
	"nested batch ends");
    CHECK(state.calls.empty(), "inner batch end does not dispatch");
    CHECK(ged_event_batch_end(gedp, result) == GED_EVENT_OK,
	"outer batch ends");

    CHECK(state.calls.size() == 1, "batch dispatches once");
    CHECK(state.names.size() == 2 && state.names[0] == "alpha",
	"publisher payloads are copied before dispatch");
    CHECK(ged_event_result_event_count(result) == 3,
	"result records valid input count");
    CHECK(ged_event_result_coalesced_count(result) == 2,
	"exact duplicate events are coalesced");
    bool found_path = false;
    for (size_t i = 0; i < ged_event_result_path_count(result); i++) {
	if (BU_STR_EQUAL(ged_event_result_path_at(result, i), "comb/alpha")) {
	    found_path = true;
	    break;
	}
    }
    CHECK(found_path, "result reports affected paths");

    CHECK(ged_event_observer_remove(gedp, observer) == 1,
	"observer is removed");
    ged_event_result_destroy(result);
}

static void
test_phase_order_and_removal(struct ged *gedp)
{
    observation phases;
    ged_event_observer_token internal = ged_event_observer_add(gedp,
	GED_EVENT_OBSERVER_INTERNAL, record_internal, &phases);
    ged_event_observer_token post = ged_event_observer_add(gedp,
	GED_EVENT_OBSERVER_POST_RECONCILE, record_post, &phases);

    struct ged_event event;
    ged_event_init(&event);
    event.kind = GED_EVENT_OBJECT_ADDED;
    event.name = "ordered";
    CHECK(ged_event_publish(gedp, &event, 1, nullptr) == GED_EVENT_OK,
	"ordered event publishes");
    CHECK(phases.calls.size() == 2 && phases.calls[0] == "internal" &&
	phases.calls[1] == "post", "internal observers precede post observers");
    ged_event_observer_remove(gedp, internal);
    ged_event_observer_remove(gedp, post);

    observation removal;
    removal.self = ged_event_observer_add(gedp,
	GED_EVENT_OBSERVER_INTERNAL, remove_observers, &removal);
    removal.other = ged_event_observer_add(gedp,
	GED_EVENT_OBSERVER_INTERNAL, record_internal, &removal);
    CHECK(ged_event_publish(gedp, &event, 1, nullptr) == GED_EVENT_OK,
	"removal event publishes");
    CHECK(removal.calls.size() == 1 && removal.calls[0] == "remove",
	"removed later observer is skipped in current dispatch");
    CHECK(ged_event_publish(gedp, &event, 1, nullptr) == GED_EVENT_OK,
	"event publishes after observer removal");
    CHECK(removal.calls.size() == 1,
	"self-removed observer does not run again");
}

static void
test_reentrancy(struct ged *gedp)
{
    observation reentrant;
    reentrant.publish_followup = true;
    ged_event_observer_token observer = ged_event_observer_add(gedp,
	GED_EVENT_OBSERVER_POST_RECONCILE, publish_followup, &reentrant);
    struct ged_event event;
    ged_event_init(&event);
    event.kind = GED_EVENT_OBJECT_MODIFIED;
    event.name = "first";
    CHECK(ged_event_publish(gedp, &event, 1, nullptr) == GED_EVENT_OK,
	"initial reentrant event publishes");
    CHECK(reentrant.calls.size() == 2 && reentrant.calls[0] == "first" &&
	reentrant.calls[1] == "second",
	"reentrant event is deferred until current dispatch completes");
    ged_event_observer_remove(gedp, observer);
}

static void
test_librt_bridge(void)
{
    struct db_i *dbip = db_open_inmem();
    CHECK(dbip != DBI_NULL, "in-memory database opens");
    if (dbip == DBI_NULL)
	return;

    struct bu_vls address = BU_VLS_INIT_ZERO;
    bu_vls_sprintf(&address, "%p", (void *)dbip);
    struct ged *gedp = ged_open("inmem", bu_vls_cstr(&address), 0);
    bu_vls_free(&address);
    CHECK(gedp != GED_NULL, "GED context adopts in-memory database");
    if (gedp == GED_NULL) {
	db_close(dbip);
	return;
    }

    observation state;
    ged_event_observer_token observer = ged_event_observer_add(gedp,
	GED_EVENT_OBSERVER_POST_RECONCILE, record_events, &state);
    int type = ID_SPH;
    struct directory *dp = db_diradd(gedp->dbip, "direct.s",
	RT_DIR_PHONY_ADDR, 0, RT_DIR_SOLID, &type);
    CHECK(dp != RT_DIR_NULL, "direct librt directory addition succeeds");
    CHECK(state.kinds.size() == 1 &&
	state.kinds[0] == GED_EVENT_OBJECT_ADDED &&
	state.names[0] == "direct.s",
	"direct librt addition publishes a GED event");

    if (dp != RT_DIR_NULL)
	CHECK(db_dirdelete(gedp->dbip, dp) == 0,
	    "direct librt directory removal succeeds");
    CHECK(state.kinds.size() == 2 &&
	state.kinds[1] == GED_EVENT_OBJECT_REMOVED &&
	state.names[1] == "direct.s",
	"direct librt removal publishes a GED event");

    ged_event_observer_remove(gedp, observer);
    ged_close(gedp);
}

static void
test_concurrent_publication(struct ged *gedp)
{
    constexpr size_t thread_count = 4;
    constexpr size_t events_per_thread = 25;
    std::atomic<size_t> publish_failures(0);
    observation state;
    ged_event_observer_token observer = ged_event_observer_add(gedp,
	GED_EVENT_OBSERVER_POST_RECONCILE, record_events, &state);

    std::vector<std::thread> publishers;
    for (size_t thread_index = 0; thread_index < thread_count;
	 thread_index++) {
	publishers.emplace_back([gedp, thread_index, events_per_thread,
		&publish_failures]() {
	    for (size_t event_index = 0; event_index < events_per_thread;
		 event_index++) {
		std::string name = "parallel_" +
		    std::to_string(thread_index) + "_" +
		    std::to_string(event_index);
		struct ged_event event;
		ged_event_init(&event);
		event.kind = GED_EVENT_OBJECT_MODIFIED;
		event.name = name.c_str();
		if (ged_event_publish(gedp, &event, 1, nullptr) != GED_EVENT_OK)
		    publish_failures++;
	    }
	});
    }
    for (std::thread &publisher : publishers)
	publisher.join();

    CHECK(publish_failures == 0, "concurrent events publish successfully");
    CHECK(state.kinds.size() == thread_count * events_per_thread,
	"concurrent publications are serialized without loss");
    ged_event_observer_remove(gedp, observer);
}

int
main(int argc, char **argv)
{
    bu_setprogname(argv[0]);
    CHECK(argc == 1, "test takes no arguments");

    struct ged context;
    ged_init(&context);
    test_batching_and_payloads(&context);
    test_phase_order_and_removal(&context);
    test_reentrancy(&context);
    test_concurrent_publication(&context);
    test_librt_bridge();

    CHECK(ged_event_batch_end(&context, nullptr) == GED_EVENT_INVALID,
	"ending a nonexistent batch is rejected");
    CHECK(ged_event_publish(&context, nullptr, 0, nullptr) ==
	GED_EVENT_INVALID, "empty publication is rejected");
    ged_free(&context);

    if (failures)
	bu_log("%d event test failure(s)\n", failures);
    return failures ? 1 : 0;
}

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
