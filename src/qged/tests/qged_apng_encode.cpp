/*               Q G E D _ A P N G _ E N C O D E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "bu/app.h"
#include "icv.h"

#include <algorithm>
#include <charconv>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

struct qged_raw_frame {
    std::string path;
    size_t width = 0;
    size_t height = 0;
    uint64_t elapsedUsec = 0;
};

static bool
qged_parse_size(const std::string &value, size_t *result)
{
    if (!result || value.empty())
	return false;
    uint64_t parsed = 0;
    const std::from_chars_result converted = std::from_chars(
	value.data(), value.data() + value.size(), parsed);
    if (converted.ec != std::errc() ||
	converted.ptr != value.data() + value.size() ||
	parsed > static_cast<uint64_t>(SIZE_MAX))
	return false;
    *result = static_cast<size_t>(parsed);
    return true;
}

static bool
qged_parse_usec(const std::string &value, uint64_t *result)
{
    if (!result || value.empty())
	return false;
    const std::from_chars_result converted = std::from_chars(
	value.data(), value.data() + value.size(), *result);
    return converted.ec == std::errc() &&
	converted.ptr == value.data() + value.size();
}

static bool
qged_read_manifest(const char *path, std::vector<qged_raw_frame> &frames)
{
    std::ifstream input(path);
    if (!input.is_open())
	return false;
    std::string line;
    while (std::getline(input, line)) {
	const size_t first = line.find('\t');
	const size_t second = first == std::string::npos ?
	    first : line.find('\t', first + 1);
	const size_t third = second == std::string::npos ?
	    second : line.find('\t', second + 1);
	if (first == std::string::npos || second == std::string::npos ||
	    third == std::string::npos)
	    return false;
	qged_raw_frame frame;
	frame.path = line.substr(0, first);
	if (!qged_parse_size(line.substr(first + 1, second - first - 1),
		&frame.width) ||
	    !qged_parse_size(line.substr(second + 1, third - second - 1),
		&frame.height) ||
	    !qged_parse_usec(line.substr(third + 1), &frame.elapsedUsec) ||
	    !frame.width || !frame.height)
	    return false;
	frames.push_back(std::move(frame));
    }
    return !frames.empty();
}

static bool
qged_read_rgba(const qged_raw_frame &frame,
    std::vector<unsigned char> &pixels)
{
    if (frame.width > SIZE_MAX / frame.height ||
	frame.width * frame.height > SIZE_MAX / 4)
	return false;
    const size_t bytes = frame.width * frame.height * 4;
    std::ifstream input(frame.path, std::ios::binary | std::ios::ate);
    if (!input.is_open() || input.tellg() < 0 ||
	static_cast<uint64_t>(input.tellg()) !=
	    static_cast<uint64_t>(bytes))
	return false;
    input.seekg(0, std::ios::beg);
    pixels.resize(bytes);
    return static_cast<bool>(input.read(
	reinterpret_cast<char *>(pixels.data()), bytes));
}

int
main(int argc, char **argv)
{
    bu_setprogname(argv[0]);
    if (argc < 3 || argc > 4 ||
	(argc == 4 && std::string(argv[3]) != "--remove-inputs")) {
	std::fprintf(stderr,
	    "Usage: %s frames.tsv output.apng [--remove-inputs]\n",
	    argv[0]);
	return 2;
    }

    std::vector<qged_raw_frame> frames;
    if (!qged_read_manifest(argv[1], frames)) {
	std::fprintf(stderr, "Unable to read frame manifest: %s\n", argv[1]);
	return 1;
    }
    const size_t width = frames.front().width;
    const size_t height = frames.front().height;
    if (width > std::numeric_limits<uint32_t>::max() ||
	height > std::numeric_limits<uint32_t>::max())
	return 1;
    for (const qged_raw_frame &frame : frames) {
	if (frame.width != width || frame.height != height) {
	    std::fprintf(stderr,
		"Frame dimensions changed within capture: %s\n",
		frame.path.c_str());
	    return 1;
	}
    }

    icv_anim_t *animation = icv_anim_create(ICV_ANIM_APNG,
	static_cast<uint32_t>(width), static_cast<uint32_t>(height), 60);
    if (!animation)
	return 1;
    std::vector<unsigned char> pixels;
    for (size_t frameIndex = 0; frameIndex < frames.size(); ++frameIndex) {
	const qged_raw_frame &frame = frames[frameIndex];
	if (!qged_read_rgba(frame, pixels)) {
	    std::fprintf(stderr, "Unable to read captured frame: %s\n",
		frame.path.c_str());
	    icv_anim_destroy(animation);
	    return 1;
	}
	icv_image_t *image = icv_create_with_channels(
	    width, height, ICV_COLOR_SPACE_RGB, 4);
	if (!image) {
	    icv_anim_destroy(animation);
	    return 1;
	}
	for (size_t y = 0; y < height; ++y) {
	    const size_t targetY = height - y - 1;
	    for (size_t x = 0; x < width; ++x) {
		const size_t source = (y * width + x) * 4;
		const size_t target = (targetY * width + x) * 4;
		for (size_t channel = 0; channel < 4; ++channel)
		    image->data[target + channel] =
			ICV_CONV_8BIT(pixels[source + channel]);
	    }
	}
	if (icv_anim_add_frame(animation, image) != 0) {
	    icv_destroy(image);
	    icv_anim_destroy(animation);
	    return 1;
	}
	icv_destroy(image);
	uint64_t delay = 16667;
	if (frameIndex + 1 < frames.size() &&
	    frames[frameIndex + 1].elapsedUsec > frame.elapsedUsec)
	    delay = frames[frameIndex + 1].elapsedUsec - frame.elapsedUsec;
	delay = std::max<uint64_t>(1, std::min<uint64_t>(delay, 1000000));
	if (icv_anim_set_frame_delay(animation, frameIndex,
		static_cast<uint32_t>(delay)) != 0) {
	    icv_anim_destroy(animation);
	    return 1;
	}
    }

    const int status = icv_anim_write(animation, argv[2]);
    icv_anim_destroy(animation);
    if (status != 0) {
	std::fprintf(stderr, "Unable to encode APNG: %s\n", argv[2]);
	return 1;
    }
    if (argc == 4) {
	for (const qged_raw_frame &frame : frames)
	    (void)std::remove(frame.path.c_str());
    }
    return 0;
}

/*
 * Local Variables:
 * tab-width: 8
 * mode: C++
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
