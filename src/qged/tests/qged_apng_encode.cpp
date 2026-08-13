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
    bool removeInputs = false;
    bool padToMaximum = false;
    if (argc < 3 || argc > 5) {
	std::fprintf(stderr,
	    "Usage: %s frames.tsv output.apng [--pad-to-max] [--remove-inputs]\n",
	    argv[0]);
	return 2;
    }
    for (int i = 3; i < argc; ++i) {
	const std::string option(argv[i]);
	if (option == "--remove-inputs")
	    removeInputs = true;
	else if (option == "--pad-to-max")
	    padToMaximum = true;
	else {
	    std::fprintf(stderr, "Unknown option: %s\n", argv[i]);
	    return 2;
	}
    }

    std::vector<qged_raw_frame> frames;
    if (!qged_read_manifest(argv[1], frames)) {
	std::fprintf(stderr, "Unable to read frame manifest: %s\n", argv[1]);
	return 1;
    }
    size_t width = frames.front().width;
    size_t height = frames.front().height;
    if (padToMaximum) {
	for (const qged_raw_frame &frame : frames) {
	    width = std::max(width, frame.width);
	    height = std::max(height, frame.height);
	}
    }
    if (width > std::numeric_limits<uint32_t>::max() ||
	height > std::numeric_limits<uint32_t>::max() ||
	width > SIZE_MAX / height || width * height > SIZE_MAX / 4)
	return 1;
    for (const qged_raw_frame &frame : frames) {
	if (!padToMaximum &&
	    (frame.width != width || frame.height != height)) {
	    std::fprintf(stderr,
		"Frame dimensions changed within capture: %s\n",
		frame.path.c_str());
	    return 1;
	}
    }

    icv_anim_writer_t *animation = icv_anim_writer_create(argv[2],
	static_cast<uint32_t>(width), static_cast<uint32_t>(height),
	frames.size(), 60);
    if (!animation)
	return 1;
    std::vector<unsigned char> pixels;
    std::vector<unsigned char> canvas(width * height * 4, 0);
    for (size_t frameIndex = 0; frameIndex < frames.size(); ++frameIndex) {
	const qged_raw_frame &frame = frames[frameIndex];
	if (!qged_read_rgba(frame, pixels)) {
	    std::fprintf(stderr, "Unable to read captured frame: %s\n",
		frame.path.c_str());
	    (void)icv_anim_writer_close(animation);
	    return 1;
	}
	/* APNG requires a fixed canvas.  A resize capture deliberately changes
	 * frame dimensions, so center each actual frame without scaling on an
	 * opaque black maximum-size canvas.  This preserves pixel evidence and
	 * makes the widget extent change visible rather than distorting it. */
	std::fill(canvas.begin(), canvas.end(), 0);
	for (size_t pixel = 0; pixel < width * height; ++pixel)
	    canvas[pixel * 4 + 3] = 255;
	const size_t xOffset = (width - frame.width) / 2;
	const size_t yOffset = (height - frame.height) / 2;
	for (size_t y = 0; y < frame.height; ++y) {
	    const size_t targetY = yOffset + y;
	    for (size_t x = 0; x < frame.width; ++x) {
		const size_t source = (y * frame.width + x) * 4;
		const size_t target =
		    (targetY * width + xOffset + x) * 4;
		for (size_t channel = 0; channel < 4; ++channel)
		    canvas[target + channel] = pixels[source + channel];
	    }
	}
	uint64_t delay = 16667;
	if (frameIndex + 1 < frames.size() &&
	    frames[frameIndex + 1].elapsedUsec > frame.elapsedUsec)
	    delay = frames[frameIndex + 1].elapsedUsec - frame.elapsedUsec;
	delay = std::max<uint64_t>(1, std::min<uint64_t>(delay, 1000000));
	if (icv_anim_writer_add_rgba8(animation, canvas.data(),
		static_cast<uint32_t>(width), static_cast<uint32_t>(height),
		static_cast<uint32_t>(delay), 0, 0) != 0) {
	    (void)icv_anim_writer_close(animation);
	    return 1;
	}
    }

    const int status = icv_anim_writer_close(animation);
    if (status != 0) {
	std::fprintf(stderr, "Unable to encode APNG: %s\n", argv[2]);
	return 1;
    }
    if (removeInputs) {
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
