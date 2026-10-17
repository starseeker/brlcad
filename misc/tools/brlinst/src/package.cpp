// SPDX-License-Identifier: BSD-2-Clause
#include "installer/core.hpp"
#include <zlib.h>
#include <algorithm>
#include <array>
#include <fstream>
#include <limits>
#include <map>
#include <set>
#ifdef _WIN32
#include <windows.h>
#endif

namespace installer {
namespace {
constexpr std::array<char, 16> magic{
    {'C', 'A', 'D', 'I', 'N', 'S', 'T', 'A', 'L', 'L', 'E', 'R', '0', '0', '1', '\0'}};
constexpr std::uint64_t footer_size = 48;
constexpr std::uint64_t max_file_size = 64ULL * 1024 * 1024 * 1024;
constexpr std::uint64_t max_total_size = 1024ULL * 1024 * 1024 * 1024;
std::uint64_t number(std::istream &in, unsigned bytes) {
    std::uint64_t n = 0;
    for (unsigned i = 0; i < bytes; ++i) {
        int c = in.get();
        if (c == EOF)
            throw Error("Truncated package header");
        n |= static_cast<std::uint64_t>(c) << (i * 8);
    }
    return n;
}
void number(std::ostream &out, std::uint64_t value, unsigned bytes) {
    for (unsigned i = 0; i < bytes; ++i)
        out.put(static_cast<char>(value >> (i * 8)));
}
std::uint64_t pe_end(std::istream &in, std::uint64_t length) {
    if (length < 256)
        throw Error("Runtime is not a Windows x64 executable");
    in.seekg(0);
    if (number(in, 2) != 0x5a4d)
        throw Error("Runtime has no MZ header");
    in.seekg(0x3c);
    auto pe = number(in, 4);
    if (pe > length - 200)
        throw Error("Invalid PE header offset");
    in.seekg(static_cast<std::streamoff>(pe));
    if (number(in, 4) != 0x4550 || number(in, 2) != 0x8664)
        throw Error("Expected Windows x64 PE runtime");
    in.seekg(static_cast<std::streamoff>(pe + 24));
    if (number(in, 2) != 0x20b)
        throw Error("Expected PE32+ runtime");
    in.seekg(static_cast<std::streamoff>(pe + 24 + 108));
    if (number(in, 4) < 5)
        return length;
    in.seekg(static_cast<std::streamoff>(pe + 24 + 112 + 4 * 8));
    auto certificates = number(in, 4);
    auto size = number(in, 4);
    if (certificates == 0 && size == 0)
        return length;
    if (certificates < pe + 200 || certificates > length || size != length - certificates || certificates % 8)
        throw Error("Invalid PE certificate table");
    return certificates;
}
std::string folded(const std::string &path) {
#ifdef _WIN32
    auto w = path_from_utf8(path).wstring();
    CharUpperBuffW(w.data(), static_cast<DWORD>(w.size()));
    return utf8(fs::path(w));
#else
    std::string result = path;
    std::transform(result.begin(), result.end(), result.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    return result;
#endif
}
void check_entries(const std::vector<FileEntry> &entries, std::uint64_t available) {
    std::set<std::string> files, dirs;
    std::uint64_t offset = 0, total = 0;
    bool has_uninstaller = false;
    if (entries.empty() || entries.size() > 200000)
        throw Error("Invalid package file count");
    for (const auto &entry : entries) {
        validate_relative(entry.path, true);
        const bool internal = entry.path.rfind(std::string(state_directory) + "/", 0) == 0;
        if (internal && entry.path != std::string(state_directory) + "/" + uninstall_filename &&
            entry.path != std::string(state_directory) + "/" + notices_filename)
            throw Error("Unexpected internal payload file");
        if (!internal)
            validate_relative(entry.path);
        if (entry.path == std::string(state_directory) + "/" + uninstall_filename)
            has_uninstaller = true;
        auto key = folded(entry.path);
        if (!files.insert(key).second)
            throw Error("Duplicate Windows filename: " + entry.path);
        auto parent = path_from_utf8(key).parent_path();
        while (!parent.empty()) {
            dirs.insert(parent.generic_u8string());
            parent = parent.parent_path();
        }
        if (entry.offset != offset || offset > available || entry.compressed_size > available - offset ||
            entry.size > max_file_size)
            throw Error("Invalid payload file bounds");
        offset += entry.compressed_size;
        if (entry.size > max_total_size - total)
            throw Error("Expanded payload exceeds limit");
        total += entry.size;
    }
    for (const auto &dir : dirs)
        if (files.count(dir))
            throw Error("File/directory collision: " + dir);
    if (!has_uninstaller)
        throw Error("Package has no uninstaller");
    if (available < offset || available - offset > 7)
        throw Error("Unexpected bytes between payload and manifest");
}
void inflate_entry(const Package &package, const FileEntry &entry, std::ostream *out, const Cancel &cancel) {
    std::ifstream in(package.executable, std::ios::binary);
    if (!in)
        throw Error("Cannot reopen package");
    in.seekg(static_cast<std::streamoff>(package.payload_offset + entry.offset));
    z_stream stream{};
    if (inflateInit(&stream) != Z_OK)
        throw Error("Cannot initialize decompressor");
    struct End {
        z_stream *s;
        ~End() {
            inflateEnd(s);
        }
    } end{&stream};
    std::array<unsigned char, buffer_size> input{}, output{};
    std::uint64_t left = entry.compressed_size, produced = 0;
    uLong checksum = 0;
    int result = Z_OK;
    while (result != Z_STREAM_END) {
        check_cancel(cancel);
        if (stream.avail_in == 0 && left) {
            auto n = static_cast<std::size_t>(std::min<std::uint64_t>(left, input.size()));
            if (!in.read(reinterpret_cast<char *>(input.data()), static_cast<std::streamsize>(n)))
                throw Error("Truncated compressed file");
            left -= n;
            stream.next_in = input.data();
            stream.avail_in = static_cast<uInt>(n);
        }
        stream.next_out = output.data();
        stream.avail_out = static_cast<uInt>(output.size());
        result = inflate(&stream, Z_NO_FLUSH);
        auto n = output.size() - stream.avail_out;
        if (result != Z_OK && result != Z_STREAM_END)
            throw Error("Corrupt compressed file: " + entry.path);
        if (n > entry.size - produced)
            throw Error("Expanded file exceeds declared size");
        produced += n;
        checksum = crc32(checksum, output.data(), static_cast<uInt>(n));
        if (out && !out->write(reinterpret_cast<char *>(output.data()), static_cast<std::streamsize>(n)))
            throw Error("Cannot write " + entry.path);
        if (result != Z_STREAM_END && !left && !stream.avail_in && !n)
            throw Error("Truncated deflate stream");
    }
    if (left || stream.avail_in || produced != entry.size || checksum != entry.crc)
        throw Error("Payload integrity check failed: " + entry.path);
}
FileEntry deflate_file(const fs::path &source, const std::string &name, std::ostream &out,
                       std::uint64_t origin) {
    check_no_links(source);
    if (!fs::is_regular_file(source))
        throw Error("Not a regular file: " + utf8(source));
    FileEntry entry;
    entry.path = name;
    entry.offset = static_cast<std::uint64_t>(out.tellp()) - origin;
    std::ifstream in(source, std::ios::binary);
    if (!in)
        throw Error("Cannot read " + utf8(source));
    z_stream stream{};
    if (deflateInit(&stream, Z_BEST_COMPRESSION) != Z_OK)
        throw Error("Cannot initialize compressor");
    struct End {
        z_stream *s;
        ~End() {
            deflateEnd(s);
        }
    } end{&stream};
    std::array<unsigned char, buffer_size> input{}, output{};
    int result = Z_OK;
    do {
        in.read(reinterpret_cast<char *>(input.data()), input.size());
        auto n = static_cast<std::size_t>(in.gcount());
        if (!in && !in.eof())
            throw Error("Payload read error");
        entry.size += n;
        if (entry.size > max_file_size)
            throw Error("File exceeds package size limit");
        entry.crc = static_cast<std::uint32_t>(crc32(entry.crc, input.data(), static_cast<uInt>(n)));
        stream.next_in = input.data();
        stream.avail_in = static_cast<uInt>(n);
        const int flush = in.eof() ? Z_FINISH : Z_NO_FLUSH;
        do {
            stream.next_out = output.data();
            stream.avail_out = static_cast<uInt>(output.size());
            result = deflate(&stream, flush);
            if (result != Z_OK && result != Z_STREAM_END)
                throw Error("Compression failed");
            auto used = output.size() - stream.avail_out;
            if (!out.write(reinterpret_cast<char *>(output.data()), static_cast<std::streamsize>(used)))
                throw Error("Package write failed");
        } while (stream.avail_out == 0);
    } while (result != Z_STREAM_END);
    entry.compressed_size = static_cast<std::uint64_t>(out.tellp()) - origin - entry.offset;
    return entry;
}
} // namespace

Package open_package(const fs::path &executable) {
    Package package;
    package.executable = executable;
    std::ifstream in(executable, std::ios::binary);
    if (!in)
        throw Error("Cannot read package");
    auto length = fs::file_size(executable);
    if (length > max_package_size)
        throw Error("Package exceeds supported executable size");
    auto end = pe_end(in, length);
    if (end < footer_size)
        throw Error("Missing package footer");
    in.seekg(static_cast<std::streamoff>(end - footer_size));
    std::array<char, 16> actual{};
    in.read(actual.data(), actual.size());
    if (actual != magic)
        throw Error("No installer payload found");
    auto payload = number(in, 8), manifest = number(in, 8), size = number(in, 8);
    auto crc = number(in, 4), version = number(in, 4);
    if (version != 1 || payload < 256 || manifest < payload || manifest > end - footer_size ||
        size != end - footer_size - manifest || size > max_manifest_size)
        throw Error("Invalid package manifest bounds");
    in.seekg(static_cast<std::streamoff>(manifest));
    std::string text(static_cast<std::size_t>(size), '\0');
    if (!in.read(text.data(), static_cast<std::streamsize>(size)))
        throw Error("Cannot read package manifest");
    if (crc32(0, reinterpret_cast<const Bytef *>(text.data()), static_cast<uInt>(text.size())) != crc)
        throw Error("Corrupt package manifest");
    auto data = parse_json(text);
    if (data.at("formatVersion") != 1)
        throw Error("Unsupported package version");
    package.product = data.at("product");
    validate_product(package.product, true);
    package.payload_offset = payload;
    if (!data.at("files").is_array() || data.at("files").size() > 200000)
        throw Error("Invalid package file list");
    for (const auto &item : data.at("files")) {
        for (auto key : {"size", "compressedSize", "offset", "crc"})
            if (!item.at(key).is_number_unsigned())
                throw Error("Invalid file metadata");
        if (item.at("crc").get<std::uint64_t>() > UINT32_MAX)
            throw Error("Invalid file checksum");
        package.files.push_back({item.at("path").get<std::string>(), item.at("size"),
                                 item.at("compressedSize"), item.at("offset"), item.at("crc")});
    }
    check_entries(package.files, manifest - payload);
    return package;
}
void verify_package(const Package &p, const Progress &progress, const Cancel &cancel) {
    std::uint64_t total = 0, done = 0;
    for (const auto &e : p.files)
        total += e.size;
    for (const auto &e : p.files) {
        if (progress)
            progress("Verifying " + e.path, done, total);
        inflate_entry(p, e, nullptr, cancel);
        done += e.size;
    }
}
void extract_file(const Package &p, const FileEntry &e, const fs::path &output, const Cancel &cancel) {
    check_no_links(output);
    if (fs::exists(output))
        throw Error("Refusing to replace existing file: " + utf8(output));
    std::ofstream out(output, std::ios::binary);
    if (!out)
        throw Error("Cannot create " + utf8(output));
    inflate_entry(p, e, &out, cancel);
    out.flush();
    if (!out)
        throw Error("Cannot flush " + utf8(output));
    out.close();
    flush_file(output);
}
void create_package(const fs::path &config, const fs::path &source, const fs::path &stub,
                    const fs::path &uninstaller, const fs::path &output) {
    auto input = prepare_input(config, source, output);
    auto &product = input.product;
    auto &icon = input.icon;
    auto inputs = std::move(input.files);
    inputs.insert(inputs.begin(), {std::string(state_directory) + "/" + uninstall_filename, uninstaller});
    inputs.emplace_back(std::string(state_directory) + "/" + notices_filename, fs::path{});
    std::sort(inputs.begin() + 1, inputs.end());
    // Reject collisions before creating output (including Windows case aliases).
    std::vector<FileEntry> preflight;
    for (const auto &i : inputs)
        preflight.push_back({i.first, 0, 0, 0, 0});
    check_entries(preflight, 0);
    {
        std::ifstream in(stub, std::ios::binary);
        if (pe_end(in, fs::file_size(stub)) != fs::file_size(stub))
            throw Error("Use an unsigned setup stub; sign the assembled installer");
    }
    {
        std::ifstream in(uninstaller, std::ios::binary);
        pe_end(in, fs::file_size(uninstaller));
    }
    auto temporary = output;
    temporary += ".building";
    check_no_links(temporary);
    if (fs::exists(temporary))
        throw Error("Remove stale package output: " + utf8(temporary));
    try {
        fs::copy_file(stub, temporary);
#ifdef _WIN32
        update_resources(temporary, product, icon);
#else
        if (!icon.empty())
            throw Error("Branding executable resources requires Windows");
#endif
        std::ofstream out(temporary, std::ios::binary | std::ios::app);
        const auto payload = fs::file_size(temporary);
        Json files = Json::array();
        std::vector<FileEntry> entries;
        for (const auto &input : inputs) {
            FileEntry entry;
            if (input.second.empty()) {
                const auto text = license_notices();
                uLongf length = compressBound(static_cast<uLong>(text.size()));
                std::vector<Bytef> compressed(length);
                if (compress2(compressed.data(), &length, reinterpret_cast<const Bytef *>(text.data()),
                              static_cast<uLong>(text.size()), Z_BEST_COMPRESSION) != Z_OK)
                    throw Error("Cannot compress dependency notices");
                entry = {input.first, text.size(), length, static_cast<std::uint64_t>(out.tellp()) - payload,
                         static_cast<std::uint32_t>(crc32(0, reinterpret_cast<const Bytef *>(text.data()),
                                                          static_cast<uInt>(text.size())))};
                if (!out.write(reinterpret_cast<const char *>(compressed.data()),
                               static_cast<std::streamsize>(length)))
                    throw Error("Cannot write dependency notices");
            } else
                entry = deflate_file(input.second, input.first, out, payload);
            files.push_back({{"path", entry.path},
                             {"size", entry.size},
                             {"compressedSize", entry.compressed_size},
                             {"offset", entry.offset},
                             {"crc", entry.crc}});
            entries.push_back(entry);
        }
        std::string manifest = Json{{"formatVersion", 1}, {"product", product}, {"files", files}}.dump();
        if (manifest.size() > max_manifest_size)
            throw Error("Manifest exceeds size limit");
        // SignTool aligns certificate data to eight bytes. Align our footer now
        // so its position is identical before and after Authenticode signing.
        while ((static_cast<std::uint64_t>(out.tellp()) + manifest.size() + footer_size) % 8)
            out.put('\0');
        const auto manifest_offset = static_cast<std::uint64_t>(out.tellp());
        check_entries(entries, manifest_offset - payload);
        out.write(manifest.data(), static_cast<std::streamsize>(manifest.size()));
        out.write(magic.data(), magic.size());
        number(out, payload, 8);
        number(out, manifest_offset, 8);
        number(out, manifest.size(), 8);
        number(out,
               crc32(0, reinterpret_cast<const Bytef *>(manifest.data()), static_cast<uInt>(manifest.size())),
               4);
        number(out, 1, 4);
        out.flush();
        if (!out || static_cast<std::uint64_t>(out.tellp()) > max_package_size)
            throw Error("Package too large or write failed");
        out.close();
        verify_package(open_package(temporary));
        flush_file(temporary);
        publish_file(temporary, output);
    } catch (...) {
        std::error_code ec;
        fs::remove(temporary, ec);
        throw;
    }
}
} // namespace installer
