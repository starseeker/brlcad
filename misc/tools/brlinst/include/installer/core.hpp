// SPDX-License-Identifier: BSD-2-Clause
#pragma once
#include <json.hpp>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>
namespace installer {
namespace fs = std::filesystem;
using Json = nlohmann::json;
constexpr std::size_t buffer_size = 64 * 1024;
constexpr std::size_t max_manifest_size = 16 * 1024 * 1024;
constexpr std::uint64_t max_package_size = 0xff000000ULL;
constexpr const char *state_directory = ".brlinst";
constexpr const char *state_filename = "state.json";
constexpr const char *notices_filename = "NOTICE.txt";
constexpr const char *uninstall_filename = "uninstall.exe";
enum class ExitCode : int {
    success = 0,
    arguments = 2,
    existing = 3,
    elevation = 4,
    cancelled = 5,
    failure = 6
};
struct Error : std::runtime_error {
    ExitCode code;
    explicit Error(const std::string &message, ExitCode c = ExitCode::failure)
        : std::runtime_error(message), code(c) {}
};
struct Options {
    bool quiet = false, machine = false, start_menu = true, desktop = false, path = false;
    bool apply = false, help = false, recover = false;
    fs::path directory, log;
};
struct FileEntry {
    std::string path;
    std::uint64_t size = 0, compressed_size = 0, offset = 0;
    std::uint32_t crc = 0;
};
struct Package {
    fs::path executable;
    Json product;
    std::vector<FileEntry> files;
    std::uint64_t payload_offset = 0;
};
// Validated input shared by both output formats. Internal runtime files are
// supplied by each backend, never by the staged application.
struct PackageInput {
    Json product;
    fs::path icon;
    std::vector<std::pair<std::string, fs::path>> files;
};
PackageInput prepare_input(const fs::path &config, const fs::path &source, const fs::path &output);
std::string license_notices();
struct Fingerprint {
    std::uint64_t size;
    std::uint32_t crc;
};
using Progress = std::function<void(const std::string &, std::uint64_t, std::uint64_t)>;
using Cancel = std::function<bool()>;
std::string utf8(const fs::path &path);
fs::path path_from_utf8(const std::string &path);
std::string read_text(const fs::path &path, std::size_t limit = max_manifest_size);
Json parse_json(const std::string &text);
void flush_file(const fs::path &path);
void publish_file(const fs::path &source, const fs::path &destination);
void durable_json(const fs::path &path, const Json &data);
void validate_relative(const std::string &path, bool internal = false);
void validate_product(const Json &product, bool packaged = false);
void check_no_links(const fs::path &path);
std::string instance_id(const Json &product);
std::string display_name(const Json &product, bool machine);
Options parse_options(const std::vector<std::string> &args, bool uninstall = false);
std::string option_help(bool uninstall);
Fingerprint fingerprint(const fs::path &path);
void check_cancel(const Cancel &cancel);
Package open_package(const fs::path &executable);
void verify_package(const Package &package, const Progress &progress = {}, const Cancel &cancel = {});
void extract_file(const Package &package, const FileEntry &entry, const fs::path &output,
                  const Cancel &cancel = {});
void create_package(const fs::path &config, const fs::path &source, const fs::path &stub,
                    const fs::path &uninstaller, const fs::path &output);
struct Integration {
    virtual ~Integration() = default;
    virtual void install(Json &state, const fs::path &journal) = 0;
    virtual void remove(Json &state, const fs::path &journal) = 0;
};
Json load_state(const fs::path &root);
void install(const Package &package, const Options &options, Integration &integration,
             const Progress &progress = {}, const Cancel &cancel = {});
std::vector<std::string> uninstall(const fs::path &root, Integration &integration,
                                   const Progress &progress = {}, bool defer_uninstaller = false);
#ifdef _WIN32
void update_resources(const fs::path &executable, const Json &product, const fs::path &icon);
#endif
} // namespace installer
