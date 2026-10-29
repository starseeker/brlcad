// SPDX-License-Identifier: BSD-2-Clause
#include "installer/core.hpp"
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
using namespace installer;
namespace {
constexpr std::uint16_t pe_machine_x64 = 0x8664;
constexpr std::uint16_t pe_machine_arm64 = 0xaa64;
void require(bool ok, const char *message) {
    if (!ok)
        throw std::runtime_error(message);
}
template <class F> void fails(F action, const char *message) {
    bool failed = false;
    try {
        action();
    } catch (const std::exception &) {
        failed = true;
    }
    require(failed, message);
}
void write(const fs::path &p, const std::string &data) {
    std::ofstream f(p, std::ios::binary);
    f << data;
    if (!f)
        throw std::runtime_error("test write failed");
}
void put(std::string &s, std::size_t at, std::uint64_t value, unsigned bytes) {
    for (unsigned i = 0; i < bytes; ++i)
        s.at(at + i) = static_cast<char>(value >> (i * 8));
}
std::string stub(std::uint16_t machine = pe_machine_x64) {
    std::string s(512, '\0');
    put(s, 0, 0x5a4d, 2);
    put(s, 0x3c, 64, 4);
    put(s, 64, 0x4550, 4);
    put(s, 68, machine, 2);
    put(s, 88, 0x20b, 2);
    put(s, 88 + 108, 16, 4);
    return s;
}
struct Backend : Integration {
    bool fail = false, fail_remove = false;
    void install(Json &s, const fs::path &journal) override {
        s["integration"]["test"] = true;
        durable_json(journal, s);
        if (fail)
            throw Error("simulated registration failure");
    }
    void remove(Json &s, const fs::path &journal) override {
        if (fail_remove)
            throw Error("simulated removal failure");
        s["integration"] = Json::object();
        durable_json(journal, s);
    }
};
} // namespace
int main(int argc, char **argv) {
    if (argc == 3 && std::string(argv[1]) == "--fixture-stubs") {
        const auto destination = path_from_utf8(argv[2]);
        fs::create_directories(destination);
#ifdef _WIN32
        const auto binaries = fs::absolute(argv[0]).parent_path();
        fs::copy_file(binaries / "brlinst-stub.exe", destination / "brlinst-stub.exe");
        fs::copy_file(binaries / "brlinst-uninstall.exe", destination / "brlinst-uninstall.exe");
#else
        write(destination / "brlinst-stub.exe", stub());
        write(destination / "brlinst-uninstall.exe", stub());
#endif
        return 0;
    }
    auto root = fs::temp_directory_path() / ("installer-test-" + std::to_string(std::random_device{}()));
    fs::create_directories(root / "payload" / "bin");
    try {
        for (auto path : {"../x", "/tmp/x", "a/../b", "a\\b", "C:/x", "a:stream", "a//b", "NUL.txt", "a/b.",
                          ".brlinst/state.json"})
            fails([&] { validate_relative(path); }, "unsafe filename accepted");
        validate_relative("bin/program.exe");
        validate_relative("share/模型.txt");
        fails([] { parse_json("{\"a\":1,\"a\":2}"); }, "duplicate JSON accepted");
        require(parse_options({"--quiet", "--scope", "machine", "--path", "on"}).machine, "scope parsing");
        fails([] { parse_options({"--quiet", "--quiet"}); }, "duplicate switch accepted");
        fails([] { parse_options({"--path", "yes"}); }, "invalid boolean accepted");
        fails([] { parse_options({"--scope", "user"}, true); }, "uninstaller accepts scope override");
#ifdef _WIN32
        auto binaries = fs::absolute(argv[0]).parent_path();
        fs::copy_file(binaries / "brlinst-stub.exe", root / "stub.exe");
        fs::copy_file(binaries / "brlinst-uninstall.exe", root / "uninstall.exe");
#else
        write(root / "stub.exe", stub());
        write(root / "uninstall.exe", stub());
#endif
        write(root / "license.txt", "Test license\n");
        const std::string contents(200000, 'x');
        write(root / "payload/bin/app.exe", contents);
        write(root / "payload/empty", "");
        fs::create_directories(root / "payload/share/doc");
        write(root / "payload/share/doc/manual.txt", "Nested payload");
        Json config{{"formatVersion", 1},
                    {"id", "test.product"},
                    {"name", "Test Product"},
                    {"version", "1.2.3"},
                    {"architecture", "x64"},
                    {"publisher", "Tests"},
                    {"license", "license.txt"},
                    {"shortcuts", Json::array({{{"name", "Test"}, {"target", "bin/app.exe"}}})}};
        auto branded_config = config;
        branded_config["branding"] = {
            {"style", "classic"},
            {"exe", {{"header", "header.bmp"}, {"sidebar", "sidebar.bmp"}}},
            {"msi",
             {{"icon", "product.ico"},
              {"banner", "banner.bmp"},
              {"dialog", "dialog.bmp"},
              {"license", "license.rtf"}}}};
        validate_product(branded_config);
        require(product_architecture(config) == Architecture::x64, "x64 architecture not parsed");
        auto arm_config = config;
        arm_config["architecture"] = "arm64";
        validate_product(arm_config);
        require(product_architecture(arm_config) == Architecture::arm64, "arm64 architecture not parsed");
        require(instance_id(config) != instance_id(arm_config), "Architectures share an instance identity");
        auto missing_architecture = config;
        missing_architecture.erase("architecture");
        fails([&] { validate_product(missing_architecture); }, "Missing architecture accepted");
        auto invalid_architecture = config;
        invalid_architecture["architecture"] = "amd64";
        fails([&] { validate_product(invalid_architecture); }, "Invalid architecture accepted");
        auto incomplete_branding = branded_config;
        incomplete_branding["branding"]["exe"].erase("sidebar");
        fails([&] { validate_product(incomplete_branding); }, "Incomplete EXE branding accepted");
        auto unknown_branding = branded_config;
        unknown_branding["branding"]["style"] = "custom";
        fails([&] { validate_product(unknown_branding); }, "Unknown branding style accepted");
        write(root / "product.json", config.dump());
        create_package(root / "product.json", root / "payload", root / "stub.exe", root / "uninstall.exe",
                       root / "setup.exe");
        create_package(root / "product.json", root / "payload", root / "stub.exe", root / "uninstall.exe",
                       root / "repeat.exe");
        require(read_text(root / "setup.exe") == read_text(root / "repeat.exe"),
                "package generation is not deterministic");
#ifndef _WIN32
        write(root / "product-arm64.json", arm_config.dump());
        write(root / "stub-arm64.exe", stub(pe_machine_arm64));
        write(root / "uninstall-arm64.exe", stub(pe_machine_arm64));
        create_package(root / "product-arm64.json", root / "payload", root / "stub-arm64.exe",
                       root / "uninstall-arm64.exe", root / "setup-arm64.exe");
        auto arm_package = open_package(root / "setup-arm64.exe");
        verify_package(arm_package);
        require(product_architecture(arm_package.product) == Architecture::arm64,
                "ARM64 package architecture was not retained");
        fails(
            [&] {
                create_package(root / "product-arm64.json", root / "payload", root / "stub.exe",
                               root / "uninstall-arm64.exe", root / "mismatched-arm64.exe");
            },
            "Mismatched setup architecture accepted");
#endif
        auto package = open_package(root / "setup.exe");
        verify_package(package);
        require(package.files.size() == 5, "wrong packaged files");
        auto signed_data = read_text(root / "setup.exe");
        std::uint32_t pe = 0;
        for (unsigned i = 0; i < 4; ++i)
            pe |= static_cast<std::uint32_t>(static_cast<unsigned char>(signed_data[0x3c + i])) << (8 * i);
        put(signed_data, pe + 24 + 112 + 32, signed_data.size(), 4);
        put(signed_data, pe + 24 + 112 + 36, 16, 4);
        signed_data.append(16, 's');
        write(root / "signed.exe", signed_data);
        verify_package(open_package(root / "signed.exe"));
        auto corrupt = read_text(root / "setup.exe");
        corrupt.at(static_cast<std::size_t>(package.payload_offset) + 8) ^= 0x3f;
        write(root / "corrupt.exe", corrupt);
        fails([&] { verify_package(open_package(root / "corrupt.exe")); }, "corrupt payload accepted");
        auto truncated = read_text(root / "setup.exe");
        truncated.resize(truncated.size() - 2);
        write(root / "truncated.exe", truncated);
        fails([&] { open_package(root / "truncated.exe"); }, "truncated package accepted");
        Backend backend;
        Options options;
        options.directory = root / "install";
        install(package, options, backend);
        require(read_text(root / "install/bin/app.exe") == contents, "payload differs after install");
        require(load_state(options.directory)["phase"] == "installed", "install state not committed");
        fails([&] { install(package, options, backend); }, "nonempty destination accepted");
        write(root / "install/user-file", "Keep me");
        write(root / "install/bin/app.exe", "Modified");
        auto kept = uninstall(options.directory, backend);
        require(kept.size() == 1, "modified payload not reported");
        require(fs::exists(root / "install/user-file") && fs::exists(root / "install/bin/app.exe"),
                "user data removed");
        require(!fs::exists(root / "install/empty"), "owned file not removed");
        require(!fs::exists(root / "install/share"), "nested owned directories not removed");
        options.directory = root / "rollback";
        backend.fail = true;
        fails([&] { install(package, options, backend); }, "backend failure ignored");
        require(!fs::exists(options.directory), "rollback left owned payload");
        backend.fail = false;
        options.directory = root / "cancel";
        bool writing = false;
        fails(
            [&] {
                install(
                    package, options, backend,
                    [&](const std::string &m, std::uint64_t, std::uint64_t) {
                        if (m.rfind("Installing ", 0) == 0)
                            writing = true;
                    },
                    [&] { return writing; });
            },
            "cancel ignored");
        require(!fs::exists(options.directory), "cancel left destination");
        options.directory = root / "recovery";
        install(package, options, backend);
        backend.fail_remove = true;
        fails([&] { uninstall(options.directory, backend); }, "removal failure ignored");
        require(fs::exists(options.directory / state_directory / uninstall_filename),
                "failed removal lost uninstaller");
        require(fs::exists(options.directory / "bin/app.exe"), "failed integration removal deleted payload");
        backend.fail_remove = false;
        auto state = load_state(options.directory);
        state["phase"] = "installing";
        durable_json(options.directory / state_directory / state_filename, state);
        uninstall(options.directory, backend);
        require(!fs::exists(options.directory), "interrupted install recovery failed");
        // On ordinary Windows volumes this spelling aliases the existing file;
        // creating it would overwrite that file, not create a collision fixture.
        if (!fs::exists(root / "payload/bin/APP.EXE")) {
            write(root / "payload/bin/APP.EXE", "collision");
            fails(
                [&] {
                    create_package(root / "product.json", root / "payload", root / "stub.exe",
                                   root / "uninstall.exe", root / "collision.exe");
                },
                "case collision accepted");
            fs::remove(root / "payload/bin/APP.EXE");
        }
#ifndef _WIN32
        auto outside = root / "outside";
        write(outside, "untouched");
        fs::create_symlink(outside, root / "unsafe.exe.building");
        fails(
            [&] {
                create_package(root / "product.json", root / "payload", root / "stub.exe",
                               root / "uninstall.exe", root / "unsafe.exe");
            },
            "output symlink accepted");
        require(read_text(outside) == "untouched", "output symlink changed an unrelated file");
        fs::create_directory_symlink(root / "install", root / "link");
        fails([&] { check_no_links(root / "link/file"); }, "symlink destination accepted");
#endif
        fs::remove_all(root);
        std::cout << "All core scenarios passed\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << "\nTest directory: " << root << '\n';
        return 1;
    }
}
