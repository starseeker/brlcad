// SPDX-License-Identifier: BSD-2-Clause
#include "installer/core.hpp"
#include "notices.hpp"
#include <algorithm>
#include <cctype>
#include <set>
#ifdef _WIN32
#include <windows.h>
#endif
namespace installer {
std::string license_notices() {
    return std::string(reinterpret_cast<const char *>(distribution_notices));
}
PackageInput prepare_input(const fs::path &config, const fs::path &source, const fs::path &output) {
    PackageInput result;
    result.product = parse_json(read_text(config));
    auto &product = result.product;
    validate_product(product);
    auto asset = [&](const std::string &name) {
        auto p = path_from_utf8(name);
        return p.is_absolute() ? p : config.parent_path() / p;
    };
    product["licenseText"] = read_text(asset(product.at("license").get<std::string>()), 1024 * 1024);
    product.erase("license");
    if (product.contains("icon")) {
        result.icon = asset(product["icon"].get<std::string>());
        check_no_links(result.icon);
        if (!fs::is_regular_file(result.icon))
            throw Error("Missing product icon");
        product.erase("icon");
    }
    validate_product(product, true);
    check_no_links(source);
    if (!fs::is_directory(source))
        throw Error("Source must be a staged directory");
    for (const auto &s : product["shortcuts"]) {
        for (auto key : {"target", "icon"})
            if (s.contains(key)) {
                auto path = source / path_from_utf8(s[key]);
                check_no_links(path);
                if (!fs::is_regular_file(path))
                    throw Error("Missing shortcut file: " + utf8(path));
            }
        auto dir = s.value("workingDirectory", "@profile");
        if (dir != "@profile" && dir != ".") {
            check_no_links(source / path_from_utf8(dir));
            if (!fs::is_directory(source / path_from_utf8(dir)))
                throw Error("Missing shortcut working directory");
        }
    }
    check_no_links(output);
    if (fs::exists(output))
        throw Error("Output already exists: " + utf8(output));
    auto relative =
        fs::absolute(output).lexically_normal().lexically_relative(fs::absolute(source).lexically_normal());
    if (!relative.empty() && *relative.begin() != "..")
        throw Error("Output must be outside source tree");
    // Include directories in the collision set, so differently spelled parent
    // directories cannot merge silently when a Unix tree is installed on Windows.
    std::set<std::string> names;
    for (const auto &entry : fs::recursive_directory_iterator(source)) {
        check_no_links(entry.path());
        auto name = entry.path().lexically_relative(source).generic_u8string();
        validate_relative(name);
        auto key = name;
#ifdef _WIN32
        auto wide = path_from_utf8(key).wstring();
        CharUpperBuffW(wide.data(), static_cast<DWORD>(wide.size()));
        key = utf8(fs::path(wide));
#else
        std::transform(key.begin(), key.end(), key.begin(),
                       [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
#endif
        if (!names.insert(key).second)
            throw Error("Duplicate Windows filename: " + name);
        if (entry.is_directory())
            continue;
        if (!entry.is_regular_file())
            throw Error("Unsupported payload file type");
        result.files.emplace_back(name, entry.path());
        if (result.files.size() > 199998)
            throw Error("Too many payload files");
    }
    std::sort(result.files.begin(), result.files.end());
    return result;
}
} // namespace installer
