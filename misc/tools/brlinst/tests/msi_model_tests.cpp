// SPDX-License-Identifier: BSD-2-Clause
#include "installer/msi.hpp"
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <set>
using namespace installer;
using namespace installer::msi;
namespace {
void require(bool condition, const std::string &message) {
    if (!condition)
        throw Error(message);
}
template <class F> void rejects(F f) {
    bool rejected = false;
    try {
        f();
    } catch (const std::exception &) {
        rejected = true;
    }
    require(rejected, "Invalid MSI input was accepted");
}
void write(const fs::path &p, const std::string &text) {
    std::ofstream(p, std::ios::binary) << text;
}
const Json &by_key(Model &m, const char *table, const std::string &key) {
    for (const auto &row : m.table(table).rows)
        if (row[0] == key)
            return row;
    throw Error("Missing MSI row: " + key);
}
// Check relational integrity independently of the generator. This catches
// broken control navigation, components, file keys and schema changes on Unix.
void relational(Model &m) {
    for (const auto &table : m.tables) {
        std::set<std::string> keys;
        for (const auto &row : table.rows) {
            require(row.size() == table.columns.size(), "Row width: " + table.name);
            Json key = Json::array();
            for (std::size_t i = 0; i < table.columns.size(); ++i) {
                const auto &column = table.columns[i];
                const auto &value = row[i];
                if (column.key)
                    key.push_back(value);
                require(column.nullable || !value.is_null(),
                        "Missing field: " + table.name + "." + column.name);
                if (value.is_null())
                    continue;
                if (column.type == "SHORT" || column.type == "LONG") {
                    require(value.is_number_integer(), "Non-integer MSI field");
                    if (column.type == "SHORT")
                        require(value.get<int>() >= -32767 && value.get<int>() <= 32767, "SHORT overflow");
                } else if (column.type == "OBJECT")
                    require(value.is_object() && value.contains("stream"), "Invalid binary stream");
                else {
                    require(value.is_string(), "Non-string MSI field");
                    auto length = std::stoi(column.type.substr(5));
                    if (length)
                        require(value.get<std::string>().size() <= static_cast<std::size_t>(length),
                                "Oversized field: " + table.name + "." + column.name);
                }
                if (!column.foreign_table.empty()) {
                    const auto &foreign = m.table(column.foreign_table);
                    bool found = false;
                    for (const auto &candidate : foreign.rows) {
                        bool matches = candidate[column.foreign_column - 1] == value;
                        for (int prefix = 0; prefix < column.foreign_column - 1; ++prefix)
                            matches = matches && candidate[prefix] == row[prefix];
                        found = found || matches;
                    }
                    require(found, "Broken foreign key: " + table.name + "." + column.name);
                }
            }
            require(keys.insert(key.dump()).second, "Duplicate primary key: " + table.name);
        }
    }
    for (const auto &row : m.table("Component").rows) {
        const auto table = row[3].get<int>() & 4 ? "Registry" : "File";
        by_key(m, table, row[5]);
    }
    for (const auto &row : m.table("ControlEvent").rows) {
        if (row[2] == "NewDialog" || row[2] == "SpawnDialog")
            by_key(m, "Dialog", row[3]);
        if (row[2] == "DoAction")
            by_key(m, "CustomAction", row[3]);
    }
    for (const auto &d : m.table("Dialog").rows) {
        std::map<std::string, std::string> next;
        for (const auto &c : m.table("Control").rows)
            if (c[0] == d[0] && !c[10].is_null())
                next[c[1]] = c[10];
        for (const auto &c : next) {
            std::set<std::string> visited;
            auto cursor = c.first;
            while (visited.insert(cursor).second) {
                require(next.count(cursor), "Tab navigation leads to an unfocusable control");
                cursor = next.at(cursor);
            }
            require(cursor == c.first, "Tab navigation is not a closed cycle");
        }
        for (int field : {7, 8, 9})
            if (!d[field].is_null()) {
                bool found = false;
                for (const auto &c : m.table("Control").rows)
                    found = found || (c[0] == d[0] && c[1] == d[field]);
                require(found, "Missing dialog focus/default/cancel control");
            }
    }
}
} // namespace
int main() {
    auto root = fs::temp_directory_path() / ("brlinst-msi-model-" + std::to_string(std::random_device{}()));
    fs::create_directories(root / "payload/bin");
    try {
        write(root / "license.txt", "Test license {with} [brackets] and Unicode: 模型\n");
        write(root / "payload/bin/app.exe", "Test executable");
        write(root / "payload/data.txt", "Test data");
        fs::create_directories(root / "payload/share/doc");
        write(root / "payload/share/doc/manual.txt", "Nested documentation");
        write(root / "NOTICE.txt", "Notices");
        write(root / "actions.dll", "Test DLL");
        Json config{{"formatVersion", 1},
                    {"id", "test.msi"},
                    {"name", "MSI Test"},
                    {"version", "1.2.3"},
                    {"publisher", "Test"},
                    {"license", "license.txt"},
                    {"shortcuts", Json::array({{{"name", "Application"},
                                                {"target", "bin/app.exe"},
                                                {"arguments", "--title [literal]"}},
                                               {{"name", "Data"},
                                                {"target", "data.txt"},
                                                {"desktop", false},
                                                {"workingDirectory", "."}}})}};
        write(root / "product.json", config.dump());
        auto input = prepare_input(root / "product.json", root / "payload", root / "test.msi");
        // A deterministic stand-in lets tests check the complete identity inputs
        // without requiring Windows CNG. Windows tests check real GUIDs as well.
        std::map<std::string, std::string> ids;
        auto guid = [&](const std::string &key) {
            auto inserted = ids.emplace(key, "");
            if (inserted.second) {
                auto value = std::to_string(ids.size());
                inserted.first->second =
                    "{00000000-0000-8000-8000-" + std::string(12 - value.size(), '0') + value + "}";
            }
            return inserted.first->second;
        };
        auto model = make_model(input, root / "NOTICE.txt", root / "actions.dll", guid);
        relational(model);
        for (const auto &file : model.files)
            require(file.path.find('\\') == std::string::npos,
                    "Windows separators leaked into payload manifest");
        require(model.files.size() == 4, "Missing payload or notices");
        require(model.table("Shortcut").rows.size() == 3, "Desktop eligibility ignored");
        require(by_key(model, "Property", "ALLUSERS")[1] == "2", "Dual-purpose package missing");
        require(by_key(model, "Property", "MSIINSTALLPERUSER")[1] == "1", "Per-user default missing");
        bool automatic_upgrade = false;
        for (const auto &r : model.table("InstallExecuteSequence").rows)
            automatic_upgrade = automatic_upgrade || r[0] == "RemoveExistingProducts";
        require(!automatic_upgrade, "Automatic upgrading was authored");
        require(by_key(model, "InstallExecuteSequence", "BrlinstCheck")[2].get<int>() <
                    by_key(model, "InstallExecuteSequence", "InstallInitialize")[2].get<int>(),
                "Destination preflight runs after mutation starts");
        auto old = model;
        input.product["version"] = "2.0.0";
        auto newer = make_model(input, root / "NOTICE.txt", root / "actions.dll", guid);
        require(newer.product_code != old.product_code, "Versions share a ProductCode");
        std::set<std::string> old_components;
        for (const auto &r : old.table("Component").rows)
            old_components.insert(r[1]);
        for (const auto &r : newer.table("Component").rows)
            require(!old_components.count(r[1]), "Versions share a component GUID");
        input.product["version"] = "1.2.3";
        input.files.insert(input.files.begin(), {"aardvark.txt", root / "payload/data.txt"});
        auto changed = make_model(input, root / "NOTICE.txt", root / "actions.dll", guid);
        require(changed.product_code == old.product_code, "File list changes ProductCode");
        auto component_guid = [](Model &m, const std::string &path) {
            for (const auto &f : m.files)
                if (f.path == path) {
                    auto component = by_key(m, "File", f.id)[1];
                    return by_key(m, "Component", component)[1];
                }
            throw Error("Missing test file");
        };
        require(component_guid(changed, "bin/app.exe") == component_guid(old, "bin/app.exe"),
                "Component identity depends on enumeration order");
        for (auto v : {"1", "1.2", "1.2.3.4", "01.2.3", "256.1.1", "1.256.1", "1.1.65536", "1.2.3-beta"}) {
            input.product["version"] = v;
            rejects([&] { make_model(input, root / "NOTICE.txt", root / "actions.dll", guid); });
        }
        require(formatted_literal("a[b]c") == "a[\\[]b[\\]]c", "Literal MSI property escaping failed");
        fs::remove_all(root);
        std::cout << "MSI model tests passed\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        fs::remove_all(root);
        return 1;
    }
}
