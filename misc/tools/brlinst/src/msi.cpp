// SPDX-License-Identifier: BSD-2-Clause
#include "installer/msi.hpp"
#include <windows.h>
#include <msi.h>
#include <msiquery.h>
#include <bcrypt.h>
#include <fci.h>
#include <objbase.h>
#include <algorithm>
#include <array>
#include <cstdlib>
#include <cerrno>
#include <climits>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
namespace installer::msi {
namespace {
struct Handle {
    MSIHANDLE value = 0;
    Handle() = default;
    explicit Handle(MSIHANDLE h) : value(h) {}
    Handle(const Handle &) = delete;
    Handle &operator=(const Handle &) = delete;
    ~Handle() {
        if (value)
            MsiCloseHandle(value);
    }
    operator MSIHANDLE() const {
        return value;
    }
};
std::wstring wide(const std::string &s) {
    return path_from_utf8(s).wstring();
}
void check(UINT result, const std::string &what) {
    if (result == ERROR_SUCCESS)
        return;
    std::string detail;
    Handle error(MsiGetLastErrorRecord());
    if (error.value) {
        DWORD length = 0;
        wchar_t empty = 0;
        if (MsiFormatRecordW(0, error, &empty, &length) == ERROR_MORE_DATA) {
            std::wstring text(length + 1, L'\0');
            if (MsiFormatRecordW(0, error, text.data(), &length) == ERROR_SUCCESS) {
                text.resize(length);
                detail = ": " + utf8(fs::path(text));
            }
        }
    }
    throw Error(what + " (Windows Installer " + std::to_string(result) + ")" + detail);
}
std::wstring record_string(MSIHANDLE record, UINT field) {
    DWORD length = 0;
    wchar_t empty = 0;
    auto result = MsiRecordGetStringW(record, field, &empty, &length);
    if (result != ERROR_MORE_DATA && result != ERROR_SUCCESS)
        check(result, "Read MSI field");
    std::wstring text(length + 1, L'\0');
    check(MsiRecordGetStringW(record, field, text.data(), &length), "Read MSI field");
    text.resize(length);
    return text;
}
void query(MSIHANDLE db, const std::string &sql, MSIHANDLE parameters = 0) {
    Handle view;
    check(MsiDatabaseOpenViewW(db, wide(sql).c_str(), &view.value), "Prepare " + sql);
    check(MsiViewExecute(view, parameters), "Execute " + sql);
}
std::string guid_text(const std::array<unsigned char, 16> &bytes) {
    std::ostringstream out;
    out << '{' << std::uppercase << std::hex << std::setfill('0');
    for (std::size_t i = 0; i < bytes.size(); ++i) {
        if (i == 4 || i == 6 || i == 8 || i == 10)
            out << '-';
        out << std::setw(2) << static_cast<unsigned>(bytes[i]);
    }
    return out.str() + '}';
}
std::string stable_guid(const std::string &name) {
    // UUIDv8: first 128 SHA-256 bits with the version and variant bits set.
    // This identity scheme is part of the backend's v1 contract; do not change it.
    BCRYPT_ALG_HANDLE alg = nullptr;
    if (BCryptOpenAlgorithmProvider(&alg, BCRYPT_SHA256_ALGORITHM, nullptr, 0) < 0)
        throw Error("Cannot initialize SHA-256 identity generator");
    struct Close {
        BCRYPT_ALG_HANDLE h;
        ~Close() {
            BCryptCloseAlgorithmProvider(h, 0);
        }
    } close{alg};
    std::array<unsigned char, 32> hash{};
    if (BCryptHash(alg, nullptr, 0, reinterpret_cast<PUCHAR>(const_cast<char *>(name.data())),
                   static_cast<ULONG>(name.size()), hash.data(), static_cast<ULONG>(hash.size())) < 0)
        throw Error("Cannot generate MSI identity");
    std::array<unsigned char, 16> id{};
    std::copy_n(hash.begin(), id.size(), id.begin());
    id[6] = static_cast<unsigned char>((id[6] & 15) | 0x80);
    id[8] = static_cast<unsigned char>((id[8] & 63) | 0x80);
    return guid_text(id);
}
std::string random_guid() {
    GUID guid{};
    if (FAILED(CoCreateGuid(&guid)))
        throw Error("Cannot create package GUID");
    wchar_t text[40]{};
    if (!StringFromGUID2(guid, text, 40))
        throw Error("Cannot format package GUID");
    return utf8(fs::path(text));
}
void representable(const std::wstring &value) {
    if (value.empty())
        return;
    BOOL substituted = FALSE;
    if (!WideCharToMultiByte(1252, WC_NO_BEST_FIT_CHARS, value.data(), static_cast<int>(value.size()),
                             nullptr, 0, nullptr, &substituted) ||
        substituted)
        throw Error(
            "MSI v1 requires Windows-1252 product text and payload names; EXE supports broader Unicode");
}
void insert(MSIHANDLE db, const Table &table, const Json &row) {
    if (row.size() != table.columns.size())
        throw Error("Invalid MSI row width: " + table.name);
    Handle record(MsiCreateRecord(static_cast<UINT>(row.size())));
    if (!record.value)
        throw Error("Cannot allocate MSI record");
    std::string columns, placeholders;
    for (std::size_t i = 0; i < row.size(); ++i) {
        if (i) {
            columns += ',';
            placeholders += ',';
        }
        columns += '`' + table.columns[i].name + '`';
        placeholders += '?';
        auto field = static_cast<UINT>(i + 1);
        const auto &value = row[i];
        if (value.is_null())
            continue;
        if (value.is_number_integer()) {
            auto integer = value.get<std::int64_t>();
            if (integer <= INT_MIN || integer > INT_MAX)
                throw Error("MSI integer out of range");
            check(MsiRecordSetInteger(record, field, static_cast<int>(integer)), "Set MSI integer");
        } else if (value.is_object()) {
            auto path = path_from_utf8(value.at("stream").get<std::string>());
            check_no_links(path);
            check(MsiRecordSetStreamW(record, field, path.c_str()), "Set MSI stream");
        } else {
            auto text = wide(value.get<std::string>());
            representable(text);
            check(MsiRecordSetStringW(record, field, text.c_str()), "Set MSI string");
        }
    }
    query(db, "INSERT INTO `" + table.name + "` (" + columns + ") VALUES (" + placeholders + ")", record);
}
void write(const fs::path &path, const std::string &text) {
    std::ofstream out(path, std::ios::binary);
    if (!out.write(text.data(), static_cast<std::streamsize>(text.size())))
        throw Error("Cannot write " + utf8(path));
}
struct Cabinet {
    fs::path work;
    std::map<std::string, fs::path> sources;
    unsigned temporary = 0;
    fs::path path(const char *name) {
        auto it = sources.find(name);
        if (it != sources.end())
            return it->second;
        auto relative = path_from_utf8(name);
        if (relative.has_parent_path() || relative.empty())
            throw Error("Unexpected cabinet work filename");
        return work / relative;
    }
};
void *DIAMONDAPI cab_alloc(ULONG size) {
    return std::malloc(size);
}
void DIAMONDAPI cab_free(void *p) {
    std::free(p);
}
INT_PTR DIAMONDAPI cab_open(char *name, int flags, int, int *error, void *context) {
    try {
        auto &cab = *static_cast<Cabinet *>(context);
        auto path = cab.path(name);
        check_no_links(path);
        const bool writing = (flags & (_O_WRONLY | _O_RDWR)) != 0;
        // Source payloads are immutable inputs even if FCI asks unexpectedly.
        if (writing && cab.sources.count(name)) {
            *error = EACCES;
            return -1;
        }
        auto handle = CreateFileW(path.c_str(), writing ? GENERIC_READ | GENERIC_WRITE : GENERIC_READ,
                                  FILE_SHARE_READ, nullptr, flags & _O_CREAT ? CREATE_ALWAYS : OPEN_EXISTING,
                                  FILE_ATTRIBUTE_NORMAL, nullptr);
        if (handle == INVALID_HANDLE_VALUE) {
            *error = static_cast<int>(GetLastError());
            return -1;
        }
        return reinterpret_cast<INT_PTR>(handle);
    } catch (...) {
        *error = EIO;
        return -1;
    }
}
UINT DIAMONDAPI cab_read(INT_PTR file, void *data, UINT count, int *error, void *) {
    DWORD done = 0;
    if (!ReadFile(reinterpret_cast<HANDLE>(file), data, count, &done, nullptr)) {
        *error = static_cast<int>(GetLastError());
        return static_cast<UINT>(-1);
    }
    return done;
}
UINT DIAMONDAPI cab_write(INT_PTR file, void *data, UINT count, int *error, void *) {
    DWORD done = 0;
    if (!WriteFile(reinterpret_cast<HANDLE>(file), data, count, &done, nullptr) || done != count) {
        *error = static_cast<int>(GetLastError());
        return static_cast<UINT>(-1);
    }
    return done;
}
int DIAMONDAPI cab_close(INT_PTR file, int *error, void *) {
    if (CloseHandle(reinterpret_cast<HANDLE>(file)))
        return 0;
    *error = static_cast<int>(GetLastError());
    return -1;
}
long DIAMONDAPI cab_seek(INT_PTR file, long distance, int origin, int *error, void *) {
    LARGE_INTEGER offset{}, result{};
    offset.QuadPart = distance;
    if (!SetFilePointerEx(reinterpret_cast<HANDLE>(file), offset, &result, static_cast<DWORD>(origin)) ||
        result.QuadPart > LONG_MAX) {
        *error = EIO;
        return -1;
    }
    return static_cast<long>(result.QuadPart);
}
int DIAMONDAPI cab_delete(char *name, int *error, void *context) {
    try {
        auto &cab = *static_cast<Cabinet *>(context);
        if (cab.sources.count(name)) {
            *error = EACCES;
            return -1;
        }
        if (DeleteFileW(cab.path(name).c_str()))
            return 0;
    } catch (...) {
    }
    *error = EIO;
    return -1;
}
BOOL DIAMONDAPI cab_temp(char *name, int capacity, void *context) {
    try {
        auto &cab = *static_cast<Cabinet *>(context);
        auto text = "scratch" + std::to_string(++cab.temporary);
        if (text.size() + 1 > static_cast<std::size_t>(capacity))
            return FALSE;
        std::memcpy(name, text.c_str(), text.size() + 1);
        return TRUE;
    } catch (...) {
        return FALSE;
    }
}
int DIAMONDAPI cab_placed(PCCAB, char *, long, BOOL, void *) {
    return 0;
}
BOOL DIAMONDAPI cab_next(PCCAB, ULONG, void *) {
    return FALSE;
}
long DIAMONDAPI cab_status(UINT, ULONG, ULONG, void *) {
    return 0;
}
INT_PTR DIAMONDAPI cab_info(char *name, USHORT *date, USHORT *time, USHORT *attributes, int *error,
                            void *context) {
    *date = 33; // 1980-01-01, deterministic DOS timestamp
    *time = 0;
    *attributes = FILE_ATTRIBUTE_ARCHIVE;
    return cab_open(name, _O_RDONLY, 0, error, context);
}
void cabinet(const Model &model, const fs::path &work) {
    Cabinet context{work, {}, 0};
    CCAB parameters{};
    parameters.cb = 0x7fffffff;
    parameters.cbFolderThresh = 64 * 1024 * 1024;
    parameters.setID = 1;
    std::strcpy(parameters.szCab, "payload.cab");
    ERF error{};
    auto handle = FCICreate(&error, cab_placed, cab_alloc, cab_free, cab_open, cab_read, cab_write, cab_close,
                            cab_seek, cab_delete, cab_temp, &parameters, &context);
    if (!handle)
        throw Error("Cannot create cabinet context");
    struct Destroy {
        HFCI h;
        ~Destroy() {
            FCIDestroy(h);
        }
    } destroy{handle};
    for (const auto &file : model.files) {
        auto token = "source" + file.id;
        auto member = file.id;
        context.sources.emplace(token, file.source);
        if (!FCIAddFile(handle, token.data(), member.data(), FALSE, cab_next, cab_status, cab_info,
                        tcompTYPE_MSZIP))
            throw Error("Cannot compress " + file.path + " (FCI " + std::to_string(error.erfOper) + ", " +
                        std::to_string(error.erfType) + ")");
    }
    if (!FCIFlushCabinet(handle, FALSE, cab_next, cab_status))
        throw Error("Cannot finish cabinet");
}
void file_metadata(Model &model) {
    auto &rows = model.table("File").rows;
    for (std::size_t i = 0; i < model.files.size(); ++i) {
        const auto &file = model.files[i];
        wchar_t version[256]{}, language[256]{};
        DWORD vn = 256, ln = 256;
        auto status = MsiGetFileVersionW(file.source.c_str(), version, &vn, language, &ln);
        if (status == ERROR_SUCCESS && vn) {
            rows[i][4] = utf8(fs::path(version));
            if (ln)
                rows[i][5] = utf8(fs::path(language));
        } else if (status == ERROR_FILE_INVALID || status == ERROR_FILE_NOT_FOUND ||
                   (status == ERROR_SUCCESS && !vn)) {
            MSIFILEHASHINFO hash{};
            hash.dwFileHashInfoSize = sizeof(hash);
            check(MsiGetFileHashW(file.source.c_str(), 0, &hash), "Hash " + file.path);
            model.table("MsiFileHash")
                .rows.push_back({file.id, 0, static_cast<std::int32_t>(hash.dwData[0]),
                                 static_cast<std::int32_t>(hash.dwData[1]),
                                 static_cast<std::int32_t>(hash.dwData[2]),
                                 static_cast<std::int32_t>(hash.dwData[3])});
        } else
            check(status, "Read file version: " + file.path);
    }
}
void summary(MSIHANDLE db, const Model &model, const Json &product) {
    Handle info;
    check(MsiGetSummaryInformationW(db, nullptr, 20, &info.value), "Open MSI summary");
    auto text = [&](UINT id, const std::string &value) {
        auto w = wide(value);
        representable(w);
        check(MsiSummaryInfoSetPropertyW(info, id, VT_LPSTR, 0, nullptr, w.c_str()),
              "Set MSI summary string");
    };
    auto number = [&](UINT id, int value, UINT type = VT_I4) {
        check(MsiSummaryInfoSetPropertyW(info, id, type, value, nullptr, nullptr), "Set MSI summary number");
    };
    number(1, 1252, VT_I2); // PID_CODEPAGE
    text(2, "Installation Database");
    text(3, model.name);
    text(4, product.at("publisher"));
    text(5, "Installer");
    text(7, "x64;1033");
    text(9, random_guid());
    number(14, 500); // Windows Installer 5.0, needed for dual-purpose authoring
    number(15, 10);  // Compressed payload and LUA (per-user-capable) package
    text(18, "brlinst");
    number(19, 2); // Read-only recommended
    check(MsiSummaryInfoPersist(info), "Persist MSI summary");
}
} // namespace
void verify(const fs::path &package) {
    check(MsiVerifyPackageW(package.c_str()), "Verify MSI storage");
    Handle db, view;
    check(MsiOpenDatabaseW(package.c_str(), MSIDBOPEN_READONLY, &db.value), "Open generated MSI");
    check(MsiDatabaseOpenViewW(db, L"SELECT `Name` FROM `_Tables`", &view.value), "Read MSI table list");
    check(MsiViewExecute(view, 0), "Query MSI tables");
    for (;;) {
        Handle record;
        auto result = MsiViewFetch(view, &record.value);
        if (result == ERROR_NO_MORE_ITEMS)
            break;
        check(result, "Fetch MSI table");
        auto name = record_string(record, 1);
        if (name.empty() || name.front() == L'_')
            continue;
        Handle table;
        check(MsiDatabaseOpenViewW(db, (L"SELECT * FROM `" + name + L"`").c_str(), &table.value),
              "Open MSI table");
        check(MsiViewExecute(table, 0), "Read MSI table");
        for (;;) {
            Handle item;
            result = MsiViewFetch(table, &item.value);
            if (result == ERROR_NO_MORE_ITEMS)
                break;
            check(result, "Fetch MSI row");
            auto validated = MsiViewModify(table, MSIMODIFY_VALIDATE, item);
            if (validated != ERROR_SUCCESS) {
                wchar_t column[256]{};
                DWORD length = 256;
                auto error = MsiViewGetErrorW(table, column, &length);
                throw Error("Invalid MSI row in " + utf8(fs::path(name)) + "." + utf8(fs::path(column)) +
                            " (database validation " + std::to_string(static_cast<int>(error)) + ")");
            }
        }
    }
    Handle streams;
    check(MsiDatabaseOpenViewW(db, L"SELECT `Data` FROM `_Streams` WHERE `Name` = 'payload.cab'",
                               &streams.value),
          "Find embedded cabinet");
    check(MsiViewExecute(streams, 0), "Read embedded cabinet");
    Handle payload;
    check(MsiViewFetch(streams, &payload.value), "Fetch embedded cabinet");
    char magic[4]{};
    DWORD count = 4;
    check(MsiRecordReadStream(payload, 1, magic, &count), "Read cabinet header");
    if (count != 4 || std::memcmp(magic, "MSCF", 4))
        throw Error("Missing MSI cabinet signature");
}
void create(const fs::path &config, const fs::path &source, const fs::path &output, const fs::path &actions) {
    auto input = prepare_input(config, source, output);
    check_no_links(actions);
    if (!fs::is_regular_file(actions))
        throw Error("Missing brlinst-msi-actions.dll beside the generator");
    auto work = output;
    work += ".building";
    check_no_links(work);
    if (!fs::create_directory(work))
        throw Error("Remove stale MSI work directory: " + utf8(work));
    struct Cleanup {
        fs::path p;
        ~Cleanup() {
            std::error_code ec;
            fs::remove_all(p, ec);
        }
    } cleanup{work};
    auto notices = work / "NOTICE.txt";
    write(notices, license_notices());
    auto model = make_model(input, notices, actions, stable_guid);
    // Hold read handles without write/delete sharing until both the metadata and
    // CAB have been written. A changing build tree must not produce a package
    // whose version/hash records describe different bytes from the payload.
    struct Locks {
        std::vector<HANDLE> files;
        ~Locks() {
            for (auto h : files)
                CloseHandle(h);
        }
    } locks;
    for (const auto &file : model.files) {
        auto handle = CreateFileW(file.source.c_str(), GENERIC_READ, FILE_SHARE_READ, nullptr, OPEN_EXISTING,
                                  FILE_ATTRIBUTE_NORMAL, nullptr);
        if (handle == INVALID_HANDLE_VALUE)
            throw Error("Cannot lock staged file: " + file.path);
        locks.files.push_back(handle);
        LARGE_INTEGER size{};
        if (!GetFileSizeEx(handle, &size) || static_cast<std::uint64_t>(size.QuadPart) != file.size)
            throw Error("Staged file changed during packaging: " + file.path);
    }
    file_metadata(model);
    cabinet(model, work);
    auto temporary = work / "package.msi";
    {
        Handle db;
        check(MsiOpenDatabaseW(temporary.c_str(), MSIDBOPEN_CREATE, &db.value), "Create MSI database");
        write(work / "codepage.idt", "\r\n\r\n1252\t_ForceCodepage\r\n");
        check(MsiDatabaseImportW(db, work.c_str(), L"codepage.idt"), "Set MSI database code page");
        for (const auto &t : model.tables)
            query(db, sql_create(t));
        for (const auto &t : model.tables)
            for (const auto &r : t.rows)
                insert(db, t, r);
        Table streams{"_Streams",
                      {{"Name", "CHAR(62)", "", false, true, false, "", 0},
                       {"Data", "OBJECT", "", false, false, false, "", 0}},
                      {}};
        insert(db, streams, Json{"payload.cab", {{"stream", utf8(work / "payload.cab")}}});
        summary(db, model, input.product);
        check(MsiDatabaseCommit(db), "Commit MSI database");
    }
    if (fs::file_size(temporary) >= 0x7fffffffULL)
        throw Error("MSI exceeds the 2 GiB package limit");
    verify(temporary);
    flush_file(temporary);
    publish_file(temporary, output);
}
} // namespace installer::msi
