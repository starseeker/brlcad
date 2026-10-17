// SPDX-License-Identifier: BSD-2-Clause
#include "installer/windows.hpp"
#include <commctrl.h>
#include <shellapi.h>
#include <array>
#include <sstream>

namespace installer {
namespace {
constexpr int id_back = 100, id_next = 101, id_cancel = 102, id_user = 103, id_machine = 104;
constexpr int id_directory = 105, id_browse = 106, id_start = 107, id_desktop = 108, id_path = 109;
constexpr int id_content = 110, id_progress = 111;
std::wstring control_text(HWND control) {
    int n = GetWindowTextLengthW(control);
    std::wstring text(static_cast<std::size_t>(n) + 1, L'\0');
    GetWindowTextW(control, text.data(), n + 1);
    text.resize(n);
    return text;
}
} // namespace
struct Wizard::Impl {
    Json product;
    Options &options;
    bool removing = false, cancel = false, ready = false, closed = false, complete = false, running = false;
    int page = 0;
    HWND hwnd = nullptr;
    HFONT font = nullptr, title_font = nullptr;
    HBRUSH background = nullptr;
    double scale = 1;
    std::vector<HWND> controls;
    std::wstring finished;
    Impl(const Json &p, Options &o, bool r) : product(p), options(o), removing(r) {
        INITCOMMONCONTROLSEX common{sizeof(common), ICC_PROGRESS_CLASS | ICC_STANDARD_CLASSES};
        InitCommonControlsEx(&common);
        WNDCLASSW wc{};
        wc.lpfnWndProc = procedure;
        wc.hInstance = GetModuleHandleW(nullptr);
        wc.lpszClassName = L"BrlinstWizard";
        wc.hCursor = LoadCursorW(nullptr, IDC_ARROW);
        wc.hbrBackground = reinterpret_cast<HBRUSH>(COLOR_WINDOW + 1);
        wc.hIcon = LoadIconW(wc.hInstance, MAKEINTRESOURCEW(1));
        if (!wc.hIcon)
            wc.hIcon = LoadIconW(nullptr, IDI_APPLICATION);
        RegisterClassW(&wc);
        auto screen = GetDC(nullptr);
        scale = static_cast<double>(GetDeviceCaps(screen, LOGPIXELSX)) / 96.0;
        ReleaseDC(nullptr, screen);
        font =
            CreateFontW(-px(14), 0, 0, 0, FW_NORMAL, FALSE, FALSE, FALSE, DEFAULT_CHARSET, OUT_DEFAULT_PRECIS,
                        CLIP_DEFAULT_PRECIS, CLEARTYPE_QUALITY, DEFAULT_PITCH, L"Segoe UI");
        title_font = CreateFontW(-px(23), 0, 0, 0, FW_SEMIBOLD, FALSE, FALSE, FALSE, DEFAULT_CHARSET,
                                 OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS, CLEARTYPE_QUALITY, DEFAULT_PITCH,
                                 L"Segoe UI");
        RECT rect{0, 0, px(700), px(500)};
        AdjustWindowRectEx(&rect, WS_OVERLAPPED | WS_CAPTION | WS_SYSMENU | WS_MINIMIZEBOX, FALSE, 0);
        auto title = wide(product.at("name")) + (r ? L" Uninstall" : L" Setup");
        hwnd = CreateWindowExW(0, wc.lpszClassName, title.c_str(),
                               WS_OVERLAPPED | WS_CAPTION | WS_SYSMENU | WS_MINIMIZEBOX, CW_USEDEFAULT,
                               CW_USEDEFAULT, rect.right - rect.left, rect.bottom - rect.top, nullptr,
                               nullptr, wc.hInstance, this);
        win_check(hwnd != nullptr, "Cannot create installer window");
        show_page();
        ShowWindow(hwnd, SW_SHOW);
        UpdateWindow(hwnd);
    }
    ~Impl() {
        if (hwnd && IsWindow(hwnd))
            DestroyWindow(hwnd);
        if (font)
            DeleteObject(font);
        if (title_font)
            DeleteObject(title_font);
        if (background)
            DeleteObject(background);
    }
    int px(int value) const {
        return static_cast<int>(value * scale);
    }
    HWND add(const wchar_t *cls, const std::wstring &text, int id, int x, int y, int w, int h,
             DWORD style = 0) {
        HWND child = CreateWindowExW(cls == std::wstring(L"EDIT") ? WS_EX_CLIENTEDGE : 0, cls, text.c_str(),
                                     WS_CHILD | WS_VISIBLE | style, px(x), px(y), px(w), px(h), hwnd,
                                     reinterpret_cast<HMENU>(static_cast<INT_PTR>(id)),
                                     GetModuleHandleW(nullptr), nullptr);
        win_check(child != nullptr, "Cannot create wizard control");
        SendMessageW(child, WM_SETFONT, reinterpret_cast<WPARAM>(font), TRUE);
        controls.push_back(child);
        return child;
    }
    void label(const std::wstring &text, int y, int h = 65) {
        add(L"STATIC", text, 0, 28, y, 640, h);
    }
    void show_page() {
        for (auto c : controls)
            DestroyWindow(c);
        controls.clear();
        auto heading = page == 0   ? (removing ? L"Remove this version" : L"Welcome")
                       : page == 1 ? L"License information"
                       : page == 2 ? L"Installation location"
                       : page == 3 ? L"Shortcuts and command line"
                       : page == 4 ? L"Ready to install"
                       : page == 5 ? (removing ? L"Removing files" : L"Installing")
                                   : L"Finished";
        auto title = add(L"STATIC", heading, 0, 28, 22, 640, 38);
        SendMessageW(title, WM_SETFONT, reinterpret_cast<WPARAM>(title_font), TRUE);
        if (page == 0) {
            auto text = wide(display_name(product, options.machine));
            text +=
                removing
                    ? L"\r\n\r\nRemove this installation and its shortcuts. User-created and modified files "
                      L"will be preserved."
                    : L"\r\n\r\n" +
                          wide(product.value("welcome", "This wizard will guide you through installation."));
            if (product.contains("website"))
                text += L"\r\n\r\n" + wide(product.at("website"));
            if (removing)
                text += L"\r\n\r\nLocation: " + options.directory.wstring();
            label(text, 90, 300);
        } else if (page == 1) {
            auto text = wide(product.at("licenseText"));
            std::wstring normalized;
            for (auto c : text) {
                if (c == L'\n' && (normalized.empty() || normalized.back() != L'\r'))
                    normalized += L'\r';
                normalized += c;
            }
            add(L"EDIT", normalized, id_content, 28, 80, 640, 335,
                ES_MULTILINE | ES_READONLY | WS_VSCROLL | WS_TABSTOP);
        } else if (page == 2) {
            label(L"Install for:", 82, 25);
            auto user = add(L"BUTTON", L"Only me (no administrator permission needed)", id_user, 30, 115, 630,
                            28, BS_AUTORADIOBUTTON | WS_GROUP | WS_TABSTOP);
            auto machine = add(L"BUTTON", L"All users (administrator permission required)", id_machine, 30,
                               150, 630, 28, BS_AUTORADIOBUTTON | WS_TABSTOP);
            SendMessageW(options.machine ? machine : user, BM_SETCHECK, BST_CHECKED, 0);
            label(L"Installation folder:", 210, 25);
            add(L"EDIT", options.directory.wstring(), id_directory, 28, 243, 530, 30,
                ES_AUTOHSCROLL | WS_TABSTOP);
            add(L"BUTTON", L"Browse...", id_browse, 570, 243, 100, 30, WS_TABSTOP);
            label(L"Choose a new or empty folder. Different versions can be installed separately. All-users "
                  L"installations use a folder inside Program Files.",
                  295, 100);
        } else if (page == 3) {
            auto a = add(L"BUTTON", L"Create Start Menu shortcuts", id_start, 28, 90, 640, 30,
                         BS_AUTOCHECKBOX | WS_TABSTOP);
            auto b = add(L"BUTTON", L"Create desktop shortcuts", id_desktop, 28, 135, 640, 30,
                         BS_AUTOCHECKBOX | WS_TABSTOP);
            auto c = add(L"BUTTON", L"Add this version's bin directory to PATH", id_path, 28, 180, 640, 30,
                         BS_AUTOCHECKBOX | WS_TABSTOP);
            SendMessageW(a, BM_SETCHECK, options.start_menu ? BST_CHECKED : BST_UNCHECKED, 0);
            SendMessageW(b, BM_SETCHECK, options.desktop ? BST_CHECKED : BST_UNCHECKED, 0);
            SendMessageW(c, BM_SETCHECK, options.path ? BST_CHECKED : BST_UNCHECKED, 0);
            label(L"Shortcuts identify their version and installation scope.\r\n\r\nPATH changes apply to "
                  L"newly started applications. This version is appended; an earlier matching command "
                  L"elsewhere in PATH still takes precedence.",
                  245, 160);
        } else if (page == 4) {
            auto text = wide(display_name(product, options.machine)) + L"\r\n\r\nDestination: " +
                        options.directory.wstring();
            text += L"\r\n\r\nStart Menu: " + std::wstring(options.start_menu ? L"yes" : L"no");
            text += L"\r\nDesktop: " + std::wstring(options.desktop ? L"yes" : L"no");
            text += L"\r\nAdd to PATH: " + std::wstring(options.path ? L"yes" : L"no");
            if (options.machine && !elevated())
                text += L"\r\n\r\nWindows will ask for administrator permission.";
            label(text, 90, 320);
        } else if (page == 5) {
            label(L"Please wait. The operation log records detailed progress.", 85, 55);
            add(L"STATIC", L"Preparing...", id_content, 28, 155, 640, 130);
            add(PROGRESS_CLASSW, L"", id_progress, 28, 300, 640, 25);
            SendMessageW(GetDlgItem(hwnd, id_progress), PBM_SETRANGE32, 0, 1000);
        } else
            add(L"EDIT", finished, id_content, 28, 85, 640, 325,
                ES_MULTILINE | ES_READONLY | WS_VSCROLL | WS_TABSTOP);
        auto back = add(L"BUTTON", L"< Back", id_back, 352, 450, 95, 30, WS_TABSTOP);
        auto next =
            add(L"BUTTON",
                page == 6 ? L"Close" : (page == 4 ? L"Install" : (removing ? L"Uninstall" : L"Next >")),
                id_next, 457, 450, 95, 30, WS_TABSTOP | BS_DEFPUSHBUTTON);
        auto stop = add(L"BUTTON", L"Cancel", id_cancel, 570, 450, 100, 30, WS_TABSTOP);
        EnableWindow(back, page > 0 && page < 5 && !removing);
        EnableWindow(next, page != 5);
        EnableWindow(stop, page != 6 && !(page == 5 && removing));
        SendMessageW(hwnd, DM_SETDEFID, id_next, 0);
        if (page != 5)
            SetFocus(next);
    }
    void capture() {
        if (page == 2)
            options.directory = fs::path(control_text(GetDlgItem(hwnd, id_directory)));
        if (page == 3) {
            options.start_menu = SendMessageW(GetDlgItem(hwnd, id_start), BM_GETCHECK, 0, 0) == BST_CHECKED;
            options.desktop = SendMessageW(GetDlgItem(hwnd, id_desktop), BM_GETCHECK, 0, 0) == BST_CHECKED;
            options.path = SendMessageW(GetDlgItem(hwnd, id_path), BM_GETCHECK, 0, 0) == BST_CHECKED;
        }
    }
    void browse() {
        BROWSEINFOW bi{};
        bi.hwndOwner = hwnd;
        bi.lpszTitle = L"Choose an empty installation folder";
        bi.ulFlags = BIF_RETURNONLYFSDIRS | BIF_NEWDIALOGSTYLE;
        auto item = SHBrowseForFolderW(&bi);
        if (item) {
            std::array<wchar_t, MAX_PATH> path{};
            if (SHGetPathFromIDListW(item, path.data()))
                SetWindowTextW(GetDlgItem(hwnd, id_directory), path.data());
            CoTaskMemFree(item);
        }
    }
    void command(int id) {
        if (id == id_cancel) {
            cancel = true;
            if (!running)
                closed = true;
        } else if (id == id_back) {
            capture();
            --page;
            show_page();
        } else if (id == id_next) {
            if (page == 6) {
                closed = true;
                return;
            }
            capture();
            if (page == 2) {
                validate_destination(options.directory, options.machine);
                if (!fs::is_directory(options.directory.parent_path()))
                    throw Error("Choose a folder whose parent already exists.");
                if (fs::exists(options.directory) && !fs::is_empty(options.directory))
                    throw Error(
                        "Choose a new or empty folder. Uninstall an existing copy before reinstalling.");
                check_registered(product, options.machine);
            }
            if (page == 4 || removing) {
                ready = true;
                running = true;
                page = 5;
                show_page();
            } else {
                ++page;
                show_page();
            }
        } else if (id == id_user || id == id_machine) {
            capture();
            auto old = default_directory(product, options.machine);
            options.machine = id == id_machine;
            if (options.directory == old)
                options.directory = default_directory(product, options.machine);
            SetWindowTextW(GetDlgItem(hwnd, id_directory), options.directory.c_str());
        } else if (id == id_browse)
            browse();
    }
    void pump() {
        MSG msg;
        while (PeekMessageW(&msg, nullptr, 0, 0, PM_REMOVE)) {
            if (msg.message == WM_QUIT) {
                cancel = true;
                closed = true;
                continue;
            }
            if (!IsDialogMessageW(hwnd, &msg)) {
                TranslateMessage(&msg);
                DispatchMessageW(&msg);
            }
        }
    }
    static LRESULT CALLBACK procedure(HWND window, UINT message, WPARAM wp, LPARAM lp) {
        auto self = reinterpret_cast<Impl *>(GetWindowLongPtrW(window, GWLP_USERDATA));
        if (message == WM_NCCREATE) {
            self = static_cast<Impl *>(reinterpret_cast<CREATESTRUCTW *>(lp)->lpCreateParams);
            SetWindowLongPtrW(window, GWLP_USERDATA, reinterpret_cast<LONG_PTR>(self));
        }
        if (self)
            try {
                if (message == WM_COMMAND && HIWORD(wp) == BN_CLICKED) {
                    self->command(LOWORD(wp));
                    return 0;
                }
                if (message == WM_CLOSE) {
                    if (!self->running || !self->removing) {
                        self->cancel = true;
                        if (!self->running)
                            self->closed = true;
                    }
                    return 0;
                }
            } catch (const std::exception &e) {
                MessageBoxW(window, wide(e.what()).c_str(), L"brlinst", MB_OK | MB_ICONERROR);
                return 0;
            }
        return DefWindowProcW(window, message, wp, lp);
    }
};
Wizard::Wizard(const Json &p, Options &o, bool r) : impl(std::make_unique<Impl>(p, o, r)) {}
Wizard::~Wizard() = default;
bool Wizard::configure() {
    while (!impl->ready && !impl->closed) {
        impl->pump();
        if (!impl->ready && !impl->closed)
            WaitMessage();
    }
    return impl->ready;
}
void Wizard::progress(const std::string &text, std::uint64_t done, std::uint64_t total) {
    if (impl->page != 5) {
        impl->page = 5;
        impl->running = true;
        impl->show_page();
    }
    SetWindowTextW(GetDlgItem(impl->hwnd, id_content), wide(text).c_str());
    SendMessageW(GetDlgItem(impl->hwnd, id_progress), PBM_SETPOS,
                 total ? static_cast<WPARAM>(done * 1000 / total) : 0, 0);
    impl->pump();
}
bool Wizard::cancelled() {
    impl->pump();
    return impl->cancel;
}
void Wizard::finish(const std::string &text) {
    impl->running = false;
    impl->page = 6;
    impl->finished = wide(text);
    impl->show_page();
    while (!impl->closed) {
        impl->pump();
        if (!impl->closed)
            WaitMessage();
    }
}
HWND Wizard::window() const {
    return impl->hwnd;
}
} // namespace installer
