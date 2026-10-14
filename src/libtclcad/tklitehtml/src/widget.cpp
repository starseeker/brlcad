/*
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above
 * copyright notice, this list of conditions and the following
 * disclaimer in the documentation and/or other materials provided
 * with the distribution.
 *
 * 3. The name of the author may not be used to endorse or promote
 * products derived from this software without specific prior written
 * permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 * GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifdef BRLCADBUILD
#include "common.h"
#endif

#include "widget.h"
#include "litehtml_compat.h"

#include <litehtml.h>
#include <litehtml/url.h>
#include <litehtml/render_item.h>
#include <tk.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <functional>
#include <initializer_list>
#include <limits>
#include <map>
#include <memory>
#include <new>
#include <numeric>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace tklitehtml
{
namespace
{
constexpr int    kDefaultWidth       = 640;
constexpr int    kDefaultHeight      = 480;
constexpr int    kDefaultFontSize    = 16;
constexpr int    kScrollUnitPixels   = 40;
constexpr float  kCssPixelsPerInch   = 96.0F;
constexpr float  kPointsPerInch      = 72.0F;
constexpr double kCssPixelsPerPoint  = kCssPixelsPerInch / kPointsPerInch;
constexpr int    kBoldWeightMinimum  = 600;
constexpr int    kMediaColorBits     = 8;
constexpr int    kMediaPaletteSize   = 256;
constexpr double kOutlineInnerRadiusSquared = 0.55;
constexpr char   kFileScheme[]       = "file://";
constexpr char   kDefaultBackground[] = "white";
constexpr char   kFindBackground[]    = "#fff29a";
constexpr char   kCurrentBackground[] = "#ffbf69";
constexpr char   kSelectionBackground[] = "#80c0ff";
constexpr int    kDragThresholdPixels = 4;
constexpr char   kFrameImageType[]    = "tklitehtmlframe";
constexpr std::size_t kDefaultMaxDocumentBytes = 16U * 1024U * 1024U;
constexpr std::size_t kDefaultMaxImagePixels    = 16U * 1024U * 1024U;
constexpr std::size_t kDefaultMaxViewportPixels = 16U * 1024U * 1024U;
constexpr std::array<const char*, 10> kOptionNames = {
    "-width",           "-height",         "-background",      "-maxdocumentbytes",
    "-maximagepixels",  "-maxviewportpixels", "-xscrollcommand", "-yscrollcommand",
    "-linkcommand",     "-resourcecommand"};

enum PixelChannel
{
    Red,
    Green,
    Blue,
    Alpha,
    PixelChannelCount
};

constexpr unsigned int kColorComponentMaximum = 255U;

int to_int(litehtml::pixel_t value)
{
    const double number = std::round(static_cast<double>(value));
    if(!std::isfinite(number))
    {
        return 0;
    }
    return static_cast<int>(std::clamp(number, static_cast<double>(std::numeric_limits<int>::min()),
                                       static_cast<double>(std::numeric_limits<int>::max())));
}

int exception_result(Tcl_Interp* interp, const char* operation)
{
    try
    {
        throw;
    } catch(const std::bad_alloc&)
    {
        Tcl_SetObjResult(interp, Tcl_ObjPrintf("%s failed: memory limit exceeded", operation));
    } catch(const std::exception& error)
    {
        Tcl_SetObjResult(interp, Tcl_ObjPrintf("%s failed: %s", operation, error.what()));
    } catch(...)
    {
        Tcl_SetObjResult(interp, Tcl_ObjPrintf("%s failed: unknown error", operation));
    }
    return TCL_ERROR;
}

bool dimensions_within_limit(int width, int height, std::size_t pixel_limit)
{
    if(width < 0 || height < 0)
    {
        return false;
    }
    if(width == 0 || height == 0)
    {
        return true;
    }
    return static_cast<std::size_t>(width) <= pixel_limit / static_cast<std::size_t>(height);
}

int saturating_add(int left, int right)
{
    const auto sum = static_cast<long long>(left) + right;
    return static_cast<int>(std::clamp(sum, static_cast<long long>(std::numeric_limits<int>::min()),
                                       static_cast<long long>(std::numeric_limits<int>::max())));
}

std::filesystem::path path_from_utf8(std::string_view path)
{
    return std::filesystem::u8path(path.begin(), path.end());
}

void replace_obj(Tcl_Obj*& destination, Tcl_Obj* value)
{
    if(value != nullptr)
    {
        Tcl_IncrRefCount(value);
    }
    if(destination != nullptr)
    {
        Tcl_DecrRefCount(destination);
    }
    destination = value;
}

int eval_command(Tcl_Interp* interp, std::initializer_list<Tcl_Obj*> arguments)
{
    for(auto* argument : arguments)
    {
        Tcl_IncrRefCount(argument);
    }
    const int result = Tcl_EvalObjv(interp, static_cast<int>(arguments.size()), arguments.begin(), TCL_EVAL_GLOBAL);
    for(auto* argument : arguments)
    {
        Tcl_DecrRefCount(argument);
    }
    return result;
}

std::string trim_font_family(std::string family)
{
    const auto comma = family.find(',');
    if(comma != std::string::npos)
    {
        family.resize(comma);
    }
    const auto first = family.find_first_not_of(" \t\r\n\"'");
    const auto last  = family.find_last_not_of(" \t\r\n\"'");
    if(first == std::string::npos)
    {
        return "Helvetica";
    }
    family = family.substr(first, last - first + 1);
    if(family == "serif")
    {
        return "Times";
    }
    if(family == "sans-serif" || family == "system-ui")
    {
        return "Helvetica";
    }
    if(family == "monospace")
    {
        return "Courier";
    }
    return family;
}

std::vector<Tcl_UniChar> unicode_chars(std::string_view text)
{
    std::vector<Tcl_UniChar> result;
    const char*              cursor = text.data();
    const char*              end    = cursor + text.size();
    while(cursor < end)
    {
        Tcl_UniChar character = 0;
        const int   length    = Tcl_UtfToUniChar(cursor, &character);
        if(length <= 0)
        {
            break;
        }
        result.push_back(character);
        cursor += length;
    }
    return result;
}

const char* utf_at(std::string_view text, int character_index)
{
    const char* cursor = text.data();
    const char* end    = cursor + text.size();
    for(int index = 0; index < character_index && cursor < end; ++index)
    {
        cursor = Tcl_UtfNext(cursor);
    }
    return cursor;
}

struct TextRange
{
    int first = 0;
    int last  = 0;
};

struct Font
{
    Tk_Font                    tk_font = nullptr;
    Tk_FontMetrics             metrics{};
    litehtml::font_description description;
};

struct TextRun
{
    TextRange             range;
    std::string           text;
    litehtml::position    position;
    Font*                 font = nullptr;
    litehtml::web_color   color;
};

std::pair<int, int> text_run_location(const TextRun& run)
{
    return {to_int(run.position.x), to_int(run.position.y)};
}

struct Highlight
{
    litehtml::web_color     background;
    std::vector<TextRange> ranges;
};

struct PaintText
{
    std::string         text;
    Font*               font = nullptr;
    litehtml::web_color color;
    int                 x = 0;
    int                 y = 0;
};

struct TextColor
{
    XColor* color = nullptr;
    GC      gc    = None;
};

struct Image
{
    std::string name;
    Tk_Image    handle   = nullptr;
    int         width    = 0;
    int         height   = 0;
    std::size_t pixels   = 0;
    bool        internal = false;
};

class HtmlWidget final : public litehtml::document_container
{
  public:
    static int create(Tcl_Interp* interp, int objc, Tcl_Obj* const objv[]);
    static void register_frame_image();

    HtmlWidget(Tcl_Interp* interp, Tk_Window window, std::string canvas_command) :
        interp_(interp),
        window_(window),
        display_(Tk_Display(window)),
        canvas_command_(std::move(canvas_command))
    {
    }

    ~HtmlWidget() override
    {
        if(redraw_pending_)
        {
            Tcl_CancelIdleCall(redraw_idle, this);
        }
        destroy_frame_image();
        document_.reset();
        frame_document_.reset();
        clear_images();
        clear_highlights();
        clear_text_colors();
        free_cursor();
        destroy_backing_image();
        if(background_ != nullptr)
        {
            Tk_FreeColor(background_);
        }
        replace_obj(xscroll_command_, nullptr);
        replace_obj(yscroll_command_, nullptr);
        replace_obj(link_command_, nullptr);
        replace_obj(resource_command_, nullptr);
    }

    HtmlWidget(const HtmlWidget&)            = delete;
    HtmlWidget& operator=(const HtmlWidget&) = delete;

  private:
    static int  widget_command(ClientData client_data, Tcl_Interp* interp, int objc, Tcl_Obj* const objv[]);
    static void command_deleted(ClientData client_data);
    static void free_widget(char* memory);
    static void event_handler(ClientData client_data, XEvent* event);
    static void redraw_idle(ClientData client_data);
    static void image_changed(ClientData client_data, int, int, int, int, int, int);
    static void backing_image_changed(ClientData, int, int, int, int, int, int);
    static int frame_created(Tcl_Interp* interp, const char*, int objc, Tcl_Obj* const objv[],
                             const Tk_ImageType*, Tk_ImageMaster master, ClientData* client_data);
    static ClientData frame_acquired(Tk_Window, ClientData client_data);
    static void frame_displayed(ClientData client_data, Display*, Drawable drawable,
                                 int image_x, int image_y, int width, int height, int drawable_x, int drawable_y);
    static void frame_released(ClientData client_data, Display*);
    static void frame_deleted(ClientData client_data);

    int command(int objc, Tcl_Obj* const objv[]);
    int configure(int objc, Tcl_Obj* const objv[]);
    int cget(Tcl_Obj* option);
    int load_document(int objc, Tcl_Obj* const objv[]);
    int search_command(int objc, Tcl_Obj* const objv[]);
    int highlight_command(int objc, Tcl_Obj* const objv[]);
    int find_command(int objc, Tcl_Obj* const objv[]);
    int fragment_command(Tcl_Obj* fragment_obj);
    int view_command(bool horizontal, int objc, Tcl_Obj* const objv[]);
    int bbox_command(Tcl_Obj* index_obj);
    int see_command(Tcl_Obj* index_obj);
    int selection_command(int objc, Tcl_Obj* const objv[]);
    int set_selection(int first, int last);
    int notify(const char* event);
    int wrong_args(const char* usage) const;

    int set_option(Tcl_Obj* option, Tcl_Obj* value);
    Tcl_Obj* option_value(std::string_view option) const;
    Tcl_Obj* configuration_info() const;

    void handle_event(XEvent* event);
    void schedule_redraw(bool layout = false);
    void redraw();
    void paint_frame(Drawable drawable, int image_x, int image_y, int width, int height,
                      int drawable_x, int drawable_y);
    int  ensure_backing_image();
    int  ensure_frame_image();
    void destroy_frame_image();
    void destroy_backing_image();
    int  resize_canvas();
    int  upload_canvas();
    void layout_document();
    int  ensure_text_index();
    void build_text_index();
    void index_element(const litehtml::element::ptr& element, char& separator);
    void invalidate_text_index();
    bool has_active_highlights() const;
    void update_scroll_commands();
    void invoke_scroll_command(Tcl_Obj* command, int offset, int total, int viewport);
    void clamp_scroll();
    void set_scroll(bool horizontal, int value);

    std::vector<TextRange> search(std::string_view pattern, bool no_case) const;
    int                    ranges_result(const std::vector<TextRange>& ranges, bool all) const;
    int                    parse_ranges(Tcl_Obj* object, std::vector<TextRange>& ranges) const;
    TextRun*               run_at(int index);
    const TextRun*         run_at(int index) const;
    int                    index_from_obj(Tcl_Obj* object, int& index) const;
    int                    index_at(int x, int y) const;

    Highlight& ensure_highlight(const std::string& name);
    int        set_highlight_color(Highlight& highlight, const char* color_name);
    void       clear_highlights();
    void       paint_highlights(const TextRun& run, int screen_x, int screen_y);

    TextColor* text_color_for(const litehtml::web_color& color);
    void       clear_text_colors();
    void fill_rectangle(const litehtml::web_color& color, int x, int y, int width, int height);
    void fill_ellipse(const litehtml::web_color& color, int x, int y, int width, int height, bool outline);
    void blend_pixel(int x, int y, const litehtml::web_color& color);
    litehtml::position active_clip() const;
    void reset_clip();

    std::string resolve_url(const char* source, const char* base) const;
    std::string internal_image_name(std::string_view purpose, unsigned int serial) const;
    std::string resource(const char* kind, const std::string& url);
    Image*      image_for(const std::string& url);
    void        clear_images();
    void        clear_images(std::map<std::string, Image>& images);
    void        delete_internal_image(const Image& image);

    void invoke_link(const std::string& url, const litehtml::element::ptr& element);
    void show_fragment(std::string_view fragment);
    void free_cursor();

    litehtml::uint_ptr create_font(const litehtml::font_description& description, const litehtml::document* document,
                                   litehtml::font_metrics* metrics) override;
    void              delete_font(litehtml::uint_ptr font) override;
    litehtml::pixel_t text_width(const char* text, litehtml::uint_ptr font) override;
    void draw_text(litehtml::uint_ptr drawable, const char* text, litehtml::uint_ptr font,
                   litehtml::web_color color, const litehtml::position& position) override;
    litehtml::pixel_t pt_to_px(float points) const override;
    litehtml::pixel_t get_default_font_size() const override;
    const char*       get_default_font_name() const override;
    void draw_list_marker(litehtml::uint_ptr drawable, const litehtml::list_marker& marker) override;
    void load_image(const char* source, const char* base, bool redraw_on_ready) override;
    void get_image_size(const char* source, const char* base, litehtml::size& size) override;
    void draw_image(litehtml::uint_ptr drawable, const litehtml::background_layer& layer, const std::string& url,
                    const std::string& base) override;
    void draw_solid_fill(litehtml::uint_ptr drawable, const litehtml::background_layer& layer,
                         const litehtml::web_color& color) override;
    void draw_linear_gradient(litehtml::uint_ptr drawable, const litehtml::background_layer& layer,
                              const litehtml::background_layer::linear_gradient& gradient) override;
    void draw_radial_gradient(litehtml::uint_ptr drawable, const litehtml::background_layer& layer,
                              const litehtml::background_layer::radial_gradient& gradient) override;
    void draw_conic_gradient(litehtml::uint_ptr drawable, const litehtml::background_layer& layer,
                             const litehtml::background_layer::conic_gradient& gradient) override;
    void draw_borders(litehtml::uint_ptr drawable, const litehtml::borders& borders,
                      const litehtml::position& position, bool root) override;
    void set_caption(const char* caption) override;
    void set_base_url(const char* url) override;
    void link(const litehtml::document::ptr& document, const litehtml::element::ptr& element) override;
    void on_anchor_click(const char* url, const litehtml::element::ptr& element) override;
    void on_mouse_event(const litehtml::element::ptr& element, litehtml::mouse_event event) override;
    void set_cursor(const char* cursor) override;
    void transform_text(std::string& text, litehtml::text_transform transform) override;
    void import_css(std::string& text, const std::string& url, std::string& base) override;
    void set_clip(const litehtml::position& position, const litehtml::border_radiuses& radii) override;
    void del_clip() override;
    void get_viewport(litehtml::position& viewport) const override;
    litehtml::element::ptr create_element(const char* tag_name, const litehtml::string_map& attributes,
                                          const litehtml::document::ptr& document) override;
    void get_media_features(litehtml::media_features& media) const override;
    void get_language(std::string& language, std::string& culture) const override;

    Tcl_Interp* interp_  = nullptr;
    Tk_Window   window_  = nullptr;
    Display*    display_ = nullptr;
    Tcl_Command command_token_ = nullptr;

    int         requested_width_        = kDefaultWidth;
    int         requested_height_       = kDefaultHeight;
    int         scroll_x_               = 0;
    int         scroll_y_               = 0;
    int         document_width_         = 0;
    int         document_height_        = 0;
    std::size_t max_document_bytes_     = kDefaultMaxDocumentBytes;
    std::size_t max_image_pixels_       = kDefaultMaxImagePixels;
    std::size_t max_viewport_pixels_    = kDefaultMaxViewportPixels;
    std::size_t document_resource_bytes_ = 0;
    std::size_t image_pixels_            = 0;

    Tcl_Obj* xscroll_command_ = nullptr;
    Tcl_Obj* yscroll_command_ = nullptr;
    Tcl_Obj* link_command_    = nullptr;
    Tcl_Obj* resource_command_ = nullptr;
    XColor*  background_       = nullptr;

    std::string base_url_;
    std::string canvas_command_;
    std::string document_base_url_;
    std::string html_source_;
    std::string user_stylesheet_;
    std::string plain_text_;
    std::string title_;

    litehtml::document::ptr  document_;
    litehtml::document::ptr  frame_document_;
    std::vector<TextRun>     text_runs_;
    std::vector<std::size_t> text_run_lookup_;
    std::map<std::string, Highlight> highlights_;
    std::vector<TextRange>  find_ranges_;
    int                     current_find_         = -1;
    int                     text_character_count_ = 0;
    TextRange               selection_;
    int                     button_x_ = 0;
    int                     button_y_ = 0;
    bool                    button_down_ = false;
    bool                    suppress_link_ = false;

    std::map<std::uint32_t, TextColor> text_colors_;
    std::map<std::string, Image>     images_;
    unsigned int                    image_serial_ = 0;

    std::vector<unsigned char> pixels_;
    std::vector<PaintText>     paint_text_;
    int                        canvas_width_  = 0;
    int                        canvas_height_ = 0;
    std::string                backing_image_name_;
    Tk_PhotoHandle             backing_photo_ = nullptr;
    Tk_Image                   backing_image_ = nullptr;
    std::string                frame_image_name_;
    Tk_ImageMaster             frame_master_ = nullptr;

    std::vector<litehtml::position> clip_stack_;
    Tk_Cursor                      cursor_ = None;

    bool redraw_pending_ = false;
    bool layout_pending_ = false;
    bool indexing_       = false;
    bool text_index_valid_ = false;
    bool loading_        = false;
    bool destroying_     = false;
};

int HtmlWidget::create(Tcl_Interp* interp, int objc, Tcl_Obj* const objv[])
{
    if(objc < 2)
    {
        Tcl_WrongNumArgs(interp, 1, objv, "pathName ?option value ...?");
        return TCL_ERROR;
    }

    Tk_Window main_window = Tk_MainWindow(interp);
    if(main_window == nullptr)
    {
        Tcl_SetObjResult(interp, Tcl_NewStringObj("Tk is not initialized", -1));
        return TCL_ERROR;
    }

    const char* path = Tcl_GetString(objv[1]);
    std::string canvas_command = std::string(path) + ".__canvas";
    // A canvas owns the platform buffering; our frame image paints its
    // background and text into the same drawable through public Tk APIs.
    if(eval_command(interp, {Tcl_NewStringObj("canvas", -1), objv[1],
                             Tcl_NewStringObj("-borderwidth", -1), Tcl_NewIntObj(0),
                             Tcl_NewStringObj("-highlightthickness", -1), Tcl_NewIntObj(0)}) != TCL_OK)
    {
        return TCL_ERROR;
    }
    Tk_Window window = Tk_NameToWindow(interp, path, main_window);
    if(window == nullptr)
    {
        return TCL_ERROR;
    }
    if(eval_command(interp, {Tcl_NewStringObj("rename", -1), objv[1],
                             Tcl_NewStringObj(canvas_command.c_str(), -1)}) != TCL_OK)
    {
        Tk_DestroyWindow(window);
        return TCL_ERROR;
    }

    try
    {
        auto widget = std::make_unique<HtmlWidget>(interp, window, std::move(canvas_command));
        Tk_SetClass(window, "TkLiteHtml");
        widget->background_ = Tk_GetColor(interp, window, Tk_GetUid(kDefaultBackground));
        if(widget->background_ == nullptr)
        {
            Tk_DestroyWindow(window);
            return TCL_ERROR;
        }
        Tk_SetWindowBackground(window, widget->background_->pixel);

        if(widget->configure(objc - 2, objv + 2) != TCL_OK)
        {
            Tk_DestroyWindow(window);
            return TCL_ERROR;
        }

        Tk_CreateEventHandler(
            window,
            ExposureMask | StructureNotifyMask | PointerMotionMask | LeaveWindowMask | ButtonPressMask |
                ButtonReleaseMask,
            event_handler, widget.get());
        widget->command_token_ =
            Tcl_CreateObjCommand(interp, Tk_PathName(window), widget_command, widget.get(), command_deleted);
        Tk_GeometryRequest(window, widget->requested_width_, widget->requested_height_);
        Tcl_SetObjResult(interp, Tcl_NewStringObj(Tk_PathName(window), -1));
        widget.release();
        return TCL_OK;
    } catch(...)
    {
        Tk_DestroyWindow(window);
        return exception_result(interp, "widget creation");
    }
}

int HtmlWidget::widget_command(ClientData client_data, Tcl_Interp* interp, int objc, Tcl_Obj* const objv[])
{
    auto* widget = static_cast<HtmlWidget*>(client_data);
    Tcl_Preserve(widget);
    int result = TCL_ERROR;
    try
    {
        result = widget->command(objc, objv);
    } catch(...)
    {
        result = exception_result(widget->interp_, "widget command");
    }
    Tcl_Obj* error_result = nullptr;
    if(result == TCL_ERROR)
    {
        if(Tcl_GetCharLength(Tcl_GetObjResult(interp)) == 0)
        {
            Tcl_SetObjResult(interp, Tcl_NewStringObj("widget command failed after the widget was destroyed", -1));
        }
        error_result = Tcl_GetObjResult(interp);
        Tcl_IncrRefCount(error_result);
    }
    Tcl_Release(widget);
    if(error_result != nullptr)
    {
        Tcl_SetObjResult(interp, error_result);
        Tcl_DecrRefCount(error_result);
    }
    return result;
}

void HtmlWidget::command_deleted(ClientData client_data)
{
    auto* widget = static_cast<HtmlWidget*>(client_data);
    widget->command_token_ = nullptr;
    // Detach image instances before widget deletion can release their
    // last reference from inside Tk's image cleanup callback.
    widget->destroy_frame_image();
    if(widget->window_ != nullptr)
    {
        widget->destroying_ = true;
        Tk_DestroyWindow(widget->window_);
    }
    Tcl_EventuallyFree(widget, free_widget);
}

void HtmlWidget::free_widget(char* memory)
{
    delete reinterpret_cast<HtmlWidget*>(memory);
}

void HtmlWidget::event_handler(ClientData client_data, XEvent* event)
{
    auto* widget = static_cast<HtmlWidget*>(client_data);
    Tcl_Preserve(widget);
    try
    {
        widget->handle_event(event);
    } catch(...)
    {
        exception_result(widget->interp_, "widget event");
        Tcl_BackgroundException(widget->interp_, TCL_ERROR);
    }
    Tcl_Release(widget);
}

void HtmlWidget::redraw_idle(ClientData client_data)
{
    auto* widget = static_cast<HtmlWidget*>(client_data);
    Tcl_Preserve(widget);
    try
    {
        widget->redraw_pending_ = false;
        widget->redraw();
    } catch(...)
    {
        exception_result(widget->interp_, "widget redraw");
        Tcl_BackgroundException(widget->interp_, TCL_ERROR);
    }
    Tcl_Release(widget);
}

void HtmlWidget::image_changed(ClientData client_data, int, int, int, int, int, int)
{
    auto* widget = static_cast<HtmlWidget*>(client_data);
    Tcl_Preserve(widget);
    try
    {
        for(auto& [url, image] : widget->images_)
        {
            const std::size_t other_pixels = widget->image_pixels_ - image.pixels;
            Tk_SizeOfImage(image.handle, &image.width, &image.height);
            if(!dimensions_within_limit(image.width, image.height,
                                        widget->max_image_pixels_ - other_pixels))
            {
                image.width = 0;
                image.height = 0;
                image.pixels = 0;
                widget->image_pixels_ = other_pixels;
                const std::string message = "image '" + url + "' exceeds the " +
                    std::to_string(widget->max_image_pixels_) + " pixel limit";
                Tcl_SetObjResult(widget->interp_, Tcl_NewStringObj(message.c_str(), -1));
                Tcl_BackgroundException(widget->interp_, TCL_ERROR);
            } else
            {
                image.pixels = static_cast<std::size_t>(image.width) * image.height;
                widget->image_pixels_ = other_pixels + image.pixels;
            }
        }
        widget->schedule_redraw(true);
    } catch(...)
    {
        exception_result(widget->interp_, "image update");
        Tcl_BackgroundException(widget->interp_, TCL_ERROR);
    }
    Tcl_Release(widget);
}

void HtmlWidget::backing_image_changed(ClientData, int, int, int, int, int, int) {}

void HtmlWidget::register_frame_image()
{
    // Tk retains image type records for the lifetime of the process.
    static Tk_ImageType type = [] {
        Tk_ImageType image_type{};
        image_type.name = kFrameImageType;
        image_type.createProc = frame_created;
        image_type.getProc = frame_acquired;
        image_type.displayProc = frame_displayed;
        image_type.freeProc = frame_released;
        image_type.deleteProc = frame_deleted;
        return image_type;
    }();
    Tk_CreateImageType(&type);
}

int HtmlWidget::frame_created(Tcl_Interp* interp, const char*, int objc, Tcl_Obj* const objv[],
                               const Tk_ImageType*, Tk_ImageMaster master, ClientData* client_data)
{
    Tcl_CmdInfo command{};
    if(objc != 1 || !Tcl_GetCommandInfo(interp, Tcl_GetString(objv[0]), &command) ||
       command.objProc != widget_command)
    {
        Tcl_SetObjResult(interp, Tcl_NewStringObj("a frame image requires an HTML widget", -1));
        return TCL_ERROR;
    }
    auto* widget = static_cast<HtmlWidget*>(command.objClientData);
    widget->frame_master_ = master;
    *client_data = widget;
    return TCL_OK;
}

ClientData HtmlWidget::frame_acquired(Tk_Window, ClientData client_data)
{
    // Canvas image instances must release the widget before it is freed.
    Tcl_Preserve(client_data);
    return client_data;
}

void HtmlWidget::frame_displayed(ClientData client_data, Display*, Drawable drawable,
                                  int image_x, int image_y, int width, int height, int drawable_x, int drawable_y)
{
    auto* widget = static_cast<HtmlWidget*>(client_data);
    Tcl_Preserve(widget);
    try
    {
        widget->paint_frame(drawable, image_x, image_y, width, height, drawable_x, drawable_y);
    } catch(...)
    {
        exception_result(widget->interp_, "frame painting");
        Tcl_BackgroundException(widget->interp_, TCL_ERROR);
    }
    Tcl_Release(widget);
}

void HtmlWidget::frame_released(ClientData client_data, Display*)
{
    Tcl_Release(client_data);
}

void HtmlWidget::frame_deleted(ClientData client_data)
{
    static_cast<HtmlWidget*>(client_data)->frame_master_ = nullptr;
}

int HtmlWidget::command(int objc, Tcl_Obj* const objv[])
{
    static const char* const commands[] = {"bbox",      "cget", "clear", "configure", "find", "fragment",
                                            "highlight", "index", "load", "search", "see", "selection", "text", "title",
                                            "xview",     "yview", nullptr};
    enum Command
    {
        Bbox,
        Cget,
        Clear,
        Configure,
        Find,
        Fragment,
        HighlightCommand,
        Index,
        Load,
        Search,
        See,
        SelectionCommand,
        Text,
        Title,
        Xview,
        Yview
    };

    if(objc < 2)
    {
        Tcl_WrongNumArgs(interp_, 1, objv, "subcommand ?arg ...?");
        return TCL_ERROR;
    }
    int index = 0;
    if(Tcl_GetIndexFromObj(interp_, objv[1], commands, "subcommand", 0, &index) != TCL_OK)
    {
        return TCL_ERROR;
    }
    if(loading_ && static_cast<Command>(index) != Cget)
    {
        Tcl_SetObjResult(interp_, Tcl_NewStringObj("widget command is unavailable during document loading", -1));
        return TCL_ERROR;
    }

    switch(static_cast<Command>(index))
    {
    case Bbox:
        if(objc != 3)
        {
            Tcl_WrongNumArgs(interp_, 2, objv, "index");
            return TCL_ERROR;
        }
        return bbox_command(objv[2]);
    case Cget:
        if(objc != 3)
        {
            Tcl_WrongNumArgs(interp_, 2, objv, "option");
            return TCL_ERROR;
        }
        return cget(objv[2]);
    case Clear:
        if(objc != 2)
        {
            Tcl_WrongNumArgs(interp_, 2, objv, nullptr);
            return TCL_ERROR;
        }
        document_.reset();
        html_source_.clear();
        title_.clear();
        invalidate_text_index();
        find_ranges_.clear();
        current_find_ = -1;
        document_width_ = document_height_ = scroll_x_ = scroll_y_ = 0;
        document_resource_bytes_ = 0;
        clear_images();
        for(auto& [name, highlight] : highlights_)
        {
            highlight.ranges.clear();
        }
        if(set_selection(0, 0) != TCL_OK)
        {
            return TCL_ERROR;
        }
        schedule_redraw();
        update_scroll_commands();
        return notify("<<DocumentChanged>>");
    case Configure:
        if(objc == 2)
        {
            Tcl_SetObjResult(interp_, configuration_info());
            return TCL_OK;
        }
        if(objc == 3)
        {
            return cget(objv[2]);
        }
        return configure(objc - 2, objv + 2);
    case Find:
        return find_command(objc - 2, objv + 2);
    case Fragment:
        if(objc != 3)
        {
            Tcl_WrongNumArgs(interp_, 2, objv, "identifier");
            return TCL_ERROR;
        }
        return fragment_command(objv[2]);
    case HighlightCommand:
        return highlight_command(objc - 2, objv + 2);
    case Index:
    {
        if(objc != 3)
        {
            return wrong_args("index index");
        }
        if(ensure_text_index() != TCL_OK)
        {
            return TCL_ERROR;
        }
        int text_index = 0;
        if(index_from_obj(objv[2], text_index) != TCL_OK)
        {
            return TCL_ERROR;
        }
        Tcl_SetObjResult(interp_, Tcl_NewIntObj(text_index));
        return TCL_OK;
    }
    case Load:
        return load_document(objc - 2, objv + 2);
    case Search:
        return search_command(objc - 2, objv + 2);
    case See:
        if(objc != 3)
        {
            Tcl_WrongNumArgs(interp_, 2, objv, "index");
            return TCL_ERROR;
        }
        return see_command(objv[2]);
    case SelectionCommand:
        return selection_command(objc - 2, objv + 2);
    case Text:
        if(objc != 2)
        {
            Tcl_WrongNumArgs(interp_, 2, objv, nullptr);
            return TCL_ERROR;
        }
        if(ensure_text_index() != TCL_OK)
        {
            return TCL_ERROR;
        }
        Tcl_SetObjResult(interp_, Tcl_NewStringObj(plain_text_.data(), static_cast<int>(plain_text_.size())));
        return TCL_OK;
    case Title:
        if(objc != 2)
        {
            Tcl_WrongNumArgs(interp_, 2, objv, nullptr);
            return TCL_ERROR;
        }
        Tcl_SetObjResult(interp_, Tcl_NewStringObj(title_.data(), static_cast<int>(title_.size())));
        return TCL_OK;
    case Xview:
        return view_command(true, objc - 2, objv + 2);
    case Yview:
        return view_command(false, objc - 2, objv + 2);
    }
    return TCL_ERROR;
}

int HtmlWidget::configure(int objc, Tcl_Obj* const objv[])
{
    if(objc % 2 != 0)
    {
        Tcl_SetObjResult(interp_, Tcl_NewStringObj("options must be specified as option-value pairs", -1));
        return TCL_ERROR;
    }
    for(int index = 0; index < objc; index += 2)
    {
        if(set_option(objv[index], objv[index + 1]) != TCL_OK)
        {
            return TCL_ERROR;
        }
    }
    if(window_ != nullptr)
    {
        Tk_GeometryRequest(window_, requested_width_, requested_height_);
    }
    schedule_redraw(true);
    return TCL_OK;
}

int HtmlWidget::set_option(Tcl_Obj* option, Tcl_Obj* value)
{
    const std::string name = Tcl_GetString(option);
    if(name == "-width" || name == "-height")
    {
        int pixels = 0;
        if(Tk_GetPixelsFromObj(interp_, window_, value, &pixels) != TCL_OK)
        {
            return TCL_ERROR;
        }
        if(pixels <= 0)
        {
            Tcl_SetObjResult(interp_, Tcl_NewStringObj("width and height must be positive", -1));
            return TCL_ERROR;
        }
        (name == "-width" ? requested_width_ : requested_height_) = pixels;
        return TCL_OK;
    }
    if(name == "-background")
    {
        XColor* color = Tk_GetColor(interp_, window_, Tk_GetUid(Tcl_GetString(value)));
        if(color == nullptr)
        {
            return TCL_ERROR;
        }
        if(background_ != nullptr)
        {
            Tk_FreeColor(background_);
        }
        background_ = color;
        Tk_SetWindowBackground(window_, background_->pixel);
        return TCL_OK;
    }
    if(name == "-maxdocumentbytes" || name == "-maximagepixels" || name == "-maxviewportpixels")
    {
        Tcl_WideInt limit = 0;
        if(Tcl_GetWideIntFromObj(interp_, value, &limit) != TCL_OK)
        {
            return TCL_ERROR;
        }
        const auto maximum = name == "-maxdocumentbytes" ?
            static_cast<unsigned long long>(std::numeric_limits<int>::max()) :
            static_cast<unsigned long long>(std::numeric_limits<std::size_t>::max());
        if(limit <= 0 || static_cast<unsigned long long>(limit) > maximum)
        {
            const std::string message = name + " must be between 1 and " + std::to_string(maximum);
            Tcl_SetObjResult(interp_, Tcl_NewStringObj(message.c_str(), -1));
            return TCL_ERROR;
        }
        const auto requested_limit = static_cast<std::size_t>(limit);
        const auto viewport_pixels = static_cast<std::size_t>(canvas_width_) * canvas_height_;
        if((name == "-maxdocumentbytes" && requested_limit < document_resource_bytes_) ||
           (name == "-maximagepixels" && requested_limit < image_pixels_) ||
           (name == "-maxviewportpixels" && requested_limit < viewport_pixels))
        {
            Tcl_SetObjResult(interp_, Tcl_ObjPrintf("%s is smaller than the current document", name.c_str()));
            return TCL_ERROR;
        }
        auto& destination = name == "-maxdocumentbytes" ? max_document_bytes_ :
            name == "-maximagepixels" ? max_image_pixels_ : max_viewport_pixels_;
        destination = requested_limit;
        return TCL_OK;
    }
    if(name == "-xscrollcommand")
    {
        replace_obj(xscroll_command_, value);
        return TCL_OK;
    }
    if(name == "-yscrollcommand")
    {
        replace_obj(yscroll_command_, value);
        return TCL_OK;
    }
    if(name == "-linkcommand")
    {
        replace_obj(link_command_, value);
        return TCL_OK;
    }
    if(name == "-resourcecommand")
    {
        replace_obj(resource_command_, value);
        return TCL_OK;
    }
    Tcl_SetObjResult(interp_, Tcl_ObjPrintf("unknown option \"%s\"", name.c_str()));
    return TCL_ERROR;
}

Tcl_Obj* HtmlWidget::option_value(std::string_view option) const
{
    if(option == "-width")
    {
        return Tcl_NewIntObj(requested_width_);
    }
    if(option == "-height")
    {
        return Tcl_NewIntObj(requested_height_);
    }
    if(option == "-background")
    {
        return Tcl_NewStringObj(background_ != nullptr ? Tk_NameOfColor(background_) : "", -1);
    }
    if(option == "-maxdocumentbytes")
    {
        return Tcl_NewWideIntObj(static_cast<Tcl_WideInt>(max_document_bytes_));
    }
    if(option == "-maximagepixels")
    {
        return Tcl_NewWideIntObj(static_cast<Tcl_WideInt>(max_image_pixels_));
    }
    if(option == "-maxviewportpixels")
    {
        return Tcl_NewWideIntObj(static_cast<Tcl_WideInt>(max_viewport_pixels_));
    }
    Tcl_Obj* object = nullptr;
    if(option == "-xscrollcommand")
    {
        object = xscroll_command_;
    } else if(option == "-yscrollcommand")
    {
        object = yscroll_command_;
    } else if(option == "-linkcommand")
    {
        object = link_command_;
    } else if(option == "-resourcecommand")
    {
        object = resource_command_;
    }
    return object != nullptr ? Tcl_DuplicateObj(object) : Tcl_NewObj();
}

Tcl_Obj* HtmlWidget::configuration_info() const
{
    Tcl_Obj* result = Tcl_NewListObj(0, nullptr);
    for(const char* name : kOptionNames)
    {
        Tcl_Obj* pair[2] = {Tcl_NewStringObj(name, -1), option_value(name)};
        Tcl_ListObjAppendElement(interp_, result, Tcl_NewListObj(2, pair));
    }
    return result;
}

int HtmlWidget::cget(Tcl_Obj* option)
{
    const std::string name = Tcl_GetString(option);
    const auto known = std::find_if(kOptionNames.begin(), kOptionNames.end(), [&name](const char* option_name) {
        return name == option_name;
    });
    if(known == kOptionNames.end())
    {
        Tcl_SetObjResult(interp_, Tcl_ObjPrintf("unknown option \"%s\"", name.c_str()));
        return TCL_ERROR;
    }
    Tcl_SetObjResult(interp_, option_value(name));
    return TCL_OK;
}

int HtmlWidget::wrong_args(const char* usage) const
{
    const char* path = window_ != nullptr ? Tk_PathName(window_) : "::tklitehtml::widget";
    Tcl_SetObjResult(interp_, Tcl_ObjPrintf("wrong # args: should be \"%s %s\"", path, usage));
    return TCL_ERROR;
}

int HtmlWidget::load_document(int objc, Tcl_Obj* const objv[])
{
    std::string base  = base_url_;
    std::string style;
    int         index = 0;
    while(index < objc - 1)
    {
        const std::string option = Tcl_GetString(objv[index]);
        if(option != "-baseurl" && option != "-userstyle")
        {
            break;
        }
        if(index + 1 >= objc - 1)
        {
            Tcl_SetObjResult(interp_, Tcl_ObjPrintf("missing value for \"%s\"", option.c_str()));
            return TCL_ERROR;
        }
        if(option == "-baseurl")
        {
            base = Tcl_GetString(objv[index + 1]);
        } else
        {
            style = Tcl_GetString(objv[index + 1]);
        }
        index += 2;
    }
    if(index != objc - 1)
    {
        return wrong_args("load ?-baseurl uri? ?-userstyle css? html");
    }

    int         length = 0;
    const char* bytes  = Tcl_GetStringFromObj(objv[index], &length);
    const auto source_size = static_cast<std::size_t>(length);
    if(source_size > max_document_bytes_ || style.size() > max_document_bytes_ - source_size)
    {
        const std::string message = "HTML document exceeds the " + std::to_string(max_document_bytes_) +
            " byte limit";
        Tcl_SetObjResult(interp_, Tcl_NewStringObj(message.c_str(), -1));
        return TCL_ERROR;
    }
    auto previous_document = std::move(document_);
    auto previous_images = std::move(images_);
    auto previous_text_runs = std::move(text_runs_);
    auto previous_text_run_lookup = std::move(text_run_lookup_);
    auto previous_plain_text = std::move(plain_text_);
    auto previous_html_source = std::move(html_source_);
    auto previous_base_url = std::move(base_url_);
    auto previous_document_base_url = std::move(document_base_url_);
    auto previous_user_stylesheet = std::move(user_stylesheet_);
    auto previous_title = std::move(title_);
    const std::size_t previous_document_resource_bytes = document_resource_bytes_;
    const std::size_t previous_image_pixels = image_pixels_;
    const int previous_text_character_count = text_character_count_;
    const int previous_document_width = document_width_;
    const int previous_document_height = document_height_;
    const int previous_scroll_x = scroll_x_;
    const int previous_scroll_y = scroll_y_;
    const bool previous_text_index_valid = text_index_valid_;

    const auto restore_previous_document = [&]() {
        document_.reset();
        clear_images();
        document_ = std::move(previous_document);
        images_ = std::move(previous_images);
        text_runs_ = std::move(previous_text_runs);
        text_run_lookup_ = std::move(previous_text_run_lookup);
        plain_text_ = std::move(previous_plain_text);
        html_source_ = std::move(previous_html_source);
        base_url_ = std::move(previous_base_url);
        document_base_url_ = std::move(previous_document_base_url);
        user_stylesheet_ = std::move(previous_user_stylesheet);
        title_ = std::move(previous_title);
        document_resource_bytes_ = previous_document_resource_bytes;
        image_pixels_ = previous_image_pixels;
        text_character_count_ = previous_text_character_count;
        document_width_ = previous_document_width;
        document_height_ = previous_document_height;
        scroll_x_ = previous_scroll_x;
        scroll_y_ = previous_scroll_y;
        text_index_valid_ = previous_text_index_valid;
    };

    loading_ = true;
    try
    {
        html_source_.assign(bytes, source_size);
        base_url_ = std::move(base);
        document_base_url_ = base_url_;
        user_stylesheet_ = std::move(style);
        document_resource_bytes_ = source_size + user_stylesheet_.size();
        image_pixels_ = 0;
        title_.clear();
        scroll_x_ = 0;
        scroll_y_ = 0;
        invalidate_text_index();
        document_ = litehtml::document::createFromString(html_source_.c_str(), this, litehtml::master_css,
                                                         user_stylesheet_);
        if(document_ == nullptr || window_ == nullptr)
        {
            throw std::runtime_error("HTML parser did not produce a document");
        }
        layout_document();
        if(window_ == nullptr)
        {
            throw std::runtime_error("widget was destroyed while loading the document");
        }
    } catch(...)
    {
        loading_ = false;
        restore_previous_document();
        return exception_result(interp_, "document load");
    }
    loading_ = false;
    previous_document.reset();
    clear_images(previous_images);
    clear_highlights();
    find_ranges_.clear();
    current_find_ = -1;
    if(set_selection(0, 0) != TCL_OK)
    {
        return TCL_ERROR;
    }
    button_down_ = false;
    suppress_link_ = false;
    update_scroll_commands();
    schedule_redraw();
    return notify("<<DocumentChanged>>");
}

int HtmlWidget::fragment_command(Tcl_Obj* fragment_obj)
{
    const char* fragment = Tcl_GetString(fragment_obj);
    if(fragment[0] == '#')
    {
        ++fragment;
    }
    show_fragment(fragment);
    return TCL_OK;
}

void HtmlWidget::handle_event(XEvent* event)
{
    if(event->type == DestroyNotify)
    {
        window_ = nullptr;
        if(command_token_ != nullptr && !destroying_)
        {
            Tcl_DeleteCommandFromToken(interp_, command_token_);
        }
        return;
    }
    if(event->type == ConfigureNotify)
    {
        schedule_redraw(true);
        return;
    }
    if(event->type == Expose)
    {
        schedule_redraw();
        return;
    }
    if(document_ == nullptr)
    {
        return;
    }

#if TKLITEHTML_REDRAW_CALLBACK
    auto redraw = [this](const litehtml::position&) { schedule_redraw(); };
#else
    litehtml::position::vector redraw;
#endif
    switch(event->type)
    {
    case MotionNotify:
        if(button_down_ &&
           (std::abs(static_cast<long long>(event->xmotion.x) - button_x_) +
            std::abs(static_cast<long long>(event->xmotion.y) - button_y_) >= kDragThresholdPixels))
        {
            suppress_link_ = true;
        }
        document_->on_mouse_over(event->xmotion.x + scroll_x_, event->xmotion.y + scroll_y_, event->xmotion.x,
                                 event->xmotion.y, redraw);
        break;
    case LeaveNotify:
        document_->on_mouse_leave(redraw);
        break;
    case ButtonPress:
        if(event->xbutton.button == Button1)
        {
            button_down_ = true;
            suppress_link_ = (event->xbutton.state & (ControlMask | ShiftMask)) != 0;
            button_x_ = event->xbutton.x;
            button_y_ = event->xbutton.y;
            document_->on_lbutton_down(event->xbutton.x + scroll_x_, event->xbutton.y + scroll_y_, event->xbutton.x,
                                       event->xbutton.y, redraw);
        }
        break;
    case ButtonRelease:
        if(event->xbutton.button == Button1)
        {
            button_down_ = false;
            document_->on_lbutton_up(event->xbutton.x + scroll_x_, event->xbutton.y + scroll_y_, event->xbutton.x,
                                     event->xbutton.y, redraw);
        }
        break;
    default:
        break;
    }
#if !TKLITEHTML_REDRAW_CALLBACK
    if(!redraw.empty())
    {
        schedule_redraw();
    }
#endif
}

void HtmlWidget::schedule_redraw(bool layout)
{
    layout_pending_ = layout_pending_ || layout;
    if(!redraw_pending_ && window_ != nullptr)
    {
        redraw_pending_ = true;
        Tcl_DoWhenIdle(redraw_idle, this);
    }
}

void HtmlWidget::redraw()
{
    if(window_ == nullptr || !Tk_IsMapped(window_))
    {
        return;
    }
    if(layout_pending_)
    {
        layout_pending_ = false;
        layout_document();
    }
    if(ensure_backing_image() != TCL_OK)
    {
        Tcl_BackgroundException(interp_, TCL_ERROR);
        return;
    }

    if(resize_canvas() != TCL_OK)
    {
        Tcl_BackgroundException(interp_, TCL_ERROR);
        return;
    }
    if(has_active_highlights() && ensure_text_index() != TCL_OK)
    {
        Tcl_BackgroundException(interp_, TCL_ERROR);
        return;
    }
    reset_clip();
    paint_text_.clear();
    // Tk may repaint the previous frame before the next redraw idle.
    // Keep its fonts alive when a page is replaced or cleared.
    frame_document_ = document_;
    const litehtml::web_color background(static_cast<litehtml::byte>(background_->red >> 8),
                                         static_cast<litehtml::byte>(background_->green >> 8),
                                         static_cast<litehtml::byte>(background_->blue >> 8));
    fill_rectangle(background, 0, 0, canvas_width_, canvas_height_);
    if(document_ != nullptr)
    {
        litehtml::position clip(0, 0, canvas_width_, canvas_height_);
        document_->draw(0, -scroll_x_, -scroll_y_, &clip);
    }
    reset_clip();
    if(upload_canvas() != TCL_OK)
    {
        Tcl_BackgroundException(interp_, TCL_ERROR);
        return;
    }
    if(ensure_frame_image() != TCL_OK)
    {
        Tcl_BackgroundException(interp_, TCL_ERROR);
        return;
    }
    Tk_ImageChanged(frame_master_, 0, 0, canvas_width_, canvas_height_, canvas_width_, canvas_height_);
}

void HtmlWidget::paint_frame(Drawable drawable, int image_x, int image_y, int width, int height,
                              int drawable_x, int drawable_y)
{
    if(backing_image_ == nullptr)
    {
        return;
    }
    Tk_RedrawImage(backing_image_, image_x, image_y, width, height, drawable, drawable_x, drawable_y);
    const int offset_x = drawable_x - image_x;
    const int offset_y = drawable_y - image_y;
    for(const auto& text : paint_text_)
    {
        TextColor* color = text_color_for(text.color);
        if(color == nullptr)
        {
            continue;
        }
        Tk_DrawChars(display_, drawable, color->gc, text.font->tk_font, text.text.c_str(),
                     static_cast<int>(text.text.size()), saturating_add(text.x, offset_x),
                     saturating_add(saturating_add(text.y, offset_y), text.font->metrics.ascent));
    }
}

int HtmlWidget::ensure_backing_image()
{
    if(backing_image_ != nullptr)
    {
        return TCL_OK;
    }
    backing_image_name_ = internal_image_name("backing", 0);
    const int result = eval_command(interp_, {Tcl_NewStringObj("image", -1), Tcl_NewStringObj("create", -1),
                                             Tcl_NewStringObj("photo", -1),
                                             Tcl_NewStringObj(backing_image_name_.c_str(), -1)});
    if(result != TCL_OK)
    {
        return TCL_ERROR;
    }
    backing_photo_ = Tk_FindPhoto(interp_, backing_image_name_.c_str());
    backing_image_ = Tk_GetImage(interp_, window_, backing_image_name_.c_str(), backing_image_changed, this);
    if(backing_photo_ == nullptr || backing_image_ == nullptr)
    {
        Tcl_SetObjResult(interp_, Tcl_NewStringObj("unable to create widget backing image", -1));
        destroy_backing_image();
        return TCL_ERROR;
    }
    Tcl_ResetResult(interp_);
    return TCL_OK;
}

int HtmlWidget::ensure_frame_image()
{
    if(frame_master_ != nullptr)
    {
        return TCL_OK;
    }
    frame_image_name_ = internal_image_name("frame", 0);
    if(eval_command(interp_, {Tcl_NewStringObj("image", -1), Tcl_NewStringObj("create", -1),
                              Tcl_NewStringObj(kFrameImageType, -1), Tcl_NewStringObj(frame_image_name_.c_str(), -1),
                              Tcl_NewStringObj(Tk_PathName(window_), -1)}) != TCL_OK)
    {
        return TCL_ERROR;
    }
    const int result = eval_command(interp_, {Tcl_NewStringObj(canvas_command_.c_str(), -1),
                                             Tcl_NewStringObj("create", -1), Tcl_NewStringObj("image", -1),
                                             Tcl_NewIntObj(0), Tcl_NewIntObj(0),
                                             Tcl_NewStringObj("-anchor", -1), Tcl_NewStringObj("nw", -1),
                                             Tcl_NewStringObj("-image", -1),
                                             Tcl_NewStringObj(frame_image_name_.c_str(), -1)});
    if(result != TCL_OK)
    {
        destroy_frame_image();
    }
    return result;
}

void HtmlWidget::destroy_frame_image()
{
    if(!frame_image_name_.empty())
    {
        Tk_DeleteImage(interp_, frame_image_name_.c_str());
        frame_image_name_.clear();
    }
}

void HtmlWidget::destroy_backing_image()
{
    if(backing_image_ != nullptr)
    {
        Tk_FreeImage(backing_image_);
        backing_image_ = nullptr;
    }
    backing_photo_ = nullptr;
    if(backing_image_name_.empty())
    {
        return;
    }
    Image image;
    image.name     = backing_image_name_;
    image.internal = true;
    delete_internal_image(image);
    backing_image_name_.clear();
}

int HtmlWidget::resize_canvas()
{
    canvas_width_  = std::max(Tk_Width(window_), 1);
    canvas_height_ = std::max(Tk_Height(window_), 1);
    const auto width = static_cast<std::size_t>(canvas_width_);
    const auto height = static_cast<std::size_t>(canvas_height_);
    if(width > max_viewport_pixels_ / height)
    {
        const std::string message = "widget viewport exceeds the " + std::to_string(max_viewport_pixels_) +
            " pixel limit";
        Tcl_SetObjResult(interp_, Tcl_NewStringObj(message.c_str(), -1));
        return TCL_ERROR;
    }
    pixels_.assign(width * height * PixelChannelCount, 0);
    return TCL_OK;
}

int HtmlWidget::upload_canvas()
{
    Tk_PhotoImageBlock block;
    block.pixelPtr  = pixels_.data();
    block.width     = canvas_width_;
    block.height    = canvas_height_;
    block.pitch     = canvas_width_ * PixelChannelCount;
    block.pixelSize = PixelChannelCount;
    block.offset[Red]   = Red;
    block.offset[Green] = Green;
    block.offset[Blue]  = Blue;
    block.offset[Alpha] = Alpha;
    return Tk_PhotoPutBlock(interp_, backing_photo_, &block, 0, 0, canvas_width_, canvas_height_,
                            TK_PHOTO_COMPOSITE_SET);
}

void HtmlWidget::layout_document()
{
    if(document_ == nullptr || window_ == nullptr)
    {
        return;
    }
    const int width = std::max(Tk_Width(window_), 1);
    document_->render(width);
    document_width_  = std::max(to_int(document_->width()), width);
    document_height_ = std::max(to_int(document_->height()), Tk_Height(window_));
    clamp_scroll();
    invalidate_text_index();
    if(!loading_)
    {
        update_scroll_commands();
    }
}

int HtmlWidget::ensure_text_index()
{
    if(layout_pending_ && document_ != nullptr && window_ != nullptr)
    {
        layout_pending_ = false;
        layout_document();
    }
    if(text_index_valid_ || document_ == nullptr)
    {
        return TCL_OK;
    }
    try
    {
        build_text_index();
    } catch(...)
    {
        invalidate_text_index();
        return exception_result(interp_, "text indexing");
    }
    return TCL_OK;
}

void HtmlWidget::build_text_index()
{
    text_runs_.clear();
    plain_text_.clear();
    text_character_count_ = 0;
    if(document_ == nullptr)
    {
        text_index_valid_ = true;
        return;
    }
    indexing_ = true;
    try
    {
        char separator = '\0';
        index_element(document_->root(), separator);
        indexing_ = false;
        text_run_lookup_.resize(text_runs_.size());
        std::iota(text_run_lookup_.begin(), text_run_lookup_.end(), 0);
        std::sort(text_run_lookup_.begin(), text_run_lookup_.end(), [this](std::size_t left, std::size_t right) {
            return text_run_location(text_runs_[left]) < text_run_location(text_runs_[right]);
        });
        text_index_valid_ = true;
    } catch(...)
    {
        indexing_ = false;
        throw;
    }
}

void HtmlWidget::index_element(const litehtml::element::ptr& element, char& separator)
{
    if(element->css().get_display() == litehtml::display_none || element->is_comment())
    {
        return;
    }
    if(element->is_break())
    {
        plain_text_ += '\n';
        ++text_character_count_;
        separator = '\0';
        return;
    }
    if(element->is_text())
    {
        const auto parent = element->parent();
        if(parent == nullptr || parent->css().get_visibility() != litehtml::visibility_visible)
        {
            return;
        }
        if(element->is_white_space())
        {
            if(separator == '\0')
            {
                separator = ' ';
            }
            return;
        }
        const auto render = element->get_render_item();
        if(render == nullptr)
        {
            return;
        }
        if(separator != '\0' && !plain_text_.empty() && plain_text_.back() != '\n')
        {
            plain_text_ += separator;
            ++text_character_count_;
        }
        separator = '\0';
        // Traverse in document order, but let litehtml supply the exact
        // displayed text (including transforms) and current geometry.
        const auto placement = render->get_placement();
        element->draw(0, placement.x - render->pos().x, placement.y - render->pos().y,
                      nullptr, render);
        return;
    }
    const bool block = !element->is_inline();
    if(block)
    {
        separator = '\n';
    }
    for(const auto& child : element->children())
    {
        index_element(child, separator);
    }
    if(block)
    {
        separator = '\n';
    }
}

void HtmlWidget::invalidate_text_index()
{
    text_runs_.clear();
    text_run_lookup_.clear();
    plain_text_.clear();
    text_character_count_ = 0;
    text_index_valid_ = false;
}

bool HtmlWidget::has_active_highlights() const
{
    return std::any_of(highlights_.begin(), highlights_.end(), [](const auto& item) {
        return !item.second.ranges.empty();
    });
}

void HtmlWidget::clamp_scroll()
{
    const int viewport_width  = window_ != nullptr ? Tk_Width(window_) : requested_width_;
    const int viewport_height = window_ != nullptr ? Tk_Height(window_) : requested_height_;
    scroll_x_ = std::clamp(scroll_x_, 0, std::max(document_width_ - viewport_width, 0));
    scroll_y_ = std::clamp(scroll_y_, 0, std::max(document_height_ - viewport_height, 0));
}

void HtmlWidget::set_scroll(bool horizontal, int value)
{
    (horizontal ? scroll_x_ : scroll_y_) = value;
    clamp_scroll();
    update_scroll_commands();
    schedule_redraw();
}

void HtmlWidget::update_scroll_commands()
{
    if(window_ == nullptr)
    {
        return;
    }
    invoke_scroll_command(xscroll_command_, scroll_x_, document_width_, Tk_Width(window_));
    invoke_scroll_command(yscroll_command_, scroll_y_, document_height_, Tk_Height(window_));
}

void HtmlWidget::invoke_scroll_command(Tcl_Obj* command, int offset, int total, int viewport)
{
    if(command == nullptr || Tcl_GetCharLength(command) == 0)
    {
        return;
    }
    const double first = total > 0 ? static_cast<double>(offset) / total : 0.0;
    const double last  = total > 0 ? static_cast<double>(std::min(offset + viewport, total)) / total : 1.0;
    Tcl_Obj* call = Tcl_DuplicateObj(command);
    Tcl_IncrRefCount(call);
    if(Tcl_ListObjAppendElement(interp_, call, Tcl_NewDoubleObj(first)) == TCL_OK &&
       Tcl_ListObjAppendElement(interp_, call, Tcl_NewDoubleObj(last)) == TCL_OK)
    {
        const int code = Tcl_EvalObjEx(interp_, call, TCL_EVAL_GLOBAL);
        if(code != TCL_OK)
        {
            Tcl_BackgroundException(interp_, code);
        }
    }
    Tcl_DecrRefCount(call);
}

int HtmlWidget::view_command(bool horizontal, int objc, Tcl_Obj* const objv[])
{
    const int total    = horizontal ? document_width_ : document_height_;
    const int viewport = window_ != nullptr ? (horizontal ? Tk_Width(window_) : Tk_Height(window_)) : 1;
    const int offset   = horizontal ? scroll_x_ : scroll_y_;
    if(objc == 0)
    {
        Tcl_Obj* fractions[2] = {
            Tcl_NewDoubleObj(total > 0 ? static_cast<double>(offset) / total : 0.0),
            Tcl_NewDoubleObj(total > 0 ? static_cast<double>(std::min(offset + viewport, total)) / total : 1.0)};
        Tcl_SetObjResult(interp_, Tcl_NewListObj(2, fractions));
        return TCL_OK;
    }

    static const char* const operations[] = {"moveto", "scroll", nullptr};
    int                      operation    = 0;
    if(Tcl_GetIndexFromObj(interp_, objv[0], operations, "operation", 0, &operation) != TCL_OK)
    {
        return TCL_ERROR;
    }
    if(operation == 0)
    {
        if(objc != 2)
        {
            return wrong_args(horizontal ? "xview moveto fraction" : "yview moveto fraction");
        }
        double fraction = 0;
        if(Tcl_GetDoubleFromObj(interp_, objv[1], &fraction) != TCL_OK)
        {
            return TCL_ERROR;
        }
        set_scroll(horizontal, static_cast<int>(std::lround(fraction * total)));
        return TCL_OK;
    }
    if(objc != 3)
    {
        return wrong_args(horizontal ? "xview scroll number units|pages" : "yview scroll number units|pages");
    }
    int count = 0;
    if(Tcl_GetIntFromObj(interp_, objv[1], &count) != TCL_OK)
    {
        return TCL_ERROR;
    }
    const std::string units = Tcl_GetString(objv[2]);
    int               step  = 0;
    if(units == "units")
    {
        step = kScrollUnitPixels;
    } else if(units == "pages")
    {
        step = std::max(viewport - kScrollUnitPixels, 1);
    } else
    {
        Tcl_SetObjResult(interp_, Tcl_NewStringObj("bad argument: must be units or pages", -1));
        return TCL_ERROR;
    }
    set_scroll(horizontal, offset + count * step);
    return TCL_OK;
}

std::vector<TextRange> HtmlWidget::search(std::string_view pattern, bool no_case) const
{
    std::vector<TextRange> matches;
    auto                   haystack = unicode_chars(plain_text_);
    auto                   needle   = unicode_chars(pattern);
    if(needle.empty())
    {
        return matches;
    }
    if(no_case)
    {
        std::transform(haystack.begin(), haystack.end(), haystack.begin(), Tcl_UniCharToLower);
        std::transform(needle.begin(), needle.end(), needle.begin(), Tcl_UniCharToLower);
    }
    const auto searcher = std::boyer_moore_searcher(needle.begin(), needle.end());
    for(auto start = haystack.begin(); start != haystack.end();)
    {
        const auto found = std::search(start, haystack.end(), searcher);
        if(found == haystack.end())
        {
            break;
        }
        const auto first = static_cast<int>(std::distance(haystack.begin(), found));
        matches.push_back({first, first + static_cast<int>(needle.size())});
        start = found + static_cast<std::ptrdiff_t>(needle.size());
    }
    return matches;
}

int HtmlWidget::search_command(int objc, Tcl_Obj* const objv[])
{
    bool all     = false;
    bool no_case = false;
    int  index   = 0;
    for(; index < objc; ++index)
    {
        const std::string option = Tcl_GetString(objv[index]);
        if(option == "--")
        {
            ++index;
            break;
        }
        if(option == "-all")
        {
            all = true;
        } else if(option == "-nocase")
        {
            no_case = true;
        } else if(!option.empty() && option.front() == '-')
        {
            Tcl_SetObjResult(interp_, Tcl_ObjPrintf("unknown search option \"%s\"", option.c_str()));
            return TCL_ERROR;
        } else
        {
            break;
        }
    }
    if(index != objc - 1)
    {
        return wrong_args("search ?-all? ?-nocase? ?--? pattern");
    }
    if(ensure_text_index() != TCL_OK)
    {
        return TCL_ERROR;
    }
    return ranges_result(search(Tcl_GetString(objv[index]), no_case), all);
}

int HtmlWidget::ranges_result(const std::vector<TextRange>& ranges, bool all) const
{
    Tcl_Obj* result = Tcl_NewListObj(0, nullptr);
    const auto append = [this, result](const TextRange& range) {
        Tcl_Obj* pair[2] = {Tcl_NewIntObj(range.first), Tcl_NewIntObj(range.last)};
        Tcl_ListObjAppendElement(interp_, result, Tcl_NewListObj(2, pair));
    };
    if(all)
    {
        for(const auto& range : ranges)
        {
            append(range);
        }
    } else if(!ranges.empty())
    {
        append(ranges.front());
    }
    Tcl_SetObjResult(interp_, result);
    return TCL_OK;
}

int HtmlWidget::parse_ranges(Tcl_Obj* object, std::vector<TextRange>& ranges) const
{
    int count = 0;
    if(Tcl_ListObjLength(interp_, object, &count) != TCL_OK)
    {
        return TCL_ERROR;
    }
    const int text_length = text_character_count_;
    for(int index = 0; index < count; ++index)
    {
        Tcl_Obj* pair = nullptr;
        int      pair_length = 0;
        Tcl_ListObjIndex(interp_, object, index, &pair);
        if(pair == nullptr || Tcl_ListObjLength(interp_, pair, &pair_length) != TCL_OK || pair_length != 2)
        {
            Tcl_SetObjResult(interp_, Tcl_NewStringObj("each highlight range must contain a first and last index", -1));
            return TCL_ERROR;
        }
        Tcl_Obj* first_obj = nullptr;
        Tcl_Obj* last_obj  = nullptr;
        int      first     = 0;
        int      last      = 0;
        Tcl_ListObjIndex(interp_, pair, 0, &first_obj);
        Tcl_ListObjIndex(interp_, pair, 1, &last_obj);
        if(Tcl_GetIntFromObj(interp_, first_obj, &first) != TCL_OK ||
           Tcl_GetIntFromObj(interp_, last_obj, &last) != TCL_OK)
        {
            return TCL_ERROR;
        }
        if(first < 0 || last < first || last > text_length)
        {
            Tcl_SetObjResult(interp_, Tcl_NewStringObj("highlight range is outside the document text", -1));
            return TCL_ERROR;
        }
        ranges.push_back({first, last});
    }
    std::sort(ranges.begin(), ranges.end(), [](const TextRange& left, const TextRange& right) {
        return left.first < right.first || (left.first == right.first && left.last < right.last);
    });
    std::vector<TextRange> normalized;
    normalized.reserve(ranges.size());
    for(const auto& range : ranges)
    {
        if(range.first == range.last)
        {
            continue;
        }
        if(normalized.empty() || range.first > normalized.back().last)
        {
            normalized.push_back(range);
        } else
        {
            normalized.back().last = std::max(normalized.back().last, range.last);
        }
    }
    ranges = std::move(normalized);
    return TCL_OK;
}

Highlight& HtmlWidget::ensure_highlight(const std::string& name)
{
    auto [iterator, inserted] = highlights_.try_emplace(name);
    if(inserted)
    {
        const char* color = name == "selection" ? kSelectionBackground :
            name == "find-current" ? kCurrentBackground : kFindBackground;
        set_highlight_color(iterator->second, color);
    }
    return iterator->second;
}

int HtmlWidget::set_highlight_color(Highlight& highlight, const char* color_name)
{
    XColor* color = Tk_GetColor(interp_, window_, Tk_GetUid(color_name));
    if(color == nullptr)
    {
        return TCL_ERROR;
    }
    highlight.background = litehtml::web_color(static_cast<litehtml::byte>(color->red >> 8),
                                                static_cast<litehtml::byte>(color->green >> 8),
                                                static_cast<litehtml::byte>(color->blue >> 8));
    Tk_FreeColor(color);
    return TCL_OK;
}

int HtmlWidget::highlight_command(int objc, Tcl_Obj* const objv[])
{
    static const char* const commands[] = {"clear", "configure", "set", nullptr};
    if(objc < 1)
    {
        return wrong_args("highlight clear|configure|set ?arg ...?");
    }
    int command = 0;
    if(Tcl_GetIndexFromObj(interp_, objv[0], commands, "highlight operation", 0, &command) != TCL_OK)
    {
        return TCL_ERROR;
    }
    if(command == 0)
    {
        if(objc > 2)
        {
            return wrong_args("highlight clear ?name?");
        }
        if(objc == 1)
        {
            for(auto& [name, highlight] : highlights_)
            {
                highlight.ranges.clear();
            }
        } else
        {
            const auto iterator = highlights_.find(Tcl_GetString(objv[1]));
            if(iterator != highlights_.end())
            {
                iterator->second.ranges.clear();
            }
        }
        schedule_redraw();
        return TCL_OK;
    }
    if(command == 1)
    {
        if(objc != 4 || std::string(Tcl_GetString(objv[2])) != "-background")
        {
            return wrong_args("highlight configure name -background color");
        }
        auto& highlight = ensure_highlight(Tcl_GetString(objv[1]));
        if(set_highlight_color(highlight, Tcl_GetString(objv[3])) != TCL_OK)
        {
            return TCL_ERROR;
        }
        schedule_redraw();
        return TCL_OK;
    }
    if(objc != 3)
    {
        return wrong_args("highlight set name ranges");
    }
    if(ensure_text_index() != TCL_OK)
    {
        return TCL_ERROR;
    }
    std::vector<TextRange> ranges;
    if(parse_ranges(objv[2], ranges) != TCL_OK)
    {
        return TCL_ERROR;
    }
    ensure_highlight(Tcl_GetString(objv[1])).ranges = std::move(ranges);
    schedule_redraw();
    return TCL_OK;
}

int HtmlWidget::find_command(int objc, Tcl_Obj* const objv[])
{
    static const char* const commands[] = {"clear", "current", "next", "previous", "set", nullptr};
    if(objc < 1)
    {
        return wrong_args("find clear|current|next|previous|set ?arg ...?");
    }
    int command = 0;
    if(Tcl_GetIndexFromObj(interp_, objv[0], commands, "find operation", 0, &command) != TCL_OK)
    {
        return TCL_ERROR;
    }
    if(command == 0)
    {
        if(objc != 1)
        {
            return wrong_args("find clear");
        }
        find_ranges_.clear();
        current_find_ = -1;
        ensure_highlight("find").ranges.clear();
        ensure_highlight("find-current").ranges.clear();
        schedule_redraw();
        return TCL_OK;
    }
    if(ensure_text_index() != TCL_OK)
    {
        return TCL_ERROR;
    }
    if(command == 4)
    {
        bool no_case = true;
        int  index   = 1;
        if(index < objc - 1 && std::string(Tcl_GetString(objv[index])) == "-case")
        {
            no_case = false;
            ++index;
        }
        if(index != objc - 1)
        {
            return wrong_args("find set ?-case? pattern");
        }
        find_ranges_ = search(Tcl_GetString(objv[index]), no_case);
        ensure_highlight("find").ranges = find_ranges_;
        current_find_ = find_ranges_.empty() ? -1 : 0;
    } else if((command == 1 || command == 2 || command == 3) && objc != 1)
    {
        const char* operation = command == 1 ? "current" : command == 2 ? "next" : "previous";
        return wrong_args((std::string("find ") + operation).c_str());
    } else if(command == 2 && !find_ranges_.empty())
    {
        current_find_ = (current_find_ + 1) % static_cast<int>(find_ranges_.size());
    } else if(command == 3 && !find_ranges_.empty())
    {
        current_find_ = (current_find_ - 1 + static_cast<int>(find_ranges_.size())) %
                        static_cast<int>(find_ranges_.size());
    } else if(command == 1)
    {
        if(current_find_ < 0)
        {
            Tcl_ResetResult(interp_);
        } else
        {
            return ranges_result({find_ranges_[static_cast<std::size_t>(current_find_)]}, false);
        }
        return TCL_OK;
    }

    auto& current = ensure_highlight("find-current");
    current.ranges.clear();
    if(current_find_ >= 0)
    {
        current.ranges.push_back(find_ranges_[static_cast<std::size_t>(current_find_)]);
        const auto* run = run_at(current.ranges.front().first);
        if(run != nullptr)
        {
            set_scroll(false, to_int(run->position.y) - kScrollUnitPixels);
        }
    }
    schedule_redraw();
    Tcl_SetObjResult(interp_, Tcl_NewIntObj(static_cast<int>(find_ranges_.size())));
    return TCL_OK;
}

TextRun* HtmlWidget::run_at(int index)
{
    return const_cast<TextRun*>(std::as_const(*this).run_at(index));
}

const TextRun* HtmlWidget::run_at(int index) const
{
    const auto iterator = std::lower_bound(text_runs_.begin(), text_runs_.end(), index,
                                           [](const TextRun& run, int value) { return run.range.last <= value; });
    return iterator != text_runs_.end() ? &*iterator : nullptr;
}

int HtmlWidget::index_from_obj(Tcl_Obj* object, int& index) const
{
    const std::string value = Tcl_GetString(object);
    const int text_length = text_character_count_;
    if(!value.empty() && value.front() == '@')
    {
        const auto comma = value.find(',');
        int x = 0;
        int y = 0;
        if(comma != std::string::npos)
        {
            const char* begin = value.data();
            const char* end = begin + value.size();
            const auto horizontal = std::from_chars(begin + 1, begin + comma, x);
            const auto vertical = std::from_chars(begin + comma + 1, end, y);
            if(horizontal.ec == std::errc{} && horizontal.ptr == begin + comma &&
               vertical.ec == std::errc{} && vertical.ptr == end)
            {
                index = index_at(x, y);
                return TCL_OK;
            }
        }
        Tcl_SetObjResult(interp_, Tcl_NewStringObj("bad coordinate index: must be @x,y", -1));
        return TCL_ERROR;
    }
    if(value == "end")
    {
        index = text_length;
        return TCL_OK;
    }
    if(Tcl_GetIntFromObj(interp_, object, &index) != TCL_OK)
    {
        return TCL_ERROR;
    }
    if(index < 0 || index > text_length)
    {
        Tcl_SetObjResult(interp_, Tcl_NewStringObj("text index is outside the document", -1));
        return TCL_ERROR;
    }
    return TCL_OK;
}

int HtmlWidget::index_at(int x, int y) const
{
    const long long document_x = static_cast<long long>(x) + scroll_x_;
    const long long document_y = static_cast<long long>(y) + scroll_y_;
    const TextRun* nearest = nullptr;
    auto best_distance = std::make_pair(std::numeric_limits<long long>::max(),
                                       std::numeric_limits<long long>::max());
    for(const auto& run : text_runs_)
    {
        const long long left = to_int(run.position.x);
        const long long top = to_int(run.position.y);
        const long long right = left + Tk_TextWidth(run.font->tk_font, run.text.data(),
                                                    static_cast<int>(run.text.size()));
        const long long bottom = top + to_int(run.position.height);
        const auto distance = std::make_pair(
            std::max({top - document_y, document_y - bottom, 0LL}),
            std::max({left - document_x, document_x - right, 0LL}));
        if(distance < best_distance)
        {
            nearest = &run;
            best_distance = distance;
        }
    }
    if(nearest == nullptr)
    {
        return 0;
    }
    const long long offset = document_x - to_int(nearest->position.x);
    if(offset <= 0)
    {
        return nearest->range.first;
    }
    if(offset >= std::numeric_limits<int>::max())
    {
        return nearest->range.last;
    }
    int width = 0;
    const int bytes = Tk_MeasureChars(nearest->font->tk_font, nearest->text.data(),
                                      static_cast<int>(nearest->text.size()),
                                      static_cast<int>(offset), 0, &width);
    int index = nearest->range.first + Tcl_NumUtfChars(nearest->text.data(), bytes);
    if(index < nearest->range.last)
    {
        const char* character = nearest->text.data() + bytes;
        const int character_width = Tk_TextWidth(nearest->font->tk_font, character,
                                                  static_cast<int>(Tcl_UtfNext(character) - character));
        if(offset - width >= (character_width + 1) / 2)
        {
            ++index;
        }
    }
    return index;
}

int HtmlWidget::selection_command(int objc, Tcl_Obj* const objv[])
{
    static const char* const operations[] = {"clear", "get", "range", "set", nullptr};
    enum Operation { ClearSelection, GetSelection, SelectionRange, SetSelection };
    int operation = 0;
    if(objc < 1)
    {
        return wrong_args("selection clear|get|range|set ?first last?");
    }
    if(Tcl_GetIndexFromObj(interp_, objv[0], operations, "selection operation", 0, &operation) != TCL_OK)
    {
        return TCL_ERROR;
    }
    if(objc != (operation == SetSelection ? 3 : 1))
    {
        return wrong_args(operation == SetSelection ? "selection set first last" : "selection clear|get|range");
    }
    if(ensure_text_index() != TCL_OK)
    {
        return TCL_ERROR;
    }
    switch(static_cast<Operation>(operation))
    {
    case ClearSelection:
        return set_selection(0, 0);
    case GetSelection:
    {
        const char* first = utf_at(plain_text_, selection_.first);
        const char* last = utf_at(plain_text_, selection_.last);
        Tcl_SetObjResult(interp_, Tcl_NewStringObj(first, static_cast<int>(last - first)));
        return TCL_OK;
    }
    case SelectionRange:
        if(selection_.first != selection_.last)
        {
            Tcl_Obj* range[] = {Tcl_NewIntObj(selection_.first), Tcl_NewIntObj(selection_.last)};
            Tcl_SetObjResult(interp_, Tcl_NewListObj(2, range));
        }
        return TCL_OK;
    case SetSelection:
    {
        int first = 0;
        int last = 0;
        if(index_from_obj(objv[1], first) != TCL_OK || index_from_obj(objv[2], last) != TCL_OK)
        {
            return TCL_ERROR;
        }
        return set_selection(std::min(first, last), std::max(first, last));
    }
    }
    return TCL_ERROR;
}

int HtmlWidget::set_selection(int first, int last)
{
    if(first == last)
    {
        first = last = 0;
    }
    if(selection_.first == first && selection_.last == last)
    {
        return TCL_OK;
    }
    selection_ = {first, last};
    if(first != last)
    {
        ensure_highlight("selection").ranges = {selection_};
    } else
    {
        const auto highlight = highlights_.find("selection");
        if(highlight != highlights_.end())
        {
            highlight->second.ranges.clear();
        }
    }
    schedule_redraw();
    return notify("<<Selection>>");
}

int HtmlWidget::notify(const char* event)
{
    if(window_ == nullptr)
    {
        return TCL_ERROR;
    }
    return eval_command(interp_, {Tcl_NewStringObj("event", -1), Tcl_NewStringObj("generate", -1),
                                  Tcl_NewStringObj(Tk_PathName(window_), -1), Tcl_NewStringObj(event, -1),
                                  Tcl_NewStringObj("-when", -1), Tcl_NewStringObj("tail", -1)});
}

int HtmlWidget::bbox_command(Tcl_Obj* index_obj)
{
    if(ensure_text_index() != TCL_OK)
    {
        return TCL_ERROR;
    }
    int index = 0;
    if(index_from_obj(index_obj, index) != TCL_OK)
    {
        return TCL_ERROR;
    }
    const TextRun* run = run_at(index);
    if(run == nullptr || index < run->range.first || index >= run->range.last)
    {
        Tcl_ResetResult(interp_);
        return TCL_OK;
    }
    const int relative = index - run->range.first;
    const char* start  = run->text.data();
    const char* middle = utf_at(run->text, relative);
    const char* next   = Tcl_UtfNext(middle);
    const int x = to_int(run->position.x) + Tk_TextWidth(run->font->tk_font, start, static_cast<int>(middle - start));
    const int width = Tk_TextWidth(run->font->tk_font, middle, static_cast<int>(next - middle));
    Tcl_Obj* values[4] = {Tcl_NewIntObj(x - scroll_x_), Tcl_NewIntObj(to_int(run->position.y) - scroll_y_),
                          Tcl_NewIntObj(width), Tcl_NewIntObj(to_int(run->position.height))};
    Tcl_SetObjResult(interp_, Tcl_NewListObj(4, values));
    return TCL_OK;
}

int HtmlWidget::see_command(Tcl_Obj* index_obj)
{
    if(ensure_text_index() != TCL_OK)
    {
        return TCL_ERROR;
    }
    int index = 0;
    if(index_from_obj(index_obj, index) != TCL_OK)
    {
        return TCL_ERROR;
    }
    const TextRun* run = run_at(index);
    if(run == nullptr || index < run->range.first)
    {
        return TCL_OK;
    }
    const int top    = to_int(run->position.y);
    const int bottom = to_int(run->position.bottom());
    if(top < scroll_y_)
    {
        set_scroll(false, top);
    } else if(bottom > scroll_y_ + Tk_Height(window_))
    {
        set_scroll(false, bottom - Tk_Height(window_));
    }
    return TCL_OK;
}

void HtmlWidget::paint_highlights(const TextRun& run, int screen_x, int screen_y)
{
    for(const auto& [name, highlight] : highlights_)
    {
        auto range = std::lower_bound(highlight.ranges.begin(), highlight.ranges.end(), run.range.first,
                                      [](const TextRange& candidate, int first) {
                                          return candidate.last <= first;
                                      });
        for(; range != highlight.ranges.end() && range->first < run.range.last; ++range)
        {
            const int first = std::max(range->first, run.range.first);
            const int last  = std::min(range->last, run.range.last);
            if(first >= last)
            {
                continue;
            }
            const char* begin = run.text.data();
            const char* left  = utf_at(run.text, first - run.range.first);
            const char* right = utf_at(run.text, last - run.range.first);
            const int   x = screen_x + Tk_TextWidth(run.font->tk_font, begin, static_cast<int>(left - begin));
            const int   width = Tk_TextWidth(run.font->tk_font, left, static_cast<int>(right - left));
            fill_rectangle(highlight.background, x, screen_y, width, to_int(run.position.height));
        }
    }
}

void HtmlWidget::clear_highlights()
{
    highlights_.clear();
}

TextColor* HtmlWidget::text_color_for(const litehtml::web_color& color)
{
    const std::uint32_t key = (static_cast<std::uint32_t>(color.red) << 16U) |
                              (static_cast<std::uint32_t>(color.green) << 8U) | color.blue;
    const auto existing = text_colors_.find(key);
    if(existing != text_colors_.end())
    {
        return &existing->second;
    }
    char name[8];
    std::snprintf(name, sizeof(name), "#%02x%02x%02x", color.red, color.green, color.blue);
    XColor* allocated = Tk_GetColor(interp_, window_, Tk_GetUid(name));
    if(allocated == nullptr)
    {
        Tcl_ResetResult(interp_);
        return nullptr;
    }
    XGCValues values{};
    values.foreground = allocated->pixel;
    TextColor text_color;
    text_color.color = allocated;
    text_color.gc    = Tk_GetGC(window_, GCForeground, &values);
    auto [iterator, inserted] = text_colors_.emplace(key, text_color);
    return inserted ? &iterator->second : nullptr;
}

void HtmlWidget::clear_text_colors()
{
    for(const auto& [key, color] : text_colors_)
    {
        if(color.gc != None)
        {
            Tk_FreeGC(display_, color.gc);
        }
        if(color.color != nullptr)
        {
            Tk_FreeColor(color.color);
        }
    }
    text_colors_.clear();
}

litehtml::position HtmlWidget::active_clip() const
{
    litehtml::position clip(0, 0, canvas_width_, canvas_height_);
    for(const auto& item : clip_stack_)
    {
        clip = clip.intersect(item);
    }
    return clip;
}

void HtmlWidget::blend_pixel(int x, int y, const litehtml::web_color& color)
{
    if(x < 0 || y < 0 || x >= canvas_width_ || y >= canvas_height_ || color.alpha == 0)
    {
        return;
    }
    const std::size_t offset =
        (static_cast<std::size_t>(y) * canvas_width_ + x) * PixelChannelCount;
    const unsigned int alpha = color.alpha;
    const unsigned int inverse = kColorComponentMaximum - alpha;
    pixels_[offset + Red] = static_cast<unsigned char>(
        (color.red * alpha + pixels_[offset + Red] * inverse) / kColorComponentMaximum);
    pixels_[offset + Green] = static_cast<unsigned char>(
        (color.green * alpha + pixels_[offset + Green] * inverse) / kColorComponentMaximum);
    pixels_[offset + Blue] = static_cast<unsigned char>(
        (color.blue * alpha + pixels_[offset + Blue] * inverse) / kColorComponentMaximum);
    pixels_[offset + Alpha] = static_cast<unsigned char>(
        alpha + pixels_[offset + Alpha] * inverse / kColorComponentMaximum);
}

void HtmlWidget::fill_rectangle(const litehtml::web_color& color, int x, int y, int width, int height)
{
    if(width <= 0 || height <= 0)
    {
        return;
    }
    const auto clip = active_clip();
    const int left   = std::max(x, to_int(clip.left()));
    const int top    = std::max(y, to_int(clip.top()));
    const int right  = std::min(saturating_add(x, width), to_int(clip.right()));
    const int bottom = std::min(saturating_add(y, height), to_int(clip.bottom()));
    for(int row = top; row < bottom; ++row)
    {
        for(int column = left; column < right; ++column)
        {
            blend_pixel(column, row, color);
        }
    }
}

void HtmlWidget::fill_ellipse(const litehtml::web_color& color, int x, int y, int width, int height, bool outline)
{
    if(width <= 0 || height <= 0)
    {
        return;
    }
    const double rx = width / 2.0;
    const double ry = height / 2.0;
    const double cx = x + rx;
    const double cy = y + ry;
    const auto clip = active_clip();
    const int top = std::max(y, to_int(clip.top()));
    const int bottom = std::min(saturating_add(y, height), to_int(clip.bottom()));
    const int left = std::max(x, to_int(clip.left()));
    const int right = std::min(saturating_add(x, width), to_int(clip.right()));
    for(int row = top; row < bottom; ++row)
    {
        for(int column = left; column < right; ++column)
        {
            const double dx = (column + 0.5 - cx) / rx;
            const double dy = (row + 0.5 - cy) / ry;
            const double distance = dx * dx + dy * dy;
            if((!outline && distance <= 1.0) ||
               (outline && distance <= 1.0 && distance >= kOutlineInnerRadiusSquared))
            {
                blend_pixel(column, row, color);
            }
        }
    }
}

void HtmlWidget::reset_clip()
{
    clip_stack_.clear();
}

std::string HtmlWidget::resolve_url(const char* source, const char* base) const
{
    if(source == nullptr || *source == '\0')
    {
        return {};
    }
    std::string url(source);
    if(url.front() == '#')
    {
        return url;
    }
    const litehtml::url reference(url);
#if defined(_WIN32)
    const bool windows_drive_path = url.size() >= 3 && std::isalpha(static_cast<unsigned char>(url[0])) != 0 &&
                                    url[1] == ':' && (url[2] == '/' || url[2] == '\\');
#else
    constexpr bool windows_drive_path = false;
#endif
    if(reference.has_scheme() && reference.scheme() != "file" && !windows_drive_path)
    {
        return url;
    }
    if(reference.scheme() == "file")
    {
        url = litehtml::url::decode(reference.path());
#if defined(_WIN32)
        if(url.size() >= 4 && url[0] == '/' && std::isalpha(static_cast<unsigned char>(url[1])) != 0 &&
           url[2] == ':')
        {
            url.erase(0, 1);
        }
#endif
    }
    std::string base_string = base != nullptr && *base != '\0' ? base : document_base_url_;
    const litehtml::url base_reference(base_string);
#if defined(_WIN32)
    const bool windows_drive_base = base_string.size() >= 3 &&
                                    std::isalpha(static_cast<unsigned char>(base_string[0])) != 0 &&
                                    base_string[1] == ':' && (base_string[2] == '/' || base_string[2] == '\\');
#else
    constexpr bool windows_drive_base = false;
#endif
    if(base_reference.has_scheme() && base_reference.scheme() != "file" && !windows_drive_base)
    {
        return litehtml::resolve(base_reference, reference).str();
    }
    std::filesystem::path path = path_from_utf8(litehtml::url::decode(url));
    if(path.is_absolute())
    {
        return path.lexically_normal().u8string();
    }
    if(base_string.rfind(kFileScheme, 0) == 0)
    {
        base_string.erase(0, std::strlen(kFileScheme));
        base_string = litehtml::url::decode(base_string);
    }
    std::filesystem::path base_path = path_from_utf8(base_string);
    std::error_code       filesystem_error;
    const bool            base_is_directory = std::filesystem::is_directory(base_path, filesystem_error);
    if(!base_string.empty() && base_string.back() != '/' && !base_is_directory)
    {
        base_path = base_path.parent_path();
    }
    return (base_path / path).lexically_normal().u8string();
}

std::string HtmlWidget::internal_image_name(std::string_view purpose, unsigned int serial) const
{
    std::string name = "::tklitehtml::";
    name.append(purpose);
    name.push_back('_');
    const std::string path = window_ != nullptr ? Tk_PathName(window_) : "widget";
    for(unsigned char character : path)
    {
        name.push_back(std::isalnum(character) != 0 ? static_cast<char>(character) : '_');
    }
    name.push_back('_');
    name += std::to_string(serial);
    return name;
}

std::string HtmlWidget::resource(const char* kind, const std::string& url)
{
    if(resource_command_ == nullptr || Tcl_GetCharLength(resource_command_) == 0)
    {
        return {};
    }
    Tcl_Obj* call = Tcl_DuplicateObj(resource_command_);
    Tcl_IncrRefCount(call);
    Tcl_ListObjAppendElement(interp_, call, Tcl_NewStringObj(Tk_PathName(window_), -1));
    Tcl_ListObjAppendElement(interp_, call, Tcl_NewStringObj(kind, -1));
    Tcl_ListObjAppendElement(interp_, call, Tcl_NewStringObj(url.c_str(), static_cast<int>(url.size())));
    std::string result;
    const int code = Tcl_EvalObjEx(interp_, call, TCL_EVAL_GLOBAL);
    if(code == TCL_OK)
    {
        int         length = 0;
        const char* bytes  = Tcl_GetStringFromObj(Tcl_GetObjResult(interp_), &length);
        result.assign(bytes, static_cast<std::size_t>(length));
    } else
    {
        const std::string message = Tcl_GetStringResult(interp_);
        Tcl_DecrRefCount(call);
        Tcl_ResetResult(interp_);
        throw std::runtime_error(std::string(kind) + " resource callback failed for '" + url + "': " + message);
    }
    Tcl_DecrRefCount(call);
    if(window_ == nullptr)
    {
        Tcl_ResetResult(interp_);
        throw std::runtime_error("widget was destroyed by a resource callback");
    }
    Tcl_ResetResult(interp_);
    return result;
}

Image* HtmlWidget::image_for(const std::string& url)
{
    const auto found = images_.find(url);
    if(found != images_.end())
    {
        return &found->second;
    }

    Image image;
    image.name = resource("image", url);
    std::error_code filesystem_error;
    std::size_t file_bytes = 0;
    if(image.name.empty() && std::filesystem::is_regular_file(path_from_utf8(url), filesystem_error))
    {
        const auto size = std::filesystem::file_size(path_from_utf8(url), filesystem_error);
        if(!filesystem_error)
        {
            if(size > max_document_bytes_ - document_resource_bytes_)
            {
                throw std::runtime_error("image resource exceeds the document byte limit");
            }
            file_bytes = static_cast<std::size_t>(size);
        }
        image.name     = internal_image_name("resource", ++image_serial_);
        image.internal = true;
        Tcl_Obj* arguments[6] = {Tcl_NewStringObj("image", -1), Tcl_NewStringObj("create", -1),
                                 Tcl_NewStringObj("photo", -1), Tcl_NewStringObj(image.name.c_str(), -1),
                                 Tcl_NewStringObj("-file", -1), Tcl_NewStringObj(url.c_str(), -1)};
        for(auto* argument : arguments)
        {
            Tcl_IncrRefCount(argument);
        }
        const int code = Tcl_EvalObjv(interp_, 6, arguments, TCL_EVAL_GLOBAL);
        for(auto* argument : arguments)
        {
            Tcl_DecrRefCount(argument);
        }
        if(code != TCL_OK)
        {
            image.name.clear();
            image.internal = false;
            Tcl_ResetResult(interp_);
        }
    }
    if(image.name.empty())
    {
        return nullptr;
    }
    image.handle = Tk_GetImage(interp_, window_, image.name.c_str(), image_changed, this);
    if(image.handle == nullptr)
    {
        delete_internal_image(image);
        Tcl_ResetResult(interp_);
        return nullptr;
    }
    Tk_SizeOfImage(image.handle, &image.width, &image.height);
    if(!dimensions_within_limit(image.width, image.height, max_image_pixels_ - image_pixels_))
    {
        Tk_FreeImage(image.handle);
        image.handle = nullptr;
        delete_internal_image(image);
        const std::string message = "image '" + url + "' exceeds the " + std::to_string(max_image_pixels_) +
            " pixel limit";
        Tcl_SetObjResult(interp_, Tcl_NewStringObj(message.c_str(), -1));
        throw std::runtime_error(Tcl_GetStringResult(interp_));
    }
    image.pixels = static_cast<std::size_t>(image.width) * image.height;
    auto [iterator, inserted] = images_.emplace(url, std::move(image));
    if(inserted)
    {
        image_pixels_ += iterator->second.pixels;
        document_resource_bytes_ += file_bytes;
    }
    return inserted ? &iterator->second : nullptr;
}

void HtmlWidget::delete_internal_image(const Image& image)
{
    if(!image.internal || image.name.empty())
    {
        return;
    }
    Tcl_Obj* arguments[3] = {Tcl_NewStringObj("image", -1), Tcl_NewStringObj("delete", -1),
                             Tcl_NewStringObj(image.name.c_str(), -1)};
    for(auto* argument : arguments)
    {
        Tcl_IncrRefCount(argument);
    }
    Tcl_EvalObjv(interp_, 3, arguments, TCL_EVAL_GLOBAL);
    for(auto* argument : arguments)
    {
        Tcl_DecrRefCount(argument);
    }
    Tcl_ResetResult(interp_);
}

void HtmlWidget::clear_images()
{
    clear_images(images_);
    image_pixels_ = 0;
}

void HtmlWidget::clear_images(std::map<std::string, Image>& images)
{
    for(const auto& [url, image] : images)
    {
        if(image.handle != nullptr)
        {
            Tk_FreeImage(image.handle);
        }
        delete_internal_image(image);
    }
    images.clear();
}

void HtmlWidget::invoke_link(const std::string& url, const litehtml::element::ptr& element)
{
    if(link_command_ == nullptr || Tcl_GetCharLength(link_command_) == 0)
    {
        return;
    }
    Tcl_Obj* attributes = Tcl_NewDictObj();
    for(const char* name : {"href", "title", "id", "class"})
    {
        const char* value = element != nullptr ? element->get_attr(name) : nullptr;
        if(value != nullptr)
        {
            Tcl_DictObjPut(interp_, attributes, Tcl_NewStringObj(name, -1), Tcl_NewStringObj(value, -1));
        }
    }
    Tcl_Obj* call = Tcl_DuplicateObj(link_command_);
    Tcl_IncrRefCount(call);
    Tcl_ListObjAppendElement(interp_, call, Tcl_NewStringObj(Tk_PathName(window_), -1));
    Tcl_ListObjAppendElement(interp_, call, Tcl_NewStringObj(url.c_str(), static_cast<int>(url.size())));
    Tcl_ListObjAppendElement(interp_, call, attributes);
    const int code = Tcl_EvalObjEx(interp_, call, TCL_EVAL_GLOBAL);
    if(code != TCL_OK)
    {
        Tcl_BackgroundException(interp_, code);
    }
    Tcl_DecrRefCount(call);
}

void HtmlWidget::show_fragment(std::string_view fragment)
{
    if(document_ == nullptr || fragment.empty())
    {
        return;
    }
    auto element = document_->root()->select_one("#" + std::string(fragment));
    if(element != nullptr)
    {
        set_scroll(false, to_int(element->get_placement().y));
    }
}

litehtml::uint_ptr HtmlWidget::create_font(const litehtml::font_description& description,
                                           const litehtml::document*, litehtml::font_metrics* metrics)
{
    Tcl_Obj* specification = Tcl_NewListObj(0, nullptr);
    Tcl_IncrRefCount(specification);
    const std::string family = trim_font_family(description.family);
    const auto append = [this, specification](const char* option, Tcl_Obj* value) {
        Tcl_ListObjAppendElement(interp_, specification, Tcl_NewStringObj(option, -1));
        Tcl_ListObjAppendElement(interp_, specification, value);
    };
    append("-family", Tcl_NewStringObj(family.c_str(), static_cast<int>(family.size())));
    append("-size", Tcl_NewIntObj(-std::max(to_int(description.size), 1)));
    append("-weight", Tcl_NewStringObj(description.weight >= kBoldWeightMinimum ? "bold" : "normal", -1));
    append("-slant", Tcl_NewStringObj(description.style == litehtml::font_style_italic ? "italic" : "roman", -1));
    append("-underline",
           Tcl_NewBooleanObj((description.decoration_line & litehtml::text_decoration_line_underline) != 0));
    append("-overstrike",
           Tcl_NewBooleanObj((description.decoration_line & litehtml::text_decoration_line_line_through) != 0));

    Tk_Font tk_font = Tk_AllocFontFromObj(interp_, window_, specification);
    Tcl_DecrRefCount(specification);
    if(tk_font == nullptr)
    {
        Tcl_ResetResult(interp_);
        return 0;
    }

    auto font         = std::make_unique<Font>();
    font->tk_font     = tk_font;
    font->description = description;
    Tk_GetFontMetrics(tk_font, &font->metrics);
    if(metrics != nullptr)
    {
        metrics->font_size  = description.size;
        metrics->height     = font->metrics.linespace;
        metrics->ascent     = font->metrics.ascent;
        metrics->descent    = font->metrics.descent;
        metrics->x_height   = std::max(font->metrics.ascent * 2 / 3, 1);
        metrics->ch_width   = Tk_TextWidth(tk_font, "0", 1);
        metrics->draw_spaces = true;
        metrics->sub_shift   = description.size / 5;
        metrics->super_shift = description.size / 3;
    }
    return reinterpret_cast<litehtml::uint_ptr>(font.release());
}

void HtmlWidget::delete_font(litehtml::uint_ptr handle)
{
    auto* font = reinterpret_cast<Font*>(handle);
    if(font != nullptr)
    {
        Tk_FreeFont(font->tk_font);
        delete font;
    }
}

litehtml::pixel_t HtmlWidget::text_width(const char* text, litehtml::uint_ptr handle)
{
    const auto* font = reinterpret_cast<Font*>(handle);
    return font != nullptr && text != nullptr ? Tk_TextWidth(font->tk_font, text, static_cast<int>(std::strlen(text)))
                                               : 0;
}

void HtmlWidget::draw_text(litehtml::uint_ptr, const char* text, litehtml::uint_ptr handle,
                           litehtml::web_color color, const litehtml::position& position)
{
    auto* font = reinterpret_cast<Font*>(handle);
    if(font == nullptr || text == nullptr)
    {
        return;
    }
    if(indexing_)
    {
        const int first = text_character_count_;
        const int count = Tcl_NumUtfChars(text, static_cast<int>(std::strlen(text)));
        text_runs_.push_back({{first, first + count}, text, position, font, color});
        plain_text_ += text;
        text_character_count_ += count;
        return;
    }
    if(color.alpha == 0)
    {
        return;
    }

    const int screen_x = to_int(position.x);
    const int screen_y = to_int(position.y);
    const int doc_x    = saturating_add(screen_x, scroll_x_);
    const int doc_y    = saturating_add(screen_y, scroll_y_);
    if(has_active_highlights())
    {
        const std::pair<int, int> location{doc_x, doc_y};
        auto candidate = std::lower_bound(
            text_run_lookup_.begin(), text_run_lookup_.end(), location,
            [this](std::size_t index, const std::pair<int, int>& value) {
                return text_run_location(text_runs_[index]) < value;
            });
        for(; candidate != text_run_lookup_.end(); ++candidate)
        {
            const auto& run = text_runs_[*candidate];
            if(text_run_location(run) != location)
            {
                break;
            }
            if(run.font == font && run.text == text)
            {
                paint_highlights(run, screen_x, screen_y);
                break;
            }
        }
    }
    paint_text_.push_back({text, font, color, screen_x, screen_y});
}

litehtml::pixel_t HtmlWidget::pt_to_px(float points) const
{
    return static_cast<float>(points * kCssPixelsPerPoint);
}

litehtml::pixel_t HtmlWidget::get_default_font_size() const
{
    return kDefaultFontSize;
}

const char* HtmlWidget::get_default_font_name() const
{
    return "sans-serif";
}

void HtmlWidget::draw_list_marker(litehtml::uint_ptr, const litehtml::list_marker& marker)
{
    if(indexing_)
    {
        return;
    }
    const int x      = to_int(marker.pos.x);
    const int y      = to_int(marker.pos.y);
    const int width  = std::max(to_int(marker.pos.width), 1);
    const int height = std::max(to_int(marker.pos.height), 1);
    switch(marker.marker_type)
    {
    case litehtml::list_style_type_circle:
        fill_ellipse(marker.color, x, y, width, height, true);
        break;
    case litehtml::list_style_type_disc:
        fill_ellipse(marker.color, x, y, width, height, false);
        break;
    case litehtml::list_style_type_square:
        fill_rectangle(marker.color, x, y, width, height);
        break;
    default:
        break;
    }
}

void HtmlWidget::load_image(const char* source, const char* base, bool)
{
    image_for(resolve_url(source, base));
}

void HtmlWidget::get_image_size(const char* source, const char* base, litehtml::size& size)
{
    if(auto* image = image_for(resolve_url(source, base)); image != nullptr)
    {
        size.width  = image->width;
        size.height = image->height;
    } else
    {
        size.width = size.height = 0;
    }
}

void HtmlWidget::draw_image(litehtml::uint_ptr, const litehtml::background_layer& layer,
                            const std::string& url, const std::string& base)
{
    if(indexing_)
    {
        return;
    }
    auto* image = image_for(resolve_url(url.c_str(), base.c_str()));
    if(image == nullptr)
    {
        return;
    }
    Tk_PhotoHandle photo = Tk_FindPhoto(interp_, image->name.c_str());
    Tk_PhotoImageBlock source{};
    if(photo == nullptr || Tk_PhotoGetImage(photo, &source) == 0 || source.width <= 0 || source.height <= 0)
    {
        return;
    }
    const int target_width  = std::max(to_int(layer.origin_box.width), 1);
    const int target_height = std::max(to_int(layer.origin_box.height), 1);
    const int origin_x      = to_int(layer.origin_box.x);
    const int origin_y      = to_int(layer.origin_box.y);
    const int clip_left     = to_int(layer.clip_box.left());
    const int clip_top      = to_int(layer.clip_box.top());
    const int clip_right    = to_int(layer.clip_box.right());
    const int clip_bottom   = to_int(layer.clip_box.bottom());
    const auto active = active_clip();
    const int visible_left = std::max(clip_left, to_int(active.left()));
    const int visible_top = std::max(clip_top, to_int(active.top()));
    const int visible_right = std::min(clip_right, to_int(active.right()));
    const int visible_bottom = std::min(clip_bottom, to_int(active.bottom()));
    const auto channel = [&source](const unsigned char* pixel, int component, unsigned char fallback) {
        const int offset = source.offset[component];
        return offset >= 0 && offset < source.pixelSize ? pixel[offset] : fallback;
    };
    for(int y = std::max(origin_y, visible_top);
        y < std::min(saturating_add(origin_y, target_height), visible_bottom); ++y)
    {
        const auto source_y_offset =
            (static_cast<long long>(y) - origin_y) * source.height / target_height;
        const int source_y = static_cast<int>(
            std::clamp(source_y_offset, 0LL, static_cast<long long>(source.height - 1)));
        for(int x = std::max(origin_x, visible_left);
            x < std::min(saturating_add(origin_x, target_width), visible_right); ++x)
        {
            const auto source_x_offset =
                (static_cast<long long>(x) - origin_x) * source.width / target_width;
            const int source_x = static_cast<int>(
                std::clamp(source_x_offset, 0LL, static_cast<long long>(source.width - 1)));
            const unsigned char* pixel = source.pixelPtr + source_y * source.pitch + source_x * source.pixelSize;
            const litehtml::web_color color(channel(pixel, Red, 0), channel(pixel, Green, 0),
                                             channel(pixel, Blue, 0),
                                             channel(pixel, Alpha, kColorComponentMaximum));
            blend_pixel(x, y, color);
        }
    }
}

void HtmlWidget::draw_solid_fill(litehtml::uint_ptr, const litehtml::background_layer& layer,
                                 const litehtml::web_color& color)
{
    if(indexing_ || color.alpha == 0)
    {
        return;
    }
    fill_rectangle(color, to_int(layer.border_box.x), to_int(layer.border_box.y), to_int(layer.border_box.width),
                   to_int(layer.border_box.height));
}

void HtmlWidget::draw_linear_gradient(litehtml::uint_ptr, const litehtml::background_layer& layer,
                                      const litehtml::background_layer::linear_gradient& gradient)
{
    if(indexing_ || gradient.color_points.empty())
    {
        return;
    }
    const auto& first = gradient.color_points.front().color;
    const auto& last  = gradient.color_points.back().color;
    const bool horizontal = std::abs(gradient.end.x - gradient.start.x) > std::abs(gradient.end.y - gradient.start.y);
    const int length = std::max(horizontal ? to_int(layer.border_box.width) : to_int(layer.border_box.height), 1);
    const int origin = horizontal ? to_int(layer.border_box.x) : to_int(layer.border_box.y);
    const auto clip = active_clip();
    const int clip_start = horizontal ? to_int(clip.left()) : to_int(clip.top());
    const int clip_end = horizontal ? to_int(clip.right()) : to_int(clip.bottom());
    const auto first_distance = static_cast<long long>(clip_start) - origin;
    const auto last_distance = static_cast<long long>(clip_end) - origin;
    const int first_offset = static_cast<int>(std::clamp(first_distance, 0LL, static_cast<long long>(length)));
    const int last_offset = static_cast<int>(std::clamp(last_distance, 0LL, static_cast<long long>(length)));
    for(int offset = first_offset; offset < last_offset; ++offset)
    {
        const double fraction = length == 1 ? 0.0 : static_cast<double>(offset) / (length - 1);
        const auto mix = [fraction](int from, int to) {
            return static_cast<litehtml::byte>(std::lround(from + (to - from) * fraction));
        };
        const litehtml::web_color color(mix(first.red, last.red), mix(first.green, last.green),
                                        mix(first.blue, last.blue), mix(first.alpha, last.alpha));
        const int x = saturating_add(to_int(layer.border_box.x), horizontal ? offset : 0);
        const int y = saturating_add(to_int(layer.border_box.y), horizontal ? 0 : offset);
        fill_rectangle(color, x, y, horizontal ? 1 : to_int(layer.border_box.width),
                       horizontal ? to_int(layer.border_box.height) : 1);
    }
}

void HtmlWidget::draw_radial_gradient(litehtml::uint_ptr drawable, const litehtml::background_layer& layer,
                                      const litehtml::background_layer::radial_gradient& gradient)
{
    if(!gradient.color_points.empty())
    {
        draw_solid_fill(drawable, layer, gradient.color_points.back().color);
    }
}

void HtmlWidget::draw_conic_gradient(litehtml::uint_ptr drawable, const litehtml::background_layer& layer,
                                     const litehtml::background_layer::conic_gradient& gradient)
{
    if(!gradient.color_points.empty())
    {
        draw_solid_fill(drawable, layer, gradient.color_points.back().color);
    }
}

void HtmlWidget::draw_borders(litehtml::uint_ptr, const litehtml::borders& borders,
                              const litehtml::position& position, bool)
{
    if(indexing_)
    {
        return;
    }
    const int x      = to_int(position.x);
    const int y      = to_int(position.y);
    const int width  = to_int(position.width);
    const int height = to_int(position.height);
    const auto draw_border = [&](const litehtml::border& border, int bx, int by, int bw, int bh) {
        if(static_cast<float>(border.width) > 0.0F && border.style != litehtml::border_style_none &&
           border.style != litehtml::border_style_hidden)
        {
            fill_rectangle(border.color, bx, by, bw, bh);
        }
    };
    const int left   = to_int(borders.left.width);
    const int right  = to_int(borders.right.width);
    const int top    = to_int(borders.top.width);
    const int bottom = to_int(borders.bottom.width);
    draw_border(borders.top, x, y, width, top);
    draw_border(borders.right, saturating_add(x, width - right), y, right, height);
    draw_border(borders.bottom, x, saturating_add(y, height - bottom), width, bottom);
    draw_border(borders.left, x, y, left, height);
}

void HtmlWidget::set_caption(const char* caption)
{
    title_ = caption != nullptr ? caption : "";
}

void HtmlWidget::set_base_url(const char* url)
{
    document_base_url_ = url != nullptr ? resolve_url(url, base_url_.c_str()) : base_url_;
}

void HtmlWidget::link(const litehtml::document::ptr&, const litehtml::element::ptr&) {}

void HtmlWidget::on_anchor_click(const char* url, const litehtml::element::ptr& element)
{
    if(url == nullptr || suppress_link_ || selection_.first != selection_.last)
    {
        return;
    }
    const std::string target = resolve_url(url, document_base_url_.c_str());
    if(url[0] == '#')
    {
        show_fragment(std::string_view(url + 1));
    }
    invoke_link(target, element);
}

void HtmlWidget::on_mouse_event(const litehtml::element::ptr&, litehtml::mouse_event) {}

void HtmlWidget::free_cursor()
{
    if(cursor_ != None)
    {
        Tk_FreeCursor(display_, cursor_);
        cursor_ = None;
    }
}

void HtmlWidget::set_cursor(const char* cursor)
{
    if(window_ == nullptr)
    {
        return;
    }
    std::string name = cursor != nullptr ? cursor : "";
    if(name == "pointer")
    {
        name = "hand2";
    } else if(name == "text")
    {
        name = "xterm";
    } else if(name == "auto" || name == "default")
    {
        name.clear();
    }
    free_cursor();
    if(name.empty())
    {
        Tk_UndefineCursor(window_);
        return;
    }
    cursor_ = Tk_GetCursor(interp_, window_, Tk_GetUid(name.c_str()));
    if(cursor_ != None)
    {
        Tk_DefineCursor(window_, cursor_);
    } else
    {
        Tcl_ResetResult(interp_);
    }
}

void HtmlWidget::transform_text(std::string& text, litehtml::text_transform transform)
{
    if(transform == litehtml::text_transform_none)
    {
        return;
    }
    const auto input = unicode_chars(text);
    std::string output;
    output.reserve(text.size());
    bool capitalize = transform == litehtml::text_transform_capitalize;
    for(Tcl_UniChar character : input)
    {
        if(transform == litehtml::text_transform_uppercase || capitalize)
        {
            character = Tcl_UniCharToUpper(character);
        } else if(transform == litehtml::text_transform_lowercase)
        {
            character = Tcl_UniCharToLower(character);
        }
        char buffer[TCL_UTF_MAX];
        const int length = Tcl_UniCharToUtf(character, buffer);
        output.append(buffer, static_cast<std::size_t>(length));
        capitalize = false;
    }
    text = std::move(output);
}

void HtmlWidget::import_css(std::string& text, const std::string& url, std::string& base)
{
    const std::string resolved = resolve_url(url.c_str(), base.c_str());
    text = resource("stylesheet", resolved);
    if(text.size() > max_document_bytes_ - document_resource_bytes_)
    {
        throw std::runtime_error("stylesheet exceeds the document byte limit");
    }
    if(text.empty())
    {
        const auto path = path_from_utf8(resolved);
        std::error_code filesystem_error;
        if(!std::filesystem::is_regular_file(path, filesystem_error))
        {
            return;
        }
        const auto size = std::filesystem::file_size(path, filesystem_error);
        if(!filesystem_error && size > max_document_bytes_ - document_resource_bytes_)
        {
            throw std::runtime_error("stylesheet exceeds the document byte limit");
        }
        std::ifstream stream(path, std::ios::binary);
        if(stream)
        {
            std::ostringstream contents;
            contents << stream.rdbuf();
            text = contents.str();
            if(text.size() > max_document_bytes_ - document_resource_bytes_)
            {
                throw std::runtime_error("stylesheet exceeds the document byte limit");
            }
        }
    }
    if(!text.empty())
    {
        document_resource_bytes_ += text.size();
        base = resolved;
    }
}

void HtmlWidget::set_clip(const litehtml::position& position, const litehtml::border_radiuses&)
{
    if(indexing_)
    {
        return;
    }
    clip_stack_.push_back(position);
}

void HtmlWidget::del_clip()
{
    if(indexing_)
    {
        return;
    }
    if(!clip_stack_.empty())
    {
        clip_stack_.pop_back();
    }
}

void HtmlWidget::get_viewport(litehtml::position& viewport) const
{
    const int width  = window_ != nullptr ? Tk_Width(window_) : requested_width_;
    const int height = window_ != nullptr ? Tk_Height(window_) : requested_height_;
    viewport = litehtml::position(scroll_x_, scroll_y_, width, height);
}

litehtml::element::ptr HtmlWidget::create_element(const char*, const litehtml::string_map&,
                                                  const litehtml::document::ptr&)
{
    return nullptr;
}

void HtmlWidget::get_media_features(litehtml::media_features& media) const
{
    litehtml::position viewport;
    get_viewport(viewport);
    media.type          = litehtml::media_type_screen;
    media.width         = viewport.width;
    media.height        = viewport.height;
    media.device_width  = viewport.width;
    media.device_height = viewport.height;
    media.color         = kMediaColorBits;
    media.monochrome    = 0;
    media.color_index   = kMediaPaletteSize;
    media.resolution    = kCssPixelsPerInch;
}

void HtmlWidget::get_language(std::string& language, std::string& culture) const
{
    language = "en";
    culture.clear();
}

int widget_factory(ClientData, Tcl_Interp* interp, int objc, Tcl_Obj* const objv[])
{
    try
    {
        return HtmlWidget::create(interp, objc, objv);
    } catch(...)
    {
        return exception_result(interp, "widget creation");
    }
}
} // namespace

int register_widget_command(Tcl_Interp* interp)
{
    HtmlWidget::register_frame_image();
    if(Tcl_CreateObjCommand(interp, "::tklitehtml::widget", widget_factory, nullptr, nullptr) == nullptr)
    {
        return TCL_ERROR;
    }
    return TCL_OK;
}
} // namespace tklitehtml
