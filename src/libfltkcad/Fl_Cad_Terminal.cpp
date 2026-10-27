/* BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this file; see the file named COPYING for more
 * information.
 */

#include "common.h"

#include "fltkcad/Fl_Cad_Terminal.H"

#include "fltkcad/CommandHistory.H"

#include <FL/Fl.H>
#include <FL/Fl_Hold_Browser.H>
#include <FL/Fl_Text_Display.H>
#include <FL/Fl_Text_Buffer.H>
#include <FL/fl_draw.H>
#include <FL/fl_utf8.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <utility>

namespace {

constexpr int kCompletionRows = 6;
constexpr int kCompletionRowHeight = 22;
constexpr char kTruncationNotice[] = "\n[output truncated by terminal limit]\n";
// Fl_Text_Display reserves the printable style codes 'A' through '~'.
constexpr std::size_t kStyleCount = '~' - 'A' + 1;
constexpr std::size_t kEscapeLimit = 64;
constexpr int kMaximumEscapeParameter = std::numeric_limits<unsigned short>::max();
constexpr int kAnsiDefault = -1;

struct ConsoleStyle
{
  int foreground = kAnsiDefault;
  int background = kAnsiDefault;
  Fl_Font font = FL_COURIER;
  unsigned attributes = 0;

  bool operator==(const ConsoleStyle& other) const
  {
    return foreground == other.foreground && background == other.background &&
      font == other.font && attributes == other.attributes;
  }
};

Fl_Color rgb_color(int rgb)
{
  return fl_rgb_color(static_cast<unsigned char>(rgb >> 16),
                      static_cast<unsigned char>(rgb >> 8),
                      static_cast<unsigned char>(rgb));
}

int ansi_color(int index)
{
  // The xterm palette consists of 16 base colors, a 6-cube, and 24 grays.
  constexpr std::array<int, 16> base{{
    0x000000, 0xcd0000, 0x00cd00, 0xcdcd00,
    0x0000ee, 0xcd00cd, 0x00cdcd, 0xe5e5e5,
    0x7f7f7f, 0xff0000, 0x00ff00, 0xffff00,
    0x5c5cff, 0xff00ff, 0x00ffff, 0xffffff}};
  if (index < 16) return base[static_cast<std::size_t>(index)];
  if (index >= 232) {
    const int gray = 8 + (index - 232) * 10;
    return (gray << 16) | (gray << 8) | gray;
  }
  const int cube = index - 16;
  const auto level = [](int value) { return value == 0 ? 0 : 55 + value * 40; };
  return (level(cube / 36) << 16) | (level((cube / 6) % 6) << 8) | level(cube % 6);
}


std::size_t previous_character(const std::string& text, std::size_t position)
{
  if (position > 0) {
    --position;
    while (position > 0 && (static_cast<unsigned char>(text[position]) & 0xc0) == 0x80) {
      --position;
    }
  }
  return position;
}

std::size_t next_character(const std::string& text, std::size_t position)
{
  if (position < text.size()) {
    ++position;
    while (position < text.size() &&
           (static_cast<unsigned char>(text[position]) & 0xc0) == 0x80) {
      ++position;
    }
  }
  return position;
}

void bound_utf8(std::string& text, std::size_t limit)
{
  if (text.size() <= limit) {
    return;
  }
  while (limit > 0 && (static_cast<unsigned char>(text[limit]) & 0xc0) == 0x80) {
    --limit;
  }
  text.resize(limit);
}

class CommandConsole : public Fl_Text_Display
{
public:
  using Fl_Text_Display::Fl_Text_Display;

  int columns() const
  {
    if (display_needs_recalc_) const_cast<CommandConsole*>(this)->recalc_display();
    fl_font(textfont(), textsize());
    const int cell = std::max(1, static_cast<int>(fl_width("M")));
    return std::max(1, text_area.w / cell);
  }

  void wrapping(bool enabled)
  {
    wrap_mode(enabled ? WRAP_AT_BOUNDS : WRAP_NONE, 0);
  }

  int first_visible_position() const { return mFirstChar; }
  bool cursor_visible() const
  {
    int x = 0, y = 0;
    return position_to_xy(insert_position(), &x, &y) != 0;
  }

  std::function<int(int)> input_handler;

  int handle(int event) override
  {
    if (input_handler && input_handler(event)) return 1;
    return Fl_Text_Display::handle(event);
  }
};

std::string
common_prefix(const std::vector<std::string>& strings)
{
  if (strings.empty()) {
    return {};
  }

  std::string prefix = strings.front();
  for (auto candidate = strings.begin() + 1; candidate != strings.end();
       ++candidate) {
    const auto mismatch = std::mismatch(prefix.begin(), prefix.end(),
                                        candidate->begin(), candidate->end());
    prefix.erase(mismatch.first, prefix.end());
    if (prefix.empty()) {
      break;
    }
  }
  return prefix;
}

} // namespace

class Fl_Cad_Terminal::Impl
{
public:
  explicit Impl(Fl_Cad_Terminal& owner) : owner_(owner)
  {
    owner_.begin();
    terminal_ = new CommandConsole(owner_.x(), owner_.y(), owner_.w(), owner_.h());
    terminal_->buffer(buffer_);
    terminal_->textfont(FL_COURIER);
    terminal_->textcolor(FL_BLACK);
    terminal_->color(FL_WHITE);
    terminal_->cursor_color(FL_BLACK);
    terminal_->wrapping(wrap_history_);
    terminal_->show_cursor();
    styles_used_.push_back(ConsoleStyle{});
    update_palette();
    terminal_->highlight_data(&styles_, style_table_.data(),
                               static_cast<int>(style_table_.size()), 0, nullptr, nullptr);
    buffer_.add_modify_callback(buffer_changed, this);
    terminal_->input_handler = [this](int event) { return handle(event); };
    completion_ = new Fl_Hold_Browser(owner_.x(), owner_.y(), owner_.w(), popup_height());
    completion_->hide();
    completion_->callback(completion_selected, this);
    owner_.end();
    set_prompt("mged> ");
    layout(owner_.x(), owner_.y(), owner_.w(), owner_.h());
  }

  void layout(int x, int y, int width, int height)
  {
    terminal_->resize(x, y, width, height);
    completion_->resize(x, y + std::max(0, height - popup_height()), width,
                        std::min(height, popup_height()));
  }

  ~Impl()
  {
    terminal_->input_handler = nullptr;
    terminal_->highlight_data(nullptr, nullptr, 0, 0, nullptr, nullptr);
    terminal_->buffer(nullptr);
    buffer_.remove_modify_callback(buffer_changed, this);
  }

  void erase_editor()
  {
    if (dispatching_) return;
    buffer_.remove(output_end_, buffer_.length());
  }

  void render_editor(bool follow_cursor = true)
  {
    if (dispatching_) return;
    std::string line;
    if (output_end_ > 0 && buffer_.byte_at(output_end_ - 1) != '\n') line += '\n';
    line += prompt_;
    input_start_ = output_end_ + static_cast<int>(line.size());
    line += input_;
    const char saved_style = style_byte_;
    style_byte_ = 'A';
    char* existing = buffer_.text_range(output_end_, buffer_.length());
    if (line != existing) buffer_.replace(output_end_, buffer_.length(), line.c_str());
    std::free(existing);
    style_byte_ = saved_style;
    terminal_->insert_position(input_start_ + static_cast<int>(cursor_));
    if (follow_cursor) terminal_->show_insert_position();
  }

  void trim_output()
  {
    const int lines = buffer_.count_lines(0, output_end_);
    int remove = lines > limits_.scrollback_lines ?
      buffer_.skip_lines(0, lines - limits_.scrollback_lines) : 0;
    const int excess = output_end_ - static_cast<int>(limits_.max_scrollback_bytes);
    if (excess > remove) {
      remove = buffer_.utf8_align(excess);
      if (remove < excess) remove = buffer_.next_char(remove);
    }
    if (remove > 0) {
      buffer_.remove(0, remove);
      output_end_ -= remove;
      stream_position_ = std::max(0, stream_position_ - remove);
    }
  }

  void set_prompt(std::string text)
  {
    if (text.find('\0') != std::string::npos) {
      report_error("Terminal prompt contains a NUL byte");
      return;
    }
    if (text.size() > limits_.max_prompt_bytes) {
      bound_utf8(text, limits_.max_prompt_bytes);
      report_error("Terminal prompt exceeded byte limit");
    }
    prompt_ = std::move(text);
    render_editor();
  }

  void set_input(std::string text, bool reset_history = true)
  {
    if (text.find('\0') != std::string::npos) {
      report_error("Terminal input contains a NUL byte");
      return;
    }
    if (text.size() > limits_.max_input_bytes) {
      bound_utf8(text, limits_.max_input_bytes);
      report_error("Terminal input exceeded byte limit");
    }
    input_ = std::move(text);
    cursor_ = input_.size();
    if (reset_history) history_.reset_navigation();
    render_editor();
  }

  bool erase_input_selection()
  {
    int start = 0, end = 0;
    if (!buffer_.selection_position(&start, &end) || start < input_start_) return false;
    cursor_ = static_cast<std::size_t>(start - input_start_);
    input_.erase(cursor_, static_cast<std::size_t>(end - start));
    buffer_.unselect();
    return true;
  }

  void insert(std::string text)
  {
    text.erase(std::remove_if(text.begin(), text.end(), [](unsigned char c) {
      return c < ' ' && c != '\n' && c != '\t' && c != '\r';
    }), text.end());
    std::replace(text.begin(), text.end(), '\r', '\n');
    erase_input_selection();
    const std::size_t available = limits_.max_input_bytes - input_.size();
    if (text.size() > available) {
      bound_utf8(text, available);
      report_error("Terminal input exceeded byte limit");
    }
    input_.insert(cursor_, text);
    cursor_ += text.size();
    history_.reset_navigation();
    hide_completion();
    render_editor();
  }

  int handle(int event)
  {
    if (event == FL_KEYDOWN || event == FL_PASTE) {
      const int position = terminal_->insert_position();
      if (position >= input_start_ && position <= buffer_.length())
        cursor_ = static_cast<std::size_t>(position - input_start_);
    }
    if (event == FL_PASTE) {
      insert(std::string(Fl::event_text(), static_cast<std::size_t>(Fl::event_length())));
      return 1;
    }
    if (event == FL_RELEASE && Fl::event_button() == FL_MIDDLE_MOUSE) {
      Fl::paste(*terminal_, 0);
      return 1;
    }
    if (event != FL_KEYDOWN) return 0;
    const int key = Fl::event_key();
    const int state = Fl::event_state();
    if ((state & FL_SHIFT) && (key == FL_Page_Up || key == FL_Page_Down)) return 0;
    if ((state & FL_CTRL) && (state & FL_SHIFT) && key == 'c') {
      return terminal_->Fl_Text_Display::handle(event);
    }
    if (((state & FL_CTRL) && key == 'v') ||
        ((state & FL_SHIFT) && key == FL_Insert)) {
      Fl::paste(*terminal_, 1);
      return 1;
    }
    if (state & FL_CTRL) {
      switch (key) {
        case 'a': move_cursor(0); return 1;
        case 'e': move_cursor(input_.size()); return 1;
        case 'b': move_cursor(previous_character(input_, cursor_)); return 1;
        case 'f': move_cursor(next_character(input_, cursor_)); return 1;
        case 'p': recall(history_.previous(input_)); return 1;
        case 'n': recall(history_.next()); return 1;
        case 'x': {
          char* selection = buffer_.selection_text();
          if (*selection) Fl::copy(selection, static_cast<int>(std::char_traits<char>::length(selection)), 1);
          std::free(selection);
          if (erase_input_selection()) render_editor();
          return 1;
        }
        case 'u': input_.erase(0, cursor_); cursor_ = 0; break;
        case 'k': input_.erase(cursor_); break;
        case 'w': {
          const std::size_t end = cursor_;
          while (cursor_ > 0 && std::isspace(static_cast<unsigned char>(input_[cursor_ - 1])))
            cursor_ = previous_character(input_, cursor_);
          while (cursor_ > 0 && !std::isspace(static_cast<unsigned char>(input_[cursor_ - 1])))
            cursor_ = previous_character(input_, cursor_);
          input_.erase(cursor_, end - cursor_);
          break;
        }
        case 'c':
          input_.clear(); cursor_ = 0;
          append_output("^C\n");
          if (interrupt_handler_) {
            try { interrupt_handler_(); }
            catch (const std::exception& error) {
              report_error(std::string("Terminal interrupt callback failed: ") + error.what());
            } catch (...) { report_error("Terminal interrupt callback failed"); }
          }
          break;
        case 'l': clear(); return 1;
        case 'd': input_.erase(cursor_, next_character(input_, cursor_) - cursor_); break;
        default: return 0;
      }
      history_.reset_navigation();
      hide_completion();
      render_editor();
      return 1;
    }
    if (handle_key(key)) return 1;
    if (state & (FL_ALT | FL_META)) return 0;
    int deleted = 0;
    if (Fl::compose(deleted)) {
      while (deleted-- > 0 && cursor_ > 0) {
        const std::size_t previous = previous_character(input_, cursor_);
        input_.erase(previous, cursor_ - previous);
        cursor_ = previous;
      }
      insert(std::string(Fl::event_text(), static_cast<std::size_t>(Fl::event_length())));
      return 1;
    }
    return 0;
  }

  void move_cursor(std::size_t position)
  {
    hide_completion();
    int anchor = input_start_ + static_cast<int>(cursor_);
    int start = 0, end = 0;
    if (buffer_.selection_position(&start, &end) && start >= input_start_)
      anchor = anchor == start ? end : start;
    cursor_ = position;
    render_editor();
    if (Fl::event_state() & FL_SHIFT) {
      const int insertion = input_start_ + static_cast<int>(cursor_);
      buffer_.select(std::min(anchor, insertion), std::max(anchor, insertion));
    } else buffer_.unselect();
  }

  int handle_key(int key)
  {
    if (completion_->Fl_Widget::visible()) {
      if (key == FL_Escape) {
        hide_completion();
        return 1;
      }
      if (key == FL_Up || key == FL_Down) {
        move_completion(key == FL_Up ? -1 : 1);
        return 1;
      }
      if (key == FL_Enter || key == FL_KP_Enter) {
        apply_selected_completion();
        return 1;
      }
    }

    switch (key) {
      case FL_Enter:
      case FL_KP_Enter:
        submit();
        return 1;
      case FL_Up:
        recall(history_.previous(input_));
        return 1;
      case FL_Down:
        recall(history_.next());
        return 1;
      case FL_Tab:
        complete();
        return 1;
      case FL_Left: move_cursor(previous_character(input_, cursor_)); return 1;
      case FL_Right: move_cursor(next_character(input_, cursor_)); return 1;
      case FL_Home: move_cursor(0); return 1;
      case FL_End: move_cursor(input_.size()); return 1;
      case FL_BackSpace: {
        if (erase_input_selection()) break;
        const std::size_t previous = previous_character(input_, cursor_);
        input_.erase(previous, cursor_ - previous);
        cursor_ = previous;
        history_.reset_navigation();
        break;
      }
      case FL_Delete:
        if (erase_input_selection()) break;
        input_.erase(cursor_, next_character(input_, cursor_) - cursor_);
        history_.reset_navigation();
        break;
      default:
        hide_completion();
        return 0;
    }
    hide_completion();
    render_editor();
    return 1;
  }

  static void buffer_changed(int position, int inserted, int deleted,
                              int, const char*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    std::string styles(static_cast<std::size_t>(inserted), self.style_byte_);
    for (int offset = 0; offset < inserted; ++offset) {
      const unsigned char byte = self.buffer_.byte_at(position + offset);
      if ((byte & 0xc0) == 0x80) styles[static_cast<std::size_t>(offset)] = '\0';
    }
    self.styles_.replace(position, position + deleted, styles.data(), inserted);
  }

  void update_palette()
  {
    for (std::size_t index = 0; index < style_table_.size(); ++index) {
      const ConsoleStyle& style = styles_used_[index < styles_used_.size() ? index : 0];
      style_table_[index] = {
        style.foreground == kAnsiDefault ? terminal_->textcolor() : rgb_color(style.foreground),
        style.font, terminal_->textsize(),
        style.attributes | Fl_Text_Display::ATTR_BGCOLOR,
        style.background == kAnsiDefault ? terminal_->color() : rgb_color(style.background)};
    }
    terminal_->redisplay_range(0, buffer_.length());
  }

  void select_style()
  {
    auto found = std::find(styles_used_.begin(), styles_used_.end(), stream_style_);
    if (found == styles_used_.end()) {
      if (styles_used_.size() == style_table_.size()) {
        // Retain existing styles; excessive palette changes fall back to defaults.
        style_byte_ = 'A';
        return;
      }
      styles_used_.push_back(stream_style_);
      update_palette();
      found = styles_used_.end() - 1;
    }
    style_byte_ = static_cast<char>('A' + (found - styles_used_.begin()));
  }

  void rendition(const std::vector<int>& codes)
  {
    for (std::size_t index = 0; index < codes.size(); ++index) {
      const int code = codes[index];
      if (code == 0) stream_style_ = ConsoleStyle{};
      else if (code == 1) stream_style_.font |= FL_BOLD;
      else if (code == 3) stream_style_.font |= FL_ITALIC;
      else if (code == 22) stream_style_.font &= ~FL_BOLD;
      else if (code == 23) stream_style_.font &= ~FL_ITALIC;
      else if (code == 4) stream_style_.attributes |= Fl_Text_Display::ATTR_UNDERLINE;
      else if (code == 24) stream_style_.attributes &= ~Fl_Text_Display::ATTR_UNDERLINE;
      else if (code == 9) stream_style_.attributes |= Fl_Text_Display::ATTR_STRIKE_THROUGH;
      else if (code == 29) stream_style_.attributes &= ~Fl_Text_Display::ATTR_STRIKE_THROUGH;
      else if (code == 39) stream_style_.foreground = kAnsiDefault;
      else if (code == 49) stream_style_.background = kAnsiDefault;
      else if (code >= 30 && code <= 37) stream_style_.foreground = ansi_color(code - 30);
      else if (code >= 40 && code <= 47) stream_style_.background = ansi_color(code - 40);
      else if (code >= 90 && code <= 97) stream_style_.foreground = ansi_color(code - 90 + 8);
      else if (code >= 100 && code <= 107) stream_style_.background = ansi_color(code - 100 + 8);
      else if ((code == 38 || code == 48) && index + 1 < codes.size()) {
        int color = kAnsiDefault;
        if (codes[index + 1] == 5 && index + 2 < codes.size()) {
          color = ansi_color(std::clamp(codes[index + 2], 0, 255));
          index += 2;
        } else if (codes[index + 1] == 2 && index + 4 < codes.size()) {
          color = (std::clamp(codes[index + 2], 0, 255) << 16) |
            (std::clamp(codes[index + 3], 0, 255) << 8) |
            std::clamp(codes[index + 4], 0, 255);
          index += 4;
        }
        if (code == 38) stream_style_.foreground = color;
        else stream_style_.background = color;
      }
    }
    select_style();
  }

  void finish_escape()
  {
    if (escape_.size() < 3 || escape_[1] != '[') return;
    std::vector<int> codes(1, 0);
    for (std::size_t index = 2; index + 1 < escape_.size(); ++index) {
      const unsigned char byte = escape_[index];
      if (byte == ';') codes.push_back(0);
      else if (byte >= '0' && byte <= '9') codes.back() =
        std::min(kMaximumEscapeParameter, codes.back() * 10 + byte - '0');
      else return;
    }
    if (escape_.back() == 'm') rendition(codes);
    else if (escape_.back() == 'K') {
      const int begin = buffer_.line_start(stream_position_);
      const int end = buffer_.line_end(stream_position_);
      if (codes.front() == 0) buffer_.remove(stream_position_, end);
      else if (codes.front() == 2) {
        buffer_.remove(begin, end);
        stream_position_ = begin;
      }
    }
  }

  void append_stream(const std::string& text)
  {
    std::string run = std::move(pending_utf8_);
    pending_utf8_.clear();
    const auto flush = [this, &run]() {
      if (run.empty()) return;
      int end = stream_position_;
      for (std::size_t byte = 0; byte < run.size() && end < buffer_.length();
           byte = next_character(run, byte)) end = buffer_.next_char(end);
      buffer_.replace(stream_position_, end, run.c_str());
      stream_position_ += static_cast<int>(run.size());
      run.clear();
    };
    for (unsigned char byte : text) {
      if (!escape_.empty()) {
        const bool overflow = escape_.size() >= kEscapeLimit;
        if (!overflow) escape_ += static_cast<char>(byte);
        const bool csi_end = escape_.size() > 2 && escape_[1] == '[' && byte >= '@' && byte <= '~';
        const bool osc_end = escape_[1] == ']' && (byte == '\a' ||
          (byte == '\\' && osc_escape_));
        osc_escape_ = byte == '\x1b';
        if (csi_end || osc_end || (escape_.size() == 2 && byte != '[' && byte != ']')) {
          if (!overflow) finish_escape();
          escape_.clear();
          osc_escape_ = false;
        }
      } else if (byte == '\x1b') {
        flush();
        escape_ += static_cast<char>(byte);
      } else if (byte == '\n' && stream_position_ == buffer_.length()) {
        run += '\n';
      } else if (byte == '\r' || byte == '\b' || byte == '\n') {
        flush();
        if (byte == '\r') stream_position_ = buffer_.line_start(stream_position_);
        else if (byte == '\b') stream_position_ = std::max(
          buffer_.line_start(stream_position_), buffer_.prev_char(stream_position_));
        else {
          stream_position_ = buffer_.line_end(stream_position_);
          buffer_.insert(stream_position_++, "\n");
        }
      } else if (byte >= ' ' || byte == '\t') run += static_cast<char>(byte);
    }
    // A subprocess may split a UTF-8 character across appends.
    std::size_t last = run.size();
    while (last > 0 && (static_cast<unsigned char>(run[last - 1]) & 0xc0) == 0x80) --last;
    if (last > 0) {
      --last;
      const int length = fl_utf8len(static_cast<unsigned char>(run[last]));
      if (length > 0 && run.size() - last < static_cast<std::size_t>(length)) {
        pending_utf8_ = run.substr(last);
        run.erase(last);
      }
    }
    flush();
  }

  void append_output(const std::string& text)
  {
    const bool selected = buffer_.selected();
    int selection_start = 0, selection_end = 0;
    buffer_.selection_position(&selection_start, &selection_end);
    const bool historical_selection = selected && selection_end <= output_end_;
    const bool input_selection = selected && selection_start >= input_start_;
    const int input_selection_start = selection_start - input_start_;
    const int input_selection_end = selection_end - input_start_;
    const int insertion = terminal_->insert_position();
    if (!dispatching_ && insertion >= input_start_ && insertion <= buffer_.length())
      cursor_ = static_cast<std::size_t>(insertion - input_start_);
    const int top = terminal_->first_visible_position();
    const int horizontal = terminal_->scroll_col();
    const bool reading_history = historical_selection || insertion < input_start_ ||
      !terminal_->cursor_visible();
    erase_editor();
    std::string accepted = text;
    bound_utf8(accepted, limits_.max_append_bytes);
    append_stream(accepted);
    if (accepted.size() != text.size()) {
      escape_.clear();
      pending_utf8_.clear();
      osc_escape_ = false;
      buffer_.append(kTruncationNotice);
      stream_position_ = buffer_.length();
      report_error("Terminal append exceeded byte limit");
    }
    output_end_ = buffer_.length();
    const int before_trim = output_end_;
    trim_output();
    const int removed = before_trim - output_end_;
    render_editor(!reading_history);
    if (historical_selection) {
      buffer_.select(std::max(0, selection_start - removed),
                      std::max(0, selection_end - removed));
    } else if (input_selection) {
      buffer_.select(input_start_ + input_selection_start, input_start_ + input_selection_end);
    }
    if (reading_history) terminal_->scroll(
      terminal_->count_lines(0, std::max(0, top - removed), true) + 1, horizontal);
  }

  void clear()
  {
    buffer_.text("");
    output_end_ = stream_position_ = 0;
    escape_.clear();
    pending_utf8_.clear();
    osc_escape_ = false;
    render_editor();
  }

  void submit()
  {
    hide_completion();
    erase_editor();
    if (output_end_ > 0 && buffer_.byte_at(output_end_ - 1) != '\n') buffer_.append("\n");
    const std::string command = input_;
    const char saved_style = style_byte_;
    style_byte_ = 'A';
    buffer_.append((prompt_ + command + "\n").c_str());
    style_byte_ = saved_style;
    output_end_ = stream_position_ = buffer_.length();
    trim_output();
    if (!command.empty()) history_.add(command);
    input_.clear();
    cursor_ = 0;
    dispatching_ = true;
    dispatch_command(command);
    dispatching_ = false;
    render_editor();
  }

  void dispatch_command(const std::string& command)
  {
    if (command_handler_) {
      try {
        command_handler_(command);
      } catch (const std::exception& error) {
        report_error(std::string("Terminal command callback failed: ") +
                     error.what());
      } catch (...) {
        report_error("Terminal command callback failed");
      }
    }
  }

  void complete()
  {
    hide_completion();
    if (!completion_handler_) {
      return;
    }

    const std::string current = input_;
    const auto cursor = cursor_;
    Completion result;
    try {
      result = completion_handler_(current, cursor);
    } catch (const std::exception& error) {
      report_error(std::string("Terminal completion callback failed: ") +
                   error.what());
      return;
    } catch (...) {
      report_error("Terminal completion callback failed");
      return;
    }
    if (result.begin > result.end || result.end > current.size() ||
        result.candidates.empty()) {
      return;
    }

    active_completion_ = std::move(result);
    if (active_completion_.candidates.size() >
        limits_.max_completion_candidates) {
      active_completion_.candidates.resize(limits_.max_completion_candidates);
      report_error("Terminal completion candidate limit reached");
    }
    const std::size_t unchanged_bytes =
      current.size() - (active_completion_.end - active_completion_.begin);
    if (unchanged_bytes > limits_.max_input_bytes) {
      report_error("Terminal input exceeds the configured byte limit");
      return;
    }
    const std::size_t candidate_limit = limits_.max_input_bytes - unchanged_bytes;
    const auto new_end = std::remove_if(
      active_completion_.candidates.begin(), active_completion_.candidates.end(),
      [candidate_limit](const std::string& candidate) {
        return candidate.size() > candidate_limit;
      });
    if (new_end != active_completion_.candidates.end()) {
      active_completion_.candidates.erase(new_end,
                                          active_completion_.candidates.end());
      report_error("Terminal completion exceeded input byte limit");
    }
    if (active_completion_.candidates.empty()) {
      return;
    }
    if (active_completion_.candidates.size() == 1) {
      apply_completion(active_completion_.candidates.front());
      return;
    }

    const std::string prefix = common_prefix(active_completion_.candidates);
    const std::string selected = current.substr(
      active_completion_.begin,
      active_completion_.end - active_completion_.begin);
    if (prefix.size() > selected.size()) {
      apply_completion(prefix);
      active_completion_.end = active_completion_.begin + prefix.size();
    }

    completion_->clear();
    for (const auto& candidate : active_completion_.candidates) {
      completion_->add(candidate.c_str());
    }
    completion_->value(1);
    completion_->show();
    completion_->redraw();
  }

  void apply_completion(const std::string& candidate)
  {
    input_.replace(active_completion_.begin,
                   active_completion_.end - active_completion_.begin, candidate);
    cursor_ = active_completion_.begin + candidate.size();
    history_.reset_navigation();
    render_editor();
  }

  void apply_selected_completion()
  {
    const int selected = completion_->value();
    const char* candidate = selected > 0 ? completion_->text(selected) : nullptr;
    if (candidate) {
      apply_completion(candidate);
    }
    hide_completion();
  }

  static void completion_selected(Fl_Widget*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    if (Fl::event_clicks()) {
      self.apply_selected_completion();
      self.terminal_->take_focus();
    }
  }

  void recall(const std::optional<std::string>& value)
  {
    if (!value) {
      return;
    }
    set_input(*value, false);
  }

  void move_completion(int delta)
  {
    const int count = completion_->size();
    if (count == 0) {
      return;
    }
    const int next = std::clamp(completion_->value() + delta, 1, count);
    completion_->value(next);
  }

  void hide_completion()
  {
    if (completion_->Fl_Widget::visible()) {
      completion_->hide();
      owner_.redraw();
    }
  }

  static constexpr int popup_height()
  {
    return kCompletionRows * kCompletionRowHeight;
  }

  void report_error(std::string message)
  {
    last_error_ = std::move(message);
    if (error_handler_) {
      try {
        error_handler_(last_error_);
      } catch (...) {
        // Error reporting must not turn a recoverable widget failure into an
        // exception escaping through FLTK's event loop.
      }
    }
  }

  Fl_Cad_Terminal& owner_;
  Fl_Text_Buffer buffer_;
  Fl_Text_Buffer styles_;
  std::array<Fl_Text_Display::Style_Table_Entry, kStyleCount> style_table_{};
  std::vector<ConsoleStyle> styles_used_;
  ConsoleStyle stream_style_;
  char style_byte_ = 'A';
  std::string escape_;
  bool osc_escape_ = false;
  std::string pending_utf8_;
  int stream_position_ = 0;
  CommandConsole* terminal_ = nullptr;
  std::string input_;
  std::size_t cursor_ = 0;
  int output_end_ = 0;
  int input_start_ = 0;
  bool wrap_history_ = true;
  bool dispatching_ = false;
  Fl_Hold_Browser* completion_ = nullptr;
  fltkcad::CommandHistory history_;
  std::string prompt_;
  CommandHandler command_handler_;
  CompletionHandler completion_handler_;
  ErrorHandler error_handler_;
  InterruptHandler interrupt_handler_;
  Completion active_completion_;
  Limits limits_;
  std::string last_error_;
};

Fl_Cad_Terminal::Fl_Cad_Terminal(int x, int y, int width, int height,
                                 const char* label)
  : Fl_Group(x, y, width, height, label), impl_(std::make_unique<Impl>(*this))
{
  resizable(&output());
}

Fl_Cad_Terminal::~Fl_Cad_Terminal() = default;

void
Fl_Cad_Terminal::append(const std::string& text)
{
  impl_->last_error_.clear();
  impl_->append_output(text);
}

void
Fl_Cad_Terminal::clear()
{
  impl_->clear();
}

void
Fl_Cad_Terminal::prompt(std::string text)
{
  impl_->set_prompt(std::move(text));
}

const std::string&
Fl_Cad_Terminal::prompt() const noexcept
{
  return impl_->prompt_;
}

void
Fl_Cad_Terminal::input(std::string text)
{
  impl_->last_error_.clear();
  impl_->set_input(std::move(text));
}

std::string
Fl_Cad_Terminal::input() const
{
  return impl_->input_;
}

void
Fl_Cad_Terminal::submit()
{
  impl_->submit();
}

void
Fl_Cad_Terminal::complete()
{
  impl_->complete();
}

void
Fl_Cad_Terminal::command_handler(CommandHandler handler)
{
  impl_->command_handler_ = std::move(handler);
}

void
Fl_Cad_Terminal::completion_handler(CompletionHandler handler)
{
  impl_->completion_handler_ = std::move(handler);
}

void
Fl_Cad_Terminal::interrupt_handler(InterruptHandler handler)
{
  impl_->interrupt_handler_ = std::move(handler);
}

void
Fl_Cad_Terminal::error_handler(ErrorHandler handler)
{
  impl_->error_handler_ = std::move(handler);
}

void
Fl_Cad_Terminal::limits(Limits value)
{
  const std::size_t maximum_buffer =
    static_cast<std::size_t>(std::numeric_limits<int>::max()) - sizeof(kTruncationNotice);
  if (value.command_history == 0 || value.scrollback_lines <= 0 ||
      value.max_input_bytes == 0 || value.max_append_bytes == 0 ||
      value.max_prompt_bytes == 0 || value.max_completion_candidates == 0 ||
      value.max_scrollback_bytes == 0 ||
      value.max_input_bytes >= maximum_buffer ||
      value.max_prompt_bytes >= maximum_buffer - value.max_input_bytes ||
      value.max_scrollback_bytes >= maximum_buffer - value.max_input_bytes - value.max_prompt_bytes ||
      value.max_append_bytes > maximum_buffer - value.max_scrollback_bytes) {
    impl_->report_error("Terminal limits must be positive and supported by FLTK");
    return;
  }
  impl_->limits_ = value;
  impl_->history_.limit(value.command_history);
  impl_->erase_editor();
  impl_->trim_output();
  impl_->render_editor();
  if (impl_->prompt_.size() > value.max_prompt_bytes) {
    impl_->set_prompt(impl_->prompt_);
  }
  if (impl_->input_.size() > value.max_input_bytes) {
    impl_->set_input(impl_->input_);
  }
}

Fl_Cad_Terminal::Limits
Fl_Cad_Terminal::limits() const noexcept
{
  return impl_->limits_;
}

const std::string&
Fl_Cad_Terminal::last_error() const noexcept
{
  return impl_->last_error_;
}

std::size_t
Fl_Cad_Terminal::history_size() const
{
  return impl_->history_.size();
}

std::string
Fl_Cad_Terminal::history_at(std::size_t index) const
{
  return impl_->history_.at(index);
}

bool
Fl_Cad_Terminal::consolidate_history(std::size_t begin, std::size_t end)
{
  return impl_->history_.consolidate(begin, end);
}

bool
Fl_Cad_Terminal::consolidate_history(std::size_t begin, std::size_t end, const std::string& command)
{
  return impl_->history_.consolidate(begin, end, command);
}

Fl_Text_Display&
Fl_Cad_Terminal::output() noexcept
{
  return *impl_->terminal_;
}

const Fl_Text_Display&
Fl_Cad_Terminal::output() const noexcept
{
  return *impl_->terminal_;
}

int
Fl_Cad_Terminal::columns() const noexcept
{
  return impl_->terminal_->columns();
}

void
Fl_Cad_Terminal::colors(Fl_Color foreground, Fl_Color background)
{
  impl_->terminal_->textcolor(foreground);
  impl_->terminal_->cursor_color(foreground);
  impl_->terminal_->color(background);
  impl_->update_palette();
  impl_->terminal_->redraw();
}

std::string
Fl_Cad_Terminal::text() const
{
  char* contents = impl_->buffer_.text();
  std::string result(contents);
  std::free(contents);
  return result;
}

void
Fl_Cad_Terminal::wrap_history(bool enabled)
{
  if (impl_->wrap_history_ == enabled) return;
  impl_->wrap_history_ = enabled;
  impl_->terminal_->wrapping(enabled);
}

bool
Fl_Cad_Terminal::wrap_history() const noexcept
{
  return impl_->wrap_history_;
}

void
Fl_Cad_Terminal::resize(int x, int y, int width, int height)
{
  Fl_Group::resize(x, y, width, height);
  impl_->layout(x, y, width, height);
}

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
