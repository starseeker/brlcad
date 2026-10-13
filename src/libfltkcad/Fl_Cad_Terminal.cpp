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
#include <FL/Fl_Terminal.H>
#include <cctype>

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <utility>

namespace {

constexpr int kCompletionRows = 6;
constexpr int kCompletionRowHeight = 22;
constexpr int kTabColumns = 8;
constexpr std::size_t kAppendChunkBytes = 64U * 1024U;
constexpr const char* kTruncationNotice = "\n[output truncated by terminal limit]\n";

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

class CommandTerminal : public Fl_Terminal
{
public:
  using Fl_Terminal::Fl_Terminal;
  using Fl_Terminal::cursor_col;
  using Fl_Terminal::cursor_row;
  using Fl_Terminal::clear_eod;
  using Fl_Terminal::scroll;

  int columns() const
  {
    const int pixels = w() - Fl::box_dw(box()) - margin_left() -
      margin_right() - scrollbar_actual_size();
    return std::max(1, w_to_col(std::max(0, pixels)));
  }

  void fit_columns()
  {
    const int available = columns();
    if (available > display_columns()) display_columns(available);
    wrap_columns(available);
    hscrollbar->value(0);
  }

  void colors(Fl_Color foreground, Fl_Color background)
  {
    const Fl_Color old_foreground = textfgcolor_default();
    const Fl_Color old_background = textbgcolor_default();
    if (foreground == old_foreground && background == old_background) return;
    for (int row = 0; row < ring_rows(); ++row) {
      Utf8Char* cells = u8c_ring_row(row);
      for (int col = 0; col < ring_cols(); ++col) {
        Utf8Char& cell = cells[col];
        CharStyle style = current_style();
        const std::string text(cell.text_utf8(), cell.length());
        style.attrib(cell.attrib());
        style.fgcolor(cell.fgcolor() == old_foreground ? foreground : cell.fgcolor());
        style.bgcolor(cell.bgcolor() == old_background ? background : cell.bgcolor());
        style.charflags(cell.charflags());
        cell.text_utf8(text.data(), static_cast<int>(text.size()), style);
      }
    }
    textfgcolor_default(foreground);
    textbgcolor_default(background);
    if (textfgcolor() == old_foreground) textfgcolor(foreground);
    if (textbgcolor() == old_background) textbgcolor(background);
    color(background);
    cursorfgcolor(background);
    cursorbgcolor(foreground);
    redraw();
  }

  std::function<int(int)> input_handler;
  std::function<void()> before_resize;
  std::function<void()> after_resize;

  int handle(int event) override
  {
    if (input_handler && input_handler(event)) {
      return 1;
    }
    return Fl_Terminal::handle(event);
  }

  void resize(int x, int y, int width, int height) override
  {
    if (before_resize) before_resize();
    Fl_Terminal::resize(x, y, width, height);
    fit_columns();
    if (after_resize) after_resize();
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
    terminal_ = new CommandTerminal(owner_.x(), owner_.y(), owner_.w(), owner_.h());
    terminal_->history_lines(limits_.scrollback_lines);
    terminal_->textfont(FL_COURIER);
    terminal_->colors(FL_BLACK, FL_WHITE);
    terminal_->input_handler = [this](int event) { return handle(event); };
    terminal_->before_resize = [this]() { erase_editor(); };
    terminal_->after_resize = [this]() {
      output_row_ = terminal_->cursor_row();
      output_col_ = terminal_->cursor_col();
      render_editor();
    };
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

  void erase_editor()
  {
    if (dispatching_) return;
    terminal_->cursor_row(output_row_);
    terminal_->cursor_col(output_col_);
    terminal_->clear_eod();
  }

  void render_editor()
  {
    if (dispatching_) return;
    erase_editor();
    const std::string line = prompt_ + input_;
    const std::size_t insertion = prompt_.size() + cursor_;
    const int columns = std::max(1, terminal_->columns());
    const int rows = std::max(1, terminal_->display_rows());
    struct Cell { std::size_t byte; std::size_t length; int row; int col; };
    std::vector<Cell> cells;
    int row = 0, col = 0, cursor_row = 0, cursor_col = 0;
    for (std::size_t byte = 0; byte <= line.size();) {
      if (byte == insertion) { cursor_row = row; cursor_col = col; }
      if (byte == line.size()) break;
      const std::size_t next = next_character(line, byte);
      if (line[byte] == '\n') { ++row; col = 0; }
      else if (line[byte] == '\t') {
        const int spaces = kTabColumns - col % kTabColumns;
        row += (col + spaces) / columns;
        col = (col + spaces) % columns;
      } else {
        cells.push_back({byte, next - byte, row, col});
        if (++col == columns) { ++row; col = 0; }
      }
      byte = next;
    }
    const int visible_rows = std::min(rows, row + 1);
    int start = output_row_ + (output_col_ > 0 ? 1 : 0);
    const int overflow = std::max(0, start + visible_rows - rows);
    if (overflow > 0) {
      terminal_->scroll(overflow);
      output_row_ -= overflow;
      start -= overflow;
      if (output_row_ < 0) { output_row_ = 0; output_col_ = 0; }
    }
    start = std::max(0, start);
    // Very long commands keep the insertion point visible without adding
    // unfinished input to the transcript or consuming unbounded screen rows.
    const int first_row = std::clamp(cursor_row - visible_rows + 1, 0,
                                      std::max(0, row - visible_rows + 1));
    for (const Cell& cell : cells) {
      const int screen_row = start + cell.row - first_row;
      if (screen_row >= start && screen_row < rows) {
        terminal_->plot_char(line.data() + cell.byte, static_cast<int>(cell.length),
                              screen_row, cell.col);
      }
    }
    terminal_->cursor_row(start + cursor_row - first_row);
    terminal_->cursor_col(cursor_col);
    terminal_->redraw();
  }

  void set_prompt(std::string text)
  {
    if (text.size() > limits_.max_prompt_bytes) {
      bound_utf8(text, limits_.max_prompt_bytes);
      report_error("Terminal prompt exceeded byte limit");
    }
    prompt_ = std::move(text);
    render_editor();
  }

  void set_input(std::string text, bool reset_history = true)
  {
    if (text.size() > limits_.max_input_bytes) {
      bound_utf8(text, limits_.max_input_bytes);
      report_error("Terminal input exceeded byte limit");
    }
    input_ = std::move(text);
    cursor_ = input_.size();
    if (reset_history) history_.reset_navigation();
    render_editor();
  }

  void insert(std::string text)
  {
    text.erase(std::remove_if(text.begin(), text.end(), [](unsigned char c) {
      return c < ' ' && c != '\n' && c != '\t' && c != '\r';
    }), text.end());
    std::replace(text.begin(), text.end(), '\r', '\n');
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
      // Fl_Terminal owns selection and clipboard copying.
      return terminal_->Fl_Terminal::handle(event);
    }
    if (((state & FL_CTRL) && key == 'v') ||
        ((state & FL_SHIFT) && key == FL_Insert)) {
      Fl::paste(*terminal_, 1);
      return 1;
    }
    if (state & FL_CTRL) {
      switch (key) {
        case 'a': cursor_ = 0; break;
        case 'e': cursor_ = input_.size(); break;
        case 'b': cursor_ = previous_character(input_, cursor_); break;
        case 'f': cursor_ = next_character(input_, cursor_); break;
        case 'p': recall(history_.previous(input_)); return 1;
        case 'n': recall(history_.next()); return 1;
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
      case FL_Left: cursor_ = previous_character(input_, cursor_); break;
      case FL_Right: cursor_ = next_character(input_, cursor_); break;
      case FL_Home: cursor_ = 0; break;
      case FL_End: cursor_ = input_.size(); break;
      case FL_BackSpace: {
        const std::size_t previous = previous_character(input_, cursor_);
        input_.erase(previous, cursor_ - previous);
        cursor_ = previous;
        history_.reset_navigation();
        break;
      }
      case FL_Delete:
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

  void append_output(const std::string& text)
  {
    erase_editor();
    const std::size_t accepted = std::min(text.size(), limits_.max_append_bytes);
    for (std::size_t offset = 0; offset < accepted; offset += kAppendChunkBytes) {
      const std::size_t length = std::min(kAppendChunkBytes, accepted - offset);
      terminal_->append(text.data() + offset, static_cast<int>(length));
    }
    if (accepted != text.size()) {
      terminal_->append(kTruncationNotice);
      report_error("Terminal append exceeded byte limit");
    }
    output_row_ = terminal_->cursor_row();
    output_col_ = terminal_->cursor_col();
    render_editor();
  }

  void clear()
  {
    erase_editor();
    terminal_->clear_screen_home();
    output_row_ = output_col_ = 0;
    render_editor();
  }

  void submit()
  {
    hide_completion();
    erase_editor();
    if (terminal_->cursor_col() > 0) terminal_->append("\n");
    const std::string command = input_;
    const bool ansi = terminal_->ansi();
    terminal_->ansi(false);
    terminal_->append((prompt_ + command + "\n").c_str());
    terminal_->ansi(ansi);
    if (!command.empty()) history_.add(command);
    input_.clear();
    cursor_ = 0;
    dispatching_ = true;
    dispatch_command(command);
    dispatching_ = false;
    output_row_ = terminal_->cursor_row();
    output_col_ = terminal_->cursor_col();
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
  CommandTerminal* terminal_ = nullptr;
  std::string input_;
  std::size_t cursor_ = 0;
  int output_row_ = 0;
  int output_col_ = 0;
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
  if (value.command_history == 0 || value.scrollback_lines <= 0 ||
      value.max_input_bytes == 0 || value.max_append_bytes == 0 ||
      value.max_prompt_bytes == 0 || value.max_completion_candidates == 0 ||
      value.max_input_bytes > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    impl_->report_error("Terminal limits must be positive and supported by FLTK");
    return;
  }
  impl_->limits_ = value;
  impl_->history_.limit(value.command_history);
  impl_->terminal_->history_lines(value.scrollback_lines);
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

Fl_Terminal&
Fl_Cad_Terminal::output() noexcept
{
  return *impl_->terminal_;
}

const Fl_Terminal&
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
  impl_->erase_editor();
  impl_->terminal_->colors(foreground, background);
  impl_->render_editor();
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
