/*                   H T M L _ T E X T . C P P
 * BRL-CAD
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

#include <algorithm>
#include <climits>
#include <cstdlib>
#include <stdexcept>
#include <wchar.h>

#include <gumbo.h>

#ifdef HAVE_WINDOWS_H
#  include "bio.h"
#endif

#include "html_text.h"

namespace {

constexpr size_t list_indent = 2;
constexpr size_t block_indent = 4;
constexpr size_t tab_width = 8;
constexpr size_t table_gap = 3;

struct glyph {
    std::string text;
    size_t columns;
    unsigned style;
};

/* Decode by scalar value, not by byte, for both wrapping and clipping. */
glyph next_glyph(const std::string &text, size_t &pos, unsigned style)
{
    size_t start = pos;
    unsigned char first = (unsigned char)text[pos++];
    unsigned cp = first;
    size_t extra = 0;
    if (first >= 0xc2 && first <= 0xdf) {
	cp = first & 0x1f;
	extra = 1;
    } else if (first >= 0xe0 && first <= 0xef) {
	cp = first & 0x0f;
	extra = 2;
    } else if (first >= 0xf0 && first <= 0xf4) {
	cp = first & 0x07;
	extra = 3;
    } else if (first >= 0x80) {
	return {"?", 1, style};
    }
    for (size_t i = 0; i < extra; i++) {
	if (pos == text.size() || ((unsigned char)text[pos] & 0xc0) != 0x80)
	    return {"?", 1, style};
	cp = (cp << 6) | ((unsigned char)text[pos++] & 0x3f);
    }
    if ((extra == 1 && cp < 0x80) || (extra == 2 && cp < 0x800) ||
	(extra == 3 && cp < 0x10000) || cp > 0x10ffff ||
	(cp >= 0xd800 && cp <= 0xdfff))
	return {"?", 1, style};

    if (cp < 0x20 || (cp >= 0x7f && cp < 0xa0))
	return {"", 0, style};

    int columns = 1;
#ifdef BRLMAN_HAVE_WCWIDTH
    columns = wcwidth((wchar_t)cp);
    if (columns >= 0)
	return {text.substr(start, pos - start), (size_t)columns, style};
    columns = 1;
#endif
    /* A deterministic fallback also covers non-UTF-8 POSIX locales. */
#ifdef HAVE_WINDOWS_H
    if (cp <= 0xffff) {
	WCHAR wc = (WCHAR)cp;
	WORD type = 0;
	if (GetStringTypeW(CT_CTYPE3, &wc, 1, &type) &&
	    (type & C3_NONSPACING))
	    columns = 0;
    }
#endif
    if ((cp >= 0x0300 && cp <= 0x036f) || (cp >= 0x1ab0 && cp <= 0x1aff) ||
	(cp >= 0x1dc0 && cp <= 0x1dff) || (cp >= 0x20d0 && cp <= 0x20ff) ||
	(cp >= 0xfe00 && cp <= 0xfe0f) || (cp >= 0xfe20 && cp <= 0xfe2f) ||
	(cp >= 0xe0100 && cp <= 0xe01ef) || cp == 0x200b || cp == 0x200c ||
	cp == 0x200d || cp == 0xfeff)
	columns = 0;
    /* Conventional wide terminal ranges; ambiguous characters stay narrow. */
    if ((cp >= 0x1100 && cp <= 0x115f) || cp == 0x2329 || cp == 0x232a ||
	(cp >= 0x2e80 && cp <= 0xa4cf && cp != 0x303f) ||
	(cp >= 0xac00 && cp <= 0xd7a3) || (cp >= 0xf900 && cp <= 0xfaff) ||
	(cp >= 0xfe10 && cp <= 0xfe19) || (cp >= 0xfe30 && cp <= 0xfe6f) ||
	(cp >= 0xff00 && cp <= 0xff60) || (cp >= 0xffe0 && cp <= 0xffe6) ||
	(cp >= 0x1f300 && cp <= 0x1faff) || (cp >= 0x20000 && cp <= 0x3fffd))
	columns = 2;
    return {text.substr(start, pos - start), (size_t)columns, style};
}

void add_run(brlman::text_line &line, const std::string &text, unsigned style)
{
    if (text.empty())
	return;
    if (!line.runs.empty() && line.runs.back().style == style)
	line.runs.back().text += text;
    else
	line.runs.push_back({text, style});
}

const GumboVector *children(const GumboNode *node)
{
    if (node->type == GUMBO_NODE_DOCUMENT)
	return &node->v.document.children;
    if (node->type == GUMBO_NODE_ELEMENT)
	return &node->v.element.children;
    return nullptr;
}

std::string attribute(const GumboNode *node, const char *name)
{
    if (node->type != GUMBO_NODE_ELEMENT)
	return {};
    const GumboAttribute *attr = gumbo_get_attribute(&node->v.element.attributes, name);
    return attr ? attr->value : "";
}

bool ignored(const GumboNode *node)
{
    if (node->type != GUMBO_NODE_ELEMENT)
	return node->type == GUMBO_NODE_COMMENT || node->type == GUMBO_NODE_TEMPLATE;
    switch (node->v.element.tag) {
	case GUMBO_TAG_HEAD:
	case GUMBO_TAG_STYLE:
	case GUMBO_TAG_SCRIPT:
	case GUMBO_TAG_NOSCRIPT:
	    return true;
	default:
	    break;
    }
    return gumbo_get_attribute(&node->v.element.attributes, "hidden") ||
	attribute(node, "id") == "footer";
}

bool block_tag(GumboTag tag)
{
    switch (tag) {
	case GUMBO_TAG_DIV: case GUMBO_TAG_SECTION: case GUMBO_TAG_ARTICLE:
	case GUMBO_TAG_HEADER: case GUMBO_TAG_FOOTER: case GUMBO_TAG_MAIN:
	case GUMBO_TAG_P: case GUMBO_TAG_H1: case GUMBO_TAG_H2:
	case GUMBO_TAG_H3: case GUMBO_TAG_H4: case GUMBO_TAG_H5:
	case GUMBO_TAG_H6: case GUMBO_TAG_DT: case GUMBO_TAG_DD:
	case GUMBO_TAG_DL: case GUMBO_TAG_HR:
	    return true;
	default:
	    return false;
    }
}

class formatter {
public:
    explicit formatter(size_t columns) : width(std::max(columns, size_t(1))) {}

    std::vector<brlman::text_line> format(const GumboNode *node)
    {
	walk(node, 0, brlman::normal);
	flush();
	while (!lines.empty() && lines.back().runs.empty())
	    lines.pop_back();
	return lines;
    }

private:
    size_t width;
    size_t indent = 0;
    size_t hanging = 0;
    size_t block = 0;
    bool space = false;
    std::vector<glyph> pending;
    std::vector<brlman::text_line> lines;

    void blank()
    {
	if (!lines.empty() && !lines.back().runs.empty())
	    lines.push_back({{}, block});
    }

    void append(const std::string &text, unsigned style)
    {
	for (size_t pos = 0; pos < text.size();) {
	    unsigned char ch = (unsigned char)text[pos];
	    if (ch == ' ' || ch == '\t' || ch == '\r' || ch == '\n' || ch == '\f') {
		space = !pending.empty() && pending.back().text != " ";
		pos++;
		continue;
	    }
	    glyph g = next_glyph(text, pos, style);
	    if (g.text.empty())
		continue;
	    if (space) {
		pending.push_back({" ", 1, brlman::normal});
		space = false;
	    }
	    pending.push_back(g);
	}
    }

    void flush()
    {
	if (pending.empty()) {
	    hanging = 0;
	    space = false;
	    return;
	}
	size_t start = 0;
	size_t margin = std::min(indent, width - 1);
	while (start < pending.size()) {
	    size_t end = start;
	    size_t used = margin;
	    size_t last_space = start;
	    while (end < pending.size()) {
		if (used + pending[end].columns > width && end > start)
		    break;
		used += pending[end].columns;
		if (pending[end].text == " ")
		    last_space = end;
		end++;
	    }
	    size_t next = end;
	    if (end < pending.size() && last_space > start) {
		end = last_space;
		next = last_space + 1;
	    }
	    while (end > start && pending[end - 1].text == " ")
		end--;
	    brlman::text_line line{{}, block};
	    add_run(line, std::string(margin, ' '), brlman::normal);
	    for (size_t i = start; i < end; i++)
		add_run(line, pending[i].text, pending[i].style);
	    lines.push_back(line);
	    start = next;
	    while (start < pending.size() && pending[start].text == " ")
		start++;
	    margin = std::min(indent + hanging, width - 1);
	}
	pending.clear();
	hanging = 0;
	space = false;
	block++;
    }

    void walk_children(const GumboNode *node, size_t depth, unsigned style)
    {
	const GumboVector *nodes = children(node);
	if (!nodes)
	    return;
	for (size_t i = 0; i < nodes->length; i++)
	    walk((const GumboNode *)nodes->data[i], depth, style);
    }

    void literal_text(const GumboNode *node, std::string &out)
    {
	if (ignored(node))
	    return;
	if (node->type == GUMBO_NODE_TEXT || node->type == GUMBO_NODE_WHITESPACE || node->type == GUMBO_NODE_CDATA) {
	    out += node->v.text.text;
	    return;
	}
	if (node->type == GUMBO_NODE_ELEMENT && node->v.element.tag == GUMBO_TAG_BR) {
	    out += '\n';
	    return;
	}
	const GumboVector *nodes = children(node);
	if (nodes)
	    for (size_t i = 0; i < nodes->length; i++)
		literal_text((const GumboNode *)nodes->data[i], out);
    }

    void literal(const GumboNode *node, size_t depth, unsigned style)
    {
	std::string text;
	literal_text(node, text);
	brlman::text_line line{{}, block, true};
	size_t column = depth;
	add_run(line, std::string(depth, ' '), brlman::normal);
	for (size_t pos = 0; pos < text.size();) {
	    char ch = text[pos];
	    if (ch == '\n') {
		lines.push_back(line);
		line = {{}, block, true};
		add_run(line, std::string(depth, ' '), brlman::normal);
		column = depth;
		pos++;
	    } else if (ch == '\t') {
		size_t count = tab_width - ((column - depth) % tab_width);
		add_run(line, std::string(count, ' '), style);
		column += count;
		pos++;
	    } else {
		glyph g = next_glyph(text, pos, style);
		add_run(line, g.text, style);
		column += g.columns;
	    }
	}
	if (column > depth)
	    lines.push_back(line);
	block++;
    }

    void list(const GumboNode *node, size_t depth, unsigned style)
    {
	const GumboVector *nodes = children(node);
	long number = 1;
	std::string start = attribute(node, "start");
	if (!start.empty()) {
	    char *end = nullptr;
	    long value = strtol(start.c_str(), &end, 10);
	    if (end && !*end && value > 0)
		number = value;
	}
	for (size_t i = 0; nodes && i < nodes->length; i++) {
	    const GumboNode *item = (const GumboNode *)nodes->data[i];
	    if (item->type != GUMBO_NODE_ELEMENT || item->v.element.tag != GUMBO_TAG_LI)
		continue;
	    std::string prefix = node->v.element.tag == GUMBO_TAG_OL ?
		std::to_string(number) + ". " : "* ";
	    if (number < LONG_MAX)
		number++;
	    indent = depth;
	    append(prefix, style);
	    /* Materialize the marker's trailing space before its first block. */
	    if (space) {
		pending.push_back({" ", 1, brlman::normal});
		space = false;
	    }
	    hanging = prefix.size();
	    const GumboVector *body = children(item);
	    for (size_t j = 0; body && j < body->length; j++) {
		const GumboNode *child = (const GumboNode *)body->data[j];
		if (child->type == GUMBO_NODE_ELEMENT && child->v.element.tag == GUMBO_TAG_P) {
		    walk_children(child, depth + prefix.size(), style);
		    flush();
		} else {
		    walk(child, depth + prefix.size(), style);
		}
	    }
	    flush();
	}
    }

    void table_rows(const GumboNode *node, std::vector<const GumboNode *> &rows)
    {
	const GumboVector *nodes = children(node);
	for (size_t i = 0; nodes && i < nodes->length; i++) {
	    const GumboNode *child = (const GumboNode *)nodes->data[i];
	    if (child->type != GUMBO_NODE_ELEMENT)
		continue;
	    if (child->v.element.tag == GUMBO_TAG_TR)
		rows.push_back(child);
	    else if (child->v.element.tag != GUMBO_TAG_TABLE)
		table_rows(child, rows);
	}
    }

    void table(const GumboNode *node, size_t depth, unsigned style)
    {
	using cell = std::vector<brlman::text_line>;
	std::vector<std::vector<cell>> rows;
	std::vector<const GumboNode *> row_nodes;
	std::vector<size_t> widths;
	bool header = false;
	bool spans = false;
	const GumboVector *body = children(node);
	for (size_t i = 0; body && i < body->length; i++) {
	    const GumboNode *child = (const GumboNode *)body->data[i];
	    if (child->type == GUMBO_NODE_ELEMENT && child->v.element.tag == GUMBO_TAG_CAPTION) {
		walk_children(child, depth, style | brlman::bold);
		flush();
	    }
	}
	table_rows(node, row_nodes);
	for (const GumboNode *row_node : row_nodes) {
	    std::vector<cell> row;
	    const GumboVector *nodes = children(row_node);
	    for (size_t i = 0; nodes && i < nodes->length; i++) {
		const GumboNode *entry = (const GumboNode *)nodes->data[i];
		if (entry->type != GUMBO_NODE_ELEMENT ||
		    (entry->v.element.tag != GUMBO_TAG_TD && entry->v.element.tag != GUMBO_TAG_TH))
		    continue;
		if (rows.empty() && entry->v.element.tag == GUMBO_TAG_TH)
		    header = true;
		for (const char *name : {"colspan", "rowspan"}) {
		    std::string span = attribute(entry, name);
		    spans = spans || (!span.empty() && span != "1");
		}
		formatter cell_format(width);
		row.push_back(cell_format.format(entry));
		if (widths.size() < row.size())
		    widths.resize(row.size());
		for (const auto &line : row.back())
		    widths[row.size() - 1] = std::max(widths[row.size() - 1],
			brlman::text_width(brlman::line_text(line)));
	    }
	    if (!row.empty())
		rows.push_back(row);
	}
	size_t total = depth;
	for (size_t w : widths)
	    total += w;
	if (!widths.empty())
	    total += table_gap * (widths.size() - 1);
	if (total <= width && !spans) {
	    for (const auto &row : rows) {
		size_t height = 0;
		for (const auto &entry : row)
		    height = std::max(height, entry.size());
		for (size_t y = 0; y < height; y++) {
		    brlman::text_line line{{}, block};
		    add_run(line, std::string(depth, ' '), brlman::normal);
		    for (size_t x = 0; x < row.size(); x++) {
			size_t used = 0;
			if (y < row[x].size()) {
			    for (const auto &run : row[x][y].runs) {
				add_run(line, run.text, run.style | style);
				used += brlman::text_width(run.text);
			    }
			}
			if (x + 1 < row.size())
			    add_run(line, std::string(widths[x] - used + table_gap, ' '), brlman::normal);
		    }
		    lines.push_back(line);
		}
		block++;
	    }
	} else {
	    /* Labeled cells retain reading order at widths unsuitable for columns.
	     * Spanning tables use positional cells rather than misleading labels. */
	    bool labels = header && !spans && rows.size() > 1;
	    if (labels)
		block++;
	    for (size_t y = labels ? 1 : 0; y < rows.size(); y++) {
		size_t row_block = block;
		for (size_t x = 0; x < rows[y].size(); x++) {
		    indent = depth;
		    if (labels && x < rows[0].size() && !rows[0][x].empty()) {
			for (size_t label_line = 0; label_line < rows[0][x].size(); label_line++) {
			    if (label_line)
				append(" ", brlman::bold);
			    append(brlman::line_text(rows[0][x][label_line]), brlman::bold);
			}
			append(": ", brlman::bold);
		    }
		    for (const auto &line : rows[y][x]) {
			if (line.literal) {
			    flush();
			    brlman::text_line literal_line{{}, block, true};
			    add_run(literal_line, std::string(depth, ' '), brlman::normal);
			    for (const auto &run : line.runs)
				add_run(literal_line, run.text, run.style | style);
			    lines.push_back(literal_line);
			    continue;
			}
			for (const auto &run : line.runs)
			    append(run.text, run.style | style);
			append(" ", brlman::normal);
		    }
		    flush();
		}
		/* A row is one resize anchor in either table presentation. */
		for (auto it = lines.rbegin(); it != lines.rend() && it->block >= row_block; ++it)
		    it->block = row_block;
		block = row_block + 1;
		blank();
	    }
	}
    }

    void walk(const GumboNode *node, size_t depth, unsigned style)
    {
	if (ignored(node))
	    return;
	if (node->type == GUMBO_NODE_TEXT || node->type == GUMBO_NODE_WHITESPACE || node->type == GUMBO_NODE_CDATA) {
	    if (pending.empty())
		indent = depth;
	    append(node->v.text.text, style);
	    return;
	}
	if (node->type != GUMBO_NODE_ELEMENT) {
	    walk_children(node, depth, style);
	    return;
	}
	GumboTag tag = node->v.element.tag;
	if (tag == GUMBO_TAG_B || tag == GUMBO_TAG_STRONG || tag == GUMBO_TAG_TH ||
	    (tag >= GUMBO_TAG_H1 && tag <= GUMBO_TAG_H6) || tag == GUMBO_TAG_DT)
	    style |= brlman::bold;
	if (tag == GUMBO_TAG_I || tag == GUMBO_TAG_EM || tag == GUMBO_TAG_U)
	    style |= brlman::underline;
	if (tag == GUMBO_TAG_BR) {
	    flush();
	    return;
	}
	if (tag == GUMBO_TAG_IMG) {
	    std::string alt = attribute(node, "alt");
	    if (!alt.empty())
		append("[" + alt + "]", style);
	    return;
	}
	bool boundary = block_tag(tag) || tag == GUMBO_TAG_PRE || tag == GUMBO_TAG_TABLE ||
	    tag == GUMBO_TAG_UL || tag == GUMBO_TAG_OL || tag == GUMBO_TAG_BLOCKQUOTE;
	if (boundary) {
	    flush();
	    blank();
	    indent = depth;
	}
	if (tag == GUMBO_TAG_DD || tag == GUMBO_TAG_BLOCKQUOTE)
	    depth += block_indent;
	if (boundary)
	    indent = depth;
	if (tag == GUMBO_TAG_PRE)
	    literal(node, depth, style);
	else if (tag == GUMBO_TAG_UL || tag == GUMBO_TAG_OL)
	    list(node, depth + list_indent, style);
	else if (tag == GUMBO_TAG_TABLE)
	    table(node, depth, style);
	else {
	    size_t first = pending.size();
	    walk_children(node, depth, style);
	    if (tag == GUMBO_TAG_A) {
		std::string href = attribute(node, "href");
		std::string label;
		for (size_t i = first; i < pending.size(); i++)
		    label += pending[i].text;
		label.erase(0, label.find_first_not_of(' '));
		if (!label.empty())
		    label.erase(label.find_last_not_of(' ') + 1);
		if (!href.empty() && href[0] != '#' && href != label &&
		    (href.compare(0, 7, "http://") == 0 || href.compare(0, 8, "https://") == 0 ||
		     href.compare(0, 7, "mailto:") == 0))
		    append(" <" + href + ">", style);
	    }
	}
	if (boundary) {
	    flush();
	    blank();
	}
    }
};

}

namespace brlman {

html_document::html_document(const std::string &html) : source(html),
    tree(gumbo_parse_with_options(&kGumboDefaultOptions, source.data(), source.size()))
{
    if (!tree)
	throw std::runtime_error("Unable to parse HTML manual page");
}

html_document::~html_document()
{
    gumbo_destroy_output(&kGumboDefaultOptions, tree);
}

std::vector<text_line>
html_document::format(size_t width) const
{
    return formatter(width).format(tree->root);
}

size_t
text_width(const std::string &text)
{
    size_t columns = 0;
    for (size_t pos = 0; pos < text.size();)
	columns += next_glyph(text, pos, normal).columns;
    return columns;
}

std::string
line_text(const text_line &line, size_t left, size_t width, bool styled)
{
    std::string out;
    size_t column = 0;
    unsigned active = normal;
    bool base_visible = false;
    for (const auto &run : line.runs) {
	for (size_t pos = 0; pos < run.text.size();) {
	    glyph g = next_glyph(run.text, pos, run.style);
	    if (!g.columns && !base_visible)
		continue;
	    size_t start = column;
	    column += g.columns;
	    if (column <= left && g.columns)
		continue;
	    if (start < left) {
		if (column > left)
		    out.append(std::min(column - left, width), ' ');
		continue;
	    }
	    if (start - left > width || g.columns > width - (start - left))
		return styled && active ? out + "\033[0m" : out;
	    if (styled && active != run.style) {
		out += "\033[0m";
		if (run.style & bold)
		    out += "\033[1m";
		if (run.style & underline)
		    out += "\033[4m";
		active = run.style;
	    }
	    out += g.text;
	    if (g.columns)
		base_visible = true;
	}
    }
    if (styled && active)
	out += "\033[0m";
    return out;
}

}

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
