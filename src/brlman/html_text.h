/*                     H T M L _ T E X T . H
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

#ifndef BRLMAN_HTML_TEXT_H
#define BRLMAN_HTML_TEXT_H

#include <string>
#include <vector>

struct GumboInternalOutput;

namespace brlman {

constexpr size_t default_width = 80;

enum text_style {
    normal = 0,
    bold = 1,
    underline = 2
};

struct text_run {
    std::string text;
    unsigned style;
};

struct text_line {
    std::vector<text_run> runs;
    size_t block;
    bool literal = false;
};

/* Keep the input alive: Gumbo's tree contains views into its source. */
class html_document {
public:
    explicit html_document(const std::string &html);
    ~html_document();
    html_document(const html_document &) = delete;
    html_document &operator=(const html_document &) = delete;
    std::vector<text_line> format(size_t width) const;

private:
    std::string source;
    GumboInternalOutput *tree;
};

/* Column-based clipping preserves UTF-8 and never emits document controls. */
std::string line_text(const text_line &line, size_t left = 0,
	    size_t width = std::string::npos, bool styled = false);
size_t text_width(const std::string &text);

bool console_interactive();
bool console_dump(const std::vector<text_line> &lines);
bool console_view(const html_document &document, size_t width);

}

#endif

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
