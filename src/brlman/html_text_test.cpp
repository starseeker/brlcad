/*                   H T M L _ T E X T _ T E S T . C P P
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

#include <clocale>
#include <iostream>
#include <string>

#include "bu/app.h"

#include "html_text.h"

namespace {

std::string render(const std::string &html, size_t width)
{
    brlman::html_document document(html);
    std::string text;
    for (const auto &line : document.format(width))
	text += brlman::line_text(line) + "\n";
    return text;
}

bool expect(const std::string &name, const std::string &actual, const std::string &expected)
{
    if (actual == expected)
	return true;
    std::cerr << name << ": expected:\n" << expected << "got:\n" << actual;
    return false;
}

}

int main(int, char **argv)
{
    bu_setprogname(argv[0]);
    std::setlocale(LC_CTYPE, "");
    bool ok = true;
    ok &= expect("inline boundaries and entities",
	render("<head><title>omit</title><style>omit</style></head>"
	    "<h1>Manual</h1><p>one<b>two</b>three &amp; &lt;four&gt;</p>"
	    "<script>omit</script><div hidden>omit</div><div id='footer'>omit</div>", 80),
	"Manual\n\nonetwothree & <four>\n");
    ok &= expect("paragraph reflow", render("<p>one two three four five</p>", 12),
	"one two\nthree four\nfive\n");
    ok &= expect("description indentation",
	render("<dl><dt>--option</dt><dd>one two three four five</dd></dl>", 18),
	"--option\n\n    one two three\n    four five\n");
    ok &= expect("lists and continuation",
	render("<ul><li><p>first item with enough words</p><ul><li>nested</li></ul>"
	    "<p>continued</p></li><li>second</li></ul>", 25),
	"  * first item with\n    enough words\n\n      * nested\n\n    continued\n  * second\n");
    ok &= expect("ordered list start", render("<ol start='9'><li>nine</li><li>ten</li></ol>", 80),
	"  9. nine\n  10. ten\n");
    ok &= expect("list source whitespace", render("<ul><li>\n<p>item</p>\n</li></ul>", 80),
	"  * item\n");
    ok &= expect("literal whitespace",
	render("<pre>  a\tb\n\n  an intentionally long literal line\n</pre>", 12),
	"  a     b\n\n  an intentionally long literal line\n");
    ok &= expect("links and images",
	render("<p><a href='https://example.org'>reference</a> "
	    "<a href='#local'>local</a> <img alt='diagram' src='missing.png'></p>", 80),
	"reference <https://example.org> local [diagram]\n");
    const std::string table = "<table><tr><th>Key</th><th>Value</th></tr>"
	"<tr><td>A</td><td><p>two<b>three</b>four</p></td></tr></table>";
    ok &= expect("table columns", render(table, 40),
	"Key   Value\nA     twothreefour\n");
    ok &= expect("narrow table", render(table, 16),
	"Key: A\nValue:\ntwothreefour\n");
    ok &= expect("spanning table content",
	render("<table><tr><th colspan='2'>Heading</th></tr><tr><td>A</td><td>B</td></tr></table>", 80),
	"Heading\n\nA\nB\n");
    ok &= expect("table caption", render("<table><caption>Caption</caption><tr><td>A</td></tr></table>", 80),
	"Caption\nA\n");
    ok &= expect("header-only narrow table", render("<table><tr><th>Long heading</th></tr></table>", 6),
	"Long\nheadin\ng\n");
    ok &= expect("document controls",
	render(std::string("<p>safe") + char(27) + "[31m" + char(8) + "end</p>", 80),
	"safe\xef\xbf\xbd[31m\xef\xbf\xbd" "end\n");
    ok &= expect("Unicode entities", render("<p>caf&#233; &#x4e2d; e&#x301;</p>", 80),
	"caf\xc3\xa9 \xe4\xb8\xad e\xcc\x81\n");
    brlman::html_document styled("<p>a<b>b</b>c</p>");
    auto lines = styled.format(80);
    ok &= expect("style spans", brlman::line_text(lines[0], 0, 80, true),
	"a\033[0m\033[1mb\033[0mc");
    ok &= expect("clipping", brlman::line_text(lines[0], 1, 1), "b");
    brlman::html_document wide("<p>abcdefgh</p>");
    ok &= expect("repeatable layout", brlman::line_text(wide.format(3)[0]), "abc");
    ok &= expect("repeatable wider layout", brlman::line_text(wide.format(80)[0]), "abcdefgh");
    brlman::text_line unicode_line{{{"\xe4\xb8\xad" "e\xcc\x81x", brlman::normal}}, 0};
    ok &= expect("wide-character clipping", brlman::line_text(unicode_line, 1, 3), " e\xcc\x81x");
    ok &= expect("clipped combining base", brlman::line_text(unicode_line, 3, 1), "x");
    std::setlocale(LC_CTYPE, "C");
    ok &= expect("non-UTF-8 locale", render("<p>caf&#233; word</p>", 5),
	"caf\xc3\xa9\nword\n");
    std::setlocale(LC_CTYPE, "");
    brlman::html_document malformed("<p>start<b>bold<p>end");
    if (malformed.format(20).empty()) {
	std::cerr << "HTML recovery discarded content\n";
	ok = false;
    }
    return ok ? 0 : 1;
}

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
