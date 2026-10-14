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
#include "library_script.h"

#include <tk.h>

#if defined(_WIN32)
#define TKLITEHTML_EXPORT __declspec(dllexport)
#elif defined(__GNUC__)
#define TKLITEHTML_EXPORT __attribute__((visibility("default")))
#else
#define TKLITEHTML_EXPORT
#endif

extern "C" TKLITEHTML_EXPORT int Tklitehtml_Init(Tcl_Interp* interp)
{
    if(Tcl_InitStubs(interp, "8.6", 0) == nullptr || Tk_InitStubs(interp, "8.6", 0) == nullptr)
    {
        return TCL_ERROR;
    }

    if(Tcl_PkgPresent(interp, "tklitehtml", "0.1.0", 1) != nullptr)
    {
        return TCL_OK;
    }
    Tcl_ResetResult(interp);

    if(Tcl_FindNamespace(interp, "::tklitehtml", nullptr, 0) == nullptr &&
       Tcl_CreateNamespace(interp, "::tklitehtml", nullptr, nullptr) == nullptr)
    {
        return TCL_ERROR;
    }

    if(tklitehtml::register_widget_command(interp) != TCL_OK ||
       Tcl_EvalEx(interp, tklitehtml_library_script,
                  static_cast<int>(sizeof(tklitehtml_library_script) - 1), TCL_EVAL_GLOBAL) != TCL_OK)
    {
        return TCL_ERROR;
    }

    return Tcl_PkgProvide(interp, "tklitehtml", "0.1.0");
}
