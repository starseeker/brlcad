/*                   F A S T _ F A C E . H
 * BRL-CAD
 *
 * Copyright (c) 2007-2026 United States Government as represented by
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
/** @file cdt/fast_face.h
 * Private face-local display tessellation results.
 */
#ifndef LIBBREP_CDT_FAST_FACE_H
#define LIBBREP_CDT_FAST_FACE_H

#include "common.h"

#include <string>
#include <vector>

#include "brep/cdt.h"
#include "fast_budget.h"

struct fast_face_mesh {
    std::vector<int> faces;
    std::vector<fastf_t> normals;
    std::vector<fastf_t> points;
    std::vector<const void *> sources;

    void clear()
    {
	faces.clear();
	normals.clear();
	points.clear();
	sources.clear();
    }

    fast_output_size storage() const
    {
	return {faces.size() * sizeof(int) +
	    (normals.size() + points.size()) * sizeof(fastf_t) +
	    sources.size() * sizeof(const void *), points.size() / 3};
    }
};

enum fast_face_outcome {
    FAST_FACE_UNATTEMPTED,
    FAST_FACE_FAILED,
    FAST_FACE_MEMORY_LIMIT,
    FAST_FACE_EXCEPTION,
    FAST_FACE_COMPLETED,
    FAST_FACE_SKIPPED_DEGENERATE,
    FAST_FACE_SKIPPED_TOLERANCE
};

struct fast_cdt_face_result {
    fast_face_mesh mesh;
    fast_face_outcome outcome = FAST_FACE_UNATTEMPTED;
    int diagnostic_result = BREP_CDT_RESULT_UNATTEMPTED;
    int diagnostic_stage = BREP_CDT_STAGE_NONE;
    std::string diagnostic_message;

    bool completed() const
    {
	return outcome == FAST_FACE_COMPLETED || skipped_degenerate() ||
	    skipped_tolerance();
    }
    bool failed() const
    {
	return outcome == FAST_FACE_FAILED || outcome == FAST_FACE_MEMORY_LIMIT ||
	    outcome == FAST_FACE_EXCEPTION;
    }
    bool skipped_degenerate() const { return outcome == FAST_FACE_SKIPPED_DEGENERATE; }
    bool skipped_tolerance() const { return outcome == FAST_FACE_SKIPPED_TOLERANCE; }
};


double fast_brep_diagonal(const ON_Brep *brep);

fast_cdt_face_result fast_evaluate_face(const ON_BrepFace &face,
    const struct bg_tess_tol *ttol, const struct bn_tol *tol,
    double model_diagonal, const struct brep_cdt_fast_options *options);

#endif

/*
 * Local Variables:
 * tab-width: 8
 * mode: C++
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
