/*                   B R E P _ C D T . C P P
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
/** @file cdt/fast.cpp
 * Bounded scheduling, adaptive display quality, and ordered mesh assembly.
 */

#include "common.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <mutex>
#include <new>
#include <vector>

#include "bu/datetime.h"
#include "bu/env.h"
#include "bu/malloc.h"
#include "bu/parallel.h"
#include "bu/log.h"
#include "brep/cdt.h"
#include "fast_face.h"
#include "heal.h"

static const size_t FAST_CDT_DEFAULT_WORKING_BYTES =
    (size_t)1024 * 1024 * 1024;
static const size_t FAST_CDT_BASE_FACE_WORKING_BYTES =
    (size_t)16 * 1024 * 1024;
static const size_t FAST_CDT_NURBS_CV_WORKING_BYTES =
    (size_t)64 * 1024;
static const size_t FAST_CDT_TRIM_SPAN_WORKING_BYTES =
    (size_t)1024;
static const size_t FAST_CDT_DEFAULT_MAX_TRIANGLES =
    (size_t)256 * 1024;
static const double FAST_CDT_DEFAULT_RELATIVE_TOLERANCE = 0.01;
static const double FAST_CDT_MAX_RELATIVE_TOLERANCE = 0.5;
static const double FAST_CDT_DEFAULT_COARSE_RELATIVE_TOLERANCE =
    FAST_CDT_MAX_RELATIVE_TOLERANCE;
static const double FAST_CDT_RELATIVE_PROGRESS_ROUNDOFF = 1.0e-12;
static const int FAST_CDT_MAX_REFINEMENT_PASSES = 12;
static const double FAST_CDT_DEFAULT_AREA_CHANGE_TOLERANCE = 0.0025;


static size_t
fast_saturating_add(size_t first, size_t second, size_t limit)
{
    if (first >= limit || second >= limit - first)
	return limit;
    return first + second;
}


static size_t
fast_saturating_multiply(size_t first, size_t second, size_t limit)
{
    if (!first || !second)
	return 0;
    if (first >= limit / second)
	return limit;
    return first * second;
}


static size_t
fast_face_working_estimate(const ON_BrepFace &face, size_t limit)
{
    if (!limit)
	return 0;
    size_t estimate = std::min(FAST_CDT_BASE_FACE_WORKING_BYTES, limit);
    const ON_NurbsSurface *nurbs = ON_NurbsSurface::Cast(face.SurfaceOf());
    if (nurbs) {
	ON_BoundingBox trim_box = ON_BoundingBox::EmptyBoundingBox;
	bool have_trim_box = false;
	bool trim_box_failed = false;
	for (int loop_index = 0; loop_index < face.LoopCount(); ++loop_index) {
	    const ON_BrepLoop *loop = face.Loop(loop_index);
	    if (!loop)
		continue;
	    for (int trim_index = 0; trim_index < loop->TrimCount();
		    ++trim_index) {
		const ON_BrepTrim *trim = loop->Trim(trim_index);
		if (!trim || !trim->GetBoundingBox(trim_box,
			trim_box.IsValid())) {
		    trim_box_failed = true;
		    continue;
		}
		have_trim_box = trim_box.IsValid();
	    }
	}
	size_t relevant_control_points[2] = {
	    (size_t)std::max(0, nurbs->CVCount(0)),
	    (size_t)std::max(0, nurbs->CVCount(1))
	};
	if (have_trim_box && !trim_box_failed) {
	    for (int direction = 0; direction < 2; ++direction) {
		const int span_count = nurbs->SpanCount(direction);
		std::vector<double> spans((size_t)std::max(0,
		    span_count) + 1);
		if (span_count <= 0 ||
			!nurbs->GetSpanVector(direction, spans.data()))
		    continue;
		const double minimum = direction ? trim_box.Min().y :
		    trim_box.Min().x;
		const double maximum = direction ? trim_box.Max().y :
		    trim_box.Max().x;
		if (!std::isfinite(minimum) || !std::isfinite(maximum) ||
			minimum > maximum)
		    continue;
		size_t overlapping_spans = 0;
		for (int span = 0; span < span_count; ++span) {
		    if (spans[(size_t)span + 1] >= minimum &&
			    spans[(size_t)span] <= maximum)
			overlapping_spans++;
		}
		if (overlapping_spans) {
		    relevant_control_points[direction] = std::min(
			relevant_control_points[direction],
			overlapping_spans +
			(size_t)std::max(0, nurbs->Order(direction) - 1));
		}
	    }
	}
	const size_t control_points = fast_saturating_multiply(
	    relevant_control_points[0], relevant_control_points[1], limit);
	estimate = fast_saturating_add(estimate,
	    fast_saturating_multiply(control_points,
		FAST_CDT_NURBS_CV_WORKING_BYTES, limit), limit);
    }
    size_t trim_spans = 0;
    for (int loop_index = 0; loop_index < face.LoopCount(); ++loop_index) {
	const ON_BrepLoop *loop = face.Loop(loop_index);
	if (!loop)
	    continue;
	for (int trim_index = 0; trim_index < loop->TrimCount(); ++trim_index) {
	    const ON_BrepTrim *trim = loop->Trim(trim_index);
	    if (!trim)
		continue;
	    trim_spans = fast_saturating_add(trim_spans,
		(size_t)std::max(1, trim->SpanCount()), limit);
	}
    }
    estimate = fast_saturating_add(estimate,
	fast_saturating_multiply(trim_spans,
	    FAST_CDT_TRIM_SPAN_WORKING_BYTES, limit), limit);
    return std::max((size_t)1, std::min(estimate, limit));
}


struct fast_work_budget {
    explicit fast_work_budget(size_t byte_limit) : limit(byte_limit) {}

    bool acquire(size_t amount, std::atomic<bool> &stop,
	    std::atomic<bool> &hit_time_limit, int64_t deadline)
    {
	amount = std::min(std::max((size_t)1, amount), limit);
	std::unique_lock<std::mutex> lock(mutex);
	while (!stop.load() && in_use > limit - amount) {
	    if (deadline > 0 && bu_gettime() >= deadline) {
		hit_time_limit = true;
		stop = true;
		condition.notify_all();
		return false;
	    }
	    condition.wait_for(lock, std::chrono::milliseconds(10));
	}
	if (stop.load())
	    return false;
	in_use += amount;
	peak = std::max(peak, in_use);
	return true;
    }

    void release(size_t amount)
    {
	amount = std::min(std::max((size_t)1, amount), limit);
	std::lock_guard<std::mutex> lock(mutex);
	in_use = amount > in_use ? 0 : in_use - amount;
	condition.notify_all();
    }

    size_t peak_bytes()
    {
	std::lock_guard<std::mutex> lock(mutex);
	return peak;
    }

    size_t limit;
    size_t in_use = 0;
    size_t peak = 0;
    std::mutex mutex;
    std::condition_variable condition;
};


struct fast_work_reservation {
    fast_work_reservation(fast_work_budget *budget, size_t amount) :
	work_budget(budget), reserved(amount) {}
    ~fast_work_reservation()
    {
	if (work_budget && reserved)
	    work_budget->release(reserved);
    }

    fast_work_budget *work_budget;
    size_t reserved;
};


struct fast_cdt_parallel_state {
    const ON_Brep *brep;
    const struct bg_tess_tol *ttol;
    const struct bg_tess_tol *planar_ttol;
    const struct bn_tol *tol;
    const struct brep_cdt_fast_options *options;
    double model_diagonal;
    std::vector<fast_cdt_face_result> *results;
    const std::vector<size_t> *working_estimates;
    std::atomic<int> next_face;
    int end_face;
    fast_output_budget *output_budget;
    std::atomic<bool> stop;
    std::atomic<bool> hit_time_limit;
    std::atomic<bool> hit_memory_limit;
    std::atomic<bool> hit_point_limit;
    int64_t deadline;
    bool allow_partial;
    fast_work_budget *work_budget;
};

static void
fast_cdt_face_worker(int UNUSED(cpu), void *data)
{
    fast_cdt_parallel_state *state = (fast_cdt_parallel_state *)data;
    const int face_cnt = state->end_face;

    for (;;) {
	if (state->stop.load())
	    return;
	const int face_index = state->next_face.fetch_add(1);
	if (face_index >= face_cnt)
	    return;
	if (state->deadline > 0 && bu_gettime() >= state->deadline) {
	    state->hit_time_limit = true;
	    state->stop = true;
	    return;
	}

	const ON_BrepFace& face = state->brep->m_F[face_index];
	const struct bg_tess_tol *face_ttol = state->ttol;
	const ON_Surface *surface = face.SurfaceOf();
	if (state->planar_ttol && surface &&
		surface->IsPlanar(NULL, state->tol->dist))
	    face_ttol = state->planar_ttol;
	const bool requested_tolerance = !state->planar_ttol ||
	    face_ttol == state->planar_ttol;
	const size_t working_estimate = state->working_estimates ?
	    (*state->working_estimates)[(size_t)face_index] :
	    fast_face_working_estimate(face, state->work_budget->limit);
	if (!state->work_budget->acquire(working_estimate, state->stop,
		state->hit_time_limit, state->deadline))
	    return;
	fast_work_reservation reservation(state->work_budget,
	    working_estimate);
	fast_cdt_face_result &result = (*state->results)[face_index];
	result = fast_evaluate_face(face, face_ttol, state->tol,
	    state->model_diagonal, state->options);
	if (result.outcome == FAST_FACE_MEMORY_LIMIT) {
	    state->hit_memory_limit = true;
	    state->stop = true;
	}
	if (result.failed()) {
	    if (!state->allow_partial)
		state->stop = true;
	    continue;
	}
	if (result.skipped_degenerate())
	    continue;
	if (result.skipped_tolerance()) {
	    /* Coarse sampling can collapse a face that the requested tolerance
	     * resolves.  Leave it incomplete for the recovery pass. */
	    if (!requested_tolerance)
		result.outcome = FAST_FACE_FAILED;
	    continue;
	}

	const auto admission = state->output_budget->replace({},
	    result.mesh.storage());
	if (admission != fast_output_budget::accepted) {
	    if (admission == fast_output_budget::byte_limit)
		state->hit_memory_limit = true;
	    else
		state->hit_point_limit = true;
	    state->stop = true;
	    result.mesh.clear();
	    result.outcome = FAST_FACE_FAILED;
	    return;
	}
    }
}


static double
fast_face_mesh_area(const fast_face_mesh &output)
{
    double area = 0.0;
    for (size_t triangle = 0; triangle + 2 < output.faces.size();
	    triangle += 3) {
	ON_3dPoint points[3];
	bool valid = true;
	for (size_t corner = 0; corner < 3; ++corner) {
	    const int index = output.faces[triangle + corner];
	    if (index < 0 || (size_t)index * 3 + 2 >= output.points.size()) {
		valid = false;
		break;
	    }
	    points[corner] = ON_3dPoint(output.points[(size_t)index * 3],
		output.points[(size_t)index * 3 + 1],
		output.points[(size_t)index * 3 + 2]);
	}
	if (!valid)
	    continue;
	const double triangle_area = 0.5 * ON_CrossProduct(
	    points[1] - points[0], points[2] - points[0]).Length();
	if (std::isfinite(triangle_area))
	    area += triangle_area;
    }
    return area;
}


static void
fast_bbox_include(ON_BoundingBox &bbox, const ON_3dPoint &point)
{
    if (!point.IsValid() || !std::isfinite(point.x) ||
	    !std::isfinite(point.y) || !std::isfinite(point.z))
	return;
    if (bbox.IsValid())
	bbox.Set(point, true);
    else
	bbox = ON_BoundingBox(point, point);
}


static ON_BoundingBox
fast_edge_boundary_box(const ON_BrepEdge &edge)
{
    ON_BoundingBox bbox = ON_BoundingBox::EmptyBoundingBox;
    const ON_Interval domain = edge.Domain();
    if (!domain.IsIncreasing())
	return bbox;

    /* Proxy edges are sampled in their native spans so extrema represented
     * by the actual curve, rather than potentially remote NURBS control
     * points, define the guard.  Degree-one spans need only their endpoints;
     * avoid a separate IsLinear query because it may be more expensive than
     * the bounded sampling itself on large imported models. */
    int span_count = edge.SpanCount();
    std::vector<double> spans;
    if (span_count > 0 && span_count <= 4096) {
	spans.resize((size_t)span_count + 1);
	if (!edge.GetSpanVector(spans.data()))
	    spans.clear();
    }
    if (spans.empty()) {
	span_count = 1;
	spans.push_back(domain.Min());
	spans.push_back(domain.Max());
    }
    const int samples_per_span = edge.Degree() <= 1 ? 1 :
	std::max(1, std::min(4, 128 / span_count));
    for (int span = 0; span < span_count; ++span) {
	const ON_Interval interval(spans[(size_t)span],
	    spans[(size_t)span + 1]);
	for (int sample = 0; sample <= samples_per_span; ++sample)
	    fast_bbox_include(bbox, edge.PointAt(interval.ParameterAt(
		(double)sample / (double)samples_per_span)));
    }
    return bbox;
}


static std::vector<ON_BoundingBox>
fast_face_boundary_boxes(const ON_Brep *brep)
{
    std::vector<ON_BoundingBox> edge_boxes((size_t)brep->m_E.Count(),
	ON_BoundingBox::EmptyBoundingBox);
    for (int edge_index = 0; edge_index < brep->m_E.Count(); ++edge_index)
	edge_boxes[(size_t)edge_index] = fast_edge_boundary_box(
	    brep->m_E[edge_index]);

    std::vector<ON_BoundingBox> face_boxes((size_t)brep->m_F.Count(),
	ON_BoundingBox::EmptyBoundingBox);
    for (int face_index = 0; face_index < brep->m_F.Count(); ++face_index) {
	ON_BoundingBox &bbox = face_boxes[(size_t)face_index];
	const ON_BrepFace &face = brep->m_F[face_index];
	for (int loop_index = 0; loop_index < face.LoopCount(); ++loop_index) {
	    const ON_BrepLoop *loop = face.Loop(loop_index);
	    if (!loop)
		continue;
	    for (int trim_index = 0; trim_index < loop->TrimCount();
		    ++trim_index) {
		const ON_BrepTrim *trim = loop->Trim(trim_index);
		if (!trim)
		    continue;
		const ON_BrepVertex *first = trim->Vertex(0);
		const ON_BrepVertex *second = trim->Vertex(1);
		if (first)
		    fast_bbox_include(bbox, first->Point());
		if (second)
		    fast_bbox_include(bbox, second->Point());
		if (trim->m_ei < 0 ||
			trim->m_ei >= (int)edge_boxes.size())
		    continue;
		const ON_BoundingBox &edge_box =
		    edge_boxes[(size_t)trim->m_ei];
		if (!edge_box.IsValid())
		    continue;
		if (bbox.IsValid())
		    bbox.Union(edge_box);
		else
		    bbox = edge_box;
	    }
	}
    }
    return face_boxes;
}


struct fast_boundary_coverage {
    double fraction = 1.0;
    bool covered = true;
};


static fast_boundary_coverage
fast_face_boundary_coverage(const fast_face_mesh &output,
	const ON_BoundingBox &reference, double tolerance)
{
    fast_boundary_coverage coverage;
    if (!reference.IsValid())
	return coverage;

    ON_BoundingBox mesh = ON_BoundingBox::EmptyBoundingBox;
    /* Output construction emits only triangle vertices, so scanning the
     * compact point array avoids revisiting shared vertices for every face
     * corner on large meshes. */
    for (size_t point = 0; point + 2 < output.points.size(); point += 3)
	fast_bbox_include(mesh, ON_3dPoint(output.points[point],
	    output.points[point + 1], output.points[point + 2]));
    if (!mesh.IsValid()) {
	coverage.fraction = 0.0;
	coverage.covered = false;
	return coverage;
    }

    coverage.fraction = 1.0;
    for (int axis = 0; axis < 3; ++axis) {
	const double extent = reference.m_max[axis] -
	    reference.m_min[axis];
	if (!std::isfinite(extent) || extent <= 2.0 * tolerance)
	    continue;
	const double low = std::max(0.0, mesh.m_min[axis] -
	    reference.m_min[axis] - tolerance);
	const double high = std::max(0.0, reference.m_max[axis] -
	    mesh.m_max[axis] - tolerance);
	const double missing = low + high;
	coverage.fraction = std::min(coverage.fraction,
	    std::max(0.0, 1.0 - missing / extent));
	coverage.covered = coverage.covered && missing <= 0.0;
    }
    return coverage;
}


static struct bg_tess_tol
fast_scaled_tolerance(const struct bg_tess_tol *requested,
    double relative_tolerance)
{
    struct bg_tess_tol scaled = *requested;
    const double requested_relative = requested->rel > ON_ZERO_TOLERANCE ?
	requested->rel : FAST_CDT_DEFAULT_RELATIVE_TOLERANCE;
    relative_tolerance = std::max(requested_relative,
	std::min(FAST_CDT_MAX_RELATIVE_TOLERANCE, relative_tolerance));
    if (requested->rel > ON_ZERO_TOLERANCE ||
	    requested->abs <= ON_ZERO_TOLERANCE)
	scaled.rel = relative_tolerance;
    else
	scaled.abs = requested->abs * relative_tolerance /
	    requested_relative;
    if (requested->norm > ON_ZERO_TOLERANCE)
	scaled.norm = std::min(ON_PI / 2.0, requested->norm *
	    relative_tolerance / requested_relative);
    return scaled;
}


static size_t
fast_object_triangle_budget(const ON_Brep *brep,
    const struct brep_cdt_fast_options &options)
{
    const size_t hard_by_bytes = std::max((size_t)1,
	options.max_result_bytes / 96);
    const size_t hard_limit = std::max((size_t)1,
	std::min(options.max_triangles, hard_by_bytes));
    size_t target = std::min((size_t)4096, hard_limit);
    target = fast_saturating_add(target,
	fast_saturating_multiply((size_t)brep->m_F.Count(), 32,
	    hard_limit), hard_limit);
    target = fast_saturating_add(target,
	fast_saturating_multiply((size_t)brep->m_E.Count(), 4,
	    hard_limit), hard_limit);
    return std::max((size_t)1, std::min(target, hard_limit));
}


static std::vector<size_t>
fast_face_triangle_budgets(const std::vector<size_t> &working_estimates,
    size_t object_budget)
{
    const size_t face_count = working_estimates.size();
    std::vector<size_t> budgets(face_count, 1);
    if (!face_count)
	return budgets;
    const size_t minimum = std::max((size_t)1,
	std::min((size_t)16, object_budget / face_count));
    const size_t reserved = std::min(object_budget, minimum * face_count);
    const size_t distributable = object_budget - reserved;
    std::vector<double> weights(face_count, 1.0);
    double weight_sum = 0.0;
    for (size_t face_index = 0; face_index < face_count; ++face_index) {
	weights[face_index] = std::sqrt(1.0 +
	    (double)working_estimates[face_index] /
	    (1024.0 * 1024.0));
	weight_sum += weights[face_index];
    }
    size_t assigned = 0;
    for (size_t face_index = 0; face_index < face_count; ++face_index) {
	const size_t share = weight_sum > 0.0 ? (size_t)std::floor(
	    distributable * weights[face_index] / weight_sum) : 0;
	budgets[face_index] = minimum + share;
	assigned += budgets[face_index];
    }
    for (size_t face_index = 0; assigned < object_budget; ++assigned)
	budgets[face_index++ % face_count]++;
    return budgets;
}


struct fast_face_refinement {
    double area = 0.0;
    double boundary_coverage = 1.0;
    double accepted_relative_tolerance = 0.0;
    int stable_passes = 0;
    bool boundary_covered = true;
    bool area_converged = false;
    bool budget_limited = false;
};


struct fast_cdt_refine_state {
    const ON_Brep *brep;
    const struct bg_tess_tol *ttol;
    const struct bn_tol *tol;
    const struct brep_cdt_fast_options *options;
    double model_diagonal;
    double relative_tolerance;
    double area_change_tolerance;
    std::vector<fast_cdt_face_result> *results;
    std::vector<fast_face_refinement> *refinement;
    const std::vector<ON_BoundingBox> *boundary_boxes;
    const std::vector<double> *boundary_tolerances;
    const std::vector<size_t> *triangle_budgets;
    const std::vector<size_t> *working_estimates;
    std::atomic<int> next_face;
    std::atomic<bool> stop;
    std::atomic<bool> time_limited;
    fast_output_budget *output_budget;
    bool coarsening;
    bool incomplete_only;
    int64_t deadline;
    fast_work_budget *work_budget;
};


static void
fast_cdt_refine_worker(int UNUSED(cpu), void *data)
{
    fast_cdt_refine_state *state = (fast_cdt_refine_state *)data;
    const int face_count = state->brep->m_F.Count();
    for (;;) {
	if (state->stop.load())
	    return;
	const int face_index = state->next_face.fetch_add(1);
	if (face_index >= face_count)
	    return;
	if (state->deadline > 0 && bu_gettime() >= state->deadline) {
	    state->time_limited = true;
	    state->stop = true;
	    return;
	}

	fast_face_refinement &quality =
	    (*state->refinement)[(size_t)face_index];
	fast_cdt_face_result &current =
	    (*state->results)[(size_t)face_index];
	if (state->incomplete_only && current.completed())
	    continue;
	if (current.skipped_degenerate() || current.skipped_tolerance() ||
		(state->coarsening && (!quality.budget_limited ||
		quality.area_converged || !quality.boundary_covered)) ||
		(!state->coarsening &&
		(quality.area_converged || quality.budget_limited) &&
		quality.boundary_covered))
	    continue;
	if (!state->coarsening && current.completed() &&
		state->relative_tolerance >=
		quality.accepted_relative_tolerance * (1.0 - FAST_CDT_RELATIVE_PROGRESS_ROUNDOFF))
	    continue;

	const ON_BrepFace &face = state->brep->m_F[face_index];
	const size_t working_estimate =
	    (*state->working_estimates)[(size_t)face_index];
	std::atomic<bool> hard_time_limit(false);
	if (!state->work_budget->acquire(working_estimate, state->stop,
		hard_time_limit, state->deadline)) {
	    if (hard_time_limit.load())
		state->time_limited = true;
	    return;
	}
	fast_work_reservation reservation(state->work_budget,
	    working_estimate);

	fast_cdt_face_result candidate = fast_evaluate_face(face, state->ttol,
	    state->tol, state->model_diagonal, state->options);
	if (candidate.skipped_degenerate()) {
	    if (!current.completed())
		current = std::move(candidate);
	    continue;
	}
	if (candidate.skipped_tolerance()) {
	    if (state->incomplete_only && !current.completed())
		current = std::move(candidate);
	    continue;
	}
	if (!candidate.completed())
	    continue;

	const size_t triangles = candidate.mesh.faces.size() / 3;
	const size_t triangle_budget =
	    (*state->triangle_budgets)[(size_t)face_index];
	const fast_boundary_coverage candidate_coverage =
	    fast_face_boundary_coverage(candidate.mesh,
		(*state->boundary_boxes)[(size_t)face_index],
		(*state->boundary_tolerances)[(size_t)face_index]);
	if (state->coarsening) {
	    if (!candidate_coverage.covered || !current.completed() ||
		triangles >= current.mesh.faces.size() / 3)
		continue;
	} else if (current.completed()) {
	    if (quality.boundary_covered && !candidate_coverage.covered)
		continue;
	    if (!quality.boundary_covered &&
		    !candidate_coverage.covered &&
		    candidate_coverage.fraction <=
		    quality.boundary_coverage + 1.0e-12)
		continue;
	    if (quality.boundary_covered && triangles > triangle_budget) {
		quality.budget_limited = true;
		continue;
	    }
	}

	const fast_output_size previous = current.completed() ?
	    current.mesh.storage() : fast_output_size{};
	if (state->output_budget->replace(previous,
		candidate.mesh.storage()) != fast_output_budget::accepted) {
	    quality.budget_limited = true;
	    continue;
	}
	const double area = fast_face_mesh_area(candidate.mesh);
	if (!state->coarsening && candidate_coverage.covered &&
		current.completed() &&
		!current.skipped_degenerate() &&
		quality.area > ON_ZERO_TOLERANCE && area > 0.0) {
	    const double change = fabs(area - quality.area) /
		std::max(area, quality.area);
	    quality.stable_passes =
		change <= state->area_change_tolerance ?
		quality.stable_passes + 1 : 0;
	    quality.area_converged = quality.stable_passes >= 2;
	} else if (!candidate_coverage.covered) {
	    quality.stable_passes = 0;
	    quality.area_converged = false;
	}
	current = std::move(candidate);
	quality.area = area;
	quality.boundary_coverage = candidate_coverage.fraction;
	quality.boundary_covered = candidate_coverage.covered;
	quality.accepted_relative_tolerance = state->relative_tolerance;
	quality.budget_limited = triangles > triangle_budget;
	if (state->coarsening) {
	    quality.stable_passes = 0;
	    quality.area_converged = false;
	}
    }
}

struct fast_refinement_report {
    int passes = 0;
    bool time_limited = false;
};

static fast_refinement_report
fast_refine_faces(fast_cdt_refine_state &state, double coarse_relative,
	double requested_relative, int64_t start_time)
{
    fast_refinement_report report;
    const ON_Brep *brep = state.brep;
    const int brep_face_count = brep->m_F.Count();
    const auto *ttol = state.ttol;
    const auto *tol = state.tol;
    const auto &options = *state.options;
    const int64_t deadline = state.deadline;
    auto &face_results = *state.results;
    auto &refinement = *state.refinement;
    const auto &face_boundary_boxes = *state.boundary_boxes;
    const auto &face_boundary_tolerances = *state.boundary_tolerances;
    const auto &face_triangle_budgets = *state.triangle_budgets;
    for (int face_index = 0; face_index < brep_face_count;
	    face_index++) {
	fast_cdt_face_result &result =
	    face_results[(size_t)face_index];
	if (!result.completed() || result.skipped_degenerate() ||
		result.skipped_tolerance())
	    continue;
	fast_face_refinement &quality =
	    refinement[(size_t)face_index];
	quality.area = fast_face_mesh_area(result.mesh);
	const fast_boundary_coverage coverage =
	    fast_face_boundary_coverage(result.mesh,
		face_boundary_boxes[(size_t)face_index],
		face_boundary_tolerances[(size_t)face_index]);
	quality.boundary_coverage = coverage.fraction;
	quality.boundary_covered = coverage.covered;
	const ON_Surface *surface =
	    brep->m_F[face_index].SurfaceOf();
	const bool planar = surface &&
	    surface->IsPlanar(NULL, tol->dist);
	quality.accepted_relative_tolerance = planar ?
	    requested_relative : coarse_relative;
	quality.area_converged = planar && coverage.covered;
	quality.budget_limited = result.mesh.faces.size() / 3 >
	    face_triangle_budgets[(size_t)face_index];
    }
    const int64_t initial_elapsed = std::max((int64_t)1,
	bu_gettime() - start_time);
    int64_t previous_pass_elapsed = initial_elapsed;

    auto run_refinement_pass = [&](double relative, bool coarsening,
	    bool incomplete_only = false) {
	const int64_t now = bu_gettime();
	/* A recovery pass schedules only missing faces, so a full-object
	 * elapsed-time estimate would prevent the retry it is meant to make.
	 * Its workers still honor the absolute deadline before each face. */
	const int64_t expected = incomplete_only ? 1 :
	    std::max(initial_elapsed, previous_pass_elapsed);
	if (deadline > 0 && (now >= deadline ||
	    expected >= deadline - now)) {
	    report.time_limited = true;
	    return false;
	}
	const int64_t pass_start = now;
	struct bg_tess_tol pass_tolerance =
	    fast_scaled_tolerance(ttol, relative);
	state.ttol = &pass_tolerance;
	state.relative_tolerance = relative;
	state.next_face = 0;
	state.stop = false;
	state.time_limited = false;
	state.coarsening = coarsening;
	state.incomplete_only = incomplete_only;
	bu_parallel(fast_cdt_refine_worker, options.max_workers,
	    &state);
	state.ttol = ttol;
	report.passes++;
	previous_pass_elapsed = std::max((int64_t)1,
	    bu_gettime() - pass_start);
	if (state.time_limited.load())
	    report.time_limited = true;
	return !report.time_limited;
    };

    /* Tessellation success is not monotonic with tolerance: coarse sampling
     * can collapse a narrow boundary which is resolved at the requested
     * tolerance.  Recover missing faces there directly rather than spending
     * several passes walking toward the only known sufficient setting. */
    bool have_incomplete_face = false;
    for (const fast_cdt_face_result &result : face_results)
	have_incomplete_face = have_incomplete_face || !result.completed();
    if (have_incomplete_face)
	run_refinement_pass(requested_relative, false, true);

    /* A topology-derived share is a target, not permission to discard
     * a valid boundary.  First try progressively coarser surface samples
     * for oversized faces; an irreducible trim boundary remains tagged. */
    for (double relative = std::min(FAST_CDT_MAX_RELATIVE_TOLERANCE, coarse_relative * 2.0);
	    !report.time_limited &&
	    relative > coarse_relative * (1.0 + FAST_CDT_RELATIVE_PROGRESS_ROUNDOFF);
	    relative = std::min(FAST_CDT_MAX_RELATIVE_TOLERANCE, relative * 2.0)) {
	bool have_oversized_face = false;
	for (const fast_face_refinement &quality : refinement)
	    have_oversized_face = have_oversized_face ||
		(quality.budget_limited && !quality.area_converged &&
		quality.boundary_covered);
	if (!have_oversized_face)
	    break;
	if (!run_refinement_pass(relative, true))
	    break;
	if (relative >= FAST_CDT_MAX_RELATIVE_TOLERANCE)
	    break;
    }

    /* Walk toward the requested quality.  Area stability and the soft
     * triangle share may stop a face only after its realized triangles
     * cover the independently sampled authoritative boundary envelope. */
    double relative = coarse_relative * 0.5;
    for (int pass = 0; !report.time_limited && pass < FAST_CDT_MAX_REFINEMENT_PASSES; ++pass) {
	bool have_refinable_face = false;
	for (int face_index = 0; face_index < brep_face_count;
		++face_index) {
	    const fast_cdt_face_result &result =
		face_results[(size_t)face_index];
	    const fast_face_refinement &quality =
		refinement[(size_t)face_index];
	    const bool needs_refinement = !result.completed() ||
		(quality.accepted_relative_tolerance >
		requested_relative * (1.0 + FAST_CDT_RELATIVE_PROGRESS_ROUNDOFF) &&
		!quality.area_converged &&
		(!quality.budget_limited ||
		!quality.boundary_covered));
	    have_refinable_face = have_refinable_face ||
		(!result.skipped_degenerate() &&
		!result.skipped_tolerance() && needs_refinement);
	}
	if (!have_refinable_face)
	    break;
	relative = std::max(requested_relative, relative);
	if (relative >= coarse_relative * (1.0 - FAST_CDT_RELATIVE_PROGRESS_ROUNDOFF))
	    break;
	if (!run_refinement_pass(relative, false))
	    break;
	if (relative <= requested_relative * (1.0 + FAST_CDT_RELATIVE_PROGRESS_ROUNDOFF))
	    break;
	relative *= 0.5;
    }
    return report;
}

/* Healing is an optional recovery step even when drawing has no deadline.
 * Bound its topology search so a malformed import cannot monopolize drawing. */
static constexpr long FAST_CDT_HEALING_TIME_MS = 1000;

static void
fast_heal_failed_faces(const ON_Brep &source, int first_face, int end_face,
	const struct brep_cdt_fast_options &options,
	const struct bg_tess_tol *ttol, const struct bn_tol *tol,
	double model_diagonal, int64_t deadline,
	std::vector<fast_cdt_face_result> &results, fast_work_budget &work_budget,
	fast_output_budget &output_budget, std::vector<int> &recovered)
{
    /* Authoritative trim samples and source tokens belong to the original
     * topology.  An owned copy cannot silently change their identities. */
    if (options.trim_sample_count || options.trim_sample ||
	options.trim_sample_source || options.point_source ||
	options.preserve_pullback_samples)
	return;
    bool needed = false;
    for (int fi = first_face; fi < end_face; ++fi) {
	const fast_cdt_face_result &result = results[(size_t)fi];
	needed = needed || (result.failed() && !result.completed() &&
	    result.diagnostic_result == BREP_CDT_RESULT_INVALID_PSLG);
    }
    if (!needed)
	return;
    long healing_time = FAST_CDT_HEALING_TIME_MS;
    if (deadline > 0)
	healing_time = (long)std::min((int64_t)healing_time,
	    (deadline - bu_gettime()) / 1000);
    if (healing_time <= 0)
	return;

    /* Keep a face-sized scratch allowance beyond the B-Rep copies, while
     * leaving at least half the working budget for the subsequent mesher. */
    const size_t source_bytes = source.SizeOf();
    const size_t healing_bytes = std::min(options.max_working_bytes / 2,
	fast_saturating_add(fast_saturating_multiply(source_bytes, 2,
	    options.max_working_bytes), FAST_CDT_BASE_FACE_WORKING_BYTES,
	    options.max_working_bytes));
    if (!healing_bytes || source_bytes > healing_bytes / 2)
	return;
    std::atomic<bool> stop(false), time_limited(false);
    if (!work_budget.acquire(healing_bytes, stop, time_limited, deadline))
	return;
    fast_work_reservation healing_reservation(&work_budget, healing_bytes);
    cdt_healing healing;
    if (!cdt_heal_topology(source, tol->dist, options.max_points,
	    healing_bytes, healing_time, healing, false,
	    CDT_HEAL_DRAWABLE_COMPONENTS) || !healing.brep ||
	    healing.brep->m_F.Count() != source.m_F.Count())
	return;

    for (int fi : healing.faces) {
	if (fi < first_face || fi >= end_face)
	    continue;
	fast_cdt_face_result &current = results[(size_t)fi];
	if (current.completed() || !current.failed() ||
	    current.diagnostic_result != BREP_CDT_RESULT_INVALID_PSLG)
	    continue;
	if (deadline > 0 && bu_gettime() >= deadline)
	    break;
	const ON_BrepFace &face = healing.brep->m_F[fi];
	const size_t working_bytes = fast_face_working_estimate(face,
	    options.max_working_bytes - healing_bytes);
	if (!work_budget.acquire(working_bytes, stop, time_limited, deadline))
	    break;
	fast_work_reservation face_reservation(&work_budget, working_bytes);
	fast_cdt_face_result candidate = fast_evaluate_face(face, ttol, tol,
	    model_diagonal, &options);
	if (candidate.outcome == FAST_FACE_MEMORY_LIMIT)
	    throw std::bad_alloc();
	if (candidate.outcome == FAST_FACE_EXCEPTION) {
	    bu_log("BRep drawing topology recovery threw an unexpected exception\n");
	    break;
	}
	if (candidate.outcome != FAST_FACE_COMPLETED)
	    continue;
	if (output_budget.replace({}, candidate.mesh.storage()) !=
		fast_output_budget::accepted)
	    continue;
	recovered.push_back(fi);
	current = std::move(candidate);
    }
}

void
brep_cdt_fast_options_default(struct brep_cdt_fast_options *options)
{
    if (!options)
	return;
    options->max_workers = std::min((size_t)8, bu_avail_cpus());
    if (!options->max_workers)
	options->max_workers = 1;
    options->max_result_bytes = (size_t)512 * 1024 * 1024;
    options->max_points = (size_t)16 * 1024 * 1024;
    options->max_time_ms = 0;
    options->allow_partial = 1;
    options->face_status = NULL;
    options->face_status_data = NULL;
    options->face_diagnostic = NULL;
    options->face_diagnostic_data = NULL;
    options->face_output = NULL;
    options->face_output_data = NULL;
    options->trim_sample_count = NULL;
    options->trim_sample = NULL;
    options->trim_sample_data = NULL;
    options->trim_sample_source = NULL;
    options->point_source = NULL;
    options->point_source_data = NULL;
    options->max_working_bytes = FAST_CDT_DEFAULT_WORKING_BYTES;
    options->max_triangles = FAST_CDT_DEFAULT_MAX_TRIANGLES;
    options->adaptive_quality = 1;
    options->coarse_relative_tolerance =
	FAST_CDT_DEFAULT_COARSE_RELATIVE_TOLERANCE;
    options->area_change_tolerance =
	FAST_CDT_DEFAULT_AREA_CHANGE_TOLERANCE;
    options->preserve_pullback_samples = 0;
}

static struct brep_cdt_fast_options
fast_resolve_options(const struct brep_cdt_fast_options *user_options,
	int face_count)
{
    struct brep_cdt_fast_options defaults;
    brep_cdt_fast_options_default(&defaults);
    struct brep_cdt_fast_options options = user_options ? *user_options : defaults;
    if (!options.max_workers)
	options.max_workers = defaults.max_workers;
    if (!options.max_result_bytes)
	options.max_result_bytes = defaults.max_result_bytes;
    if (!options.max_points)
	options.max_points = defaults.max_points;
    if (!options.max_working_bytes)
	options.max_working_bytes = defaults.max_working_bytes;
    if (!options.max_triangles)
	options.max_triangles = defaults.max_triangles;
    if (!std::isfinite(options.coarse_relative_tolerance) ||
	options.coarse_relative_tolerance <= 0.0)
	options.coarse_relative_tolerance = defaults.coarse_relative_tolerance;
    if (!std::isfinite(options.area_change_tolerance) ||
	options.area_change_tolerance <= 0.0)
	options.area_change_tolerance = defaults.area_change_tolerance;
    options.max_workers = std::max((size_t)1,
	std::min(options.max_workers, (size_t)face_count));

    ssize_t available = bu_mem(BU_MEM_AVAIL, NULL);
    if (available > 0)
	options.max_result_bytes = std::min(options.max_result_bytes,
	    (size_t)available / 4);
    if (available > 0)
	options.max_working_bytes = std::min(options.max_working_bytes,
	    std::max((size_t)1, (size_t)available / 4));

    return options;
}

static void
fast_assemble_faces(fast_face_mesh &output, struct brep_cdt_fast_report &summary,
	const std::vector<fast_cdt_face_result> &face_results,
	const std::vector<fast_face_refinement> &refinement,
	int first_face, int end_face, const struct brep_cdt_fast_options &options,
	bool adaptive_quality, double requested_relative)
{
    auto &all_faces = output.faces;
    auto &all_norms = output.normals;
    auto &all_pnts = output.points;
    for (int fi = first_face; fi < end_face; fi++) {
	const fast_cdt_face_result &result = face_results[(size_t)fi];
	if (!result.completed()) {
	    summary.failed_faces++;
	    if (result.failed() && options.face_diagnostic)
		options.face_diagnostic(fi, result.diagnostic_result,
		    result.diagnostic_stage,
		    result.diagnostic_message.c_str(),
		    options.face_diagnostic_data);
	    if (options.face_status)
		options.face_status(fi, result.failed() ?
		    BREP_CDT_FAST_FACE_FAILED :
		    BREP_CDT_FAST_FACE_NOT_PROCESSED,
		    options.face_status_data);
	    continue;
	}
	summary.completed_faces++;
	if (result.skipped_degenerate()) {
	    summary.skipped_degenerate_faces++;
	    if (options.face_status)
		options.face_status(fi,
		    BREP_CDT_FAST_FACE_SKIPPED_DEGENERATE,
		    options.face_status_data);
	    continue;
	}
	if (result.skipped_tolerance()) {
	    summary.skipped_tolerance_faces++;
	    if (options.face_status)
		options.face_status(fi,
		    BREP_CDT_FAST_FACE_SKIPPED_TOLERANCE,
		    options.face_status_data);
	    continue;
	}
	const fast_face_refinement &quality = refinement[(size_t)fi];
	const bool approximated = adaptive_quality &&
	    (quality.budget_limited ||
	    !quality.boundary_covered ||
	    quality.accepted_relative_tolerance >
	    requested_relative * (1.0 + FAST_CDT_RELATIVE_PROGRESS_ROUNDOFF));
	if (approximated)
	    summary.approximated_faces++;
	if (quality.area_converged)
	    summary.area_converged_faces++;
	if (quality.budget_limited)
	    summary.triangle_budget_limited_faces++;
	if (!quality.boundary_covered)
	    summary.boundary_envelope_incomplete_faces++;
	if (options.face_status)
	    options.face_status(fi, approximated ?
		BREP_CDT_FAST_FACE_APPROXIMATED :
		BREP_CDT_FAST_FACE_COMPLETED,
		options.face_status_data);
	const size_t point_offset = all_pnts.size() / 3;
	if (options.face_output)
	    options.face_output(fi, all_faces.size() / 3,
		result.mesh.faces.size() / 3, point_offset,
		result.mesh.points.size() / 3, options.face_output_data);
	if (options.point_source) {
	    for (size_t point = 0; point < result.mesh.points.size() / 3; ++point) {
		const void *source = point < result.mesh.sources.size() ?
		    result.mesh.sources[point] : NULL;
		options.point_source(fi, point_offset + point, source,
		    options.point_source_data);
	    }
	}
	for (size_t i = 0; i < result.mesh.faces.size(); i++)
	    all_faces.push_back((int)point_offset + result.mesh.faces[i]);
	all_norms.insert(all_norms.end(), result.mesh.normals.begin(),
	    result.mesh.normals.end());
	all_pnts.insert(all_pnts.end(), result.mesh.points.begin(), result.mesh.points.end());
    }

    summary.result_bytes = output.storage().bytes;
}

int
brep_cdt_fast_ex(int **faces, int *face_cnt, vect_t **pnt_norms,
	point_t **pnts, int *pntcnt, const ON_Brep *brep, int index,
	const struct bg_tess_tol *ttol, const struct bn_tol *tol,
	const struct brep_cdt_fast_options *user_options,
	struct brep_cdt_fast_report *report)
{
    if (!faces || !face_cnt || !pnt_norms || !pnts || !pntcnt)
	return BREP_CDT_FAST_ERROR;

    *faces = NULL;
    *face_cnt = 0;
    *pnt_norms = NULL;
    *pnts = NULL;
    *pntcnt = 0;
    if (report)
	memset(report, 0, sizeof(*report));

    if (!brep || !ttol || !tol)
	return BREP_CDT_FAST_ERROR;

    const int brep_face_count = brep->m_F.Count();

    if (index < -1 || (index >= 0 && index >= brep_face_count))
	return BREP_CDT_FAST_ERROR;
    if (!brep_face_count)
	return BREP_CDT_FAST_OK;

    const struct brep_cdt_fast_options options =
	fast_resolve_options(user_options, brep_face_count);

    const int64_t start_time = bu_gettime();
    const int64_t deadline = options.max_time_ms > 0 ?
	start_time + (int64_t)options.max_time_ms * 1000 : 0;
    const double model_diagonal = fast_brep_diagonal(brep);
    const double requested_relative = ttol->rel > ON_ZERO_TOLERANCE ?
	ttol->rel : FAST_CDT_DEFAULT_RELATIVE_TOLERANCE;
    const bool adaptive_quality = index == -1 &&
	options.adaptive_quality != 0;
    const double coarse_relative = adaptive_quality ?
	std::max(requested_relative, options.coarse_relative_tolerance) :
	requested_relative;
    struct bg_tess_tol initial_tolerance = adaptive_quality ?
	fast_scaled_tolerance(ttol, coarse_relative) : *ttol;

    fast_face_mesh output;
    auto &all_faces = output.faces;
    auto &all_norms = output.normals;
    auto &all_pnts = output.points;

    const int first_face = (index == -1) ? 0 : index;
    const int end_face = (index == -1) ? brep_face_count : index + 1;
    std::vector<fast_cdt_face_result> face_results((size_t)brep_face_count);
    bool hit_time_limit = false;
    bool hit_memory_limit = false;
    bool hit_point_limit = false;
    fast_work_budget work_budget(options.max_working_bytes);
    fast_output_budget output_budget({options.max_result_bytes, options.max_points});
    const size_t triangle_budget = adaptive_quality ?
	fast_object_triangle_budget(brep, options) : options.max_triangles;
    std::vector<size_t> face_working_estimates((size_t)brep_face_count, 0);
    if (adaptive_quality) {
	for (int face_index = 0; face_index < brep_face_count; ++face_index)
	    face_working_estimates[(size_t)face_index] =
		fast_face_working_estimate(brep->m_F[face_index],
		    options.max_working_bytes);
    }
    std::vector<size_t> face_triangle_budgets = adaptive_quality ?
	fast_face_triangle_budgets(face_working_estimates, triangle_budget) :
	std::vector<size_t>((size_t)brep_face_count, options.max_triangles);
    std::vector<ON_BoundingBox> face_boundary_boxes = adaptive_quality ?
	fast_face_boundary_boxes(brep) :
	std::vector<ON_BoundingBox>((size_t)brep_face_count,
	    ON_BoundingBox::EmptyBoundingBox);
    std::vector<double> face_boundary_tolerances(
	(size_t)brep_face_count, tol->dist);
    if (adaptive_quality) {
	for (int face_index = 0; face_index < brep_face_count; ++face_index) {
	    const ON_BoundingBox &bbox =
		face_boundary_boxes[(size_t)face_index];
	    if (!bbox.IsValid())
		continue;
	    double allowance = tol->dist;
	    if (ttol->abs > ON_ZERO_TOLERANCE)
		allowance = std::max(allowance, ttol->abs);
	    if (ttol->rel > ON_ZERO_TOLERANCE)
		allowance = std::max(allowance,
		    ttol->rel * bbox.Diagonal().Length());
	    face_boundary_tolerances[(size_t)face_index] = allowance;
	}
    }
    std::vector<fast_face_refinement> refinement(
	(size_t)brep_face_count);

    /* Workers own face-local output.  Assemble it in face order only after
     * the pass finishes, preserving deterministic indices and callbacks. */
    fast_cdt_parallel_state initial_state = {
	brep, &initial_tolerance, adaptive_quality ? ttol : NULL,
	tol, &options, model_diagonal,
	&face_results, adaptive_quality ? &face_working_estimates : NULL,
	first_face, end_face, &output_budget,
	false, false, false, false,
	deadline, (bool)options.allow_partial, &work_budget
    };
    if (index == -1)
	bu_parallel(fast_cdt_face_worker, options.max_workers, &initial_state);
    else
	fast_cdt_face_worker(0, &initial_state);
    hit_time_limit = initial_state.hit_time_limit.load();
    hit_memory_limit = initial_state.hit_memory_limit.load();
    hit_point_limit = initial_state.hit_point_limit.load();
    if (index >= 0 && hit_time_limit &&
	face_results[(size_t)index].outcome == FAST_FACE_UNATTEMPTED)
	face_results[(size_t)index].outcome = FAST_FACE_FAILED;

    fast_refinement_report refinement_report;
    if (adaptive_quality && !hit_time_limit && !hit_memory_limit &&
	    !hit_point_limit) {
	fast_cdt_refine_state refine_state = {
	    brep, ttol, tol, &options, model_diagonal,
	    coarse_relative, options.area_change_tolerance, &face_results,
	    &refinement, &face_boundary_boxes, &face_boundary_tolerances,
	    &face_triangle_budgets, &face_working_estimates,
	    0, false, false, &output_budget, false, false, deadline, &work_budget
	};
	refinement_report = fast_refine_faces(refine_state, coarse_relative,
	    requested_relative, start_time);
    }

    if (!hit_time_limit && !hit_memory_limit && !hit_point_limit) {
	std::vector<int> recovered;
	try {
	    fast_heal_failed_faces(*brep,
		first_face, end_face, options, ttol, tol, model_diagonal,
		deadline, face_results, work_budget, output_budget, recovered);
	} catch (const std::bad_alloc &) {
	    hit_memory_limit = true;
	} catch (...) {
	    /* Optional healing must leave the original per-face diagnostics
	     * and any completed drawing available when recovery throws. */
	    bu_log("BRep drawing topology recovery threw an unexpected exception\n");
	}
	for (int fi : recovered) {
	    const fast_cdt_face_result &result = face_results[(size_t)fi];
	    fast_face_refinement &quality = refinement[(size_t)fi];
	    const fast_boundary_coverage coverage = fast_face_boundary_coverage(
		result.mesh, face_boundary_boxes[(size_t)fi],
		face_boundary_tolerances[(size_t)fi]);
	    quality.area = fast_face_mesh_area(result.mesh);
	    quality.boundary_coverage = coverage.fraction;
	    quality.boundary_covered = coverage.covered;
	    quality.accepted_relative_tolerance = requested_relative;
	    quality.area_converged = coverage.covered &&
		brep->m_F[fi].SurfaceOf()->IsPlanar(NULL, tol->dist);
	    quality.budget_limited = result.mesh.faces.size() / 3 >
		face_triangle_budgets[(size_t)fi];
	}
	if (deadline > 0 && bu_gettime() >= deadline)
	    hit_time_limit = true;
    }

    struct brep_cdt_fast_report summary = {};
    fast_assemble_faces(output, summary, face_results, refinement,
	first_face, end_face, options, adaptive_quality, requested_relative);

    summary.requested_faces = end_face - first_face;
    summary.hit_time_limit = hit_time_limit;
    summary.hit_memory_limit = hit_memory_limit;
    summary.hit_point_limit = hit_point_limit;
    summary.peak_working_bytes = work_budget.peak_bytes();
    summary.triangle_budget = triangle_budget;
    summary.refinement_passes = refinement_report.passes;
    summary.refinement_time_limited = refinement_report.time_limited;
    if (report)
	*report = summary;

    const bool hit_limit = hit_time_limit || hit_memory_limit ||
	hit_point_limit;
    if (summary.failed_faces && (!options.allow_partial || !summary.completed_faces))
	return hit_limit ? BREP_CDT_FAST_LIMIT : BREP_CDT_FAST_ERROR;

    /* A BRep containing only exact or tolerance-level empty faces has no
     * drawable output, but the request was still handled completely. */
    if (!all_faces.size() || !all_pnts.size() || !all_norms.size())
	return (!summary.failed_faces && summary.completed_faces == summary.skipped_degenerate_faces +
	    summary.skipped_tolerance_faces) ?
	    BREP_CDT_FAST_OK : BREP_CDT_FAST_ERROR;

    if (all_faces.size() / 3 > INT_MAX || all_pnts.size() / 3 > INT_MAX)
	return BREP_CDT_FAST_LIMIT;

    (*face_cnt) = (int)all_faces.size()/3;
    (*pntcnt) = (int)all_pnts.size()/3;
    (*faces) = (int *)bu_calloc(all_faces.size(), sizeof(int), "faces");
    (*pnt_norms) = (vect_t *)bu_calloc(all_norms.size(), sizeof(fastf_t), "normals");
    (*pnts) = (point_t *)bu_calloc(all_pnts.size()/3, sizeof(point_t), "pnts");
    std::copy(all_faces.begin(), all_faces.end(), *faces);
    std::copy(all_norms.begin(), all_norms.end(), (fastf_t *)*pnt_norms);
    std::copy(all_pnts.begin(), all_pnts.end(), (fastf_t *)*pnts);
    if (hit_limit)
	return BREP_CDT_FAST_LIMIT;
    return summary.failed_faces ? BREP_CDT_FAST_PARTIAL : BREP_CDT_FAST_OK;
}

int
brep_cdt_fast(int **faces, int *face_cnt, vect_t **pnt_norms,
	point_t **pnts, int *pntcnt, const ON_Brep *brep, int index,
	const struct bg_tess_tol *ttol, const struct bn_tol *tol)
{
    struct brep_cdt_fast_options options;
    struct brep_cdt_fast_report report;
    brep_cdt_fast_options_default(&options);
    options.allow_partial = 1;
    int ret = brep_cdt_fast_ex(faces, face_cnt, pnt_norms, pnts, pntcnt,
	brep, index, ttol, tol, &options, &report);
    return (ret == BREP_CDT_FAST_OK || ret == BREP_CDT_FAST_PARTIAL ||
	(ret == BREP_CDT_FAST_LIMIT && *face_cnt > 0)) ?
	BRLCAD_OK : BRLCAD_ERROR;
}

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
