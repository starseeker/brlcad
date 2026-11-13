# Run with jq -e --arg mode shaded -f lod_lucy_contract.jq report.json.
# Emit every failed condition before the exit-status Boolean so an observation
# defect cannot hide a separate quality failure in the same replay.
def checkpoint($name):
  first(.samples[] | select((.checkpoint? // "") |
    endswith("/" + $name + ".png"))) // null;

def check($name; $sample; $passed):
  {condition: $name, passed: ($passed == true),
   event_index: $sample.event_index, checkpoint: $sample.checkpoint};

def mesh_presented:
  .active_lod_cad_payloads > 0 and
  .active_progressive_cad_faces > 0 and
  .presented_cad_work_exact == true and
  (if $mode == "shaded" then .presented_cad_faces > 0
   else .presented_cad_lines > 0 end);

def no_proxy:
  .active_cad_subpixel_proxy_points == 0 and
  .active_cad_subpixel_proxy_draw_points == 0;

def quality_checkpoints:
  ["zoom-in-stable", "smooth-zoom-close-stable", "smooth-zoom-out-stable",
   "smooth-zoom-return", "background-cache-complete", "final-stable"];

. as $report |
checkpoint("first-coverage-ready") as $coverage |
checkpoint("first-cad-mesh-ready") as $mesh |
checkpoint("background-cache-complete") as $stable |
.samples[-1] as $final |
# Draw return can precede both autoview and the first affordable mesh.  The
# bounded first-useful-image and cold/warm coverage deadlines remain separate
# matrix gates.  Once any mesh is observed, even an unsettled exact frame may
# not regress to a point.  Readiness independently closes the startup window.
([.samples[] | select(mesh_presented) | .event_index] | min) as $mesh_event |
[
  check("mode"; null; $mode == "shaded" or $mode == "wire"),
  check("coverage.checkpoint"; $coverage; $coverage != null),
  check("coverage.payload"; $coverage; $coverage.active_lod_cad_payloads > 0),
  # Coverage publication can precede System GL frame certification.  The
  # separate first-mesh checkpoint and first-useful image gate prove drawing.
  check("coverage.geometry"; $coverage;
    $coverage.active_progressive_cad_faces > 0 or
    $coverage.presented_cad_faces > 0 or $coverage.presented_cad_lines > 0),
  check("coverage.structural_boxes"; $coverage;
    $coverage.visible_structural_fallback_boxes == 0),
  check("first_mesh.checkpoint"; $mesh; $mesh != null),
  check("first_mesh.presentation"; $mesh; $mesh | mesh_presented),
  check("registry.availability"; $final;
    ($final.deep_lod_diagnostics | type) == "boolean"),
  # Compact storage is optional observation, not the semantic definition of
  # drawing a bare primitive or a combination root.  Enforce it when measured.
  (if $final.deep_lod_diagnostics == true then
    check("registry.entries"; $final; $final.compact_lod_entries > 0),
    check("registry.payload_entries"; $final;
      $final.compact_lod_entries_with_payload > 0 and
      $final.compact_lod_entries_with_payload <= $final.compact_lod_entries and
      $final.compact_lod_entries_with_payload == $final.active_lod_cad_payloads)
   else empty end),
  check("final.presentation"; $final; $final | mesh_presented),
  check("final.resident_assets"; $final; $final.lod_service_resident_assets > 0),
  check("final.resident_bytes"; $final; $final.lod_service_resident_bytes > 0),
  (.samples[] | select(.presented_cad_work_exact == true) |
    select(.lod_convergence_view_ready == true or
      ($mesh_event != null and .event_index >= $mesh_event) or
      ($coverage != null and .event_index >= $coverage.event_index)) |
    check("presentation.no_proxy_after_coverage"; .; no_proxy)),
  check("settled.checkpoint"; $stable; $stable != null),
  (if $mode == "shaded" then
    check("settled.prominent_mesh"; $stable; $stable.lod_prominent_cad_payloads > 0),
    check("settled.pixel_target"; $stable;
      $stable.lod_convergence_performance_limited == true or
      (($stable.lod_max_cad_normalized_error | type) == "number" and
       # Existing normalized-error limit includes discrete-cut hysteresis.
       $stable.lod_max_cad_normalized_error <= 1.251))
   else empty end),
  # These quiet camera endpoints are mandatory, even if a malformed report
  # omits them or incorrectly labels them as still moving.  Capacity limits
  # cannot excuse prominent-floor debt at a settled Lucy view.
  (quality_checkpoints[] |
    . as $name | ($report | checkpoint($name)) as $sample |
    check("quality.checkpoint." + $name; $sample; $sample != null)),
  (.samples[] | select(.checkpoint != null) |
    select(.lod_convergence_view_ready == true or
      ((.checkpoint | split("/") | last | rtrimstr(".png")) as $name |
        quality_checkpoints | index($name) != null)) |
    select($mode == "shaded" or (.checkpoint | contains("/smooth-zoom-"))) |
    check("quality.prominent_floor"; .;
      .lod_prominent_cad_quality_floor_violations == 0) +
      {violations: .lod_prominent_cad_quality_floor_violations,
       normalized_error: .lod_max_cad_normalized_error,
       performance_limited: .lod_convergence_performance_limited})
] as $checks |
{contract: "lucy_coverage_presentation",
 registry_diagnostics: (if $final.deep_lod_diagnostics == true then "measured"
   elif $final.deep_lod_diagnostics == false then "not-collected" else "missing" end),
 failures: [$checks[] | select(.passed | not)]},
all($checks[]; .passed)
