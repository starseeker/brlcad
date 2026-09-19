# Completed frames and deadline interruptions are separate observations. A fast
# recovery frame cannot erase the missed deadline which selected its ceiling.
def checkpoint($name):
  first(.samples[] | select((.checkpoint? // "") |
    endswith("/" + $name + ".png"))) // null;

def nonnegative_number:
  type == "number" and . >= 0;

def counter:
  nonnegative_number and floor == .;

def check($name; $passed):
  {condition: $name, passed: ($passed == true)};

1.05 as $completed_deadline_slack |
1.01 as $responsive_pixel_limit |
0.95 as $retained_population_fraction |
checkpoint("zoom-return-stable") as $before |
checkpoint("zoom-out-motion") as $zoom_out |
checkpoint("rotate-held-end") as $held |
checkpoint("rotate-motion") as $motion |
(first(.samples[] | select(.action == "wheel")) // null) as $zoom_in |
(first(.samples[] | select(.action == "mouse_press")) // null) as $press |
(first(.samples[] | select(.action == "mouse_release" and
  .event_index > $press.event_index)) // null) as $release |
# Evidence after the held checkpoint cannot justify a ceiling already visible
# there. The release sample belongs to the later quiet handoff.
[.samples[] | select(.event_index >= $press.event_index and
  .event_index <= $held.event_index)] as $rotation |
($before.lod_interactive_target_fps // 0) as $fps |
(if ($fps | type) == "number" and $fps > 0 then
  1000.0 * $completed_deadline_slack / $fps else null end) as $soft_deadline |
($rotation | map(.last_render_ms) | max) as $completed_peak |
[$rotation | to_entries[] | select(.key > 0) |
  .value as $current | $rotation[.key - 1] as $prior |
  select($current.presentation_interrupted_frames >
    $prior.presentation_interrupted_frames) |
  {event_index: $current.event_index,
   elapsed_ms: $current.presentation_last_interrupted_ms,
   deadline_ms: $current.presentation_deadline_interactive_ms}] as $interruptions |
($held.lod_target_pixel_error > $responsive_pixel_limit or
  $held.lod_interactive_progressive_ceiling >= 0) as $corrected |
(if $mode == "wire" then "presented_cad_lines"
 else "presented_cad_faces" end) as $population |
[
  check("mode"; $mode == "shaded" or $mode == "wire"),
  check("checkpoints"; all([$before, $zoom_out, $held, $motion,
    $zoom_in, $press, $release][]; . != null)),
  check("event_order"; $before.event_index < $press.event_index and
    $press.event_index < $held.event_index and
    $held.event_index < $release.event_index),
  check("zoom_scale"; $zoom_in.lod_scale_changing_interaction == true and
    $zoom_out.lod_scale_changing_interaction == true),
  check("pose_scale"; $held.lod_scale_changing_interaction == false and
    $motion.lod_scale_changing_interaction == false),
  (if .backend == "system_gl" then
    check("timing"; $soft_deadline != null and
      ($before.last_render_ms | nonnegative_number) and
      ($rotation | length) > 1 and
      all($rotation[]; .last_render_ms | nonnegative_number)),
    check("interruption_counters";
      all($rotation[]; .presentation_interrupted_frames | counter) and
      all(range(1; $rotation | length);
        $rotation[.].presentation_interrupted_frames >=
        $rotation[. - 1].presentation_interrupted_frames)),
    check("interruption_deadlines";
      all($interruptions[];
        (.elapsed_ms | nonnegative_number) and
        (.deadline_ms | nonnegative_number) and .deadline_ms > 0 and
        .elapsed_ms >= .deadline_ms)),
    (if $soft_deadline != null and $before.last_render_ms <= $soft_deadline then
      if $completed_peak > $soft_deadline or ($interruptions | length) > 0 then
        check("pressure_response"; $corrected)
      else
        check("responsive_pixel_error";
          ($held.lod_target_pixel_error | nonnegative_number) and
          $held.lod_target_pixel_error <= $responsive_pixel_limit),
        check("responsive_ceiling"; $held.lod_interactive_progressive_ceiling == -1),
        check("responsive_population";
          $before.presented_cad_work_exact == true and
          $held.presented_cad_work_exact == true and
          ($before[$population] | nonnegative_number) and $before[$population] > 0 and
          ($held[$population] | nonnegative_number) and
          $held[$population] >= $before[$population] * $retained_population_fraction)
      end
     else empty end)
   else empty end)
] as $checks |
{contract: "zoom_pose_retained_cut", completed_peak_ms: $completed_peak,
 completed_deadline_ms: $soft_deadline, interruptions: $interruptions,
 failures: [$checks[] | select(.passed | not)]},
all($checks[]; .passed)
