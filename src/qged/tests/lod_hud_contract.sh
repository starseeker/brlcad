#!/bin/sh

set -eu

if test "$#" -ne 4; then
    echo "Usage: $0 JQ IDENTIFY CONVERT REPORT" >&2
    exit 2
fi

jq_executable="$1"
identify_executable="$2"
convert_executable="$3"
report="$4"
minimum_fill_pixels=100
idle_phase=0
crop_width=40
crop_top=20
crop_bottom=80

# A checkpoint reads presented pixels, while telemetry describes current
# retained features. An interrupted traversal may retire its render latch while
# a bounded recovery is pending, so exact CAD work and an empty latch do not
# identify the displayed HUD. Authenticate its feature revision against the
# completed framebuffer before applying the captured-geometry visibility floor.
# Interior rows exclude fractional endpoints without assuming GL rounding.
samples=$("$jq_executable" -r --argjson minimum "$minimum_fill_pixels" \
    --argjson idle "$idle_phase" '
    .samples[] |
        select((.checkpoint? // "") != "" and
            ((.presented_cad_faces // 0) > 0 or
                (.presented_cad_lines // 0) > 0) and
            (.lod_convergence_visible_targets // 0) > 0 and
            (.lod_convergence_phase // $idle) > $idle) |
        ((.lod_convergence_terminal_error == true) and
            (.view_lod_policy // 0) != 0 and
            (.view_lod_mesh_enabled == true or
                .view_lod_csg_enabled == true)) as $terminal_error |
        select(.lod_progress_fill_present == true or $terminal_error) |
        if (.render_requested | type) != "boolean" or
            (.cad_aggregate_work_exact | type) != "boolean" then
            error("HUD checkpoint has no valid presentation state")
        else . end |
        select(.render_requested == false and .cad_aggregate_work_exact) |
        if (.feature_presentation_revision | type) != "string" or
            (.feature_presentation_revision | test("^[0-9]+$")) != true or
            (has("presented_feature_revision") | not) or
            (.presented_feature_revision != null and
                ((.presented_feature_revision | type) != "string" or
                 (.presented_feature_revision | test("^[0-9]+$")) != true)) then
            error("HUD checkpoint has no valid feature presentation identity")
        else . end |
        if $terminal_error and
            .presented_feature_revision != .feature_presentation_revision then
            error("Terminal error HUD has not reached the completed framebuffer")
        else . end |
        select(.presented_feature_revision == .feature_presentation_revision) |
        if $terminal_error and
            (.lod_progress_fill_present != true or
             .lod_progress_label_present != true or
             (.lod_progress_label_text | startswith("View incomplete")) != true or
             .lod_progress_fill_color != [255, 90, 80] or
             .lod_progress_track_present != true or
             .lod_progress_track_minimum_y != .lod_progress_track_maximum_y or
             .lod_progress_fill_maximum_y != .lod_progress_track_maximum_y) then
            error("Terminal error HUD has no complete red fill and error label")
        else . end |
        if (.lod_progress_fill_minimum_y | type) != "number" or
            (.lod_progress_fill_maximum_y | type) != "number" or
            (.lod_progress_fill_line_width | type) != "number" or
            .lod_progress_fill_maximum_y < .lod_progress_fill_minimum_y or
            .lod_progress_fill_line_width <= 0 then
            error("HUD fill has invalid geometry")
        else . end |
        ((.lod_progress_fill_maximum_y | floor) -
            (.lod_progress_fill_minimum_y | ceil)) as $interior_rows |
        select($interior_rows * (.lod_progress_fill_line_width | floor) >=
            $minimum) |
        .checkpoint
' "$report")

if test -z "$samples"; then
    echo "SKIP: no current HUD checkpoint with $minimum_fill_pixels interior fill pixels"
    exit 0
fi

printf '%s\n' "$samples" | while IFS= read -r sample; do
    if ! test -f "$sample"; then
        echo "LoD HUD checkpoint is missing: $sample" >&2
        exit 1
    fi

    dimensions=$("$identify_executable" -format '%wx%h' "$sample")
    width="${dimensions%x*}"
    height="${dimensions#*x}"
    case "$width:$height" in
        *[!0-9:]* | :* | *:)
            echo "LoD HUD checkpoint has invalid dimensions: $dimensions" >&2
            exit 1
            ;;
    esac
    if test "$width" -le "$crop_width" ||
        test "$height" -le "$((crop_top + crop_bottom))"; then
        echo "LoD HUD checkpoint is too small: $dimensions" >&2
        exit 1
    fi

    # Keep the established visibility floor and phase-color tolerance.  A missing
    # fill behind its gray track must fail even when diagnostic progress advances.
    crop_height=$((height - crop_top - crop_bottom))
    fill_pixels=$("$convert_executable" "$sample" \
        -crop "${crop_width}x${crop_height}+$((width - crop_width))+$crop_top" \
        -fuzz 20% -fill white -opaque '#60dcff' \
        -opaque '#70eb87' -opaque '#ffcd48' -opaque '#ffaa40' \
        -opaque '#ff5a50' \
        -fill black +opaque white -format '%[fx:mean*w*h]' info:)
    if ! awk -v pixels="$fill_pixels" -v minimum="$minimum_fill_pixels" \
        'BEGIN { exit !(pixels ~ /^[0-9]+([.][0-9]+)?$/ && pixels + 0 >= minimum) }'; then
        echo "LoD convergence HUD track has no visible progress fill ($fill_pixels pixels): $sample" >&2
        exit 1
    fi
    echo "LoD HUD visible fill: $fill_pixels pixels in $sample"
done
