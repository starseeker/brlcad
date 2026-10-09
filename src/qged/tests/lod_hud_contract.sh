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
minimum_accent_pixels=40
idle_phase=0
crop_width=450
crop_height=42

# A checkpoint reads presented pixels, while telemetry describes current
# retained features.  Authenticate the feature revision against the completed
# framebuffer before checking the card pixels.
samples=$("$jq_executable" -r --argjson idle "$idle_phase" '
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
        select(.lod_progress_card_present == true or $terminal_error) |
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
            (.lod_progress_card_present != true or
             (.lod_progress_card_title |
                startswith("Geometry preparation failed")) != true or
             .lod_progress_card_color != [255, 89, 79] or
             .lod_progress_card_terminal != true or
             .lod_progress_card_ready != false) then
            error("Terminal error HUD has no red failure card")
        else . end |
        .checkpoint
' "$report")

if test -z "$samples"; then
    echo "SKIP: no current LoD card checkpoint"
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
    if test "$width" -le "$crop_width" || test "$height" -le "$crop_height"; then
        echo "LoD HUD checkpoint is too small: $dimensions" >&2
        exit 1
    fi

    accent_pixels=$("$convert_executable" "$sample" \
        -crop "${crop_width}x${crop_height}+$((width - crop_width))+$((height - crop_height))" \
        -fuzz 20% -fill white -opaque '#61bfff' \
        -opaque '#70eb87' -opaque '#ffbf47' -opaque '#ff594f' \
        -fill black +opaque white -format '%[fx:mean*w*h]' info:)
    if ! awk -v pixels="$accent_pixels" -v minimum="$minimum_accent_pixels" \
        'BEGIN { exit !(pixels ~ /^[0-9]+([.][0-9]+)?$/ && pixels + 0 >= minimum) }'; then
        echo "LoD status card has no visible accent ($accent_pixels pixels): $sample" >&2
        exit 1
    fi
    echo "LoD HUD visible accent: $accent_pixels pixels in $sample"
done
