#!/usr/bin/env bash
#
# Release qualification for the complete installed NIST BRep corpus.  Every
# renderer/LoD pair owns one initially empty cache and one immediate warm
# replay.  Large model files are referenced in place; no database is copied.

set -uo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source_root="$(cd "$script_dir/../../.." && pwd)"
build_dir="${1:-$source_root/.build}"
artifact_dir="${2:-}"
run_timeout="${BOBOL_NIST_MATRIX_TIMEOUT:-180}"

build_dir="$(realpath -m "$build_dir")"
test_bin="$build_dir/bin/test_qtcad_obol_real_models"
if [[ ! -x "$test_bin" ]]; then
    echo "ERROR: missing NIST real-model test executable: $test_bin" >&2
    exit 2
fi
for tool in identify timeout; do
    command -v "$tool" >/dev/null 2>&1 || {
	echo "ERROR: qtcad NIST matrix requires $tool" >&2
	exit 2
    }
done
if [[ -z "$artifact_dir" ]]; then
    artifact_dir="$(mktemp -d /tmp/qtcad-obol-nist.XXXXXX)"
elif [[ -e "$artifact_dir" ]]; then
    echo "ERROR: artifact directory already exists: $artifact_dir" >&2
    exit 2
else
    mkdir -p "$artifact_dir"
fi
artifact_dir="$(realpath -m "$artifact_dir")"
mkdir -p "$artifact_dir"/{caches,runs}
results="$artifact_dir/results.csv"
printf 'status,case,backend,lod,phase,seconds\n' > "$results"

cases=(
    nist_pmi1_shaded
    nist_pmi2_shaded
    nist_pmi3_shaded
    nist_pmi4_shaded
    nist_pmi5_shaded
    nist_pmi6_shaded
    nist_pmi11_shaded
    nist_pmi7_10_shaded
    nist_pmi7_10_wire
)
backends=(system osmesa)
lod_modes=(auto off)
phases=(cold warm)
[[ -n "${BOBOL_NIST_MATRIX_CASES:-}" ]] && read -r -a cases <<<"$BOBOL_NIST_MATRIX_CASES"
[[ -n "${BOBOL_NIST_MATRIX_BACKENDS:-}" ]] && read -r -a backends <<<"$BOBOL_NIST_MATRIX_BACKENDS"
[[ -n "${BOBOL_NIST_MATRIX_LOD_MODES:-}" ]] && read -r -a lod_modes <<<"$BOBOL_NIST_MATRIX_LOD_MODES"
[[ -n "${BOBOL_NIST_MATRIX_PHASES:-}" ]] && read -r -a phases <<<"$BOBOL_NIST_MATRIX_PHASES"
failures=0
passes=0

for case_name in "${cases[@]}"; do
    for backend in "${backends[@]}"; do
	if [[ "$backend" == system && -z "${DISPLAY:-}" ]]; then
	    echo "ERROR: System-GL NIST matrix requires DISPLAY" >&2
	    exit 125
	fi
	expected_renderer="$backend"
	[[ "$backend" == system ]] && expected_renderer="system-gl"
	for lod in "${lod_modes[@]}"; do
	    cache="$artifact_dir/caches/${case_name}-${backend}-${lod}"
	    mkdir -p "$cache"
	    for phase in "${phases[@]}"; do
		run="${case_name}-${backend}-${lod}-${phase}"
		out="$artifact_dir/runs/$run"
		mkdir -p "$out"
		capture="$out/frame.png"
		log="$out/test.log"
		started=$SECONDS
		env_args=(
		    "BRLCAD_ROOT=$build_dir"
		    "BU_DIR_CACHE=$cache"
		    "BOBOL_QTCAD_REAL_MODEL_LOD=$lod"
		    "BOBOL_QTCAD_REAL_MODEL_TIMING=1"
		    "BOBOL_QTCAD_REAL_MODEL_CAPTURE=$capture"
		)
		if [[ "$backend" == system ]]; then
		    env_args+=("BOBOL_QTCAD_REAL_MODEL_GL=1" "QT_QPA_PLATFORM=xcb")
		else
		    env_args+=("BOBOL_QTCAD_REAL_MODEL_GL=0" "QT_QPA_PLATFORM=offscreen"
			"DISPLAY=" "XAUTHORITY=")
		fi
		if timeout --signal=TERM "$run_timeout" env "${env_args[@]}" \
			"$test_bin" "$case_name" >"$log" 2>&1; then
		    terminal="$capture.redraw2.png"
		    [[ "$lod" == auto ]] && terminal="$capture.settled.png"
		    if [[ -s "$capture" && -s "$terminal" ]] &&
			identify "$capture" "$terminal" >>"$log" 2>&1 &&
			grep -q "CONFIG $case_name renderer=$expected_renderer lod=$lod " "$log"; then
			printf 'PASS,%s,%s,%s,%s,%s\n' "$case_name" "$backend" \
			    "$lod" "$phase" "$((SECONDS-started))" >> "$results"
			passes=$((passes + 1))
			continue
		    fi
		fi
		printf 'FAIL,%s,%s,%s,%s,%s\n' "$case_name" "$backend" \
		    "$lod" "$phase" "$((SECONDS-started))" >> "$results"
		failures=$((failures + 1))
	    done
	done
    done
done

if [[ "$failures" -ne 0 ]]; then
    echo "qtcad NIST matrix FAIL passes=$passes failures=$failures artifacts=$artifact_dir" >&2
    cat "$results" >&2
    exit 1
fi

echo "qtcad NIST matrix PASS rows=$passes cases=${#cases[@]} artifacts=$artifact_dir"
if [[ "${BRLCAD_KEEP_TEST_ARTIFACTS:-0}" == 1 ]]; then
    exit 0
fi
rm -rf "$artifact_dir"
exit 0
