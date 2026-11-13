#!/usr/bin/env bash

set -u

usage()
{
    echo "Usage: $0 <build-dir> <output-dir> [repeat] [quick|full] [data-build-dir]" >&2
}

if [ "$#" -lt 2 ] || [ "$#" -gt 5 ]; then
    usage
    exit 1
fi

build_dir=$1
output_dir=$2
repeat=${3:-3}
tier=${4:-quick}
data_build_dir=${5:-$build_dir}
perf_bin="$build_dir/bin/ged_draw_perf"

if [ ! -x "$perf_bin" ]; then
    echo "obol_draw_perf: missing executable $perf_bin" >&2
    exit 1
fi
if [ "$tier" != "quick" ] && [ "$tier" != "full" ]; then
    usage
    exit 1
fi
case "$repeat" in
    ''|*[!0-9]*|0)
	usage
	exit 1
	;;
esac

mkdir -p "$output_dir"
summary="$output_dir/results.txt"
failures=0
cases_run=0
: > "$summary"

run_case()
{
    label=$1
    database=$2
    root=$3
    mode=$4
    log="$output_dir/$label.log"
    cache="$output_dir/cache-$label"

    if [ ! -f "$database" ]; then
	echo "skip=$label reason=missing_database path=$database" | tee -a "$summary"
	return
    fi

    cases_run=$((cases_run + 1))
    echo "case=$label database=$database root=$root mode=${mode:-wire}" | tee -a "$summary"
    if [ -n "$mode" ]; then
	"$perf_bin" --profile --autoview --render --repeat "$repeat" \
	    --cache-dir "$cache" --clear-cache -- "$database" "$mode" "$root" \
	    > "$log" 2>&1
    else
	"$perf_bin" --profile --autoview --render --repeat "$repeat" \
	    --cache-dir "$cache" --clear-cache -- "$database" "$root" \
	    > "$log" 2>&1
    fi
    status=$?
    iteration_count=$(grep -c '^iter=' "$log")
    grep '^iter=' "$log" | tee -a "$summary"
    if [ "$status" -ne 0 ]; then
	echo "case=$label status=failed log=$log" | tee -a "$summary"
	failures=$((failures + 1))
	return
    fi
    if [ "$iteration_count" -lt "$repeat" ]; then
	echo "case=$label status=failed reason=incomplete_profile expected=$repeat actual=$iteration_count log=$log" | tee -a "$summary"
	failures=$((failures + 1))
	return
    fi
    echo "case=$label status=ok log=$log" | tee -a "$summary"
}

run_case moss_wire "$data_build_dir/share/db/moss.g" all.g ""
run_case moss_evaluated "$data_build_dir/share/db/moss.g" all.g -m3
run_case moss_shaded "$data_build_dir/share/db/moss.g" all.g -m2
run_case moss_hidden "$data_build_dir/share/db/moss.g" all.g -m4

run_case m35_wire "$data_build_dir/share/db/m35.g" all.g ""
run_case havoc_wire "$data_build_dir/share/db/havoc.g" havoc ""
run_case generic_twin_wire "$data_build_dir/share/db/faa/Generic_Twin.g" all ""

if [ "$tier" = "full" ]; then
    run_case m35_evaluated "$data_build_dir/share/db/m35.g" all.g -m3
    run_case havoc_evaluated "$data_build_dir/share/db/havoc.g" havoc -m3
    run_case generic_twin_evaluated \
	"$data_build_dir/share/db/faa/Generic_Twin.g" all -m3
    run_case m35_shaded "$data_build_dir/share/db/m35.g" all.g -m2
    run_case havoc_shaded "$data_build_dir/share/db/havoc.g" havoc -m2
    run_case generic_twin_shaded \
	"$data_build_dir/share/db/faa/Generic_Twin.g" all -m2
    run_case m35_hidden "$data_build_dir/share/db/m35.g" all.g -m4
    run_case havoc_hidden "$data_build_dir/share/db/havoc.g" havoc -m4
    run_case generic_twin_hidden \
	"$data_build_dir/share/db/faa/Generic_Twin.g" all -m4
fi

if [ "$cases_run" -eq 0 ]; then
    echo "status=failed reason=no_performance_inputs" | tee -a "$summary"
    exit 1
fi

if [ "$failures" -ne 0 ]; then
    exit 1
fi
