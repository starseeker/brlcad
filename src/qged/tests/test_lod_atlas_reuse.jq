include "lod_atlas_reuse";

{lod_gpu_atlas_suffix_upload_bytes: 100,
 lod_gpu_atlas_lineage_reuses: 4,
 lod_gpu_ordinary_lineage_reuses: 2,
 lod_gpu_ordinary_part_buffer_bytes: 200,
 lod_gpu_ordinary_full_upload_bytes: 300,
 lod_gpu_ordinary_lineage_replacements: 0} as $before |
[
  {name: "unchanged", after: $before, expected: false},
  {name: "same_generation_suffix", expected: true,
   after: ($before | .lod_gpu_atlas_suffix_upload_bytes += 100)},
  {name: "atlas_generation_reuse", expected: true,
   after: ($before | .lod_gpu_atlas_lineage_reuses += 1)},
  {name: "ordinary_generation_reuse", expected: true,
   after: ($before | .lod_gpu_ordinary_lineage_reuses += 1)},
  {name: "ordinary_first_upload", expected: true,
   after: ($before | .lod_gpu_ordinary_part_buffer_bytes += 100 |
     .lod_gpu_ordinary_full_upload_bytes += 100)},
  {name: "ordinary_replacement", expected: false,
   after: ($before | .lod_gpu_ordinary_part_buffer_bytes += 100 |
     .lod_gpu_ordinary_full_upload_bytes += 100 |
     .lod_gpu_ordinary_lineage_replacements += 1)},
  {name: "atlas_full_upload_only", expected: false,
   after: ($before | .lod_gpu_atlas_full_upload_bytes = 1000)},
  {name: "ordinary_reservation_only", expected: false,
   after: ($before | .lod_gpu_ordinary_part_buffer_bytes += 100)},
  {name: "suffix_counter_rewound", expected: false,
   after: ($before | .lod_gpu_atlas_suffix_upload_bytes = 0)},
  {name: "missing_baseline", before: {}, expected: false,
   after: ($before | .lod_gpu_atlas_suffix_upload_bytes += 100)},
  {name: "missing_sample", after: {}, expected: false},
  {name: "string_counter", expected: false,
   after: ($before | .lod_gpu_atlas_suffix_upload_bytes = "200")},
  {name: "fractional_counter", expected: false,
   after: ($before | .lod_gpu_atlas_suffix_upload_bytes = 200.5)}
] |
map(. + {passed: (lod_atlas_append_or_reuse(.before // $before; .after) == .expected)}) |
{contract: "atlas_incremental_realization", cases: length,
 failures: map(select(.passed | not))},
all(.[]; .passed)
