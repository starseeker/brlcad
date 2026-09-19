# A cut can expose more of one immutable generation.  Its suffix upload is
# incremental realization even when the generation-swap counter stays fixed.
def lod_gpu_counter:
  if type == "number" then . >= 0 and . == floor else false end;

def lod_gpu_counter_advanced($before; $after):
  ($before | lod_gpu_counter) and ($after | lod_gpu_counter) and
  $after > $before;

def lod_atlas_append_or_reuse($before; $after):
  lod_gpu_counter_advanced($before.lod_gpu_atlas_suffix_upload_bytes;
    $after.lod_gpu_atlas_suffix_upload_bytes) or
  lod_gpu_counter_advanced($before.lod_gpu_atlas_lineage_reuses;
    $after.lod_gpu_atlas_lineage_reuses) or
  lod_gpu_counter_advanced($before.lod_gpu_ordinary_lineage_reuses;
    $after.lod_gpu_ordinary_lineage_reuses) or
  # A newly exposed ordinary part has no preceding generation to reuse.
  (lod_gpu_counter_advanced($before.lod_gpu_ordinary_part_buffer_bytes;
     $after.lod_gpu_ordinary_part_buffer_bytes) and
   lod_gpu_counter_advanced($before.lod_gpu_ordinary_full_upload_bytes;
     $after.lod_gpu_ordinary_full_upload_bytes) and
   ($before.lod_gpu_ordinary_lineage_replacements | lod_gpu_counter) and
   ($after.lod_gpu_ordinary_lineage_replacements | lod_gpu_counter) and
   $after.lod_gpu_ordinary_lineage_replacements ==
     $before.lod_gpu_ordinary_lineage_replacements);
