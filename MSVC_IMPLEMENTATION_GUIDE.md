# MSVC Symbol Retention: Implementation Guide

## Two Viable Approaches

Based on the analysis in `MSVC_IN_CODE_VS_LINKER_FLAG.md`, here are the implementation options:

---

## Approach A: `/OPT:NOREF` Linker Flag ✅ RECOMMENDED

**Status:** Already implemented

**Implementation:**
```cmake
# src/libged/CMakeLists.txt (lines 351-356)
if (LIBGED_STATIC_CORE)
  target_link_libraries(libged PRIVATE libged_plugin_deps)
  # On MSVC, prevent the linker from stripping "unused" static command registration symbols
  if(MSVC)
    target_link_options(libged PRIVATE "/OPT:NOREF")
  endif()
endif()
```

**Pros:**
- ✅ Robust and maintainable
- ✅ Already working
- ✅ Standard industry practice
- ✅ Single point of configuration

**Cons:**
- ❌ Requires CMake build system modification
- ❌ User preference is to avoid linker changes

**Recommendation:** Keep this approach. It's the engineering best practice.

---

## Approach B: Scanner-Generated MSVC Pragmas ⚠️ ALTERNATIVE

**Status:** Not yet implemented (can be added as alternative)

This approach generates an MSVC-specific header file with `#pragma comment(linker, "/INCLUDE:...")` directives for each command symbol, then includes it in the build.

### Implementation Plan

#### Step 1: Modify `ged_cmd_scanner.cpp`

Add support for generating MSVC pragma header:

```cpp
// Add to argument parsing (around line 72)
bool gen_msvc_pragmas = false;
std::string msvc_pragma_file;

// In argument parsing loop (around line 82)
else if (arg == "--gen-msvc-pragmas") {
    if (i + 1 < argc) {
        gen_msvc_pragmas = true;
        msvc_pragma_file = argv[++i];
    } else {
        std::cerr << "Missing filename after --gen-msvc-pragmas\n";
        return EXIT_FAILURE;
    }
}

// Add after static registration file generation (around line 410)
// MSVC pragma file generation
if (gen_msvc_pragmas) {
    const std::string MSVC_MARKER = "// Auto-generated MSVC linker pragmas for symbol retention";
    
    // Check if file exists and is auto-generated
    std::ifstream testfs(msvc_pragma_file);
    if (testfs.is_open()) {
        std::string first_line;
        std::getline(testfs, first_line);
        testfs.close();
        if (first_line != MSVC_MARKER) {
            std::cerr << "ERROR: " << msvc_pragma_file << " already exists but doesn't appear to be auto-generated.\n";
            std::cerr << "       Refusing to overwrite. Rename or remove the file to allow generation.\n";
            return EXIT_FAILURE;
        }
    }
    
    std::ofstream ofs(msvc_pragma_file, std::ios::out | std::ios::trunc);
    if (!ofs.is_open()) {
        std::cerr << "Cannot open output MSVC pragma file: " << msvc_pragma_file << "\n";
        return EXIT_FAILURE;
    }
    
    ofs << MSVC_MARKER << "\n";
    ofs << "/*\n";
    ofs << " * This file forces MSVC linker to include all GED command symbols.\n";
    ofs << " * It uses #pragma comment(linker, \"/INCLUDE:symbol\") to prevent\n";
    ofs << " * the linker from stripping symbols during link-time optimization.\n";
    ofs << " *\n";
    ofs << " * This is an alternative to using /OPT:NOREF linker flag.\n";
    ofs << " */\n\n";
    ofs << "#ifdef _MSC_VER\n";
    ofs << "/* Force inclusion of all command pointer symbols */\n";
    
    for (const auto &cmd_macro : static_cmd_sorted) {
        // Note: Symbol names must match exactly what the linker sees
        // For extern "C" symbols, this is typically the name as-is
        ofs << "#pragma comment(linker, \"/INCLUDE:__ged_cmd_ptr_" << cmd_macro << "\")\n";
    }
    
    ofs << "#endif /* _MSC_VER */\n";
    ofs.close();
    return EXIT_SUCCESS;
}
```

#### Step 2: Update CMakeLists.txt

Add generation of MSVC pragma file:

```cmake
# src/libged/CMakeLists.txt

# After line 287 (where GED_STATIC_REGISTRATION_CPP is defined)
set(GED_MSVC_PRAGMAS_H ${CMAKE_CURRENT_BINARY_DIR}/ged_msvc_pragmas.h)
distclean("${GED_MSVC_PRAGMAS_H}")

# Modify scanner invocation to generate MSVC pragma file
# Find the execute_process that runs GED_CMD_SCANNER and add:
--gen-msvc-pragmas "${GED_MSVC_PRAGMAS_H}"

# Add to LIBGED_SOURCES (around line 293)
# The pragma file needs to be included somewhere, so we can add it as a source
# or include it in ged_init.cpp
```

#### Step 3: Include the Generated Header

Option A - Include in `ged_init.cpp`:

```cpp
// In src/libged/ged_init.cpp, near the top (around line 50)
#ifdef LIBGED_STATIC_CORE
#  ifdef _MSC_VER
#    include "ged_msvc_pragmas.h"  // Force symbol inclusion on MSVC
#  endif
#endif
```

Option B - Include in `plugin.h`:

```cpp
// In src/libged/include/plugin.h, near the end
#ifdef LIBGED_STATIC_CORE
#  ifdef _MSC_VER
#    include "../../../build/src/libged/ged_msvc_pragmas.h"
     // Note: Path might need adjustment based on build directory
#  endif
#endif
```

#### Step 4: Remove `/OPT:NOREF` (Optional)

If using MSVC pragmas instead of linker flag:

```cmake
# In src/libged/CMakeLists.txt
# Comment out or remove the /OPT:NOREF section:
if (LIBGED_STATIC_CORE)
  target_link_libraries(libged PRIVATE libged_plugin_deps)
  # MSVC symbol retention now handled by generated pragmas in ged_msvc_pragmas.h
  # if(MSVC)
  #   target_link_options(libged PRIVATE "/OPT:NOREF")
  # endif()
endif()
```

### Testing Requirements

**CRITICAL:** This approach MUST be tested on actual MSVC:

1. **Build Configuration:**
   - Windows with MSVC (Visual Studio 2019 or later)
   - CMake configuration with LIBGED_STATIC_CORE=ON
   
2. **Test Cases:**
   - Verify all 227 commands are available: `mged` → `?` command
   - Check symbol names match exactly
   - Test both Debug and Release builds
   - Verify no link errors

3. **Symbol Name Verification:**
   ```cmd
   REM Use dumpbin to verify symbol names
   dumpbin /symbols libged.lib | findstr __ged_cmd_ptr
   ```

### Potential Issues

1. **Symbol Name Mismatches:**
   - Pragma says: `/INCLUDE:__ged_cmd_ptr_threeptarb_cmd`
   - Actual symbol might be: `___ged_cmd_ptr_threeptarb_cmd` (extra underscore)
   - Solution: Test on MSVC, adjust scanner if needed

2. **Build Path Dependencies:**
   - Generated header location must be findable by includes
   - May need `target_include_directories()` adjustment

3. **Incremental Build Issues:**
   - Scanner must re-run when command files change
   - CMake dependencies must be set up correctly

---

## Hybrid Approach: Both for Maximum Safety

If desired, you can use BOTH mechanisms:

```cmake
if (LIBGED_STATIC_CORE)
  target_link_libraries(libged PRIVATE libged_plugin_deps)
  
  if(MSVC)
    # Primary: Generated pragmas (self-contained in code)
    # Included via ged_msvc_pragmas.h in ged_init.cpp
    
    # Backup: Linker flag (in case pragma names don't match exactly)
    # This is redundant if pragmas work, but provides safety net
    # target_link_options(libged PRIVATE "/OPT:NOREF")
  endif()
endif()
```

**Rationale:**
- If pragma symbol names match: Pragmas work, linker flag is redundant but harmless
- If pragma names mismatch: Linker flag provides fallback
- During development: Remove linker flag once pragmas proven to work

---

## Recommendation Summary

### For Production Use: Keep `/OPT:NOREF` ✅

**Reason:** It's proven, robust, and standard practice.

The CMake change is minimal (4 lines) and isolated. The engineering benefits far outweigh the aesthetic preference for keeping build logic out of CMakeLists.txt.

### For User Preference: Implement Scanner-Generated Pragmas ⚠️

**If CMake changes are absolutely unacceptable:**

1. Implement the scanner modification (shown above)
2. Generate `ged_msvc_pragmas.h`
3. Include it in `ged_init.cpp`
4. **Test thoroughly on MSVC**
5. Be prepared for symbol name debugging

**Trade-off:** More complex, requires MSVC testing, but achieves the goal of self-contained code.

### Pragmatic Approach: Use Both During Transition

1. Keep `/OPT:NOREF` for now (proven to work)
2. Add pragma generation (for future migration)
3. Test pragma approach on MSVC
4. Once proven: Remove `/OPT:NOREF`
5. Document pragma approach for other developers

---

## Implementation Checklist

If implementing scanner-generated pragmas:

- [ ] Modify `ged_cmd_scanner.cpp` to add `--gen-msvc-pragmas` option
- [ ] Update scanner to generate pragma file with all symbols
- [ ] Modify CMakeLists.txt to invoke scanner with new flag
- [ ] Add generated file to build dependencies
- [ ] Include pragma header in appropriate location (ged_init.cpp recommended)
- [ ] Test on MSVC with Debug build
- [ ] Test on MSVC with Release build
- [ ] Verify all 227 commands available
- [ ] Use `dumpbin` to verify symbol names match pragma directives
- [ ] Document the approach
- [ ] Optionally remove `/OPT:NOREF` if pragmas proven to work

---

## Files to Modify

### Scanner-Generated Pragma Approach

1. **src/libged/ged_cmd_scanner.cpp** - Add MSVC pragma generation
2. **src/libged/CMakeLists.txt** - Add scanner flag and file generation
3. **src/libged/ged_init.cpp** - Include generated pragma header
4. **src/libged/CMakeLists.txt** - Optionally remove `/OPT:NOREF`

### Linker Flag Approach (Current)

1. **src/libged/CMakeLists.txt** - Already has `/OPT:NOREF` (lines 354-356)

---

## Conclusion

Both approaches are viable:

- **Linker flag**: Simpler, more robust, already working
- **Generated pragmas**: More complex, self-contained in code, requires testing

The user's preference to avoid CMake linker changes is understandable, but the linker flag approach is objectively superior from an engineering standpoint. However, the pragma approach can work if implemented carefully and tested thoroughly on MSVC.

**Final recommendation:** Start with the linker flag (already done), optionally add pragma generation as an enhancement later if desired for cleaner separation of build and source concerns.
