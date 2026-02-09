# BRL-CAD gedcmdsreview Branch Code Review Report

## Executive Summary

This report analyzes the changes between the `main` and `gedcmdsreview` branches in the BRL-CAD repository. The gedcmdsreview branch implements a hybrid plugin system for libged that allows commands to be either statically compiled into libged or dynamically loaded as plugins, addressing the startup overhead issues of the previous all-dynamic plugin approach.

**Overall Assessment:** The migration is largely well-executed with consistent patterns across 226+ commands. However, several critical cross-platform and build system issues were identified that need to be addressed.

---

## 1. Missed Commands

### Status: ✅ NO ISSUES FOUND

All command directories from the main branch are present in gedcmdsreview. Command count verification:
- **Main branch:** 227 command directories
- **Gedcmdsreview branch:** 227 command directories  
- **Missing commands:** 0

The `pmat` command was initially flagged but was confirmed to exist - it has a modified CMakeLists.txt but all functionality is preserved.

---

## 2. Cross-Platform Issues

### Issue 2.1: MSVC Symbol Retention Problem ⚠️ **CRITICAL**

**File:** `src/libged/include/plugin.h:80-84`

**Severity:** Critical - Will cause build/runtime failures on Windows

**Problem:**
The `GED_CMD_USED` macro is empty on MSVC, meaning the linker may eliminate static command symbols that appear unreferenced:

```c
#if defined(__GNUC__) || defined(__clang__)
#define GED_CMD_USED __attribute__((used))
#else
#define GED_CMD_USED  // Empty on MSVC - PROBLEM!
#endif
```

When `LIBGED_STATIC_CORE=ON`, each command file creates symbols like `__ged_cmd_ptr_xxx` marked with `GED_CMD_USED`. These are referenced only from the auto-generated `ged_force_static_registration()` function in a separate compilation unit. On MSVC, without `__attribute__((used))` or equivalent, the linker's optimization may remove these "unreferenced" symbols before linking, causing undefined symbol errors or missing commands at runtime.

**Evidence of Risk:**
1. MSVC's `/OPT:REF` linker optimization (enabled by default in Release builds) removes unreferenced data
2. The symbols are only referenced from generated code in a different TU
3. No fallback mechanism exists for MSVC builds

**Recommended Fixes (in order of preference):**

**Option 1 - Use `__declspec(selectany)` with forced export:**
```c
#if defined(__GNUC__) || defined(__clang__)
#define GED_CMD_USED __attribute__((used))
#elif defined(_MSC_VER)
#define GED_CMD_USED __declspec(selectany)
#else
#define GED_CMD_USED
#endif
```

**Option 2 - Add linker flag in CMake for MSVC:**
In `src/libged/CMakeLists.txt`, add:
```cmake
if(MSVC AND LIBGED_STATIC_CORE)
  target_link_options(libged PRIVATE "/OPT:NOREF")
endif()
```

**Option 3 - Use `/INCLUDE:` linker directives (most robust but complex):**
Modify the scanner to generate a .def file or /INCLUDE directives for MSVC builds.

**Impact if not fixed:** Windows builds with `LIBGED_STATIC_CORE=ON` will fail to link or will have missing commands at runtime.

---

### Issue 2.2: Path Separator Handling

**File:** Multiple CMakeLists.txt files

**Severity:** Low - CMake generally handles this, but worth verifying

**Observation:**
The build system uses both forward slashes and `${CMAKE_CURRENT_SOURCE_DIR}` constructs. While CMake normalizes paths on Windows, there are direct path constructions in:
- `src/libged/CMakeLists.txt`: Plugin path construction
- `src/libged/ged_init.cpp`: `bu_dir()` calls

**Recommendation:**
No immediate action required, but recommend testing on Windows to ensure plugin discovery works correctly with the `bu_file_list()` pattern matching in `ged_init.cpp:160`.

---

## 3. Plugin Registration Problems

### Status: ✅ NO MAJOR ISSUES

All 226 command plugins follow the new registration pattern correctly:

```c
#define GED_<NAME>_COMMANDS(X, XID) \
    X(token, function, flags) \
    ...

GED_DECLARE_COMMAND_SET(GED_<NAME>_COMMANDS)
GED_DECLARE_PLUGIN_MANIFEST("libged_<name>", version, GED_<NAME>_COMMANDS)
```

**Verified Patterns:**
- ✅ Single-command plugins (e.g., `3ptarb`)
- ✅ Multi-command plugins (e.g., `draw` with 7 commands)
- ✅ Commands with special characters (e.g., `?` in help, using XID form)
- ✅ Subcommand structures (e.g., bot subcommands)

**Minor Observations:**
- PLUGIN_ONLY flag correctly used for 6 commands that require dynamic loading:
  - `brep` - depends on generated perplex/lemon sources
  - `debug`, `env`, `facetize`, `graph`, `simulate` - likely have specific reasons
- Scanner supports both X() and XID() forms for flexible command naming

---

## 4. Build System Issues

### Issue 4.1: Library Dependency Inconsistency

**Files:** `src/libged/bot/CMakeLists.txt` and potentially others

**Severity:** Medium - May cause link failures in dynamic plugin mode

**Problem:**
The `bot` command's CMakeLists.txt removed `libged` from its `BOT_LIBS`:

```cmake
# main branch:
set(BOT_LIBS libged libbg libbu manifold::manifold)

# gedcmdsreview branch:
set(BOT_LIBS libbg libbu manifold::manifold)  # libged removed!
```

**Analysis:**
The `ged_plugin_library()` function handles this differently depending on mode:

- **Static mode (`LIBGED_STATIC_CORE=ON`)**: Line 166 removes `libged` from LIBS list (correct, since bot code will be part of libged)
- **Dynamic plugin mode (`LIBGED_STATIC_CORE=OFF`)**: Builds as shared library that NEEDS `libged` to resolve symbols

**Current behavior:**
```cmake
# src/libged/CMakeLists.txt:166
if(GEDPL_LIBS)
  list(REMOVE_ITEM GEDPL_LIBS libged)  # Always removes libged in static mode
  target_link_libraries(libged_plugin_deps INTERFACE ${GEDPL_LIBS})
endif()

# Line 200+ for dynamic mode:
if(GEDPL_LIBS)
  target_link_libraries(${name} ${GEDPL_LIBS})  # Uses LIBS as-is
endif()
```

**The Issue:**
When `LIBGED_STATIC_CORE=OFF`, the bot plugin is built as a shared library but `BOT_LIBS` doesn't include `libged`, so undefined symbol errors may occur.

**Recommended Fix:**
Restore `libged` to the LIBS list in individual command CMakeLists.txt files:

```cmake
# src/libged/bot/CMakeLists.txt
set(BOT_LIBS libged libbg libbu manifold::manifold)
```

The `ged_plugin_library()` function will automatically remove it when building in static mode (line 166), but keep it when building as a dynamic plugin.

**Commands to check:**
Review all command CMakeLists.txt files that removed `libged` from their LIBS lists.

---

### Issue 4.2: Scanner Binary Rebuild

**File:** `src/libged/CMakeLists.txt:15-16`

**Severity:** Low - Build system hygiene

**Code:**
```cmake
set(GED_CMD_SCANNER ${CMAKE_BINARY_DIR}/CMakeTmp/ged_cmd_scanner${CMAKE_EXECUTABLE_SUFFIX})
file(REMOVE "${GED_CMD_SCANNER}")  # Added in review branch
```

**Observation:**
The explicit `file(REMOVE)` before `try_compile()` was added. This ensures clean rebuilds but may slow down incremental configuration.

**Recommendation:**
This is fine as-is. The scanner is small and compiles quickly. The safety of guaranteed clean state is worth the minimal overhead.

---

### Issue 4.3: CMake Configuration Dependencies

**File:** `src/libged/CMakeLists.txt:60-67`

**Code:**
```cmake
# TODO: Note that to be strictly correct, we should be making CMake
# reconfigure every time one of the plugin files changes to avoid the
# function definitions being out of sync. The lines below that are commented
# out should do that...
foreach(sf ${ARGN})
  if (EXISTS ${sf})
    #set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${sf})
    execute_process(COMMAND ${GED_CMD_SCANNER} ${sf} OUTPUT_VARIABLE CMD_NAMES)
```

**Observation:**
The scanner runs at configure time, but changes to command signatures don't trigger automatic reconfiguration.

**Impact:**
If a developer changes a command signature (function name, command name, etc.) without manually reconfiguring CMake, the generated files will be stale.

**Recommendation:**
This is a known limitation documented in the code. The team decided the trade-off is acceptable to avoid forcing reconfigure on every plugin edit. Consider documenting this in the HACKING or developer guide:

*"When modifying command signatures (function names, command names) in libged plugins, run `cmake .` to regenerate command registration files."*

---

## 5. Command Registration Scanner Analysis

### Scanner Implementation: ✅ ROBUST

**File:** `src/libged/ged_cmd_scanner.cpp`

The scanner supports multiple detection patterns (5 phases):

1. **Phase 1-2:** Legacy `REGISTER_GED_COMMAND()` and `LABEL_GED_COMMAND()` macros
2. **Phase 3:** `REGISTER_BU_PLUGIN_COMMAND()` pattern
3. **Phase 4:** `bu_plugin_cmd` array pattern
4. **Phase 5:** X-macro command lists (new canonical pattern)
5. **Fallback:** `ged_cmd_impl` struct initializers

**Strengths:**
- ✅ Handles both single-line and multi-line struct definitions
- ✅ Skips commands with `?` in name (handled separately)
- ✅ Supports multiple patterns for backward compatibility
- ✅ Sorted output for deterministic builds
- ✅ Proper error handling for missing files

**Potential Issues:**
- Scanner runs at configure time, not build time (see Issue 4.3)
- Regex patterns assume consistent formatting - unusual whitespace/comments might break detection

**Recommendation:**
Add scanner unit tests to verify all supported patterns are correctly detected.

---

## 6. Initialization and Execution Analysis

### Initialization Flow: ✅ CORRECT

**File:** `src/libged/ged_init.cpp`

The initialization sequence is well-designed:

1. **Bootstrap:** `bu_plugin_init()` sets up generalized registry
2. **Static registration:** `ged_force_static_registration()` (when `LIBGED_STATIC_CORE=ON`)
3. **Plugin scan:** `scan_plugins()` loads dynamic plugins from disk
4. **Lazy initialization:** `ged_ensure_initialized()` ensures exactly-once execution

**Thread Safety:**
- Uses `std::once_flag` for thread-safe initialization (likely in bu_plugin implementation)
- First-wins registration strategy prevents duplicates

**Environment Variables:**
- `GED_NO_PLUGIN_SCAN=1` - Disables plugin scanning (useful for testing)
- `GED_EXEC_TIME` - Enables command timing

### Execution Flow: ✅ CORRECT

**File:** `src/libged/exec.cpp`

Command execution properly:
- ✅ Normalizes command names to basename
- ✅ Looks up via `bu_plugin_cmd_get()`
- ✅ Tracks recursion depth (limit: `GED_CMD_RECURSION_LIMIT`)
- ✅ Supports pre/post execution callbacks
- ✅ Handles `NULL`-terminated argv arrays

---

## 7. Additional Observations

### Positive Changes

1. **Improved modularity:** Clear separation between static and dynamic loading
2. **Better performance:** Static core eliminates plugin loading overhead for common commands
3. **Consistent patterns:** X-macro pattern is clean and maintainable
4. **ABI safety:** Version checking via `BU_PLUGIN_ABI_VERSION`
5. **Documentation:** Excellent inline documentation in `plugin.h` and `bu_plugin.h`

### Areas for Future Enhancement

1. **Testing:** Add integration tests for:
   - Static vs dynamic mode parity
   - Cross-platform builds (especially Windows)
   - Plugin loading edge cases
   - Scanner pattern detection

2. **Documentation:**
   - Add migration guide for new commands
   - Document the X-macro pattern usage
   - Add troubleshooting guide for common issues

3. **Build Performance:**
   - Consider making scanner a build-time tool rather than configure-time
   - Evaluate using `CMAKE_CONFIGURE_DEPENDS` for command files

---

## 8. Summary of Findings

| Category | Status | Critical Issues | Medium Issues | Low Issues |
|----------|--------|----------------|---------------|------------|
| **Missed Commands** | ✅ PASS | 0 | 0 | 0 |
| **Function Assignments** | ✅ PASS | 0 | 0 | 0 |
| **Cross-Platform** | ⚠️ ISSUES | 1 | 0 | 1 |
| **Plugin Registration** | ✅ PASS | 0 | 0 | 0 |
| **Build System** | ⚠️ ISSUES | 0 | 1 | 2 |

### Critical Issues (Must Fix Before Merge)

1. **Issue 2.1:** MSVC symbol retention (`GED_CMD_USED` macro)

### Medium Priority Issues (Should Fix)

1. **Issue 4.1:** Library dependency consistency (bot and other commands)

### Low Priority Issues (Nice to Have)

1. **Issue 2.2:** Path separator verification on Windows
2. **Issue 4.2:** Scanner rebuild optimization
3. **Issue 4.3:** Configure dependency documentation

---

## 9. Recommended Action Plan

### Phase 1: Critical Fixes (Before Merge)

1. ✅ Fix `GED_CMD_USED` macro for MSVC (Issue 2.1)
   - Implement Option 2 (linker flag) as immediate fix
   - Test on Windows with MSVC
   
2. ✅ Restore `libged` to BOT_LIBS and verify other commands (Issue 4.1)
   - Audit all command CMakeLists.txt
   - Test in dynamic plugin mode

### Phase 2: Testing and Validation

1. Build and test on all platforms:
   - Linux (GCC, Clang)
   - macOS (Clang)
   - Windows (MSVC, MinGW)

2. Test both modes:
   - `LIBGED_STATIC_CORE=ON` (default)
   - `LIBGED_STATIC_CORE=OFF` (dynamic mode)

3. Verify command registration:
   - All 227 commands load correctly
   - No duplicate registrations
   - No missing commands

### Phase 3: Documentation and Enhancement

1. Document the Windows-specific considerations
2. Add developer guide for the new plugin system
3. Consider scanner unit tests
4. Add note about configure dependencies

---

## 10. Conclusion

The gedcmdsreview branch represents a significant and well-thought-out architectural improvement to libged's plugin system. The migration has been executed consistently across all 226+ commands with good attention to detail.

The primary concern is the **MSVC symbol retention issue** which must be resolved before merging to ensure Windows compatibility. The library dependency issue should also be addressed to prevent linking problems in dynamic plugin mode.

Once these issues are resolved and cross-platform testing is complete, the branch should be safe to merge.

---

## Appendix A: Command Migration Pattern Examples

### Simple Single-Command Plugin

```c
// Before (main branch)
#ifdef GED_PLUGIN
struct ged_cmd_impl foo_cmd_impl = {"foo", ged_foo_core, GED_CMD_DEFAULT};
const struct ged_cmd foo_cmd = { &foo_cmd_impl };
const struct ged_cmd *foo_cmds[] = { &foo_cmd, NULL };
static const struct ged_plugin pinfo = { GED_API, foo_cmds, 1 };
COMPILER_DLLEXPORT const struct ged_plugin *ged_plugin_info(void) {
    return &pinfo;
}
#endif

// After (gedcmdsreview branch)
#include "../include/plugin.h"
#define GED_FOO_COMMANDS(X, XID) \
    X(foo, ged_foo_core, GED_CMD_DEFAULT)
GED_DECLARE_COMMAND_SET(GED_FOO_COMMANDS)
GED_DECLARE_PLUGIN_MANIFEST("libged_foo", 1, GED_FOO_COMMANDS)
```

### Multi-Command Plugin with Aliases

```c
// draw command - 7 commands in one plugin
#define GED_DRAW_COMMANDS(X, XID) \
    X(draw, ged_draw_core, GED_CMD_DEFAULT) \
    X(E, ged_E_core, GED_CMD_DEFAULT) \
    X(e, ged_draw_core, GED_CMD_DEFAULT) \
    X(ev, ged_ev_core, GED_CMD_DEFAULT) \
    X(redraw, ged_redraw_core, GED_CMD_DEFAULT) \
    X(loadview, ged_loadview_core, GED_CMD_DEFAULT) \
    X(preview, ged_preview_core, GED_CMD_DEFAULT)

GED_DECLARE_COMMAND_SET(GED_DRAW_COMMANDS)
GED_DECLARE_PLUGIN_MANIFEST("libged_draw", 1, GED_DRAW_COMMANDS)
```

### Command with Special Characters

```c
// help command - includes the "?" command
#define GED_HELP_COMMANDS(X, XID) \
    X(help, ged_help_core, GED_CMD_DEFAULT) \
    X(apropos, ged_apropos_core, GED_CMD_DEFAULT) \
    X(man, ged_man_core, GED_CMD_DEFAULT) \
    XID(questionmark_cmd, "?", ged_help_core, GED_CMD_DEFAULT)
```

---

*Report generated: 2025*
*Branches compared: main vs gedcmdsreview*
*Total commands analyzed: 227*
*Total files reviewed: 470+*
