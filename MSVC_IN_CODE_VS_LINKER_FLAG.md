# MSVC Symbol Retention: In-Code vs Linker Flag Trade-offs

## Executive Summary

**User Preference:** Avoid modifying the linker line; prefer a self-contained solution in the `GED_CMD_USED` macro.

**Recommendation:** Use `#pragma comment(linker, "/INCLUDE:...")` approach with macro generation, OR accept `/OPT:NOREF` as the most robust solution.

**Bottom Line:** Each approach has trade-offs. The in-code solutions are more fragile and complex than `/OPT:NOREF`, but can work if carefully implemented.

---

## Option Analysis

### Option 1: Current Approach - `/OPT:NOREF` Linker Flag ✅ SIMPLEST

**Implementation:**
```cmake
# In CMakeLists.txt
if(MSVC AND LIBGED_STATIC_CORE)
  target_link_options(libged PRIVATE "/OPT:NOREF")
endif()
```

**Pros:**
- ✅ **Simple and robust** - Single CMake line
- ✅ **Guaranteed to work** - Disables entire optimization pass
- ✅ **Low maintenance** - No per-symbol intervention needed
- ✅ **No fragility** - Cannot be broken by code changes
- ✅ **Clean separation** - Build config separate from source code
- ✅ **Standard practice** - Similar to how Linux kernel handles similar issues

**Cons:**
- ❌ Requires CMake build system modification
- ❌ User prefers avoiding linker line changes
- ❌ Disables optimization broadly (though only for libged, not whole program)
- ❌ Slightly increases binary size (keeps all symbols, not just needed ones)

**Fragility:** **Low** - Very robust, hard to break

---

### Option 2: `__declspec(selectany)` ❌ DOESN'T SOLVE THE PROBLEM

**Implementation:**
```c
#if defined(_MSC_VER)
#define GED_CMD_USED __declspec(selectany)
#else
#define GED_CMD_USED __attribute__((used))
#endif
```

**What it does:**
- Marks symbol as COMDAT (Common Data)
- Allows multiple definitions, linker picks one
- Used for template instantiations and inline variables

**Pros:**
- ✅ Self-contained in macro
- ✅ No CMake changes needed

**Cons:**
- ❌ **Does NOT prevent symbol elimination** - `selectany` allows duplicates but doesn't force retention
- ❌ Not designed for this use case
- ❌ Linker can still remove the symbol if unreferenced
- ❌ **Will NOT fix the problem**

**Verdict:** **REJECTED** - Doesn't solve the core issue

---

### Option 3: `#pragma comment(linker, "/INCLUDE:symbol")` ⚠️ COMPLEX BUT WORKABLE

**Implementation Approach A - Manual per-symbol pragmas:**
```c
#ifdef _MSC_VER
#define GED_CMD_USED
#define REGISTER_GED_COMMAND(cmd_symbol)                                          \
    const struct ged_cmd cmd_symbol##_rcmd = { &cmd_symbol##_impl };             \
    extern "C" const struct ged_cmd * const __ged_cmd_ptr_##cmd_symbol =         \
        &cmd_symbol##_rcmd;                                                       \
    __pragma(comment(linker, "/INCLUDE:__ged_cmd_ptr_" #cmd_symbol))
#else
#define GED_CMD_USED __attribute__((used))
// ... normal definition
#endif
```

**How it works:**
- `#pragma comment(linker, "/INCLUDE:symbol")` forces linker to include specific symbol
- Embedded directly in the REGISTER_GED_COMMAND macro
- Self-contained, no CMake changes

**Pros:**
- ✅ Self-contained in source code
- ✅ No CMake modifications needed
- ✅ Precise - only includes exactly the symbols needed
- ✅ Standard MSVC feature for this exact purpose

**Cons:**
- ❌ **Symbol name must be known at preprocessing time** - requires macro stringification
- ❌ **Name mangling issues** - Must match exact linker symbol name (with underscores, etc.)
- ❌ **Fragile** - If REGISTER_GED_COMMAND macro changes, pragma must change too
- ❌ **Harder to debug** - Pragma errors are cryptic
- ❌ **C++ name mangling** - `extern "C"` helps but still fragile
- ❌ **Platform-specific complexity** in macro definitions

**Fragility:** **Medium-High** - Symbol naming must be exact

---

### Option 4: Scanner-Generated Pragma File ⚠️ MODERATELY COMPLEX

**Implementation:**
Modify `ged_cmd_scanner.cpp` to generate an MSVC-specific header with pragmas:

```cpp
// Generated: ged_static_includes_msvc.h (MSVC only)
#ifdef _MSC_VER
#pragma comment(linker, "/INCLUDE:__ged_cmd_ptr_threeptarb_cmd")
#pragma comment(linker, "/INCLUDE:__ged_cmd_ptr_adc_cmd")
// ... for all 227 commands
#endif
```

Then include in plugin.h or ged_init.cpp.

**Pros:**
- ✅ Automated generation - scanner knows all command names
- ✅ Exact symbol names - no manual typing
- ✅ Self-contained in code (header file)
- ✅ No per-command macro complexity

**Cons:**
- ❌ Requires scanner modifications
- ❌ Generates another intermediate file
- ❌ Still fragile to symbol naming changes
- ❌ Build system complexity (where to include the generated header?)
- ❌ Not truly "self-contained" - depends on generated file

**Fragility:** **Medium** - Automated but still dependent on exact symbol names

---

### Option 5: `__declspec(dllexport)` ❌ POLLUTES EXPORT TABLE

**Implementation:**
```c
#ifdef _MSC_VER
#define GED_CMD_USED __declspec(dllexport)
#else
#define GED_CMD_USED __attribute__((used))
#endif
```

**Pros:**
- ✅ Forces symbol retention
- ✅ Self-contained in macro

**Cons:**
- ❌ **Pollutes DLL export table** with 227+ internal symbols
- ❌ Exposes internal implementation details
- ❌ Not appropriate for static library (libged is typically static)
- ❌ Performance impact (export table lookups)
- ❌ **Wrong tool for the job**

**Verdict:** **REJECTED** - Inappropriate for static builds

---

## Detailed Comparison

| Aspect | `/OPT:NOREF` | `#pragma comment` (per-symbol) | `#pragma comment` (generated) |
|--------|--------------|-------------------------------|------------------------------|
| **Simplicity** | ⭐⭐⭐⭐⭐ Trivial | ⭐⭐ Complex macros | ⭐⭐⭐ Medium |
| **Robustness** | ⭐⭐⭐⭐⭐ Very robust | ⭐⭐ Fragile | ⭐⭐⭐ Moderately robust |
| **Maintenance** | ⭐⭐⭐⭐⭐ None needed | ⭐⭐ High (name changes break it) | ⭐⭐⭐ Medium |
| **Self-contained** | ⭐ CMake changes | ⭐⭐⭐⭐⭐ Fully in code | ⭐⭐⭐⭐ Header file |
| **Debugging** | ⭐⭐⭐⭐⭐ Clear | ⭐⭐ Cryptic pragma errors | ⭐⭐⭐ Generated file issues |
| **Binary size** | ⭐⭐⭐ Keeps all symbols | ⭐⭐⭐⭐⭐ Only needed symbols | ⭐⭐⭐⭐⭐ Only needed symbols |
| **Correctness** | ⭐⭐⭐⭐⭐ Guaranteed | ⭐⭐⭐ Must match names exactly | ⭐⭐⭐⭐ Automated name matching |

---

## Fragility Analysis

### `/OPT:NOREF` Fragility: **VERY LOW** ✅

**Failure modes:**
- CMake syntax error (caught immediately at configure time)
- Wrong target name (caught at configure time)

**Resilience:**
- Works regardless of symbol names
- Works regardless of code structure
- Works with any number of commands
- Immune to refactoring

### `#pragma comment` (per-symbol) Fragility: **HIGH** ⚠️

**Failure modes:**
- Symbol name typo → Link error
- Name mangling mismatch → Link error
- C vs C++ linkage confusion → Link error
- Macro expansion issues → Compile error
- Underscore prefix changes → Link error

**Example fragility:**
```c
// Macro creates: __ged_cmd_ptr_threeptarb_cmd
#pragma comment(linker, "/INCLUDE:__ged_cmd_ptr_threeptarb_cmd")  // Must match exactly!

// If symbol is actually: ___ged_cmd_ptr_threeptarb_cmd (3 underscores)
// → Link error!
```

**Resilience:**
- Breaks if symbol naming changes
- Breaks if linker behavior changes
- Breaks if macro refactored
- Must be tested on actual MSVC builds

### Scanner-Generated Pragmas Fragility: **MEDIUM** ⚠️

**Failure modes:**
- Scanner symbol detection bug → Missing pragmas
- Symbol name extraction error → Wrong names in pragmas
- Generated file not included properly → No effect
- Build system doesn't regenerate → Stale pragmas

**Resilience:**
- Automated, less manual error
- Still dependent on exact symbol names
- Requires scanner maintenance
- Build dependency complexity

---

## Real-World Experience

### Linux Kernel Approach
The Linux kernel uses linker scripts and flags (similar to `/OPT:NOREF`) for section-based registration, not per-symbol pragmas. This is considered the **robust, maintainable approach**.

### GCC Plugin System
GCC's plugin system uses `__attribute__((used))` on supported platforms and **accepts** that some platforms need linker flags for full functionality.

### Industry Practice
Most large cross-platform projects (LLVM, Chromium, etc.) use:
- Attributes where available (GCC/Clang)
- Linker flags where necessary (MSVC)
- **NOT** per-symbol pragmas (too fragile)

---

## Recommendations

### Recommended: Keep `/OPT:NOREF` (Option 1) ✅

**Rationale:**
- Most robust and maintainable
- Industry standard approach
- Already implemented and tested
- Single point of configuration

**Trade-off:**
- Requires one CMake conditional
- User prefers avoiding this, but it's the best engineering choice

### Alternative: Scanner-Generated Pragmas (Option 4) ⚠️

**If user strongly objects to CMake changes:**

**Implementation plan:**
1. Modify `ged_cmd_scanner.cpp` to generate `ged_msvc_includes.h`
2. Include this header in `ged_init.cpp` or one central location
3. Test thoroughly on MSVC

**Pros:**
- No CMake linker changes
- Automated generation reduces manual error

**Cons:**
- More complex build process
- Another generated file to manage
- Still fragile to symbol naming
- Requires scanner modification and testing

### Not Recommended: Per-Symbol Pragmas (Option 3) ❌

**Rationale:**
- Too fragile
- High maintenance burden
- Error-prone
- Not used in major projects for good reason

---

## Specific Symbol Name Concerns

### Testing Required
The `#pragma comment(linker, "/INCLUDE:symbol")` approach requires **exact** symbol names:

```c
// What the preprocessor sees:
#pragma comment(linker, "/INCLUDE:__ged_cmd_ptr_threeptarb_cmd")

// What the linker sees (might be):
__ged_cmd_ptr_threeptarb_cmd     // C symbol
_ged_cmd_ptr_threeptarb_cmd      // Different calling convention
___ged_cmd_ptr_threeptarb_cmd    // With name prefix
```

The **exact format depends on:**
- Calling convention (`__cdecl`, `__stdcall`, etc.)
- Name mangling rules (C vs C++)
- Target architecture (x86 vs x64)
- MSVC version

**This is why `/OPT:NOREF` is superior** - it doesn't care about symbol names.

---

## Conclusion

### Engineering Recommendation: `/OPT:NOREF` ✅

**Best choice because:**
1. **Robust** - Can't break with code changes
2. **Simple** - One line in CMake
3. **Standard** - Used by major projects
4. **Maintainable** - No ongoing upkeep
5. **Debuggable** - Clear failure modes

### Compromise: Scanner-Generated Pragmas ⚠️

**Acceptable if CMake changes are truly unacceptable:**
- More complex but automated
- Requires careful implementation
- Needs thorough MSVC testing
- Still more fragile than `/OPT:NOREF`

### Not Recommended: Manual Per-Symbol Pragmas ❌

**Avoid because:**
- Too fragile
- Too error-prone
- Not maintainable
- No major project uses this approach

---

## User's Question Answered

> What are the drawbacks of the in-code solution as opposed to the linker line solution?

**Drawbacks of in-code `#pragma comment` approach:**

1. **Fragility**: Symbol names must match exactly; breaks if naming changes
2. **Platform-specific complexity**: MSVC-specific macro logic clutters source
3. **Debugging difficulty**: Pragma link errors are cryptic and hard to diagnose
4. **Name mangling risk**: C vs C++ linkage, calling conventions, architecture differences
5. **Maintenance burden**: Each command addition requires careful pragma management
6. **Testing requirements**: Must be tested on actual MSVC, can't verify on Linux/Mac
7. **Not industry standard**: Major projects don't use this pattern for good reason

> Is the former more fragile?

**Yes, significantly more fragile.**

- `/OPT:NOREF`: Changes to code structure, symbol names, or macro definitions don't break it
- `#pragma comment`: Any mismatch between pragma string and actual symbol name causes link failure

**Bottom line:** The linker flag is the engineering best practice. The in-code solution can work but is objectively more fragile and complex.
