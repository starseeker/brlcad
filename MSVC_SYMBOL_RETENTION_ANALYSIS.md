# MSVC Symbol Retention - Why Both Mechanisms Are Needed

## Question

The code review identified an MSVC symbol retention issue and recommended adding `/OPT:NOREF` linker flag. However, there's already a scanner-generated `ged_static_registration.cpp` file that explicitly references all command symbols. Does this existing mechanism already solve the problem, making `/OPT:NOREF` unnecessary?

## Short Answer

**No - both mechanisms are needed.** The `/OPT:NOREF` linker flag is still necessary even with the generated registration code.

## Detailed Explanation

### How the System Works

1. **Symbol Definition (in command source files):**
   ```c++
   // In each command file (e.g., 3ptarb.c)
   REGISTER_GED_COMMAND(threeptarb_cmd)
   // Expands to:
   GED_CMD_USED const struct ged_cmd threeptarb_cmd_rcmd = { &threeptarb_cmd_impl };
   extern "C" GED_CMD_USED const struct ged_cmd * const __ged_cmd_ptr_threeptarb_cmd = &threeptarb_cmd_rcmd;
   ```

2. **Scanner-Generated Code (ged_static_registration.cpp):**
   ```c++
   extern "C" const struct ged_cmd * const __ged_cmd_ptr_threeptarb_cmd;
   // ... for all 227 commands ...
   
   extern "C" void ged_force_static_registration(void)
   {
       (void)ged_register_command(__ged_cmd_ptr_threeptarb_cmd);
       // ... for all 227 commands ...
   }
   ```

3. **Initialization (ged_init.cpp):**
   ```c++
   #if defined(LIBGED_STATIC_CORE)
       ged_force_static_registration();  // Called at startup
   #endif
   ```

### The Problem on MSVC

The issue is with MSVC's **link-time code generation (LTCG)** and optimization:

1. **GCC/Clang Solution:**
   - `__attribute__((used))` marks symbols to prevent elimination
   - Compiler guarantees symbol retention even if no direct references found

2. **MSVC Challenge:**
   - No equivalent to `__attribute__((used))` (GED_CMD_USED is empty on MSVC)
   - Even with `ged_force_static_registration()` referencing symbols via `extern`
   - The linker's `/OPT:REF` optimization runs in multiple passes:
     - **Pass 1:** Identifies "unreferenced" symbols in object files
     - **Pass 2:** Strips those symbols
     - **Pass 3:** Links remaining symbols
   
3. **The Chicken-and-Egg Problem:**
   - Command source files create `__ged_cmd_ptr_xxx` symbols
   - These appear "unreferenced" in their own translation units
   - MSVC's linker may strip them **before** seeing the `extern` references in `ged_static_registration.cpp`
   - Result: Undefined symbol errors or missing commands

### Why Both Mechanisms Are Necessary

**Generated Registration Code Alone:**
- ✅ Creates explicit references to symbols
- ✅ Ensures symbols are used if they exist
- ❌ **Cannot prevent MSVC from stripping symbols before linking**
- The `extern` declarations don't force retention of the actual definitions

**/OPT:NOREF Flag:**
- ✅ Tells MSVC linker to skip the "remove unreferenced symbols" optimization
- ✅ Ensures all symbol definitions survive to final linking
- ✅ Works together with generated code to guarantee symbols are both present and used

### Analogy

Think of it like a relay race:

- **Generated code** = "The runner is ready to receive the baton"
- **/OPT:NOREF** = "Make sure the previous runner actually shows up to the handoff"

Without `/OPT:NOREF`, MSVC might eliminate the "previous runner" (symbol definitions) before the handoff (linking) even happens.

## Testing the Theory

To verify this is a real issue, you could:

1. **Test WITHOUT /OPT:NOREF on MSVC:**
   - Build with `LIBGED_STATIC_CORE=ON`
   - Remove the `/OPT:NOREF` flag
   - Expect: Link errors for undefined `__ged_cmd_ptr_xxx` symbols, OR
   - Expect: Successful link but commands missing at runtime

2. **Test WITH /OPT:NOREF on MSVC:**
   - Build with `LIBGED_STATIC_CORE=ON`
   - Keep the `/OPT:NOREF` flag
   - Expect: Clean build, all 227 commands available

## Alternative Solutions (Not Recommended)

Other theoretical approaches that would be more complex:

1. **Explicit `/INCLUDE:` directives:**
   ```cmake
   # Would need to list all 227 symbols
   target_link_options(libged PRIVATE 
     "/INCLUDE:__ged_cmd_ptr_threeptarb_cmd"
     "/INCLUDE:__ged_cmd_ptr_adc_cmd"
     # ... 225 more ...
   )
   ```
   - Too verbose and fragile

2. **Force all command symbols to be exported:**
   ```c
   __declspec(dllexport) const struct ged_cmd * const __ged_cmd_ptr_xxx
   ```
   - Pollutes export table
   - Not needed for static library

3. **Put all commands in a single .cpp file:**
   - Defeats modularity
   - Makes incremental compilation much slower

## Conclusion

The `/OPT:NOREF` linker flag is the **simplest and most reliable** solution for MSVC. It works in conjunction with the scanner-generated registration code:

- **Scanner code:** Ensures symbols are actually used during initialization
- **/OPT:NOREF:** Ensures symbols survive to be available for use

Both mechanisms are complementary and necessary for a robust cross-platform solution.

## Implementation Status

✅ **ALREADY IMPLEMENTED** - The fix is correct and should be kept:

```cmake
# src/libged/CMakeLists.txt
if (LIBGED_STATIC_CORE)
  target_link_libraries(libged PRIVATE libged_plugin_deps)
  # On MSVC, prevent the linker from stripping "unused" static command registration symbols
  if(MSVC)
    target_link_options(libged PRIVATE "/OPT:NOREF")
  endif()
endif()
```

## References

- MSVC `/OPT:REF` documentation: https://docs.microsoft.com/en-us/cpp/build/reference/opt-optimizations
- Link-time code generation: https://docs.microsoft.com/en-us/cpp/build/reference/ltcg-link-time-code-generation
