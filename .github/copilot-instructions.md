# BRL-CAD Coding Agent Instructions

## Repository Overview

BRL-CAD is a powerful cross-platform open source combinatorial solid modeling system with over 1 million lines of code comprising 400+ tools and utilities. It includes an interactive 3D solid geometry editor, high-performance ray-tracer, geometric analysis tools, and extensive libraries for solid modeling operations. The project has been actively developed since 1979 and supports a wide variety of geometric representations including CSG primitives, NURBS, and mesh geometry.

**Key Repository Facts:**
- **Size:** ~1M+ lines of code, 400+ tools
- **Languages:** Primarily C/C++ with Tcl/Tk, some Python
- **Build System:** CMake 3.20+ required
- **Target Platforms:** Cross-platform (Linux, Windows, macOS, BSD, etc.)
- **External Dependencies:** Managed via separate 'bext' repository
- **License:** LGPL (libraries) and BSD (build system)

## Build System and Configuration

### Quick Build Commands

**CRITICAL:** Always use the cached bext approach for significantly faster builds (can reduce configure time from hours to minutes).

#### Standard Configuration (Recommended):
```bash
# If using cached bext (strongly recommended for agents):
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DBRLCAD_ENABLE_STRICT=OFF -DBRLCAD_EXT_DIR=/path/to/bext_output

# For fresh build without cached bext (very slow):
mkdir build && cd build  
cmake .. -DCMAKE_BUILD_TYPE=Release -DBRLCAD_ENABLE_STRICT=OFF -DBRLCAD_BUNDLED_LIBS=ON
```

#### Minimal/Fast Component Build:
For development targeting specific functionality, use BRLCAD_COMPONENTS to build only required parts:
```bash
# Example: Only build libbu, libbn, and related dependencies
cmake .. -DCMAKE_BUILD_TYPE=Release -DBRLCAD_ENABLE_STRICT=OFF -DBRLCAD_COMPONENTS="libbu:libbn"

# Example: Build geometry tools only
cmake .. -DCMAKE_BUILD_TYPE=Release -DBRLCAD_ENABLE_STRICT=OFF -DBRLCAD_COMPONENTS="gtools"
```

### Build Process
```bash
# Build (parallel, then serial fallback if needed)
make -j$(nproc) || make

# Run tests 
make test
# OR
ctest

# Install
make install
```

### Critical Configuration Options

- **BRLCAD_ENABLE_STRICT=OFF**: Always use OFF to avoid strict compiler warnings that can break builds
- **CMAKE_BUILD_TYPE=Release**: Use Release for performance; Debug only when debugging
- **BRLCAD_EXT_DIR**: Specify precompiled external dependencies directory (CRITICAL for speed)
- **BRLCAD_BUNDLED_LIBS=ON**: Use if not using BRLCAD_EXT_DIR (much slower)
- **BRLCAD_COMPONENTS**: Limit build scope to specific libraries/tools (dramatically faster)

### External Dependencies (bext) Strategy

The CI system uses sophisticated caching for external dependencies. If agents have access to cached bext builds:

1. **Use cached bext**: Set `BRLCAD_EXT_DIR` to cached bext output directory
2. **Cache key logic**: Based on bext commit SHA + OS + compiler version
3. **Build time impact**: Can reduce configure from hours to minutes

### Build Issues and Workarounds

**Common Problems:**
- **X11/OpenGL not found**: Normal on headless systems; build continues with reduced functionality  
- **Configure timeout**: External dependency download/build can take 1-6 hours on first run
- **Parallel build failures**: Fall back to serial builds (`make` without `-j`)
- **Memory issues**: Large builds may need `make -j1` or `-j2` instead of full parallelism
- **CMake configuration interruption**: If configure stops unexpectedly, clear build directory and retry
- **Target not found errors**: Ensure full CMake configure completed before attempting builds

**Build Timing:**
- Full configure (first time): 2-6 hours depending on dependencies
- Full build: 30-90 minutes  
- With cached bext: Configure 5-15 minutes, build 20-60 minutes
- Component builds: Can be 5-20 minutes total
- **Critical**: Always allow extra time for bext dependency builds on first configure

**Environment Requirements:**
- **Disk space**: Minimum 3GB, recommend 10GB+ for full builds
- **Memory**: 4GB+ recommended, 8GB+ for parallel builds
- **Time**: Plan for several hours on first build without cached dependencies

## Project Architecture and Layout

### Source Organization
```
src/
├── libbu/          # Basic utilities (foundation library)
├── libbn/          # Basic numerics (depends on libbu)
├── libbg/          # Basic geometry (depends on libbn, libbu)
├── libbv/          # Basic viewing (depends on libbg)
├── librt/          # Ray-tracing core (depends on libbrep, libnmg, libbv)
├── libged/         # Geometry editor (high-level interface)
├── mged/           # Interactive geometry editor
├── archer/         # Tk-based GUI
├── qged/           # Qt-based GUI
├── conv/           # Geometry converters
├── gtools/         # Geometry tools
├── util/           # Utility programs
└── ...
```

### Key Configuration Files
- **CMakeLists.txt**: Main build configuration
- **src/source_dirs.cmake**: Component dependency definitions
- **misc/CMake/**: Extensive build system configurations
- **misc/astyle.opt**: Code formatting rules (K&R style)
- **.github/workflows/push.yml**: CI configuration (reference for caching)

### Library Dependencies (in build order)
```
libbu (foundation)
├── libbn (numerics)
    ├── libbg (geometry)
        ├── libbv (viewing)
        ├── libnmg (n-manifold geometry)
        ├── libbrep (boundary representation)
            ├── librt (ray-tracing)
                ├── libwdb (write database)
                ├── libanalyze (analysis)
                ├── liboptical (optics)
                ├── libdm (display manager)
                    ├── libged (geometry editor)
```

## Testing and Validation

### Test Execution
```bash
# Run all tests
make test
# OR
ctest

# Run specific test categories
ctest -L Benchmark    # Benchmark tests
ctest -L Unit         # Unit tests

# Run benchmark suite
make benchmark
```

### Test Structure
- **Unit tests**: Located in `src/*/tests/` directories
- **Integration tests**: In `regress/` directory  
- **Benchmark suite**: In `bench/` directory
- **CTest framework**: Integrated with CMake build

### Code Quality Tools
```bash
# Format code (follow K&R style)
astyle --options=misc/astyle.opt src/path/to/file.c

# Lint geometry files (built-in tool)
bin/lint -g filename.g

# Repository integrity check
make distcheck-repo-verify
```

## Continuous Integration

The CI system (`.github/workflows/push.yml`) provides the definitive reference for working build configurations:

### CI Build Process
1. **Cache bext**: Uses SHA + OS + compiler for cache key
2. **Configure**: Uses Ninja generator with Release build type
3. **Build**: Attempts parallel (`-j2`), falls back to serial
4. **Test**: Runs `check` target for validation
5. **Package**: Creates distribution packages

### CI Environment
- **Generator**: Ninja (faster than Make)
- **Parallelism**: `-j2` with serial fallback
- **Dependencies**: Cached bext with `BRLCAD_EXT_DIR`
- **Features**: Enables Qt, OpenGL when available

## Development Workflow

### Making Changes
1. **Scope changes**: Use BRLCAD_COMPONENTS for targeted builds
2. **Build incrementally**: `make target_name` for specific targets
3. **Test changes**: Use `ctest -R pattern` for relevant tests
4. **Format code**: Follow K&R style with astyle
5. **Validate**: Run `make check` before submitting

### Key Validation Steps
```bash
# After making changes:
make clean                    # Clean previous build
make target_name             # Build specific target  
ctest -R relevant_pattern    # Run relevant tests
make check                   # Full validation

# Quick validation without full build:
cmake --build . --target check  # Alternative to make check
ctest --output-on-failure       # See test failures immediately

# Benchmark validation (if needed):
make benchmark                   # Run benchmark suite (requires completed build)
```

### Performance Tips
- **Use Release builds** for performance testing
- **Leverage BRLCAD_COMPONENTS** for faster iteration
- **Cache bext outputs** when possible
- **Use parallel builds** but have serial fallback ready

## Important Notes for Agents

1. **Trust these instructions**: Only search for additional information if these instructions are incomplete or incorrect
2. **Always use BRLCAD_ENABLE_STRICT=OFF** unless specifically testing strict compilation  
3. **Leverage BRLCAD_COMPONENTS** for targeted development to avoid unnecessary compilation time
4. **Expect long configure times** on first build without cached bext (2-6 hours) - this is normal
5. **Use cached bext when available** - can reduce build times by 10x or more
6. **Fall back to serial builds** if parallel builds fail (`make` instead of `make -j`)
7. **Build system is mature but complex** - follow established patterns rather than trying to simplify
8. **Cross-platform considerations**: Code must work on Windows, Linux, macOS, BSD
9. **Plan for significant time investment** on initial setup - this is a large, complex system
10. **Monitor disk space and memory** - builds can be resource intensive

This system has extensive functionality and build complexity, but following these instructions should enable efficient development without the trial-and-error typically required for such a large codebase.