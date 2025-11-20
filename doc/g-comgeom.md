# g-comgeom - BRL-CAD to COMGEOM Exporter

## NAME
g-comgeom - export BRL-CAD database geometry to COMGEOM format

## SYNOPSIS
```
g-comgeom [options] input.g [objects...]
```

## DESCRIPTION
**g-comgeom** exports BRL-CAD database geometry to COMGEOM (COMbinatorial GEOMetry) format, producing card deck files suitable for use with ballistic analysis tools. It is the modern successor to the **vdeck** command with enhanced primitive support, multiple format versions, and JSON reporting capabilities.

The exporter generates three output files:
- **prefix.sol** - Solid table (primitive definitions)
- **prefix.reg** - Region table (boolean tree expressions)
- **prefix.id** - Ident table (region attributes)

## OPTIONS

### Version and Format
- **-V** *version*  
  Specify COMGEOM version (0, 1, 4, or 5). Default: 5  
  - v1: Fixed-column legacy format
  - v4: Token-based format with improved precision
  - v5: Enhanced format with extended attributes

### Output Control
- **-o** *prefix*  
  Output file prefix. Default: "output"  
  Generates: prefix.sol, prefix.reg, prefix.id

- **-s** *num*  
  Starting solid number. Default: 1

- **-r** *num*  
  Starting region number. Default: 1

### Reporting
- **-j** *path*  
  Generate JSON report at specified path  
  Report includes: solid/region counts, primitive breakdown, errors, warnings

- **-v**  
  Verbose output (can be repeated for increased verbosity)

- **-h**  
  Display help message

## SUPPORTED PRIMITIVES

### Fully Supported
- **TOR** - Torus
- **TGC** - Truncated General Cone (including RCC, TRC, REC, TEC variants)
- **ELL** - Ellipsoid (including SPH, ELL1, ELLG variants)
- **ARB8** - Arbitrary polyhedron (8 vertices, reduced to ARB4-8 as appropriate)
- **ARBN** - Arbitrary convex polyhedron (N planes)
- **ARS** - Arbitrary faceted solid
- **HALF** - Half-space

### Unsupported (Warning Generated)
- **BOT** - Bag of Triangles (future: tessellation fallback)
- **BREP** - Boundary Representation (future: triangulation fallback)
- Other primitives not in COMGEOM specification

## PRIMITIVE CLASSIFICATION

The exporter automatically classifies primitives into their most specific COMGEOM types based on geometric properties (Requicha & Voelcker 1982):

### TGC Classification
- **RCC** - Right Circular Cylinder (a = b, c = d, |a| = |c|)
- **TRC** - Truncated Right Cone (a = b, c = d, |a| ≠ |c|)
- **REC** - Right Elliptical Cylinder (c = d, |a| ≠ |b|)
- **TEC** - Truncated Elliptical Cone (c = d = 0, |a| ≠ |b|)
- **TGC** - General case

### ELL Classification
- **SPH** - Sphere (|a| = |b| = |c|)
- **ELL1** - Ellipsoid of Revolution (two axes equal)
- **ELLG** - General ellipsoid

## COMGEOM FORMAT VERSIONS

### Version 1 (v1) - Legacy Format
Fixed-column format with 10-character fields:
```
SSSS TYPE   PPPPPPPPPPPPPPPPPPPPPPPPPPPPPP...
```
Where:
- SSSS = Solid number (5 columns)
- TYPE = Primitive type (5 columns)
- P = Parameter fields (10 columns each)

### Version 4 (v4) - Token Format
Token-based format with improved precision:
```
solid_num type param1 param2 ... paramN
```

### Version 5 (v5) - Extended Format
Enhanced format with region attributes:
```
Solid table: same as v4
Region table: parenthesized boolean expressions
Ident table: reg_num ident air material LOS color_rgb
```

## BOOLEAN TREE SERIALIZATION

Region boolean trees are serialized using standard operators (Foley & Van Dam 1996):
- **u** - Union
- **+** - Intersection
- **-** - Subtraction

Format example (v5):
```
region_num u solid1 - solid2 + solid3
```

## REGION ATTRIBUTES (v5)

The ident file contains region attributes:
- **Region number** - Unique identifier
- **Ident** - Region identification code
- **Air code** - Air material code (0 = not air)
- **Material code** - Material identifier
- **LOS** - Line-of-sight thickness (percentage)
- **Color** - RGB color values (optional)

## EXIT STATUS
- **0** - Success (all primitives exported)
- **1** - Error (missing input, unsupported primitives, write failures)

## EXAMPLES

### Basic Export
```bash
g-comgeom input.g
```
Exports all primitives from input.g to output.sol, output.reg, output.id (v5 format)

### Specific Objects with Custom Prefix
```bash
g-comgeom -o mymodel -V 4 input.g tank.r turret.r
```
Exports specified regions using v4 format to mymodel.sol, mymodel.reg, mymodel.id

### With JSON Report
```bash
g-comgeom -j report.json -v input.g
```
Exports with verbose output and generates report.json containing statistics

### Custom Starting Numbers
```bash
g-comgeom -s 1000 -r 2000 -o export input.g
```
Starts solid numbering at 1000 and region numbering at 2000

## JSON REPORT FORMAT

When **-j** is specified, a JSON report is generated:

```json
{
  "version": 5,
  "output_prefix": "output",
  "statistics": {
    "solids": 42,
    "regions": 8,
    "unsupported": 2,
    "errors": 0,
    "warnings": 1
  },
  "primitives": {
    "TOR": 3,
    "TGC": 12,
    "ELL": 8,
    "ARB": 15,
    "ARBN": 2,
    "ARS": 0,
    "HALF": 2,
    "BOT": 2,
    "BREP": 0,
    "Other": 0
  }
}
```

## ENVIRONMENT VARIABLES

Currently, no environment variables are defined. Future versions may support:
- **BRLCAD_COMGEOM_PRECISION** - Output precision control
- **BRLCAD_COMGEOM_UNITS** - Unit conversion options

## ACADEMIC REFERENCES

The algorithmic approaches used in g-comgeom are based on:

1. **Requicha, A. A. G., & Voelcker, H. B. (1982)**  
   "Solid Modeling: A Historical Summary and Contemporary Assessment"  
   IEEE Computer Graphics and Applications, 2(2), 9-22.  
   *Basis for CSG primitive classification*

2. **Foley, J. D., & Van Dam, A. (1996)**  
   "Computer Graphics: Principles and Practice" (2nd ed.)  
   Addison-Wesley.  
   *Scene graph traversal and boolean tree serialization*

3. **Hoffmann, C. M. (1989)**  
   "Geometric and Solid Modeling: An Introduction"  
   Morgan Kaufmann.  
   *Format conversion strategies and data exchange*

## SEE ALSO
- **vdeck**(1) - Legacy COMGEOM exporter
- **comgeom-g**(1) - COMGEOM to BRL-CAD importer
- **g-obj**(1), **g-stl**(1) - Other geometry exporters

## DIAGNOSTICS

Warning messages are generated for:
- Unsupported primitive types (BOT, BREP, etc.)
- Missing objects in database
- Write failures

Error messages indicate:
- Cannot open database
- Cannot open output files
- Invalid command-line options

## BUGS

- Boolean tree serialization is currently simplified (stub implementation)
- ARB reduction to minimal forms (ARB4-ARB7) not yet implemented
- No pattern matching for object filters

Report bugs to: BRL-CAD bug tracker

## HISTORY
g-comgeom was developed in 2025 as a modern replacement for vdeck, incorporating advances in geometric data exchange and providing enhanced primitive support.

## AUTHORS
BRL-CAD Development Team
