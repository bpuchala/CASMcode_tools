---
name: casm
description: Using CASM for crystallography, symmetry analysis, group theory, structure enumeration, high-throughput DFT calculation, machine-learn potential training data construction, cluster expansion, Monte Carlo, kinetic Monte Carlo, free energy integration, phase diagram construction
---

## Key Concepts

- **Lattice**: A 3-dimensional lattice
- **Structure**: A crystal structure
- **Prim**: Primitive crystal structure — lattice vectors, basis sites, and allowed degrees of freedom (DoF) (occupants, displacements, magnetic spin, strain).
- **Supercell**: Superlattice of the prim defined by a 3×3 integer transformation matrix. Named e.g. `SCEL4_2_2_1_0_0_0`.
- **Configuration**: A specific DoF assignment within a supercell. Occupation stored as integer array indexing into each site's `occ_dof` list.
- **DoF**: Degrees of Freedom (DoF) include site occupants (discrete, may include magnetic spin), site displacements (continuous), site magnetic spin (continuous), and homogeneous lattice strain (continuous).
- **Clexulator**: A CLuster EXpansion calcULATOR (clexulator), used to evaluate cluster expansion basis set functions. Based on the symmetry of a prim, CASM generates code for efficient calculation of basis functions for a configuration representing one crystal structure consistent with the prim. This source code may be compiled, linked, and used at runtime via Clexulator.
- **Cluster expansion**: basis-function expansion of energy (or other crystal properties) for a configuration; coefficients (ECIs) fit to DFT data.
- **Canonical form**: symmetry-equivalent representative of a configuration, supercell, lattice, etc., used for comparison

---

## Packages

**Use MCP server `superchivo` to search Python package docs.**

Two Python namespaces: **`libcasm`** (C++ extension packages) and **`casm`** (pure Python packages).

### libcasm-global
- **Import**: `libcasm.casmglobal`, `libcasm.counter`
- Foundation: global constants, generic utilities (Counter), no CASM deps.

### libcasm-xtal
- **Import**: 
  - `import libcasm.xtal as xtal` - `Lattice`, `Structure`, `Prim`, `SymOp`, etc.
  - `import libcasm.xtal.prims as xtal_prims`  — `BCC`, `FCC`, `HCP`, `cubic` 
  - `import libcasm.xtal.lattices as xtal_lattices`  — `BCC`, `FCC`, `HCP`, `cubic`, `from_lattice_parameters`, `hexagonal`, etc. 
  - `import libcasm.xtal.structures as xtal_structures`  — `BCC`, `FCC`, `HCP` 
- Crystallography building blocks, superlattice enumeration, primitive/reduced cell finding, site degrees of freedom (occupant, displacement, strain).

Key constructors:
```python
import libcasm.xtal as xtal
import numpy as np

# Lattice vectors as columns of a 3x3 array
lattice = xtal.Lattice(np.array([[a,0,0],[0,b,0],[0,0,c]]).T)

# Occupant with per-atom properties (e.g. Ising spin)
occ = xtal.Occupant(
    name="A.up",
    atoms=[xtal.AtomComponent(name="A", coordinate=np.zeros(3),
                               properties={"Cmagspin": np.array([1.0])})],
)

prim = xtal.Prim(
    lattice=lattice,
    coordinate_frac=np.array([[0., 0., 0.]]).T,
    occ_dof=[["Si", "Ge"]],          # per-site allowed occupants
    title="SiGe_occ",
)
```

Prim can also be defined from a JSON dict with keys `lattice_vectors`, `basis`, `coordinate_mode`.


### libcasm-composition
- **Import**: `from libcasm.composition import (...)`
- Composition axes (`CompositionConverter`), mol↔parametric coordinate conversion, sublattice composition.

### libcasm-mapping
- **Import**: 
  - `import libcasm.mapping.methods as mapmethods`
  - `import libcasm.mapping.info as mapinfo`
  - `import libcasm.mapping.mapsearch as mapsearch`
- Structure mapping (lattice mapping, atom mapping, full structure mapping) with symmetry. Generates symmetrically equivalent mappings, interpolated structures. Used to map relaxed DFT structures back to ideal configurations.

### libcasm-configuration
- **Import**:
  - `import libcasm.configuration as casmconfig` — `Supercell`, `Configuration`, `ConfigurationSet`
  - `import libcasm.clusterography as casmclust` — `Cluster`, `ClusterSpecs`, cluster comparison/enumeration
  - `from libcasm.enumerate import (...)` — configuration enumeration (by supercell, by property, distinct supercells)
  - `import libcasm.irreps as casmirreps` — irreducible space decomposition, order parameter symmetry
  - `import libcasm.sym_info as sym_info` — symmetry groups (`SymGroup`), factor groups, representations
  - `from libcasm.occ_events import (...)` — occupation events for kinetic Monte Carlo
  - `from libcasm.local_configuration import (...)` — local environment classification

A `Configuration` is specified by a `Supercell` describing the perdiodicity of crystal (determined by a 3×3 transformation matrix of the prim) and an integer occupation array (one entry per site, indexing into the site's `occ_dof` list).

### libcasm-clexulator
- **Import**: `libcasm.clexulator`
- Evaluates cluster expansion basis functions using pre-generated C++ clexulator code. Neighbor list construction, correlation calculation, order parameter evaluation. Does **not** generate basis functions — that is casm-bset's role.
- Requires a C++17 compiler (GCC ≥ 10).

### casm-bset
- **Install**: `pip install casm-bset`
- **Import**: `casm.bset`
- Constructs cluster expansion basis sets: enumerate clusters, build symmetry-adapted functions, generate C++ clexulator source code for use by libcasm-clexulator.

Typical usage via `casm-project`:
```python
bset = project.bset.get(id="default")
bset.make_bspecs()        # build cluster specs (max_length per cluster size)
bset.update()             # generate + compile clexulator C++ code
bset.print_orbits()       # inspect symmetrically equivalent clusters
bset.display_functions()  # show basis function formulas
```

### libcasm-monte
- **Import**: `libcasm.monte`
- Building blocks for Monte Carlo simulations: sampling classes, equilibration/convergence checking, event definitions. Independent of the xtal/configuration stack.

### libcasm-clexmonte
- **Import**: `libcasm.clexmonte`
- Cluster expansion Monte Carlo simulations (canonical, semi-grand canonical, kinetic MC). Combines clexulator, monte, and configuration layers.

### casm-project
- **Import**: `from casm.project import Project`
- High-level project workflow. Reads/writes CASM project files; manages prim, supercells, configurations, calculations, and results. Each sub-workflow is accessed via a command object on the project.

Key command objects:
| Attribute | Type | Purpose |
|---|---|---|
| `project.sym` | `SymCommand` | Print/inspect symmetry |
| `project.enum` | `EnumCommand` | Enumeration workflow |
| `project.bset` | `BsetCommand` | Basis set construction |
| `project.calc` | `CalcCommand` | Calculation setup/status |

### casm-tools
- **Import**: `casm.tools.shared`
- CLI tools: `casm-calc` (set up/run DFT calculations), `casm-map` (structure mapping). ASE integration for I/O.

---

## casm-project Directory Structure

Standardized using casm.project.DirectoryStructure:

```
<project>/
├── .casm/                            # Project settings, prim, etc.
├── calculation_settings/
│   └── calctype.<id>/                # Input files for this calculation type (typically DFT)
│       ├── INCAR, KPOINTS            # Template DFT input files
│       ├── calc.json                 # ASE VASP calculator args
│       ├── meta.json                 # Description of the calculation type
│       └── vasp_relax_greatlakes.sh  # SLURM submit script template
└── enumerations/
    └── enum.<id>/
        ├── scel_set.json           # enumerated supercells
        ├── config_set.json         # enumerated configurations, ConfigurationSet
        ├── config_list.json        # enumerated configuration, list[Configuration]
        ├── config_selection.<name>.json
        └── training_data/
            └── calctype.<id>/<config_name>/
                ├── run.0/          # DFT run
                ├── run.1/          # DFT run
                ├── ...             # ...
                ├── run.final/      # DFT final run
                ├── config.json     # Input CASM configuration
                ├── structure.json  # Input CASM structure
                ├── structure_with_properties.json # Final CASM structure with calculated properties
                └── status.json
```

---

## Common Workflow Patterns

### 1. Project initialization
```python
from casm.project import Project

# From a JSON dict or xtal.Prim object
project = Project(path="MyProject", prim=prim_data)
project.sym.print_factor_group(coord="frac")
```

### 2. Configuration enumeration
```python
enum = project.enum.get(id="main")

# Enumerate symmetrically distinct configurations in supercells of volume 1–4
enum.occ_by_supercell(max=4, min=1)

# For 2D systems, restrict to in-plane supercells
enum.occ_by_supercell(max=10, min=1, dirs="ab", verbose=True)
```

### 3. Iterating configurations
```python
for record in enum.configuration_set:
    name = record.configuration_name          # e.g. "SCEL2_2_1_1_0_0_0/1"
    occ  = record.configuration.occupation    # integer array
    struct = record.configuration.to_structure()
    atom_types = struct.atom_type()
```

### 4. ConfigSelection — filtering and selecting
```python
sel = enum.config_selection   # or load a named selection

for record in sel:            # iterates selected records only
    ...
for record in sel.all:        # iterates all records
    ...

record.select()
record.deselect()
sel.commit()                  # write selection to disk
```

### 5. Basis set construction (via casm-project)
```python
bset = project.bset.get(id="default")
bset.make_bspecs()   # define cluster specs with max_length per orbit size
bset.update()        # generates C++ clexulator and compiles
```

### 6. Correlations and ECI fitting
```python
corr_calc = bset.make_corr_calculator()
correlations = corr_calc.per_unitcell(config_list)  # shape: (n_configs, n_basis)
# fit ECIs externally (e.g. sklearn, lasso) then store in project
```

---

## Common Workflow

1. **Define prim** — lattice, basis sites, allowed occupants/DoF (libcasm-xtal or JSON)
2. **Initialize project** — `Project()` (casm-project)
3. **Enumerate configurations** — `enum.occ_by_supercell(max=N)` (casm-project)
4. **Set up DFT calculations** — write POSCAR/input files, submit jobs (casm-tools / manual)
5. **Map relaxed structures** — map DFT outputs back to configurations (libcasm-mapping / casm-tools)
6. **Build basis set** — `bset.make_bspecs()` + `bset.update()` → compiles clexulator (casm-bset)
7. **Calculate correlations** — `corr_calc.per_unitcell(configs)` (libcasm-clexulator)
8. **Fit ECIs** — linear regression of formation energies vs. correlations
9. **Monte Carlo** — semi-grand canonical or canonical simulations (libcasm-clexmonte)
10. **Analyze** — composition axes, free energy, phase boundaries (libcasm-composition, casm-project)

---
