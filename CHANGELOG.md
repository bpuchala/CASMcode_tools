# Changelog

All notable changes to `casm-tools` will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `casm-convert`: New CLI tool for converting crystal structure files between
  formats. Formats are inferred from file extensions (VASP POSCAR, CASM JSON,
  or any format supported by ASE). Supports `-i`/`--input-format`,
  `-o`/`--output-format`, and `-f`/`--force` options.
- `casm.tools.map.methods.make_chain_info`: Added an `is_primitive` option, to skip
  making the structures primitive when the chain is already primitive, as from
  `make_primitive_chain`.

### Changed

- `casm.tools.map.methods.make_chain_info`: Determine the volume and symmetry of each
  structure in a chain from the primitive equivalent structure, using
  `libcasm.xtal.make_primitive_structure`, so that `volume`, `factor_group_size`, and
  `spacegroup_type` describe the crystal and do not depend on which superstructure the
  interpolated structure is expressed in. An equivalent chain of primitive structures,
  as from `make_primitive_chain`, now gives the same results. Previously the structure
  as given was used, which reported the symmetry restricted to that superstructure. For
  example, `C2/m` with a 32 element factor group for a bcc structure in a det(T)=8
  supercell, rather than `Im-3m` with the 48 element prim factor group.
- `casm.tools.map.methods`: Use `libcasm.xtal.make_primitive_structure` instead of a
  local helper that made structures primitive via `libcasm.xtal.Prim`.
- Update libcasm-xtal dependency to `>=3.0a1`
- Update libcasm-mapping dependency to `>=3.0a1`
- Update libcasm-configuration dependency to `>=3.0a1`

### Fixed

- `casm.tools.map.methods.make_supercell_info`: Fixed the supercell `volume`, which
  truncated instead of rounding the determinant of the transformation matrix. For
  example, a det(T)=8 supercell was reported as volume 7, because `numpy.linalg.det`
  evaluates the determinant of that integer matrix as 7.999999999999998.


## [2.0a3] - 2026-02-26

### Changed

- Set requires-python to ">=3.10,<3.15" to match the wheels being built for CASM C++ extensions and distributed on PyPI.
- Changed libcasm package dependencies to latest versions, which will avoid potential compatibility issues between older versions of libcasm packages and the planned libcasm-xtal>=3. 


## [2.0a2] - 2024-08-07

### Changed

- Use `libcasm.mapping.mapsearch.make_atom_to_site_cost_future` for making 
  atom to site cost functions.


## [2.0a1] - 2025-08-04

This release creates the casm-tools package, which provides pure Python CLI tools and 
helper functions. This includes:

- casm-calc: Setup, run, and report results of structure calculations
- casm-map: Structure mapping and import
- casm.tools.shard: Helper functions for I/O, integrating with ASE, and context 
  managers.

