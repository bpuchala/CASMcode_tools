# Changelog

All notable changes to `casm-tools` will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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

