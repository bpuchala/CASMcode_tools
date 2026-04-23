---
name: casm-dev
description: Developer's guide for CASM (libcasm packages and casm packages)
---

## General notes

- Do not include "Co-Authored-By" in commit messages.
- Assume `CASM_PREFIX` is already set in the environment.
- Use MCP server `superchivo` to search Python package docs.

---

## CASM packages

### C++ extension packages

- `CASMcode_global` (libcasm-global) — Constants, definitions, and generically useful tools used by CASM
- `CASMcode_crystallography` (libcasm-xtal) — CASM crystallography
- `CASMcode_composition` (libcasm-composition) — CASM composition axes, conversions, and calculations
- `CASMcode_mapping` (libcasm-mapping) — CASM structure mapping
- `CASMcode_clexulator` (libcasm-clexulator) — CASM clexulator
- `CASMcode_configuration` (libcasm-configuration) — CASM configuration comparison and enumeration
- `CASMcode_monte` (libcasm-monte) — CASM building blocks for Monte Carlo simulations
- `CASMcode_clexmonte` (libcasm-clexmonte) — CASM cluster expansion Monte Carlo

### Pure Python packages

- `CASMcode_bset` (casm-bset) — CASM cluster expansion basis set generation
- `CASMcode_project` (casm-project) — CASM project interface
- `CASMcode_tools` (casm-tools) — CASM command line tools

---

## Package structure

### C++ extension packages

- `.github/` — GitHub workflow files
- `cmake/` — CMake modules
- `doc/` — Doxygen documentation source
- `include/` — C++ headers
- `python/doc/` — Sphinx docs; `python/libcasm/` — Python source; `python/src/` — pybind11 bindings; `python/tests/` — tests
- `src/` — C++ source; `tests/unit/` — C++ unit tests
- Key files: `CMakeLists.txt.in`, `make_CMakeLists.py`, `pyproject.toml`, `setup.py`, `stylize.sh`, `CHANGELOG.md`
- `python/pyproject.toml` and `python/setup.py` are **editable-mode only**

### Pure Python packages

- `.github/` — workflows
- `casm/` — Python source
- `doc/` — Sphinx docs
- `tests/` — Python tests
- Key files: `pyproject.toml`, `setup.py`, `CHANGELOG.md`, `requirements.txt`

---

## CASM Package Dependencies

```
libcasm-global           [none]
libcasm-xtal             [global]
libcasm-composition      [global]
libcasm-mapping          [global, xtal]
libcasm-clexulator       [global, xtal]
libcasm-configuration    [global, xtal, clexulator]
casm-bset                [configuration]
libcasm-monte            [global, xtal, composition]
libcasm-clexmonte        [global, xtal, composition, clexulator, configuration, monte]
casm-tools               [xtal, mapping, configuration]
casm-project             [mapping, clexmonte, bset, tools]
```

---

## Development

Each repo's `agents/AGENTS.md` contains repo-specific dev workflow commands.

For versioning procedures, use the `/casm-version` skill.
For release procedures, use the `/casm-release` skill.

Release and workflow scripts live in `CASMcode_global/dev/`:

- `dev/release.py` — full release workflow (run from package root as `python ../CASMcode_global/dev/release.py`)
- `dev/download_release.py` — download build artifacts from GitHub Actions
- `dev/update_workflow_versions.py` — update CASM dependency versions in `.github/workflows/*.yml`
