---
name: casm-dev
description: Developer's guide for CASM (libcasm packages and casm packages)
---

## General notes

- Do not include "Co-Authored-By" in commit messages.
- Keep commit messages concise: a title plus at most a short paragraph. Detailed rationale/root-cause goes in code comments or the chat response, not the commit body.
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

libcasm (C++ extension) packages:

- Add methods to pybind11-bound classes (Configuration, Prim, Supercell, etc.) in the C++ binding source (`python/src/*.cpp`, via `.def()` / `.def_property_readonly()`), not by monkey-patching in `_methods.py` or `__init__.py`.

Documentation guidelines:

- For writing voice, phrasing, plain-language rules, and preserving user prose, see `documentation-style.md` in this skill directory.
- For CASM vocabulary, notation, and math symbols, use the `casm-notation` skill.
- The numpydoc/formatting conventions below apply on top of those.

- **Use numpydoc style** — Parameters, Returns, Raises, Notes sections with section name on one line and a dashes underline (-------).
- **Show array shapes in type annotations** — Write numpy.ndarray[numpy.float64[3, n]] not just np.ndarray; include shape hints like shape=(3,3) in the parameter description line when the full type is too verbose.
- **Lead with a one-line imperative summary** — First line should be a single sentence fragment starting with a verb (e.g., "Convert fractional coordinates to Cartesian coordinates"), followed by a blank line before any extended description.
- **Use :class:, :func:, :data: cross-reference roles** — Always link to related types and functions (e.g., :func:`make_primitive`, :class:`~libcasm.xtal.Lattice`). Use the ~ prefix to suppress the module path in rendered link text.
- **Use italics and code block consistently** - Refer to programs, packages, files, paths, arguments, and variables using italics. For code snippets, commands, command line options, and defining environment variables, use inline code (``x=3``) or  code blocks (".. code-block:: Python"). 
- **Use ".. rubric::" for informal sections** — Non-numpydoc sections like ".. rubric:: Constructor", ".. rubric:: Special Methods", or ".. rubric:: Example" usage keep the structure consistent with Sphinx rendering without triggering numpydoc section parsing. To document any special members of a class, such as comparison operators (*<*, *<=*, *>*, *>=*, etc.) or arithmetic operators (`*`, `*=`, `+`, `+=`, `-`, `-=`, etc.), use `.. rubric:: Special Methods` to create a section in a class docstring.
- **Put math in ":math:" / ".. math::" directives** — Inline math uses :math:`\vec{r}`; block equations use ".. math::".Never use plain text or ASCII art for formulas.
- **Show defaults in the parameter type field** — Write int, default=1 or Optional[np.ndarray] = None on the type line, not buried in the description   
  prose.
- **Include ".. code-block:: Python" examples for non-trivial behavior** — When the semantics are best explained by showing the equivalent Python expression (e.g., r_after = op.matrix() @ r_before + op.translation()), put it in a code block right after the summary, before Parameters.
- **Put "Equivalent to" cross-references in the extended summary** — When a function wraps or dispatches to another, add a line like "Equivalent to :func:make_primitive_prim when obj is a :class:Prim." immediately after the summary, with no section header (per existing memory rule).
- **Document Optional parameters using "Optional[Type] = None" in the type field** - Mark other optional-with-default parameters as "Type, default=value**". Distinguish "can be None" from "has a non-None default" explicitly.


Testing:

- Add Python tests for *casm.<subpackage>* in *python/tests/<subpackage>*, using pytest.
- If data files are needed for testing, they can be placed in *python/tests/<subpackage>/data/*.
- To access data files, use the *shared_datadir* fixture available from the `pytest-datadir <https://pypi.org/project/pytest-datadir/>`_ plugin.
- For tests reading and writing files, use temporary directories created by the `tmpdir` and `tmpdir_factory` pytest fixtures.
- For tests that involve an expensive setup process, such as compiling Clexulators, a session-length shared datadir can be constructed once and re-used as done in `CASMcode_clexulator/python/tests/clexulator/conftest.py`.
- Expensive tests can also be set to run optionally using flags as demonstrated in CASMcode_clexulator.


Adding dependencies:

- Do not add dependencies unless the user specifically requests or approves it
- If a dependency is only required for building, testing, or documentation, add it to *build_requirements.txt*, *test_requirements.txt*, or *doc_requirements.txt* instead of adding it as dependency.


Formatting and linting:

- For C++ extensions, use clang-format via the `stylize.sh` script
- For `libcasm` packages Python, use `black python` and `ruff check --fix python`
- for `casm` packages, use `black casm tests` for ``ruff check --fix casm tests`

For versioning procedures, use the `/casm-version` skill.
For release procedures, use the `/casm-release` skill.

Release and workflow scripts live in `CASMcode_global/dev/`:

- `dev/release.py` — full release workflow (run from package root as `python ../CASMcode_global/dev/release.py`)
- `dev/download_release.py` — download build artifacts from GitHub Actions
- `dev/update_workflow_versions.py` — update CASM dependency versions in `.github/workflows/*.yml`
