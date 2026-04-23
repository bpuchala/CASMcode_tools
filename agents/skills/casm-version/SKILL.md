---
name: casm-version
description: Guides updating the version across all required files in a CASM package.
---

## Usage

Invoke with `/casm-version` when preparing a new version in any CASM package.

## Versioning checklist

Update version in all that apply (MAJOR.MINOR.PATCH only where noted):

1. `CMakeLists.txt.in` (M.N.P only)
2. `tests/CMakeLists.txt.in` (M.N.P only)
3. `tests/unit/<package>/version_test.cpp` (if exists)
4. `pyproject.toml`
5. `setup.py`
6. `python/tests/casmglobal/test_casmglobal.py` (if exists)
7. `python/doc/conf.py`
8. `python/setup.py`
9. `doc/doxygen_config` (PROJECT_NUMBER, M.N.P only)
10. `CHANGELOG.md` (add entry, set date)

After updating `CMakeLists.txt.in`, run:

    python make_CMakeLists.py

---

## Updating CASM dependency version constraints

Update all that apply:

- `pyproject.toml` — both `[build-system]` requires and `[project]` dependencies
- `build_requirements.txt`
- `python/pyproject.toml` — editable-mode build requires

---

## GitHub Actions workflow cache keys

Linux only; macOS uses pip and needs no cache updates.

- `test-linux-dependencies.yml` — update cache key **and** `ref:` checkout tag
- `test-linux-build.yml`, `test-linux.yml`, `test-linux-cxx-only.yml` — update cache key only
