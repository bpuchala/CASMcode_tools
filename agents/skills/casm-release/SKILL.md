---
name: casm-release
description: Guides the release process for a CASM package.
---

## Usage

Invoke with `/casm-release` when releasing a CASM package. Run from the package root after all changes are committed and GitHub Actions have passed.

A helper script automates the steps below:

    python ../CASMcode_global/dev/release.py <version> --dev-branch 2.X --pydocs ../CASMcode_pydocs
    # Use --start-step N or --steps LIST to resume

Replace `<package>`, `<Package>`, `<version>`, `<Package_>` as appropriate in the commands below.

## Step 1 — Build wheels and upload to PyPI

```bash
python ../CASMcode_global/dev/download_release.py <version>
python label_wheels.py <version>
python -m twine upload dist/<version>/*
```

## Step 2 — Install from PyPI and run tests

```bash
pip install --upgrade <Package>
pytest -rsap -x python/tests
```

## Step 3 — Create and push version tag

```bash
git branch --show-current  # should be <version>
git tag v<version>
git push public v<version>
```

## Step 4 — Create GitHub release

Use CHANGELOG.md entry as release notes.

```bash
gh release create v<version> --title "v<version>" --notes "<notes>" dist/<version>/<Package_>-<version>.tar.gz
```

## Step 5 — Merge tag into main, then main into 2.X

```bash
git fetch --all --tags
# Verify: git log --oneline -5 public/main && public/2.X && v<version>
git checkout main && git pull && git merge --ff-only v<version> && git push public main
git checkout 2.X && git pull public 2.X && git merge main --no-edit && git push public 2.X
```

## Step 6 — Build docs

```bash
rm -rf python/doc/reference/libcasm/_autosummary
sphinx-build -b html python/doc $LIBCASM_PYDOCS/<package>/2.0
```

`<package>`: `clexmonte`, `clexulator`, `composition`, `configuration`, `global`, `mapping`, `monte`, or `xtal`

## Step 7 — Update overview page

Edit `CASMcode_pydocs/src/overview/latest/index.rst` to update the version, then:

```bash
cd ../CASMcode_pydocs
sphinx-build -b html src/overview/latest docs/overview/latest/
```

## Step 8 — Commit and push docs

```bash
git add docs/libcasm/<package>/2.0/ docs/overview/latest/ src/overview/latest/index.rst
git commit -m "<Package> v<version>"
git push
```
