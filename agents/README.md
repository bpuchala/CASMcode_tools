# agents/

Agent and AI assistant configuration for CASM. These files are used to generate `AGENTS.md` and `CLAUDE.md` files for the various CASMcode repositories, assuming they are all in the same parent directory.

## Files

| File | Purpose |
|------|---------|
| `AGENTS_libcasm.md.template` | Template for C++ extension package repos |
| `AGENTS_casm.md.template` | Template for pure Python package repos |
| `generate_libcasm_agents.py` | Generate `AGENTS.md` or `CLAUDE.md` for libcasm repos |
| `generate_casm_agents.py` | Generate `AGENTS.md` or `CLAUDE.md` for casm repos |
| `skills/casm/` | On-demand skill: CASM concepts, APIs, and workflow patterns |
| `skills/casm-dev/` | On-demand skill: CASM developer guide |
| `skills/casm-version/` | On-demand skill: versioning checklist |
| `skills/casm-release/` | On-demand skill: release procedure |

Release and workflow scripts live in `dev/` of the CASMcode_global repository:

| Script | Purpose |
|--------|---------|
| `dev/release.py` | Full release workflow (run from package root as `python ../CASMcode_global/dev/release.py`) |
| `dev/download_release.py` | Download build artifacts from GitHub Actions |
| `dev/update_workflow_versions.py` | Update CASM dependency versions in `.github/workflows/*.yml` |

## Generating repo files

Run from `CASMcode_tools/agents/` (or anywhere using the full path):

```bash
# Generate CLAUDE.md for all repos (Claude Code)
python generate_libcasm_agents.py --filename CLAUDE.md
python generate_casm_agents.py --filename CLAUDE.md

# Generate AGENTS.md for all repos (Codex, Jules, etc.)
python generate_libcasm_agents.py
python generate_casm_agents.py
```

## Skills setup

Copy skills to `~/.claude/skills/` (invoke with `/casm`, `/casm-dev`, `/casm-version`, `/casm-release`):

```bash
mkdir -p ~/.claude/skills
cp -r agents/skills/casm ~/.claude/skills/casm
cp -r agents/skills/casm-dev ~/.claude/skills/casm-dev
cp -r agents/skills/casm-version ~/.claude/skills/casm-version
cp -r agents/skills/casm-release ~/.claude/skills/casm-release
```

## Notes

- `AGENTS.md` is also read by some other agents (e.g. OpenAI Codex, Google Jules).
- These files are maintained in `CASMcode_tools` as the canonical source.
