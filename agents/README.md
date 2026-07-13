# agents/

Agent and AI assistant configuration for CASM. These files are used to generate `AGENTS.md` and `CLAUDE.md` files for the various CASMcode repositories, assuming they are all in the same parent directory.

## Files

| File | Purpose |
|------|---------|
| `AGENTS_libcasm.md.template` | Template for C++ extension package repos |
| `AGENTS_casm.md.template` | Template for pure Python package repos |
| `generate_libcasm_agents.py` | Generate `AGENTS.md` or `CLAUDE.md` for libcasm repos |
| `generate_casm_agents.py` | Generate `AGENTS.md` or `CLAUDE.md` for casm repos |
| `install_skills.py` | Install `skills/` into an AI tool's skills directory |
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

Install the skills with `install_skills.py` (invoke with `/casm`, `/casm-dev`, `/casm-version`, `/casm-release`). Every `skills/<name>/` folder containing a `SKILL.md` is copied to the destination, replacing any existing copy.

Run from `CASMcode_tools/agents/` (or anywhere using the full path):

```bash
# Shortcut options for common CLI AI tools
python install_skills.py --claude     # ~/.claude/skills
python install_skills.py --codex      # ~/.codex/skills
python install_skills.py --gemini     # ~/.gemini/skills
python install_skills.py --opencode   # ~/.config/opencode/skills

# Any other tool / arbitrary directory
python install_skills.py --dest ~/.config/other-tool/skills
```

A destination is required (there is no default).

## Notes

- `AGENTS.md` is also read by some other agents (e.g. OpenAI Codex, Google Jules).
- These files are maintained in `CASMcode_tools` as the canonical source.
