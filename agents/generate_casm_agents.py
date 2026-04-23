"""Generate AGENTS.md (or CLAUDE.md) for each casm-X repo from a shared template.

Run from CASMcode_tools/agents/:
    python generate_agents.py
    python generate_agents.py --filename CLAUDE.md

Or from anywhere:
    python ../CASMcode_tools/agents/generate_agents.py
"""

import argparse
import pathlib

REPOS = [
    # (repo_name, package_name, docs_name)
    ("CASMcode_bset", "casm-bset", "bset"),
    ("CASMcode_project", "casm-project", "project"),
    ("CASMcode_tools", "casm-tools", "tools"),
]

parser = argparse.ArgumentParser()
parser.add_argument("--filename", default="AGENTS.md", choices=["AGENTS.md", "CLAUDE.md"])
args = parser.parse_args()

agents_dir = pathlib.Path(__file__).parent
template_path = agents_dir / "AGENTS_casm.md.template"
root = agents_dir.parent.parent  # .../CASMcode_tools/agents/ -> .../

template = template_path.read_text()

for repo_name, package_name, docs_name in REPOS:
    content = template.format(
        repo_name=repo_name,
        package_name=package_name,
        docs_name=docs_name,
    )
    out_path = root / repo_name / args.filename
    out_path.write_text(content)
    print(f"wrote {out_path}")
