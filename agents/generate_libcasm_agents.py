"""Generate AGENTS.md (or CLAUDE.md) for each libcasm-X repo from a shared template.

Run from CASMcode_tools/agents/:
    python generate_libcasm_agents.py
    python generate_libcasm_agents.py --filename CLAUDE.md

Or from anywhere:
    python ../CASMcode_tools/agents/generate_libcasm_agents.py
"""

import argparse
import pathlib

REPOS = [
    # (repo_name, package_name, docs_name)
    ("CASMcode_global", "libcasm-global", "global"),
    ("CASMcode_crystallography", "libcasm-xtal", "xtal"),
    ("CASMcode_composition", "libcasm-composition", "composition"),
    ("CASMcode_mapping", "libcasm-mapping", "mapping"),
    ("CASMcode_clexulator", "libcasm-clexulator", "clexulator"),
    ("CASMcode_configuration", "libcasm-configuration", "configuration"),
    ("CASMcode_monte", "libcasm-monte", "monte"),
    ("CASMcode_clexmonte", "libcasm-clexmonte", "clexmonte"),
]

parser = argparse.ArgumentParser()
parser.add_argument(
    "--filename", default="AGENTS.md", choices=["AGENTS.md", "CLAUDE.md"]
)
args = parser.parse_args()

agents_dir = pathlib.Path(__file__).parent
template_path = agents_dir / "AGENTS_libcasm.md.template"
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
