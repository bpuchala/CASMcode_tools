"""Install CASM skills into an AI assistant's skills directory.

Each subdirectory of ``skills/`` holding a ``SKILL.md`` is copied to the
destination. This works for any AI tool that loads skills from a directory
of ``<name>/SKILL.md`` folders; use a shortcut option for a common CLI tool
or ``--dest`` for an arbitrary directory.

Run from CASMcode_tools/agents/:
    python install_skills.py --claude
    python install_skills.py --dest ~/.config/other-tool/skills

Or from anywhere:
    python ../CASMcode_tools/agents/install_skills.py --claude

Use ``--remove`` to uninstall the skills from the destination instead:
    python install_skills.py --claude --remove
"""

import argparse
import pathlib
import shutil

# (option_name, destination, tool_name)
SHORTCUTS = [
    ("claude", "~/.claude/skills", "Claude Code"),
    ("codex", "~/.codex/skills", "OpenAI Codex CLI"),
    ("gemini", "~/.gemini/skills", "Google Gemini CLI"),
    ("opencode", "~/.config/opencode/skills", "opencode"),
]

parser = argparse.ArgumentParser(description=__doc__)
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument(
    "--dest",
    help="destination skills directory",
)
for option, destination, tool_name in SHORTCUTS:
    group.add_argument(
        f"--{option}",
        dest="dest",
        action="store_const",
        const=destination,
        help=f"install to the {tool_name} skills directory ({destination})",
    )
parser.add_argument(
    "--remove",
    action="store_true",
    help="remove the skills from the destination instead of installing them",
)
args = parser.parse_args()

agents_dir = pathlib.Path(__file__).parent
skills_dir = agents_dir / "skills"
dest_dir = pathlib.Path(args.dest).expanduser()

if not args.remove:
    dest_dir.mkdir(parents=True, exist_ok=True)

for skill_path in sorted(skills_dir.iterdir()):
    if not (skill_path / "SKILL.md").is_file():
        continue
    out_path = dest_dir / skill_path.name
    if args.remove:
        if out_path.exists():
            shutil.rmtree(out_path)
            print(f"removed {out_path}")
        else:
            print(f"not installed {out_path}")
        continue
    if out_path.exists():
        shutil.rmtree(out_path)
    shutil.copytree(skill_path, out_path)
    print(f"installed {skill_path.name} -> {out_path}")
