"""Implements ``casm-convert ...``"""

import argparse
import os
import pathlib
import sys


# <-- max width = 80 characters                            --> #
################################################################
def print_desc():
    desc = """
# The `casm-convert` command:

## Method

Convert a crystal structure file from one format to another.

The input and output formats are inferred from the file
extensions when possible. To override the inferred formats, or
when files have no extension, use the `--input-format` and
`--output-format` options.

Supported formats:

- 'vasp': VASP POSCAR format. Used when the file has no
  extension or a '.vasp' extension.
- 'casm': CASM Structure JSON format. Used when the file has a
  '.json' or '.casm' extension.
- Any format recognized by `ase.io` if ASE is installed.


## Parameters

Positional arguments:

input: pathlib.Path
    Input structure file.
output: pathlib.Path
    Output structure file.

Options:

-i, --input-format: Optional[str]=None
    Format for reading the input file. If not specified, the
    format is inferred from the file extension. Supported
    formats include 'vasp', 'casm', and any format recognized
    by `ase.io.read` if ASE is installed.
-o, --output-format: Optional[str]=None
    Format for writing the output file. If not specified, the
    format is inferred from the file extension. Supported
    formats include 'vasp', 'casm', and any format recognized
    by `ase.io.write` if ASE is installed.
-f, --force: bool=False
    If given, overwrite the output file if it already exists.
"""
    print(desc)


def run_convert(args):
    """Implements ``casm-convert ...``

    Parameters
    ----------
    args : argparse.Namespace
        The parsed arguments from the command line.

    Returns
    -------
    code: int
        A return code indicating success (0) or failure (non-zero).

    """
    from casm.tools.shared.structure_io import read_structure, write_structure

    if args.desc:
        print_desc()
        return 0

    try:
        structure = read_structure(path=args.input, format=args.input_format)
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return 1

    try:
        write_structure(
            path=args.output,
            casm_structure=structure,
            format=args.output_format,
            force=args.force,
        )
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return 1

    return 0


def make_parser():
    """Make the ``casm-convert ...`` argument parser

    Returns
    -------
    parser: argparse.ArgumentParser
        The argument parser for the `casm-convert` program.

    """
    parser = argparse.ArgumentParser(
        description="Convert a crystal structure file from one format to another.",
    )

    ### Positional arguments
    positional = parser.add_argument_group("Positional arguments")
    positional.add_argument(
        "input",
        type=pathlib.Path,
        help="Input structure file.",
    )
    positional.add_argument(
        "output",
        type=pathlib.Path,
        help="Output structure file.",
    )

    ### Options
    options = parser.add_argument_group("Options")
    options.add_argument(
        "-i",
        "--input-format",
        type=str,
        default=None,
        help=(
            "Input file format (default= inferred from file extension ). "
            "Supported formats include 'vasp', 'casm', and any format "
            "recognized by ase.io.read if ASE is installed."
        ),
    )
    options.add_argument(
        "-o",
        "--output-format",
        type=str,
        default=None,
        help=(
            "Output file format (default= inferred from file extension ). "
            "Supported formats include 'vasp', 'casm', and any format "
            "recognized by ase.io.write if ASE is installed."
        ),
    )
    options.add_argument(
        "-f",
        "--force",
        action="store_true",
        default=False,
        help="Overwrite output file if it already exists.",
    )

    ### Other options
    other = parser.add_argument_group("Other options")
    other.add_argument(
        "--desc",
        action="store_true",
        help="Print an extended description of the method and parameters.",
    )

    parser.set_defaults(func=run_convert)
    return parser


def main(argv=None, working_dir=None):
    """Implements ``casm-convert ...``

    Parameters
    ----------
    argv : list of str, optional
        The command line arguments to parse. If None, uses `sys.argv`.
    working_dir : str, optional
        The working directory to use. If None, uses the current working directory.

    Returns
    -------
    code: int
        A return code indicating success (0) or failure (non-zero).

    """
    from casm.tools.shared import contexts

    if argv is None:
        argv = sys.argv
    if working_dir is None:
        working_dir = os.getcwd()

    parser = make_parser()

    if "--desc" in argv:
        print_desc()
        return 0

    if len(argv) < 2:
        parser.print_help()
        return 1

    args = parser.parse_args(argv[1:])

    with contexts.working_dir(working_dir):
        code = args.func(args)

    return code
