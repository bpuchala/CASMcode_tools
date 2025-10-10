"""Implements ``casm-calc submit ...``"""

# Note:
# The command `casm-calc submit` is implemented by `run_submit` which requires
# `casm.project` to be installed. This is a bit odd, because `casm.tools` is a
# dependency of `casm.project`. At some point `casm-calc submit` will be refactored
# to be a plugin provided by `casm.project`.


def run_submit(args):
    """Implements ``casm-calc submit ...``

    Parameters
    ----------
    args : argparse.Namespace
        The parsed arguments from the command line. Uses:

        - `args.enum`: str, Enumeration ID
        - `args.calctype`: str, A calctype ID to specify the calculation type.
        - `args.clex`: str, Cluster expansion ID to specify the calctype. This is the
          fallback, used if `args.calctype` is not specified. If this is also not
          specified, the default cluster expansion is used.
        - `args.selection`: str, optional, Name of a selection in the enumeration with
          selected configurations to submit calculations for. Select all configurations
          if not specified.
        - `args.dry_run`: If set, prints configurations that would be submitted and
          their current calc status, but does not submit.
        - `args.cancel`: If set, calculations with status="started" or "submitted" are
          canceled. No jobs are submitted.

    Returns
    -------
    code: int
        A return code indicating success (0) or failure (non-zero).

    """
    from .status import (
        _validate_casm_project,
        list_options,
        make_config_selection,
    )

    _validate_casm_project()

    from casm.project import Project, project_path
    from casm.tools.shared.json_io import safe_dump

    # --- Load project ---
    path = project_path()
    if path is None:
        print("No CASM project found.")
        return 1
    project = Project(path=path)

    # --- list options ---
    if args.list:
        list_options(
            project=project,
            args=args,
        )
        return 0
    elif args.enum is None:
        print("No enumeration ID provided. Use --list to see available options.")
        return 1

    # --- Variables ---
    if args.dry_run:
        submit_jobs = False
    else:
        submit_jobs = True

    # --- Select configurations ---
    config_selection = make_config_selection(
        project=project,
        args=args,
    )
    if isinstance(config_selection, int):
        # An error occurred in make_config_selection
        return config_selection

    # --- Submit jobs or print status ---
    dry_run_msg = "(dry-run)"
    print(f"{'Name':36}{'Status':12}{'Job ID':12}{'Message':18}")
    print("-" * 78)

    n_jobs_not_ready = 0

    for record in config_selection:

        if args.cancel:
            if record.calc_status == "started" or record.calc_status == "submitted":
                # --- Update status.json ---
                status_data = record.calc_status_data
                status_data["status"] = "canceled"
                safe_dump(
                    data=status_data,
                    path=record.calc_dir / "status.json",
                    quiet=True,
                    force=True,
                )
                jobid = record.calc_jobid
                msg = "canceled"

                # --- Cancel job ---
                completed_processes = record.run_subprocess(
                    args=["scancel", record.calc_jobid],
                    capture_output=True,
                    text=True,
                )

                if completed_processes[0].returncode != 0:
                    print(
                        f"Error canceling job for {record.name}: "
                        f"{completed_processes[0].stderr}"
                    )
                    return 1

                # --- Print status ---
                print(f"{record.name:36}{record.calc_status:12}{jobid:12}{msg:18}")
            continue

        if not submit_jobs:
            if record.calc_status != "setup" or record.calc_jobid != "none":
                n_jobs_not_ready += 1
                msg = "skipping"
            else:
                msg = dry_run_msg
            jobid = record.calc_jobid

            # --- Printstatus ---
            print(f"{record.name:36}{record.calc_status:12}{jobid:12}{msg:18}")
        elif record.calc_status != "setup" or record.calc_jobid != "none":
            n_jobs_not_ready += 1
            msg = "skipping"
            jobid = record.calc_jobid

            # --- Print status ---
            print(f"{record.name:36}{record.calc_status:12}{jobid:12}{msg:18}")

        else:

            # --- Submit job w/ hold ---
            completed_processes = record.run_subprocess(
                args=["sbatch", "--hold", "submit.sh"],
                capture_output=True,
                text=True,
            )

            if completed_processes[0].returncode != 0:
                print(
                    f"Error submitting job for {record.name}: "
                    f"{completed_processes[0].stderr}"
                )
                return 1

            # --- Extract job ID from output ---
            stdout = completed_processes[0].stdout
            if "Submitted batch job" in stdout:
                jobid = stdout.strip().split()[-1]

                # --- Update status.json ---
                status_data = record.calc_status_data
                status_data["status"] = "submitted"
                status_data["jobid"] = jobid
                safe_dump(
                    data=status_data,
                    path=record.calc_dir / "status.json",
                    quiet=True,
                    force=True,
                )

                # --- Release job ---
                record.run_subprocess(
                    args=["scontrol", "release", f"{jobid}"],
                    capture_output=True,
                    text=True,
                )
            else:
                print("~" * 40)
                print(f"Error submitting job for {record.name}: ")
                print("---")
                print(stdout)
                print("---")
                print(completed_processes[0].stderr)
                print("~" * 40)

            # --- Print status ---
            msg = ""
            print(f"{record.name:36}{record.calc_status:12}{jobid:12}{msg:18}")

    if n_jobs_not_ready:
        print()
        print(
            f"Warning: {n_jobs_not_ready} jobs have status != 'setup' or "
            "jobid != 'none', so they were not submitted."
        )
        print()

    return 0


def run_submitd(args):
    """Implements ``casm-calc submitd ...``

    Acts on all calculation subdirectories of a specified top level directory.
    Calculation subdirectories are identified by the presence of a `status.json` file.
    Each calculation directory is expected to have a `submit.sh` script that can be
    submitted with `sbatch`.

    Parameters
    ----------
    args : argparse.Namespace
        The parsed arguments from the command line. Uses:

        - `args.dir`: pathlib.Path, The top level directory. Subdirectories with a
          `status.json` file are treated as calculation directories.
        - `args.dry_run`: If set, prints configurations that would be submitted and
          their current calc status, but does not submit.
        - `args.cancel`: If set, calculations with status="started" or "submitted" are
          canceled. No jobs are submitted.
        - `args.name_parts`: int, The number of path parts to include in the name.
          For example, if the calculation directory is `a/b/c/d` and `name_parts`
          is 2, the name will be `c/d`. If `name_parts` is greater than the number
          of parts in the path, the entire path is used.

    Returns
    -------
    code: int
        0 for success, non-zero for failure
    """
    import os

    from casm.tools.shared.calc_dir import CalcDir

    top_dir = args.dir
    if not top_dir.is_dir():
        print(f"Directory not found: {args.dir}")
        return 1

    calc_dirs = [
        CalcDir(path=p.parent, name_parts=args.name_parts)
        for p in top_dir.rglob("status.json")
    ]
    if not calc_dirs:
        print("No calculation directories found.")
        return 1

    print(f"{'Name':36}{'Status':12}{'Job ID':12}{'Message':18}")
    print("-" * 78)
    n_jobs_not_ready = 0

    for calc_dir in calc_dirs:

        name = calc_dir.name
        calc_status = calc_dir.status
        jobid = calc_dir.jobid

        if args.cancel:
            if calc_status in ("started", "submitted"):
                calc_dir.status_data["status"] = "canceled"
                calc_dir.commit_status(quiet=True)
                # Attempt to cancel job
                completed_process = os.system(f"scancel {jobid}")
                msg = "canceled" if completed_process == 0 else "cancel error"
                print(f"{name:36}{calc_status:12}{jobid:12}{msg:18}")
            continue

        if args.dry_run:
            if calc_status != "setup" or jobid != "none":
                n_jobs_not_ready += 1
                msg = "skipping"
            else:
                msg = "(dry-run)"
            print(f"{name:36}{calc_status:12}{jobid:12}{msg:18}")
        elif calc_status != "setup" or jobid != "none":
            n_jobs_not_ready += 1
            msg = "skipping"
            print(f"{name:36}{calc_status:12}{jobid:12}{msg:18}")
        else:
            # --- Submit job w/ hold ---
            import subprocess

            completed_process = subprocess.run(
                args=["sbatch", "--hold", "submit.sh"],
                cwd=calc_dir.path,
                capture_output=True,
                text=True,
            )

            if completed_process.returncode != 0:
                print(
                    f"Error submitting job for {name}: " f"{completed_process.stderr}"
                )
                return 1

            # --- Extract job ID from output ---
            stdout = completed_process.stdout
            if "Submitted batch job" in stdout:
                jobid = stdout.strip().split()[-1]

                # --- Update status.json ---
                calc_dir.status_data["status"] = "submitted"
                calc_dir.status_data["jobid"] = jobid
                calc_dir.commit_status(quiet=True)

                # --- Release job ---
                subprocess.run(
                    args=["scontrol", "release", f"{jobid}"],
                    cwd=calc_dir.path,
                    capture_output=True,
                    text=True,
                )
            else:
                print("~" * 40)
                print(f"Error submitting job for {name}: ")
                print("---")
                print(stdout)
                print("---")
                print(completed_process.stderr)
                print("~" * 40)

            # --- Print status ---
            msg = ""
            print(f"{name:36}{calc_dir.status:12}{jobid:12}{msg:18}")

    if n_jobs_not_ready:
        print()
        print(
            f"Warning: {n_jobs_not_ready} jobs have status != 'setup' or "
            "jobid != 'none', so they were not submitted."
        )
        print()

    return 0


################################################################################


def print_desc(argv=None):
    print("No extended description available.")


def make_submit_subparser(c):
    """Constructs the ``casm-calc submit ...`` argument parser, and attaches the methods
    for running the subcommands.

    Parameters
    ----------
    c: argparse._SubParsersAction
        The output from ``parser.add_subparsers`` to which ``casm-calc submit``
        arguments are added.

    Returns
    -------
    code: int
        A return code indicating success (0) or failure (non-zero).

    """
    submit = c.add_parser(
        "submit",
        help="Submit CASM project calculations",
        description="Submit CASM project calculations",
    )

    ### casm-calc submit ....
    submit.set_defaults(func=run_submit)

    submit.add_argument(
        "enum",
        type=str,
        help=("Enumeration ID"),
        nargs="?",
    )
    submit.add_argument(
        "-c",
        "--selection",
        type=str,
        help=(
            "Name of a selection in the enumeration with selected configurations to "
            "submit calculations for. Select all configurations if not specified. "
        ),
    )
    submit.add_argument(
        "--calctype",
        type=str,
        help=(
            "A calctype ID, used to indicate which calctype to submit"
            "calculations for. Uses --clex if not specified. "
        ),
    )
    submit.add_argument(
        "--clex",
        type=str,
        help=(
            "A cluster expansion key, used to indicate which calctype to submit "
            "calculations for. Uses the default clex if not specified. "
        ),
    )
    submit.add_argument(
        "--dry-run",
        action="store_true",
        help=(
            "If given, print configurations that would be submitted and their current "
            "calc status, but do not submit."
        ),
    )
    submit.add_argument(
        "--cancel",
        action="store_true",
        help=(
            'If given, calculations with status="started" or "submitted" are canceled. '
            "No jobs are submitted."
        ),
    )
    submit.add_argument(
        "-l",
        "--list",
        action="store_true",
        help=("If given, list enumeration, calctype, clex, and selection options."),
    )


def make_submitd_subparser(c):
    """Constructs the ``casm-calc submitd ...`` argument parser, and attaches the
    methods for running the subcommands.

    Parameters
    ----------
    c: argparse._SubParsersAction
        The output from ``parser.add_subparsers`` to which ``casm-calc submitd``
        arguments are added.

    Returns
    -------
    code: int
        A return code indicating success (0) or failure (non-zero).

    """
    import pathlib

    submitd = c.add_parser(
        "submitd",
        help=("Submit calculations in a directory, recursively"),
        description=(
            "Submit calculations in all calculation subdirectories of a specified top "
            "level directory. Calculation subdirectories are identified by the "
            "presence of a `status.json` file. Each calculation directory is expected "
            "to have a `submit.sh` script that can be submitted with `sbatch`. "
        ),
    )

    ### casm-calc submitd ....
    submitd.set_defaults(func=run_submitd)

    submitd.add_argument(
        "dir",
        type=pathlib.Path,
        help=("The top level directory."),
    )
    submitd.add_argument(
        "--name-parts",
        type=int,
        help=(
            "The number of path parts to include in the name. For "
            "example, if the calculation directory is `a/b/c/d` and `name_parts` is 2, "
            "the name will be `c/d`. If `name_parts` is greater than the number of "
            "parts in the path, the entire path is used. Default is 2."
        ),
        default=2,
    )
    submitd.add_argument(
        "--dry-run",
        action="store_true",
        help=(
            "If given, print configurations that would be submitted and their current "
            "calc status, but do not submit."
        ),
    )
    submitd.add_argument(
        "--cancel",
        action="store_true",
        help=(
            'If given, calculations with status="started" or "submitted" are canceled. '
            "No jobs are submitted."
        ),
    )
