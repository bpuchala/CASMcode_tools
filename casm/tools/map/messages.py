import pathlib
import sys


def results_dir_exists_error(results_dir: pathlib.Path) -> None:
    """Print an error message if the results directory already exists."""

    error = f"""
################################################################################
# Error: Results directory already exists                                      #
#                                                                              #
# A directory already exists at the specified path.                            #
#                                                                              #
# To merge new results, use --merge. Otherwise, delete the existing directory  #
# or specify a new one.                                                        #

--results-dir={results_dir}

# Stopping...                                                                  #
################################################################################
"""
    print(error)
    sys.exit(1)


def atom_types_mismatch_error(parent_atom_types, child_atom_types) -> None:
    """Print an error message if the parent and child atom types do not match."""
    error = f"""
    ################################################################################
    # Error: Child cannot map to parent due to atom types mismatch                 #
    #                                                                              #

    - Parent atom types: {parent_atom_types}
    - Child atom types: {child_atom_types}

    # Stopping...                                                                  #
    ################################################################################
    """
    print(error)
    sys.exit(1)


def atom_fraction_mismatch_error(parent_atom_frac, child_atom_frac) -> None:
    """Print an error message if the parent and child atom types do not match."""
    error = f"""
    ################################################################################
    # Error: Parent and child structures have different atom fractions             #
    #                                                                              #

    - Parent atom fraction: {parent_atom_frac}
    - Child atom fraction: {child_atom_frac}

    # Stopping...                                                                  #
    ################################################################################
    """
    print(error)
    sys.exit(1)


def invalid_forced_on_values_error(
    parent_site_index: int,
    child_atom_index: int,
    child_type: str,
    allowed_types: list[str],
) -> None:
    """Print an error message if the `--forced-on` option is used with invalid
    values."""

    error = f"""
################################################################################
# Error: Invalid --forced-on values                                            #
#                                                                              #
# The `--forced-on` option requires that the child atom type is allowed on the #
# parent site.                                                                 #

child_atom_index={child_atom_index}
child_type={child_type}
parent_site_index={parent_site_index}
allowed_types={allowed_types}

# Stopping...                                                                  #
################################################################################
"""
    print(error)
    sys.exit(1)


def invalid_fix_parent_error() -> None:
    """Print an error message if the `--fix-parent` option is used with a parent and
    child structure that have different numbers of atoms."""

    error = """
################################################################################
# Error: --fix-parent requires parent and child w/ same number of atoms.       #
#                                                                              #
# The `--fix-parent` option is used to map to the parent structure as          #
# provided, without searching over parent superstructures and lattice          #
# reorientations. It is only allowed if the number of atoms in the parent      #
# structure is the same as the number of atoms in the child structure.         #
#                                                                              #
# Stopping...                                                                  #
################################################################################
"""
    print(error)
    sys.exit(1)


def different_parent_error(results_dir: pathlib.Path) -> None:
    """When merging, print an error message if the parent has changed."""

    error = """
################################################################################
# Error: parent structure has changed                                          #
#                                                                              #
# When using the --merge option, the parent structure must remain the same.    #
#                                                                              #
# Stopping...                                                                  #
################################################################################
"""
    print(error)
    sys.exit(1)


def different_child_error(results_dir: pathlib.Path) -> None:
    """When merging, print an error message if the child has changed."""

    error = """
################################################################################
# Error: child structure has changed                                           #
#                                                                              #
# When using the --merge option, the child structure must remain the same.     #
#                                                                              #
# Stopping...                                                                  #
################################################################################
"""
    print(error)
    sys.exit(1)


def different_lattice_mapping_cost_method_error() -> None:
    """When merging, print an error message if the lattice mapping cost method has
    changed."""

    error = """
################################################################################
# Error: lattice mapping cost method has changed                               #
#                                                                              #
# When using the --merge option, the lattice mapping cost method must remain   #
# the same.                                                                    #
#                                                                              #
# Stopping...                                                                  #
################################################################################
"""
    print(error)
    sys.exit(1)


def different_atom_mapping_cost_method_error() -> None:
    """When merging, print an error message if the atom mapping cost method has
    changed."""

    error = """
################################################################################
# Error: atom mapping cost method has changed                                  #
#                                                                              #
# When using the --merge option, the atom mapping cost method must remain      #
# the same.                                                                    #
#                                                                              #
# Stopping...                                                                  #
################################################################################
"""
    print(error)
    sys.exit(1)


def different_lattice_cost_weight_error() -> None:
    """When merging, print an error message if the lattice cost weight has changed."""

    error = """
################################################################################
# Error: lattice cost weight has changed                                       #
#                                                                              #
# When using the --merge option, the lattice cost weight must remain the same. #
#                                                                              #
# Stopping...                                                                  #
################################################################################
"""
    print(error)
    sys.exit(1)


def primitive_parent_notice() -> None:
    """Write a notice to the console that the parent is not primitive."""

    notice = """
################################################################################
# Notice: parent is not primitive                                              #
# Writing primitive parent: parent.primitive.json                              #
#                                                                              #
# The parent is not primitive, and the search will continue with the           #
# non-primitive parent structure. If you want to use the primitive parent,     #
# please use the file `parent.primitive.json` instead.                         #
################################################################################
"""
    print(notice)
    sys.stdout.flush()


def primitive_child_notice() -> None:
    """Write a notice to the console that the child structure is not primitive."""

    notice = """
################################################################################
# Notice: child is not primitive                                               #
# Writing primitive child structure: child.primitive.json                      #
#                                                                              #
# The child structure is not primitive, and the search will continue with the  #
# non-primitive child structure. If you want to use the primitive child,       #
# please use the file `child.primitive.json` instead.                           #
################################################################################
"""
    print(notice)
    sys.stdout.flush()


def invalid_min_n_atoms_error(min_n_atoms: int):
    """Print an error message for invalid min_n_atoms."""

    error = f"""
################################################################################
# Error: Invalid min_n_atoms                                                   #
#                                                                              #
# The value of min_n_atoms must be at least 1.                                 #
#                                                                              #

--min-n-atoms={min_n_atoms}

# Stopping...                                                                  #
################################################################################
"""
    print(error)
    sys.exit(1)


def invalid_max_n_atoms_error(
    min_n_atoms: int,
    max_n_atoms: int,
    computed_msg: str,
):
    """Print an error message for invalid max_n_atoms."""

    error = f"""
################################################################################
# Error: Invalid max_n_atoms                                                   #
#                                                                              #
# The value of max_n_atoms must be greater than or equalt to min_n_atoms.      #
# equal to the minimum.                                                        #

--min-n-atoms={min_n_atoms}
--max-n-atoms={max_n_atoms} {computed_msg}

# Stopping...                                                                  #
################################################################################
"""
    print(error)
    sys.exit(1)


def invalid_lattice_mapping_cost_method_error(method: str):
    """Print an error message for invalid lattice mapping cost method."""

    error = f"""
################################################################################
# Error: Invalid lattice mapping cost method                                   #
#                                                                              #
# The lattice mapping cost method must be one of:                              #
# - 'isotropic_strain_cost'                                                    #
# - 'symmetry_breaking_strain_cost'                                            #
#                                                                              #

--lattice-cost-method={method}

# Stopping...                                                                  #
################################################################################
"""
    print(error)
    sys.exit(1)


def invalid_atom_mapping_cost_method_error(method: str):
    """Print an error message for invalid atom mapping cost method."""

    error = f"""
################################################################################
# Error: Invalid atom mapping cost method                                      #
#                                                                              #
# The atom mapping cost method must be one of:                                 #
# - 'isotropic_disp_cost'                                                      #
# - 'symmetry_breaking_disp_cost'                                              #
#                                                                              #

--atom-cost-method={method}

# Stopping...                                                                  #
################################################################################
"""
    print(error)
    sys.exit(1)


def invalid_deduplication_interpolation_factors_error(dedup_factors):
    """Print an error message for invalid deduplication interpolation factors."""

    error = f"""
################################################################################
# Error: Invalid deduplication interpolation factors                           #
#                                                                              #
# The deduplication interpolation factors must be a list of floats.            #
#                                                                              #

--dedup-interp-factors={dedup_factors}

# Stopping...                                                                  #
################################################################################
"""
    print(error)
    sys.exit(1)
