import math
import typing

import numpy as np

import libcasm.xtal as xtal


# Take a CASM structure and target kpra (k-points per reciprocal atom) and
# generate a k-point mesh that meets or exceeds the target kpra.
def generate_kpoints_mesh_kpra(
    structure: xtal.Structure,
    kpra: float,
) -> list[int]:
    """Generate k-point mesh for a given structure and target k-points per
    reciprocal atom (kpra).

    Parameters
    ----------
    structure : xtal.Structure
        The structure to generate the k-point mesh for.
    kpra : float
        The target k-points per reciprocal atom.

    Returns
    -------
    list[int]
        The k-point mesh as a list of three integers.
    """
    # Get the reciprocal lattice vectors
    rec_lattice = structure.lattice().reciprocal()

    # Calculate the lengths of the reciprocal lattice vectors
    rec_lengths = rec_lattice.lengths_and_angles()[:3]

    # Calculate the number of atoms in the unit cell
    num_atoms = len(structure.atom_type())

    # Calculate the k-point mesh
    # The total number of k-points should be approximately kpra * num_atoms
    # We distribute k-points proportionally to reciprocal lattice vector lengths
    # k1 * k2 * k3 ~ kpra * num_atoms
    # k_i ~ (kpra * num_atoms)^(1/3) * (rec_length_i / geometric_mean(rec_lengths))

    target_total_kpoints = kpra * num_atoms
    geometric_mean = (rec_lengths[0] * rec_lengths[1] * rec_lengths[2]) ** (1.0 / 3.0)

    kpoint_mesh = []
    for length in rec_lengths:
        k_i = (target_total_kpoints ** (1.0 / 3.0)) * (length / geometric_mean)
        kpoint_mesh.append(max(1, int(math.ceil(k_i))))

    return kpoint_mesh


def generate_kpoints_mesh_rk(
    structure: xtal.Structure,
    rk: float,
) -> list[int]:
    """Generate k-point mesh for a given structure and target Rk (real-space
    distance between k-points).

    Parameters
    ----------
    structure : xtal.Structure
        The structure to generate the k-point mesh for.
    rk : float
        The target Rk (real-space distance between k-points).

    Returns
    -------
    list[int]
        The k-point mesh as a list of three integers.
    """
    # Get the reciprocal lattice vectors
    rec_lattice = structure.lattice().reciprocal()

    # Calculate the lengths of the reciprocal lattice vectors
    rec_lengths = np.array(rec_lattice.lengths_and_angles()[:3])

    # Calculate the k-point mesh
    # kspacing_i = 2*pi / Rk
    # n_i = rec_length_i * 2*pi / kspacing_i = rec_length_i * Rk

    kpoint_mesh = []
    for length in rec_lengths:
        k_i = length * rk / (2.0 * math.pi)
        kpoint_mesh.append(max(1, int(math.ceil(k_i))))

    return kpoint_mesh


def generate_kpoints_mesh_kspacing(
    structure: xtal.Structure,
    kspacing: float,
) -> list[int]:
    """Generate k-point mesh for a given structure and target k-spacing.

    Parameters
    ----------
    structure : xtal.Structure
        The structure to generate the k-point mesh for.
    kspacing : float
        The target k-spacing.

    Returns
    -------
    list[int]
        The k-point mesh as a list of three integers.
    """
    # Get the reciprocal lattice vectors
    rec_lattice = structure.lattice().reciprocal()

    # Calculate the lengths of the reciprocal lattice vectors
    rec_lengths = np.array(rec_lattice.lengths_and_angles()[:3])

    # Calculate the k-point mesh
    # n_i = rec_length_i * kspacing

    kpoint_mesh = []
    for length in rec_lengths:
        k_i = length / kspacing
        kpoint_mesh.append(max(1, int(math.ceil(k_i))))

    return kpoint_mesh


def generate_kpoints_mesh(
    structure: xtal.Structure,
    kpra: typing.Optional[float] = None,
    rk: typing.Optional[float] = None,
    kspacing: typing.Optional[float] = None,
) -> list[int]:
    """Generate k-point mesh for a given structure and target kpra, rk, or kspacing.

    Notes
    -----

    Exactly one of kpra, rk, or kspacing must be provided.

    Parameters
    ----------
    structure : xtal.Structure
        The structure to generate the k-point mesh for.
    kpra : Optional[float] = None
        The k-points per reciprocal atom.
    rk : Optional[float] = None
        The k-points length, :math:`R_k`.
    kspacing : Optional[float] = None
        The KSPACING value.

    Returns
    -------
    list[int]
        The k-point mesh as a list of three integers.
    """
    # If >1 of kpra, rk, kspacing are provided, raise an error:
    num_provided = sum(x is not None for x in [kpra, rk, kspacing])
    if num_provided != 1:
        raise ValueError("Exactly one of kpra, rk, or kspacing must be provided.")

    if kpra is not None:
        return generate_kpoints_mesh_kpra(structure, kpra)
    elif rk is not None:
        return generate_kpoints_mesh_rk(structure, rk)
    elif kspacing is not None:
        return generate_kpoints_mesh_kspacing(structure, kspacing)
    else:
        raise ValueError("One of kpra, rk, or kspacing must be provided.")


def kpoints_info(
    structure: xtal.Structure,
    kpoints_mesh: list[int],
) -> dict[str, float]:
    """Calculate k-points information from a given k-point mesh.

    Parameters
    ----------
    structure : xtal.Structure
        The structure to calculate k-points information for.
    kpoints_mesh : list[int]
        The k-point mesh as a list of three integers.

    Returns
    -------
    info: dict
        A dictionary containing:

        .. code-block:: python

            {
                "total_kpoints": int,
                "kpra": float,
                "kspacing": list[float],
                "rk": list[float],
            }
    """
    rec_lattice = structure.lattice().reciprocal()
    rec_lengths = np.array(rec_lattice.lengths_and_angles()[:3])
    num_atoms = len(structure.atom_type())

    total_kpoints = kpoints_mesh[0] * kpoints_mesh[1] * kpoints_mesh[2]
    kpra = total_kpoints / num_atoms

    # kspacing: reciprocal space distance between k-points in each direction
    kspacing = rec_lengths / kpoints_mesh

    # rk: real-space distance between k-points in each direction
    rk = (2.0 * math.pi) / kspacing

    return {
        "total_kpoints": total_kpoints,
        "kpra": kpra,
        "kspacing": kspacing.tolist(),
        "rk": rk.tolist(),
    }
