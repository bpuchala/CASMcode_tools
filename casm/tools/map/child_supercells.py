import numpy as np

import libcasm.xtal as xtal

from . import methods as mthds


class FromList:
    """A generator for child supercells based on a list of transformation matrices.

    This generator yields transformation matrices for child superstructures that are
    provided in the `child_T_list`.
    """

    def __init__(self, child_T_list: list[np.ndarray]):
        self.child_T_list = child_T_list

    def __call__(self, child: xtal.Structure):
        for child_T in self.child_T_list:
            yield child_T


class ByNumberOfAtoms:
    """A generator for child supercells based on the number of atoms in the child
    structure.

    This generator yields transformation matrices for child superstructures that have
    a number of atoms between `min_n_atoms` and `max_n_atoms`, inclusive.
    """

    def __init__(self, min_n_atoms: int, max_n_atoms: int):
        self.min_n_atoms = min_n_atoms
        self.max_n_atoms = max_n_atoms

    def __call__(self, child: xtal.Structure):
        # Parameters
        child_crystal_point_group = xtal.make_structure_crystal_point_group(child)
        child_n_atoms = len(child.atom_type())

        self._child_T_list = []

        child_superlattices = xtal.enumerate_superlattices(
            unit_lattice=child.lattice(),
            point_group=child_crystal_point_group,
            max_volume=mthds.floordiv(self.max_n_atoms, child_n_atoms),
            min_volume=mthds.ceildiv(self.min_n_atoms, child_n_atoms),
        )
        for child_superlattice in child_superlattices:
            yield xtal.make_transformation_matrix_to_super(
                unit_lattice=child.lattice(),
                superlattice=child_superlattice,
            )
