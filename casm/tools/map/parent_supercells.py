import math
from typing import Optional

import numpy as np

import libcasm.configuration as casmconfig
import libcasm.xtal as xtal


class FromList:
    """A generator for parent supercells based on a list of transformation matrices.

    This generator yields transformation matrices for parent_T_list superstructures
    that are provided in the `parent_T_list`.
    """

    def __init__(self, parent_T_list: list[np.ndarray]):
        self.parent_T_list = parent_T_list

    def __call__(
        self,
        child_T: np.ndarray,
    ):
        for parent_T in self.parent_T_list:
            yield parent_T


class ByVolume:
    """A generator for parent supercells based on a range of volumes.

    This generator yields transformation matrices for superstructures with volume in
    the range `[min_parent_vol, max_parent_vol]`, inclusive, where the volumes are
    given as integer multiples of the parent prim volume.
    """

    def __init__(self, min_parent_vol: int, max_parent_vol: int):
        self.min_parent_vol = min_parent_vol
        self.max_parent_vol = max_parent_vol

    def __call__(
        self,
        child_T: np.ndarray,
    ):
        parent_superlattices = xtal.enumerate_superlattices(
            unit_lattice=self.parent_prim.xtal_prim.lattice(),
            point_group=self.parent_prim.crystal_point_group.elements,
            max_volume=self.max_parent_vol,
            min_volume=self.min_parent_vol,
        )
        for parent_superlattice in parent_superlattices:
            yield xtal.make_transformation_matrix_to_super(
                unit_lattice=self.parent_prim.xtal_prim.lattice(),
                superlattice=parent_superlattice,
            )


class ByAtomsPerUnitcell:
    """A generator for parent supercells based on the number of atoms per parent unit
    cell.

    This generator yields transformation matrices for parent superstructures that have
    a number of atoms per parent unit cell between `min_atoms_per_unitcell` and
    `max_atoms_per_unitcell`, inclusive.

    .. rubric:: Special methods

    - `__call__(child_T: np.ndarray)`: Given a transformation matrix for the child
      structure, yields transformation matrices for parent superstructures that have
      a number of atoms per parent unit cell between `min_atoms_per_unitcell` and
      `max_atoms_per_unitcell`, inclusive.

    """

    def __init__(
        self,
        child: xtal.Structure,
        parent_prim: casmconfig.Prim,
        min_atoms_per_unitcell: float,
        max_atoms_per_unitcell: float,
    ):
        """

        .. rubric:: Constructor


        Parameters
        ----------
        child: xtal.Structure
            The child structure.
        parent_prim: casmconfig.Prim
            The parent prim.
        min_atoms_per_unitcell: float
            The minimum number of atoms per parent unit cell the parent superstructure
            should have.
        max_atoms_per_unitcell: float
            The maximum number of atoms per parent unit cell the parent superstructure
            should have.

        """

        self.child = child
        self.parent_prim = parent_prim
        self.min_atoms_per_unitcell = min_atoms_per_unitcell
        self.max_atoms_per_unitcell = max_atoms_per_unitcell

        if (
            self.min_atoms_per_unitcell <= 0.0
            or self.max_atoms_per_unitcell <= 0.0
            or self.min_atoms_per_unitcell > self.max_atoms_per_unitcell
        ):
            raise ValueError(
                "Error in ParentVolumeSearchOptions: invalid range: "
                f"({self.min_atoms_per_unitcell}, {self.min_atoms_per_unitcell})'."
            )

    def __call__(
        self,
        child_T: np.ndarray,
    ):
        child_vol = int(round(np.linalg.det(child_T)))
        n_child_atoms = len(self.child.atom_type()) * child_vol

        min_parent_vol = int(math.floor(n_child_atoms / self.max_atoms_per_unitcell))
        max_parent_vol = int(math.ceil(n_child_atoms / self.min_atoms_per_unitcell))

        parent_superlattices = xtal.enumerate_superlattices(
            unit_lattice=self.parent_prim.xtal_prim.lattice(),
            point_group=self.parent_prim.crystal_point_group.elements,
            max_volume=max_parent_vol,
            min_volume=min_parent_vol,
        )
        for parent_superlattice in parent_superlattices:
            yield xtal.make_transformation_matrix_to_super(
                unit_lattice=self.parent_prim.xtal_prim.lattice(),
                superlattice=parent_superlattice,
            )


class ByPointDefectCount:
    """A generator for parent supercells based on the number of vacancies and
    interstitials expected to be found.

    Calculates the parent supercell volume as:

    .. math::

        V_{parent} = \\frac{N_{child} - N_{interstitial} +
        N_{vacancy}}{N_{expected\\_per}}

    where :math:`N_{child}` is the number of atoms in the child supercell,
    :math:`N_{interstitial}` is the number of interstitials expected,
    :math:`N_{vacancy}` is the number of vacancies expected, and
    :math:`N_{expected\\_per}` is the number of atoms expected in the parent unit cell
    if there are no point defects. By default, :math:`N_{expected\\_per}` is the
    number of sites in the `parent_prim` which do not have a vacancy as the first
    occupant.

    .. rubric:: Special methods

    - `__call__(child_T: np.ndarray)`: Given a transformation matrix for the child
      structure, yields transformation matrices for parent superstructures that have
      the requested number of vacancies and interstitials, assuming that atoms fill
      sublattices as expected.

    """

    def __init__(
        self,
        child: xtal.Structure,
        parent_prim: casmconfig.Prim,
        expected_n_vacancy: Optional[int] = None,
        expected_n_interstitial: Optional[int] = None,
        n_expected_atoms_per_parent_unitcell: Optional[int] = None,
    ):
        """

        .. rubric:: Constructor


        Parameters
        ----------
        child: xtal.Structure
            The child structure.
        parent_prim: casmconfig.Prim
            The parent prim.
        expected_n_vacancy: Optional[int] = None
            The number of vacancies expected after mapping, if the `method` is set to
            "point-defect-count". If not provided, no vacancies are expected.
        expected_n_interstitial: Optional[int] = None
            The number of interstitials expected after mapping, if the `method` is set
            to "point-defect-count". If not provided, no interstitials are expected.
        n_expected_atoms_per_parent_unitcell: Optional[int] = None
            The number of atoms expected in the parent unit cell if there are no point
            defects, used if `method` is set to "point-defect-count". By default, the
            number of sites in the `parent_prim` which do not have a vacancy as the
            first occupant is used. If this is provided, it overrides the default value.


        """
        self.child = child
        self.parent_prim = parent_prim

        self.expected_n_vacancy: Optional[int] = expected_n_vacancy
        """Optional[int]: The number of vacancies expected after mapping, if the 
        `method` is set to "point-defect-count"."""

        self.expected_n_interstitial: Optional[int] = expected_n_interstitial
        """Optional[int]: The number of interstitials expected after mapping, if the 
        `method` is set to "point-defect-count"."""

        # Make default n_expected_atoms_per_parent_unitcell if necessary:
        if n_expected_atoms_per_parent_unitcell is None:
            vacancy_names = ["Va", "va", "VA"]
            _expected_per = 0
            occ_dof = parent_prim.xtal_prim.occ_dof()
            for site_dof in occ_dof:
                if len(site_dof):
                    if site_dof[0] not in vacancy_names:
                        _expected_per += 1
            n_expected_atoms_per_parent_unitcell = _expected_per

        self.n_expected_atoms_per_parent_unitcell: Optional[int] = (
            n_expected_atoms_per_parent_unitcell
        )
        """Optional[int]: The number of atoms expected in the parent unit cell if there
        are no point defects, used if `method` is set to "point-defect-count".

        By default, the number of sites in the `parent_prim` which do not have a 
        vacancy as the first occupant is used. If this is provided, it overrides the
        default value."""

    def __call__(
        self,
        child_T: np.ndarray,
    ):
        child_vol = int(round(np.linalg.det(child_T)))
        n_child_atoms = len(self.child.atom_type()) * child_vol

        n_atoms = n_child_atoms
        n_atoms -= self.expected_n_interstitial
        n_atoms += self.expected_n_vacancy

        parent_vol = n_atoms / self.n_expected_atoms_per_parent_unitcell

        # check if parent_vol is approximately integer:
        if not np.isclose(parent_vol, round(parent_vol), atol=1e-5):
            raise ValueError(
                "Error in ByPointDefectCount: "
                "the expected number of vacancies and interstitials does not result in "
                "an integer parent supercell volume."
            )

        parent_vol = int(round(parent_vol))

        parent_superlattices = xtal.enumerate_superlattices(
            unit_lattice=self.parent_prim.xtal_prim.lattice(),
            point_group=self.parent_prim.crystal_point_group.elements,
            max_volume=parent_vol,
            min_volume=parent_vol,
        )
        for parent_superlattice in parent_superlattices:
            yield xtal.make_transformation_matrix_to_super(
                unit_lattice=self.parent_prim.xtal_prim.lattice(),
                superlattice=parent_superlattice,
            )
