from typing import Optional

import numpy as np


class TotalMappingOptions:
    """Total structure mappping options."""

    def __init__(
        self,
        min_cost: float = 0.0,
        max_cost: float = 0.3,
        k_best: int = 1,
        lattice_cost_weight: float = 0.5,
        cost_tol: float = 1e-5,
        deduplication_interpolation_factors: Optional[list[float]] = None,
    ):
        """

        .. rubric:: Constructor

        Parameters
        ----------
        tmin_cost : float = 0.0
            The minimum total cost mapping to include in search results.
        max_cost : float = 0.3
            The maximum total cost mapping to include in search results.
        k_best : int = 1
            Keep the `k_best` mappings with the lowest total cost that also
            satisfy the min/max cost criteria. Approximate ties with the
            current `k_best`-ranked result are also kept.
        lattice_cost_weight : float = 0.5
            The weight of the lattice cost in the total structure mapping cost.
        cost_tol : float = 1e-5
            Tolerance for checking if mapping costs are approximately equal.
        deduplication_interpolation_factors : Optional[list[float]] = None
            Interpolation factors to use for deduplication. If None, the default value
            ``[0.5, 1.0]`` is used.

        """

        self.min_cost: float = min_cost
        self.max_cost: float = max_cost
        self.k_best: float = k_best
        self.lattice_cost_weight: float = lattice_cost_weight
        self.cost_tol: float = cost_tol

        if deduplication_interpolation_factors is None:
            deduplication_interpolation_factors = [0.5, 1.0]

        self.deduplication_interpolation_factors: Optional[list[float]] = (
            deduplication_interpolation_factors
        )


class LatticeMappingStepOptions:
    """Options for the lattice mapping step in structure mapping search."""

    def __init__(
        self,
        cost_method: str = "symmetry_breaking_strain_cost",
        min_cost: float = 0.0,
        max_cost: float = 1e20,
        k_best: int = 10,
        reorientation_range: int = 1,
        fix_parent: bool = False,
        cost_tol: float = 1e-5,
    ):
        """
        .. rubric:: Constructor

        Parameters
        ----------
        cost_method : str = 'symmetry_breaking_strain_cost'
            Selects the method used to calculate lattice mapping costs. Used when
            `map_lattices_with_reorientation` is True. One of
            "isotropic_strain_cost" or "symmetry_breaking_strain_cost".
        min_cost : float = 0.0
            Keep lattice mappings with cost >= min_cost. Used when
            `map_lattices_with_reorientation` is True.
        max_cost : float = 1e20
            Keep results with cost <= max_cost. Used when
            `map_lattices_with_reorientation` is True.
        k_best : int = 10
            If not None, then only keep the k-best results (i.e. k lattice mappings
            with minimum cost) satisfying the min_cost and max_cost constraints.
            If there are approximate ties, those will also be kept. Used when
            `map_lattices_with_reorientation` is True.
        reorientation_range : int = 1
            The absolute value of the maximum element in the lattice mapping
            reorientation matrix, :math:`N`. This determines how many equivalent
            lattice vector reorientations are checked. Increasing the value results in
            more checks. The value 1 is generally expected to be sufficient because
            reduced cell lattices are compared internally.
        fix_parent : bool = False
            If True, map to the parent structure as provided and skip searching over
            parent superstructures and lattice reorientations. The deformation
            gradient is still calculated and atom mapping is still performed. Only
            allowed if the number of atoms in the parent structure is the same as
            the number of atoms in the child structure.
        cost_tol : float = 1e-5
            Tolerance for checking if mapping costs are approximately equal.

        """
        self.cost_method = cost_method
        self.min_cost = min_cost
        self.max_cost = max_cost
        self.k_best = k_best
        self.reorientation_range = reorientation_range
        self.fix_parent = fix_parent
        self.cost_tol = cost_tol


class AtomMappingStepOptions:

    def __init__(
        self,
        cost_method: str = "symmetry_breaking_disp_cost",
        remove_mean_displacement: bool = True,
        forced_on: Optional[dict[int, int]] = None,
        forced_off: Optional[list[tuple[int, int]]] = None,
    ):
        """
        .. rubric:: Constructor

        Parameters
        ----------
        cost_method : str = 'symmetry_breaking_disp_cost'
            Selects the method used to calculate atom mapping costs. One of
            "isotropic_disp_cost" or "symmetry_breaking_disp_cost".
        remove_mean_displacement : bool = True
            If True, the mean displacement is removed from atom mappings. Otherwise,
            the mean displacement is not removed and at least one atom is mapped without
            any displacement.
        forced_on : Optional[dict[int, int]] = None
            A map of assignments `parent_atom_index: child_atom_index` that are forced
            on. Indices begin at 0. Requires that `fix_parent` is True for the lattice
            mapping step. If `fix_parent` is False, this option is ignored.
        forced_off : Optional[list[tuple[int, int]]] = None
            A list of tuples of assignments `(parent_atom_index, child_atom_index) that
            are forced off. Indices begin at 0. Requires that `fix_parent` is True for
            the lattice mapping step. If `fix_parent` is False, this option is ignored.
        """
        self.cost_method = cost_method
        self.remove_mean_displacement = remove_mean_displacement
        self.forced_on = forced_on
        self.forced_off = forced_off


class SupercellsOptions:
    """Options for supercell generation"""

    def __init__(
        self,
        method: str,
        kwargs: dict,
    ):
        """

        .. rubric:: Constructor

        Parameters
        ----------
        method: str
            The name of a supercell generation method class, i.e. "FromList",
            "ByNumberOfAtoms", etc.
        kwargs: dict
            Arguments for the supercell generation method constructor.
        """
        self.method: str = method
        self.kwargs: dict = kwargs

    def to_dict(self):
        return {
            "method": self.method,
            "kwargs": self.kwargs,
        }

    @staticmethod
    def from_dict(data: dict):
        return SupercellsOptions(
            method=data["method"],
            kwargs=data["kwargs"],
        )


class StructureMappingSearchOptions_v2:
    def __init__(
        self,
        total: Optional[TotalMappingOptions] = None,
        child_supercells: Optional[SupercellsOptions] = None,
        parent_supercells: Optional[SupercellsOptions] = None,
        lattice_step: Optional[LatticeMappingStepOptions] = None,
        atom_step: Optional[AtomMappingStepOptions] = None,
    ):
        if total is None:
            total = TotalMappingOptions()
        if lattice_step is None:
            lattice_step = LatticeMappingStepOptions()
        if atom_step is None:
            atom_step = AtomMappingStepOptions()

        self.total = total
        self.child_supercells = child_supercells
        self.parent_supercells = parent_supercells
        self.lattice_step = lattice_step
        self.atom_step = atom_step


class StructureMappingSearchOptions:
    """Options controlling the structure mapping search."""

    def __init__(
        self,
        max_n_atoms: Optional[int] = None,
        min_n_atoms: int = 1,
        child_transformation_matrix_to_super_list: Optional[list[np.ndarray]] = None,
        parent_transformation_matrix_to_super_list: Optional[list[np.ndarray]] = None,
        total_min_cost: float = 0.0,
        total_max_cost: float = 0.3,
        total_k_best: int = 1,
        no_remove_mean_displacement: bool = False,
        fix_parent: bool = False,
        lattice_mapping_min_cost: Optional[float] = 0.0,
        lattice_mapping_max_cost: Optional[float] = 1e20,
        lattice_mapping_k_best: Optional[int] = 10,
        lattice_mapping_reorientation_range: Optional[int] = 1,
        lattice_mapping_cost_method: str = "symmetry_breaking_strain_cost",
        atom_mapping_cost_method: str = "symmetry_breaking_disp_cost",
        forced_on: Optional[dict[int, int]] = None,
        forced_off: Optional[list[tuple[int, int]]] = None,
        lattice_cost_weight: float = 0.5,
        cost_tol: Optional[float] = 1e-5,
        deduplication_interpolation_factors: Optional[list[float]] = None,
    ):
        """

        .. rubric:: Constructor

        Parameters
        ----------
        max_n_atoms : Optional[int] = None
            The maximum number of atoms in the superstructures that should be included
            in the search. If None, the least common multiple of the number of atoms
            in the child and parent structures.
        min_n_atoms : int = 1
            The minimum number of atoms in the superstructures that should be included
            in the search.
        child_transformation_matrix_to_super_list : Optional[list[np.ndarray]] = None
            If provided, overrides the `min_n_atoms` and `max_n_atoms` options to
            directly specify the transformation matrices to use for the child
            superstructures.
        parent_transformation_matrix_to_super_list : Optional[list[np.ndarray]] = None
            If provided, only use the specified transformation matrices to create
            parent superstructures. If None, the parent superstructures are
            enumerated.
        total_min_cost : float = 0.0
            The minimum total cost mapping to include in search results.
        total_max_cost : float = 0.3
            The maximum total cost mapping to include in search results.
        total_k_best : int = 1
            Keep the `k_best` mappings with the lowest total cost that also
            satisfy the min/max cost criteria. Approximate ties with the
            current `k_best`-ranked result are also kept.
        no_remove_mean_displacement : bool = False
            If True, do not remove the mean displacement from the atom mapping.
        fix_parent : bool = False
            If True, map to the parent structure as provided and skip searching over
            parent superstructures and lattice reorientations. The deformation
            gradient is still calculated and atom mapping is still performed. Only
            allowed if the number of atoms in the parent structure is the same as
            the number of atoms in the child structure.
        lattice_mapping_min_cost : float = 0.0
            Keep lattice mappings with cost >= min_cost. Used when
            `map_lattices_with_reorientation` is True.
        lattice_mapping_max_cost : float = 1e20
            Keep results with cost <= max_cost. Used when
            `map_lattices_with_reorientation` is True.
        lattice_mapping_k_best : int = 10
            If not None, then only keep the k-best results (i.e. k lattice mappings
            with minimum cost) satisfying the min_cost and max_cost constraints.
            If there are approximate ties, those will also be kept. Used when
            `map_lattices_with_reorientation` is True.
        lattice_mapping_reorientation_range : int = 1
            The absolute value of the maximum element in the lattice mapping
            reorientation matrix, :math:`N`. This determines how many equivalent
            lattice vector reorientations are checked. Increasing the value results in
            more checks. The value 1 is generally expected to be sufficient because
            reduced cell lattices are compared internally.
        lattice_mapping_cost_method : str = 'symmetry_breaking_strain_cost'
            Selects the method used to calculate lattice mapping costs. Used when
            `map_lattices_with_reorientation` is True. One of
            "isotropic_strain_cost" or "symmetry_breaking_strain_cost".
        atom_mapping_cost_method : str = 'symmetry_breaking_disp_cost'
            Selects the method used to calculate atom mapping costs. One of
            "isotropic_disp_cost" or "symmetry_breaking_disp_cost".
        forced_on : Optional[dict[int, int]] = None
            A map of assignments `parent_atom_index: child_atom_index` that are forced
            on. Indices begin at 0. Requires that `fix_parent` is True.
        forced_off : Optional[list[tuple[int, int]]] = None
            A list of tuples of assignments `(parent_atom_index, child_atom_index) that
            are forced off. Indices begin at 0. Requires that `fix_parent` is True.
        lattice_cost_weight : float = 0.5
            The weight of the lattice cost in the total structure mapping cost.
        cost_tol : float = 1e-5
            Tolerance for checking if mapping costs are approximately equal.
        deduplication_interpolation_factors : Optional[list[float]] = None
            Interpolation factors to use for deduplication. If None, the default value
            ``[0.5, 1.0]`` is used.

        """
        self.min_n_atoms = min_n_atoms
        self.max_n_atoms = max_n_atoms
        self.child_transformation_matrix_to_super_list = (
            child_transformation_matrix_to_super_list
        )
        self.parent_transformation_matrix_to_super_list = (
            parent_transformation_matrix_to_super_list
        )

        self.total_min_cost = total_min_cost
        self.total_max_cost = total_max_cost
        self.total_k_best = total_k_best
        self.no_remove_mean_displacement = no_remove_mean_displacement
        self.fix_parent = fix_parent
        self.lattice_mapping_min_cost = lattice_mapping_min_cost
        self.lattice_mapping_max_cost = lattice_mapping_max_cost
        self.lattice_mapping_k_best = lattice_mapping_k_best
        self.lattice_mapping_reorientation_range = lattice_mapping_reorientation_range
        self.lattice_mapping_cost_method = lattice_mapping_cost_method
        self.atom_mapping_cost_method = atom_mapping_cost_method
        self.forced_on = forced_on
        self.forced_off = forced_off
        self.lattice_cost_weight = lattice_cost_weight
        self.cost_tol = cost_tol

        # Deduplication options
        if deduplication_interpolation_factors is None:
            deduplication_interpolation_factors = [0.5, 1.0]
        self.deduplication_interpolation_factors = deduplication_interpolation_factors

    def to_dict(self):
        return {
            "min_n_atoms": self.min_n_atoms,
            "max_n_atoms": self.max_n_atoms,
            "child_transformation_matrix_to_super_list": (
                [x.tolist() for x in self.child_transformation_matrix_to_super_list]
                if self.child_transformation_matrix_to_super_list is not None
                else None
            ),
            "parent_transformation_matrix_to_super_list": (
                [x.tolist() for x in self.parent_transformation_matrix_to_super_list]
                if self.parent_transformation_matrix_to_super_list is not None
                else None
            ),
            "total_min_cost": self.total_min_cost,
            "total_max_cost": self.total_max_cost,
            "total_k_best": self.total_k_best,
            "no_remove_mean_displacement": self.no_remove_mean_displacement,
            "fix_parent": self.fix_parent,
            "lattice_mapping_min_cost": self.lattice_mapping_min_cost,
            "lattice_mapping_max_cost": self.lattice_mapping_max_cost,
            "lattice_mapping_k_best": self.lattice_mapping_k_best,
            "lattice_mapping_reorientation_range": self.lattice_mapping_reorientation_range,  # noqa: E501
            "lattice_mapping_cost_method": self.lattice_mapping_cost_method,
            "atom_mapping_cost_method": self.atom_mapping_cost_method,
            "forced_on": (
                [[key, value] for key, value in self.forced_on.items()]
                if self.forced_on is not None
                else None
            ),
            "forced_off": self.forced_off,
            "lattice_cost_weight": self.lattice_cost_weight,
            "cost_tol": self.cost_tol,
            "deduplication_interpolation_factors": self.deduplication_interpolation_factors,  # noqa: E501
        }

    @staticmethod
    def from_dict(data: dict):
        return StructureMappingSearchOptions(
            max_n_atoms=data["max_n_atoms"],
            min_n_atoms=data["min_n_atoms"],
            child_transformation_matrix_to_super_list=(
                [np.array(x) for x in data["child_transformation_matrix_to_super_list"]]
                if data["child_transformation_matrix_to_super_list"] is not None
                else None
            ),
            parent_transformation_matrix_to_super_list=(
                [
                    np.array(x)
                    for x in data["parent_transformation_matrix_to_super_list"]
                ]
                if data["parent_transformation_matrix_to_super_list"] is not None
                else None
            ),
            total_min_cost=data["total_min_cost"],
            total_max_cost=data["total_max_cost"],
            total_k_best=data["total_k_best"],
            no_remove_mean_displacement=data["no_remove_mean_displacement"],
            fix_parent=data["fix_parent"],
            lattice_mapping_min_cost=data["lattice_mapping_min_cost"],
            lattice_mapping_max_cost=data["lattice_mapping_max_cost"],
            lattice_mapping_k_best=data["lattice_mapping_k_best"],
            lattice_mapping_reorientation_range=data[
                "lattice_mapping_reorientation_range"
            ],
            lattice_mapping_cost_method=data["lattice_mapping_cost_method"],
            atom_mapping_cost_method=data["atom_mapping_cost_method"],
            forced_on=(
                {x[0]: x[1] for x in data["forced_on"]}
                if data["forced_on"] is not None
                else None
            ),
            forced_off=(
                [tuple(x) for x in data["forced_off"]]
                if data["forced_off"] is not None
                else None
            ),
            lattice_cost_weight=data["lattice_cost_weight"],
            cost_tol=data["cost_tol"],
            deduplication_interpolation_factors=data[
                "deduplication_interpolation_factors"
            ],
        )
