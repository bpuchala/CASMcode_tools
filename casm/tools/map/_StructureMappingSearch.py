import math
import pathlib
import sys
import uuid
from typing import Optional

import numpy as np
from tabulate import tabulate

import libcasm.configuration as casmconfig
import libcasm.mapping.info as mapinfo
import libcasm.mapping.mapsearch as mapsearch
import libcasm.mapping.methods as mapmethods
import libcasm.xtal as xtal
from casm.tools.shared.json_io import (
    read_optional,
    read_required,
    safe_dump,
)

# from .methods import (
#     _get_max_n_atoms_for_parent_prim,
#     ceildiv,
#     chain_is_in_orbit,
#     floordiv,
#     make_chain_orbit,
#     make_child_supercell_info,
#     make_child_transformation_matrix_to_super,
#     make_parent_supercell_info,
#     make_primitive_chain,
#     make_primitive_chain_orbit,
#     parent_supercell_size,
# )
from . import messages as msgs
from . import methods as mthds
from ._StructureMappingSearchOptions import StructureMappingSearchOptions


class SearchResult:
    def __init__(self):
        self.parent_structure: Optional[xtal.Structure] = None
        """Optional[xtal.Structure]: The parent structure, with lattice 
        :math:`L_{1}`."""

        self.parent_prim: Optional[xtal.Prim] = None
        """Optional[xtal.Prim]: The `xtal.Prim` used to represent the parent structure
        in the search methods."""

        self.child_structure: Optional[xtal.Structure] = None
        """Optional[xtal.Structure]: The child structure, with lattice :math:`L_{2}`."""

        self.child_T: Optional[np.ndarray] = None
        """Optional[np.ndarray]: The transformation matrix :math:`T_{2}` used to create
        a superstructure of the child for input to the search methods."""

        self.child_superstructure: Optional[xtal.Structure] = None
        """Optional[xtal.Structure]: The child superstructure, with lattice
        :math:`L_{2} T_{2}`."""

        self.scored_structure_mapping: Optional[mapinfo.ScoredStructureMapping] = None
        """Optional[mapinfo.ScoredStructureMapping]: The scored structure mapping
        between the parent structure and the child superstructure."""

        self.dedup_chain_orbit: Optional[list[list[xtal.Structure]]] = None
        """Optional[list[list[xtal.Structure]]]: A list of lists of structures, where
        `chain_orbit[i][0]` is the parent and `chain_orbit[i][-1]` is the child in the
        :math:`i`-th chain of interpolated structures in the orbit of equivalent 
        chains, put into primitive, canonical form for use in deduplication."""


def update_options_to_next_n_atoms(
    self,
    options: StructureMappingSearchOptions,
    results_dir: pathlib.Path,
):
    n_atoms_parent = len(self.parent.atom_type())
    n_atoms_child = len(self.child.atom_type())
    n_atoms_lcm = math.lcm(n_atoms_parent, n_atoms_child)

    last_max_n_atoms = None
    if len(self.options_history):
        if self.options_history[-1].max_n_atoms is None:
            last_max_n_atoms = n_atoms_lcm
        else:
            last_max_n_atoms = self.options_history[-1].max_n_atoms

    if last_max_n_atoms is None:
        next_max_n_atoms = n_atoms_lcm
    else:
        next_max_n_atoms = 0
        while next_max_n_atoms <= last_max_n_atoms:
            next_max_n_atoms += n_atoms_lcm

    options.max_n_atoms = next_max_n_atoms
    options.min_n_atoms = next_max_n_atoms


class StructureMappingSearch:
    """Search for mappings between superstructures of parent and child structures.

    Find structure mappings of the type:

    - parent and child structures have matching atom types and fractions
    - no vacancies, no parent sites with >1 allowed atom type
    - enable mean displacement removal


    To do so, this method finds lattice mappings of the type:

    .. math::

        F L_1 T_1 N = L_2 T_2

    where:

    - :math:`T_2` is a shape=(3,3) integer transformation matrix that generates a
      superlattice of the child lattice :math:`L_2`
    - other variables are defined as in the class
      :class:`libcasm.mapping.info.LatticeMapping`, using :math:`T_1` for :math:`T`.


    Notes
    -----

    This mapping search is limited to the case where the parent and child structures:

    - have the same atom types
    - have the same atom fractions

    Constraints that can be applied to the search:

    - min / max number of atoms
    - parent / child supercells used
    - min / max total cost of the mapping
    - min / max cost of the lattice mapping
    - lattice mapping reorientation range
    - k-best mappings to keep

    Other options include:

    - Choice of mapping cost methods:

      - Lattice mapping cost: "isotropic_strain_cost" or "symmetry_breaking_strain_cost"
      - Atom mapping cost: "isotropic_disp_cost" or "symmetry_breaking_disp_cost"
      - Lattice cost weight: The fraction of the total cost that is due to the lattice
        mapping cost. The remainder is due to the atom mapping cost.

    - Choice of interpolation factors used for deduplication

    """

    def __init__(
        self,
        opt: StructureMappingSearchOptions,
    ):
        self.opt: StructureMappingSearchOptions = opt
        """StructureMappingSearchOptions: Options for the search."""

    def _get_max_n_atoms(self, parent: xtal.Structure, child: xtal.Structure):
        """Get the maximum supercell size of the child structure based on the parent
        structure.

        If `child_max_supercell_size` is not set, the maximum supercell size is set to
        the least common multiple of the number of atoms in the child and parent
        structures.
        """
        if self.opt.max_n_atoms is not None:
            return self.opt.max_n_atoms

        n_atoms_parent = len(parent.atom_type())
        n_atoms_child = len(child.atom_type())
        return math.lcm(n_atoms_parent, n_atoms_child)

    def _get_child_to_parent_vol(
        self,
        parent: xtal.Structure,
        child: xtal.Structure,
    ):
        max_n_atoms = self._get_max_n_atoms(parent, child)
        child_n_atoms = len(child.atom_type())
        parent_n_atoms = len(parent.atom_type())

        child_to_parent_vol = {}
        child_vol = 1
        while child_vol * child_n_atoms <= max_n_atoms:
            child_superstructure_n_atoms = child_n_atoms * child_vol
            _vol = child_superstructure_n_atoms / parent_n_atoms

            # if parent_vol is integer, then it is a valid supercell size:
            if _vol.is_integer():
                child_to_parent_vol[child_vol] = int(_vol)

            child_vol += 1

        return child_to_parent_vol

    def _enable_symmetry_breaking_atom_cost(self):
        """Check if symmetry breaking atom cost is enabled based on the options."""
        return self.opt.atom_mapping_cost_method == "symmetry_breaking_disp_cost"

    def _atom_cost_f(self):
        """Get the atom cost function based on the options."""
        if self.opt.atom_mapping_cost_method == "isotropic_disp_cost":
            return mapsearch.IsotropicAtomCost()
        elif self.opt.atom_mapping_cost_method == "symmetry_breaking_disp_cost":
            return mapsearch.SymmetryBreakingAtomCost()
        else:
            raise ValueError(
                f"Unknown atom mapping cost method: {self.opt.atom_mapping_cost_method}"
            )

    def _total_cost_f(self):
        """Get the total cost function based on the options."""
        return mapsearch.WeightedTotalCost(
            lattice_cost_weight=self.opt.lattice_cost_weight
        )

    def validate(
        self,
        parent: xtal.Structure,
        child: xtal.Structure,
    ) -> None:
        """Raise if atom types or fractions differ between parent and child."""

        # Check atom types and stoichiometry
        parent_atom_types, parent_counts = np.unique(
            parent.atom_type(), return_counts=True
        )
        total_atoms = np.sum(parent_counts)
        parent_atom_frac = parent_counts / total_atoms

        child_atom_types, child_counts = np.unique(
            child.atom_type(), return_counts=True
        )
        total_atoms = np.sum(child_counts)
        child_atom_frac = child_counts / total_atoms

        if (parent_atom_types != child_atom_types).any():
            print("Error: Parent atom types differs from child atom types")
            print(f"- Parent atom types: {parent_atom_types}")
            print(f"- Child atom types: {child_atom_types}")
            print()
            print("Stopping")
            sys.exit(1)

        if not np.allclose(parent_atom_frac, child_atom_frac):
            print("Error: Parent and child structures have different atom fractions")
            print(f"- Atom types: {parent_atom_types}")
            print(f"- Parent atom fraction: {parent_atom_frac}")
            print(f"- Child atom fraction: {child_atom_frac}")
            print()
            print("Stopping")
            sys.exit(1)

        if self.opt.forced_on is not None:
            _allowed = [list([x]) for x in parent.atom_type()]
            _child_types = child.atom_type()
            for parent_site_index, child_atom_index in self.opt.forced_on.items():
                child_type = _child_types[child_atom_index]
                if child_type not in _allowed[parent_site_index]:
                    msgs.invalid_forced_on_values_error(
                        parent_site_index=parent_site_index,
                        child_atom_index=child_atom_index,
                        child_type=child_type,
                        allowed_types=_allowed[parent_site_index],
                    )

        if self.opt.fix_parent:
            child_n_atoms = len(child.atom_type())
            parent_n_atoms = len(parent.atom_type())
            if child_n_atoms != parent_n_atoms:
                msgs.invalid_fix_parent_error()

        else:
            # Print notice if parent or child are not primitive, and write the
            # primitive structures
            primitive_parent = xtal.make_primitive_structure(parent)
            if len(primitive_parent.atom_type()) != len(parent.atom_type()):
                safe_dump(
                    xtal.pretty_json(primitive_parent.to_dict()),
                    path="parent.primitive.json",
                    force=True,
                    quiet=True,
                )
                msgs.primitive_parent_notice()

            primitive_child = xtal.make_primitive_structure(child)
            if len(primitive_child.atom_type()) != len(child.atom_type()):
                safe_dump(
                    xtal.pretty_json(primitive_child.to_dict()),
                    path="child.primitive.json",
                    force=True,
                    quiet=True,
                )
                msgs.primitive_child_notice()

        if self.opt.child_transformation_matrix_to_super_list is None:
            # Validate the min/max number of atoms
            if self.opt.min_n_atoms < 1:
                msgs.invalid_min_n_atoms_error(min_n_atoms=self.opt.min_n_atoms)

            _max_n_atoms = self._get_max_n_atoms(parent, child)
            if _max_n_atoms < self.opt.min_n_atoms:
                computed_msg = (
                    "(computed from lcm of atom counts)"
                    if self.opt.max_n_atoms is None
                    else ""
                )
                msgs.invalid_max_n_atoms_error(
                    min_n_atoms=self.opt.min_n_atoms,
                    max_n_atoms=_max_n_atoms,
                    computed_msg=computed_msg,
                )

        # Validate lattice mapping cost method
        if self.opt.lattice_mapping_cost_method not in [
            "isotropic_strain_cost",
            "symmetry_breaking_strain_cost",
        ]:
            msgs.invalid_lattice_mapping_cost_method_error(
                self.opt.lattice_mapping_cost_method
            )

        # Validate atom mapping cost method
        if self.opt.atom_mapping_cost_method not in [
            "isotropic_disp_cost",
            "symmetry_breaking_disp_cost",
        ]:
            msgs.invalid_atom_mapping_cost_method_error(
                self.opt.atom_mapping_cost_method
            )

        # Validate that deduplication_interpolation_factors is a list of floats:
        dedup_factors = self.opt.deduplication_interpolation_factors
        if not isinstance(dedup_factors, list) or not all(
            isinstance(factor, float) for factor in dedup_factors
        ):
            msgs.invalid_deduplication_interpolation_factors_error(dedup_factors)

    def _make_T_pairs(
        self,
        parent: xtal.Structure,
        child: xtal.Structure,
        parent_prim: casmconfig.Prim,
        min_n_atoms: int,
        max_n_atoms: int,
        child_T_list: Optional[list[np.ndarray]] = None,
        parent_T_list: Optional[list[np.ndarray]] = None,
    ):
        """Make a list of (T_child, T_parent) pairs for the search.

        Parameters
        ----------
        parent : xtal.Structure
            The parent structure.
        child : xtal.Structure
            The child structure.
        parent_prim : casmconfig.Prim
            The primitive parent structure.
        min_n_atoms : int
            The minimum number of atoms in the superstructures that should be included
            in the search.
        max_n_atoms : int
            The maximum number of atoms in the superstructures that should be included
            in the search.
        child_T_list : Optional[list[np.ndarray]] = None
            For the child superstructures, a list of transformation matrices
            :math:`T_{2}` to use. If None, the child superstructures are enumerated
            based on the `min_n_atoms` and `max_n_atoms` options.
        parent_T_list : Optional[list[np.ndarray]] = None
            For the parent superstructures, a list of transformation matrices
            :math:`T_{1}` to use. If None, the parent superstructures are enumerated
            based on the `min_n_atoms` and `max_n_atoms` options.

        Returns
        -------
        T_pairs: list[tuple[np.ndarray, np.ndarray]]
            List of (T_child, T_parent) pairs.

        """
        # Results, list of (T_child, T_parent) pairs
        T_pairs = []

        # Parameters
        child_crystal_point_group = xtal.make_structure_crystal_point_group(child)
        child_n_atoms = len(child.atom_type())
        child_to_parent_vol = self._get_child_to_parent_vol(
            parent=parent,
            child=child,
        )

        # If child_T_list is not provided, enumerate the child supercells
        if child_T_list is None:
            child_T_list = []

            child_superlattices = xtal.enumerate_superlattices(
                unit_lattice=child.lattice(),
                point_group=child_crystal_point_group,
                max_volume=mthds.floordiv(max_n_atoms, child_n_atoms),
                min_volume=mthds.ceildiv(min_n_atoms, child_n_atoms),
            )
            for child_superlattice in child_superlattices:
                child_T_list.append(
                    xtal.make_transformation_matrix_to_super(
                        unit_lattice=child.lattice(),
                        superlattice=child_superlattice,
                    )
                )

        # For each child superstructure...
        for child_T in child_T_list:
            child_vol = int(round(np.linalg.det(child_T)))

            # If no valid parent volume, continue
            if child_vol not in child_to_parent_vol:
                continue
            parent_vol = child_to_parent_vol[child_vol]

            # Get the list of valid parent supercells
            restricted_parent_T_list = []

            # If parent_T_list is not provided, enumerate the parent supercells
            if parent_T_list is None:
                parent_superlattices = xtal.enumerate_superlattices(
                    unit_lattice=parent.lattice(),
                    point_group=parent_prim.crystal_point_group.elements,
                    max_volume=parent_vol,
                    min_volume=parent_vol,
                )
                for parent_superlattice in parent_superlattices:
                    restricted_parent_T_list.append(
                        xtal.make_transformation_matrix_to_super(
                            unit_lattice=parent.lattice(),
                            superlattice=parent_superlattice,
                        )
                    )

            # If parent_T_list is provided, filter the parent supercells
            else:
                for parent_T in parent_T_list:
                    if int(round(np.linalg.det(parent_T))) == parent_vol:
                        restricted_parent_T_list.append(parent_T)

            # Add the (child_T, parent_T) pairs
            for parent_T in restricted_parent_T_list:
                T_pairs.append((child_T, parent_T))

        return T_pairs

    def __call__(
        self,
        parent: xtal.Structure,
        parent_prim: Optional[casmconfig.Prim],
        child: xtal.Structure,
        results_dir: pathlib.Path,
        merge: bool = False,
    ):
        """Perform the structure mapping search.

        Parameters
        ----------
        parent : xtal.Structure
            The parent structure.
        parent_prim : Optional[casmconfig.Prim]
            The parent primitive structure. If both `parent` and `parent_prim` are
            provided, the `parent` structure is used to create the fix the supercell
            being mapped to with the "fix_parent" option.
        child : xtal.Structure
            The child structure.
        results_dir : pathlib.Path
            The directory to write the results to. If the directory already exists,
            the program exits with an error.
        merge: bool = False
            If True, merge the results with existing results in the directory. If False,
            exit with error if the directory already exists.
        """
        alloy = False
        if parent_prim is None:
            parent_prim = casmconfig.Prim(
                xtal.Prim.from_atom_coordinates(structure=parent)
            )
        else:
            alloy = True
            if not self.opt.fix_parent:
                raise NotImplementedError(
                    "Mapping to a prim is only supported with the --fix-parent option."
                )

        self._supercell_set = casmconfig.SupercellSet(prim=parent_prim)

        search_results = []
        uuids = []
        chain_orbits = []

        if results_dir.exists():
            if merge is False:
                msgs.results_dir_exists_error(results_dir=results_dir)
                sys.exit(1)
            else:
                data = read_required(results_dir / "mappings.json")

                # Validate same parent and child structures:
                if alloy is False:
                    _last_parent = xtal.Structure.from_dict(data.get("parent"))
                    if not parent.is_equivalent_to(_last_parent):
                        msgs.different_parent_error()
                    _last_child = xtal.Structure.from_dict(data.get("child"))
                    if not child.is_equivalent_to(_last_child):
                        msgs.different_child_error()
                else:
                    # TODO validation
                    pass

                search_results = [
                    mapinfo.ScoredStructureMapping.from_dict(
                        data=x, prim=parent_prim.xtal_prim
                    )
                    for x in data["mappings"]
                ]
                uuids = data.get("uuids", [])

                options_data = read_required(results_dir / "options_history.json")
                last_options = StructureMappingSearchOptions.from_dict(options_data[-1])

                if (
                    self.opt.lattice_mapping_cost_method
                    != last_options.lattice_mapping_cost_method
                ):
                    msgs.different_lattice_mapping_cost_method_error()
                if (
                    self.opt.atom_mapping_cost_method
                    != last_options.atom_mapping_cost_method
                ):
                    msgs.different_atom_mapping_cost_method_error()
                if not math.isclose(
                    self.opt.lattice_cost_weight,
                    last_options.lattice_cost_weight,
                    abs_tol=1e-5,
                ):
                    msgs.different_lattice_cost_weight_error()

        if alloy is False:
            self.validate(parent, child)
        else:
            # TODO
            pass

        ## Parameters
        if alloy is False:
            _max_n_atoms = self._get_max_n_atoms(parent, child)
        else:
            _max_n_atoms = mthds.get_max_n_atoms_for_parent_prim(
                max_n_atoms=self.opt.max_n_atoms,
                child=child,
            )
        _min_n_atoms = self.opt.min_n_atoms
        _child_T_list = self.opt.child_transformation_matrix_to_super_list
        _parent_T_list = self.opt.parent_transformation_matrix_to_super_list
        _enable_symmetry_breaking_atom_cost = self._enable_symmetry_breaking_atom_cost()
        _total_min_cost = self.opt.total_min_cost
        _total_max_cost = self.opt.total_max_cost
        _total_k_best = self.opt.total_k_best
        _no_remove_mean_displacement = self.opt.no_remove_mean_displacement
        _enable_remove_mean_displacement = not self.opt.no_remove_mean_displacement
        _lattice_cost_weight = self.opt.lattice_cost_weight
        _lattice_mapping_min_cost = self.opt.lattice_mapping_min_cost
        _lattice_mapping_max_cost = self.opt.lattice_mapping_max_cost
        _lattice_mapping_cost_method = self.opt.lattice_mapping_cost_method
        _lattice_mapping_k_best = self.opt.lattice_mapping_k_best
        _lattice_mapping_reorientation_range = (
            self.opt.lattice_mapping_reorientation_range
        )
        _atom_cost_f = self._atom_cost_f()
        _forced_on = self.opt.forced_on if self.opt.forced_on is not None else {}
        _forced_off = self.opt.forced_off if self.opt.forced_off is not None else []
        _total_cost_f = self._total_cost_f()
        _cost_tol = self.opt.cost_tol

        ## Fixed parameters
        _infinity = 1e20
        _atom_to_site_cost_f = mapsearch.make_atom_to_site_cost

        ## Create a parent structure search data object.
        prim_search_data = mapsearch.PrimSearchData(
            prim=parent_prim.xtal_prim,
            enable_symmetry_breaking_atom_cost=_enable_symmetry_breaking_atom_cost,
        )
        init_child_structure_data = mapsearch.StructureSearchData(
            lattice=child.lattice(),
            atom_coordinate_cart=child.atom_coordinate_cart(),
            atom_type=child.atom_type(),
            override_structure_factor_group=None,
        )

        if self.opt.fix_parent:
            if alloy is False:
                I_matrix = np.eye(3, dtype="int")
                T_pairs = [(I_matrix, I_matrix)]
            else:
                # If fixing the parent, we only need one pair of transformation matrices
                # (identity for both child and parent).
                I_matrix = np.eye(3, dtype="int")
                T_parent = xtal.make_transformation_matrix_to_super(
                    superlattice=parent.lattice(),
                    unit_lattice=parent_prim.xtal_prim.lattice(),
                )
                T_pairs = [(I_matrix, T_parent)]
        else:
            if alloy is False:
                # Get a list of (T_child, T_parent) pairs
                T_pairs = self._make_T_pairs(
                    parent=parent,
                    child=child,
                    parent_prim=parent_prim,
                    min_n_atoms=_min_n_atoms,
                    max_n_atoms=_max_n_atoms,
                    child_T_list=_child_T_list,
                    parent_T_list=_parent_T_list,
                )
            else:
                raise NotImplementedError(
                    "Alloy structures are only supported with the --fix-parent option."
                )

        total = len(T_pairs)
        print(f"Beginning search over {total} parent / child superstructure pairs...")
        print()
        print("Search results so far:")
        print()
        print()
        sys.stdout.flush()

        last_child_T = None
        child_structure_data = None

        n_atoms = 0
        if len(T_pairs):
            # Get the number of atoms in the child structure
            child_vol = int(round(np.linalg.det(T_pairs[0][0])))
            child_n_atoms = len(child.atom_type())
            n_atoms = child_n_atoms * child_vol
        min_total_cost = 0.0
        max_total_cost = 0.0

        for i_pair, _pair in enumerate(T_pairs):

            child_T, parent_T = _pair
            if last_child_T is None or not np.allclose(child_T, last_child_T):
                child_structure_data = mapsearch.make_superstructure_data(
                    prim_structure_data=init_child_structure_data,
                    transformation_matrix_to_super=child_T,
                )

            child_vol = int(round(np.linalg.det(child_T)))
            child_n_atoms = len(child.atom_type())
            n_atoms = child_n_atoms * child_vol

            # Create a MappingSearch object.
            # This will hold a queue of possible mappings,
            # sorted by cost, as we generate them.
            search = mapsearch.MappingSearch(
                min_cost=_total_min_cost,
                max_cost=_total_max_cost,
                k_best=_total_k_best,
                atom_cost_f=_atom_cost_f,
                total_cost_f=_total_cost_f,
                atom_to_site_cost_f=_atom_to_site_cost_f,
                enable_remove_mean_displacement=_enable_remove_mean_displacement,
                infinity=_infinity,
                cost_tol=_cost_tol,
            )
            atom_mapping_step = AtomMappingStep(
                search=search,
                prim_search_data=prim_search_data,
                child_structure_data=child_structure_data,
                opt=AtomMappingOptions(
                    no_remove_mean_displacement=_no_remove_mean_displacement,
                    forced_on=_forced_on,
                    forced_off=_forced_off,
                ),
            )

            lattice_mapping_step = LatticeMappingStep(
                prim_search_data=prim_search_data,
                child_structure_data=child_structure_data,
                lattice_cost_weight=_lattice_cost_weight,
                opt=LatticeMappingOptions(
                    lattice_mapping_cost_method=_lattice_mapping_cost_method,
                    lattice_mapping_min_cost=_lattice_mapping_min_cost,
                    lattice_mapping_max_cost=_lattice_mapping_max_cost,
                    lattice_mapping_k_best=_lattice_mapping_k_best,
                    lattice_mapping_reorientation_range=_lattice_mapping_reorientation_range,
                    fix_parent=self.opt.fix_parent,
                    cost_tol=_cost_tol,
                ),
            )

            # Might be able to tighten lattice max cost limit:
            curr_total_max_cost = _total_max_cost
            if len(search_results) >= _total_k_best:
                curr_total_max_cost = search_results[-1].total_cost()

            lattice_mappings = lattice_mapping_step(
                parent_T=parent_T,
                curr_total_max_cost=curr_total_max_cost,
            )

            for scored_lattice_mapping in lattice_mappings:
                atom_mapping_step(scored_lattice_mapping)

            while search.size():
                search.partition()

            search_results, uuids, chain_orbits = self.add_new_results(
                new_results=search.results().data(),
                existing_results=search_results,
                uuids=uuids,
                chain_orbits=chain_orbits,
                parent=parent,
                child=child,
                parent_prim=parent_prim,
                k_best=_total_k_best,
                cost_tol=_cost_tol,
            )

            self.write_results(
                search_results=search_results,
                uuids=uuids,
                parent=parent,
                child=child,
                parent_prim=parent_prim,
                results_dir=results_dir,
            )

            if len(search_results) > 0:
                min_total_cost = search_results[0].total_cost()
                max_total_cost = search_results[-1].total_cost()

            # Delete the last line
            sys.stdout.write("\033[F")  # Move cursor up one line
            sys.stdout.write("\033[K")  # Clear the line
            sys.stdout.flush()

            print(
                (
                    f"Pair: {i_pair + 1} / {total} (#atoms: {n_atoms}), "
                    f"MinTotalCost: {min_total_cost:.5f}, "
                    f"MaxTotalCost: {max_total_cost:.5f}, "
                    f"#mappings: {len(search_results)}"
                )
            )
            sys.stdout.flush()
            # pbar.update(1)

        print("DONE")
        print()
        sys.stdout.flush()
        print(f"# Results: {len(search_results)}\n")
        sys.stdout.flush()

        self.tabulate_results(
            search_results=search_results,
            uuids=uuids,
            parent=parent,
            child=child,
            parent_prim=parent_prim,
        )

        # Write the options history
        self.write_options_history(results_dir=results_dir)

        return 0

    def add_new_results(
        self,
        new_results: list[mapinfo.ScoredStructureMapping],
        existing_results: list[mapinfo.ScoredStructureMapping],
        uuids: list[str],
        chain_orbits: list[list[xtal.Structure]],
        parent: xtal.Structure,
        child: xtal.Structure,
        parent_prim: casmconfig.Prim,
        k_best: int,
        cost_tol: float,
    ) -> tuple[
        list[mapinfo.ScoredStructureMapping],
        list[str],
        list[list[xtal.Structure]],
    ]:
        """Add new results to the existing search results, deduplicating them.

        Parameters
        ----------
        new_results : list[libcasm.mapping.info.ScoredStructureMapping]
            The new results to add to the existing search results.
        existing_results : list[libcasm.mapping.info.ScoredStructureMapping]
            The existing search results to which the new results will be added.
        uuids : list[str]
            The UUIDs of the existing search results.
        chain_orbits : list[list[xtal.Structure]]
            The chain orbits of the existing search results.
        parent : xtal.Structure
            The parent structure.
        child : xtal.Structure
            The child structure.
        parent_prim : casmconfig.Prim
            The parent structure, as a Prim.
        k_best : int
            The number of best results to keep after deduplication. Any approximate ties
            will also be kept.
        cost_tol : float
            The tolerance for comparing costs.

        Returns
        -------
        search_results : list[libcasm.mapping.info.ScoredStructureMapping]
            The updated list of search results after deduplication.
        uuids : list[str]
            The updated list of UUIDs corresponding to the search results.
        chain_orbits : list[list[xtal.Structure]]
            The updated list of chain orbits corresponding to the search results.

        """
        search_results = existing_results

        # Deduplicate the new results
        f_chain = self.opt.deduplication_interpolation_factors

        def make_chain(structure_mapping):
            return mthds.make_primitive_chain(
                parent_lattice=parent_prim.xtal_prim.lattice(),
                child=child,
                structure_mapping=structure_mapping,
                f_chain=f_chain,
            )

        def make_orbit(chain_prototype):
            return mthds.make_chain_orbit(
                chain_prototype=chain_prototype,
                parent_prim=parent_prim,
            )

        while len(chain_orbits) < len(search_results):
            smap = search_results[len(chain_orbits)]
            chain_orbits.append(make_orbit(make_chain(smap)))
            uuids.append(str(uuid.uuid4()))

        if len(new_results) == 0:
            return search_results, uuids, chain_orbits

        for i, smap_new in enumerate(new_results):
            primitive_chain = make_chain(smap_new)

            # Check for duplicates:
            found_duplicate = False
            i_duplicate = 0
            for smap_existing, chain_orbit_existing in zip(
                search_results, chain_orbits
            ):
                if mthds.chain_is_in_orbit(primitive_chain, chain_orbit_existing):
                    found_duplicate = True
                    break
                i_duplicate += 1

            if found_duplicate:
                smap_existing = search_results[i_duplicate]
                scel_size_new = mthds.parent_supercell_size(smap_new)
                scel_size_existing = mthds.parent_supercell_size(smap_existing)

                prefer_new = False
                if scel_size_new < scel_size_existing:
                    prefer_new = True

                # prefer smaller volume mappings
                if prefer_new:
                    search_results[i_duplicate] = smap_new
                    uuids[i_duplicate] = str(uuid.uuid4())
                    chain_orbits[i_duplicate] = make_orbit(primitive_chain)
                else:
                    continue
            else:
                search_results.append(smap_new)
                uuids.append(str(uuid.uuid4()))
                chain_orbits.append(make_orbit(primitive_chain))

        # Sort the search results and chain orbits, by total cost
        sys.stdout.flush()
        isorted = [
            x[0]
            for x in sorted(enumerate(search_results), key=lambda x: x[1].total_cost())
        ]
        search_results = [search_results[i] for i in isorted]
        uuids = [uuids[i] for i in isorted]
        chain_orbits = [chain_orbits[i] for i in isorted]

        # Keep only the k-best results
        if len(search_results) > k_best:
            next_index = k_best
            while next_index < len(search_results):
                next_cost = search_results[next_index].total_cost()
                if math.isclose(
                    search_results[k_best - 1].total_cost(), next_cost, abs_tol=cost_tol
                ):
                    next_index += 1
                else:
                    break

            search_results = search_results[:(next_index)]
            uuids = uuids[:(next_index)]
            chain_orbits = chain_orbits[:(next_index)]

        return search_results, uuids, chain_orbits

    def write_results(
        self,
        search_results: list[mapinfo.ScoredStructureMapping],
        uuids: list[str],
        parent: xtal.Structure,
        child: xtal.Structure,
        parent_prim: casmconfig.Prim,
        results_dir: pathlib.Path,
    ) -> None:
        """Write the results of the search."""
        data = {
            "parent": parent.to_dict(),
            "child": child.to_dict(),
            "parent_prim": parent_prim.to_dict(),
            "mappings": [smap.to_dict() for smap in search_results],
            "uuids": [x for x in uuids],
        }
        safe_dump(
            data,
            path=results_dir / "mappings.json",
            force=True,
            quiet=True,
        )

    def write_options_history(
        self,
        results_dir: pathlib.Path,
    ) -> None:
        options = read_optional(results_dir / "options_history.json", default=[])
        options.append(self.opt.to_dict())
        safe_dump(
            options,
            path=results_dir / "options_history.json",
            force=True,
            quiet=True,
        )

    def tabulate_results(
        self,
        search_results: list[mapinfo.ScoredStructureMapping],
        uuids: list[str],
        parent: xtal.Structure,
        child: xtal.Structure,
        parent_prim: casmconfig.Prim,
    ) -> str:
        """Tabulate the results of the search."""

        prec = 5
        headers = [
            "Index",
            "TotCost",
            "LatCost",
            "AtmCost",
            "Parent Vol., Grp., #Ops",
            "Child Vol., Grp., #Ops",
            "Mult.",
            "UUID",
        ]
        f_chain = self.opt.deduplication_interpolation_factors
        child_prim = casmconfig.Prim(xtal.Prim.from_atom_coordinates(structure=child))

        data = []
        for i, scored_structure_mapping in enumerate(search_results):
            smap = scored_structure_mapping

            latmap = smap.lattice_mapping()
            T_parent = latmap.transformation_matrix_to_super()
            parent_volume = abs(int(round(np.linalg.det(T_parent))))
            T_child = mthds.make_child_transformation_matrix_to_super(
                parent_lattice=parent_prim.xtal_prim.lattice(),
                child_lattice=child.lattice(),
                structure_mapping=scored_structure_mapping,
            )
            child_volume = abs(int(round(np.linalg.det(T_child))))

            total_cost = f"{smap.total_cost():.{prec}f}"
            lattice_cost = f"{smap.lattice_cost():.{prec}f}"
            atom_cost = f"{smap.atom_cost():.{prec}f}"

            chain_orbit = mthds.make_primitive_chain_orbit(
                parent_prim=parent_prim,
                child=child,
                structure_mapping=smap,
                f_chain=f_chain,
            )
            mult = len(chain_orbit)

            parent_info = mthds.make_parent_supercell_info(
                structure_mapping=smap,
                parent_prim=parent_prim,
            )
            parent_grp = parent_info["spacegroup_type"]["international_short"]
            fg_size = parent_info["factor_group_size"]

            child_info = mthds.make_child_supercell_info(
                T_child=T_child,
                child_prim=child_prim,
            )
            child_grp = child_info["spacegroup_type"]["international_short"]
            child_fg_size = child_info["factor_group_size"]

            data.append(
                [
                    i,
                    total_cost,
                    lattice_cost,
                    atom_cost,
                    str(parent_volume) + ", " + parent_grp + ", " + str(fg_size),
                    str(child_volume) + ", " + child_grp + ", " + str(child_fg_size),
                    mult,
                    uuids[i],
                ]
            )

        print("Lattice cost method:", self.opt.lattice_mapping_cost_method)
        print("Atom cost method:", self.opt.atom_mapping_cost_method)
        print("Lattice cost weight:", self.opt.lattice_cost_weight)
        print(tabulate(data, headers=headers, tablefmt="grid"))
        print()


class AtomMappingOptions:
    def __init__(
        self,
        no_remove_mean_displacement: bool = False,
        forced_on: Optional[dict[int, int]] = None,
        forced_off: Optional[list[tuple[int, int]]] = None,
    ):
        self.no_remove_mean_displacement = no_remove_mean_displacement
        self.forced_on = forced_on
        self.forced_off = forced_off


# A LatticeMappingOptions class,
# for use by LatticeMappingStep.
class LatticeMappingOptions:

    def __init__(
        self,
        lattice_mapping_cost_method: str = "symmetry_breaking_strain_cost",
        lattice_mapping_min_cost: float = 0.0,
        lattice_mapping_max_cost: float = 1e20,
        lattice_mapping_k_best: int = 10,
        lattice_mapping_reorientation_range: int = 1,
        fix_parent: bool = False,
        cost_tol: float = 1e-5,
    ):
        self.lattice_mapping_cost_method = lattice_mapping_cost_method
        self.lattice_mapping_min_cost = lattice_mapping_min_cost
        self.lattice_mapping_max_cost = lattice_mapping_max_cost
        self.lattice_mapping_k_best = lattice_mapping_k_best
        self.lattice_mapping_reorientation_range = lattice_mapping_reorientation_range
        self.fix_parent = fix_parent
        self.cost_tol = cost_tol


class LatticeMappingStep:
    def __init__(
        self,
        prim_search_data: mapsearch.PrimSearchData,
        child_structure_data: mapsearch.StructureSearchData,
        lattice_cost_weight: float,
        opt: LatticeMappingOptions,
    ):
        self.prim_search_data = prim_search_data
        self.child_structure_data = child_structure_data
        self.lattice_cost_weight = lattice_cost_weight
        self.opt = opt

    def run_without_reorientation(
        self,
        parent_T,
    ):
        prim_search_data = self.prim_search_data
        child_structure_data = self.child_structure_data
        lattice_mapping_cost_method = self.opt.lattice_mapping_cost_method

        ###
        lattice_mapping = mapmethods.map_lattices_without_reorientation(
            lattice1=prim_search_data.prim_lattice(),
            lattice2=child_structure_data.lattice(),
            transformation_matrix_to_super=parent_T,
        )
        F = lattice_mapping.deformation_gradient()
        if lattice_mapping_cost_method == "isotropic_strain_cost":
            lattice_cost = mapinfo.isotropic_strain_cost(
                deformation_gradient=F,
            )
        elif lattice_mapping_cost_method == "symmetry_breaking_strain_cost":
            lattice_cost = mapinfo.symmetry_breaking_strain_cost(
                deformation_gradient=F,
                lattice1_point_group=prim_search_data.prim_crystal_point_group(),
            )
        else:
            raise ValueError(
                f"Unknown lattice mapping cost method: "
                f"{lattice_mapping_cost_method}"
            )
        return [
            mapinfo.ScoredLatticeMapping(
                lattice_cost=lattice_cost,
                lattice_mapping=lattice_mapping,
            )
        ]

    def run_with_reorientation(
        self,
        parent_T,
        curr_total_max_cost: float,
    ):
        prim_search_data = self.prim_search_data
        child_structure_data = self.child_structure_data
        lattice_cost_weight = self.lattice_cost_weight
        lattice_mapping_cost_method = self.opt.lattice_mapping_cost_method
        lattice_mapping_min_cost = self.opt.lattice_mapping_min_cost
        lattice_mapping_max_cost = self.opt.lattice_mapping_max_cost
        lattice_mapping_k_best = self.opt.lattice_mapping_k_best
        lattice_mapping_reorientation_range = (
            self.opt.lattice_mapping_reorientation_range
        )
        cost_tol = self.opt.cost_tol

        ###

        # # Might be able to tighten lattice max cost limit:
        # _curr_search_max = total_max_cost
        # if len(search_results) >= _total_k_best:
        #     _curr_search_max = search_results[-1].total_cost()
        #
        _curr_lattice_max = min(
            lattice_mapping_max_cost,
            curr_total_max_cost / lattice_cost_weight,
        )

        return mapmethods.map_lattices(
            lattice1=prim_search_data.prim_lattice(),
            lattice2=child_structure_data.lattice(),
            transformation_matrix_to_super=parent_T,
            lattice1_point_group=prim_search_data.prim_crystal_point_group(),
            lattice2_point_group=child_structure_data.structure_crystal_point_group(),
            min_cost=lattice_mapping_min_cost,
            max_cost=_curr_lattice_max,
            cost_method=lattice_mapping_cost_method,
            k_best=lattice_mapping_k_best,
            reorientation_range=lattice_mapping_reorientation_range,
            cost_tol=cost_tol,
        )

    def __call__(
        self,
        parent_T: np.ndarray,
        curr_total_max_cost: float = 1e20,
    ):
        if self.opt.fix_parent:
            return self.run_without_reorientation(parent_T=parent_T)

        else:
            return self.run_with_reorientation(
                parent_T=parent_T,
                curr_total_max_cost=curr_total_max_cost,
            )


class AtomMappingStep:
    def __init__(
        self,
        search: mapsearch.MappingSearch,
        prim_search_data: mapsearch.PrimSearchData,
        child_structure_data: mapsearch.StructureSearchData,
        opt: AtomMappingOptions,
    ):
        self.search = search
        self.prim_search_data = prim_search_data
        self.child_structure_data = child_structure_data
        self.opt = opt

    def __call__(
        self,
        scored_lattice_mapping,
    ):
        search = self.search
        prim_search_data = self.prim_search_data
        child_structure_data = self.child_structure_data
        no_remove_mean_displacement = self.opt.no_remove_mean_displacement
        forced_on = self.opt.forced_on
        forced_off = self.opt.forced_off

        # Make lattice mapping data
        lattice_mapping_data = mapsearch.LatticeMappingSearchData(
            prim_data=prim_search_data,
            structure_data=child_structure_data,
            lattice_mapping=scored_lattice_mapping,
        )

        # Check if 'forced_on' values are valid.
        if len(forced_on) > 0:
            _allowed = lattice_mapping_data.supercell_allowed_atom_types()
            _child_types = child_structure_data.atom_type()
            for parent_site_index, child_atom_index in forced_on.items():
                child_type = _child_types[child_atom_index]
                if child_type not in _allowed[parent_site_index]:
                    raise ValueError(
                        f"Invalid --forced-on values: "
                        f"child atom {child_atom_index} (type={child_type}) "
                        f"is not allowed to map to "
                        f"parent site {parent_site_index} "
                        f"(allowed types: {_allowed[parent_site_index]})"
                    )

        # Generate possible translations
        if no_remove_mean_displacement:
            # If mean displacement removal is disabled, then we need info
            # on which parent/atom mappings to force on. (We could also allow
            # generating every combination here.)
            if len(forced_on) == 0:
                raise ValueError(
                    "If --no-remove-mean-displacement is set, "
                    "the --forced-on option must be set."
                )
            # If forced_on is set, also use parent/child pairs to generate
            # trial translations
            trial_translations = []
            parent_cart = prim_search_data.prim_site_coordinate_cart()
            child_cart = lattice_mapping_data.atom_coordinate_cart_in_supercell()
            for parent_index, child_index in forced_on.items():
                trial_translations.append(
                    parent_cart[:, parent_index] - child_cart[:, child_index]
                )
        else:
            # Make a minimal set of trial translations
            trial_translations = mapsearch.make_trial_translations(
                lattice_mapping_data=lattice_mapping_data,
            )

        # For each combination of lattice mapping and translation,
        # make and insert a mapping solution (MappingNode)
        for trial_translation in trial_translations:
            search.make_and_insert_mapping_node(
                lattice_cost=scored_lattice_mapping.lattice_cost(),
                lattice_mapping_data=lattice_mapping_data,
                trial_translation_cart=trial_translation,
                forced_on=forced_on,
                forced_off=forced_off,
            )


#
#
# input: (T_child, T_parent)
# input: existing results
# input: parameters
#
# for each pair:
#     make MappingSearch(max_total_cost, ...)
#     make lattice mappings
#     for each lattice mapping:
#         make trial translations
#         for each trial translation:
#             make and insert mapping node
#     search
#     merge results
#     write results
#     update min/max cost
