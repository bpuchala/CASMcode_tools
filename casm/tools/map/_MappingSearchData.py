from typing import Optional, Union

import numpy as np

import libcasm.configuration as casmconfig
import libcasm.mapping.info as mapinfo
import libcasm.xtal as xtal
from casm.tools.shared.json_io import (
    safe_dump,
)

from . import messages as msgs
from . import methods as mthds
from ._StructureMappingSearchOptions import StructureMappingSearchOptions


class MappingSearchData:
    def __init__(
        self,
        parent: Union[xtal.Structure, casmconfig.Prim],
        child: xtal.Structure,
        options: Optional[StructureMappingSearchOptions] = None,
        mappings: list[mapinfo.ScoredStructureMapping] = None,
        uuids: list[str] = None,
        options_history: list[StructureMappingSearchOptions] = None,
    ):
        if mappings is None:
            mappings = []
        if uuids is None:
            uuids = []
        if options_history is None:
            options_history = []

        if isinstance(parent, xtal.Structure):
            parent_structure = parent
            # noinspection PyArgumentList
            xtal_prim = xtal.Prim.from_atom_coordinates(structure=parent)
            parent_prim = casmconfig.Prim(xtal_prim)
        elif isinstance(parent, casmconfig.Prim):
            parent_structure = None
            parent_prim = parent
        else:
            raise TypeError(
                "Error in MappingSearchData: `parent` must be either a "
                "libcasm.xtal.Structure or a libcasm.configuration.Prim."
            )

        self.parent_structure: Optional[xtal.Structure] = parent_structure
        """Optional[xtal.Structure]: The parent structure, with lattice 
        :math:`L_{1}`, if mapping to a particular structure."""

        self.parent_prim: casmconfig.Prim = parent_prim
        """casmconfig.Prim: The :class:`~libcasm.configuration.Prim` for the parent, 
        which determines allowed occupants on each basis site."""

        self.child: xtal.Structure = child
        """Optional[xtal.Structure]: The child structure, with lattice :math:`L_{2}`."""

        self.options: Optional[StructureMappingSearchOptions] = options
        """Optional[StructureMappingSearchOptions]: Options for the current search."""

        self.mappings: list[mapinfo.ScoredStructureMapping] = mappings
        """list[mapinfo.ScoredStructureMapping]: The list of scored structure mappings
        between the parent structure and the child superstructure."""

        self.uuids: list[str] = uuids
        """list[str]: A list of UUIDs for the mappings."""

        self.options_history: list[StructureMappingSearchOptions] = options_history
        """list[StructureMappingSearchOptions]: A history of options used by previous
        searches."""

    @property
    def parent_atom_types(self):
        """list[str]: The list of atom types in the parent structure, sorted."""
        if self.parent_structure is not None:
            parent_atom_types = set(self.parent_structure.atom_type())
        else:
            occ_dof = self.parent_prim.xtal_prim.occ_dof()
            parent_atom_types = {name for site_dof in occ_dof for name in site_dof}
        return sorted(list(parent_atom_types))

    @property
    def parent_atom_frac(self):
        """Optional[np.ndarray]: The fraction of each atom type in the parent, in
        order corresponding to `parent_atom_types`.

        If `parent_structure` is None, or has no atoms, the value is None."""
        if self.parent_structure is None:
            return None
        if len(self.parent_structure.atom_type()) == 0:
            return None
        _atom_types = self.parent_atom_types
        _atom_count = [0] * len(_atom_types)
        for atom_type in self.parent_structure.atom_type():
            _atom_count[_atom_types.index(atom_type)] += 1
        _atom_counts = np.array(_atom_count)
        total = np.sum(_atom_counts)
        return _atom_counts / total

    @property
    def min_atom_count_per_parent_unitcell(self):
        """np.array: The minimum number of atoms per parent unit cell of each type, in
        order corresponding to `parent_atom_types`."""
        _atom_types = self.parent_atom_types
        _atom_count = [0] * len(_atom_types)
        if self.parent_structure is not None:
            for atom_type in self.parent_structure.atom_type():
                _atom_count[_atom_types.index(atom_type)] += 1
        else:
            occ_dof = self.parent_prim.xtal_prim.occ_dof()
            for site_dof in occ_dof:
                if len(site_dof) == 1:
                    _atom_count[_atom_types.index(site_dof[0])] += 1
        return np.array(_atom_count)

    @property
    def max_atom_count_per_parent_unitcell(self):
        """np.array: The maximum number of atoms per parent unit cell of each type, in
        order corresponding to `parent_atom_types`."""
        _atom_types = self.parent_atom_types
        _atom_count = [0] * len(_atom_types)
        if self.parent_structure is not None:
            for atom_type in self.parent_structure.atom_type():
                _atom_count[_atom_types.index(atom_type)] += 1
        else:
            occ_dof = self.parent_prim.xtal_prim.occ_dof()
            for site_dof in occ_dof:
                for name in site_dof:
                    _atom_count[_atom_types.index(name)] += 1
        return np.array(_atom_count)

    @property
    def child_atom_count(self):
        """np.array: The number of atoms in the child of each type, in
        order corresponding to `child_atom_types`."""
        _atom_types = self.child_atom_types
        _atom_count = [0] * len(_atom_types)
        for atom_type in self.child_structure.atom_type():
            index = _atom_types.index(atom_type)
            if index >= 0:
                _atom_count[index] += 1
        return np.array(_atom_count)

    @property
    def child_atom_count_of_parent_types(self):
        """np.array: The number of atoms in the child of each parent type, in
        order corresponding to `parent_atom_types`."""
        _atom_types = self.parent_atom_types
        _atom_count = [0] * len(_atom_types)
        for atom_type in self.child_structure.atom_type():
            index = _atom_types.index(atom_type)
            if index >= 0:
                _atom_count[index] += 1
        return np.array(_atom_count)

    @property
    def child_atom_types(self):
        """list[str]: The list of atom types in the child structure, sorted."""
        child_atom_types = set(self.child.atom_type())
        return sorted(list(child_atom_types))

    @property
    def child_atom_frac(self):
        """Optional[np.ndarray]: The fraction of each atom type in the child, in
        order corresponding to `child_atom_types`.

        If `child` has no atoms, return None."""
        if len(self.child.atom_type()) == 0:
            return None
        _atom_types = self.child_atom_types
        _atom_count = [0] * len(_atom_types)
        for atom_type in self.child.atom_type():
            _atom_count[_atom_types.index(atom_type)] += 1
        _atom_counts = np.array(_atom_count)
        total = np.sum(_atom_counts)
        return _atom_counts / total

    def validate_atom_types(self):
        """Validate parent and child atom types are consistent

        If `parent_structure` is not None, check that the atom types in the parent
        structure and child structure are the same.

        If `parent_structure` is None, check that the atom types in the child structure
        are a subset of the atom types in the parent prim.

        If the check fails, print an error message and exit.

        """
        if self.parent_structure is not None:
            if set(self.parent_atom_types) != set(self.child_atom_types):
                msgs.atom_types_mismatch_error(
                    self.parent_atom_types, self.child_atom_types
                )
        else:
            if not set(self.child_atom_types).issubset(set(self.parent_atom_types)):
                msgs.atom_types_mismatch_error(
                    self.parent_atom_types,
                    self.child_atom_types,
                )

    def validate_atom_frac(self):
        """Validate the parent and child atom fractions are consistent.

        If `parent_structure` is not None, check that the atom fractions in the parent
        structure and child structure are the same.
        """
        if self.parent_structure is not None:
            if not np.allclose(self.parent_atom_frac, self.child_atom_frac, atol=1e-5):
                msgs.atom_fraction_mismatch_error(
                    self.parent_atom_frac, self.child_atom_frac
                )

    def validate_forced_on(self):
        """Validate `forced_on` values, if the parent_structure is given

        If `parent_structure` is not None, check that the `--forced-on` option maps
        parent and child atoms of the same type.

        If `parent_structure` is None, this method currently does nothing.

        The `forced_on` values are also validated before atom mapping.

        """
        if self.parent_structure is None:
            return
        _allowed = [list([x]) for x in self.parent_structure.atom_type()]
        _child_types = self.child.atom_type()
        for parent_site_index, child_atom_index in self.opt.forced_on.items():
            child_type = _child_types[child_atom_index]
            if child_type not in _allowed[parent_site_index]:
                msgs.invalid_forced_on_values_error(
                    parent_site_index=parent_site_index,
                    child_atom_index=child_atom_index,
                    child_type=child_type,
                    allowed_types=_allowed[parent_site_index],
                )

    def validate_fix_parent(self):
        """Validate `fix_parent` option, if the parent_structure is given

        If `parent_structure` is not None, check that the `--fix-parent` option is used
        only when the number of atoms in the parent structure is the same as the number
        of atoms in the child structure.

        If `parent_structure` is None, this method currently does nothing.

        """
        if self.parent_structure is None:
            return
        if self.opt.fix_parent:
            child_n_atoms = len(self.child.atom_type())
            parent_n_atoms = len(self.parent_structure.atom_type())
            if child_n_atoms != parent_n_atoms:
                msgs.invalid_fix_parent_error()

    def notify_if_non_primitive(self):
        """Print a notice if the parent or child is not primitive, and write
        the primitive form to a file.

        If the parent is not primitive, write its primitive form to
        `parent.primitive.json`. If parent_structure is not None, it is checked.
        Otherwise, the parent_prim is checked.

        If the child is not primitive, write its primitive form to
        `child.primitive.json`.
        """
        if self.parent_structure is not None:
            parent = self.parent_structure
            primitive_parent = xtal.make_primitive_structure(parent)
            if len(primitive_parent.atom_type()) != len(parent.atom_type()):
                safe_dump(
                    xtal.pretty_json(primitive_parent.to_dict()),
                    path="parent.primitive.json",
                    force=True,
                    quiet=True,
                )
                msgs.primitive_parent_notice()
        else:
            parent = self.parent_structure
            primitive_parent = xtal.make_primitive_prim(parent)
            if len(primitive_parent.occ_dof()) != len(parent.occ_dof()):
                safe_dump(
                    xtal.pretty_json(primitive_parent.to_dict()),
                    path="parent.primitive.json",
                    force=True,
                    quiet=True,
                )
                msgs.primitive_parent_notice()

        primitive_child = xtal.make_primitive_structure(self.child)
        if len(primitive_child.atom_type()) != len(self.child.atom_type()):
            safe_dump(
                xtal.pretty_json(primitive_child.to_dict()),
                path="child.primitive.json",
                force=True,
                quiet=True,
            )
            msgs.primitive_child_notice()

    def validate_n_atoms(self):
        """Validate the min_n_atoms and max_n_atoms options.

        If `child_transformation_matrix_to_super_list` is None, this does nothing
        because the user has requested which child supercells to try mapping.

        Otherwise, it checks that the `min_n_atoms` >= 1 and that `max_n_atoms` is
        greater than or equal to `min_n_atoms`. If the user does not specify
        `max_n_atoms` explicitly, it is computed using the least common multiple of the
        number of atoms in the parent and child structures.

        """
        if self.options.child_transformation_matrix_to_super_list is None:

            min_n_atoms = self.options.min_n_atoms
            max_n_atoms = mthds.get_max_n_atoms_for_parent_structure(
                max_n_atoms=self.options.max_n_atoms,
                parent_structure=self.parent_structure,
                child=self.child,
            )

            # Validate the min/max number of atoms
            if min_n_atoms < 1:
                msgs.invalid_min_n_atoms_error(min_n_atoms=min_n_atoms)

            if self.parent_structure is not None:

                if max_n_atoms < min_n_atoms:
                    computed_msg = (
                        "(computed from lcm of atom counts)"
                        if self.options.max_n_atoms is None
                        else ""
                    )
                    msgs.invalid_max_n_atoms_error(
                        min_n_atoms=min_n_atoms,
                        max_n_atoms=max_n_atoms,
                        computed_msg=computed_msg,
                    )

    def move_options_to_history(self):
        """Move the current options to the options history."""
        if self.options is not None:
            self.options_history.append(self.options)
            self.options = None

    def to_dict(self):
        """Convert the search data to a Python dictionary.

        Notes
        -----

        This does not move the current options to the options history. Use
        :func:`move_options_to_history` before calling this method if you want to
        include the current options in the history.

        Returns
        -------
        data: dict
            A Python dict representation of the search data.

        """
        return {
            "parent_structure": (
                self.parent_structure.to_dict() if self.parent_structure else None
            ),
            "parent_prim": self.parent_prim.to_dict(),
            "child": self.child.to_dict(),
            "options": self.options.to_dict() if self.options else None,
            "mappings": [mapping.to_dict() for mapping in self.mappings],
            "uuids": self.uuids,
            "options_history": [opt.to_dict() for opt in self.options_history],
        }

    @staticmethod
    def from_dict(
        self,
        data: dict,
    ):
        """Create a MappingSearchData object from a Python dictionary.

        Parameters
        ----------
        data: dict
            A Python dict representation of the search data.

        Returns
        -------
        search_data: MappingSearchData
            The MappingSearchData object created from the dictionary.

        """
        parent_structure = (
            xtal.Structure.from_dict(data["parent_structure"])
            if data["parent_structure"] is not None
            else None
        )
        parent_prim = casmconfig.Prim.from_dict(data["parent_prim"])
        child = xtal.Structure.from_dict(data["child"])
        options = (
            StructureMappingSearchOptions.from_dict(data["options"])
            if data["options"] is not None
            else None
        )
        mappings = [
            mapinfo.ScoredStructureMapping.from_dict(data=x, prim=parent_prim.xtal_prim)
            for x in data["mappings"]
        ]
        uuids = data["uuids"]
        options_history = [
            StructureMappingSearchOptions.from_dict(data=x)
            for x in data["options_history"]
        ]

        return MappingSearchData(
            parent=parent_structure or parent_prim,
            child=child,
            options=options,
            mappings=mappings,
            uuids=uuids,
            options_history=options_history,
        )
