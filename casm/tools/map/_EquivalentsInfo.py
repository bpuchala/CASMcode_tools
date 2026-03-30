import pathlib
from typing import Union

import numpy as np

import libcasm.clexulator as casmclex
import libcasm.configuration as casmconfig
import libcasm.mapping.methods as mapmethods
import libcasm.xtal as xtal
from casm.tools.shared.misc import pretty

from ._StructureMappingSearch import read_results
from .methods import (
    get_global_dof_space,
    get_local_dof_space,
)


class EquivalentsInfo:
    def __init__(
        self,
        results_dir: Union[pathlib.Path, str],
    ):
        self.results_dir = pathlib.Path(results_dir)

        if not self.results_dir.exists():
            raise ValueError(
                "Error in EquivalentsInfo: "
                f"Results directory {self.results_dir} does not exist."
            )

        ### Read existing results ###
        parent, child, parent_prim, search_results, uuids, options_history = (
            read_results(
                results_dir=self.results_dir,
            )
        )

        self.parent = parent
        """libcasm.xtal.Structure: The parent structure used for the mapping."""

        self.child = child
        """libcasm.xtal.Structure: The child structure used for the mapping."""

        self.parent_prim = parent_prim
        """libcasm.configuration.Prim: The primitive cell of the parent structure."""

        self.search_results = search_results
        """list of libcasm.tools.map._search.SearchResult: The results of the structure 
        mapping search."""

        self.uuids = uuids
        """list of str: The UUIDs of the structures in the search results."""

        self.options_history = options_history
        """list of dict: The options used."""

        # Local continuous degrees of freedom (DoF)
        disp_dof = xtal.DoFSetBasis("disp")  # Atomic displacement

        # Global continuous degrees of freedom (DoF)
        GLstrain_dof = xtal.DoFSetBasis(
            dofname="GLstrain",
            # axis_names=["e_{1}", "e_{2}", "e_{3}", "e_{4}", "e_{5}", "e_{6}"],
            # basis=xtal.make_symmetry_adapted_strain_basis(),
        )

        n_sites = parent_prim.xtal_prim.coordinate_frac().shape[1]
        parent_prim_strain_disp = casmconfig.Prim(
            xtal_prim=xtal.Prim(
                lattice=parent_prim.xtal_prim.lattice(),
                coordinate_frac=parent_prim.xtal_prim.coordinate_frac(),
                occ_dof=parent_prim.xtal_prim.occ_dof(),
                local_dof=[[disp_dof] for x in range(n_sites)],
                global_dof=[GLstrain_dof],
                title="prim_with_strain_disp",
            )
        )

        self.parent_prim_strain_disp = parent_prim_strain_disp
        """libcasm.configuration.Prim: The primitive cell of the parent structure, 
        with local atomic displacement and global strain degrees of freedom."""

        supercell_set = casmconfig.SupercellSet(
            prim=self.parent_prim_strain_disp,
        )
        self.supercell_set = supercell_set
        """libcasm.configuration.SupercellSet: The supercell set constructed from the
        parent primitive with strain and displacement degrees of freedom."""

        ### Attributes that are set in _update() ###

        self.superdupercell = None
        """libcasm.configuration.Supercell: The fully commensurate supercell with
        strain and displacement degrees of freedom."""

        self.mapped_superchild_configuration = None
        """libcasm.configuration.Configuration: The mapped child configuration, copied
        into the superdupercell."""

        self.equiv = None
        """list[libcasm.configuration.Configuration]: The configurations equivalent to
        the mapped child configuration, in the superdupercell."""

        self.strain_irreps = None
        """list[libcasm.irreps.IrrepInfo]: The irreducible representations (irreps) of
        the global strain degrees of freedom in the superdupercell."""

        self.strain_dof_space = None
        """libcasm.clexulator.DoFSpace: The strain modes in the superdupercell."""

        self.strain_values = None
        """list[np.ndarray]: The global strain degrees of freedom for the configurations
        equivalent to the mapped child configuration in the superdupercell"""

        self.strain_order_parameters = None
        """list[np.ndarray]: The strain order parameters for the configurations
        equivalent to the mapped child configuration in the superdupercell."""

        self.unique_strain_order_parameters = None
        """list[np.ndarray]: The unique strain order parameters for the configurations
        equivalent to the mapped child configuration in the superdupercell."""

        self.disp_irreps = None
        """list[libcasm.irreps.IrrepInfo]: The irreducible representations (irreps) of
        the local displacement degrees of freedom in the superdupercell."""

        self.disp_dof_space = None
        """libcasm.clexulator.DoFSpace: The displacement modes in the superdupercell."""

        self.disp_values = None
        """list[np.ndarray]: The displacement degrees of freedom for the configurations
        equivalent to the mapped child configuration in the superdupercell"""

        self.disp_order_parameters = None
        """list[np.ndarray]: The displacement order parameters for the configurations
        equivalent to the mapped child configuration in the superdupercell."""

        self.unique_disp_order_parameters = None
        """list[np.ndarray]: The unique displacement order parameters for the
        configurations equivalent to the mapped child configuration in the
        superdupercell."""

    def _update(
        self,
        mapping: mapmethods.ScoredStructureMapping,
    ):
        """Update attributes based on a selected structure mapping."""

        ### Make mapped child structure / configuration ###

        T = mapping.lattice_mapping().transformation_matrix_to_super()
        T_int = np.round(T).astype(int)
        supercell = casmconfig.Supercell(
            prim=self.parent_prim_strain_disp,
            transformation_matrix_to_super=T_int,
        )

        mapped_child_structure = mapmethods.make_mapped_structure(
            structure_mapping=mapping.interpolated(1.0),
            unmapped_structure=self.child,
        )
        mapped_child_configuration = casmconfig.Configuration.from_structure(
            prim=self.parent_prim_strain_disp,
            structure=mapped_child_structure,
            supercells=self.supercell_set,
        )

        ### Make fully commensurate superdupercell with strain and displacement DoF ###

        superduper_lattice = xtal.make_canonical_lattice(
            xtal.make_superduperlattice(
                lattices=[supercell.superlattice],
                mode="fully_commensurate",
                point_group=self.parent_prim_strain_disp.lattice_point_group.elements,
            )
        )
        superduper_T = xtal.make_transformation_matrix_to_super(
            superlattice=superduper_lattice,
            unit_lattice=self.parent_prim_strain_disp.xtal_prim.lattice(),
        )
        superdupercell = casmconfig.Supercell(
            prim=self.parent_prim_strain_disp,
            transformation_matrix_to_super=superduper_T,
        )
        self.superdupercell = superdupercell

        ### Equivalents in the superdupercell ###

        mapped_superchild_configuration = casmconfig.copy_configuration(
            motif=mapped_child_configuration,
            supercell=superdupercell,
        )
        self.mapped_superchild_configuration = mapped_superchild_configuration

        equiv = casmconfig.make_equivalent_configurations(
            configuration=mapped_superchild_configuration,
        )
        self.equiv = equiv

        ### Strain modes in the superdupercell ###
        irreps, strain_dof_space = get_global_dof_space(
            results_dir=self.results_dir,
            supercell=superdupercell,
            dof_key="GLstrain",
        )
        self.strain_irreps = irreps

        self.strain_dof_space = strain_dof_space

        ### Calculate the strain order parameters for the mapped children ###

        f_strain = casmclex.OrderParameter(
            dof_space=self.strain_dof_space,
        )
        f_strain.update(
            transformation_matrix_to_super=self.superdupercell.transformation_matrix_to_super,
            site_index_converter=self.superdupercell.site_index_converter,
            config_dof_values=self.mapped_superchild_configuration.dof_values,
        )

        equiv_strains = []
        for config in self.equiv:
            x = pretty(config.dof_values_vector(dof_space=self.strain_dof_space))
            equiv_strains.append(x)
        self.strain_values = equiv_strains

        strain_order_parameters = []
        unique_strain_order_parameters = []
        for i, config in enumerate(self.equiv):
            f_strain.set(config.dof_values)
            x = pretty(f_strain.value())
            strain_order_parameters.append(x)
            if not any(
                np.allclose(x, op, atol=1e-8) for op in unique_strain_order_parameters
            ):
                unique_strain_order_parameters.append(x)
        self.strain_order_parameters = strain_order_parameters

        self.unique_strain_order_parameters = unique_strain_order_parameters

        ### Displacement modes in the superdupercell ###
        disp_irreps, disp_dof_space = get_local_dof_space(
            results_dir=self.results_dir,
            supercell=self.superdupercell,
            dof_key="disp",
        )
        self.disp_irreps = disp_irreps

        self.disp_dof_space = disp_dof_space

        ### Calculate the displacement order parameters for the mapped children ###
        f_disp = casmclex.OrderParameter(
            dof_space=self.disp_dof_space,
        )
        f_disp.update(
            transformation_matrix_to_super=self.superdupercell.transformation_matrix_to_super,
            site_index_converter=self.superdupercell.site_index_converter,
            config_dof_values=self.mapped_superchild_configuration.dof_values,
        )

        equiv_displacements = []
        for config in self.equiv:
            x = pretty(config.dof_values_vector(dof_space=self.disp_dof_space))
            equiv_displacements.append(x)
        self.disp_values = equiv_displacements

        disp_order_parameters = []
        unique_disp_order_parameters = []
        for i, config in enumerate(self.equiv):
            f_disp.set(config.dof_values)
            x = pretty(f_disp.value())
            disp_order_parameters.append(x)
            if not any(
                np.allclose(x, op, atol=1e-8) for op in unique_disp_order_parameters
            ):
                unique_disp_order_parameters.append(x)

        self.disp_order_parameters = disp_order_parameters

        self.unique_disp_order_parameters = unique_disp_order_parameters

    def update_by_index(self, index: int):
        """Update attributes for a particular structure mapping, specified by index in
        the search results.

        Parameters
        ----------
        index: int
            The index of the structure mapping in the search results to update for.

        """
        self._update(mapping=self.search_results[index])

    def update_by_uuid(self, uuid: str):
        """Update attributes for a particular structure mapping, specified by UUID.

        Parameters
        ----------
        uuid: str
            The UUID of the structure mapping in the search results to update for.

        """
        index = self.uuids.index(uuid)
        self.update_by_index(index=index)

    def get_coordinates(
        self,
        dof_type: str,
        irrep_index: int,
        irrep_mode_index: int,
    ):
        """Get coordinates for the equivalent configurations by projecting the DoF
        values onto a particular irrep mode.

        Parameters
        ----------
        dof_key: str
            One of "strain" or "disp".
        irrep_index: int
            The index of the irrep to project onto.
        irrep_mode_index: int
            The index of the mode within the irrep to project onto.

        Returns
        -------
        coords: np.ndarray
            The coordinates for the equivalent configurations, obtained by projecting
            the DoF values onto the specified irrep mode.
        """
        if dof_type == "strain":
            values = self.strain_values
            irreps = self.strain_irreps
        elif dof_type == "disp":
            values = self.disp_values
            irreps = self.disp_irreps
        else:
            raise ValueError(f"Invalid dof_type: {dof_type}")

        projected = []
        for i, x in enumerate(values):
            equiv_proj = []
            for j, irrep in enumerate(irreps):
                x_proj = pretty(irrep.trans_mat.real @ x)
                equiv_proj.append(x_proj)
            projected.append(equiv_proj)

        coords = [proj[irrep_index][irrep_mode_index] for proj in projected]
        return np.array(coords)
