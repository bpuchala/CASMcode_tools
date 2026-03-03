import copy
import typing

import numpy as np

import libcasm.configuration as casmconfig
import libcasm.mapping.info as mapinfo
import libcasm.mapping.methods as mapmethods
import libcasm.xtal as xtal
from casm.tools.shared.conversions import to_kelvin, voigt_to_kelvin


def make_lattice_coordinates(
    lmap_ref: mapinfo.LatticeMapping,
    lmap_i: mapinfo.LatticeMapping,
    strain_converter: xtal.StrainConverter,
):
    """Project lattice strain onto path and perpendicular directions.

    Parameters
    ----------
    lmap_ref : libcasm.mapping.info.LatticeMapping
        The reference lattice mapping defining the path direction.
    lmap_i : libcasm.mapping.info.LatticeMapping
        The lattice mapping for the current image.
    strain_converter : libcasm.xtal.StrainConverter
        The strain converter to use for converting deformation gradients to strain.

    Returns
    -------
    x_lattice_path : float
        The lattice strain coordinate along the path direction.
    x_lattice_perp : float
        The lattice strain coordinate perpendicular to the path direction.
    """

    Q_ref, U_ref = xtal.StrainConverter.F_to_QU(lmap_ref.deformation_gradient())
    v_U_ref = (U_ref - np.eye(3)).flatten()
    v_U_ref_mag = np.linalg.norm(v_U_ref)
    v_U_ref_unit = v_U_ref / v_U_ref_mag

    Q_i, U_i = xtal.StrainConverter.F_to_QU(lmap_i.deformation_gradient())
    v_U_i = (U_i - np.eye(3)).flatten()

    # Ustrain along path direction:
    v_U_i_path = np.dot(v_U_i, v_U_ref_unit) * v_U_ref_unit
    v_U_i_path_mag = np.linalg.norm(v_U_i_path)
    x_lattice_path = v_U_i_path_mag / v_U_ref_mag

    # Ustrain perpendicular to path direction:
    v_U_i_perp = v_U_i - v_U_i_path
    v_U_i_perp_mag = np.linalg.norm(v_U_i_perp)
    x_lattice_perp = v_U_i_perp_mag / v_U_ref_mag

    return (x_lattice_path, x_lattice_perp)


def make_displacement_coordinates(
    amap_ref: mapinfo.AtomMapping,
    amap_i: mapinfo.AtomMapping,
):
    """Project atomic displacements onto path and perpendicular directions.

    Parameters
    ----------
    amap_ref : libcasm.mapping.info.AtomMapping
        The reference atom mapping defining the path direction.
    amap_i : libcasm.mapping.info.AtomMapping
        The atom mapping for the current image.

    Returns
    -------
    x_disp_path : float
        The displacement coordinate along the path direction.
    x_disp_perp : float
        The displacement coordinate perpendicular to the path direction.
    """

    # a 3xN array of displacements, unrolled as a 1D array:
    d_ref = amap_ref.displacement().flatten()
    d_ref_unit = d_ref / np.linalg.norm(d_ref)
    d_ref_mag = np.linalg.norm(d_ref)

    d_i = amap_i.displacement().flatten()

    # Displacement along path direction
    d_i_path = np.dot(d_i, d_ref_unit) * d_ref_unit
    d_i_path_mag = np.linalg.norm(d_i_path)
    x_disp_path = d_i_path_mag / d_ref_mag

    # Displacement perpendicular to path direction
    d_i_perp = d_i - d_i_path
    d_i_perp_mag = np.linalg.norm(d_i_perp)
    x_disp_perp = d_i_perp_mag / d_ref_mag

    return (x_disp_path, x_disp_perp)


# def to_kelvin(value: np.ndarray) -> np.ndarray:
#     R"""Convert 3x3 symmetric tensor to Kelvin notation.
#     Parameters
#     ----------
#     value : np.ndarray
#         The 3x3 symmetric tensor.
#
#     Returns
#     -------
#     v_kelvin : np.ndarray
#         The tensor in Kelvin notation,
#         :math:`[\sigma_{xx}, \sigma_{yy}, \sigma_{zz}, \sqrt(2)\sigma_{yz},
#         \sqrt(2)\sigma_{xz}, \sqrt(2)\sigma_{xy}]`.
#     """
#     v_kelvin = np.array(
#         [
#             value[0, 0],
#             value[1, 1],
#             value[2, 2],
#             np.sqrt(2.0) * value[1, 2],
#             np.sqrt(2.0) * value[0, 2],
#             np.sqrt(2.0) * value[0, 1],
#         ]
#     )
#     return v_kelvin
#
#
# def voigt_to_kelvin(v_voigt: np.ndarray) -> np.ndarray:
#     R"""Convert tensor from Voigt notation to Kelvin notation.
#
#     Parameters
#     ----------
#     v_voigt : np.ndarray
#         The tensor in Voigt notation,
#         :math:`[\sigma_{xx}, \sigma_{yy}, \sigma_{zz}, \sigma_{yz}, \sigma_{xz},
#         \sigma_{xy}]`.
#
#     Returns
#     -------
#     v_kelvin : np.ndarray
#         The tensor in Kelvin notation,
#         :math:`[\sigma_{xx}, \sigma_{yy}, \sigma_{zz}, \sqrt(2)\sigma_{yz},
#         \sqrt(2)\sigma_{xz}, \sqrt(2)\sigma_{xy}]`.
#     """
#     v_kelvin = copy.deepcopy(v_voigt)
#     v_kelvin[3:] *= np.sqrt(2.0)
#     return v_kelvin
#
#
# def kelvin_to_voigt(v_kelvin: np.ndarray) -> np.ndarray:
#     R"""Convert tensor from Kelvin notation to Voigt notation.
#
#     Parameters
#     ----------
#     v_kelvin : np.ndarray
#         The tensor in Kelvin notation,
#         :math:`[\sigma_{xx}, \sigma_{yy}, \sigma_{zz}, \sqrt(2)\sigma_{yz},
#         \sqrt(2)\sigma_{xz}, \sqrt(2)\sigma_{xy}]`.
#
#     Returns
#     -------
#     v_voigt : np.ndarray
#         The tensor in Voigt notation,
#         :math:`[\sigma_{xx}, \sigma_{yy}, \sigma_{zz}, \sigma_{yz}, \sigma_{xz},
#         \sigma_{xy}]`.
#     """
#     v_voigt = copy.deepcopy(v_kelvin)
#     v_voigt[3:] /= np.sqrt(2.0)
#     return v_voigt


class NEBImage:
    """Image data for NEB calculations."""

    def __init__(
        self,
        structure: xtal.Structure,
        Ustrain: np.ndarray,
        stress: np.ndarray,
    ):
        R"""

        .. rubric:: Constructor

        Parameters
        ----------
        structure : libcasm.xtal.Structure
            The structure of the image. Expected to include the calculated atom
            properties "energy" and "force".

        Ustrain : np.ndarray
            The lattice strain,
            :math:`[\epsilon_{xx}, \epsilon_{yy}, \epsilon_{zz},
            \sqrt(2)}\epsilon_{yz}, \sqrt(2)}\epsilon_{xz}, \sqrt(2)}\epsilon_{xy}]`,
            relative to the parent prim.

        stress : np.ndarray
            The calculated stress,
            :math:`[\sigma_{xx}, \sigma_{yy}, \sigma_{zz}, \sqrt(2)\sigma_{yz},
            \sqrt(2)\sigma_{xz}, \sqrt(2)\sigma_{xy}]`.

        """
        self.structure = structure
        """libcasm.xtal.Structure: The structure of the image. Expected to include
        the properties "energy", "force", and "Ustrain" once calculated."""

        self.Ustrain = Ustrain
        R"""np.ndarray: The lattice strain,
        :math:`[\epsilon_{xx}, \epsilon_{yy}, \epsilon_{zz}, \sqrt(2)}\epsilon_{yz}, 
        \sqrt(2)}\epsilon_{xz}, \sqrt(2)}\epsilon_{xy}]`, relative to the parent 
        prim."""

        self.stress = stress
        R"""np.ndarray: The calculated stress, 
        :math:`[\sigma_{xx}, \sigma_{yy}, \sigma_{zz}, \sqrt(2)\sigma_{yz}, 
        \sqrt(2)\sigma_{xz}, \sqrt(2)\sigma_{xy}]`."""

    @property
    def energy(self) -> float:
        """float: The energy of the image structure."""
        return self.structure.atom_properties()["energy"][0, 0]

    @property
    def force(self) -> np.ndarray:
        """np.ndarray: The forces on the atoms in the image structure,
        shape (3, n_atoms)."""
        return self.structure.atom_properties()["force"]

    def __rmul__(self, op: xtal.SymOp):
        """Apply symmetry operation to NEBImage.

        Parameters
        ----------
        op : libcasm.xtal.SymOp
            The symmetry operation to apply. This can also be rigid rotations and
            translations.

        Returns
        -------
        new_image : NEBImage
            The new NEBImage with the symmetry operation applied to the structure and
            properties.
        """
        if not isinstance(op, xtal.SymOp):
            raise ValueError(
                "Error in NEBImage.__rmul__: `op` must be a libcasm.xtal.SymOp."
            )
        new_structure = op * self.structure

        # Same matrix_rep for stress and Ustrain:
        M = op.matrix_rep("Ustrain")
        new_Ustrain = M @ self.Ustrain
        new_stress = M @ self.stress

        return NEBImage(structure=new_structure, Ustrain=new_Ustrain, stress=new_stress)

    def copy(self):
        """Create a deep copy of the NEBImage.

        Returns
        -------
        new_image : NEBImage
            A deep copy of the NEBImage.
        """
        new_structure = self.structure.copy()
        new_Ustrain = copy.deepcopy(self.Ustrain)
        new_stress = copy.deepcopy(self.stress)
        return NEBImage(
            structure=new_structure,
            Ustrain=new_Ustrain,
            stress=new_stress,
        )

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

    def to_dict(self):
        """Convert NEBImage to a dict for serialization.

        Returns
        -------
        dict
            The NEBImage as a dictionary. Contains the following keys:

            - structure: The structure of the image.
            - Ustrain: The lattice strain.
            - stress: The calculated stress.
        """
        return {
            "structure": self.structure.to_dict(),
            "Ustrain": self.Ustrain.tolist(),
            "stress": self.stress.tolist(),
        }

    @staticmethod
    def from_dict(data: dict):
        """Create NEBImage from a dict.

        Parameters
        ----------
        data : dict
            The dictionary containing the NEBImage data. Expected to have the following
            keys:

            - structure: The structure of the image.
            - Ustrain: The lattice strain.
            - stress: The calculated stress.

        Returns
        -------
        image : NEBImage
            The NEBImage created from the dictionary.
        """
        structure = xtal.Structure.from_dict(data["structure"])
        Ustrain = np.array(data["Ustrain"])
        stress = np.array(data["stress"])
        return NEBImage(structure=structure, Ustrain=Ustrain, stress=stress)

    def __repr__(self):
        return xtal.pretty_json(self.to_dict())


class NEBPath:
    """Data for working with nudged elastic band (NEB) calculations."""

    def __init__(
        self,
        parent_prim: casmconfig.Prim,
        child: xtal.Structure,
        structure_mapping: mapinfo.StructureMapping,
    ):
        """

        .. rubric:: Constructor

        Parameters
        ----------
        parent_prim: libcasm.configuration.Prim
            The parent prim.
        child: libcasm.xtal.Structure
            The unmapped child structure.
        structure_mapping: libcasm.mapping.info.StructureMapping
            The structure mapping from parent to child, defining the initial /
            directly interpolated NEB path.
        """
        self.parent_prim = parent_prim
        """libcasm.configuration.Prim: The parent prim."""

        self.child = child
        """libcasm.xtal.Structure: The unmapped child structure."""

        self.structure_mapping = structure_mapping
        """libcasm.mapping.info.StructureMapping: The structure mapping from parent to 
        child."""

        self.symmetry_adapated_Hstrain_converter = xtal.StrainConverter(
            metric="Hstrain",
            basis=xtal.make_symmetry_adapted_strain_basis(),
        )
        """libcasm.xtal.StrainConverter: Use to convert deformation gradients to
        symmetry adapted Hencky strain."""

        self._T_parent = None
        self._T_child = None
        self._superchild = None

        initial_structure = mapmethods.make_mapped_structure(
            structure_mapping=self.structure_mapping.interpolated(0.0),
            unmapped_structure=self.superchild,
        )
        self.initial_structure = initial_structure
        """libcasm.xtal.Structure: The structure at the start of the NEB path 
        (the parent)."""

        final_structure = mapmethods.make_mapped_structure(
            structure_mapping=self.structure_mapping.interpolated(1.0),
            unmapped_structure=self.superchild,
        )
        self.final_structure = final_structure
        """libcasm.xtal.Structure: The structure at the end of the NEB path 
        (the mapped child)."""

        lmap_ref, amap_ref = mapmethods.direct_structure_mapping(
            structure1=self.initial_structure,
            structure2=self.final_structure,
            remove_mean_displacement=False,
        )
        self.lattice_mapping_ref = lmap_ref
        """libcasm.mapping.info.LatticeMapping: The reference lattice mapping
        defining the NEB path direction.
        """

        self.atom_mapping_ref = amap_ref
        """libcasm.mapping.info.AtomMapping: The reference atom mapping
        defining the NEB path direction.
        """

    @property
    def T_parent(self) -> np.ndarray:
        """np.ndarray: The transformation matrix from parent prim to the supercell
        used in the NEB calculations.
        """
        if self._T_parent is None:
            lmap = self.structure_mapping.lattice_mapping()
            self._T_parent = np.round(lmap.transformation_matrix_to_super()).astype(int)
        return self._T_parent

    @property
    def T_child(self) -> np.ndarray:
        """np.ndarray: The transformation matrix from child structure to the supercell
        used in the NEB calculations.
        """
        if self._T_child is None:
            from casm.tools.map.methods import (
                make_child_transformation_matrix_to_super,
            )

            _T_child = make_child_transformation_matrix_to_super(
                parent_lattice=self.parent_prim.xtal_prim.lattice(),
                child_lattice=self.child.lattice(),
                structure_mapping=self.structure_mapping,
            )
            self._T_child = np.round(_T_child).astype(int)
        return self._T_child

    @property
    def superchild(self) -> xtal.Structure:
        """libcasm.xtal.Structure: The supercell of the child structure that is
        mapped to parent prim."""
        if self._superchild is None:
            self._superchild = xtal.make_superstructure(
                transformation_matrix_to_super=self.T_child,
                structure=self.child,
            )
        return self._superchild

    def make_path(
        self,
        n_images: typing.Optional[int] = None,
        f_chain: typing.Optional[list[float]] = None,
    ):
        """Make a list of interpolated structures.

        Parameters
        ----------
        n_images : Optional[int] = None
            The number of equally spaced images to generate along the NEB path,
            including endpoints.
        f_chain : Optional[list[float]] = None
            A list of fractional positions along the NEB path at which to generate
            images. If provided, `n_images` is ignored.

        Returns
        -------
        path : list[libcasm.xtal.Structure]
            A list of interpolated structures between the initial and final structures.
        """
        if f_chain is None:
            if n_images is None:
                raise ValueError(
                    "If f_chain is not provided, n_images must be specified."
                )
            f_chain = np.linspace(0.0, 1.0, n_images).tolist()

        path = []
        for f in f_chain:
            structure_mapping_f = self.structure_mapping.interpolated(f)
            structure_f = mapmethods.make_mapped_structure(
                structure_mapping=structure_mapping_f,
                unmapped_structure=self._superchild,
            )
            path.append(structure_f)
        return path

    def calculate(
        self,
        structure: list[xtal.Structure],
        calculator: typing.Any,
        step_index: int,
        image_index: int,
    ) -> NEBImage:
        R"""Calculate forces for a list of structures.

        Parameters
        ----------
        structure : libcasm.xtal.Structure
            The structure to calculate energy, force, and stress for.
        calculator : libcasm.xtal.Calculator
            The calculator to use for force calculations. Calculators should be
            set to use the current working directory (".") for any file I/O. The
            working directory will be set to
            "calculations/step_{step_index}/image_{i_index}" for each calculation.
        step_index : int
            The current NEB step number, used for organizing calculation directories.
        image_index : int
            The image index, used for organizing calculation directories.

        Returns
        -------
        calculated_structure : libcasm.xtal.Structure
            The structure with calculated properties (energy, forces, stress).
        stress : np.ndarray
            The calculated stress,
            :math:`[\sigma_{xx}, \sigma_{yy}, \sigma_{zz}, \sqrt(2)\sigma_{yz},
            \sqrt(2)\sigma_{xz}, \sqrt(2)\sigma_{xy}]`.
        """
        import pathlib

        from casm.tools.shared.ase_utils import make_ase_atoms, make_casm_structure
        from casm.tools.shared.contexts import working_dir

        atoms = make_ase_atoms(casm_structure=structure)
        atoms.calc = calculator

        path = (
            pathlib.Path("calculations") / f"step_{step_index}" / f"image_{image_index}"
        )
        path.mkdir(parents=True, exist_ok=True)
        with working_dir(wd=path):
            energy = atoms.get_potential_energy()  # noqa: F841
            forces = atoms.get_forces()  # noqa: F841
            stress = voigt_to_kelvin(atoms.get_stress(voigt=True))
        calculated_structure = make_casm_structure(ase_atoms=atoms)

        lmap = mapmethods.map_lattices_without_reorientation(
            lattice1=self.parent_prim.xtal_prim.lattice(),
            lattice2=structure.lattice(),
            transformation_matrix_to_super=self.T_parent,
        )
        Ustrain = to_kelvin(lmap.right_stretch())

        return NEBImage(structure=calculated_structure, Ustrain=Ustrain, stress=stress)

    def path_coordinates(self, structure: xtal.Structure):
        """Get lattice and displacement coordinates for a given structure parallel and
        perpendicular to the direct NEB path.

        Parameters
        ----------
        structure : libcasm.xtal.Structure
            The structure to calculate lattice coordinates for.

        Returns
        -------
        x_lattice_path : float
            The lattice strain coordinate along the path direction.
        x_lattice_perp : float
            The lattice strain coordinate perpendicular to the path direction.
        x_disp_path : float
            The displacement coordinate along the path direction.
        x_disp_perp : float
            The displacement coordinate perpendicular to the path direction.
        """
        lmap_i, amap_i = mapmethods.direct_structure_mapping(
            structure1=self.initial_structure,
            structure2=structure,
            remove_mean_displacement=False,
        )
        x_lattice_path, x_lattice_perp = make_lattice_coordinates(
            lmap_ref=self.lattice_mapping_ref,
            lmap_i=lmap_i,
            strain_converter=self.symmetry_adapated_Hstrain_converter,
        )
        x_disp_path, x_disp_perp = make_displacement_coordinates(
            amap_ref=self.atom_mapping_ref,
            amap_i=amap_i,
        )
        return (x_lattice_path, x_lattice_perp, x_disp_path, x_disp_perp)

    def make_chain_orbit(self, chain_prototype: list[NEBImage]):
        """Make an orbit of NEBImages chains from a prototype.

        Parameters
        ----------
        chain_prototype : list[NEBImage]
            The list of NEBImages along the NEB path.

        Returns
        -------
        chain_orbit : list[list[NEBImage]]
            The orbit of inequivalent NEBImage chains generated by applying the
            symmetry operations of the parent prim's factor group to the prototype
            chain.
        """
        from .methods import chain_is_in_orbit

        image_chain_orbit = []
        structure_chain_orbit = []
        for op in self.parent_prim.factor_group.elements:
            image_chain = [op * image for image in chain_prototype]
            structure_chain = [img.structure for img in image_chain]
            if not chain_is_in_orbit(
                chain=structure_chain,
                orbit=structure_chain_orbit,
            ):
                image_chain_orbit.append(image_chain)
                structure_chain_orbit.append(structure_chain)
        return image_chain_orbit


class NEBPathData:
    """Data class for NEB path analysis."""

    def __init__(
        self,
        parent_prim: casmconfig.Prim,
        child: xtal.Structure,
        scored_structure_mapping: mapinfo.ScoredStructureMapping,
        relaxed_structures: list[xtal.Structure],
        energies: list[float],
        strain_converter: typing.Optional[xtal.StrainConverter] = None,
    ):
        """

        .. rubric:: Constructor

        Parameters
        ----------
        parent_prim: libcasm.configuration.Prim
            The primitive cell of the parent structure.
        child: libcasm.xtal.Structure
            The child structure.
        scored_structure_mapping: libcasm.mapping.info.ScoredStructureMapping
            The linear mapping from parent to child structure.
        relaxed_structures: list[libcasm.xtal.Structure]
            The list of relaxed structures along the NEB path.
        energies: list[float]
            The list of energies corresponding to each relaxed structure.
        strain_converter: Optional[libcasm.xtal.StrainConverter] = None
            The strain converter to use for converting deformation gradients to strain.
            By default, the symmetry adapated basis with the Hencky strain metric
            ("Hstrain") is used.
        """
        self.parent_prim = parent_prim
        """libcasm.configuration.Prim: The primitive cell of the parent structure."""

        self.child = child
        """libcasm.xtal.Structure: The child structure."""

        self.scored_structure_mapping = scored_structure_mapping
        """libcasm.mapping.info.ScoredStructureMapping: The linear mapping from
        parent to child structure."""

        self.relaxed_structures = relaxed_structures
        """list[libcasm.xtal.Structure]: The list of relaxed structures along the
        NEB path."""

        self.energies = energies
        """list[float]: The list of energies corresponding to each relaxed structure."""

        if strain_converter is None:
            from math import sqrt

            basis = np.array(
                [
                    [1.0 / sqrt(3), 1.0 / sqrt(3), 1.0 / sqrt(3), 0.0, 0.0, 0.0],
                    [1.0 / sqrt(2), -1.0 / sqrt(2), 0.0, 0.0, 0.0, 0.0],
                    [-1.0 / sqrt(6), -1.0 / sqrt(6), 2.0 / sqrt(6), 0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
                ]
            ).transpose()
            strain_converter = xtal.StrainConverter(
                metric="Hstrain",
                basis=basis,
            )

        self.strain_converter = strain_converter
        """libcasm.xtal.StrainConverter: The strain converter used for converting
        deformation gradients to strain."""

        _x_lattice_path = []
        _x_lattice_perp = []
        _x_disp_path = []
        _x_disp_perp = []

        # Create reference mappings for aligned endpoints:
        lattice_mapping_ref, atom_mapping_ref = mapmethods.direct_structure_mapping(
            structure1=self.relaxed_structures[0],
            structure2=self.relaxed_structures[-1],
            remove_mean_displacement=False,
        )
        # lmap_endpoints = scored_structure_mapping.lattice_mapping()
        # F_path = lmap_endpoints.deformation_gradient()
        # Q_path, U_path = xtal.StrainConverter.F_to_QU(F_path)
        # print("Q_path:\n", Q_path)
        # print("U_path:\n", U_path)

        _e = []
        _d = []
        for i, structure in enumerate(self.relaxed_structures):
            # print(f"NEBImage {i}: Epot = {self.energies[i]:.5f} eV")
            # print(structure)
            # x = structure.lattice().lengths_and_angles()
            # print(f"Lengths and angles: {', '.join([f"{v:.5f}" for v in x])}")
            # print(f"Volume: {structure.lattice().volume():.5f} Å³")

            lmap_i, amap_i = mapmethods.direct_structure_mapping(
                structure1=self.relaxed_structures[0],
                structure2=structure,
                remove_mean_displacement=False,
            )
            F = lmap_i.deformation_gradient()
            _e.append(self.strain_converter.from_F(F))
            x_lattice_path, x_lattice_perp = make_lattice_coordinates(
                lmap_ref=lattice_mapping_ref,
                lmap_i=lmap_i,
                strain_converter=self.strain_converter,
            )
            _x_lattice_path.append(x_lattice_path)
            _x_lattice_perp.append(x_lattice_perp)
            # print(
            #     f"  ||x_lattice_path|| = {x_lattice_path:.5f}, "
            #     f"||x_lattice_perp|| = {x_lattice_perp:.5f}"
            # )

            _d.append(amap_i.displacement())
            x_disp_path, x_disp_perp = make_displacement_coordinates(
                amap_ref=atom_mapping_ref,
                amap_i=amap_i,
            )
            _x_disp_path.append(x_disp_path)
            _x_disp_perp.append(x_disp_perp)
            # print(
            #     f"  ||x_disp_path|| = {x_disp_path:.5f}, "
            #     f"||x_disp_perp|| = {x_disp_perp:.5f}"
            # )
            #
            # print()

        self.x_lattice_path = _x_lattice_path
        """list[float]: The projection of the lattice strain vector onto the 
        parent-to-child strain vector."""

        self.x_lattice_perp = _x_lattice_perp
        """list[float]: The projection of the lattice strain vector perpendicular to the
        parent-to-child strain vector."""

        self.x_disp_path = _x_disp_path
        """list[float]: The projection of the atomic displacement vector onto the
        parent-to-child displacement vector, in the parent reference state."""

        self.x_disp_perp = _x_disp_perp
        """list[float]: The projection of the atomic displacement vector perpendicular 
        to the parent-to-child displacement vector, in the parent reference state."""

        self.e = np.array(_e).transpose()
        """np.array: Strain components for each image, shape (6, n_images)"""

        self.d = np.array(_d).transpose((1, 2, 0))
        """np.array: Displacement vectors for each image, 
        shape (3, n_atoms, n_images)"""

    def to_dict(self):
        """Convert NEBPathData to a dict for serialization.

        Returns
        -------
        dict
            The NEBPathData as a dictionary. Contains the following keys:

            - parent_prim: The primitive cell of the parent structure.
            - child: The child structure.
            - scored_structure_mapping: The linear mapping from parent to child
              structure.
            - relaxed_structures: The list of relaxed structures along the NEB path.
            - energies: The list of energies corresponding to each relaxed structure.
            - strain_metric: The strain metric used in the strain converter.
            - strain_basis: The strain basis used in the strain converter, where
              `strain_basis[i]` is the i-th basis vector.
            - x_lattice_path: The lattice path coordinates for each image.
            - x_lattice_perp: The lattice perpendicular coordinates for each image.
            - x_disp_path: The displacement path coordinates for each image.
            - x_disp_perp: The displacement perpendicular coordinates for each image.

        """
        data = {
            "parent_prim": self.parent_prim.to_dict(),
            "child": self.child.to_dict(),
            "scored_structure_mapping": self.scored_structure_mapping.to_dict(),
            "relaxed_structures": [s.to_dict() for s in self.relaxed_structures],
            "energies": self.energies,
            "strain_metric": self.strain_converter.metric(),
            "strain_basis": self.strain_converter.basis().transpose().tolist(),
            "x_lattice_path": self.x_lattice_path,
            "x_lattice_perp": self.x_lattice_perp,
            "x_disp_path": self.x_disp_path,
            "x_disp_perp": self.x_disp_perp,
        }
        return data

    @staticmethod
    def from_dict(data: dict):
        """Create NEBPathData from a dict.

        Parameters
        ----------
        data : dict
            The NEBPathData as a dictionary, as serialized by `to_dict()`.

        Returns
        -------
        path: NEBPathData
            The NEBPathData object.
        """
        parent_prim = casmconfig.Prim.from_dict(data["parent_prim"])
        child = xtal.Structure.from_dict(data["child"])
        scored_structure_mapping = mapinfo.ScoredStructureMapping.from_dict(
            data["scored_structure_mapping"],
            prim=parent_prim,
        )
        relaxed_structures = [
            xtal.Structure.from_dict(s) for s in data["relaxed_structures"]
        ]
        energies = data["energies"]
        strain_metric = data["strain_metric"]
        strain_basis = np.array(data["strain_basis"]).transpose()
        strain_converter = xtal.StrainConverter(
            metric=strain_metric,
            basis=strain_basis,
        )

        neb_path_data = NEBPathData(
            parent_prim=parent_prim,
            child=child,
            scored_structure_mapping=scored_structure_mapping,
            relaxed_structures=relaxed_structures,
            energies=energies,
            strain_converter=strain_converter,
        )
        return neb_path_data

    def plot_energy_vs_lattice_path(self, title: str):
        """Use bokeh to plot energy vs lattice path coordinate"""
        from bokeh.plotting import figure, show

        p = figure(
            title=title,
            x_axis_label="Lattice Path Coordinate",
            y_axis_label="Energy (eV)",
        )

        p.line(
            self.x_lattice_path,
            self.energies,
            line_width=2,
            legend_label="Energy vs Lattice Path",
        )

        # Include points:
        p.circle(
            self.x_lattice_path,
            self.energies,
            size=8,
            color="red",
            legend_label="Images",
        )

        show(p)  # This will open a browser window with the plot

    def plot_energy_vs_disp_path(self, title: str):
        """Use bokeh to plot energy vs displacement path coordinate"""
        from bokeh.plotting import figure, show

        p = figure(
            title=title,
            x_axis_label="Displacement Path Coordinate",
            y_axis_label="Energy (eV)",
        )

        p.line(
            self.x_disp_path,
            self.energies,
            line_width=2,
            legend_label="Energy vs Displacement Path",
        )

        # Include points:
        p.circle(
            self.x_disp_path,
            self.energies,
            size=8,
            color="red",
            legend_label="Images",
        )

        show(p)  # This will open a browser window with the plot

    def plot_path_coordinates(self, title: str):
        """2D plot of lattice and displacement path coordinates"""
        from bokeh.plotting import figure, show

        p = figure(
            title=title,
            x_axis_label="Lattice Path Coordinate",
            y_axis_label="Displacement Path Coordinate",
        )

        p.line(
            self.x_lattice_path,
            self.x_disp_path,
            line_width=2,
            legend_label="Path in Coordinate Space",
        )

        # Include points:
        p.circle(
            self.x_lattice_path,
            self.x_disp_path,
            size=8,
            color="red",
            legend_label="Images",
        )

        show(p)  # This will open a browser window with the plot

    def super_plot(self, title: str, width: int = 600, height: int = 800):
        """Plot all three plots as subfigures in one window"""
        from bokeh.layouts import column
        from bokeh.models import Div
        from bokeh.plotting import figure, show

        p1 = figure(
            title="Energy vs Lattice Path Coordinate",
            x_axis_label="Lattice Path Coordinate",
            y_axis_label="Energy (eV)",
            width=width,
            height=int(round(height / 4)),
        )
        p1.line(
            self.x_lattice_path,
            self.energies,
            line_width=2,
            # legend_label="Energy vs Lattice Path",
        )
        p1.circle(
            self.x_lattice_path,
            self.energies,
            size=8,
            color="red",
            # legend_label="Images",
        )

        p2 = figure(
            title="Energy vs Displacement Path Coordinate",
            x_axis_label="Displacement Path Coordinate",
            y_axis_label="Energy (eV)",
            width=width,
            height=int(round(height / 4)),
        )
        p2.line(
            self.x_disp_path,
            self.energies,
            line_width=2,
            # legend_label="Energy vs Displacement Path",
        )
        p2.circle(
            self.x_disp_path,
            self.energies,
            size=8,
            color="red",
            # legend_label="Images",
        )

        p3 = figure(
            title="Lattice and Displacement Path Coordinates",
            x_axis_label="Lattice Path Coordinate",
            y_axis_label="Displacement Path Coordinate",
            width=width,
            height=int(round(height / 2)),
        )
        p3.line(
            self.x_lattice_path,
            self.x_disp_path,
            line_width=2,
            # legend_label="Path in Coordinate Space",
        )
        p3.circle(
            self.x_lattice_path,
            self.x_disp_path,
            size=8,
            color="red",
            # legend_label="Images",
        )

        p4 = figure(
            title="Strain Components vs Lattice Path Coordinate",
            x_axis_label="Lattice Path Coordinate",
            y_axis_label="Strain Components",
            width=width,
            height=int(round(height / 2)),
        )
        colors = ["red", "green", "blue", "magenta", "cyan", "yellow"]
        for i in range(self.e.shape[0]):
            p4.line(
                self.x_lattice_path,
                self.e[i, :],
                line_width=2,
                color=colors[i],
                legend_label=f"e_{i+1}",
            )

        layout = column(
            Div(text=f"<h1>{title}</h1>"),
            p1,
            p2,
            p3,
            p4,
        )
        show(layout)  # This will open a browser window with the plots
