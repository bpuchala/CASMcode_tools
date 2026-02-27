import typing

import numpy as np

import libcasm.mapping.methods as mapmethods
import libcasm.xtal as xtal


class StrainDispVarTool:
    R"""Tool for working with strain and displacement variables for optimization.

    Variables are atomic displacements and/or lattice strain, defined relative
    to a reference structure:

    .. math::

        F L_1 = L_2 \\\\
        F \left(\vec{r_1}(i) + \vec{d}(i) \right) = \vec{r_2}(i)

    where:

    - :math:`L_1` is a shape=(3,3) matrix with columns containing the
      reference "parent" lattice vectors.
    - :math:`L_2` is a shape=(3,3) matrix with columns containing the
      "child" lattice vectors.
    - :math:`F` is the parent-to-child deformation gradient tensor,
      a shape=(3,3) matrix.
    - :math:`\vec{r_1}(i)` is the Cartesian coordinates of the i-th atom in
      the reference "parent" structure.
    - :math:`\vec{r_2}(i)` is the Cartesian coordinates of i-th atom in
      the "child" structure.

    Displacement variables, :math:`\vec{d}(i)`, are ordered as:

    .. math::

        [d(1)_x, d(1)_y, d(1)_z, d(2)_x, d(2)_y, d(2)_z, ...]

    Strain variables (Green-Lagrange, Hencky, or Euler-Almansi) are constructed from
    the deformation gradient, :math:`F`, using :class:`libcasm.xtal.StrainConverter`:

    - `"GLstrain"`: Green-Lagrange strain metric,
      :math:`E = \frac{1}{2}(F^{\mathsf{T}} F - I)`
    - `"Hstrain"`: Hencky strain metric, :math:`E = \frac{1}{2}\ln(F^{\mathsf{T}} F)`
    - `"EAstrain"`: Euler-Almansi strain metric,
      :math:`E = \frac{1}{2}(I−(F F^{\mathsf{T}})^{-1})`

    Strain variables may be represented in Kelvin notation using the standard basis:

    .. math::

        [E_{xx}, E_{yy}, E_{zz}, \sqrt{2}E_{yz}, \sqrt{2}E_{xz}, \sqrt{2}E_{xy}]

    or may be transformed into a symmetry-adapted or other user-specified basis,
    see :ref:`Strain DoF <sec-strain-dof>`.

    Notes
    -----

    Currently, this supports optimization of all displacement variables and/or
    all strain variables. Support for constrained optimization of particular modes
    may be added in the future.

    """

    def __init__(
        self,
        reference_structure: xtal.Structure,
        include_strain: bool = True,
        include_displacements: bool = True,
        strain_metric: str = "GLstrain",
        strain_basis: typing.Union[np.ndarray, str, None] = None,
    ):
        """
        .. rubric:: Constructor

        Parameters
        ----------
        reference_structure: xtal.Structure
            The reference structure, used for the lattice reference, number of atoms,
            and reference positions for displacements.
        include_strain: bool = True
            Whether to include strain variables.
        include_displacements: bool = True
            Whether to include displacement variables.
        strain_metric: str = "GLstrain"
            The strain metric to use. Options are "GLstrain", "Hstrain", and "EAstrain".
            See :class:`libcasm.xtal.StrainConverter` for more details. The metrics
            "Ustrain" and "Bstrain" are not supported for stress conversions.
        strain_basis: Optional[np.ndarray] = None
            The basis to use for the strain metric. If None, the identity matrix is
            used. If "symmetry_adapted", then the symmetry adapted basis is used
            (see :func:`libcasm.xtal.make_symmetry_adapted_strain_basis`).
        """
        self.reference_structure = reference_structure
        """xtal.Structure: The reference structure."""

        self.include_displacements = include_displacements
        """bool: Whether to include displacement variable."""

        self.include_strain = include_strain
        """bool: Whether to include strain variables."""

        if strain_basis is None:
            strain_basis = np.eye(6)
        elif strain_basis == "symmetry_adapted":
            strain_basis = xtal.make_symmetry_adapted_strain_basis()
        else:
            if strain_basis.shape[1] != 6 or strain_basis.shape[0] != 6:
                raise ValueError(
                    "strain_basis must be a 6x6 matrix, 'symmetry_adapted', or None"
                )

        self.strain_metric = strain_metric
        """str: The strain metric used for lattice variables."""

        self.strain_basis = strain_basis
        """np.ndarray: The strain basis used for lattice variables."""

        self.strain_basis_pinv = np.linalg.pinv(strain_basis)
        """np.ndarray: The pseudoinverse of the strain basis used for lattice 
        variables."""

        self.strain_converter = xtal.StrainConverter(
            metric=strain_metric,
            basis=strain_basis,
        )
        """xtal.StrainConverter: The strain converter for lattice variables."""

        self.n_atoms = len(self.reference_structure.atom_type())
        """int: The number of atoms in the reference structure."""

        self.n_x_disp = self.n_atoms * 3 if self.include_displacements else 0
        """int: The number of displacement variables."""

        self.n_x_strain = 6 if self.include_strain else 0
        """int: The number of strain variables."""

    def make_x(
        self,
        structure: xtal.Structure,
    ) -> np.ndarray:
        R"""Make a variable array from the given structure.

        Parameters
        ----------
        structure: xtal.Structure
            The structure to convert to variables.

        Returns
        -------
        variables: np.ndarray
            The variable array representing the structure. If displacement variables are
            included, they are unrolled as
            :math:`[d_{1x}, d_{1y}, d_{1z}, d_{2x}, ...]`. If strain variables are
            included, they are in Kelvin notation, :math:`[E_{xx}, E_{yy}, E_{zz},
            \sqrt(2)E_{yz}, \sqrt(2)E_{xz}, \sqrt(2)E_{xy}]`.
        """
        lmap_ref, amap_ref = mapmethods.direct_structure_mapping(
            structure1=self.reference_structure,
            structure2=structure,
            remove_mean_displacement=False,
        )

        x = []

        if self.include_displacements:
            D = amap_ref.displacement()  # (3, N) -> (N, 3) -> (3*N,)

            # [d1x, d1y, d1z, d2x, ...]
            x_disp = D.transpose().flatten()

            x.extend(x_disp.tolist())

        if self.include_strain:
            F = lmap_ref.deformation_gradient()  # (3, 3)
            x_strain = self.strain_converter.from_F(F)  # (6,)

            x.extend(x_strain.tolist())

        x = np.array(x)

        # print("----------------")
        # print("X=", x)
        # print("----------------")

        return x

    def make_structure(
        self,
        x: np.ndarray,
    ) -> xtal.Structure:
        """Make a structure from the given variable array.

        Parameters
        ----------
        x: np.ndarray
            See :func:`make_x` for variable array format.

        Returns
        -------
        structure: xtal.Structure
            The structure represented by the variable array.

        """
        structure = self.reference_structure.copy()

        index = 0

        if self.include_displacements:
            n_disp = self.n_x_disp
            x_disp = x[index : index + n_disp]
            index += n_disp

            # (N, 3)
            D = x_disp.reshape((self.n_atoms, 3)).transpose()

            atom_coordinate_cart = self.reference_structure.atom_coordinate_cart() + D
            atom_coordinate_frac = xtal.cartesian_to_fractional(
                lattice=structure.lattice(),
                coordinate_cart=atom_coordinate_cart,
            )
        else:
            atom_coordinate_frac = self.reference_structure.atom_coordinate_frac()

        if self.include_strain:
            x_strain = x[index : index + self.n_x_strain]
            F = self.strain_converter.to_F(x_strain)

            L_ref = self.reference_structure.lattice().column_vector_matrix()
            L = F @ L_ref
            lattice = xtal.Lattice(
                column_vector_matrix=L,
            )
        else:
            lattice = self.reference_structure.lattice()

        return xtal.Structure(
            lattice=lattice,
            atom_coordinate_frac=atom_coordinate_frac,
            atom_type=self.reference_structure.atom_type(),
        )

    def make_grad(
        self,
        x: np.ndarray,
        force: np.ndarray,
        stress: np.ndarray,
    ):
        """Create the energy gradient vector from calculated force and stress.

        Parameters
        ----------
        x: np.ndarray
            The variable array at which the gradient is evaluated.
        force: np.ndarray
            The forces on the atoms, shape (3, N_atoms).
        stress: np.ndarray
            The stress, in Kelvin notation, shape (6,).

        Returns
        -------
        grad: np.ndarray
            The gradient vector.
        """

        index = 0

        if self.include_displacements:
            index += self.n_x_disp

        if self.include_strain:
            n_strain = self.n_x_strain
            x_strain = x[index : index + n_strain]
            index += n_strain

            F = self.strain_converter.to_F(x_strain)
            F_inv = np.linalg.inv(F)
            L_ref = self.reference_structure.lattice().column_vector_matrix()
            L = F @ L_ref
            lattice = xtal.Lattice(column_vector_matrix=L)
        else:
            F = np.eye(3)
            F_inv = F
            lattice = self.reference_structure.lattice()

        volume_ref = self.reference_structure.lattice().volume()
        volume = lattice.volume()

        # Build the gradient vector
        grad = np.zeros((index,))

        # Forces are negative gradients of energy with respect to positions
        # in the deformed state, but displacements are defined in the
        # undeformed state.

        # grad_disp = dEpot/ddisp

        # force = -dEpot/dx2

        # F * (x1 + disp) = x2
        # F * dEpot/ddisp = dEpot/dx2
        # => dEpot/ddisp = F_inv * dEpot/dx2
        # => f_disp = - F_inv * force

        index = 0
        if self.include_displacements:
            grad_disp = -(F_inv @ force)  # (3, N)

            grad[index : index + self.n_x_disp] = grad_disp.transpose().flatten()
            index += self.n_x_disp

        # Stress is derivative of energy with respect to strain, but
        # depending on strain metric, conversion is needed

        # grad_strain = dEpot/dEstrain_metric

        # stress = (1/volume) * dE/dE_GL

        # For GLstrain:
        # dEpot/dEstrain_GL = volume * stress
        # => grad_strain = volume * stress
        # For Hstrain or EAstrain:
        # dEpot/dEstrain_metric = F @ (dEpot/dEstrain_GL) @ F.transpose()
        # => grad_strain = volume * F @ stress @ F.transpose()

        # If using a non-standard basis:
        # Estrain_standard = B @ Estrain_user
        # dEpot/dEstrain_user = B.pinv() @ dEpot/dEstrain_standard
        # => grad_strain = B.pinv() @ grad_strain_standard

        if self.include_strain:
            from casm.tools.shared.conversions import from_kelvin, to_kelvin

            # stress = 1/V_k * dU_k/dE^GL_kk
            # stress_pk2 = 1/V_1 * dU_k/dE^GL_1k, second Piola-Kirchhoff stress
            #     referenced to reference structure
            # stress_pk2 = (V_k / V_1) * F_inv @ stress @ F_inv.T

            stress_3x3 = from_kelvin(stress)
            stress_pk2_3x3 = (volume / volume_ref) * F_inv @ stress_3x3 @ F_inv.T

            if self.strain_metric == "GLstrain":
                # dU_k/dE^GL_1k = V_1 * stress_pk2
                grad_strain_standard = volume_ref * to_kelvin(stress_pk2_3x3)  # (6,)
            elif self.strain_metric == "Hstrain":
                # stress_tau = (V_k / V_1) * stress, kirchoff stress
                # stress_tau = 1/V_1 * dU_k/dE^H_1k
                # dU_k/dE^H_1k = V_k * stress
                grad_strain_standard = volume * stress  # (6,)

            elif self.strain_metric == "EAstrain":
                # stress = 1/V_1 * dU_k/dE^EA_1k
                # dU_k/dE^EA_1k = V_1 * stress
                grad_strain_standard = volume_ref * stress  # (6,)
            else:
                raise ValueError(
                    "Unsupported strain metric for stress conversion: "
                    f"{self.strain_metric}"
                )

            grad_strain = self.strain_basis_pinv @ grad_strain_standard  # (6,)

            grad[index : index + self.n_x_strain] = grad_strain

        # print("----------------")
        # print("X=", x)
        # print("GRAD=", grad)
        # print("----------------")

        return grad
