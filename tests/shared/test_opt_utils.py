import math

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

import casm.tools.shared.opt_utils as opt_utils
import libcasm.mapping.methods as mapmethods
import libcasm.xtal as xtal
import libcasm.xtal.lattices as xtal_lattices
import libcasm.xtal.structures as xtal_structures


def hcp_reference_structure():
    return xtal_structures.HCP(r=1.6, atom_type="Mg")


def hcp_structure_1():
    # strained only
    return xtal.Structure(
        lattice=xtal_lattices.HCP(r=1.8),
        atom_coordinate_frac=np.array(
            [
                [1 / 3.0, 2 / 3, 1 / 4],
                [2 / 3.0, 1 / 3, 3 / 4],
            ]
        ).transpose(),
        atom_type=["Mg", "Mg"],
    )


def hcp_structure_2():
    # strained, and displaced
    return xtal.Structure(
        lattice=xtal_lattices.HCP(r=1.8),
        atom_coordinate_frac=np.array(
            [
                [1 / 3.0, 2 / 3, 1 / 4 + 0.1],
                [2 / 3.0, 1 / 3, 3 / 4 - 0.1],
            ]
        ).transpose(),
        atom_type=["Mg", "Mg"],
    )


def hcp_structure_2_rotated():
    hcp_lattice = xtal_lattices.HCP(r=1.8)
    L0 = hcp_lattice.column_vector_matrix()

    # rotate L0
    rot = Rotation.from_euler("XYZ", [20, 30, 60], degrees=True)
    Q0 = rot.as_matrix()

    L = Q0 @ L0
    lattice = xtal.Lattice(column_vector_matrix=L)

    # strained, and displaced
    return xtal.Structure(
        lattice=lattice,
        atom_coordinate_frac=np.array(
            [
                [1 / 3.0, 2 / 3, 1 / 4 + 0.1],
                [2 / 3.0, 1 / 3, 3 / 4 - 0.1],
            ]
        ).transpose(),
        atom_type=["Mg", "Mg"],
    )


def hcp_structure_3():
    # displaced only
    return xtal.Structure(
        lattice=xtal_lattices.HCP(r=1.6),
        atom_coordinate_frac=np.array(
            [
                [1 / 3.0, 2 / 3, 1 / 4 + 0.1],
                [2 / 3.0, 1 / 3, 3 / 4 - 0.1],
            ]
        ).transpose(),
        atom_type=["Mg", "Mg"],
    )


def make_test_data(casm_structure):
    from ase.calculators.kim.kim import KIM

    import casm.tools.shared.ase_utils as ase_utils
    from casm.tools.shared.conversions import voigt_to_kelvin

    calculator = KIM("MEAM_LAMMPS_DickelBaskesAslam_2018_MgAlZn__MO_093637366498_002")

    ase_atoms = ase_utils.make_ase_atoms(casm_structure)
    ase_atoms.calc = calculator

    force = ase_atoms.get_forces().transpose()
    print("FORCE=\n", force.tolist())

    stress = voigt_to_kelvin(ase_atoms.get_stress(voigt=True))
    print("STRESS=\n", stress.tolist())


def hcp_structure_1_force_stress():
    f1z = -5.320824930905536
    f2z = -f1z
    force = np.array(
        [
            [0.0, 0.0],
            [0.0, 0.0],
            [f1z, f2z],
        ]
    )
    s11 = -0.08355538386475626
    s22 = s11
    s33 = -0.14039149934916226
    stress = np.array([s11, s22, s33, 0.0, 0.0, 0.0])
    return force, stress


def hcp_structure_1_disponly_expected_grad():
    f1z = -5.320824930905536
    f2z = -f1z
    return np.array([0.0, 0.0, -f1z, 0.0, 0.0, -f2z])


def hcp_structure_2_force_stress():
    f1z = -2.1854905532529636
    f2z = -f1z
    force = np.array(
        [
            [0.0, 0.0],
            [0.0, 0.0],
            [f1z, f2z],
        ]
    )
    s11 = 5.283641649499512e-05
    s22 = s11
    s33 = -0.03924824568659316
    stress = np.array([s11, s22, s33, 0.0, 0.0, 0.0])
    return force, stress


def hcp_structure_2_standard_expected_grad():
    return np.array(
        [
            0.0,
            0.0,
            1.9426582695581898,
            0.0,
            0.0,
            -1.9426582695581898,
            0.004412258007733448,
            0.004412258007733448,
            -3.2775384442766686,
            0.0,
            0.0,
            0.0,
        ]
    )


def hcp_structure_2_symmetry_adapted_expected_grad():
    return np.array(
        [
            0.0,
            0.0,
            1.9426582695581898,
            0.0,
            0.0,
            -1.9426582695581898,
            -1.887192866385484,
            0.0,
            -2.6797015271885627,
            0.0,
            0.0,
            0.0,
        ]
    )


@pytest.mark.xfail(reason="To be determined...")
def test_StrainDispVarTool_standard_rotated():
    reference_structure = hcp_reference_structure()
    structure = hcp_structure_2_rotated()

    tool = opt_utils.StrainDispVarTool(
        reference_structure=reference_structure,
        include_displacements=True,
        include_strain=True,
        strain_metric="Hstrain",
    )
    assert isinstance(tool, opt_utils.StrainDispVarTool)

    x = tool.make_x(structure=structure)
    assert isinstance(x, np.ndarray)
    assert x.shape == (12,)

    expected_x = np.zeros((12,))

    # expected displacement components:
    # F * (x1 +d) = x2
    # x2 = F * L1 * (x_frac_1 + d_frac)
    # => d = L1 * d_frac,
    #
    # d_frac = [0.0, 0.0, +/- 0.1].T

    L1 = reference_structure.lattice().column_vector_matrix()
    c = L1[2, 2]
    expected_x[0:6] = np.array([0.0, 0.0, 0.1 * c, 0.0, 0.0, -0.1 * c])

    # expected strain components:
    # Hstrain = log(F.T * F) / 2
    F11 = 1.8 / 1.6
    expected_x[6:9] = math.log(F11**2) / 2
    assert np.allclose(x, expected_x)

    # check re-created structure
    test_structure = tool.make_structure(x=x)
    test_structure_within = xtal.make_structure_within(test_structure)
    test_fg = xtal.make_structure_factor_group(test_structure_within)

    # print("STRUCTURE=\n", structure)
    # fg = xtal.make_structure_factor_group(structure)
    # print("# FG ELEMENTS=", len(fg))
    # print("TEST STRUCTURE (WITHIN)=\n", test_structure_within)
    # print("# FG ELEMENTS=", len(test_fg))

    assert len(test_fg) == 12
    assert test_structure.is_equivalent_to(structure) is False

    lmap_ref, amap_ref = mapmethods.direct_structure_mapping(
        structure1=structure,
        structure2=test_structure,
        remove_mean_displacement=False,
    )

    # assert that right_stretch and displacement are approximately zero:
    right_streth = lmap_ref.right_stretch()
    assert np.allclose(right_streth, np.eye(3))
    disp = amap_ref.displacement()
    assert np.allclose(disp, np.zeros((3, 2)))

    # assert the rotation is as expected
    rot = Rotation.from_euler("XYZ", [20, 30, 60], degrees=True)
    Q0 = rot.as_matrix()
    assert np.allclose(Q0, lmap_ref.isometry().T)

    # check gradient
    force, stress = hcp_structure_2_force_stress()
    grad = tool.make_grad(x=x, force=force, stress=stress)
    expected_grad = hcp_structure_2_standard_expected_grad()
    assert np.allclose(grad, expected_grad)


@pytest.mark.xfail(reason="To be determined...")
def test_StrainDispVarTool_standard():
    reference_structure = hcp_reference_structure()
    structure = hcp_structure_2()

    tool = opt_utils.StrainDispVarTool(
        reference_structure=reference_structure,
        include_displacements=True,
        include_strain=True,
        strain_metric="Hstrain",
    )
    assert isinstance(tool, opt_utils.StrainDispVarTool)

    x = tool.make_x(structure=structure)
    assert isinstance(x, np.ndarray)
    assert x.shape == (12,)

    expected_x = np.zeros((12,))

    # expected displacement components:
    # F * (x1 +d) = x2
    # x2 = F * L1 * (x_frac_1 + d_frac)
    # => d = L1 * d_frac,
    #
    # d_frac = [0.0, 0.0, +/- 0.1].T

    L1 = reference_structure.lattice().column_vector_matrix()
    c = L1[2, 2]
    expected_x[0:6] = np.array([0.0, 0.0, 0.1 * c, 0.0, 0.0, -0.1 * c])

    # expected strain components:
    # Hstrain = log(F.T * F) / 2
    F11 = 1.8 / 1.6
    expected_x[6:9] = math.log(F11**2) / 2
    assert np.allclose(x, expected_x)

    # check re-created structure
    test_structure = tool.make_structure(x=x)
    assert test_structure.is_equivalent_to(structure)

    # check gradient
    force, stress = hcp_structure_2_force_stress()
    grad = tool.make_grad(x=x, force=force, stress=stress)
    expected_grad = hcp_structure_2_standard_expected_grad()
    assert np.allclose(grad, expected_grad)


@pytest.mark.xfail(reason="To be determined...")
def test_StrainDispVarTool_Ustrain():
    reference_structure = hcp_reference_structure()
    structure = hcp_structure_2()

    tool = opt_utils.StrainDispVarTool(
        reference_structure=reference_structure,
        include_displacements=True,
        include_strain=True,
        strain_metric="Ustrain",
    )
    assert isinstance(tool, opt_utils.StrainDispVarTool)

    x = tool.make_x(structure=structure)
    assert isinstance(x, np.ndarray)
    assert x.shape == (12,)

    # expected_x = np.zeros((12,))
    #
    # # expected displacement components:
    # # F * (x1 +d) = x2
    # # x2 = F * L1 * (x_frac_1 + d_frac)
    # # => d = L1 * d_frac,
    # #
    # # d_frac = [0.0, 0.0, +/- 0.1].T
    #
    # L1 = reference_structure.lattice().column_vector_matrix()
    # c = L1[2, 2]
    # expected_x[0:6] = np.array([0.0, 0.0, 0.1 * c, 0.0, 0.0, -0.1 * c])
    #
    # # expected strain components:
    # # Hstrain = log(F.T * F) / 2
    # F11 = 1.8 / 1.6
    # expected_x[6:9] = math.log(F11**2) / 2
    # assert np.allclose(x, expected_x)

    # check re-created structure
    test_structure = tool.make_structure(x=x)
    assert test_structure.is_equivalent_to(structure)

    # check gradient
    force, stress = hcp_structure_2_force_stress()
    grad = tool.make_grad(x=x, force=force, stress=stress)
    expected_grad = hcp_structure_2_standard_expected_grad()
    assert np.allclose(grad, expected_grad)


@pytest.mark.xfail(reason="To be determined...")
def test_StrainDispVarTool_symmetry_adapted():
    reference_structure = hcp_reference_structure()
    structure = hcp_structure_2()

    tool = opt_utils.StrainDispVarTool(
        reference_structure=reference_structure,
        include_displacements=True,
        include_strain=True,
        strain_metric="Hstrain",
        strain_basis="symmetry_adapted",
    )
    assert isinstance(tool, opt_utils.StrainDispVarTool)

    x = tool.make_x(structure=structure)
    assert isinstance(x, np.ndarray)
    assert x.shape == (12,)

    expected_x = np.zeros((12,))

    # expected displacement components:
    L1 = reference_structure.lattice().column_vector_matrix()
    c = L1[2, 2]
    expected_x[0:6] = np.array([0.0, 0.0, 0.1 * c, 0.0, 0.0, -0.1 * c])

    # expected strain components (symmetry adapted basis):
    F11 = 1.8 / 1.6
    expected_x[6] = (math.log((F11) ** 2) / 2) * math.sqrt(3)
    assert np.allclose(x, expected_x)

    # check re-created structure
    test_structure = tool.make_structure(x=x)
    assert test_structure.is_equivalent_to(structure)

    # check gradient
    force, stress = hcp_structure_2_force_stress()
    grad = tool.make_grad(x=x, force=force, stress=stress)
    expected_grad = hcp_structure_2_symmetry_adapted_expected_grad()
    assert np.allclose(grad, expected_grad)


def test_StrainDispVarTool_disponly():
    reference_structure = hcp_reference_structure()
    structure = hcp_structure_3()  # <-- displaced only

    tool = opt_utils.StrainDispVarTool(
        reference_structure=reference_structure,
        include_displacements=True,
        include_strain=False,
    )
    assert isinstance(tool, opt_utils.StrainDispVarTool)

    x = tool.make_x(structure=structure)
    assert isinstance(x, np.ndarray)
    assert x.shape == (6,)

    # expected displacement components:
    L1 = reference_structure.lattice().column_vector_matrix()
    c = L1[2, 2]
    expected_x = np.array([0.0, 0.0, 0.1 * c, 0.0, 0.0, -0.1 * c])
    assert np.allclose(x, expected_x)

    # check re-created structure
    test_structure = tool.make_structure(x=x)
    assert test_structure.is_equivalent_to(structure)

    # check gradient
    force, stress = hcp_structure_1_force_stress()
    grad = tool.make_grad(x=x, force=force, stress=stress)
    expected_grad = hcp_structure_1_disponly_expected_grad()
    assert np.allclose(grad, expected_grad)

    # check using displaced only structure:
    structure = hcp_structure_3()
    x = tool.make_x(structure=structure)
    test_structure = tool.make_structure(x=x)
    assert test_structure.is_equivalent_to(structure)


@pytest.mark.xfail(reason="To be determined...")
def test_StrainDispVarTool_strainonly():
    reference_structure = hcp_reference_structure()
    structure = hcp_structure_1()  # <-- strained only

    tool = opt_utils.StrainDispVarTool(
        reference_structure=reference_structure,
        include_displacements=False,
        include_strain=True,
        strain_metric="Hstrain",
    )
    assert isinstance(tool, opt_utils.StrainDispVarTool)

    x = tool.make_x(structure=structure)
    assert isinstance(x, np.ndarray)
    assert x.shape == (6,)

    # expected strain components:
    expected_x = np.zeros((6,))
    F11 = 1.8 / 1.6
    expected_x[0:3] = math.log(F11**2) / 2
    assert np.allclose(x, expected_x)

    # check re-created structure
    test_structure = tool.make_structure(x=x)
    assert test_structure.is_equivalent_to(structure)

    # check gradient
    force, stress = hcp_structure_2_force_stress()
    grad = tool.make_grad(x=x, force=force, stress=stress)
    expected_grad = hcp_structure_2_standard_expected_grad()[6:]
    assert np.allclose(grad, expected_grad)

    # check using strained only structure:
    structure = hcp_structure_1()
    x = tool.make_x(structure=structure)
    test_structure = tool.make_structure(x=x)
    assert test_structure.is_equivalent_to(structure)


def test_StrainDispVarTool_strainonly_GLstrain():
    reference_structure = hcp_reference_structure()
    structure = hcp_structure_2()  # <-- strained and displaced

    tool = opt_utils.StrainDispVarTool(
        reference_structure=reference_structure,
        include_displacements=False,
        include_strain=True,
        strain_metric="GLstrain",
    )
    assert isinstance(tool, opt_utils.StrainDispVarTool)

    x = tool.make_x(structure=structure)
    assert isinstance(x, np.ndarray)
    assert x.shape == (6,)

    # expected strain components:
    expected_x = np.zeros((6,))
    F11 = 1.8 / 1.6
    expected_x[0:3] = (F11**2 - 1) / 2
    assert np.allclose(x, expected_x)

    # because disp is not included, the test structure should not be equivalent
    test_structure = tool.make_structure(x=x)
    assert not test_structure.is_equivalent_to(structure)

    # check using strained only structure:
    structure = hcp_structure_1()
    x = tool.make_x(structure=structure)
    test_structure = tool.make_structure(x=x)
    assert test_structure.is_equivalent_to(structure)
