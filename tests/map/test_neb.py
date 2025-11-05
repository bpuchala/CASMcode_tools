import pytest
from ase.calculators.kim.kim import KIM

import casm.tools.map.neb as casm_neb
import libcasm.configuration as casmconfig
import libcasm.xtal as xtal
import libcasm.xtal.structures as xtal_structures
from casm.tools.map import map_to_prim
from casm.tools.shared.ase_utils import relax


@pytest.fixture
def bcc_to_hcp_mapping():

    # Requires:
    # - ASE with KIM installed and configured
    # - KIM model: MEAM_LAMMPS_DickelBaskesAslam_2018_MgAlZn__MO_093637366498_002

    # Prepare calculator:
    calculator = KIM("MEAM_LAMMPS_DickelBaskesAslam_2018_MgAlZn__MO_093637366498_002")
    fmax = 0.01

    # Prepare unrelaxed parent and child structures:
    bcc_init = xtal_structures.BCC(r=1.6, atom_type="Mg")
    hcp_init = xtal_structures.HCP(r=1.6, atom_type="Mg")

    # Relax parent and child structures:
    bcc_relaxed, bcc_Epot = relax(
        casm_structure=bcc_init,
        calculator=calculator,
        fmax=fmax,
    )

    hcp_relaxed, hcp_Epot = relax(
        casm_structure=hcp_init,
        calculator=calculator,
        fmax=fmax,
    )

    # Construct parent prim and child structure (without calculated properties):
    parent_prim = casmconfig.Prim(
        xtal_prim=xtal.Prim(
            lattice=bcc_relaxed.lattice(),
            coordinate_frac=bcc_relaxed.atom_coordinate_frac(),
            occ_dof=[[x] for x in bcc_relaxed.atom_type()],
        )
    )
    child = xtal.Structure(
        lattice=hcp_relaxed.lattice(),
        atom_coordinate_frac=hcp_relaxed.atom_coordinate_frac(),
        atom_type=hcp_relaxed.atom_type(),
    )

    # Map child to parent prim:
    results, uuids = map_to_prim(
        child=child,
        parent_prim=parent_prim,
    )
    structure_mapping = results[0]

    return parent_prim, child, structure_mapping, calculator


def test_NEBPath(bcc_to_hcp_mapping):
    parent_prim, child, structure_mapping, calculator = bcc_to_hcp_mapping

    # Construct NEB path using parent, child, and structure mapping:
    neb_path = casm_neb.NEBPath(
        parent_prim=parent_prim,
        child=child,
        structure_mapping=structure_mapping,
    )
    assert isinstance(neb_path, casm_neb.NEBPath)
    assert isinstance(neb_path.initial_structure, xtal.Structure)
    assert isinstance(neb_path.final_structure, xtal.Structure)

    # Construct interpolated images to initialize NEB path:
    structures = neb_path.make_path(
        n_images=5,
    )
    assert len(structures) == 5
    for i, struct in enumerate(structures):
        assert isinstance(struct, xtal.Structure)
    assert neb_path.initial_structure.is_equivalent_to(structures[0])
    assert neb_path.final_structure.is_equivalent_to(structures[-1])

    # Check interpolation coordinates:
    for i, struct in enumerate(structures):
        coord = neb_path.path_coordinates(struct)

        expected_parallel_coord = i / (len(structures) - 1)
        expected_perpendicular_coord = 0.0

        # lattice strain coordinates
        assert pytest.approx(coord[0], abs=1e-6) == expected_parallel_coord
        assert pytest.approx(coord[1], abs=1e-6) == expected_perpendicular_coord

        # displacement coordinates
        assert pytest.approx(coord[2], abs=1e-6) == expected_parallel_coord
        assert pytest.approx(coord[3], abs=1e-6) == expected_perpendicular_coord


def test_NEBPath_calculate(bcc_to_hcp_mapping):
    parent_prim, child, structure_mapping, calculator = bcc_to_hcp_mapping

    # Construct NEB path using parent, child, and structure mapping:
    neb_path = casm_neb.NEBPath(
        parent_prim=parent_prim,
        child=child,
        structure_mapping=structure_mapping,
    )

    path = neb_path.make_path(n_images=5)
    step_index = 0
    for image_index, structure in enumerate(path):
        image = neb_path.calculate(
            structure=structure,
            calculator=calculator,
            step_index=step_index,
            image_index=image_index,
        )
        assert isinstance(image, casm_neb.NEBImage)
        # print(f"Image {image_index}:")
        # print(image)
        # print()


def test_NEBImage_rmul(bcc_to_hcp_mapping):
    parent_prim, child, structure_mapping, calculator = bcc_to_hcp_mapping

    # Construct NEB path using parent, child, and structure mapping:
    neb_path = casm_neb.NEBPath(
        parent_prim=parent_prim,
        child=child,
        structure_mapping=structure_mapping,
    )

    path = neb_path.make_path(n_images=5)
    step_index = 0
    images = []
    for image_index, structure in enumerate(path):
        image = neb_path.calculate(
            structure=structure,
            calculator=calculator,
            step_index=step_index,
            image_index=image_index,
        )
        assert isinstance(image, casm_neb.NEBImage)
        images.append(image)

    # Create equivalent NEBImages via __rmul__:
    orbit = []
    for op in parent_prim.factor_group.elements:
        equiv = []
        for i, image in enumerate(images):
            transformed_image = op * image
            assert isinstance(transformed_image, casm_neb.NEBImage)
            equiv.append(transformed_image)
        orbit.append(equiv)

    assert len(orbit) == len(parent_prim.factor_group.elements)
    assert len(orbit) == 48


def test_NEBPath_make_chain_orbit(bcc_to_hcp_mapping):
    parent_prim, child, structure_mapping, calculator = bcc_to_hcp_mapping

    # Construct NEB path using parent, child, and structure mapping:
    neb_path = casm_neb.NEBPath(
        parent_prim=parent_prim,
        child=child,
        structure_mapping=structure_mapping,
    )

    path = neb_path.make_path(n_images=5)
    step_index = 0
    chain = []
    for image_index, structure in enumerate(path):
        image = neb_path.calculate(
            structure=structure,
            calculator=calculator,
            step_index=step_index,
            image_index=image_index,
        )
        assert isinstance(image, casm_neb.NEBImage)
        chain.append(image)

    # Excludes equivalent chains: 12 distinct chains
    orbit = neb_path.make_chain_orbit(chain_prototype=chain)
    assert len(orbit) == 12
