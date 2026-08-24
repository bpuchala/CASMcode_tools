import numpy as np

import libcasm.configuration as casmconfig
import libcasm.xtal as xtal
import libcasm.xtal.structures as xtal_structures
from casm.tools.map import methods


def bcc_Mg_prim(a: float = 3.55):
    """Construct a bcc Mg prim, as used by the ``map_Mg_bcc_hcp`` example."""
    return casmconfig.Prim(
        xtal_prim=xtal.Prim.from_atom_coordinates(
            structure=xtal_structures.BCC(a=a, atom_type="Mg")
        )
    )


def skewed_bcc_superstructure():
    """Construct a det(T)=8 bcc Mg superstructure, in a skewed cell

    This is the ``f=0`` image (the ideal parent superstructure) of one of the
    det(T)=8 bcc->hcp mappings found by the ``map_Mg_bcc_hcp`` example. The cell is
    skewed with respect to the bcc primitive cell, which is what makes it a
    regression test: symmetry analysis of the structure as given used to raise
    "Error in make_inverse_index: multiplication table identity error", because
    `xtal::make_factor_group` returned a factor group whose first element was not
    the identity.
    """
    return xtal.Structure(
        lattice=xtal.Lattice(
            np.array(
                [
                    [-1.775, 5.325, -1.775],
                    [-1.775, -1.775, 8.875],
                    [-1.775, -1.775, -5.325],
                ]
            )
        ),
        atom_coordinate_frac=np.array(
            [
                [0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875],
                [0.0, 0.625, 0.25, 0.875, 0.5, 0.125, 0.75, 0.375],
                [0.0, 0.750, 0.50, 0.250, 0.0, 0.750, 0.50, 0.250],
            ]
        ),
        atom_type=["Mg"] * 8,
    )


def test_make_chain_info_skewed_det8_supercell():
    """Chain info can be constructed for a det(T)=8 image in a skewed cell"""
    parent_prim = bcc_Mg_prim()
    structure = skewed_bcc_superstructure()

    chain_info = methods.make_chain_info(
        chain=[structure],
        parent_prim=parent_prim,
    )

    assert len(chain_info) == 1
    info = chain_info[0]

    # The crystal is bcc: its primitive cell is the parent prim, with the 48 element
    # prim factor group, Im-3m (229) -- even though the structure was given in a
    # det(T)=8 cell
    assert np.isclose(info["volume"], 1.0)
    assert info["factor_group_size"] == 48
    assert info["spacegroup_type"]["number"] == 229
    assert info["spacegroup_type"]["international_short"] == "Im-3m"


def test_make_chain_info_is_independent_of_supercell_choice():
    """The info describes the crystal, not the supercell it is given in

    All three structures are bcc Mg, expressed in cells of different size and shape,
    so all three must give the same result.
    """
    parent_prim = bcc_Mg_prim()
    bcc = xtal_structures.BCC(a=3.55, atom_type="Mg")

    chain = [
        bcc,
        xtal.make_superstructure(np.array([[2, 0, 0], [0, 2, 0], [0, 0, 2]]), bcc),
        skewed_bcc_superstructure(),
    ]

    chain_info = methods.make_chain_info(chain=chain, parent_prim=parent_prim)

    for info in chain_info:
        assert np.isclose(info["volume"], 1.0)
        assert info["factor_group_size"] == 48
        assert info["spacegroup_type"]["number"] == 229


def test_make_chain_info_volume_is_relative_to_the_parent_prim():
    """The volume is that of the primitive cell, relative to the parent prim"""
    parent_prim = bcc_Mg_prim()
    bcc = xtal_structures.BCC(a=3.55, atom_type="Mg")
    hcp = xtal_structures.HCP(r=1.6, atom_type="Mg")

    chain_info = methods.make_chain_info(chain=[bcc, hcp], parent_prim=parent_prim)

    expected_hcp_volume = abs(hcp.lattice().volume()) / abs(bcc.lattice().volume())
    assert not np.isclose(expected_hcp_volume, 1.0)

    assert np.isclose(chain_info[0]["volume"], 1.0)
    assert np.isclose(chain_info[1]["volume"], expected_hcp_volume)

    # hcp: P6_3/mmc (194), with a 24 element prim factor group
    assert chain_info[1]["factor_group_size"] == 24
    assert chain_info[1]["spacegroup_type"]["number"] == 194


def test_make_supercell_info_volume_with_inexact_determinant():
    """The supercell volume is rounded, not truncated

    ``np.linalg.det`` of this integer matrix evaluates to 7.999999999999998. It is
    the parent supercell of one of the det(T)=8 bcc->hcp mappings found by the
    ``map_Mg_bcc_hcp`` example.
    """
    T = np.array([[1, 1, -1], [1, -1, 2], [1, -1, -2]])
    assert np.linalg.det(T) != 8.0

    supercell = casmconfig.Supercell(
        prim=bcc_Mg_prim(),
        transformation_matrix_to_super=T,
    )
    info = methods.make_supercell_info(supercell=supercell)

    assert info["volume"] == 8


def test_make_chain_info_is_primitive_flag():
    """`is_primitive=True` skips making the structures primitive"""
    parent_prim = bcc_Mg_prim()
    bcc = xtal_structures.BCC(a=3.55, atom_type="Mg")
    hcp = xtal_structures.HCP(r=1.6, atom_type="Mg")

    # Both structures are already primitive, so the flag changes nothing
    chain = [bcc, hcp]
    assert methods.make_chain_info(
        chain=chain, parent_prim=parent_prim, is_primitive=True
    ) == methods.make_chain_info(chain=chain, parent_prim=parent_prim)

    # The promise is not checked: a structure that is not primitive is used as given,
    # and the results then describe that superstructure and not the crystal
    info = methods.make_chain_info(
        chain=[skewed_bcc_superstructure()],
        parent_prim=parent_prim,
        is_primitive=True,
    )[0]
    assert np.isclose(info["volume"], 8.0)
    assert info["factor_group_size"] == 32
    assert info["spacegroup_type"]["number"] != 229
