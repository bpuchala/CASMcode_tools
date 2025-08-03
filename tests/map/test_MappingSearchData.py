import numpy as np
import pytest

import libcasm.configuration as casmconfig
import libcasm.xtal as xtal
from casm.tools.map import MappingSearchData


def test_MappingSearchData_constructor_1():
    """Test the constructor of MappingSearchData."""

    with pytest.raises(TypeError):
        data = MappingSearchData()  # noqa: F841


def test_MappingSearchData_constructor_2(get_casm_structure):
    """Test the constructor of MappingSearchData."""

    parent_structure = get_casm_structure("BCC_Li.json")
    child = get_casm_structure("HCP_Li.json")

    data = MappingSearchData(
        parent=parent_structure,
        child=child,
    )

    assert isinstance(data, MappingSearchData)
    assert isinstance(data.parent_structure, xtal.Structure)
    assert data.parent_structure is parent_structure
    assert isinstance(data.parent_prim, casmconfig.Prim)
    assert isinstance(data.child, xtal.Structure)
    assert data.options is None
    assert data.mappings == []
    assert data.uuids == []
    assert data.options_history == []

    assert data.child_atom_types == ["Li"]
    assert data.child_atom_count == np.array([2])
    assert np.allclose(data.child_atom_frac, [1])

    assert data.parent_atom_types == ["Li"]
    assert data.parent_atom_count == np.array([1])
    assert np.allclose(data.parent_atom_frac, [1])

    data.validate_atom_types()


def test_MappingSearchData_constructor_3(get_casm_structure, get_prim):
    """Test the constructor of MappingSearchData."""

    parent_prim = get_prim("BCC_LiVa_prim.json")
    child = get_casm_structure("HCP_Li.json")

    data = MappingSearchData(
        parent=parent_prim,
        child=child,
    )

    assert isinstance(data, MappingSearchData)
    assert data.parent_structure is None
    assert data.parent_prim is parent_prim
    assert isinstance(data.parent_prim, casmconfig.Prim)
    assert isinstance(data.child, xtal.Structure)
    assert data.options is None
    assert data.mappings == []
    assert data.uuids == []
    assert data.options_history == []

    assert data.child_atom_types == ["Li"]
    assert data.child_atom_count == np.array([2])
    assert np.allclose(data.child_atom_frac, [1])

    assert data.parent_atom_types == ["Li", "Va"]
    with pytest.raises(ValueError):
        parent_atom_count = data.parent_atom_count  # noqa: F841

    with pytest.raises(ValueError):
        parent_atom_frac = data.parent_atom_frac  # noqa: F841

    data.validate_atom_types()
