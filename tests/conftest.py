import pathlib

import pytest

import libcasm.configuration as casmconfig
import libcasm.xtal as xtal
from casm.tools.shared.json_io import read_required


@pytest.fixture
def data_dir():
    tests_dir = pathlib.Path(__file__).parent
    return tests_dir / "data"


@pytest.fixture
def get_structure(data_dir):
    def _get_structure(name):
        return (data_dir / "structures" / name).as_posix()

    return _get_structure


@pytest.fixture
def get_casm_structure(get_structure):
    def _get_casm_structure(name):
        return xtal.Structure.from_dict(read_required(get_structure(name)))

    return _get_casm_structure


@pytest.fixture
def get_prim(get_structure):
    def _get_prim(name):
        return casmconfig.Prim.from_dict(read_required(get_structure(name)))

    return _get_prim


@pytest.fixture
def examples_dir(data_dir):
    return data_dir / "examples"


@pytest.fixture
def settings_dir(examples_dir):
    return dir / "settings"


@pytest.fixture
def dummy_vasp_potentials_dir(examples_dir):
    return examples_dir / "dummy_vasp_potentials"
