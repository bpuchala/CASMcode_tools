<img alt="Shows the CASM logo" src="https://raw.githubusercontent.com/prisms-center/CASMcode_global/main/python/doc/_static/logo.svg" width="600" />

#### casm-tools

The casm-tools package provides pure Python CLI tools and helper functions. This includes:

- casm-calc: Setup, run, and report results of structure calculations
- casm-map: Structure mapping and import


#### Install

    pip install casm-tools


#### Usage

See the [casm docs](https://prisms-center.github.io/CASMcode_pydocs/casm/overview/latest/).


#### Testing

Some tests require [ASE](https://ase-lib.org/) and some also require [OpenKIM](https://openkim.org/).

Skip some tests using the `-m` option to select tests by marker. For example, to skip tests that require either ASE or OpenKIM, run:

    pytest -rsap -x -m "not requires_kim and not requires_ase" tests

or to skip tests that require OpenKIM but not ASE, run:

    pytest -rsap -x -m "not requires_kim" tests

See the [OpenKIM documentation](https://openkim.org/) for more information about installing and configuring OpenKIM.

Example configuration:

    source kim-api-activate
    export KIM_API_PORTABLE_MODELS_DIR="$HOME/.kim-api/2.4.1+Clang.Clang.GNU.2025-04-16-07-14-10/portable-models-dir"
    export KIM_API_MODEL_DRIVERS_DIR="$HOME/.kim-api/2.4.1+Clang.Clang.GNU.2025-04-16-07-14-10/model-drivers-dir"
    kim-api-collections-management install user MEAM_LAMMPS_KimJeonLee_2015_MgCa__MO_611309973581_002



#### Instal Zsh Completions

To install Zsh completions, cp the completions function to your custom completions directory:

    cp completions/zsh/_casm_map ~/.oh-my-zsh/custom/completions/_casm-map

Then reload your shell or run:

    compinit


#### About CASM

The casm-tools package is part of the [CASM](https://prisms-center.github.io/CASMcode_docs/) open source software package, which is designed to perform first-principles statistical mechanical studies of multi-component crystalline solids.

CASM is developed by the Van der Ven group, originally at the University of Michigan and currently at the University of California Santa Barbara.

For more information, see the [CASM homepage](https://prisms-center.github.io/CASMcode_docs/).


#### License

GNU Lesser General Public License (LGPL). Please see the file LICENSE for details.

