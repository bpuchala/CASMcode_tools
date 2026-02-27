import numpy as np
import pymatgen.core as pmg

import libcasm.xtal as xtal


def make_pmg_structure(
    casm_structure: xtal.Structure,
) -> pmg.Structure:
    """Convert a CASM Structure to a pymatgen Structure

    Notes
    -----

    This method creates a pymatgen Structure from the given CASM Structure's lattice,
    atom types, and atom fractional coordinates.

    The following convention is used to create pymatgen species from CASM atom types:

    - CASM "Va" / "VA" / "va" atom types are converted to pymatgen DummySpecie "X".
    - Other atom types are converted to either pymatgen Element, Species, or
      DummySpecies objects (in attempt order)

    Magnetic moments are included as site properties if present in the CASM Structure.
    The atom property "Cmagspin" (1 x N_sites array) is converted to a site property
    "magmom" with scalar values. The atom property "NCmagspin" (3 x N_sites array)
    or "SOmagspin" are converted to a site property "magmom" with 3-component vector
    values.

    Parameters
    ----------
    casm_structure: xtal.Structure
        The CASM Structure to convert.

    Returns
    -------
    pymatgen_structure: pmg.Structure
        The converted pymatgen Structure.
    """
    lattice = pmg.Lattice(casm_structure.lattice().column_vector_matrix().transpose())

    species = []
    for atom_type in casm_structure.atom_type():
        if atom_type.lower() == "va":
            _species = pmg.DummySpecies("X")
        else:
            try:
                _species = pmg.Element(atom_type)
            except Exception:
                try:
                    _species = pmg.Species(atom_type)
                except Exception:
                    _species = pmg.DummySpecies(atom_type)
        species.append(_species)

    if "Cmagspin" in casm_structure.atom_properties():
        magmoms = []
        Cmagspin = casm_structure.atom_properties()["Cmagspin"]
        for i in range(Cmagspin.shape[1]):
            magmoms.append(float(Cmagspin[0, i]))
        site_properties = {"magmom": magmoms}
    elif "NCmagspin" in casm_structure.atom_properties():
        magmoms = []
        NCmagspin = casm_structure.atom_properties()["NCmagspin"]
        for i in range(NCmagspin.shape[1]):
            magmoms.append(
                [
                    float(NCmagspin[0, i]),
                    float(NCmagspin[1, i]),
                    float(NCmagspin[2, i]),
                ]
            )
        site_properties = {"magmom": magmoms}
    elif "SOmagspin" in casm_structure.atom_properties():
        magmoms = []
        SOmagspin = casm_structure.atom_properties()["SOmagspin"]
        for i in range(SOmagspin.shape[1]):
            magmoms.append(
                [
                    float(SOmagspin[0, i]),
                    float(SOmagspin[1, i]),
                    float(SOmagspin[2, i]),
                ]
            )
        site_properties = {"magmom": magmoms}
    else:
        site_properties = None

    coords = casm_structure.atom_coordinate_frac().transpose()

    return pmg.Structure(
        lattice=lattice,
        species=species,
        coords=coords,
        coords_are_cartesian=False,
        site_properties=site_properties,
    )


def make_casm_structure(
    pmg_structure: pmg.Structure,
    magspin_flavor: str = "C",
):
    """Convert a pymatgen Structure to a CASM Structure

    Notes
    -----

    This method creates a CASM Structure from the given pymatgen Structure's lattice,
    atom types, and atom fractional coordinates.

    Atom types are constructed using the string representation of the pymatgen species
    at each site.

    Magnetic moments are included as site properties if present in the pymatgen
    Structure. The site property "magmom" is checked for each site. If present,
    it is included in the CASM Structure as follows:

    - If "magmom" is scalar valued (int or float) for all sites, it is included as
      an atom property "Cmagspin" (1 x N_sites array).
    - If "magmom" is vector valued (3-component list or tuple) for all sites, it is
      included as an atom property "magmom" (3 x N_sites array).
    - If "magmom" varies in type (scalar vs vector) across sites, a ValueError is
      raised.

    Parameters
    ----------
    pymatgen_structure: pmg.Structure
        The pymatgen Structure to convert.
    magspin_flavor: Optional[str] = "SO"
        The flavor prefix to use for vector-valued magnetic moments. Expected to be
        "SO" (spin-orbit coupled) or "NC" (non-collinear). The difference between
        "SO" and "NC" is that the Cartesian rotation matrix of symmetry-operations is
        applied to magnetic spins of type "SOmagspin" but not to "NCmagspin".

    Returns
    -------
    casm_structure: xtal.Structure
        The converted CASM Structure.
    """
    lattice = xtal.Lattice(
        column_vector_matrix=pmg_structure.lattice.matrix.transpose(),
    )
    atom_coordinate_frac = pmg_structure.frac_coords.transpose()
    atom_type = [str(site.specie) for site in pmg_structure.sites]
    n_sites = len(pmg_structure.sites)

    # Check if pmg_structure sites have site property "magmom" and check if it is
    # scalar or vector valued. Error if it varies from site to site.

    has_magmom = False
    has_scalar_magmom = None
    for site in pmg_structure.sites:
        if "magmom" in site.properties:
            has_magmom = True
            magmom = site.properties["magmom"]
            if isinstance(magmom, (int, float)):
                if has_scalar_magmom is None:
                    has_scalar_magmom = True
                elif has_scalar_magmom is False:
                    raise ValueError(
                        "Error in make_casm_structure: "
                        "Site property 'magmom' must be consistently scalar or vector "
                        "valued across all sites."
                    )
            elif isinstance(magmom, (list, tuple)) and len(magmom) == 3:
                if has_scalar_magmom is None:
                    has_scalar_magmom = False
                elif has_scalar_magmom is True:
                    raise ValueError(
                        "Error in make_casm_structure: "
                        "Site property 'magmom' must be consistently scalar or vector "
                        "valued across all sites."
                    )
            else:
                raise ValueError(
                    "Error in make_casm_structure: "
                    "Site property 'magmom' must be scalar or 3-component vector."
                )

    atom_properties = {}

    if has_magmom:
        if has_scalar_magmom:
            Cmagspin = np.zeros((1, n_sites))
            for i, site in enumerate(pmg_structure.sites):
                Cmagspin[0, i] = site.properties.get("magmom", 0.0)
            atom_properties["Cmagspin"] = Cmagspin
        else:
            label = magspin_flavor + "magspin"
            M = np.zeros((3, n_sites))
            for i, site in enumerate(pmg_structure.sites):
                M[:, i] = site.properties.get("magmom", [0.0, 0.0, 0.0])
            atom_properties[label] = M

    global_properties = {}

    return xtal.Structure(
        lattice=lattice,
        atom_coordinate_frac=atom_coordinate_frac,
        atom_type=atom_type,
        atom_properties=atom_properties,
        global_properties=global_properties,
    )
