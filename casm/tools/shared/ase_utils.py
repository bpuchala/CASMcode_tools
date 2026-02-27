"""Interface with `ASE <https://wiki.fysik.dtu.dk/ase/>`_"""

import pathlib
import typing

import ase
import ase.calculators.vasp
import ase.io
import numpy as np

import casm.tools.shared.json_io as json_io
import casm.tools.shared.opt_utils as opt_utils
import libcasm.configuration as casmconfig
import libcasm.mapping.info as mapinfo
import libcasm.mapping.methods as mapmethods
import libcasm.xtal as xtal


def make_ase_atoms(casm_structure: xtal.Structure) -> ase.Atoms:
    """Given a CASM Structure, convert it to an ASE Atoms

    .. attention::

        This method only works for non-magnetic atomic structures. If the structure
        contains molecular information, an error will be raised.

    Notes
    -----

    This method converts a CASM Structure object to an ASE Atoms object. It includes:

    - the lattice vectors
    - the atomic positions
    - the atomic types

    Parameters
    ----------
    casm_structure : libcasm.xtal.Structure

    Returns
    -------
    ase.Atoms

    """
    if len(casm_structure.mol_type()):
        raise ValueError(
            "Error: only non-magnetic atomic structures may be converted using "
            "to_ase_atoms"
        )

    symbols = casm_structure.atom_type()
    positions = casm_structure.atom_coordinate_cart().transpose()
    cell = casm_structure.lattice().column_vector_matrix().transpose()

    return ase.Atoms(
        symbols=symbols,
        positions=positions,
        cell=cell,
        pbc=True,
    )


def make_casm_structure(ase_atoms: ase.Atoms) -> xtal.Structure:
    """Given an ASE Atoms, convert it to a CASM Structure

    .. attention::

        This method only works for non-magnetic atomic structures.

    Notes
    -----

    This method converts an ASE Atoms object to a CASM Structure object. It includes:

    - the lattice vectors
    - atomic positions
    - atomic types, using the ASE chemical symbols.
    - Optional properties, if available from an ASE calculator:

      - Forces are added as atom properties named `"force"`.
      - Potential energy is added as the global property name `"energy"`.

    Parameters
    ----------
    ase_atoms : ase.Atoms
        A :class:`ase.Atoms` object

    Returns
    -------
    casm_structure: libcasm.xtal.Structure
        A :class:`~libcasm.xtal.Structure` object

    """

    lattice = xtal.Lattice(
        column_vector_matrix=ase_atoms.get_cell().transpose(),
    )
    atom_coordinate_frac = ase_atoms.get_scaled_positions().transpose()
    atom_type = ase_atoms.get_chemical_symbols()

    atom_properties = {}
    global_properties = {}
    if ase_atoms._calc is not None:
        try:
            forces = ase_atoms.get_forces()
            atom_properties["force"] = forces.transpose()
        except Exception:
            pass

        try:
            energy = ase_atoms.get_potential_energy()
            global_properties["energy"] = np.array([[energy]])
        except Exception:
            pass

    return xtal.Structure(
        lattice=lattice,
        atom_coordinate_frac=atom_coordinate_frac,
        atom_type=atom_type,
        atom_properties=atom_properties,
        global_properties=global_properties,
    )


def write_structure_using_ase(
    casm_structure: xtal.Structure,
    path: pathlib.Path,
    format: typing.Optional[str] = None,
    make_ase_atoms_f: typing.Optional[
        typing.Callable[[xtal.Structure], ase.Atoms]
    ] = None,
) -> None:
    """Write a structure using ASE's write function.

    .. attention::

        This method does not write magnetic moments.

    Parameters
    ----------
    casm_structure : libcasm.xtal.Structure
        The CASM Structure to write.
    path : pathlib.Path
        The path to the file where the structure will be written.
    format : Optional[str]=None
        The format to use for writing the file. If None, ASE will try to infer the
        format from the file extension.
    make_ase_atoms_f : Optional[Callable[[libcasm.xtal.Structure], ase.Atoms]] = None
        A function to convert the CASM structure to an ASE Atoms object. If None,
        the default function, :func:`make_ase_atoms` is used, which works for
        non-magnetic atomic structures.

    """
    try:
        import ase.io
    except ImportError:
        raise ImportError(
            "ASE is not installed. "
            "Please install ASE to write this structure format."
        )
    if make_ase_atoms_f is None:
        make_ase_atoms_f = make_ase_atoms

    ase_atoms = make_ase_atoms_f(casm_structure)
    ase.io.write(path.as_posix(), ase_atoms, format=format)


def write_structure_traj_using_ase(
    casm_structure_traj: list[xtal.Structure],
    path: pathlib.Path,
    format: typing.Optional[str] = None,
    make_ase_atoms_f: typing.Optional[
        typing.Callable[[xtal.Structure], ase.Atoms]
    ] = None,
) -> None:
    """Write a structure using ASE's write function.

    .. attention::

        This method does not write magnetic moments.

    Parameters
    ----------
    casm_structure_traj : list[libcasm.xtal.Structure]
        The CASM Structure trajectory to write.
    path : pathlib.Path
        The path to the file where the structure will be written.
    format : Optional[str]=None
        The format to use for writing the file. If None, ASE will try to infer the
        format from the file extension.
    make_ase_atoms_f : Optional[Callable[[libcasm.xtal.Structure], ase.Atoms]] = None
        A function to convert the CASM structure to an ASE Atoms object. If None,
        the default function, :func:`make_ase_atoms` is used, which works for
        non-magnetic atomic structures.

    """
    try:
        import ase.io
    except ImportError:
        raise ImportError(
            "ASE is not installed. "
            "Please install ASE to write this structure format."
        )
    if make_ase_atoms_f is None:
        make_ase_atoms_f = make_ase_atoms

    ase_atoms = [make_ase_atoms_f(x) for x in casm_structure_traj]
    ase.io.write(path.as_posix(), ase_atoms, format=format)


def read_structure_using_ase(
    path: pathlib.Path,
    format: typing.Optional[str] = None,
    make_casm_structure_f: typing.Optional[
        typing.Callable[[ase.Atoms], xtal.Structure]
    ] = None,
) -> xtal.Structure:
    """Read a structure using ASE's read function.

    .. attention::

        This method does not read magnetic moments.

    Parameters
    ----------
    path : pathlib.Path
        The path to the structure file.
    format : Optional[str]=None
        The format to use for reading the file. If None, ASE will try to infer the
        format from the file extension.
    make_casm_structure_f : typing.Callable[[ase.Atoms], libcasm.xtal.Structure]
        A function to convert an ASE Atoms object to a CASM structure. If None, the
        default function, :func:`make_casm_structure`, is used, which works for
        non-magnetic atomic structures.

    Returns
    -------
    casm_structure: libcasm.xtal.Structure
        A CASM Structure read from the file.

    """
    try:
        import ase.io
    except ImportError:
        raise ImportError(
            "ASE is not installed. Please install ASE to read this structure format."
        )
    if make_casm_structure_f is None:
        make_casm_structure_f = make_casm_structure

    return make_casm_structure_f(ase.io.read(path.as_posix(), format=format))


def read_structure_traj_using_ase(
    path: pathlib.Path,
    format: typing.Optional[str] = None,
    make_casm_structure_f: typing.Optional[
        typing.Callable[[ase.Atoms], xtal.Structure]
    ] = None,
) -> list[xtal.Structure]:
    """Read a structure trajectory using ASE's read function.

    .. attention::

        This method does not read magnetic moments.

    Parameters
    ----------
    path : pathlib.Path
        The path to the structure file.
    format : Optional[str]=None
        The format to use for reading the file. If None, ASE will try to infer the
        format from the file extension.
    make_casm_structure_f : typing.Callable[[ase.Atoms], libcasm.xtal.Structure]
        A function to convert an ASE Atoms object to a CASM structure. If None, the
        default function, :func:`make_casm_structure`, is used, which works for
        non-magnetic atomic structures.

    Returns
    -------
    casm_structure: libcasm.xtal.Structure
        A CASM Structure read from the file.

    """
    try:
        import ase.io
    except ImportError:
        raise ImportError(
            "ASE is not installed. Please install ASE to read this structure format."
        )
    if make_casm_structure_f is None:
        make_casm_structure_f = make_casm_structure

    return [
        make_casm_structure_f(x)
        for x in ase.io.read(path.as_posix(), format=format, index=":")
    ]


class AseVaspTool:
    def __init__(
        self,
        calctype_settings_dir: typing.Optional[pathlib.Path] = None,
        make_ase_atoms_f: typing.Optional[
            typing.Callable[[xtal.Structure], ase.Atoms]
        ] = None,
        make_casm_structure_f: typing.Optional[
            typing.Callable[[ase.Atoms], xtal.Structure]
        ] = None,
    ):
        """Setup, run, and collect VASP calculations using ASE.

        For details on the parameters and environment configuration, see the
        `ASE Vasp calculator documentation <https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html>`_.

        .. attention::

            ASE assumes that POTCAR files exist in one of `potpaw_PBE`, `potpaw`, or
            `potpaw_GGA`, located at the path specified by the environment
            variable VASP_PP_PATH.

        Parameters
        ----------
        calctype_settings_dir: Optional[pathlib.Path] = None
            Path to the directory containing settings, including a `calc`.json file for
            VASP calculator constructor arguments, and INCAR and KPOINTS
            template files. All files are optional, but if they exist, they will be
            used to set up the VASP calculator.
        make_ase_atoms_f: Optional[Callable[[libcasm.xtal.Structure], ase.Atoms]] = None
            A function to convert the CASM structure to an ASE Atoms object. The
            default function, :func:`make_ase_atoms` works for non-magnetic atomic
            structures.
        make_casm_structure_f: \
        Optional[Callable[[ase.Atoms], libcasm.xtal.Structure]] = None
            A function to convert an ASE Atoms object to a CASM structure. The
            default function, :func:`make_casm_structure` works for non-magnetic atomic
            structures.
        """
        if make_ase_atoms_f is None:
            make_ase_atoms_f = make_ase_atoms
        if make_casm_structure_f is None:
            make_casm_structure_f = make_casm_structure

        ### Read settings from the calctype_settings_dir if provided

        self.calctype_settings_dir = calctype_settings_dir
        """Optional[pathlib.Path]: Path to the directory containing settings, 
        including a `calc.json` file for 
        `ASE VASP calculator <https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html#module-ase.calculators.vasp>`_
        constructor arguments, and template INCAR and KPOINTS files. All files are 
        optional, but if they exist, they will be used to set up the VASP 
        calculations."""

        incar_path = None
        kpoints_path = None
        settings = {}
        if self.calctype_settings_dir is not None:
            settings = json_io.read_required(
                path=self.calctype_settings_dir / "calc.json"
            )

            incar_path = calctype_settings_dir / "INCAR"
            if not incar_path.exists():
                incar_path = None

            kpoints_path = calctype_settings_dir / "KPOINTS"
            if not kpoints_path.exists():
                kpoints_path = None

        self.settings = settings
        """Optional[dict]: Settings read from the calc.json file in the 
        calctype_settings_dir, which will be passed to the 
        :class:`ase.calculators.vasp.Vasp` calculator constructor."""

        self.incar_path = incar_path
        """Optional[pathlib.Path]: Path to the template INCAR file, if it exists."""

        self.kpoints_path = kpoints_path
        """Optional[pathlib.Path]: Path to the template KPOINTS file, if it exists."""

        ### Functions to convert between CASM Structure and ASE Atoms

        self._make_ase_atoms_f = make_ase_atoms_f
        """Callable[[libcasm.xtal.Structure], ase.Atoms]: Function to convert CASM 
        Structure to ASE Atoms."""

        self._make_casm_structure_f = make_casm_structure_f
        """Callable[[ase.Atoms], libcasm.xtal.Structure]: Function to convert an 
        ASE Atoms to CASM Structure."""

    def make_calculator(
        self,
        ase_atoms: ase.Atoms,
        calc_dir: pathlib.Path,
    ) -> ase.calculators.vasp.Vasp:
        """Construct an ASE VASP calculator.

        Parameters
        ----------
        ase_atoms: ase.Atoms
            The ASE Atoms object to use for the calculation.
        calc_dir: pathlib.Path
            The directory in which to store the calculation files. This directory will
            be created if it does not exist.

        Returns
        -------
        vasp_calculator: ase.calculators.vasp.Vasp
            The ASE VASP calculator object, configured with the provided settings and
            paths to INCAR and KPOINTS files if they exist.

        """
        calc_dir.mkdir(parents=True, exist_ok=True)

        vasp_calculator = ase.calculators.vasp.Vasp(
            atoms=ase_atoms,
            directory=calc_dir,
            **self.settings,
        )

        if self.incar_path is not None:
            vasp_calculator.read_incar(self.incar_path)

        if self.kpoints_path is not None:
            vasp_calculator.read_kpoints(self.kpoints_path)

        return vasp_calculator

    def setup(
        self,
        casm_structure: xtal.Structure,
        calc_dir: pathlib.Path,
        config: typing.Optional[casmconfig.Configuration] = None,
    ) -> ase.calculators.vasp.Vasp:
        """Setup a VASP calculation for a given structure.

        Parameters
        ----------
        casm_structure: libcasm.xtal.Structure
            The structure to calculate. The structure is written to the calculation
            directory as `structure.json`.
        calc_dir: pathlib.Path
            The directory in which to store the calculation files.
        config: Optional[libcasm.configuration.Configuration] = None
            If provided, the configuration object associated with the structure is
            printed to the calculation directory as `config.json`.

        Returns
        -------
        vasp_calculator: ase.calculators.vasp.Vasp
            The ASE VASP calculator object.
        """
        ase_atoms = self._make_ase_atoms_f(casm_structure)
        vasp_calculator = self.make_calculator(ase_atoms=ase_atoms, calc_dir=calc_dir)

        # Write INCAR, KPOINTS, POTCAR, POSCAR
        vasp_calculator.write_input(atoms=ase_atoms)

        return vasp_calculator

    def report(
        self,
        calc_dir: typing.Optional[pathlib.Path] = None,
        index: typing.Any = None,
    ) -> typing.Union[xtal.Structure, list[xtal.Structure]]:
        """Report the results of a VASP calculation.

        Parameters
        ----------
        calc_dir: pathlib.Path
            The directory containing the VASP calculation files.
        index: int, slice or str
            Indicates the structures to return. By default, only the last structure
            is returned. Use `":"` to return all structures.

        Returns
        -------
        results: Union[libcasm.xtal.Structure, list[libcasm.xtal.Structure]]
            A CASM Structure or a list of CASM Structures, as specified by `index`.
        """

        outcar_file = calc_dir / "OUTCAR"
        if outcar_file.exists():
            value = ase.io.read(outcar_file, format="vasp-out", index=index)
        else:
            outcar_gz_file = calc_dir / "OUTCAR.gz"
            if outcar_gz_file.exists():
                import gzip

                with gzip.open(outcar_gz_file, "rt") as f:
                    value = ase.io.read(f, format="vasp-out", index=index)
            else:
                raise FileNotFoundError(
                    f"Error in AseVaspTool.report: "
                    f"Neither OUTCAR nor OUTCAR.gz found in {calc_dir.as_posix()}"
                )

        if isinstance(value, ase.Atoms):
            results = self._make_casm_structure_f(value)
        elif isinstance(value, list):
            results = [self._make_casm_structure_f(x) for x in value]
        else:
            raise Exception(f"Unrecognized type {type(value)} from ase.io.read")

        return results


def relax(
    casm_structure: xtal.Structure,
    calculator: typing.Any,
    fmax: float,
    logfile: typing.Any = None,
):
    """Relax a CASM structure using ASE BFGS optimizer

    Notes
    -----

    Relaxes the lattice and coordinates using:

    .. code-block:: python

        from ase.filters import UnitCellFilter
        from ase.optimize import BFGS

        atoms = make_ase_atoms(casm_structure)
        atoms.calc = calculator
        optimizer = BFGS(UnitCellFilter(atoms), logfile=logfile)
        optimizer.run(fmax=fmax)


    Parameters
    ----------
    casm_structure : libcasm.xtal.Structure
        The CASM structure to relax
    calculator : typing.Any
        An ASE calculator to use for energy and force calculations
    fmax : float
        The maximum force tolerance for the relaxation
    logfile: typing.Any = None
        The ASE BFGS optimizer `logfile` parameter. If logfile is a string, a file with
        that name will be opened. Use ‘-’ for stdout. It may be a file object, Path,
        str, or None. Default is None.

    Returns
    -------
    relaxed_structure : libcasm.xtal.Structure
        The relaxed CASM structure
    potential_energy : float
        The potential energy of the relaxed structure
    """
    from ase.filters import UnitCellFilter
    from ase.optimize import BFGSLineSearch

    atoms = make_ase_atoms(casm_structure)
    atoms.calc = calculator
    # optimizer = BFGS(UnitCellFilter(atoms), logfile=logfile)
    optimizer = BFGSLineSearch(UnitCellFilter(atoms), logfile=logfile)
    optimizer.run(fmax=fmax)
    potential_energy = atoms.get_potential_energy()
    relaxed_structure = make_casm_structure(atoms)
    return (relaxed_structure, potential_energy)


def _make_alignment_symop_to_Cartesian_axes(
    lattice: xtal.Lattice,
) -> xtal.SymOp:
    # align structure so 'a' is along 'x', 'b' in 'xy' plane:
    # L_aligned = Q * L, solve for Q:
    L = lattice.column_vector_matrix()
    a = L[:, 0]
    b = L[:, 1]
    c = L[:, 2]
    a_mag = np.linalg.norm(a)
    a_norm = a / a_mag
    b_proj = b - np.dot(b, a_norm) * a_norm
    b_norm = b_proj / np.linalg.norm(b_proj)
    a_aligned = np.array([a_mag, 0.0, 0.0])
    b_aligned = np.array([np.dot(b, a_norm), np.linalg.norm(b_proj), 0.0])
    c_aligned = np.array(
        [np.dot(c, a_norm), np.dot(c, b_norm), np.dot(c, np.cross(a_norm, b_norm))]
    )
    L_aligned = np.column_stack((a_aligned, b_aligned, c_aligned))
    Q = L_aligned @ np.linalg.inv(L)

    op = xtal.SymOp(matrix=Q, translation=np.zeros(3), time_reversal=False)

    return op


def relax_NEB(
    parent_lattice: xtal.Lattice,
    child: xtal.Structure,
    scored_structure_mapping: mapinfo.ScoredStructureMapping,
    calculator: typing.Any,
    fmax: float,
    max_iterations: int = 1000,
    n_images: int = 7,
):
    """Perform a SS-NEB calculation using `tsase.neb.ssneb` with `fire_ssneb` optimizer

    Parameters
    ----------
    parent_lattice : libcasm.xtal.Lattice
        The parent lattice.
    child : libcasm.xtal.Structure
        The child structure that was mapped.
    scored_structure_mapping : libcasm.mapping.info.ScoredStructureMapping
        The structure mapping from the parent to the child.
    calculator : typing.Any
        An ASE calculator to use for energy and force calculations.
    fmax : float
        The maximum force tolerance for the relaxation.
    max_iterations : Optional[int] = 1000
        The maximum number of iterations for the relaxation. Default is 1000.
    n_images: Optional[int] = 7
        The number of images in the NEB chain, including the endpoints.

    Returns
    -------
    relaxed_structures : list[libcasm.xtal.Structure]
        The relaxed CASM structures for each image in the NEB chain. Structures
        include calculated properties as collected by the `make_casm_structure`
        function. Structures are aligned to the original parent / mapped child
        orientation.
    """

    from casm.tools.map.methods import (
        make_child_transformation_matrix_to_super,
    )

    T_child = make_child_transformation_matrix_to_super(
        parent_lattice=parent_lattice,
        child_lattice=child.lattice(),
        structure_mapping=scored_structure_mapping,
    )
    T_child = np.round(T_child).astype(int)
    superchild = xtal.make_superstructure(
        transformation_matrix_to_super=T_child,
        structure=child,
    )

    # Generate endpoints from structure mapping:
    structures = []
    for i, f in enumerate([0.0, 1.0]):
        structure = mapmethods.make_mapped_structure(
            structure_mapping=scored_structure_mapping.interpolated(f),
            unmapped_structure=superchild,
        )
        structures.append(structure)

    # Align structures to Cartesian axes for SS-NEB:
    symop_big_alignment = _make_alignment_symop_to_Cartesian_axes(
        structures[0].lattice()
    )

    aligned_structures = []
    for structure in structures:
        x = symop_big_alignment * structure
        symop_small_alignment = _make_alignment_symop_to_Cartesian_axes(x.lattice())
        aligned_structure = symop_small_alignment * x
        aligned_structures.append(aligned_structure)

    # Create ASE atoms for endpoints:
    images = []
    for i, aligned_structure in enumerate(aligned_structures):

        atoms = make_ase_atoms(casm_structure=aligned_structure)
        atoms.calc = calculator
        images.append(atoms)

    # Create reference mappings for aligned endpoints:
    lattice_mapping_ref, atom_mapping_ref = mapmethods.direct_structure_mapping(
        structure1=aligned_structures[0],
        structure2=aligned_structures[1],
    )

    # Perform SS-NEB calculation:

    # Other parameter defaults:
    # k = 5.0, tangent = "new", dneb = False, dnebOrg = False, method = 'normal',
    # onlyci = False, weight = 1, parallel = False, ss = True,
    # express = numpy.zeros((3,3)), fixstrain = numpy.ones((3,3))
    #
    # Parameters:
    #         p1.......... one endpoint of the path
    #         p2.......... the other endpoint of the path
    #         numImages... the total number of images in the path, including the
    #                      endpoints
    #         k........... the spring force constant
    #         tangent..... the tangent method to use, "new" for the new tangent,
    #                      anything else for the old tangent
    #         dneb........ set to true to use the double-nudging method
    #         dnebOrg..... set to true to use the original double-nudging method
    #         method...... "ci" for the climbing image method, anything else for
    #                      normal NEB method
    #         onlyci...... (no description)
    #         weight...... (no description)
    #         ss.......... boolean, solid-state dimer or regular dimer
    #         express..... external press, 3*3 lower triangular matrix in the
    #                      unit of GPa
    #         fixstrain... 3*3 matrix as express.
    #                      0 fixes strain at the corresponding direction

    import shutil

    import tsase

    from casm.tools.shared.contexts import captured_output

    # Make a temporary working directory for NEB output:
    tmp_dir = pathlib.Path("_relax_neb_temp_working_dir")
    if tmp_dir.exists():
        shutil.rmtree(tmp_dir)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    # Run NEB calculation:
    with captured_output(wd=tmp_dir) as (sout, serr):
        ssneb = tsase.neb.ssneb(images[0], images[1], numImages=n_images)
        opt = tsase.neb.fire_ssneb(ssneb, maxmove=0.1, dtmax=0.1, dt=0.1)
        converged = opt.minimize(forceConverged=fmax, maxIterations=max_iterations)

    # Align relaxed structures back to original parent / mapped child orientation:
    ref_parent_structure = structures[0]

    relaxed_structures = []
    for image in ssneb.path:
        structure = make_casm_structure(image)
        # Q * U * L1 * T1 * N = L2; Q.inv == Q.T; L2_mapped = Q.T * L2
        lmap = mapmethods.map_lattices_without_reorientation(
            lattice1=ref_parent_structure.lattice(),
            lattice2=structure.lattice(),
        )
        Q = lmap.isometry()
        op = xtal.SymOp(matrix=Q.T, translation=np.zeros(3), time_reversal=False)
        relaxed_structures.append(op * structure)

    assert structures[0].is_equivalent_to(relaxed_structures[0])
    assert structures[1].is_equivalent_to(relaxed_structures[-1])

    return (relaxed_structures, converged)


class StrainDispVarAseCalcTool:
    """Tool for constructing energy and gradient functions in terms of strain
    and displacement variables for use in optimization routines.

    Notes
    -----
    - Uses an ASE calculator for energy, force, and stress evaluations.
    - Caches results in memory to make it easier to track the number of unique
      evaluations.
    - Intended for use with fast calculators, but could be adapted for slower ones by
      storing results on disk instead of in memory.

    """

    def __init__(
        self,
        var_tool: opt_utils.StrainDispVarTool,
        calculator: typing.Any,
    ):
        """

        .. rubric:: Constructor

        Parameters
        ----------
        var_tool: opt_utils.StrainDispVarTool
            Tool to convert between CASM structures and strain and displacement
            variable arrays.

        calculator: typing.Any
            ASE calculator to use for energy, force, and stress evaluations.
        """

        self.var_tool = var_tool
        """opt_utils.StrainDispVarTool: Tool to convert between CASM structures and 
        variable arrays."""

        self.calculator = calculator
        """typing.Any: ASE calculator to use for energy, force, and stress 
        evaluations."""

        self.x_list = []
        """List[np.ndarray]: List of x arrays corresponding to each calculation."""

        self.potential_energy_list = []
        """List[float]: List of potential energies corresponding to each x in x_list."""

        self.force_list = []
        """List[np.ndarray]: List of forces corresponding to each x in x_list,
        stored as (3, N_atoms) arrays."""

        self.stress_list = []
        """List[np.ndarray]: List of stresses corresponding to each x in x_list,
        stored as (6,) arrays in Kelvin notation."""

    def get_index(self, x: np.ndarray) -> int:
        """Return the index of x in the history, or None if not found.

        Parameters
        ----------
        x : np.ndarray
            1D array of variables.

        Returns
        -------
        index : Optional[int]
            The index of x in the history, or None if not found. Exact matches only.
        """
        from casm.tools.shared.conversions import voigt_to_kelvin

        for i, x_hist in enumerate(self.x_list):
            if (x == x_hist).all():
                return i

        # If not found, calculate and store:
        i = len(self.x_list)

        structure = self.var_tool.make_structure(x)
        atoms = make_ase_atoms(structure)
        atoms.calc = self.calculator

        energy = atoms.get_potential_energy()
        forces = atoms.get_forces()
        force = forces.transpose()
        stress = voigt_to_kelvin(atoms.get_stress(voigt=True))

        self.x_list.append(x)
        self.potential_energy_list.append(energy)
        self.force_list.append(force)
        self.stress_list.append(stress)

        return i

    def get_potential_energy(
        self,
        x: np.ndarray,
    ) -> float:
        """Return potential energy for given x

        Notes
        -----

        - Checks if x is already in history.
        - If it is, returns stored potential energy.
        - If not, calculates the potential energy, forces, and stress, and stores them,
          then returns the potential energy.

        Parameters
        ----------
        x : np.ndarray
            1D array of variables.

        Returns
        -------
        potential_energy : float
            The potential energy corresponding to x.
        """
        i = self.get_index(x)
        potential_energy = self.potential_energy_list[i]
        # print(f"**Eval {i}: Energy = {potential_energy:.6f} eV**")
        return potential_energy

    def get_grad(
        self,
        x: np.ndarray,
        dx: typing.Union[float, np.ndarray],
    ):
        """Calculate the potential energy gradient using finite differences

        Parameters
        ----------
        x : np.ndarray
            1D array of variables.
        dx : Union[float, np.ndarray]
            The finite difference step size. If a float is provided, the same step size
            is used for all variables. If an array is provided, it must have the same
            shape as x.

        Returns
        -------
        grad : np.ndarray
            The gradient of the potential energy with respect to x, as a 1D array,
            calculated using central finite differences.
        """
        n_vars = x.shape[0]
        grad = np.zeros(n_vars)

        for i in range(n_vars):
            x_plus = x.copy()
            x_minus = x.copy()

            if isinstance(dx, float):
                delta = dx
            else:
                delta = dx[i]

            x_plus[i] += delta
            x_minus[i] -= delta

            f_plus = self.get_potential_energy(x_plus)
            f_minus = self.get_potential_energy(x_minus)

            grad[i] = (f_plus - f_minus) / (2 * delta)

        return grad

    def get_hess(
        self,
        x: np.ndarray,
        dx: typing.Union[float, np.ndarray],
    ):
        """Calculate the Hessian using finite differences

        Parameters
        ----------
        x : np.ndarray
            1D array of variables.
        dx : Union[float, np.ndarray]
            The finite difference step size. If a float is provided, the same step size
            is used for all variables. If an array is provided, it must have the same
            shape as x.

        Returns
        -------
        hess : np.ndarray
            The Hessian of the potential energy with respect to x, as a 2D array,
            calculated using central finite differences.
        """
        n_vars = x.shape[0]
        hess = np.zeros((n_vars, n_vars))

        for i in range(n_vars):
            for j in range(n_vars):

                if i == j:
                    x_plus = x.copy()
                    x_minus = x.copy()

                    if isinstance(dx, float):
                        delta = dx
                    else:
                        delta = dx[i]

                    x_plus[i] += delta
                    x_minus[i] -= delta

                    f_plus = self.get_potential_energy(x_plus)
                    f_minus = self.get_potential_energy(x_minus)
                    f_curr = self.get_potential_energy(x)

                    hess[i, i] = (f_plus - 2 * f_curr + f_minus) / (delta**2)
                    continue

                elif j < i:

                    x_pp = x.copy()
                    x_pm = x.copy()
                    x_mp = x.copy()
                    x_mm = x.copy()

                    if isinstance(dx, float):
                        delta_i = dx
                        delta_j = dx
                    else:
                        delta_i = dx[i]
                        delta_j = dx[j]

                    x_pp[i] += delta_i
                    x_pp[j] += delta_j

                    x_pm[i] += delta_i
                    x_pm[j] -= delta_j

                    x_mp[i] -= delta_i
                    x_mp[j] += delta_j

                    x_mm[i] -= delta_i
                    x_mm[j] -= delta_j

                    f_pp = self.get_potential_energy(x_pp)
                    f_pm = self.get_potential_energy(x_pm)
                    f_mp = self.get_potential_energy(x_mp)
                    f_mm = self.get_potential_energy(x_mm)

                    hess[i, j] = (f_pp - f_pm - f_mp + f_mm) / (4 * delta_i * delta_j)
                    hess[j, i] = hess[i, j]
                    continue

        return hess

    def get_force(
        self,
        x: np.ndarray,
    ) -> np.ndarray:
        """Return forces on atoms for given x

        Parameters
        ----------
        x : np.ndarray
            1D array of variables.

        Returns
        -------
        force : np.ndarray
            The force, as a (3, N) matrix, corresponding to x.
        """
        return self.force_list[self.get_index(x)]

    def get_stress(
        self,
        x: np.ndarray,
    ) -> np.ndarray:
        """Return stress for given x

        Parameters
        ----------
        x : np.ndarray
            1D array of variables.

        Returns
        -------
        stress : np.ndarray
            The stress, as a (6,) array in Kelvin notation, corresponding to x.
        """
        return self.stress_list[self.get_index(x)]

    def make_potential_energy_f(self):
        """Return a function that computes potential energy for given x

        Returns
        -------
        potential_energy_f : Callable[[np.ndarray], float]
            Function that takes a 1D array x and returns the potential energy.

        """

        def potential_energy_f(x: np.ndarray) -> float:
            return self.get_potential_energy(x)

        return potential_energy_f

    def make_grad_f(self):
        """Return a function that computes the gradient for given x

        Returns
        -------
        grad_f : Callable[[np.ndarray], np.ndarray]
            Function that takes a 1D array x and returns the gradient of the potential
            energy with respect to x as a 1D array.
        """

        def grad_f(x: np.ndarray) -> np.ndarray:
            force = self.get_force(x)
            stress = self.get_stress(x)
            return self.var_tool.make_grad(
                x=x,
                force=force,
                stress=stress,
            )

        return grad_f
