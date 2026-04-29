import subprocess

from libcasm.xtal import Structure


def test_convert_vasp_to_json(data_dir, tmp_path):
    input_file = data_dir / "structures" / "BCC_Li.vasp"
    output_file = tmp_path / "BCC_Li.json"

    result = subprocess.run(
        ["casm-convert", str(input_file), str(output_file)],
        check=True,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    assert output_file.exists()

    import json

    data = json.loads(output_file.read_text())
    structure = Structure.from_dict(data)
    assert isinstance(structure, Structure)
    assert structure.atom_type() == ["Li"]


def test_convert_json_to_vasp(data_dir, tmp_path):
    input_file = data_dir / "structures" / "BCC_Li.vasp"
    json_file = tmp_path / "BCC_Li.json"
    output_file = tmp_path / "BCC_Li_out.vasp"

    subprocess.run(
        ["casm-convert", str(input_file), str(json_file)],
        check=True,
    )

    result = subprocess.run(
        ["casm-convert", str(json_file), str(output_file)],
        check=True,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    assert output_file.exists()

    structure = Structure.from_poscar_str(output_file.read_text())
    assert isinstance(structure, Structure)
    assert structure.atom_type() == ["Li"]


def test_convert_explicit_formats(data_dir, tmp_path):
    input_file = data_dir / "structures" / "HCP_Ag.vasp"
    output_file = tmp_path / "HCP_Ag_out.json"

    result = subprocess.run(
        [
            "casm-convert",
            "--input-format",
            "vasp",
            "--output-format",
            "casm",
            str(input_file),
            str(output_file),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    assert output_file.exists()

    import json

    data = json.loads(output_file.read_text())
    structure = Structure.from_dict(data)
    assert isinstance(structure, Structure)
    assert "Ag" in structure.atom_type()


def test_convert_no_overwrite(data_dir, tmp_path):
    input_file = data_dir / "structures" / "BCC_Li.vasp"
    output_file = tmp_path / "BCC_Li.json"

    subprocess.run(
        ["casm-convert", str(input_file), str(output_file)],
        check=True,
    )

    result = subprocess.run(
        ["casm-convert", str(input_file), str(output_file)],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    assert "skipping" in result.stdout


def test_convert_force_overwrite(data_dir, tmp_path):
    input_file = data_dir / "structures" / "BCC_Li.vasp"
    output_file = tmp_path / "BCC_Li.json"

    subprocess.run(
        ["casm-convert", str(input_file), str(output_file)],
        check=True,
    )

    result = subprocess.run(
        ["casm-convert", "--force", str(input_file), str(output_file)],
        check=True,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    assert output_file.exists()


def test_convert_missing_input(tmp_path):
    result = subprocess.run(
        [
            "casm-convert",
            str(tmp_path / "nonexistent.vasp"),
            str(tmp_path / "out.json"),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0
