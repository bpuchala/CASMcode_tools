import pathlib
import shutil

# Path to the examples directory:
examples_dir = pathlib.Path(__file__).resolve().parent

# Clean up generated files
dir = examples_dir / "map_Mg_bcc_hcp"
shutil.rmtree(dir / ".ipynb_checkpoints", ignore_errors=True)
