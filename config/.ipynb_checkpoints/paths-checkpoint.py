"""
Path configuration for the TLS lung HDM project

Edit BASE_DIR to the location of local data directory
"""

from pathlib import Path

# Local data directory with Xenium outputs
DATA_DIR = Path("/home/workspace/data/temp/mouse_lung")

# Local output directory
BASE_OUTDIR = Path("/home/workspace/spatial_mouse_lung_outputs")

# Lodal directory with scRNAseq ref
REF_DIR = Path("/home/workspace/data/temp/hurskainen_ref")

# Local inputs dir, functions dir
INPUTS_DIR = Path("/home/workspace/spatial_tls_manuscript/inputs")
FUNCTIONS_DIR = Path("/home/workspace/spatial_tls_manuscript/utils")
