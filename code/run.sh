#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=240:30:00
#SBATCH --job-name=hoomd
#SBATCH --output=hoomd-%j.out

# ==========================================================
module load anaconda/anaconda3
source activate hoomd

# ==========================================================
python3 ProteinBulk.py

# ==========================================================
# MDAnalysis package is not available in Stokes's Conda env
# Need to create own env and install it for analysis

conda activate my_env
python3 long_config_analysis.py