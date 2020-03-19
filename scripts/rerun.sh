#!/bin/bash
#
# Re-run the code locally, to re-create the data and figure.
#
# Usage:
#
#   ./scripts/rerun.sh
#
#SBATCH --partition=gelifes
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --mem=10G
#SBATCH --job-name=pirex23
#SBATCH --output=example_23.log
#
rm -rf example_23
time Rscript example_23.R
zip -r pirouette_example_23.zip example_23 example_23.R scripts errors_low.png errors_mid.png errors_high.png

