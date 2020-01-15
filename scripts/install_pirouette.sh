#!/bin/bash
#
# Script to install pirouette and its dependencies
# on the Peregrine computer cluster
#
# Usage:
# 
# sbatch install_pirouette.sh
#
#
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --mem=1G
#SBATCH --job-name=install_pirouette
#SBATCH --output=install_pirouette.log
module load R
Rscript -e 'devtools::install_github("richelbilderbeek/pirouette", dependencies = TRUE)'
Rscript -e 'devtools::install_github("richelbilderbeek/becosys")'
Rscript -e 'if (!beastier::is_beast2_installed()) beastier::install_beast2()'
Rscript -e 'if (!mauricer::is_beast2_pkg_installed("NS")) mauricer::install_beast2_pkg("NS")'
