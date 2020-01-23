#!/bin/bash
#
# Re-run the code locally, to re-create the data and figure.
#
# Usage:
#
#   ./scripts/rerun.sh
#
#
rm -rf $(ls | egrep "example_23_3..")
Rscript example_23.R

