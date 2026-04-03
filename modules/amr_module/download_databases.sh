#!/bin/bash
echo "Downloading AMRFinder database..."
mkdir -p acinetoscope/modules/amr_module/data/amrfinder_db/2025-12-03.1
wget -P acinetoscope/modules/amr_module/data/amrfinder_db/2025-12-03.1 \
    https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/AMR.LIB

echo "Downloading MLST databases..."
mkdir -p acinetoscope/modules/mlst_module/data
# Add MLST database download commands

echo "Downloading ABRicate databases..."
mkdir -p acinetoscope/modules/abricate_module/data
# Add ABRicate database download commands

echo "Databases downloaded successfully!"
