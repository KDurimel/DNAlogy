#!/bin/bash
# Not really added in shell environment
echo 'alias fastGSEA="python $(pwd)'/src/fastGSEA.py'"' >> ~/.bashrc
# Create a new environment from packages.yml
conda env create -f packages.yml
# Activate the new environment
source activate gsea_env
