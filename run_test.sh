#!/bin/bash

# Create conda environment from conda.yaml
mamba env create -f conda.yaml

# Activate the environment
conda activate laura-env

# compile Laura
cd laura
mkdir build
cd build
cmake .. -DLAURA_BUILD_ROOFIT_TASK=ON -DLAURA_BUILD_EXAMPLES=ON -DLAURA_BUILD_TESTS=ON
make -j$(nproc)

# Run tests with timing


time /home/almalinux/Laura++/build/examples/GenFit3pi gen 100
time /home/almalinux/Laura++/build/examples/GenFit3pi fit 0 100