#!/bin/bash

#install miniconda and mamba
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p $HOME/miniconda
export PATH="$HOME/miniconda/bin:$PATH" 
conda init bash
source ~/.bashrc

# Create conda environment from conda.yaml
mamba env create -f conda.yaml -y

# Activate the environment
conda activate laura-env

# compile Laura
cd Laura++
mkdir build
cd build
cmake .. -DLAURA_BUILD_ROOFIT_TASK=ON -DLAURA_BUILD_EXAMPLES=ON -DLAURA_BUILD_TESTS=ON
make -j$(nproc)

# Run tests with timing

localtime=$(date +"%Y-%m-%d %H:%M:%S")
echo "Starting tests at: $localtime"

time ./examples/GenFit3pi gen 100
time ./examples/GenFit3pi fit 0 100

endtime=$(date +"%Y-%m-%d %H:%M:%S")
echo "Finished tests at: $endtime"