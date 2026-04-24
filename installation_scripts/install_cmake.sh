#!/bin/bash
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -p pi_keating
#SBATCH --mem=20000
#SBATCH -o ./pmemd24_install%j.out
#SBATCH -e ./pmemd24_install%j.err
#SBATCH --constraint="rocky8"
#SBATCH --nodelist=node[3619-3620]


module purge
module load deprecated-modules
module load cuda/12.4.0
module load openmpi/5.0.8

wget https://github.com/Kitware/CMake/releases/download/v4.0.4/cmake-4.0.4.tar.gz
tar xzf ./cmake-4.0.4.tar.gz
cd cmake-4.0.4
./bootstrap --prefix=$HOME/opt/cmake-4.0
