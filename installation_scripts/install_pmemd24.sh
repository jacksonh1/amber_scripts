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
export PATH="$HOME/opt/cmake-4.0/bin:$PATH"
module load deprecated-modules
module load cuda/12.4.0
module load openmpi/5.0.8

cd ./pmemd24_src/build
./run_cmake
make -j 8
make install
