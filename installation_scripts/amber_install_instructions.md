# installing Amber

## intro

use this website for info - https://ambermd.org/GetAmber.php

You need ambertools25 and PMEMD24. Ambertools is the free software with analysis tools and I think you can also do simulations but PMEMD24 has the software for faster simulations (like the ability to use parallel cpu and GPU). PMEMD24 is commercial but free for academic use. You'll still need ambertools25 though.


you may want to handle all the src and downloaded files in a different directory. Just copy the install scripts there, they expect the src code to be in the same directory


## Step 1 - install ambertools25 in a conda environment
instructions adapted from:
https://ambermd.org/GetAmber.php

```bash
conda create --name AmberTools25 python=3.12
conda activate AmberTools25
conda config --add channels conda-forge
conda config --set channel_priority strict
conda install dacase::ambertools-dac=25
conda install -c conda-forge mdanalysis
```

Alternatively, recreate the full environment from the provided spec (recommended for reproducibility):

```bash
conda env create -f installation_scripts/environment.yml
```

To use the code:
```bash
conda activate AmberTools25
```
there is a file that you can source to use the free version of the simulation stuff `source $CONDA_PREFIX/amber.sh`. Don't do this if you are going to use PMEMD.

## Step 2 - install PMEMD24

### step 2a - install cmake-4.0.4
1. run the install script: `sbatch ./install_cmake.sh`
hopefully that will just work


### step 2b - install pmemd24
steps:
1. Go here - https://ambermd.org/GetAmber.php and go to the section "Getting Amber24 for non-commercial use"
2. fill out the form. Download the source code. File should be called `pmemd24.tar.bz2`
3. unzip it (`tar -xvjf pmemd24.tar.bz2`)
4. run the install script `sbatch ./install_pmemd24.sh`
- note that if you need to troubleshoot the install, you should clean the build folder before trying again. navigate to `pmemd24_src/build` and run `./clean_build`
5. after compiling, there willl be a `pmemd24` folder next to the `pmemd24_src` folder. In that folder, there will be an `amber.sh` file that needs to be sourced prior to running simulations. In the `site_config.sh` file the `AMBER_SH` variable must be set to the path of this file. (`AMBER_SH="${AMBER_SH:-/orcd/pool/004/jhalpin/pmemd24/amber.sh}"`)


