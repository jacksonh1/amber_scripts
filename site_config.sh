#!/usr/bin/env bash
# site_config.sh — Site-specific configuration for the Amber REMD/MD pipeline.
#
# HOW TO USE:
#   Edit the values below once for your cluster/installation.
#   Sourced automatically by amber_scripts/*.sbatch when they run on compute nodes.
#   submit_jobs/ scripts do NOT source this — they only need AMBER_SCRIPTS_DIR set.
#
#   Any variable can be overridden by exporting it before sourcing, e.g.:
#     SCRATCH_ROOT=/tmp/myscratch sbatch ...

# === Amber environment ===
# Path to amber.sh from your Amber/pmemd installation (sources AmberTools env vars).
AMBER_SH="${AMBER_SH:-/orcd/pool/004/jhalpin/pmemd24/amber.sh}"

# Conda environment name to activate (must have ambertools/cpptraj).
CONDA_ENV="${CONDA_ENV:-ambertools25}"

# === Scratch storage ===
# Root for per-job scratch directories. Each job creates a subdirectory:
#   ${SCRATCH_ROOT}/${SLURM_JOB_ID}_${timestamp}
# This should be on fast local or parallel scratch storage.
SCRATCH_ROOT="${SCRATCH_ROOT:-/orcd/data/keating/001/jhalpin/MD}"

# === HPC environment modules ===
CUDA_MODULE="${CUDA_MODULE:-cuda/12.9.1}"
OPENMPI_MODULE="${OPENMPI_MODULE:-openmpi/4.1.4}"
