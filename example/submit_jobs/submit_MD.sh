#!/usr/bin/env bash
# Submit a standard MD production run.
# Usage: bash submit_jobs/submit_MD.sh

export AMBER_SCRIPTS_DIR="/orcd/pool/004/jhalpin/09-fragfold/RELE_simulations/amber_REMD/amber_scripts"

# === Job parameters ===
TEMP0=400
PROD_NS=50

# === System ===
PDB_IN="/orcd/pool/004/jhalpin/09-fragfold/RELE_simulations/amber_REMD/example/input_pdbs/DB2_unbound.pdb"
OUTBASE="$(basename "${PDB_IN%.*}")"

# === Output directory (auto-named from parameters) ===
OUTPUT_ROOT="/orcd/pool/004/jhalpin/09-fragfold/RELE_simulations/amber_REMD/example/outputs/output_MD"
OUTDIR="${OUTPUT_ROOT}/${OUTBASE}-${PROD_NS}ns-MD-${TEMP0}K"

sbatch --export=ALL,PROD_NS=$PROD_NS,TEMP0=$TEMP0,PDB_IN=$PDB_IN,OUTDIR=$OUTDIR,OUTBASE=$OUTBASE \
  "${AMBER_SCRIPTS_DIR}/MD.sbatch"
