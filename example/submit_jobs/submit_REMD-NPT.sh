#!/usr/bin/env bash
# Submit a fresh T-REMD NPT run.
# Usage: bash submit_jobs/submit_REMD-NPT.sh

export AMBER_SCRIPTS_DIR="/orcd/pool/004/jhalpin/09-fragfold/RELE_simulations/amber_REMD/amber_scripts"

# === Job parameters ===
REPLICAS=48
TOTAL_EXCHANGES=1000
T_MIN=300
T_MAX=450
EXCHANGE_EVERY_PS=1.0

# === System ===
PDB_IN="/orcd/pool/004/jhalpin/09-fragfold/RELE_simulations/amber_REMD/example/input_pdbs/helix_fusion.pdb"
OUTBASE="$(basename "${PDB_IN%.*}")"

# === Output directory (auto-named from parameters) ===
OUTPUT_ROOT="/orcd/pool/004/jhalpin/09-fragfold/RELE_simulations/amber_REMD/example/outputs/output_T-REMD"
NS=$(awk "BEGIN {print ($TOTAL_EXCHANGES * $EXCHANGE_EVERY_PS) / 1000}")
OUTDIR="${OUTPUT_ROOT}/${OUTBASE}-${NS}ns-REMD-${T_MIN}-${T_MAX}K-${REPLICAS}reps-NPT"

sbatch -n "$REPLICAS" \
  --export=ALL,EXCHANGE_EVERY_PS=$EXCHANGE_EVERY_PS,TOTAL_EXCHANGES=$TOTAL_EXCHANGES,T_MIN=$T_MIN,T_MAX=$T_MAX,PDB_IN=$PDB_IN,OUTDIR=$OUTDIR,OUTBASE=$OUTBASE \
  "${AMBER_SCRIPTS_DIR}/REMD-NPT-singlenode.sbatch"
