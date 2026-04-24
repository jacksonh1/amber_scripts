#!/usr/bin/env bash
# Resume a T-REMD NPT run after job time-out or maintenance.
# Safe to re-run: auto-detects completed exchanges and continues from where it left off.
# Usage: bash submit_jobs/submit_REMD-NPT-resume.sh

export AMBER_SCRIPTS_DIR="/orcd/pool/004/jhalpin/09-fragfold/RELE_simulations/amber_REMD/amber_scripts"

# === Job parameters (must match the original submission) ===
REPLICAS=48
TOTAL_EXCHANGES=250000

# === Point OUTDIR at the existing run's output directory ===
OUTBASE="helix_fusion"
OUTDIR="/orcd/pool/004/jhalpin/09-fragfold/RELE_simulations/amber_REMD/example/outputs/output_T-REMD/${OUTBASE}-250ns-REMD-300-450K-${REPLICAS}reps-NPT"

sbatch -t 13:00:00 -n "$REPLICAS" \
  --export=ALL,OUTDIR="$OUTDIR",OUTBASE="$OUTBASE",REPLICAS="$REPLICAS",TOTAL_EXCHANGES="$TOTAL_EXCHANGES" \
  "${AMBER_SCRIPTS_DIR}/REMD-NPT-resume.sbatch"
