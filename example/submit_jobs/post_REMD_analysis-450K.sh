#!/bin/bash
#SBATCH --job-name=remd_analysis
## ===== CONFIGURE FOR YOUR SITE: edit the 2 lines below =====
#SBATCH -p pi_keating     # <- your partition name
#SBATCH --gres=gpu:l40s:4 # <- your GPU type:count  (or remove if CPUs only)
## ===========================================================
#SBATCH --nodes=1
#SBATCH --ntasks=48           # Default tasks; override with sbatch -n/--ntasks
#SBATCH --cpus-per-task=1
#SBATCH --mem=250G
#SBATCH -t 4-00:00:00
#SBATCH -o ./logs/%x_%j.out
#SBATCH -e ./logs/%x_%j.err

AMBER_SCRIPTS_DIR="/orcd/pool/004/jhalpin/09-fragfold/RELE_simulations/amber_REMD/amber_scripts"
source "${AMBER_SCRIPTS_DIR}/../site_config.sh"
source activate "$CONDA_ENV"

# === Job parameters ===
REPLICAS=48
TOTAL_EXCHANGES=500000
T_MAX=450
T_ANALYSIS=450
EXCHANGE_EVERY_PS=1.0

# === Systems to analyze (set OUTDIR to the existing run's output directory) ===
OUTBASE="DB2_bound_unbound"
NS=$(awk "BEGIN {print ($TOTAL_EXCHANGES * $EXCHANGE_EVERY_PS) / 1000}")
OUTDIR="/orcd/pool/004/jhalpin/09-fragfold/RELE_simulations/amber_REMD/example/outputs/output_T-REMD/${OUTBASE}-${NS}ns-REMD-300-${T_MAX}K-${REPLICAS}reps-NVT"
bash "${AMBER_SCRIPTS_DIR}/reconstruct_300K-parallel.sh" "${OUTDIR}" "${OUTBASE}" "${OUTDIR}/${OUTBASE}_rem.log" $T_ANALYSIS $REPLICAS
# python "${AMBER_SCRIPTS_DIR}/rem_accept_cpptraj.py" "${OUTDIR}/${OUTBASE}_rem.log"

OUTBASE="DB2_unbound"
NS=$(awk "BEGIN {print ($TOTAL_EXCHANGES * $EXCHANGE_EVERY_PS) / 1000}")
OUTDIR="/orcd/pool/004/jhalpin/09-fragfold/RELE_simulations/amber_REMD/example/outputs/output_T-REMD/${OUTBASE}-${NS}ns-REMD-300-${T_MAX}K-${REPLICAS}reps-NVT"
bash "${AMBER_SCRIPTS_DIR}/reconstruct_300K-parallel.sh" "${OUTDIR}" "${OUTBASE}" "${OUTDIR}/${OUTBASE}_rem.log" $T_ANALYSIS $REPLICAS
# python "${AMBER_SCRIPTS_DIR}/rem_accept_cpptraj.py" "${OUTDIR}/${OUTBASE}_rem.log"

# python "${AMBER_SCRIPTS_DIR}/regen_REMD_export-all.py"
