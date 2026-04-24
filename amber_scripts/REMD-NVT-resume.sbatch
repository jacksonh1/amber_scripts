#!/bin/bash
#SBATCH --job-name=remd-nvt-resume
# NOTE: partition and GPU below are defaults only. run_scripts/submit_REMD-NVT-resume.sh
# overrides them from site_config.sh via --partition and --gres CLI flags.
#SBATCH -p pi_keating
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --gres=gpu:l40s:4
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH -t 4-00:00:00
#SBATCH -o ./logs/%x_%j.out
#SBATCH -e ./logs/%x_%j.err

# Usage: sbatch -t <walltime> -n <replicas> --export=ALL,OUTDIR=...,OUTBASE=... REMD-NVT-resume.sbatch
# Resumes a T-REMD NVT production run from existing *_remd.rst.* checkpoint files.
# Safe to re-run after cancellation/maintenance: auto-detects completed exchanges and
# writes each resumed segment to a new numbered trajectory (seg001, seg002, ...) in the
# original scratch directory. At completion, all segments are concatenated and post-processing
# runs unchanged.

set -e
set -o pipefail

: "${AMBER_SCRIPTS_DIR:?AMBER_SCRIPTS_DIR is not set — submit via submit_jobs/submit_*.sh or: export AMBER_SCRIPTS_DIR=/path/to/amber_scripts}"
source "${AMBER_SCRIPTS_DIR}/../site_config.sh"
CURRENT_TIME=$(date +%s)

export OMP_NUM_THREADS=1
REPLICAS=${SLURM_NTASKS:-${REPLICAS:-48}}
ENGINE="${ENGINE:-pmemd.cuda.MPI}"
TOTAL_EXCHANGES="${TOTAL_EXCHANGES:-500000}"

stop_mps() { echo quit | nvidia-cuda-mps-control >/dev/null 2>&1 || true; }
trap 'stop_mps' EXIT

echo "===== JOB START $(date) ====="

############################
# STEP 1: Validate parameters
############################
echo "===== STEP 1: Validate parameters ====="
: "${OUTDIR:?[ERROR] OUTDIR must be set via --export}"
: "${OUTBASE:?[ERROR] OUTBASE must be set via --export}"
PRMTOP="${OUTDIR}/${OUTBASE}.prmtop"
REMLOG="${OUTDIR}/${OUTBASE}_rem.log"
GROUPFILE_ORIG="${OUTDIR}/${OUTBASE}_groupfile"
[[ -f "$PRMTOP" ]]        || { echo "[ERROR] prmtop not found: $PRMTOP"; exit 1; }
[[ -f "$REMLOG" ]]        || { echo "[ERROR] rem.log not found: $REMLOG"; exit 1; }
[[ -f "$GROUPFILE_ORIG" ]] || { echo "[ERROR] Groupfile not found: $GROUPFILE_ORIG"; exit 1; }
mkdir -p "${OUTDIR}/inputs" "${OUTDIR}/outputs" "${OUTDIR}/logs" \
         "${OUTDIR}/mdinfo" "${OUTDIR}/trajectories" ./logs
echo "[OK] OUTDIR=$OUTDIR  OUTBASE=$OUTBASE"

############################
# STEP 2: Count completed exchanges
############################
echo "===== STEP 2: Count completed exchanges ====="
COMPLETED_EXCHANGES=$(grep -c "^# exchange" "$REMLOG" || true)
COMPLETED_EXCHANGES=${COMPLETED_EXCHANGES:-0}
REMAINING=$((TOTAL_EXCHANGES - COMPLETED_EXCHANGES))
echo "[INFO] Completed: $COMPLETED_EXCHANGES / $TOTAL_EXCHANGES  (remaining: $REMAINING)"
if [[ "$REMAINING" -le 0 ]]; then
  echo "[INFO] All exchanges complete. Nothing to do."
  exit 0
fi

############################
# STEP 3: Locate original scratch directory
############################
echo "===== STEP 3: Locate scratch directory ====="
FIRST_X_PATH=$(awk '{for(i=1;i<=NF;i++) if($i=="-x"){print $(i+1); exit}}' "$GROUPFILE_ORIG")
SCRATCH_DIR=$(dirname "$FIRST_X_PATH")
[[ -d "$SCRATCH_DIR" ]] || { echo "[ERROR] Scratch dir not found: $SCRATCH_DIR"; exit 1; }
SCRATCH_AVAIL=$(df -BG "$SCRATCH_DIR" | awk 'NR==2{gsub(/G/,"",$4); print $4}')
[[ "$SCRATCH_AVAIL" -lt 50 ]] && { echo "[ERROR] Low disk space: ${SCRATCH_AVAIL}GB"; exit 1; }
echo "[INFO] Scratch: $SCRATCH_DIR  (${SCRATCH_AVAIL}GB free)"

############################
# STEP 4: Determine segment tag (supports multiple restarts)
############################
echo "===== STEP 4: Determine segment tag ====="
NEXT_SEG=1
while [[ -f "${SCRATCH_DIR}/${OUTBASE}_remd_seg$(printf "%03d" $NEXT_SEG).nc.000" ]]; do
  NEXT_SEG=$((NEXT_SEG+1))
done
SEG_TAG="seg$(printf "%03d" $NEXT_SEG)"
echo "[INFO] Segment tag: $SEG_TAG"

############################
# STEP 5: Load environment
############################
echo "===== STEP 5: Loading environment ====="
source activate "$CONDA_ENV"
module purge
module load gcc
module load "$CUDA_MODULE"
module load "$OPENMPI_MODULE"
source "$AMBER_SH"
echo "[OK] Modules loaded"

############################
# STEP 6: Start MPS
############################
echo "===== STEP 6: Starting MPS ====="
export CUDA_DEVICE_ORDER=PCI_BUS_ID
export CUDA_VISIBLE_DEVICES="0,1,2,3"
export CUDA_MPS_PIPE_DIRECTORY="${SLURM_TMPDIR:-/tmp}/mps_${SLURM_JOB_ID}"
export CUDA_MPS_LOG_DIRECTORY="${SLURM_TMPDIR:-/tmp}/mps_${SLURM_JOB_ID}"
mkdir -p "${CUDA_MPS_PIPE_DIRECTORY}" "${CUDA_MPS_LOG_DIRECTORY}"
nvidia-cuda-mps-control -d
sleep 2
echo "[OK] MPS started"

############################
# STEP 7: Generate per-replica input files
############################
echo "===== STEP 7: Generating input files (numexchg=${REMAINING}) ====="
for i in $(seq 0 $((REPLICAS-1))); do
  idx=$(printf "%03d" "$i")
  ORIG_IN="${OUTDIR}/inputs/${OUTBASE}_remd.in.${idx}"
  [[ -f "$ORIG_IN" ]] || { echo "[ERROR] Original input not found: $ORIG_IN"; exit 1; }
  sed "s/numexchg=[0-9]*/numexchg=${REMAINING}/" "$ORIG_IN" \
    > "${OUTDIR}/inputs/${OUTBASE}_remd_${SEG_TAG}.in.${idx}"
done
echo "[OK] Input files written"

############################
# STEP 8: Generate groupfile
############################
echo "===== STEP 8: Generating groupfile ====="
RESUME_GROUP="${OUTDIR}/${OUTBASE}_groupfile_${SEG_TAG}"
: > "$RESUME_GROUP"
for i in $(seq 0 $((REPLICAS-1))); do
  idx=$(printf "%03d" "$i")
  RST="${OUTDIR}/restarts/${OUTBASE}_remd.rst.${idx}"
  [[ -f "$RST" ]] || { echo "[ERROR] REMD restart not found: $RST"; exit 1; }
  printf -- "-O -rem 1 -remlog %s -i %s -p %s -c %s -o %s -r %s -x %s -inf %s -l %s\n" \
    "$REMLOG" \
    "${OUTDIR}/inputs/${OUTBASE}_remd_${SEG_TAG}.in.${idx}" \
    "$PRMTOP" \
    "$RST" \
    "${OUTDIR}/outputs/${OUTBASE}_remd_${SEG_TAG}.out.${idx}" \
    "$RST" \
    "${SCRATCH_DIR}/${OUTBASE}_remd_${SEG_TAG}.nc.${idx}" \
    "${OUTDIR}/mdinfo/mdinfo_remd_${SEG_TAG}.${idx}" \
    "${OUTDIR}/logs/logfile_remd_${SEG_TAG}.${idx}" \
    >> "$RESUME_GROUP"
done
echo "[OK] Groupfile: $RESUME_GROUP"

############################
# STEP 9: Run REMD
############################
echo "===== STEP 9: Running REMD segment $SEG_TAG ($REMAINING exchanges) ====="
echo "[INFO] Progress: check $REMLOG for exchange statistics"
mpirun -np "$REPLICAS" "$ENGINE" -ng "$REPLICAS" -groupfile "$RESUME_GROUP"
echo "[OK] REMD segment $SEG_TAG done"

############################
# STEP 10: Concatenate trajectory segments
############################
echo "===== STEP 10: Concatenating trajectory segments ====="
CONCAT_FAIL=0
for i in $(seq 0 $((REPLICAS-1))); do
  idx=$(printf "%03d" "$i")

  # Collect all segments: original (seg0) then numbered resumes, in sorted order
  SEG_FILES=()
  [[ -f "${SCRATCH_DIR}/${OUTBASE}_remd.nc.${idx}" ]] \
    && SEG_FILES+=("${SCRATCH_DIR}/${OUTBASE}_remd.nc.${idx}")
  for sf in $(ls "${SCRATCH_DIR}/${OUTBASE}_remd_seg"???.nc.${idx} 2>/dev/null | sort); do
    SEG_FILES+=("$sf")
  done

  if [[ ${#SEG_FILES[@]} -eq 0 ]]; then
    echo "[WARN] No trajectory files for replica ${idx}"; CONCAT_FAIL=1; continue
  fi

  if [[ ${#SEG_FILES[@]} -eq 1 ]]; then
    ln -sf "${SEG_FILES[0]}" "${OUTDIR}/trajectories/${OUTBASE}_remd.nc.${idx}"
    continue
  fi

  COMBINED="${SCRATCH_DIR}/${OUTBASE}_remd_combined.nc.${idx}"
  CONCAT_IN=$(mktemp)
  for f in "${SEG_FILES[@]}"; do printf 'trajin %s\n' "$f"; done > "$CONCAT_IN"
  printf 'trajout %s netcdf\nrun\nquit\n' "$COMBINED" >> "$CONCAT_IN"

  if cpptraj -p "$PRMTOP" -i "$CONCAT_IN"; then
    ln -sf "$COMBINED" "${OUTDIR}/trajectories/${OUTBASE}_remd.nc.${idx}"
  else
    echo "[WARN] cpptraj concat failed for replica ${idx}"; CONCAT_FAIL=1
  fi
  rm -f "$CONCAT_IN"
done

if [[ "$CONCAT_FAIL" -eq 0 ]]; then
  echo "[OK] All 48 replicas concatenated → ${OUTDIR}/trajectories/"
else
  echo "[WARN] Some trajectory concatenations failed — check warnings above"
fi

############################
# STEP 11: Post-processing
############################
echo "===== STEP 11: Post-processing ====="
ambpdb -p "$PRMTOP" -c "${OUTDIR}/restarts/${OUTBASE}_remd.rst.000" \
  > "${OUTDIR}/${OUTBASE}_final_000.pdb" || true

bash "${AMBER_SCRIPTS_DIR}/reconstruct_300K-parallel.sh" \
  "${OUTDIR}" "${OUTBASE}" "${REMLOG}" 300 "${REPLICAS}" || true

python "${AMBER_SCRIPTS_DIR}/rem_accept_cpptraj.py" "${REMLOG}" || true
python "${AMBER_SCRIPTS_DIR}/gen_REMD_export.py" --export-folder="300K_export" --remd-folder="${OUTDIR}" --temperature=300 --overwrite

ELAPSED_S="$(( $(date +%s) - CURRENT_TIME ))"
echo ""
echo "=========================================="
echo "===== RESUME COMPLETE $(date) ====="
echo "[INFO] Elapsed: ${ELAPSED_S}s"
echo "[INFO] Scratch: $SCRATCH_DIR"
echo "[INFO] Trajectories: ${OUTDIR}/trajectories/"
echo "[INFO] Exchange log: $REMLOG"
echo "To re-run after maintenance: re-submit submit_REMD-NVT-resume.sh (safe; auto-detects progress)"
echo "=========================================="
