#!/bin/bash
#SBATCH --job-name=remd
# NOTE: partition and GPU below are defaults only. run_scripts/submit_REMD-NVT.sh overrides
# them from site_config.sh via --partition and --gres CLI flags at submission time.
#SBATCH -p pi_keating
#SBATCH --nodes=1
#SBATCH --ntasks=48           # Default tasks; override with sbatch -n/--ntasks
#SBATCH --gres=gpu:l40s:4     # 4 GPUs per node
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH -t 4-00:00:00
#SBATCH -o ./logs/%x_%j.out
#SBATCH -e ./logs/%x_%j.err

# #SBATCH --nodelist=node[3619-3620]
# #SBATCH --nodelist=node[3620,3500]
# Exit on any error
set -e
set -o pipefail


: "${AMBER_SCRIPTS_DIR:?AMBER_SCRIPTS_DIR is not set — submit via submit_jobs/submit_*.sh or: export AMBER_SCRIPTS_DIR=/path/to/amber_scripts}"
source "${AMBER_SCRIPTS_DIR}/../site_config.sh"
SCRATCH_TS="$(date +%Y%m%d-%H%M%S)"
SCRATCH_DIR="${SCRATCH_ROOT}/${SLURM_JOB_ID:-$$}_${SCRATCH_TS}"
mkdir -p "$SCRATCH_DIR"

CURRENT_TIME=$(date +%s)
# Track current step number so we can decide whether to preserve scratch on failure.
CURRENT_STEP=0
# From which step onward we PRESERVE the scratch directory on failure (configurable).
PRESERVE_FROM_STEP=${PRESERVE_FROM_STEP:-12}


export OMP_NUM_THREADS=1 #need to be the same as -c value; Not sure what this does, it might be a relic. I don't think you ever want to set -c to anything other than 1
# Use scheduler-provided task count if available; else env; else default
REPLICAS=${SLURM_NTASKS:-${REPLICAS:-48}}

# Cleanup function for errors. Preserve scratch if we have reached (or passed)
# the configured PRESERVE_FROM_STEP (default: 12) so later-stage failures retain data.
cleanup_scratch_on_error() {
  local rc=$?
  if [ -n "$SCRATCH_DIR" ] && [ -d "$SCRATCH_DIR" ]; then
    if [ "$CURRENT_STEP" -ge "$PRESERVE_FROM_STEP" ]; then
      echo "[ERROR] Job failed at step $CURRENT_STEP (rc=$rc); PRESERVING scratch: $SCRATCH_DIR"
    else
      echo "[ERROR] Job failed at step $CURRENT_STEP (rc=$rc); cleaning scratch: $SCRATCH_DIR"
      rm -rf "$SCRATCH_DIR"
    fi
  fi
  exit $rc
}
trap 'cleanup_scratch_on_error' ERR

# Check available disk space in scratch (need ~100GB for 48 replicas × 20ns)
SCRATCH_AVAIL=$(df -BG "$SCRATCH_DIR" | awk 'NR==2 {print $4}' | sed 's/G//')
if [ "$SCRATCH_AVAIL" -lt 100 ]; then
    echo "[ERROR] Low disk space in scratch: ${SCRATCH_AVAIL}GB available (need ~100GB)"
    exit 1
fi

echo "===== JOB START $(date) ====="

# module load cuda
# module load openmpi
# source ~/.bashrc
source activate "$CONDA_ENV"
module purge
module load gcc
module load "$CUDA_MODULE"
module load "$OPENMPI_MODULE"
source "$AMBER_SH"
echo "[OK] Modules loaded"


####usage of this script: sbatch --export=ALL,PARAM1=VALUE1,... REMD.sbatch
############################
# ==== USER PARAMS ==== #
############################
CURRENT_STEP=1
echo "===== STEP 1: Setting parameters ====="
PDB_IN="${PDB_IN:-IL7-2-V3-cys.pdb}"
OUTBASE="${OUTBASE:-IL7-2-V3-cys}"
OUTDIR="${OUTDIR:-remd_${OUTBASE}}"
USE_PDB4AMBER="${USE_PDB4AMBER:-1}"
FORCE="${FORCE:-0}"

BOX_TYPE="${BOX_TYPE:-TIP3PBOX}"
BOX_SHAPE="${BOX_SHAPE:-Oct}"
BOX_BUFFER="${BOX_BUFFER:-10.0}"
NEUTRALIZE="${NEUTRALIZE:-1}"
SALT_MOLAR="${SALT_MOLAR:-0.15}"

ENGINE="${ENGINE:-pmemd.cuda.MPI}"
DT_PS="${DT_PS:-0.002}"
EXCHANGE_EVERY_PS="${EXCHANGE_EVERY_PS:-0.5}"
TOTAL_EXCHANGES="${TOTAL_EXCHANGES:-40000}"
TEMPS_LIST="${TEMPS_LIST:-}"
T_MIN="${T_MIN:-300}"
T_MAX="${T_MAX:-400}"

GAMMA_LN="${GAMMA_LN:-2.0}"
CUTOFF="${CUTOFF:-9.0}"
CPU_PER_TASK="${CPU_PER_TASK:-4}"

# Create organized subdirectories
mkdir -p "${OUTDIR}/inputs" "${OUTDIR}/outputs" "${OUTDIR}/restarts" \
         "${OUTDIR}/logs" "${OUTDIR}/mdinfo" "${OUTDIR}/trajectories" ./logs
[[ -f "$PDB_IN" ]] || { echo "[ERROR] PDB not found: $PDB_IN"; exit 1; }
echo "[OK] Parameters set"

############################
# ==== TEMPS & STEPS ==== #
############################
CURRENT_STEP=2
echo "===== STEP 2: Temperature and step calculation ====="
TEMPS=$(python3 - <<PY
lst = "${TEMPS_LIST}".strip()
if lst:
    ts = [float(x) for x in lst.split(",")]
else:
    n  = int("${REPLICAS}")
    t0 = float("${T_MIN}")
    t1 = float("${T_MAX}")
    r  = (t1/t0)**(1/(n-1)) if n>1 else 1.0
    ts = [t0*(r**i) for i in range(n)]
print(",".join(f"{t:.3f}" for t in ts))
PY
)
IFS=, read -r -a TARR <<< "$TEMPS"
[[ "${#TARR[@]}" -eq "$REPLICAS" ]] || { echo "[ERROR] Temperature count mismatch"; exit 1; }
echo "[INFO] Temps(K): $TEMPS"
echo "[OK] Temperature setup done"

NSTLIM_EXCH=$(python3 - <<PY
dt=float("${DT_PS}"); ps=float("${EXCHANGE_EVERY_PS}")
print(max(1,int(round(ps/dt))))
PY
)

NSTLIM_OUT=$(python3 - <<PY
N = int("${NSTLIM_EXCH}")
print(max(1,N*10))
PY
)

SALT_NUM=$(python3 - <<PY
M=float("${SALT_MOLAR}"); N=88.67
print(max(0,int(M*N)))
PY
)

############################
# ==== SYSTEM BUILD ==== #
############################
CURRENT_STEP=3
echo "===== STEP 3: Building system with tleap ====="
PDB4_DIR="${OUTDIR}/pdb4amber"
LEAP_DIR="${OUTDIR}/leap"
MIN_DIR="${OUTDIR}/min"

mkdir -p "$PDB4_DIR" "$LEAP_DIR" "$MIN_DIR"

LEAP_IN="${LEAP_DIR}/${OUTBASE}_leap.in"
# LEAP_LOG="${LEAP_DIR}/${OUTBASE}_tleap.log"
PRMTOP="${OUTDIR}/${OUTBASE}.prmtop"
INPCRD="${OUTDIR}/${OUTBASE}.inpcrd"
INIT_PDB="${LEAP_DIR}/${OUTBASE}_init.pdb"

if [[ "$USE_PDB4AMBER" -eq 1 ]] && command -v pdb4amber >/dev/null 2>&1; then
  echo "[INFO] Running pdb4amber..."
  pdb4amber -i "$PDB_IN" -o "${PDB4_DIR}/pdb4amber.pdb" -y -p 2> "${PDB4_DIR}/pdb4amber.log"
  PDB_IN="${PDB4_DIR}/pdb4amber.pdb"
fi

[[ "$BOX_SHAPE" == "Oct" ]] && SOLV_CMD="solvateOct prot ${BOX_TYPE} ${BOX_BUFFER}" || SOLV_CMD="solvateBox prot ${BOX_TYPE} ${BOX_BUFFER}"

cat > "$LEAP_IN" <<EOF
source leaprc.protein.ff14SB
source leaprc.water.tip3p
prot = loadpdb ${PDB_IN}
${SOLV_CMD}
EOF
if [[ "$NEUTRALIZE" -eq 1 ]]; then
cat >> "$LEAP_IN" <<'EOF'
addIonsRand prot Na+ 0
addIonsRand prot Cl- 0
EOF
fi
if python3 - <<<'import sys;sys.exit(0 if float("'$SALT_MOLAR'")>0 else 1)'; then
cat >> "$LEAP_IN" <<EOF
addIonsRand prot Na+ ${SALT_NUM} Cl- ${SALT_NUM} 4.0
EOF
fi
cat >> "$LEAP_IN" <<EOF

saveAmberParm prot ${PRMTOP} ${INPCRD}
savepdb      prot ${INIT_PDB}
quit
EOF


echo "[INFO] Running tleap..."
# Run tleap in LEAP_DIR so leap.log ends up there, not in CWD
cd "${LEAP_DIR}" && tleap -f "$LEAP_IN" | tee "${OUTBASE}_tleap.log" && cd - > /dev/null
[[ -f "$PRMTOP" && -f "$INPCRD" ]] || { echo "[ERROR] tleap failed"; exit 1; }
echo "[OK] tleap done"

############################
# ==== MINIMIZATION ==== #
############################
CURRENT_STEP=4
echo "===== STEP 4: Minimization ====="
MIN_IN="${MIN_DIR}/${OUTBASE}_min.in"
MIN_RST="${MIN_DIR}/${OUTBASE}_min.rst"
MIN_OUT="${MIN_DIR}/${OUTBASE}_min.out"
cat > "$MIN_IN" <<EOF
Minimize water/ions
 &cntrl
  imin=1, maxcyc=3000, ncyc=1500,
  cut=${CUTOFF}, ntb=1, igb=0,
  ntr=1, restraint_wt=5.0,
  restraintmask='!:WAT,Na+,Cl- & !@H='
 /
EOF
pmemd.cuda -O -i "$MIN_IN" -o "$MIN_OUT" -p "$PRMTOP" -c "$INPCRD" -r "$MIN_RST" -ref "$INPCRD" -inf "${MIN_DIR}/mdinfo_min"
if [ $? -ne 0 ] || [ ! -f "$MIN_RST" ]; then
    echo "[ERROR] Minimization failed"
    exit 1
fi
echo "[OK] Minimization complete"

###############################################
# ==== SCRATCH + GROUPFILES ==== #
###############################################
# Use basename of OUTDIR for scratch to avoid double paths
# OUTDIR_BASE="$(basename "$OUTDIR")"
# SCRATCH_DIR="$HOME/orcd/scratch/${OUTDIR_BASE}_${SLURM_JOB_ID:-$$}"
# mkdir -p "$SCRATCH_DIR"

###############################################
# ==== STEP 5: Auto density equilibration ==== #
###############################################
CURRENT_STEP=5
echo "===== STEP 5: Auto density equilibration (NPT, segmented) ====="

DENS_SEG_STEPS=${DENS_SEG_STEPS:-10000}
DENS_SEG_STEPS1=${DENS_SEG_STEPS1:-5000}   # gentler first segment
DENS_MIN_SEG=${DENS_MIN_SEG:-8}
DENS_MAX_SEG=${DENS_MAX_SEG:-20}
DENS_TOL_REL=${DENS_TOL_REL:-0.005}
DENS_CPUS=${DENS_CPUS:-10}          # CPU ranks for density (pmemd.MPI)
DENS_DT_SEG1=${DENS_DT_SEG1:-0.001}
DENS_DT_OTHERS=${DENS_DT_OTHERS:-${DT_PS}}
DENS_REST_WT1=${DENS_REST_WT1:-5.0}
DENS_REST_WT=${DENS_REST_WT:-2.0}

DENS_LAST_RST="$MIN_RST"
prev_vol_avg=""
seg=0

while (( seg < DENS_MAX_SEG )); do
  seg=$((seg+1))
  echo "[INFO] Density segment $seg starting from $DENS_LAST_RST"

  if (( seg == 1 )); then
    DENS_IREST=0
    DENS_NTX=1
    # Use a gentle Berendsen barostat for the first segment (isotropic for oct boxes)
    DENS_PRES_CTL="ntp=1, barostat=1"
    DENS_DT_USE=$DENS_DT_SEG1
    DENS_STEPS_USE=$DENS_SEG_STEPS1
    DENS_REST_WT_USE=$DENS_REST_WT1
  else
    DENS_IREST=1
    DENS_NTX=5
    # Switch back to MC barostat for subsequent segments (correct NPT ensemble)
    DENS_PRES_CTL="ntp=1, barostat=2"
    DENS_DT_USE=$DENS_DT_OTHERS
    DENS_STEPS_USE=$DENS_SEG_STEPS
    DENS_REST_WT_USE=$DENS_REST_WT
  fi

  DENSITY_IN="${OUTDIR}/inputs/${OUTBASE}_density_seg${seg}.in"
  cat > "$DENSITY_IN" <<EOF
Auto density equil (restrained NPT, segmented)
 &cntrl
  imin=0, irest=${DENS_IREST}, ntx=${DENS_NTX},
  ntb=2, ${DENS_PRES_CTL}, taup=1.0,
  cut=${CUTOFF},
  ntt=3, gamma_ln=${GAMMA_LN}, temp0=${T_MIN}, tempi=${T_MIN},
  ntc=2, ntf=2,
  dt=${DENS_DT_USE}, nstlim=${DENS_STEPS_USE},
  ntpr=1000, ntwx=1000, ntwr=5000,
  ntr=1, restraint_wt=${DENS_REST_WT_USE},
  restraintmask='!:WAT,Na+,Cl- & !@H=',
 /
EOF

  DENS_RST="${OUTDIR}/restarts/${OUTBASE}_density${seg}.rst"
  DENS_OUT="${OUTDIR}/outputs/${OUTBASE}_density${seg}.out"
  DENS_NC="${SCRATCH_DIR}/${OUTBASE}_density${seg}.nc"
  VOL_DAT="${OUTDIR}/outputs/${OUTBASE}_density${seg}_volume.dat"

  mpirun -np "$DENS_CPUS" pmemd.MPI -O \
    -i "$DENSITY_IN" \
    -o "$DENS_OUT" \
    -p "$PRMTOP" \
    -c "$DENS_LAST_RST" \
    -r "$DENS_RST" \
    -x "$DENS_NC" \
    -ref "$DENS_LAST_RST"

  cpptraj <<EOF2 > /dev/null 2>&1
parm $PRMTOP
trajin $DENS_NC
volume out $VOL_DAT
run
EOF2

  vol_avg=$(awk '{sum+=$2; n++} END{if(n>0) printf "%.6f", sum/n;}' "$VOL_DAT")
  if [[ -z "$vol_avg" ]]; then
    echo "[ERROR] Could not compute average volume for density segment $seg" >&2
    exit 1
  fi
  echo "[INFO] Segment $seg average volume: $vol_avg Å^3"

  if [[ -n "$prev_vol_avg" ]]; then
    rel_change=$(python3 - <<PY
prev = float("${prev_vol_avg}")
curr = float("${vol_avg}")
print(abs(curr - prev)/prev if prev != 0.0 else 1.0)
PY
)
    echo "[INFO] Segment $seg relative volume change vs previous: $rel_change"

    if (( seg >= DENS_MIN_SEG )); then
      above=$(python3 - <<PY
thr = float("${DENS_TOL_REL}")
rc  = float("${rel_change}")
print(1 if rc > thr else 0)
PY
)
      if [[ "$above" -eq 0 ]]; then
        echo "[INFO] Density converged (relative change <= ${DENS_TOL_REL}); stopping at segment $seg."
        break
      else
        echo "[INFO] Density not converged yet; continuing."
      fi
    fi
  fi

  prev_vol_avg="$vol_avg"
  DENS_LAST_RST="$DENS_RST"
done

if (( seg == DENS_MAX_SEG )); then
  echo "[WARNING] Reached DENS_MAX_SEG=${DENS_MAX_SEG} without satisfying tolerance; using last segment."
fi

if [[ -z "$DENS_RST" ]]; then
  echo "[ERROR] Auto density failed to generate a restart" >&2
  exit 1
fi

REF_DENS_RST="$DENS_RST"
echo "[INFO] Using $REF_DENS_RST as input for NVT equilibration."


############################
# ==== MPS START ==== #
############################
CURRENT_STEP=6
echo "===== STEP 6: Start MPS on all nodes ====="

export CUDA_DEVICE_ORDER=PCI_BUS_ID
export CUDA_VISIBLE_DEVICES="0,1,2,3"  # All 4 GPUs
export CUDA_MPS_PIPE_DIRECTORY="${SLURM_TMPDIR:-/tmp}/mps_${SLURM_JOB_ID}"
export CUDA_MPS_LOG_DIRECTORY="${SLURM_TMPDIR:-/tmp}/mps_${SLURM_JOB_ID}"

# Start MPS on node
mkdir -p ${CUDA_MPS_PIPE_DIRECTORY} ${CUDA_MPS_LOG_DIRECTORY}
nvidia-cuda-mps-control -d
sleep 2  # Give MPS time to start
echo "[OK] MPS started"


###############################################
# ==== STEP 7: Per-replica NVT equil (common box) ==== #
###############################################
CURRENT_STEP=7
echo "===== STEP 7: set up Per-replica NVT equilibration ====="

EQUI_PS="${EQUI_PS:-200.0}"            # ~200 ps typical
EQUI_SAVE_EVERY_PS="${EQUI_SAVE_EVERY_PS:-5.0}"   # trajectory sampling interval (ps); set 0 to disable
read -r EQUI_NSTLIM <<<"$(python3 - <<PY
dt=float("${DT_PS}")
def steps(ps): import math; return max(1,int(round(ps/dt)))
print(steps(${EQUI_PS}))
PY
)"

# How often to write equil trajectory.
# If EQUI_SAVE_EVERY_PS=0, disable trajectory output (ntwx=0).
read -r EQUI_NTWX <<<"$(python3 - <<PY
dt=float("${DT_PS}")
save_ps=float("${EQUI_SAVE_EVERY_PS}")
if save_ps <= 0:
  print(0)
else:
  print(max(1, int(round(save_ps/dt))))
PY
)"

# Write NVT equil inputs
for i in $(seq 0 $((REPLICAS-1))); do
  idx=$(printf "%03d" "$i"); T=${TARR[$i]}
  cat > "${OUTDIR}/inputs/${OUTBASE}_equi.in.${idx}" <<EOF
NVT equil @ ${T} K (common box)
 &cntrl
  imin=0, irest=0, ntx=1,        ! create velocities
  ntb=1, ntp=0,                  ! NVT -> box stays identical for all replicas
  cut=${CUTOFF},
  ntc=2, ntf=2,
  ntt=3, gamma_ln=${GAMMA_LN}, temp0=${T}, tempi=${T}, ig=-1,
  dt=${DT_PS}, nstlim=${EQUI_NSTLIM},
  ntpr=1000, ntwx=${EQUI_NTWX}, ntwr=1000,
  ioutfm=1, ntxo=2, iwrap=1
 /
EOF
done

# Groupfile for per-replica NVT equil (start from the SAME box)
EQUI_GRP="${OUTDIR}/${OUTBASE}_equi.group"
: > "$EQUI_GRP"
for i in $(seq 0 $((REPLICAS-1))); do
  idx=$(printf "%03d" "$i")
  echo "-O -i ${OUTDIR}/inputs/${OUTBASE}_equi.in.${idx} \
-o ${OUTDIR}/outputs/${OUTBASE}_equi.out.${idx} -p ${PRMTOP} \
-c ${REF_DENS_RST} \
-r ${OUTDIR}/restarts/${OUTBASE}_equi.rst.${idx} \
-x ${SCRATCH_DIR}/${OUTBASE}_equi.nc.${idx} \
-inf ${OUTDIR}/mdinfo/mdinfo_equi.${idx} -l ${OUTDIR}/logs/logfile_equi.${idx}" >> "$EQUI_GRP"
done

############################################
#### STEP 8: Run per-replica NVT equilibration (common box) ### #
#############################################

CURRENT_STEP=8
echo "===== STEP 8: Run NVT equil ====="
mpirun -np "$REPLICAS" "$ENGINE" -ng "$REPLICAS" -groupfile "$EQUI_GRP"
EQUI_EXIT=$?
if [ $EQUI_EXIT -ne 0 ]; then
    echo "[ERROR] EQUIL failed with exit code $EQUI_EXIT"
    exit 1
fi
echo "[OK] NVT equil done"

############################
# ==== STEP 9: REMD INPUTS (NVT) ==== #
############################
for i in $(seq 0 $((REPLICAS-1))); do
  idx=$(printf "%03d" "$i"); T=${TARR[$i]}
  cat > "${OUTDIR}/inputs/${OUTBASE}_remd.in.${idx}" <<EOF
T-REMD replica ${idx} at ${T} K
 &cntrl
  imin=0, irest=1, ntx=5,        ! continue from NVT equil
  ntb=1, ntp=0,                  ! NVT for T-REMD
  cut=${CUTOFF}, ntc=2, ntf=2,
  ntt=3, gamma_ln=${GAMMA_LN}, temp0=${T},
  dt=${DT_PS}, nstlim=${NSTLIM_EXCH},
  ntpr=${NSTLIM_OUT}, ntwx=${NSTLIM_OUT}, ntwr=${NSTLIM_OUT},
  ioutfm=1, ntxo=2, iwrap=1, nscm=500,
  numexchg=${TOTAL_EXCHANGES}
 /
EOF
done

GROUP="${OUTDIR}/${OUTBASE}_groupfile"
: > "$GROUP"
for i in $(seq 0 $((REPLICAS-1))); do
  idx=$(printf "%03d" "$i")
  echo "-O -rem 1 -remlog ${OUTDIR}/${OUTBASE}_rem.log \
-i ${OUTDIR}/inputs/${OUTBASE}_remd.in.${idx} -p ${PRMTOP} \
-c ${OUTDIR}/restarts/${OUTBASE}_equi.rst.${idx} \
-o ${OUTDIR}/outputs/${OUTBASE}_remd.out.${idx} \
-r ${OUTDIR}/restarts/${OUTBASE}_remd.rst.${idx} \
-x ${SCRATCH_DIR}/${OUTBASE}_remd.nc.${idx} \
-inf ${OUTDIR}/mdinfo/mdinfo_remd.${idx} -l ${OUTDIR}/logs/logfile_remd.${idx}" >> "$GROUP"
done

CURRENT_STEP=9
echo "===== STEP 9: Running REMD ====="
echo "[INFO] Starting REMD production run (this will take a while...)"
echo "[INFO] Progress: check ${OUTDIR}/${OUTBASE}_rem.log for exchange statistics"
mpirun -np "$REPLICAS" "$ENGINE" -ng "$REPLICAS" -groupfile "$GROUP"
REMD_EXIT=$?
if [ $REMD_EXIT -ne 0 ]; then
    echo "[ERROR] REMD failed with exit code $REMD_EXIT"
    exit 1
fi
echo "[OK] REMD done"
echo "[OK] REMD done"

############################
# ==== STOP MPS ==== #
############################
CURRENT_STEP=10
echo "===== STEP 10: Stop MPS ====="
echo quit | nvidia-cuda-mps-control >/dev/null 2>&1 || true
echo "[OK] MPS stopped"

############################
# ==== FINAL OUTPUT ==== #
############################
CURRENT_STEP=11
echo "===== STEP 11: Finalizing outputs ====="
ambpdb -p "$PRMTOP" -c "${OUTDIR}/restarts/${OUTBASE}_remd.rst.000" > "${OUTDIR}/${OUTBASE}_final_000.pdb" || true
echo "[OK] Final PDB exported"

# Create symlinks to trajectories in scratch
if [[ ${seg:-0} -ge 1 ]]; then
  for density_idx in $(seq 1 "$seg"); do
    dens_src="${SCRATCH_DIR}/${OUTBASE}_density${density_idx}.nc"
    if [[ -f "$dens_src" ]]; then
      ln -sf "$dens_src" "${OUTDIR}/trajectories/${OUTBASE}_density${density_idx}.nc"
    fi
  done
fi

for i in $(seq 0 $((REPLICAS-1))); do
  idx=$(printf "%03d" "$i")
  if [[ -f "${SCRATCH_DIR}/${OUTBASE}_equi.nc.${idx}" ]]; then
    ln -sf "${SCRATCH_DIR}/${OUTBASE}_equi.nc.${idx}" "${OUTDIR}/trajectories/${OUTBASE}_equi.nc.${idx}"
  fi
  ln -sf "${SCRATCH_DIR}/${OUTBASE}_remd.nc.${idx}" "${OUTDIR}/trajectories/${OUTBASE}_remd.nc.${idx}"
done

# Summary
echo ""
CURRENT_STEP=12
echo "=========================================="
echo "===== JOB COMPLETED $(date) ====="
echo "=========================================="
echo "[INFO] Output directory: $OUTDIR"
echo "[INFO] Trajectories (symlinks): ${OUTDIR}/trajectories/*.nc.*"
echo "[INFO] Actual trajectories preserved in: $SCRATCH_DIR"
echo "[INFO] Exchange log: ${OUTDIR}/${OUTBASE}_rem.log"
echo ""
echo "Next steps:"
echo "  1. Check exchange statistics: grep 'Replica' ${OUTDIR}/${OUTBASE}_rem.log | head -50"
echo "  2. Extract temperature 0 trajectory for analysis"
echo "  3. Calculate acceptance rates (should be 20-40%)"
echo "=========================================="


NS=$(awk "BEGIN {print ($TOTAL_EXCHANGES * $EXCHANGE_EVERY_PS) / 1000}")
ELAPSED_S="$(( $(date +%s) - CURRENT_TIME ))"
PARAMS_FILE="$OUTDIR/parameters.txt"
{
  echo "Elapsed seconds: $ELAPSED_S"
  echo "Job Name: $SLURM_JOB_NAME"
  echo "Job ID: $SLURM_JOB_ID"
  echo "Replicas: $REPLICAS"
  echo "nanoseconds of simulation per replica: $NS"
  echo "Base Temperature: $T_MIN K"
  echo "Max Temperature: $T_MAX K"
  echo "Box Type: $BOX_TYPE"
  echo "Box Buffer: $BOX_BUFFER"
  echo "Neutralize: $NEUTRALIZE"
  echo "Salt Molarity: $SALT_MOLAR"
  echo "Gamma LN: $GAMMA_LN"
  echo "Cutoff: $CUTOFF"
  echo "CPU per Task: $CPU_PER_TASK"
  echo "Total Exchanges: $TOTAL_EXCHANGES"
  echo "Exchange Every PS: $EXCHANGE_EVERY_PS"
  echo "Time Step (ps): $DT_PS"
  echo "Output Directory: $OUTDIR"
  echo "Scratch directory: $SCRATCH_DIR"
  echo "number of GPUs: ${SLURM_GPUS_ON_NODE:-1}"
} > "$PARAMS_FILE"


PARAMS_FILE2="$SCRATCH_DIR/parameters.txt"
{
  echo "Elapsed seconds: $ELAPSED_S"
  echo "Job Name: $SLURM_JOB_NAME"
  echo "Job ID: $SLURM_JOB_ID"
  echo "Replicas: $REPLICAS"
  echo "nanoseconds of simulation per replica: $NS"
  echo "Base Temperature: $T_MIN K"
  echo "Max Temperature: $T_MAX K"
  echo "Box Type: $BOX_TYPE"
  echo "Box Buffer: $BOX_BUFFER"
  echo "Neutralize: $NEUTRALIZE"
  echo "Salt Molarity: $SALT_MOLAR"
  echo "Gamma LN: $GAMMA_LN"
  echo "Cutoff: $CUTOFF"
  echo "CPU per Task: $CPU_PER_TASK"
  echo "Total Exchanges: $TOTAL_EXCHANGES"
  echo "Exchange Every PS: $EXCHANGE_EVERY_PS"
  echo "Time Step (ps): $DT_PS"
  echo "Output Directory: $OUTDIR"
  echo "Scratch directory: $SCRATCH_DIR"
  echo "number of GPUs: ${SLURM_GPUS_ON_NODE:-1}"
} > "$PARAMS_FILE2"



echo "===== STEP 12: Post-processing ====="
# bash "/home/jhalpin/orcd/pool/09-fragfold/XBL_PEP_DESIGN/development/30-testMD/amber_scripts/analyze_rep0.sh" "${OUTDIR}" "${OUTBASE}"
bash "${AMBER_SCRIPTS_DIR}/reconstruct_300K-parallel.sh" "${OUTDIR}" "${OUTBASE}" "${OUTDIR}/${OUTBASE}_rem.log" $T_MIN $REPLICAS
python "${AMBER_SCRIPTS_DIR}/rem_accept_cpptraj.py" "${OUTDIR}/${OUTBASE}_rem.log"
python "${AMBER_SCRIPTS_DIR}/gen_REMD_export.py" --export-folder="300K_export" --remd-folder="${OUTDIR}" --temperature=300 --overwrite
