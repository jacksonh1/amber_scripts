#!/bin/bash
#SBATCH --job-name=MD_production
# NOTE: partition and GPU below are defaults only. run_scripts/submit_MD.sh overrides
# them from site_config.sh via --partition and --gres CLI flags at submission time.
#SBATCH -p pi_keating
#SBATCH --gres=gpu:l40s:1
#SBATCH --ntasks=1
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH -t 4-00:00:00
#SBATCH -o logs/MD_%x_%j.out
#SBATCH -e logs/MD_%x_%j.err

set -e
set -o pipefail

: "${AMBER_SCRIPTS_DIR:?AMBER_SCRIPTS_DIR is not set — submit via submit_jobs/submit_*.sh or: export AMBER_SCRIPTS_DIR=/path/to/amber_scripts}"
source "${AMBER_SCRIPTS_DIR}/../site_config.sh"

# Track current step for error handling
CURRENT_STEP=0

source activate "$CONDA_ENV"
module purge
module load gcc
module load "$CUDA_MODULE"
module load "$OPENMPI_MODULE"
source "$AMBER_SH"
echo "[OK] Modules loaded"

# Setup scratch directory
SCRATCH_TS="$(date +%Y%m%d-%H%M%S)"
SCRATCH_DIR="${SCRATCH_ROOT}/${SLURM_JOB_ID:-$$}_${SCRATCH_TS}"
mkdir -p "$SCRATCH_DIR"

# Cleanup function for errors
cleanup_scratch_on_error() {
  local rc=$?
  if [ -n "$SCRATCH_DIR" ] && [ -d "$SCRATCH_DIR" ]; then
    echo "[ERROR] Job failed at step $CURRENT_STEP (rc=$rc); cleaning scratch: $SCRATCH_DIR"
    rm -rf "$SCRATCH_DIR"
  fi
  exit $rc
}
trap 'cleanup_scratch_on_error' ERR

# Check available disk space in scratch
SCRATCH_AVAIL=$(df -BG "$SCRATCH_DIR" | awk 'NR==2 {print $4}' | sed 's/G//')
if [ "$SCRATCH_AVAIL" -lt 50 ]; then
    echo "[WARNING] Low disk space in scratch: ${SCRATCH_AVAIL}GB available"
    echo "[WARNING] Continuing anyway, but monitor disk usage!"
fi

echo "===== JOB START $(date) ====="

####usage of this script: sbatch --export=ALL,PARAM1=VALUE1,... MD.sbatch
############################
# ==== USER PARAMETERS ====#
############################
CURRENT_STEP=1
echo "===== STEP 1: Setting parameters ====="
PDB_IN="${PDB_IN:-IL7-2-V3-cys.pdb}"       # Input PDB file
OUTBASE="${OUTBASE:-IL7-2-V3-cys}"            # Output prefix
OUTDIR="${OUTDIR:-out_${OUTBASE}}"    # Output directory
USE_PDB4AMBER="${USE_PDB4AMBER:-1}"   # 1 = clean with pdb4amber, 0 = skip

# Solvation and ion setup
BOX_TYPE="${BOX_TYPE:-TIP3PBOX}"      # Water box type
BOX_SHAPE="${BOX_SHAPE:-Oct}"         # Oct / Tri / Cubic
BOX_BUFFER="${BOX_BUFFER:-10.0}"      # Solvent buffer distance (Å)
NEUTRALIZE="${NEUTRALIZE:-1}"         # 1 = neutralize charge
SALT_MOLAR="${SALT_MOLAR:-0.15}"         # Salt molarity (M), 0 = none

# MD parameters
DT_PS="${DT_PS:-0.002}"               # Time step (ps) = 2 fs
PROD_NS="${PROD_NS:-1}"               # Production length (ns)
TRAJ_EVERY_PS="${TRAJ_EVERY_PS:-5}"   # Trajectory write frequency (ps/frame)
ENERGY_EVERY_PS="${ENERGY_EVERY_PS:-5}" # Energy print frequency (ps)
RESTART_EVERY_PS="${RESTART_EVERY_PS:-100}" # Restart write frequency (ps)
TEMP0="${TEMP0:-300.0}"               # Target temperature (K)
GAMMA_LN="${GAMMA_LN:-2.0}"           # Langevin collision frequency
CUTOFF="${CUTOFF:-9.0}"               # Nonbonded cutoff (Å)
######################################

# --- Dependency checks ---
need(){ command -v "$1" >/dev/null 2>&1 || { echo "[ERROR] $1 not found"; exit 1; }; }
need tleap; need pmemd.cuda; need ambpdb
command -v pdb4amber >/dev/null 2>&1 && HAVE_PDB4=1 || HAVE_PDB4=0
[[ "$USE_PDB4AMBER" -eq 1 && $HAVE_PDB4 -eq 0 ]] && { echo "[WARN] pdb4amber not available, skipping."; USE_PDB4AMBER=0; }

# --- Setup directories and file names ---
mkdir -p "${OUTDIR}/inputs" "${OUTDIR}/outputs" "${OUTDIR}/restarts" logs
[[ -f "$PDB_IN" ]] || { echo "[ERROR] PDB not found: $PDB_IN"; exit 1; }
echo "[OK] Parameters set"

LEAP_IN="${OUTDIR}/${OUTBASE}_leap.in"
PRMTOP="${OUTDIR}/${OUTBASE}.prmtop"
INPCRD="${OUTDIR}/${OUTBASE}.inpcrd"
INIT_PDB="${OUTDIR}/${OUTBASE}_init.pdb"
FINAL_PDB="${OUTDIR}/${OUTBASE}_final.pdb"

MIN1_IN="${OUTDIR}/inputs/min1.in"; MIN2_IN="${OUTDIR}/inputs/min2.in"; HEAT_IN="${OUTDIR}/inputs/heat.in"
DENSITY_IN="${OUTDIR}/inputs/density.in"; EQ_IN="${OUTDIR}/inputs/eq.in"; PROD_IN="${OUTDIR}/inputs/prod.in"

MIN1_OUT="${OUTDIR}/outputs/${OUTBASE}_min1.out"; MIN1_RST="${OUTDIR}/restarts/${OUTBASE}_min1.rst"
MIN2_OUT="${OUTDIR}/outputs/${OUTBASE}_min2.out"; MIN2_RST="${OUTDIR}/restarts/${OUTBASE}_min2.rst"
HEAT_OUT="${OUTDIR}/outputs/${OUTBASE}_heat.out"; HEAT_RST="${OUTDIR}/restarts/${OUTBASE}_heat.rst"; HEAT_NC="${SCRATCH_DIR}/${OUTBASE}_heat.nc"
DENS_OUT="${OUTDIR}/outputs/${OUTBASE}_density.out"; DENS_RST="${OUTDIR}/restarts/${OUTBASE}_density.rst"; DENS_NC="${SCRATCH_DIR}/${OUTBASE}_density.nc"
EQ_OUT="${OUTDIR}/outputs/${OUTBASE}_eq.out"; EQ_RST="${OUTDIR}/restarts/${OUTBASE}_eq.rst"; EQ_NC="${SCRATCH_DIR}/${OUTBASE}_eq.nc"
PROD_OUT="${OUTDIR}/outputs/${OUTBASE}_prod.out"; PROD_RST="${OUTDIR}/restarts/${OUTBASE}_prod.rst"; PROD_NC="${SCRATCH_DIR}/${OUTBASE}_prod.nc"

echo "==> Output directory: $OUTDIR (prefix: $OUTBASE)"
echo "==> Input PDB: $PDB_IN"

############################
# ==== SYSTEM BUILD ==== #
############################
CURRENT_STEP=2
echo "===== STEP 2: Building system with tleap ====="

# --- Optional pdb4amber cleanup ---
if [[ "$USE_PDB4AMBER" -eq 1 ]] && command -v pdb4amber >/dev/null 2>&1; then
  echo "[INFO] Running pdb4amber cleanup..."
  pdb4amber -i "$PDB_IN" -o "${OUTDIR}/pdb4amber.pdb" -y -p 2> "${OUTDIR}/pdb4amber.log"
  PDB_IN="${OUTDIR}/pdb4amber.pdb"
  echo "    Cleaned PDB: $PDB_IN"
fi

SALT_NUM=$(python3 - <<PY
M=float("${SALT_MOLAR}"); N=88.67
print(max(0,int(M*N)))
PY
)

# --- Write MD input files ---
cat > "$MIN1_IN" <<'EOF'
Min1: heavy solute restrained, relax water/ions
 &cntrl
  imin=1, maxcyc=5000, ncyc=2500,
  cut=9.0, ntb=1, igb=0,
  ntr=1, restraint_wt=5.0,
  restraintmask='!:WAT,Na+,Cl- & !@H='
 /
EOF
cat > "$MIN2_IN" <<'EOF'
Min2: weaker restraint on solute
 &cntrl
  imin=1, maxcyc=10000, ncyc=5000,
  cut=9.0, ntb=1, igb=0,
  ntr=1, restraint_wt=1.0,
  restraintmask='!:WAT,Na+,Cl- & !@H='
 /
EOF
cat > "$HEAT_IN" <<EOF
Heat 0→${TEMP0}K (50 ps)
 &cntrl
  imin=0, irest=0, ntx=1, ntb=1, cut=${CUTOFF},
  ntc=2, ntf=2, ntt=3, gamma_ln=${GAMMA_LN},
  tempi=0.0, temp0=${TEMP0},
  nstlim=25000, dt=0.002,
  ntpr=500, ntwx=500, ntwr=5000
 /
EOF
cat > "$DENSITY_IN" <<EOF
Density equilibration (200 ps, restrained solute)
 &cntrl
  imin=0, irest=1, ntx=5,
  ntb=2, ntp=1, barostat=2, taup=1.0,
  cut=${CUTOFF},
  ntc=2, ntf=2,
  temp0=${TEMP0}, ntt=3, gamma_ln=${GAMMA_LN},
  nstlim=100000, dt=0.002,
  ntpr=1000, ntwx=1000, ntwr=10000,
  ntr=1, restraint_wt=0.5,
  restraintmask='!:WAT,Na+,Cl- & !@H='
 /
EOF
cat > "$EQ_IN" <<EOF
Equilibration NPT (500 ps)
 &cntrl
  imin=0, irest=1, ntx=5,
  ntb=2, ntp=1, barostat=2, taup=1.0,
  cut=${CUTOFF},
  ntc=2, ntf=2,
  temp0=${TEMP0}, ntt=3, gamma_ln=${GAMMA_LN},
  nstlim=250000, dt=0.002,
  ntpr=1000, ntwx=1000, ntwr=10000
 /
EOF

# --- Calculate production steps ---
NSTLIM=$(python3 - <<PY
dt=${DT_PS}; ns=float("${PROD_NS}")
print(max(1, int(round(ns*1000.0/dt))))
PY
)
ps2steps(){ python3 - "$1" "$DT_PS" <<'PY'
import sys; ps=float(sys.argv[1]); dt=float(sys.argv[2])
print(max(1, int(round(ps/dt))))
PY
}
NTWX=$(ps2steps "${TRAJ_EVERY_PS}")
NTPR=$(ps2steps "${ENERGY_EVERY_PS}")
NTWR=$(ps2steps "${RESTART_EVERY_PS}")

cat > "$PROD_IN" <<EOF
Production NPT (${PROD_NS} ns)
 &cntrl
  imin=0, irest=1, ntx=5,
  ntb=2, ntp=1, barostat=2, taup=1.0,
  cut=${CUTOFF},
  ntc=2, ntf=2,
  ntt=3, gamma_ln=${GAMMA_LN}, temp0=${TEMP0},
  nstlim=${NSTLIM}, dt=${DT_PS},
  ntpr=${NTPR}, ntwx=${NTWX}, ntwr=${NTWR},
  ioutfm=1, ntxo=2,
  iwrap=1, nscm=500
 /
EOF
echo "==> PROD_NS=${PROD_NS} ns → nstlim=${NSTLIM}"
echo "==> TRAJ every ${TRAJ_EVERY_PS} ps → ntwx=${NTWX}"
echo "==> ENERGY every ${ENERGY_EVERY_PS} ps → ntpr=${NTPR}"
echo "==> RESTART every ${RESTART_EVERY_PS} ps → ntwr=${NTWR}"

# --- Write LEaP input file ---
SOLV_LINE="solvate${BOX_SHAPE} prot ${BOX_TYPE} ${BOX_BUFFER}"
cat > "$LEAP_IN" <<EOF
source leaprc.protein.ff14SB
source leaprc.water.tip3p
prot = loadpdb ${PDB_IN}

# To manually create disulfide bonds (if pdb4amber didn't make CYX):
# bond prot.<resid1>.SG prot.<resid2>.SG

${SOLV_LINE}
EOF
if python3 - <<<'import sys;sys.exit(0 if float("'$SALT_MOLAR'")>0 else 1)'; then
cat >> "$LEAP_IN" <<EOF
addIonsRand prot Na+ ${SALT_NUM} Cl- ${SALT_NUM} 4.0
EOF
fi
if [[ "$NEUTRALIZE" -eq 1 ]]; then
cat >> "$LEAP_IN" <<'EOF'
addions prot Na+ 0
addions prot Cl- 0
EOF
fi
cat >> "$LEAP_IN" <<EOF

saveAmberParm prot ${PRMTOP} ${INPCRD}
savepdb      prot ${INIT_PDB}
quit
EOF

# --- Run tleap ---
echo "[INFO] Running tleap..."
cd "${OUTDIR}" && tleap -f "$LEAP_IN" | tee "${OUTBASE}_tleap.log" && cd - > /dev/null
[[ -f "$PRMTOP" && -f "$INPCRD" ]] || { echo "[ERROR] tleap failed"; exit 1; }
echo "[OK] tleap done"

# --- Run MD stages ---
############################
# ==== MINIMIZATION ==== #
############################
CURRENT_STEP=3
echo "===== STEP 3: Minimization ====="
echo "==> Min1 ..."
pmemd.cuda -O -i "$MIN1_IN" -o "$MIN1_OUT" -p "$PRMTOP" -c "$INPCRD"   -r "$MIN1_RST" -ref "$INPCRD" -inf "${SCRATCH_DIR}/mdinfo_min1"
echo "==> Min2 ..."
pmemd.cuda -O -i "$MIN2_IN" -o "$MIN2_OUT" -p "$PRMTOP" -c "$MIN1_RST" -r "$MIN2_RST" -ref "$INPCRD" -inf "${SCRATCH_DIR}/mdinfo_min2"
echo "[OK] Minimization complete"

############################
# ==== HEATING ==== #
############################
CURRENT_STEP=4
echo "===== STEP 4: Heating ====="
pmemd.cuda -O -i "$HEAT_IN" -o "$HEAT_OUT" -p "$PRMTOP" -c "$MIN2_RST" -r "$HEAT_RST" -x "$HEAT_NC" -inf "${SCRATCH_DIR}/mdinfo_heat"
echo "[OK] Heating complete"

############################
# ==== DENSITY ==== #
############################
CURRENT_STEP=5
echo "===== STEP 5: Density equilibration ====="
pmemd.cuda -O -i "$DENSITY_IN" -o "$DENS_OUT" -p "$PRMTOP" -c "$HEAT_RST" -r "$DENS_RST" -x "$DENS_NC" -ref "$HEAT_RST" -inf "${SCRATCH_DIR}/mdinfo_density"
echo "[OK] Density equilibration complete"

############################
# ==== EQUILIBRATION ==== #
############################
CURRENT_STEP=6
echo "===== STEP 6: Equilibration ====="
pmemd.cuda -O -i "$EQ_IN" -o "$EQ_OUT" -p "$PRMTOP" -c "$DENS_RST"  -r "$EQ_RST"  -x "$EQ_NC" -inf "${SCRATCH_DIR}/mdinfo_eq"
echo "[OK] Equilibration complete"

############################
# ==== PRODUCTION ==== #
############################
CURRENT_STEP=7
echo "===== STEP 7: Production ====="
if [[ -f "$PROD_RST" ]]; then
  START_FOR_PROD="$PROD_RST"
  echo "==> Production (RESUME) from existing restart: $PROD_RST"
else
  START_FOR_PROD="$EQ_RST"
  echo "==> Production (NEW) starting from: $EQ_RST"
fi

pmemd.cuda -O -i "$PROD_IN" -o "$PROD_OUT" -p "$PRMTOP" -c "$START_FOR_PROD" -r "$PROD_RST" -x "$PROD_NC" -inf "${SCRATCH_DIR}/mdinfo_prod"
echo "[OK] Production complete"

############################
# ==== FINAL OUTPUT ==== #
############################
CURRENT_STEP=8
echo "===== STEP 8: Finalizing outputs ====="
# --- Export final PDB ---
echo "==> Writing final PDB ..."
ambpdb -p "$PRMTOP" -c "$PROD_RST" > "$FINAL_PDB" || echo "[WARN] ambpdb failed"
echo "[OK] Final PDB exported"

# Create symlinks to trajectories in scratch
ln -sf "${SCRATCH_DIR}/${OUTBASE}_heat.nc" "${OUTDIR}/${OUTBASE}_heat.nc"
ln -sf "${SCRATCH_DIR}/${OUTBASE}_density.nc" "${OUTDIR}/${OUTBASE}_density.nc"
ln -sf "${SCRATCH_DIR}/${OUTBASE}_eq.nc" "${OUTDIR}/${OUTBASE}_eq.nc"
ln -sf "${SCRATCH_DIR}/${OUTBASE}_prod.nc" "${OUTDIR}/${OUTBASE}_prod.nc"

echo ""
echo "=========================================="
echo "===== JOB COMPLETED $(date) ====="
echo "=========================================="
echo "[INFO] Output directory: $OUTDIR"
echo "[INFO] Trajectories (symlinks): ${OUTDIR}/*.nc"
echo "[INFO] Actual trajectories preserved in: $SCRATCH_DIR"
echo ""
echo "Key outputs in $OUTDIR (prefix $OUTBASE):"
echo "  - ${OUTBASE}.prmtop / ${OUTBASE}.inpcrd / ${OUTBASE}_init.pdb"
echo "  - inputs/*: MD input files"
echo "  - outputs/*: MD output files (.out)"
echo "  - restarts/*: Restart files (.rst)"
echo "  - *.nc: Trajectory symlinks → $SCRATCH_DIR"
echo "  - ${OUTBASE}_final.pdb"
echo "=========================================="
