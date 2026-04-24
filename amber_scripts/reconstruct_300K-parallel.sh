#!/bin/bash
# Reconstruct constant-temperature (e.g. 300K) trajectory from REMD ensemble
# Then perform the same analyses as `analyze_rep0.sh` but on the demultiplexed path.
#
# PARALLEL VERSION: uses cpptraj.MPI for frame-level MPI parallelism.
# Steps 2 (process), 4 (DSSP), and 5 (clustering) are merged into a single cpptraj
# invocation so DEMUX_RAW is read only once. Step 1 (demultiplex) also uses MPI.
#
# This script builds a demultiplexed trajectory at a chosen temperature using
# all replica trajectories plus the Amber `rem.log` (exchange log), then runs:
#   - Stripping, imaging, centering
#   - RMSD / RMSF / Radius of gyration / End-to-end distance
#   - Secondary structure (DSSP)
#   - Clustering (DBSCAN) on backbone atoms
#   - Subsampled visualization trajectory + reference PDB
#   - Plot generation (Python/matplotlib)
#   - Summary report
#
# USAGE:
#   reconstruct_300K-parallel.sh <WORKDIR> <BASE> <REMLOG> [TARGET_TEMP] [NPROCS]
# EXAMPLE:
#   reconstruct_300K-parallel.sh \
#       /home/jhalpin/orcd/pool/30-testMD/output/IL7-2-0_5ns-REMD-300-450K-60reps \
#       IL7-2 \
#       IL7-2_rem.log \
#       300.0 \
#       16
#
# ARGUMENTS:
#   WORKDIR      Path containing *.prmtop, ${BASE}_remd.nc.000, ..., rem.log
#   BASE         Base system name (prefix of prmtop & trajectories)
#   REMLOG       Exchange log file (e.g. IL7-2_rem.log) relative or absolute
#   TARGET_TEMP  Temperature to extract (float). Default 300.0
#   NPROCS       Number of MPI ranks. Default: $SLURM_NTASKS if set, else 8.
#                Override MPI launcher via env: MPI_CMD=srun
#
# OUTPUT DIRECTORY:
#   analysis_demux_<TEMP>K (e.g. analysis_demux_300K)
#
# ASSUMPTIONS / NOTES:
#   - Trajectories follow pattern: ${BASE}_remd.nc.000, .001, ...
#   - AmberTools/cpptraj.MPI is available and activated (falls back to cpptraj).
#   - Water & ions stripped for processed trajectories & clustering.
#   - Cluster parameters (epsilon/minpoints) are tuned for medium-size proteins; adjust if needed.
#   - DSSP (secstruct) is computed before stripping so all residues are visible.

# set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <WORKDIR> <BASE> [REMLOG] [TARGET_TEMP] [NPROCS]" 1>&2
  exit 1
fi

WORKDIR="$1"        # Directory with REMD output files
BASE="$2"           # System base name
DEFAULTREMLOG="$WORKDIR/${BASE}_rem.log"  # Default rem.log path
REMLOG="${3:-$DEFAULTREMLOG}"         # rem.log path (exchange log)
TARGET_TEMP="${4:-300.0}"   # Temperature to reconstruct
NPROCS="${5:-${SLURM_NTASKS:-8}}"    # MPI ranks
MPI_CMD="${MPI_CMD:-mpirun}"          # MPI launcher (set MPI_CMD=srun for SLURM)
SUBDIR="trajectories"  # Subdirectory with replica trajectories

# Convert WORKDIR to absolute path
WORKDIR="$(cd "$WORKDIR" && pwd)"

# Convert REMLOG to absolute path if not already absolute
if [[ "$REMLOG" != /* ]]; then
  REMLOG="$(cd "$(dirname "$REMLOG")" && pwd)/$(basename "$REMLOG")"
fi


if [[ ! -d "$WORKDIR" ]]; then
  echo "[ERROR] WORKDIR does not exist: $WORKDIR" 1>&2
  exit 1
fi
cd "$WORKDIR"

if [[ ! -f "$BASE.prmtop" ]]; then
  echo "[ERROR] Topology not found: $BASE.prmtop" 1>&2
  exit 1
fi

# rem.log is optional for demultiplex with remdtrajvalues; warn if missing but continue
if [[ ! -f "$REMLOG" ]]; then
  echo "[WARN] rem.log not found (optional): $REMLOG" 1>&2
fi

TEMP_LABEL="${TARGET_TEMP%%.*}K"
PRMTOP="${BASE}.prmtop"

echo "====================================================================="
echo "Reconstructing ${TARGET_TEMP} trajectory from REMD ensemble (PARALLEL)"
echo "System: $BASE"
echo "Workdir: $WORKDIR"
echo "Exchange log (optional): $REMLOG"
echo "Target Temp: $TARGET_TEMP"
echo "MPI ranks: $NPROCS  (launcher: $MPI_CMD)"
echo "====================================================================="
echo

# ---------------------------------------------------------------------------
# Locate cpptraj / cpptraj.MPI
# ---------------------------------------------------------------------------
if ! command -v cpptraj.MPI >/dev/null 2>&1 && ! command -v cpptraj >/dev/null 2>&1; then
  echo "[INFO] cpptraj not in PATH; attempting to source AmberTools environment"
  [ -n "${AMBER_SH:-}" ] && source "$AMBER_SH" 2>/dev/null || true
fi

if command -v cpptraj.MPI >/dev/null 2>&1; then
  echo "[INFO] cpptraj.MPI found: $(command -v cpptraj.MPI)"
  CPPTRAJ="${MPI_CMD} -n ${NPROCS} cpptraj.MPI"
elif command -v cpptraj >/dev/null 2>&1; then
  echo "[WARN] cpptraj.MPI not found — falling back to serial cpptraj (no MPI parallelism)"
  CPPTRAJ="cpptraj"
else
  echo "[WARN] Neither cpptraj.MPI nor cpptraj found in PATH. Proceeding may fail."
  CPPTRAJ="cpptraj"
fi
echo "[INFO] cpptraj command: $CPPTRAJ"
echo "[INFO] Proceeding to collect trajectories..."

# Collect replica trajectories
shopt -s nullglob
TRAJ_FILES=( "${SUBDIR}/${BASE}_remd.nc."[0-9][0-9][0-9] )
if [[ ${#TRAJ_FILES[@]} -eq 0 ]]; then
  echo "[ERROR] No replica trajectories found matching ${BASE}_remd.nc.[000-999]" 1>&2
  exit 1
fi
NREPS=${#TRAJ_FILES[@]}
echo "[INFO] Found $NREPS replica trajectories:"; printf '  %s\n' "${TRAJ_FILES[@]}"; echo

OUTDIR_MAIN="analysis_demux_${TEMP_LABEL}"
mkdir -p "$OUTDIR_MAIN"
OUTDIR="analysis_demux_${TEMP_LABEL}/cpptraj_out"
mkdir -p "$OUTDIR"

# ---------------------------------------------------------------------------
# Step 1: Demultiplex target temperature trajectory
# ---------------------------------------------------------------------------
echo "Step 1: Demultiplexing trajectory for ${TARGET_TEMP}K"
echo "------------------------------------------------------------------------"

DEMUX_IN="$OUTDIR_MAIN/reconstruct_${TEMP_LABEL}.cpptraj"
DEMUX_RAW="$OUTDIR_MAIN/${BASE}_${TEMP_LABEL}_raw.nc"

{
  echo "# Demultiplex ${TARGET_TEMP}K trajectory using embedded REMD metadata";
  echo "parm $PRMTOP";
  # cpptraj.MPI distributes replica reading across ranks automatically
  echo "trajin ${SUBDIR}/${BASE}_remd.nc.000 remdtraj remdtrajvalues ${TARGET_TEMP}";
  echo "trajout $DEMUX_RAW netcdf";
  echo "run";
  echo "quit";
} > "$DEMUX_IN"

echo "Running cpptraj (demultiplex) with: $CPPTRAJ ..."
if ! $CPPTRAJ -i "$DEMUX_IN" > "$OUTDIR_MAIN/reconstruct_${TEMP_LABEL}.log" 2>&1; then
  echo "[ERROR] cpptraj demultiplex failed; see $OUTDIR_MAIN/reconstruct_${TEMP_LABEL}.log" 1>&2
  exit 1
fi

if [[ ! -f "$DEMUX_RAW" ]]; then
  echo "[ERROR] Demultiplexed trajectory not created: $DEMUX_RAW" 1>&2
  exit 1
fi
NFRAMES=$(cpptraj -p "$PRMTOP" -y "$DEMUX_RAW" -tl 2>/dev/null | awk '/^Frames:/ {print $2}')
if [[ -z "$NFRAMES" ]]; then
  echo "[ERROR] Unable to determine frame count for $DEMUX_RAW" 1>&2
  exit 1
fi
echo "[OK] Demultiplexed trajectory: $DEMUX_RAW"; echo "[INFO] Frames: $NFRAMES"; echo

# ---------------------------------------------------------------------------
# Steps 2 + 4 + 5 (combined): Process, DSSP, Clustering
#
# Reading DEMUX_RAW once for all three saves significant I/O on large trajectories.
# DSSP (secstruct) is placed BEFORE the strip command so it sees all protein residues.
# Strip/autoimage/center then apply for all subsequent analyses and clustering.
# ---------------------------------------------------------------------------
PROC_TRAJ="$OUTDIR_MAIN/${BASE}_${TEMP_LABEL}_stripped.nc"
STRIP_TOP="$OUTDIR_MAIN/${BASE}_${TEMP_LABEL}_stripped.parm7"
COMBINED_IN="$OUTDIR_MAIN/analyze_${TEMP_LABEL}.cpptraj"

echo "Steps 2+4+5: Process / DSSP / Clustering (single MPI cpptraj run, $NPROCS ranks)"
echo "------------------------------------------------------------------------"

{
  echo "# Combined analysis pass — reads DEMUX_RAW once for all analyses";
  echo "parm $PRMTOP [full]";
  echo "trajin $DEMUX_RAW parm [full]";
  echo "";
  echo "# Step 4: Secondary structure (DSSP) — must run before strip";
  echo "secstruct out $OUTDIR/dssp_${TEMP_LABEL}.dat sumout $OUTDIR/dssp_summary_${TEMP_LABEL}.dat";
  echo "";
  echo "# Step 2: Strip solvent/ions, image, center";
  echo "strip :WAT,Na+,Cl-";
  echo "autoimage";
  echo "center mass origin";
  echo "trajout $PROC_TRAJ netcdf";
  echo "";
  echo "# Step 2 analyses (on stripped, imaged coords)";
  echo "rms ToFirst @CA,C,N,O first out $OUTDIR/rmsd_${TEMP_LABEL}.dat mass";
  echo "atomicfluct out $OUTDIR/rmsf_${TEMP_LABEL}.dat @CA byres";
  echo "radgyr out $OUTDIR/rgyr_${TEMP_LABEL}.dat mass";
  echo "distance end2end :1@CA :*@CA out $OUTDIR/end_to_end_${TEMP_LABEL}.dat";
  echo "";
  echo "# Step 5: Clustering (DBSCAN on backbone)";
  echo "cluster C1 @CA,C,N,O dbscan minpoints 25 epsilon 2.0 \\";
  echo "  sieve 10 \\";
  echo "  summary $OUTDIR/cluster_summary_${TEMP_LABEL}.dat \\";
  echo "  info $OUTDIR/cluster_info_${TEMP_LABEL}.dat \\";
  echo "  cpopvtime $OUTDIR/cluster_cpop_${TEMP_LABEL}.agr \\";
  echo "  repout $OUTDIR/${BASE}_cluster_${TEMP_LABEL} repfmt pdb \\";
  echo "  singlerepout $OUTDIR/${BASE}_clusters_${TEMP_LABEL}.pdb singlerepfmt pdb";
  echo "";
  echo "# Write stripped topology";
  echo "parm $PRMTOP [topstrip]";
  echo "parmstrip :WAT,Na+,Cl- parm [topstrip]";
  echo "parmwrite out $STRIP_TOP parm [topstrip]";
  echo "run";
  echo "quit";
} > "$COMBINED_IN"

echo "Running cpptraj (process+DSSP+cluster) with: $CPPTRAJ ..."
$CPPTRAJ -i "$COMBINED_IN" > "$OUTDIR_MAIN/analyze_${TEMP_LABEL}.log" 2>&1
echo "[OK] Combined analysis complete (process / DSSP / clustering)"; echo

# ---------------------------------------------------------------------------
# Step 6: Plotting (Python)
# ---------------------------------------------------------------------------
echo "Step 6: Plotting analyses"
echo "------------------------------------------------------------------------"
cat > "$OUTDIR/plot_analysis_${TEMP_LABEL}.py" <<'PYTHON_EOF'
#!/usr/bin/env python3
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os
from glob import glob

def safe_load(fname):
    try:
        return np.loadtxt(fname, comments='#')
    except Exception as e:
        print(f"[WARN] Could not load {fname}: {e}")
        return None

def load_first(*names):
  for name in names:
    if name is None:
      continue
    if os.path.exists(name):
      data = safe_load(name)
      if data is not None:
        print(f"[INFO] Using data file: {name}")
        return data, name
  return None, None

print("Generating plots for demultiplexed trajectory...")
plt.style.use('seaborn-v0_8-darkgrid')

label = os.getenv('TEMP_LABEL', '300K')

# Candidate filenames for flexibility across variants
rmsd, rmsd_name = load_first(
  f'rmsd_{label}.dat', 'rmsd_300K.dat', 'rmsd.dat', 'rmsd_rep0.dat', 'rmsd_${TEMP_LABEL}.dat'
)
rmsf, rmsf_name = load_first(
  f'rmsf_{label}.dat', 'rmsf.dat', 'rmsf_rep0.dat', 'rmsf_${TEMP_LABEL}.dat'
)
rgyr, rgyr_name = load_first(
  f'rgyr_{label}.dat', 'rgyr.dat', 'rgyr_rep0.dat', 'rgyr_${TEMP_LABEL}.dat'
)
e2e, e2e_name = load_first(
  f'end_to_end_{label}.dat', 'end_to_end.dat', 'end_to_end_rep0.dat', 'end_to_end_${TEMP_LABEL}.dat'
)

if rmsd is not None:
    try:
        fig, axes = plt.subplots(2,2, figsize=(16,12))
        axes[0,0].plot(rmsd[:,0], rmsd[:,1], lw=0.5, color='steelblue')
        axes[0,0].set_xlabel('Frame'); axes[0,0].set_ylabel('RMSD (Å)'); axes[0,0].set_title('RMSD vs Frame')
        axes[0,1].hist(rmsd[:,1], bins=50, alpha=0.7, color='steelblue', edgecolor='black')
        axes[0,1].axvline(rmsd[:,1].mean(), color='red', ls='--', lw=2, label=f"Mean {rmsd[:,1].mean():.2f} Å")
        axes[0,1].legend(); axes[0,1].set_xlabel('RMSD (Å)'); axes[0,1].set_ylabel('Freq')
        win=50
        if len(rmsd[:,1]) >= win:
            avg = np.convolve(rmsd[:,1], np.ones(win)/win, mode='valid')
            axes[1,0].plot(rmsd[:,0], rmsd[:,1], lw=0.3, alpha=0.4)
            axes[1,0].plot(rmsd[win-1:,0], avg, lw=2, color='darkblue', label=f'{win}-frame avg')
            axes[1,0].legend()
        axes[1,0].set_xlabel('Frame'); axes[1,0].set_ylabel('RMSD (Å)'); axes[1,0].set_title('RMSD (Running Avg)')
        sr = np.sort(rmsd[:,1]); cum = np.arange(1,len(sr)+1)/len(sr)
        axes[1,1].plot(sr, cum*100, lw=2, color='steelblue')
        axes[1,1].axvline(np.median(rmsd[:,1]), color='red', ls='--', lw=2, label=f"Median {np.median(rmsd[:,1]):.2f} Å")
        axes[1,1].set_xlabel('RMSD (Å)'); axes[1,1].set_ylabel('Cumulative %'); axes[1,1].legend()
        plt.tight_layout(); plt.savefig(f'rmsd_analysis_{label}.png', dpi=300)
        plt.close()
        print(f'✓ rmsd_analysis_{label}.png')
    except Exception as e:
        print(f"[WARN] RMSD plot failed: {e}")

if rmsf is not None:
    try:
        plt.figure(figsize=(14,6))
        plt.plot(rmsf[:,0], rmsf[:,1], lw=1.2, marker='o', markersize=3, color='darkgreen')
        plt.fill_between(rmsf[:,0], 0, rmsf[:,1], alpha=0.25, color='lightgreen')
        plt.xlabel('Residue'); plt.ylabel('RMSF (Å)'); plt.title('Per-Residue RMSF')
        plt.tight_layout(); plt.savefig(f'rmsf_analysis_{label}.png', dpi=300); plt.close()
        print(f'✓ rmsf_analysis_{label}.png')
    except Exception as e:
        print(f"[WARN] RMSF plot failed: {e}")

if rgyr is not None:
    try:
        fig,axes=plt.subplots(1,2, figsize=(14,6))
        axes[0].plot(rgyr[:,0], rgyr[:,1], lw=0.5, color='orange')
        axes[0].axhline(rgyr[:,1].mean(), color='red', ls='--', lw=2, label=f"Mean {rgyr[:,1].mean():.2f} Å")
        axes[0].set_xlabel('Frame'); axes[0].set_ylabel('Rg (Å)'); axes[0].legend(); axes[0].set_title('Rg vs Frame')
        axes[1].hist(rgyr[:,1], bins=50, alpha=0.7, color='orange', edgecolor='black')
        axes[1].axvline(rgyr[:,1].mean(), color='red', ls='--', lw=2)
        axes[1].set_xlabel('Rg (Å)'); axes[1].set_ylabel('Freq'); axes[1].set_title('Rg Distribution')
        plt.tight_layout(); plt.savefig(f'rgyr_analysis_{label}.png', dpi=300); plt.close()
        print(f'✓ rgyr_analysis_{label}.png')
    except Exception as e:
        print(f"[WARN] Rg plot failed: {e}")

if e2e is not None:
    try:
        fig,axes=plt.subplots(1,2, figsize=(14,6))
        axes[0].plot(e2e[:,0], e2e[:,1], lw=0.5, color='purple')
        axes[0].axhline(e2e[:,1].mean(), color='red', ls='--', lw=2, label=f"Mean {e2e[:,1].mean():.2f} Å")
        axes[0].set_xlabel('Frame'); axes[0].set_ylabel('End-to-End (Å)'); axes[0].legend(); axes[0].set_title('End-to-End vs Frame')
        axes[1].hist(e2e[:,1], bins=50, alpha=0.7, color='purple', edgecolor='black')
        axes[1].axvline(e2e[:,1].mean(), color='red', ls='--', lw=2)
        axes[1].set_xlabel('End-to-End (Å)'); axes[1].set_ylabel('Freq'); axes[1].set_title('End-to-End Distribution')
        plt.tight_layout(); plt.savefig(f'end_to_end_analysis_{label}.png', dpi=300); plt.close()
        print(f'✓ end_to_end_analysis_{label}.png')
    except Exception as e:
        print(f"[WARN] End-to-End plot failed: {e}")

if rmsd is not None and rgyr is not None and e2e is not None:
    try:
        fig,axes=plt.subplots(3,1, figsize=(16,12))
        axes[0].plot(rmsd[:,0], rmsd[:,1], lw=0.5, color='steelblue'); axes[0].set_ylabel('RMSD (Å)')
        axes[1].plot(rgyr[:,0], rgyr[:,1], lw=0.5, color='orange'); axes[1].set_ylabel('Rg (Å)')
        axes[2].plot(e2e[:,0], e2e[:,1], lw=0.5, color='purple'); axes[2].set_ylabel('End-to-End (Å)'); axes[2].set_xlabel('Frame')
        for ax in axes: ax.grid(alpha=0.3)
        plt.tight_layout(); plt.savefig(f'combined_overview_{label}.png', dpi=300); plt.close()
        print(f'✓ combined_overview_{label}.png')
    except Exception as e:
        print(f"[WARN] Combined overview failed: {e}")

print("\n✓ Plot generation complete.")
PYTHON_EOF

(
  cd "$OUTDIR"
  TEMP_LABEL="$TEMP_LABEL" python3 "plot_analysis_${TEMP_LABEL}.py" || echo "[WARN] Plot script had issues"
)
