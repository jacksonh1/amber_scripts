#!/bin/bash
# Post-process a production MD trajectory.
# 1. Fix periodic boundary and strip solvent (cpptraj) → *_prod_imaged_stripped.nc on scratch
# 2. Strip solvent, align, compute RMSD/Rg/RMSF via process_productionMD.py
#
# For a full-solvent imaged trajectory (rare, e.g. to inspect water), use image_MD_solvent.sh instead.
#
# USAGE:
#   bash post_process_MD.sh <OUTDIR> <OUTBASE> <PRMTOP> <PROD_NC> <AMBER_SCRIPTS_DIR>
#
# ARGUMENTS:
#   OUTDIR             Output directory for this MD run
#   OUTBASE            System base name (file prefix)
#   PRMTOP             Path to topology file
#   PROD_NC            Path to raw production trajectory (usually on scratch)
#   AMBER_SCRIPTS_DIR  Path to amber_scripts directory (for process_productionMD.py)

set -e
set -o pipefail

if [[ $# -lt 5 ]]; then
  echo "Usage: $0 <OUTDIR> <OUTBASE> <PRMTOP> <PROD_NC> <AMBER_SCRIPTS_DIR>" >&2
  exit 1
fi

OUTDIR="$1"
OUTBASE="$2"
PRMTOP="$3"
PROD_NC="$4"
AMBER_SCRIPTS_DIR="$5"

# Store processed trajectory alongside the raw trajectory on scratch (symlink into OUTDIR)
PROD_NC_REAL="$(readlink -f "$PROD_NC")"
SCRATCH_DIR="$(dirname "$PROD_NC_REAL")"
PROCESSED_TRAJ_SCRATCH="${SCRATCH_DIR}/${OUTBASE}_prod_imaged_stripped.nc"
PROCESSED_TRAJ_LINK="${OUTDIR}/${OUTBASE}_prod_imaged_stripped.nc"

echo "====================================================================="
echo "MD Post-processing"
echo "System:       $OUTBASE"
echo "Output dir:   $OUTDIR"
echo "Raw traj:     $PROD_NC"
echo "Processed:    $PROCESSED_TRAJ_SCRATCH (symlinked to $PROCESSED_TRAJ_LINK)"
echo "====================================================================="

[[ -f "$PRMTOP" ]] || { echo "[ERROR] Topology not found: $PRMTOP" >&2; exit 1; }
[[ -f "$PROD_NC" ]] || { echo "[ERROR] Production trajectory not found: $PROD_NC" >&2; exit 1; }

# ---------------------------------------------------------------------------
# Step 1: Image, center, and strip solvent in one cpptraj pass
# ---------------------------------------------------------------------------
echo ""
echo "Step 1: Imaging and stripping trajectory (cpptraj)..."

cpptraj <<EOF
parm $PRMTOP
trajin $PROD_NC
autoimage
center mass origin
strip :WAT,Na+,Cl-
trajout $PROCESSED_TRAJ_SCRATCH netcdf
run
quit
EOF

[[ -f "$PROCESSED_TRAJ_SCRATCH" ]] || { echo "[ERROR] Processed trajectory not created: $PROCESSED_TRAJ_SCRATCH" >&2; exit 1; }
ln -sf "$PROCESSED_TRAJ_SCRATCH" "$PROCESSED_TRAJ_LINK"
echo "[OK] Processed trajectory written to scratch: $PROCESSED_TRAJ_SCRATCH"
echo "[OK] Symlink created: $PROCESSED_TRAJ_LINK"

# ---------------------------------------------------------------------------
# Step 2: Align, compute RMSD/Rg/RMSF via Python
# ---------------------------------------------------------------------------
echo ""
echo "Step 2: Python analysis (align, RMSD/Rg/RMSF)..."
python "${AMBER_SCRIPTS_DIR}/process_productionMD.py" --md_folder "$OUTDIR"
echo "[OK] Python analysis done"

echo ""
echo "===== Post-processing complete ====="
