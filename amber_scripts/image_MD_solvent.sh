#!/bin/bash
# Image a production MD trajectory retaining full solvent → *_prod_imaged.nc on scratch.
# Use this when you need to inspect water/ions (e.g. in VMD).
# For routine post-processing (stripped, analysed), use post_process_MD.sh instead.
#
# USAGE:
#   bash image_MD_solvent.sh <OUTDIR> <OUTBASE> <PRMTOP> <PROD_NC>
#
# ARGUMENTS:
#   OUTDIR    Output directory for this MD run (symlink will be placed here)
#   OUTBASE   System base name (file prefix)
#   PRMTOP    Path to topology file
#   PROD_NC   Path to raw production trajectory (usually on scratch)

set -e
set -o pipefail

if [[ $# -lt 4 ]]; then
  echo "Usage: $0 <OUTDIR> <OUTBASE> <PRMTOP> <PROD_NC>" >&2
  exit 1
fi

OUTDIR="$1"
OUTBASE="$2"
PRMTOP="$3"
PROD_NC="$4"

PROD_NC_REAL="$(readlink -f "$PROD_NC")"
SCRATCH_DIR="$(dirname "$PROD_NC_REAL")"
IMAGED_TRAJ_SCRATCH="${SCRATCH_DIR}/${OUTBASE}_prod_imaged.nc"
IMAGED_TRAJ_LINK="${OUTDIR}/${OUTBASE}_prod_imaged.nc"

echo "====================================================================="
echo "MD solvent imaging"
echo "System:      $OUTBASE"
echo "Raw traj:    $PROD_NC"
echo "Imaged traj: $IMAGED_TRAJ_SCRATCH (symlinked to $IMAGED_TRAJ_LINK)"
echo "====================================================================="

[[ -f "$PRMTOP" ]] || { echo "[ERROR] Topology not found: $PRMTOP" >&2; exit 1; }
[[ -f "$PROD_NC" ]] || { echo "[ERROR] Production trajectory not found: $PROD_NC" >&2; exit 1; }

cpptraj <<EOF
parm $PRMTOP
trajin $PROD_NC
autoimage
center mass origin
trajout $IMAGED_TRAJ_SCRATCH netcdf
run
quit
EOF

[[ -f "$IMAGED_TRAJ_SCRATCH" ]] || { echo "[ERROR] Imaged trajectory not created: $IMAGED_TRAJ_SCRATCH" >&2; exit 1; }
ln -sf "$IMAGED_TRAJ_SCRATCH" "$IMAGED_TRAJ_LINK"
echo "[OK] Imaged trajectory written to scratch: $IMAGED_TRAJ_SCRATCH"
echo "[OK] Symlink created: $IMAGED_TRAJ_LINK"
