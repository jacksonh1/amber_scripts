# Amber REMD/MD Pipeline

This readme was made by AI for an existing codebase I wrote. I did edit it.

I would recommend always doing the NPT REMD (as opposed to NVT) as the exchange rates are usually more even across the ladder even if they aren't as high

See [installation_scripts/amber_install_instructions.md](installation_scripts/amber_install_instructions.md) for instructions on installing AmberTools25 and PMEMD24.

Scripts for running Temperature Replica Exchange MD (T-REMD) and standard MD with Amber/pmemd on a SLURM cluster.

## Repository structure

```
.
├── site_config.sh        ← Edit this for your cluster (once)
├── amber_scripts/        ← Pipeline engine — do not edit
│   ├── REMD-NVT-singlenode.sbatch
│   ├── REMD-NPT-singlenode.sbatch
│   ├── REMD-NVT-resume.sbatch
│   ├── REMD-NPT-resume.sbatch
│   ├── MD.sbatch
│   ├── reconstruct_300K-parallel.sh
│   ├── reconstruct_300K.sh
│   ├── rem_accept_cpptraj.py
│   └── gen_REMD_export.py
└── example/              ← Example job for reference
    ├── input_pdbs/       ← Example input PDB files
    └── submit_jobs/      ← Example submission scripts (copy and edit for your job)
        ├── submit_REMD-NVT.sh
        ├── submit_REMD-NPT.sh
        ├── submit_MD.sh
        ├── submit_REMD-NVT-resume.sh
        ├── submit_REMD-NPT-resume.sh
        └── post_REMD_analysis.sh
```

## Pipeline overview

### NVT T-REMD

Runs replicas at constant volume. Simpler setup, no box size equilibration needed.

1. Build system with tleap (solvate, add ions)
2. Minimization
3. NPT density equilibration (convergence-checked: repeats until box volume stabilizes)
4. Per-replica NVT equilibration (each replica heated to its own target temperature from the same box)
5. T-REMD production (all replicas in parallel, exchanging temperatures every `EXCHANGE_EVERY_PS`)
6. Post-processing: reconstruct 300 K trajectory, compute acceptance rates

### NPT T-REMD

Runs replicas at constant pressure. Each replica equilibrates its own box volume before
REMD, which matters when replicas span a wide temperature range (boxes at 450 K are
meaningfully larger than at 300 K).

1. Build system with tleap (solvate, add ions)
2. Minimization
3. Per-replica NPT density equilibration (convergence-checked: repeats until box volume
   stabilizes, up to a configurable maximum number of segments)
4. Per-replica NPT equilibration (each replica starts from its own equilibrated box)
5. T-REMD production (all replicas in parallel, exchanging temperatures)
6. Post-processing: reconstruct 300 K trajectory, compute acceptance rates

**Use NPT when** your system is sensitive to pressure (e.g. membrane systems, large
conformational changes) or when you want a more physically rigorous ensemble. Use NVT
for speed or when the volume difference across the temperature range is small.

## Setup (do this once per cluster)

See [installation_scripts/amber_install_instructions.md](installation_scripts/amber_install_instructions.md) for instructions on installing AmberTools25 and PMEMD24.

The `ambertools25` conda environment includes all required Python packages (numpy, matplotlib, MDAnalysis). To create it from the provided spec:

```bash
conda env create -f installation_scripts/environment.yml
```

**1. Edit `site_config.sh`** in the repo root:

| Variable | Description |
|----------|-------------|
| `AMBER_SH` | Path to `amber.sh` from your pmemd installation |
| `CONDA_ENV` | Conda environment name (default: `ambertools25`) |
| `SCRATCH_ROOT` | Fast scratch storage root (e.g. local SSD or parallel fs). Make sure you have a lot of space free here, as the large trajectory files will be saved to this folder |
| `CUDA_MODULE` | CUDA module name on your cluster |
| `OPENMPI_MODULE` | OpenMPI module name on your cluster |

**2. Edit the `#SBATCH` headers** in the sbatch scripts you'll use (partition, GPU type):

```bash
#SBATCH -p your_partition
#SBATCH --gres=gpu:your_gpu_type:4
```

**3. Set `AMBER_SCRIPTS_DIR`** in any submit script you copy:

```bash
export AMBER_SCRIPTS_DIR="/path/to/this/repo/amber_scripts"
```

## Running a job

Copy a submit script from `example/submit_jobs/`, edit the parameters at the top, then run it:

```bash
cp example/submit_jobs/submit_REMD-NVT.sh my_run.sh
# edit my_run.sh: set AMBER_SCRIPTS_DIR, PDB_IN, OUTPUT_ROOT, and job parameters
bash my_run.sh
```

### NVT T-REMD (fresh start)

```bash
# submit_REMD-NVT.sh — key parameters:
REPLICAS=48
TOTAL_EXCHANGES=500000    # number of exchange attempts
T_MAX=450                 # max temperature (K); min is always 300 K
EXCHANGE_EVERY_PS=1.0     # exchange frequency

PDB_IN="/path/to/your/system.pdb"
OUTPUT_ROOT="/path/to/outputs"
# OUTDIR is auto-named: ${OUTPUT_ROOT}/${OUTBASE}-${NS}ns-REMD-300-${T_MAX}K-...
```

### NPT T-REMD (fresh start)

Same parameters as NVT. Use `submit_REMD-NPT.sh`.

### Standard MD

```bash
# submit_MD.sh — key parameters:
TEMP0=300      # simulation temperature (K)
PROD_NS=100    # production length (ns)
PDB_IN="/path/to/your/system.pdb"
OUTPUT_ROOT="/path/to/outputs"
```

### Resuming after a timeout or cluster maintenance

If a job is cancelled or times out, resume from where it left off — the script
auto-detects how many exchanges completed and continues from the last checkpoint:

```bash
cp example/submit_jobs/submit_REMD-NVT-resume.sh my_resume.sh
# edit: set AMBER_SCRIPTS_DIR, OUTBASE, OUTDIR (must match the original run), TOTAL_EXCHANGES
bash my_resume.sh
```

Safe to re-run multiple times. Each restart creates a new trajectory segment; segments
are concatenated automatically at the end.

There will probably be a smaller number of frames exported in the final trajectories. This is because
it writes a frame every 10 exchanges (default) and was probably interrupted in between 2 write time points. It
likely doesn't matter. You can back calculate the time corresponding to each frame if you need to (maybe from the rem log)
but you probably don't need that info

## Output directory layout

```
${OUTDIR}/
├── *.prmtop, *.inpcrd        — topology and coordinates
├── *_rem.log                 — exchange statistics
├── parameters.txt            — job parameters record
├── inputs/                   — Amber input files (equil, REMD)
├── outputs/                  — Amber output files
├── restarts/                 — checkpoint files
├── trajectories/             — symlinks to scratch trajectories
├── logs/                     — per-replica log files
├── mdinfo/                   — mdinfo files
├── density_equil/            — NPT only: density equilibration files
├── leap/                     — tleap input and log
├── min/                      — minimization files
├── analysis_demux_300K/      — demultiplexed 300 K trajectory files
└── 300K_export/              — exported 300 K trajectory files (xtc, for easy analysis)
```

There are some other files that are exported that might be useful for debugging or analysis (e.g. `rem_accept.dat` with exchange acceptance data, `parameters.txt` with a record of the job parameters, etc). The key outputs for analysis are the trajectories in `analysis_demux_300K/` and `300K_export/`.

Trajectory data lives on scratch (`SCRATCH_ROOT`) and is symlinked into `trajectories/`.
Copy or move the scratch directory before it is cleaned if you need the raw NetCDF files.

## Tips

- **Acceptance rates**: should be 20–40%. Check the `rem_accept.dat` file.
- **Scratch space**: each replica generates ~2 GB/µs. Ensure `SCRATCH_ROOT` has enough space.
- **Multiple systems**: copy/paste the system block in a submit script (or loop) to submit
  several systems as different jobs in one script call.
- **nodes**: The way that amber is installed in the install scripts, it only supports running a job on 1 node,
not splitting across multiple nodes. Splitting across multiple nodes might be slow because information will have
to transfer at each exchange step but if you know better go for it.
