# Amber REMD/MD Pipeline

Scripts for running Temperature Replica Exchange MD (T-REMD) and standard MD with Amber/pmemd on a SLURM cluster.

See [installation_scripts/amber_install_instructions.md](installation_scripts/amber_install_instructions.md) for instructions on installing AmberTools25 and PMEMD24.

See the [Recommendations](#recommendations) section for how I would use these tools.

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
conformational changes) or when you want a more physically rigorous ensemble.

### NVT T-REMD

Runs replicas at constant volume. Simpler setup, no box size equilibration needed.

1. Build system with tleap (solvate, add ions)
2. Minimization
3. NPT density equilibration (convergence-checked: repeats until box volume stabilizes)
4. Per-replica NVT equilibration (each replica heated to its own target temperature from the same box)
5. T-REMD production (all replicas in parallel, exchanging temperatures every `EXCHANGE_EVERY_PS`)
6. Post-processing: reconstruct 300 K trajectory, compute acceptance rates


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
cp example/submit_jobs/submit_REMD-NPT.sh my_run.sh
# edit my_run.sh: set AMBER_SCRIPTS_DIR, PDB_IN, OUTPUT_ROOT, and job parameters
bash my_run.sh
```

### NPT T-REMD (fresh start)

```bash
# submit_REMD-NPT.sh — key parameters:
REPLICAS=48
TOTAL_EXCHANGES=1000      # number of exchange attempts if EXCHANGE_EVERY_PS=1.0, then 1000 exchanges will be 1 ns
T_MAX=450                 # max temperature (K)
T_MIN=300                 # min temperature (K); defaults to 300
EXCHANGE_EVERY_PS=1.0     # exchange frequency

PDB_IN="/path/to/your/system.pdb"
OUTPUT_ROOT="/path/to/outputs"
# OUTDIR is auto-named: ${OUTPUT_ROOT}/${OUTBASE}-${NS}ns-REMD-300-${T_MAX}K-...
```

### NVT T-REMD (fresh start)

Same parameters as NPT. Use `submit_REMD-NVT.sh`.

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
cp example/submit_jobs/submit_REMD-NPT-resume.sh my_resume.sh
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


## Recommendations


- **Acceptance rates**: should be 20–40%. Check the `rem_accept.dat` file.
- **Scratch space**: each replica generates ~10 GB/µs. Ensure `SCRATCH_ROOT` has enough space.
- **Multiple systems**: copy/paste the system block in a submit script (or loop) to submit
  several systems as different jobs in one script call.


**Other things to look for:**
- number of round trips - each coordinate set should have multiple trips up and down the temperature ladder
- dwell time - ideally each coordinate set spends an equal amount of time at each temperature
- watch out for periodic boundary artifacts - when a protein interacts with its neighboring images across the boundary. I'm working on a check for this. If you notice it, you need to increase the BOX_BUFFER


I would recommend always doing the NPT REMD (as opposed to NVT) as the exchange rates are usually more even across the ladder even if they aren't as high

REMD is usually limited by exchange frequency. You need to get around 20-40% exchange between each replicate, or else trajectories get stuck and the results are not thermodynamically valid. If the trajectories aren't mixing it also defeats the purpose of REMD in the first place because you aren't actually getting the enhanced sampling. If the exchange is too high, it means you're just wasting compute (I think).

The exchange probability is directly related to the energy difference between 2 neighboring replicas. For systems with more atoms, you need closer temperature spacing to get good exchange. So the number of replicates you need depends on the system. Unfortunately, you need 1 cpu per replica with the amber installation I have, so there is a limit on the number of replicas you can use. My amber installation also doesn't allow you to split one job across multiple nodes. You might not want this anyway as the nodes would have to transfer information with each other at every exchange.

I would recommend running a short pilot simulation for your system (maybe ~1 ns) and see what the exchange rates look like. I usually try to span the temperature range of 300-450K, which requires ~48 replicas for my mini-protein systems. This just fits on the 3219-3620 nodes.

Here's parameters I would use for a pilot simulation:
```bash
REPLICAS=48
TOTAL_EXCHANGES=1000
T_MAX=450
EXCHANGE_EVERY_PS=1.0
```
Simulation time (ns) = ($TOTAL_EXCHANGES * $EXCHANGE_EVERY_PS) / 1000

I usually set EXCHANGE_EVERY_PS to 1ps but I've seen people use 0.5 - 2 ps.

If the exchange rates are good, I'd do a longer production simulation. I usually aim for 40-100ns for this, but it depends on the protein and timescale of the dynamics you're looking for, which vary wildly

If you can't get good exchange you can lower the max temperature. There are other ways to address this issue that I haven't gotten to work. Let me know if you figure out a solution. Here's some ideas:
- Compile amber with nvhpc or some tool for splitting jobs across nodes so you can increase the replicates to more than 48. You'll have to write new scripts and it might be slow
- Change from T-REMD to REST2. I've tried getting REST2 to work (called REAF in amber). It worked at one time, then engaging changed some modules or available software and I cannot get it to work anymore. Might be unrelated to software changes but I don't know. Message me if you want that script. I get segfaults and no clear error message. It might be worth trying to patch amber with PLUMED if that's an option. I have moved over to gromacs instead for this reason.


**Replica vs. coordinate set:** Amber uses two indexing schemes. A *coordinate set* carries a fixed set of molecular coordinates. its temperature changes as exchanges are accepted. Raw trajectory files (`_remd.nc.000`, etc.) are therefore indexed by coordinate set and follow one structure as it travels up and down the temperature ladder — they are NOT at a fixed temperature. The demux step in post-processing reconstructs a fixed-temperature trajectory by collecting frames from whichever replica happened to be at 300 K at each exchange step, which is what ends up in `analysis_demux_300K/` and `300K_export/`. For most analyses you want the demuxed fixed-temperature trajectories, not the raw coordinate set trajectories. Amber often refers to a replica as a fixed temperature, in contrast to a coordinate set, which can be confusing.




<!-- - **nodes**: The way that amber is installed in the install scripts, it only supports running a job on 1 node, not splitting across multiple nodes. Splitting across multiple nodes might be slow because information will have to transfer at each exchange step, but if you know better go for it and let me know if it works -->
