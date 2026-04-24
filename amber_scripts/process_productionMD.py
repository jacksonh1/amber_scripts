import numpy as np
from pathlib import Path
import argparse

import MDAnalysis as mda
from MDAnalysis.analysis import align, rms
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import tempfile
import os


def align_and_strip_trajectory(
    topology: Path,
    trajectory: Path,
    traj_output_file: Path,
    pdb_output_file: Path,
    trajout_slice: slice = slice(None, None, None),
    input_traj_format: str = "NCDF",
    output_rmsd: bool = False,
    output_Rg: bool = False,
    in_memory: bool = True,
) -> None:
    traj_output_file.parent.mkdir(parents=True, exist_ok=True)
    pdb_output_file.parent.mkdir(parents=True, exist_ok=True)
    universe = mda.Universe(topology, trajectory, format=input_traj_format)
    universe.trajectory[0]
    protein = universe.select_atoms("protein")
    if protein.n_atoms == 0:
        raise ValueError("Protein selection is empty; check the topology contents.")

    # Build reference on the first trajectory frame for self-alignment.
    protein.write(str(pdb_output_file))
    reference = mda.Universe(str(pdb_output_file))
    tmp_path = None
    if in_memory:
        aligner = align.AlignTraj(
            universe,
            reference,
            select="protein and backbone",
            in_memory=True,
            match_atoms=True,
            strict=True,
            weights="mass",
        ).run()
        universe.trajectory[0]  # reset to first frame after aligning
    else:
        tmp_file = tempfile.NamedTemporaryFile(suffix=".xtc", delete=False)
        tmp_file.close()
        tmp_path = Path(tmp_file.name)
        aligner = align.AlignTraj(
            universe,
            reference,
            select="protein and backbone",
            in_memory=False,
            match_atoms=True,
            strict=True,
            weights="mass",
            filename=str(tmp_path),
        ).run()
        universe = mda.Universe(str(topology), str(tmp_path))
        universe.trajectory[0]  # reset to first frame after aligning
        protein = universe.select_atoms("protein")
    with mda.Writer(str(traj_output_file), n_atoms=protein.n_atoms) as writer:
        for _ in universe.trajectory[trajout_slice]:
            writer.write(protein)
    if output_rmsd:
        rms_output = traj_output_file.parent / f"{traj_output_file.stem}_rmsd.dat"
        x = np.zeros((len(aligner.results.rmsd[trajout_slice]), 2))
        frames = np.arange(len(aligner.results.rmsd))
        x[:, 0] = frames[trajout_slice]
        x[:, 1] = aligner.results.rmsd[trajout_slice]
        np.savetxt(rms_output, x, delimiter=",")
        png_output = traj_output_file.parent / f"{traj_output_file.stem}_rmsd.png"
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(x[:, 0], x[:, 1])
        ax.set_xlabel("Frame")
        ax.set_ylabel("RMSD (Å) to first frame")
        fig.savefig(png_output)
        plt.close(fig)
    if output_Rg:
        rg_output = traj_output_file.parent / f"{traj_output_file.stem}_Rg.dat"
        rg_array = np.zeros((len(universe.trajectory[trajout_slice]), 2))
        frames = np.arange(len(universe.trajectory))
        rg_values = []
        for ts in universe.trajectory[trajout_slice]:
            rg = protein.radius_of_gyration()
            rg_values.append(rg)
        rg_array[:, 0] = frames[trajout_slice]
        rg_array[:, 1] = np.array(rg_values)
        np.savetxt(rg_output, rg_array, delimiter=",")
        png_output = traj_output_file.parent / f"{traj_output_file.stem}_Rg.png"
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(rg_array[:, 0], rg_array[:, 1])
        ax.set_xlabel("Frame")
        ax.set_ylabel("Radius of Gyration (Å)")
        fig.savefig(png_output)
        plt.close(fig)
    if tmp_path is not None and tmp_path.exists():
        os.remove(tmp_path)


def find_file(folder: Path, pattern: str) -> Path:
    found_files = list(folder.glob(pattern))
    if not found_files:
        raise FileNotFoundError(f"No files matching {pattern} found in {folder}")
    if len(found_files) > 1:
        raise ValueError(f"Multiple files matching {pattern} found in {folder}")
    return found_files[0]


def strip_trajectories(folder: Path, overwrite: bool = False):
    print(f"Processing folder: {folder}")
    output_dir = folder
    try:
        topology_input_file = find_file(folder, "*.prmtop")
        trajectory_input_file = find_file(folder, "*_prod.nc")
        trajectory_output_file = output_dir / trajectory_input_file.name.replace("_prod.nc", "_prod_stripped.xtc")
        pdb_output_file = output_dir / trajectory_input_file.name.replace("_prod.nc", "prod_f1_stripped.pdb")
    except FileNotFoundError as e:
        print(e)
        return
    if trajectory_output_file.exists() and not overwrite:
        print(f"Output file {trajectory_output_file} already exists. Skipping.")
        return
    align_and_strip_trajectory(
        topology=topology_input_file,
        trajectory=trajectory_input_file,
        traj_output_file=trajectory_output_file,
        pdb_output_file=pdb_output_file,
        output_rmsd=True,
        output_Rg=True,
        # trajout_slice=slice(-2000, None),
    )


def main(md_parent_folder: Path, overwrite: bool = False):
    for folder in md_parent_folder.iterdir():
        if folder.is_dir():
            strip_trajectories(folder, overwrite=overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process production MD trajectories.")
    parser.add_argument("--md_parent_folder", type=Path, help="Path to the folder containing the production MD trajectories.", default="./output_MD/")
    parser.add_argument("--overwrite", action="store_true", help="Whether to overwrite existing output files.")
    args = parser.parse_args()
    main(Path(args.md_parent_folder), overwrite=args.overwrite)
