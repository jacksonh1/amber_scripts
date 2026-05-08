#!/usr/bin/env python3
"""
Cluster a stripped/aligned MD trajectory by Cα RMSD.

Loads *_prod_f1_stripped.pdb (topology) + *_prod_stripped.xtc (trajectory),
runs k-means or DBSCAN on flattened Cα coordinates, and writes cluster
representative PDB files + a per-frame cluster assignment CSV.

USAGE:
    python cluster_MD.py --md_folder /path/to/run [options]

OPTIONS:
    --md_folder     Path to the MD output directory
    --method        Clustering method: kmeans (default) or dbscan
    --n_clusters    Number of clusters for k-means (default: 5)
    --epsilon       DBSCAN epsilon in Å (default: 2.0)
    --min_samples   DBSCAN min_samples (default: 10)
    --stride        Use every Nth frame (default: 1 = all frames)
    --selection     MDAnalysis atom selection for clustering (default: 'name CA')
    --outdir        Output directory for results (default: <md_folder>/clustering)
"""

import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

import MDAnalysis as mda
from MDAnalysis.analysis import align
from sklearn.cluster import KMeans, DBSCAN
from sklearn.preprocessing import StandardScaler


def restore_chain_ids(pdb_path: Path, pdb4amber_path: Path) -> None:
    """Rewrite chain IDs in pdb_path using pdb4amber.pdb as reference."""
    ref = mda.Universe(str(pdb4amber_path))
    resid_to_chain = {}
    for res in ref.residues:
        chain = res.segment.segid.strip()
        if not chain:
            chain = 'A'
        resid_to_chain[res.resid] = chain
    lines = pdb_path.read_text().splitlines()
    new_lines = []
    for line in lines:
        if line.startswith(('ATOM', 'HETATM')) and len(line) >= 26:
            try:
                resid = int(line[22:26].strip())
                chain = resid_to_chain.get(resid, 'X')
                line = line[:21] + chain + line[22:]
            except ValueError:
                pass
        new_lines.append(line)
    pdb_path.write_text('\n'.join(new_lines) + '\n')


def find_file(folder: Path, pattern: str) -> Path:
    found = list(folder.glob(pattern))
    if not found:
        raise FileNotFoundError(f"No file matching '{pattern}' in {folder}")
    if len(found) > 1:
        raise ValueError(f"Multiple files matching '{pattern}' in {folder}: {found}")
    return found[0]


def load_trajectory(md_folder: Path, stride: int) -> mda.Universe:
    topology = find_file(md_folder, "*_prod_f1_stripped.pdb")
    trajectory = find_file(md_folder, "*_prod_stripped.xtc")
    print(f"[INFO] Topology:   {topology.name}")
    print(f"[INFO] Trajectory: {trajectory.name}")
    u = mda.Universe(str(topology), str(trajectory))
    print(f"[INFO] Total frames: {len(u.trajectory)}")
    if stride > 1:
        print(f"[INFO] Stride: {stride} → using {len(u.trajectory[::stride])} frames")
    return u


def extract_coordinates(u: mda.Universe, selection: str, stride: int) -> np.ndarray:
    """Return (n_frames, n_atoms * 3) array of selected atom coordinates."""
    atoms = u.select_atoms(selection)
    if atoms.n_atoms == 0:
        raise ValueError(f"Selection '{selection}' matched no atoms.")
    print(f"[INFO] Clustering on {atoms.n_atoms} atoms (selection: '{selection}')")

    coords = []
    for _ in u.trajectory[::stride]:
        coords.append(atoms.positions.flatten())
    return np.array(coords)


def run_kmeans(coords: np.ndarray, n_clusters: int, random_state: int = 42):
    print(f"[INFO] Running k-means with k={n_clusters}...")
    km = KMeans(n_clusters=n_clusters, random_state=random_state, n_init=10)
    labels = km.fit_predict(coords)
    inertias = None
    return labels, km.cluster_centers_


def run_dbscan(coords: np.ndarray, epsilon: float, min_samples: int):
    print(f"[INFO] Running DBSCAN with epsilon={epsilon} Å, min_samples={min_samples}...")
    # Scale coordinates so epsilon is meaningful in Å units
    db = DBSCAN(eps=epsilon, min_samples=min_samples, metric='euclidean', n_jobs=-1)
    labels = db.fit_predict(coords)
    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise = (labels == -1).sum()
    print(f"[INFO] Found {n_clusters} clusters, {n_noise} noise frames ({100*n_noise/len(labels):.1f}%)")
    return labels, n_clusters


def find_representative_frame(coords: np.ndarray, labels: np.ndarray,
                               cluster_id: int) -> int:
    """Return the frame index closest to the centroid of the given cluster."""
    mask = labels == cluster_id
    cluster_coords = coords[mask]
    centroid = cluster_coords.mean(axis=0)
    dists = np.linalg.norm(cluster_coords - centroid, axis=1)
    local_idx = np.argmin(dists)
    # Map back to global frame index
    global_indices = np.where(mask)[0]
    return global_indices[local_idx]


def write_cluster_pdbs(u: mda.Universe, selection: str, stride: int,
                       labels: np.ndarray, coords: np.ndarray,
                       outdir: Path, md_folder: Path) -> None:
    all_protein = u.select_atoms("protein")
    cluster_ids = sorted(set(labels))
    pdb4amber_path = md_folder / "pdb4amber.pdb"
    for cid in cluster_ids:
        label = f"noise" if cid == -1 else f"cluster{cid:02d}"
        rep_frame = find_representative_frame(coords, labels, cid)
        actual_frame = rep_frame * stride
        u.trajectory[actual_frame]
        out_pdb = outdir / f"{label}_rep.pdb"
        all_protein.write(str(out_pdb))
        if pdb4amber_path.exists():
            restore_chain_ids(out_pdb, pdb4amber_path)
        count = (labels == cid).sum()
        print(f"[INFO] {label}: {count} frames ({100*count/len(labels):.1f}%), "
              f"rep = frame {actual_frame} → {out_pdb.name}")


def plot_cluster_populations(labels: np.ndarray, outdir: Path, method: str) -> None:
    cluster_ids, counts = np.unique(labels, return_counts=True)
    fig, ax = plt.subplots(figsize=(max(6, len(cluster_ids)), 5))
    colors = ['lightcoral' if c == -1 else 'steelblue' for c in cluster_ids]
    x_labels = [f'noise' if c == -1 else f'C{c}' for c in cluster_ids]
    ax.bar(x_labels, counts, color=colors, edgecolor='black')
    ax.set_xlabel('Cluster')
    ax.set_ylabel('Frames')
    ax.set_title(f'Cluster populations ({method})')
    for i, (xl, cnt) in enumerate(zip(x_labels, counts)):
        ax.text(i, cnt + 0.5, f'{100*cnt/len(labels):.1f}%',
                ha='center', va='bottom', fontsize=9)
    plt.tight_layout()
    fig.savefig(outdir / 'cluster_populations.png', dpi=200)
    plt.close(fig)


def plot_cluster_timeseries(labels: np.ndarray, stride: int, outdir: Path) -> None:
    frames = np.arange(len(labels)) * stride
    fig, ax = plt.subplots(figsize=(14, 4))
    ax.scatter(frames, labels, s=2, c=labels, cmap='tab10')
    ax.set_xlabel('Frame')
    ax.set_ylabel('Cluster')
    ax.set_title('Cluster assignment over time')
    plt.tight_layout()
    fig.savefig(outdir / 'cluster_timeseries.png', dpi=200)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Cluster a stripped MD trajectory.")
    parser.add_argument('--md_folder', type=Path, required=True)
    parser.add_argument('--method', choices=['kmeans', 'dbscan'], default='kmeans')
    parser.add_argument('--n_clusters', type=int, default=5,
                        help='Number of clusters (k-means only)')
    parser.add_argument('--epsilon', type=float, default=2.0,
                        help='DBSCAN epsilon in Å')
    parser.add_argument('--min_samples', type=int, default=10,
                        help='DBSCAN min_samples')
    parser.add_argument('--stride', type=int, default=1,
                        help='Use every Nth frame (default: 1 = all)')
    parser.add_argument('--selection', type=str, default='name CA',
                        help="MDAnalysis atom selection (default: 'name CA')")
    parser.add_argument('--outdir', type=Path, default=None)
    args = parser.parse_args()

    outdir = args.outdir or (args.md_folder / 'clustering')
    outdir.mkdir(parents=True, exist_ok=True)

    u = load_trajectory(args.md_folder, args.stride)
    coords = extract_coordinates(u, args.selection, args.stride)

    if args.method == 'kmeans':
        labels, centers = run_kmeans(coords, args.n_clusters)
    else:
        labels, n_clusters = run_dbscan(coords, args.epsilon, args.min_samples)

    # Save per-frame assignments
    frames = np.arange(len(labels)) * args.stride
    csv_path = outdir / 'cluster_assignments.csv'
    np.savetxt(csv_path, np.column_stack([frames, labels]),
               delimiter=',', header='frame,cluster', fmt='%d', comments='')
    print(f"[OK] Cluster assignments saved: {csv_path.name}")

    write_cluster_pdbs(u, args.selection, args.stride, labels, coords, outdir, args.md_folder)
    plot_cluster_populations(labels, outdir, args.method)
    plot_cluster_timeseries(labels, args.stride, outdir)

    print(f"\n[OK] Clustering complete. Results in: {outdir}")


if __name__ == '__main__':
    main()
