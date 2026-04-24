#!/usr/bin/env python3
import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Tuple, List

HEADER = """
Compute pairwise REMD acceptance using cpptraj
- Reads an Amber rem.log
- Runs cpptraj remlog analysis to produce acceptance stats
- Prints a clean pairwise table and optionally plots it
""".strip()


def find_cpptraj() -> str:
    exe = shutil.which("cpptraj")
    if exe:
        return exe
    # Fallback to path from environment (set by site_config.sh via CPPTRAJ_PATH)
    fallback = os.environ.get("CPPTRAJ_PATH", "cpptraj")
    if os.path.isabs(fallback) and os.path.exists(fallback):
        return fallback
    raise FileNotFoundError(
        "cpptraj not found in PATH; activate ambertools25 (source activate ambertools25)"
    )


def build_cpptraj_input(rem_log: str, accept_out: str, stats_out: str, edata_out: str, crdidx_out: str, repidx_out: str, reptime_out: str) -> str:
    return f"""
readdata {rem_log} as remlog name rem1
remlog rem1 out {crdidx_out} crdidx 
remlog rem1 out {repidx_out} repidx
remlog rem1 reptime {reptime_out} acceptout {accept_out} stats statsout {stats_out}
run
quit
""".strip()
    # edata edataout {edata_out}
# remlog rem1 edata edataout {edata_out}
    #  


def run_cpptraj(rem_log: str, workdir: str, verbose: bool = False) -> Tuple[str, str]:
    """Run cpptraj to generate acceptance/stats files.

    If cpptraj succeeds but files are missing, raise with captured output.
    Never silently suppress exceptions.
    """
    cpptraj = find_cpptraj()
    accept_out = os.path.join(workdir, "rem_accept.dat")
    stats_out = os.path.join(workdir, "rem_stats.dat")
    edata_out = os.path.join(workdir, "rem_edata.dat")
    crdidx_out = os.path.join(workdir, "rem_crdidx.dat")
    repidx_out = os.path.join(workdir, "rem_repidx.dat")
    reptime_out = os.path.join(workdir, "rem_reptime.dat")
    inp = build_cpptraj_input(rem_log, accept_out, stats_out, edata_out, crdidx_out, repidx_out, reptime_out)
    cpptraj_file = Path(rem_log).parent / "cpptraj_remlog.in"
    # with tempfile.NamedTemporaryFile(
    #     "w", delete=False, dir=workdir, prefix="cpptraj_remlog_", suffix=".in"
    # ) as f:
    with open(cpptraj_file, "w") as f:
        f.write(inp + "\n")
        inpath = str(cpptraj_file)
    proc = subprocess.run(
        [cpptraj, "-i", inpath], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    # Attempt to remove the temporary file; warn if it fails.
    # if os.path.exists(inpath):
    #     try:
    #         os.unlink(inpath)
    #     except OSError as e:
    #         sys.stderr.write(f"[WARN] Could not remove temporary file {inpath}: {e}\n")
    if verbose:
        sys.stderr.write("[cpptraj stdout]\n" + proc.stdout + "\n")
        sys.stderr.write("[cpptraj stderr]\n" + proc.stderr + "\n")
    if proc.returncode != 0:
        raise RuntimeError(
            f"cpptraj failed (exit {proc.returncode}) while processing rem.log.\nSTDERR:\n{proc.stderr}\nSTDOUT:\n{proc.stdout}"
        )
    # Validate output files
    if not os.path.exists(accept_out):
        raise FileNotFoundError(
            "cpptraj completed without creating acceptance file.\n"
            f"Expected file: {accept_out}\n"
            "Possible causes: log is malformed, insufficient permissions, or cpptraj encountered an internal parse error.\n"
            f"cpptraj stdout snippet (first 40 lines):\n{os.linesep.join(proc.stdout.splitlines()[:40])}\n"
            f"cpptraj stderr snippet (first 40 lines):\n{os.linesep.join(proc.stderr.splitlines()[:40])}"
        )
    return accept_out, stats_out


def parse_accept(accept_path: str) -> Tuple[List[float], List[float]]:
    up: List[float] = []
    down: List[float] = []
    with open(accept_path, "r") as f:
        for line in f:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 3:
                # Format: Replica  %UP  %DOWN (replica is 1-based)
                try:
                    _ = int(parts[0])
                    up.append(float(parts[1]))
                    down.append(float(parts[2]))
                except ValueError:
                    continue
    return up, down  # lists of length N, 1-based indices at i-1


def parse_first_block_temps(rem_log: str) -> List[float]:
    temps: List[float] = []
    exch_re = re.compile(r"^#\s*exchange", re.IGNORECASE)
    row_re = re.compile(
        r"^\s*(\d+)\s+([\-\d\.]+)\s+([\-\d\.]+)\s+([\-\d\.]+)\s+([\-\d\.]+)\s+([\-\d\.]+)"
    )
    in_block = False
    with open(rem_log, "r", errors="ignore") as f:
        for line in f:
            if exch_re.match(line):
                if in_block:  # already read first block
                    break
                in_block = True
                continue
            if not in_block:
                continue
            m = row_re.match(line)
            if m:
                # groups: idx, vel, T, Eptot, Temp0, NewTemp0
                temp0 = float(m.group(5))
                temps.append(temp0)
            elif line.strip().startswith("#"):
                # Next block header reached; stop.
                break
    return temps


def print_table(
    temps: List[float], up: List[float], down: List[float], csv_path: str | None = None
) -> List[tuple]:
    n = min(len(temps), len(up), len(down))
    rows = []
    print(f"{'Pair':>6} {'Temp1':>10} {'Temp2':>10} {'UP(i)':>10} {'DOWN(i+1)':>12} {'Mean%':>8}")
    for i in range(1, n):  # pairs 1..n-1
        t1 = temps[i - 1]
        t2 = temps[i]
        up_i = up[i - 1]
        down_ip1 = down[i]
        mean_acc = (up_i + down_ip1) / 2.0
        print(
            f"{i:>3}-{i+1:<3} {t1:10.2f} {t2:10.2f} {up_i:10.2f} {down_ip1:12.2f} {mean_acc:8.2f}"
        )
        rows.append((i, i + 1, t1, t2, up_i, down_ip1, mean_acc))
    if csv_path:
        import csv

        with open(csv_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["Pair", "Temp1", "Temp2", "UP_i", "DOWN_i+1", "Mean"])
            for r in rows:
                w.writerow(
                    [
                        f"{r[0]}-{r[1]}",
                        f"{r[2]:.2f}",
                        f"{r[3]:.2f}",
                        f"{r[4]:.2f}",
                        f"{r[5]:.2f}",
                        f"{r[6]:.2f}",
                    ]
                )
        print(f"\nWrote CSV: {csv_path}")
    return rows


def plot_table(temps: List[float], up: List[float], down: List[float], out_png: str) -> str:
    try:
        import matplotlib.pyplot as plt

        # Use a non-interactive backend if DISPLAY is not set
        if not os.environ.get("DISPLAY"):
            try:
                import matplotlib

                matplotlib.use("Agg")
            except Exception as e:
                sys.stderr.write(f"[WARN] Could not set matplotlib backend to Agg: {e}\n")
    except ImportError:
        print("matplotlib not installed. Skipping plot.")
        return ""
    n = min(len(temps), len(up), len(down))
    pairs = list(range(1, n))
    mean = [(up[i - 1] + down[i]) / 2.0 for i in pairs]
    midtemps = [(temps[i - 1] + temps[i]) / 2.0 for i in pairs]
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    ax1.plot(pairs, mean, "o-", label="Mean acceptance (%)")
    ax1.fill_between(pairs, 15, 40, color="green", alpha=0.1, label="Target 15-40%")
    ax1.set_xlabel("Pair (i,i+1)")
    ax1.set_ylabel("Acceptance %")
    ax1.legend()
    ax1.grid(alpha=0.3)
    ax2.plot(midtemps, mean, "o-", color="orange")
    ax2.fill_between(midtemps, 15, 40, color="green", alpha=0.1)
    ax2.set_xlabel("Avg Temp (K)")
    ax2.set_ylabel("Acceptance %")
    ax2.grid(alpha=0.3)
    plt.tight_layout()
    fig.savefig(out_png, dpi=300)
    print(f"Saved plot: {out_png}")
    # Try to show if interactive; ignore errors
    try:
        plt.show()
    except Exception as e:
        sys.stderr.write(f"[WARN] Could not display plot: {e}\n")
    plt.close(fig)
    return out_png


def main(remlog_file: Path, verbose: bool = True, output_dir: Path | None = None):

    rem_log = remlog_file.resolve()
    if not rem_log.exists():
        sys.exit(f"rem.log not found: {rem_log}")

    if output_dir is None:
        output_dir = rem_log.parent
    accept_path, stats_path = run_cpptraj(rem_log, output_dir, verbose=verbose)
    up, down = parse_accept(accept_path)
    temps = parse_first_block_temps(rem_log)
    csv_path = output_dir / f"{rem_log.stem}_acceptance.csv"
    rows = print_table(temps, up, down, csv_path=str(csv_path))
    out_png = output_dir / f"{rem_log.stem}_acceptance.png"
    plot_table(temps, up, down, str(out_png))
    # Simple parity check: mean acceptance from cpptraj vs our row means
    if rows:
        cpptraj_mean = sum(((r[4] + r[5]) / 2.0) for r in rows) / len(rows)
        print(f"\n[INFO] Mean pair acceptance (cpptraj-derived): {cpptraj_mean:.2f}%")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=HEADER)
    ap.add_argument("remlog", help="Path to IL7-2_rem.log (AMBER REMD log)")
    ap.add_argument(
        "--verbose",
        action="store_true",
        help="Print cpptraj stdout/stderr for debugging",
    )
    args = ap.parse_args()
    main(Path(args.remlog), verbose=args.verbose)
