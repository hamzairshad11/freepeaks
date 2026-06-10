from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def load_grid(path: Path):
    data = np.loadtxt(path, comments="#")
    n = int(round(np.sqrt(data.shape[0])))
    return data[:, 2].reshape(n, n), data[:, 3].reshape(n, n), data[:, 6].reshape(n, n)


def load_points(path: Path):
    if not path.exists() or path.stat().st_size == 0:
        return np.empty((0, 7))
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return data


def load_edges(path: Path):
    if not path.exists() or path.stat().st_size == 0:
        return np.empty((0, 5))
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return data


def label(suite_dir: Path):
    parts = suite_dir.name.split("_", 2)
    return f"S{parts[1]} {parts[2].replace('_', ' ')}" if len(parts) == 3 else suite_dir.name


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", default="Visualization/free_peaks_multiparty_samples")
    parser.add_argument("--result", default="Visualization/multiparty_nbn_cmaes")
    parser.add_argument("--output", default="Visualization/figures/free_peaks_multiparty_nbn_cmaes_overlay.png")
    args = parser.parse_args()

    sample_root = Path(args.samples)
    result_root = Path(args.result)
    suite_dirs = sorted(p for p in result_root.glob("suite_*") if re.match(r"suite_\d+_p\d+_", p.name))
    if not suite_dirs:
        raise RuntimeError(f"No solver suite directories found under {result_root}")

    fig, axes = plt.subplots(2, 4, figsize=(16, 8), constrained_layout=True)
    axes = axes.ravel()
    for ax, res_dir in zip(axes, suite_dirs):
        sample_dir = sample_root / res_dir.name
        X, Y, Z = load_grid(sample_dir / "grid.txt")
        cf = ax.contourf(X, Y, Z, levels=np.linspace(0, 90, 28), cmap="viridis")
        ax.contour(X, Y, Z, levels=10, colors="black", linewidths=0.2, alpha=0.25)

        for fname, color, alpha in [("party0_nbn_edges.txt", "tab:red", 0.18), ("party1_nbn_edges.txt", "tab:blue", 0.18)]:
            edges = load_edges(res_dir / fname)
            for e in edges[:250]:
                ax.plot([e[0], e[2]], [e[1], e[3]], color=color, linewidth=0.35, alpha=alpha)

        coord = load_points(res_dir / "coordination_population.txt")
        if coord.size:
            ax.scatter(coord[:, 2], coord[:, 3], s=7, c="white", edgecolors="none", alpha=0.35, label="coord")
        archive = load_points(res_dir / "archive_consensus_optima.txt")
        if archive.size:
            ax.scatter(archive[:, 2], archive[:, 3], s=32, c="red", marker="*", edgecolors="white", linewidths=0.4, label="archive")

        ax.set_title(label(res_dir), fontsize=10)
        ax.set_xlabel("x0")
        ax.set_ylabel("x1")
        ax.set_aspect("equal")
        fig.colorbar(cf, ax=ax, fraction=0.046, pad=0.03)

    fig.suptitle("MultiParty NBN-CMAES: NBN edges, coordination population, and consensus archive", fontsize=14)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=220)
    print(out)


if __name__ == "__main__":
    main()
