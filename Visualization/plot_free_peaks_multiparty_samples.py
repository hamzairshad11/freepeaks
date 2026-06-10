from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def load_grid(path: Path):
    data = np.loadtxt(path, comments="#")
    n = int(round(np.sqrt(data.shape[0])))
    if n * n != data.shape[0]:
        raise ValueError(f"{path} does not contain a square grid: {data.shape[0]} rows")
    x = data[:, 2].reshape(n, n)
    y = data[:, 3].reshape(n, n)
    f0 = data[:, 4].reshape(n, n)
    f1 = data[:, 5].reshape(n, n)
    consensus = data[:, 6].reshape(n, n)
    extent = [x.min(), x.max(), y.min(), y.max()]
    return extent, f0, f1, consensus


def load_shared(path: Path):
    if not path.exists() or path.stat().st_size == 0:
        return np.empty((0, 7))
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return data


def suite_label(suite_dir: Path) -> str:
    name = suite_dir.name
    parts = name.split("_", 2)
    if len(parts) == 3 and parts[0] == "suite":
        return f"S{parts[1]} {parts[2].replace('_', ' ')}"
    return name.replace("_", " ")


def suite_id(suite_dir: Path) -> int:
    match = re.match(r"suite_(\d+)_", suite_dir.name)
    return int(match.group(1)) if match else 10**9


def add_shared_points(ax, shared: np.ndarray):
    if shared.size == 0:
        return
    ax.scatter(
        shared[:, 2],
        shared[:, 3],
        marker="*",
        s=70,
        c="white",
        edgecolors="black",
        linewidths=0.7,
        zorder=3,
    )


def draw_image(ax, z, extent, title, cmap):
    xs = np.linspace(extent[0], extent[1], z.shape[1])
    ys = np.linspace(extent[2], extent[3], z.shape[0])
    X, Y = np.meshgrid(xs, ys)
    levels = np.linspace(0, max(90.0, float(np.nanmax(z))), 29)
    img = ax.contourf(X, Y, z, levels=levels, cmap=cmap)
    ax.contour(X, Y, z, levels=10, colors="black", linewidths=0.25, alpha=0.35)
    ax.set_title(title, fontsize=9)
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    ax.set_xlabel("x0")
    ax.set_ylabel("x1")
    ax.set_aspect("equal")
    return img


def plot_consensus(suite_dirs, output_dir: Path, show: bool):
    fig, axes = plt.subplots(2, 4, figsize=(15, 7.5), constrained_layout=True)
    axes = axes.ravel()
    for ax, suite_dir in zip(axes, suite_dirs):
        extent, _, _, consensus = load_grid(suite_dir / "grid.txt")
        shared = load_shared(suite_dir / "shared_optima.txt")
        img = draw_image(ax, consensus, extent, suite_label(suite_dir), "viridis")
        add_shared_points(ax, shared)
        fig.colorbar(img, ax=ax, fraction=0.046, pad=0.03)
    fig.suptitle("FreePeaks Multi-party 2D: consensus landscape min(f_party0, f_party1)", fontsize=13)
    out = output_dir / "free_peaks_multiparty_consensus.png"
    fig.savefig(out, dpi=220)
    if show:
        plt.show()
    plt.close(fig)
    return out


def plot_parties(suite_dirs, output_dir: Path, show: bool):
    fig, axes = plt.subplots(len(suite_dirs), 3, figsize=(11, 3.0 * len(suite_dirs)), constrained_layout=True)
    for row, suite_dir in enumerate(suite_dirs):
        extent, f0, f1, consensus = load_grid(suite_dir / "grid.txt")
        shared = load_shared(suite_dir / "shared_optima.txt")
        fields = [(f0, "party 0", "magma"), (f1, "party 1", "cividis"), (consensus, "consensus", "viridis")]
        for col, (z, label, cmap) in enumerate(fields):
            ax = axes[row, col]
            title = f"{suite_label(suite_dir)} - {label}" if col == 0 else label
            img = draw_image(ax, z, extent, title, cmap)
            add_shared_points(ax, shared)
            fig.colorbar(img, ax=ax, fraction=0.046, pad=0.03)
    out = output_dir / "free_peaks_multiparty_parties.png"
    fig.savefig(out, dpi=200)
    if show:
        plt.show()
    plt.close(fig)
    return out



def pareto_ranks_2d_max(f0: np.ndarray, f1: np.ndarray) -> np.ndarray:
    values = np.column_stack([f0.ravel(), f1.ravel()])
    rounded = np.round(values, 10)
    unique, inverse = np.unique(rounded, axis=0, return_inverse=True)
    order = np.lexsort((-unique[:, 1], -unique[:, 0]))
    layer_max_f1: list[float] = []
    unique_rank = np.empty(unique.shape[0], dtype=np.int32)

    for idx in order:
        y = unique[idx, 1]
        lo, hi = 0, len(layer_max_f1)
        while lo < hi:
            mid = (lo + hi) // 2
            if y > layer_max_f1[mid]:
                hi = mid
            else:
                lo = mid + 1
        if lo == len(layer_max_f1):
            layer_max_f1.append(float(y))
        else:
            layer_max_f1[lo] = float(y)
        unique_rank[idx] = lo + 1

    return unique_rank[inverse].reshape(f0.shape)


def plot_party_objectives(suite_dirs, output_dir: Path, show: bool):
    fig, axes = plt.subplots(len(suite_dirs), 2, figsize=(8.5, 3.0 * len(suite_dirs)), constrained_layout=True)
    for row, suite_dir in enumerate(suite_dirs):
        extent, f0, f1, _ = load_grid(suite_dir / "grid.txt")
        shared = load_shared(suite_dir / "shared_optima.txt")
        for col, (z, label, cmap) in enumerate([(f0, "party 0 objective", "magma"), (f1, "party 1 objective", "cividis")]):
            ax = axes[row, col]
            title = f"{suite_label(suite_dir)} - {label}" if col == 0 else label
            img = draw_image(ax, z, extent, title, cmap)
            add_shared_points(ax, shared)
            fig.colorbar(img, ax=ax, fraction=0.046, pad=0.03)
    out = output_dir / "free_peaks_multiparty_party_objectives.png"
    fig.savefig(out, dpi=200)
    if show:
        plt.show()
    plt.close(fig)
    return out


def plot_rank_landscape(suite_dirs, output_dir: Path, show: bool):
    fig, axes = plt.subplots(2, 4, figsize=(15, 7.5), constrained_layout=True)
    axes = axes.ravel()
    for ax, suite_dir in zip(axes, suite_dirs):
        extent, f0, f1, _ = load_grid(suite_dir / "grid.txt")
        shared = load_shared(suite_dir / "shared_optima.txt")
        ranks = pareto_ranks_2d_max(f0, f1)
        rank_cap = int(np.nanpercentile(ranks, 99.5))
        clipped = np.minimum(ranks, max(rank_cap, 2))
        z = np.log1p(clipped)
        xs = np.linspace(extent[0], extent[1], z.shape[1])
        ys = np.linspace(extent[2], extent[3], z.shape[0])
        X, Y = np.meshgrid(xs, ys)
        levels = np.linspace(np.log1p(1), np.log1p(max(rank_cap, 2)), 28)
        img = ax.contourf(X, Y, z, levels=levels, cmap="viridis_r")
        ax.contour(X, Y, z, levels=10, colors="black", linewidths=0.25, alpha=0.35)
        add_shared_points(ax, shared)
        ax.set_title(f"{suite_label(suite_dir)} rank", fontsize=9)
        ax.set_xlim(extent[0], extent[1])
        ax.set_ylim(extent[2], extent[3])
        ax.set_xlabel("x0")
        ax.set_ylabel("x1")
        ax.set_aspect("equal")
        cb = fig.colorbar(img, ax=ax, fraction=0.046, pad=0.03)
        cb.set_label("log(1 + Pareto rank)")
    fig.suptitle("FreePeaks Multi-party 2D: bi-objective Pareto rank landscape, log scale (max f0, max f1)", fontsize=13)
    out = output_dir / "free_peaks_multiparty_rank_landscape.png"
    fig.savefig(out, dpi=220)
    if show:
        plt.show()
    plt.close(fig)
    return out


def plot_single_party_2x4(suite_dirs, party: int, output_dir: Path, show: bool):
    fig, axes = plt.subplots(2, 4, figsize=(15, 7.5), constrained_layout=True)
    axes = axes.ravel()
    for ax, suite_dir in zip(axes, suite_dirs):
        extent, f0, f1, _ = load_grid(suite_dir / "grid.txt")
        shared = load_shared(suite_dir / "shared_optima.txt")
        z = f0 if party == 0 else f1
        cmap = "magma" if party == 0 else "cividis"
        img = draw_image(ax, z, extent, f"{suite_label(suite_dir)} - party {party}", cmap)
        add_shared_points(ax, shared)
        fig.colorbar(img, ax=ax, fraction=0.046, pad=0.03)
    fig.suptitle(f"FreePeaks Multi-party 2D: party {party} objective landscapes", fontsize=13)
    out = output_dir / f"free_peaks_multiparty_party{party}_objective_2x4.png"
    fig.savefig(out, dpi=220)
    if show:
        plt.show()
    plt.close(fig)
    return out


def plot_single_party_3x4(suite_dirs, party: int, output_dir: Path, show: bool):
    fig, axes = plt.subplots(3, 4, figsize=(15, 10.8), constrained_layout=True)
    axes = axes.ravel()
    for ax, suite_dir in zip(axes, suite_dirs):
        extent, f0, f1, _ = load_grid(suite_dir / "grid.txt")
        shared = load_shared(suite_dir / "shared_optima.txt")
        z = f0 if party == 0 else f1
        cmap = "magma" if party == 0 else "cividis"
        img = draw_image(ax, z, extent, f"{suite_label(suite_dir)} - party {party}", cmap)
        add_shared_points(ax, shared)
        fig.colorbar(img, ax=ax, fraction=0.046, pad=0.03)
    for ax in axes[len(suite_dirs):]:
        ax.axis("off")
    fig.suptitle(f"FreePeaks Multi-party 2D p1-p12: party {party} objective landscapes", fontsize=13)
    out = output_dir / f"free_peaks_multiparty_p1_p12_party{party}_objective_3x4.png"
    fig.savefig(out, dpi=220)
    if show:
        plt.show()
    plt.close(fig)
    return out


def plot_rank_landscape_3x4(suite_dirs, output_dir: Path, show: bool):
    fig, axes = plt.subplots(3, 4, figsize=(15, 10.8), constrained_layout=True)
    axes = axes.ravel()
    for ax, suite_dir in zip(axes, suite_dirs):
        extent, f0, f1, _ = load_grid(suite_dir / "grid.txt")
        shared = load_shared(suite_dir / "shared_optima.txt")
        ranks = pareto_ranks_2d_max(f0, f1)
        rank_cap = int(np.nanpercentile(ranks, 99.5))
        clipped = np.minimum(ranks, max(rank_cap, 2))
        z = np.log1p(clipped)
        xs = np.linspace(extent[0], extent[1], z.shape[1])
        ys = np.linspace(extent[2], extent[3], z.shape[0])
        X, Y = np.meshgrid(xs, ys)
        levels = np.linspace(np.log1p(1), np.log1p(max(rank_cap, 2)), 28)
        img = ax.contourf(X, Y, z, levels=levels, cmap="viridis_r")
        ax.contour(X, Y, z, levels=10, colors="black", linewidths=0.25, alpha=0.35)
        add_shared_points(ax, shared)
        ax.set_title(f"{suite_label(suite_dir)} rank", fontsize=9)
        ax.set_xlim(extent[0], extent[1])
        ax.set_ylim(extent[2], extent[3])
        ax.set_xlabel("x0")
        ax.set_ylabel("x1")
        ax.set_aspect("equal")
        cb = fig.colorbar(img, ax=ax, fraction=0.046, pad=0.03)
        cb.set_label("log(1 + Pareto rank)")
    for ax in axes[len(suite_dirs):]:
        ax.axis("off")
    fig.suptitle("FreePeaks Multi-party 2D p1-p12: bi-objective Pareto rank landscapes", fontsize=13)
    out = output_dir / "free_peaks_multiparty_p1_p12_rank_landscape_3x4.png"
    fig.savefig(out, dpi=220)
    if show:
        plt.show()
    plt.close(fig)
    return out


def plot_twenty_optima_suite(suite_dir: Path, output_dir: Path, show: bool):
    extent, f0, f1, _ = load_grid(suite_dir / "grid.txt")
    shared = load_shared(suite_dir / "shared_optima.txt")
    ranks = pareto_ranks_2d_max(f0, f1)
    rank_cap = int(np.nanpercentile(ranks, 99.5))
    clipped = np.minimum(ranks, max(rank_cap, 2))
    z_rank = np.log1p(clipped)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), constrained_layout=True)
    fields = [
        (f0, "party 0 objective", "magma", None),
        (f1, "party 1 objective", "cividis", None),
        (z_rank, "rank fitness landscape", "viridis_r", "log(1 + Pareto rank)"),
    ]
    for ax, (z, title, cmap, cbar_label) in zip(axes, fields):
        if cbar_label is None:
            img = draw_image(ax, z, extent, title, cmap)
        else:
            xs = np.linspace(extent[0], extent[1], z.shape[1])
            ys = np.linspace(extent[2], extent[3], z.shape[0])
            X, Y = np.meshgrid(xs, ys)
            levels = np.linspace(np.log1p(1), np.log1p(max(rank_cap, 2)), 28)
            img = ax.contourf(X, Y, z, levels=levels, cmap=cmap)
            ax.contour(X, Y, z, levels=10, colors="black", linewidths=0.25, alpha=0.35)
            ax.set_title(title, fontsize=9)
            ax.set_xlim(extent[0], extent[1])
            ax.set_ylim(extent[2], extent[3])
            ax.set_xlabel("x0")
            ax.set_ylabel("x1")
            ax.set_aspect("equal")
        add_shared_points(ax, shared)
        cb = fig.colorbar(img, ax=ax, fraction=0.046, pad=0.03)
        if cbar_label:
            cb.set_label(cbar_label)
    fig.suptitle(f"{suite_label(suite_dir)}: 20 shared optima benchmark", fontsize=13)
    out = output_dir / "free_peaks_multiparty_twenty_optima_party_rank.png"
    fig.savefig(out, dpi=220)
    if show:
        plt.show()
    plt.close(fig)
    return out


def plot_mixed_consensus_row(suite_dirs, output_dir: Path, show: bool):
    fig, axes = plt.subplots(1, 4, figsize=(18, 4.6), constrained_layout=True)
    for ax, suite_dir in zip(axes, suite_dirs):
        extent, _, _, consensus = load_grid(suite_dir / "grid.txt")
        shared = load_shared(suite_dir / "shared_optima.txt")
        img = draw_image(ax, consensus, extent, suite_label(suite_dir), "viridis")
        add_shared_points(ax, shared)
        fig.colorbar(img, ax=ax, fraction=0.046, pad=0.03)
    fig.suptitle("FreePeaks Multi-party 2D: p9-p12 mixed consensus landscapes", fontsize=13)
    out = output_dir / "free_peaks_multiparty_mixed_p9_p12_consensus_row.png"
    fig.savefig(out, dpi=220)
    if show:
        plt.show()
    plt.close(fig)
    return out


def plot_mixed_rank_row(suite_dirs, output_dir: Path, show: bool):
    fig, axes = plt.subplots(1, 4, figsize=(18, 4.6), constrained_layout=True)
    for ax, suite_dir in zip(axes, suite_dirs):
        extent, f0, f1, _ = load_grid(suite_dir / "grid.txt")
        shared = load_shared(suite_dir / "shared_optima.txt")
        ranks = pareto_ranks_2d_max(f0, f1)
        rank_cap = int(np.nanpercentile(ranks, 99.5))
        clipped = np.minimum(ranks, max(rank_cap, 2))
        z = np.log1p(clipped)
        xs = np.linspace(extent[0], extent[1], z.shape[1])
        ys = np.linspace(extent[2], extent[3], z.shape[0])
        X, Y = np.meshgrid(xs, ys)
        levels = np.linspace(np.log1p(1), np.log1p(max(rank_cap, 2)), 28)
        img = ax.contourf(X, Y, z, levels=levels, cmap="viridis_r")
        ax.contour(X, Y, z, levels=10, colors="black", linewidths=0.25, alpha=0.35)
        add_shared_points(ax, shared)
        ax.set_title(f"{suite_label(suite_dir)} rank", fontsize=9)
        ax.set_xlim(extent[0], extent[1])
        ax.set_ylim(extent[2], extent[3])
        ax.set_xlabel("x0")
        ax.set_ylabel("x1")
        ax.set_aspect("equal")
        cb = fig.colorbar(img, ax=ax, fraction=0.046, pad=0.03)
        cb.set_label("log(1 + Pareto rank)")
    fig.suptitle("FreePeaks Multi-party 2D: p9-p12 mixed Pareto rank landscapes", fontsize=13)
    out = output_dir / "free_peaks_multiparty_mixed_p9_p12_rank_row.png"
    fig.savefig(out, dpi=220)
    if show:
        plt.show()
    plt.close(fig)
    return out


def main():
    parser = argparse.ArgumentParser(description="Plot sampled FreePeaks multi-party 2D suites.")
    parser.add_argument("--input", default="Visualization/free_peaks_multiparty_samples", help="Directory containing suite_* sample folders.")
    parser.add_argument("--output", default="Visualization/figures", help="Directory for output PNG figures.")
    parser.add_argument("--show", action="store_true", help="Display figures interactively after saving.")
    args = parser.parse_args()

    input_dir = Path(args.input)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    suite_dirs = sorted(
        (p for p in input_dir.glob("suite_*") if (p / "grid.txt").exists() and re.match(r"suite_\d+_p\d+_", p.name)),
        key=suite_id,
    )
    base_dirs = [p for p in suite_dirs if 1 <= suite_id(p) <= 8]
    mixed_dirs = [p for p in suite_dirs if 9 <= suite_id(p) <= 12]
    all_dirs = [p for p in suite_dirs if 1 <= suite_id(p) <= 12]
    twenty_dirs = [p for p in suite_dirs if suite_id(p) == 9]
    if len(base_dirs) != 8:
        raise RuntimeError(f"Expected 8 base suite directories under {input_dir}, found {len(base_dirs)}")
    if len(twenty_dirs) != 1:
        raise RuntimeError(f"Expected one suite_9 twenty-optima directory under {input_dir}, found {len(twenty_dirs)}")
    if len(mixed_dirs) != 4:
        raise RuntimeError(f"Expected suite_9..suite_12 mixed directories under {input_dir}, found {len(mixed_dirs)}")
    if len(all_dirs) != 12:
        raise RuntimeError(f"Expected suite_1..suite_12 directories under {input_dir}, found {len(all_dirs)}")

    outputs = [
        plot_consensus(base_dirs, output_dir, args.show),
        plot_parties(base_dirs, output_dir, args.show),
        plot_party_objectives(base_dirs, output_dir, args.show),
        plot_single_party_2x4(base_dirs, 0, output_dir, args.show),
        plot_single_party_2x4(base_dirs, 1, output_dir, args.show),
        plot_rank_landscape(base_dirs, output_dir, args.show),
        plot_twenty_optima_suite(twenty_dirs[0], output_dir, args.show),
        plot_mixed_consensus_row(mixed_dirs, output_dir, args.show),
        plot_mixed_rank_row(mixed_dirs, output_dir, args.show),
        plot_single_party_3x4(all_dirs, 0, output_dir, args.show),
        plot_single_party_3x4(all_dirs, 1, output_dir, args.show),
        plot_rank_landscape_3x4(all_dirs, output_dir, args.show),
    ]
    for out in outputs:
        print(out)


if __name__ == "__main__":
    main()
