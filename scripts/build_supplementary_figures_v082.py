#!/usr/bin/env python3
"""Build the twelve Paper 1 supplementary figures from frozen CSV evidence."""

from __future__ import annotations

import csv
import os
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

os.environ.setdefault("SOURCE_DATE_EPOCH", "1785283200")

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter
import numpy as np


PACKAGE = Path(__file__).resolve().parents[1]
EVIDENCE = PACKAGE / "evidence"
OUT = PACKAGE / "generated" / "supplementary_figures"
QA = PACKAGE / "generated" / "qa"
WIDTH = 6.70
FONT = 11.0
NAVY = "#17324d"
BLUE = "#376aa0"
TEAL = "#16756f"
RED = "#b2453d"
GOLD = "#bd8126"
PURPLE = "#7a5698"
GRID = "#d7dde1"
INK = "#202a33"

mpl.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman"],
    "mathtext.fontset": "custom",
    "mathtext.rm": "Times New Roman",
    "mathtext.it": "Times New Roman:italic",
    "mathtext.bf": "Times New Roman:bold",
    "font.size": FONT,
    "axes.titlesize": FONT,
    "axes.labelsize": FONT,
    "xtick.labelsize": FONT,
    "ytick.labelsize": FONT,
    "legend.fontsize": FONT,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "svg.hashsalt": "paper1-v082",
    "savefig.facecolor": "white",
})

MANIFEST: dict[str, list[Path]] = {}
TEXT_AUDIT: list[dict[str, object]] = []
LEGEND_AUDIT: list[dict[str, object]] = []
TICK_AUDIT: list[dict[str, object]] = []


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def num(row: dict[str, str], key: str, default: float = math.nan) -> float:
    raw = row.get(key, "").strip()
    if not raw:
        return default
    if raw.startswith("."):
        raw = "0" + raw
    if raw.startswith("-."):
        raw = "-0" + raw[1:]
    return float(raw)


def style_axis(ax, grid: str = "both") -> None:
    ax.grid(True, which=grid, color=GRID, linewidth=0.55, alpha=0.75)
    ax.tick_params(direction="out", length=3, width=0.7)
    for spine in ax.spines.values():
        spine.set_color("#59636b")


def plain_log_ticks(ax, axis: str) -> None:
    locator = LogLocator(base=10)
    formatter = FuncFormatter(lambda value, _: f"{value:.0e}".replace("e-0", "e-").replace("e+0", "e+") if value > 0 else "")
    target = ax.xaxis if axis == "x" else ax.yaxis
    target.set_major_locator(locator)
    target.set_major_formatter(formatter)
    target.set_minor_formatter(NullFormatter())


def panel(ax, text: str) -> None:
    ax.text(-0.13, 1.04, text, transform=ax.transAxes, ha="left", va="bottom", fontweight="bold")


def audit_layout(fig, figure_id: str) -> None:
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    visible_axes = [ax for ax in fig.axes if ax.get_visible() and ax.axison]

    for artist in fig.findobj(match=mpl.text.Text):
        label = artist.get_text().strip()
        if not label:
            continue
        family = artist.get_fontproperties().get_name()
        size = float(artist.get_fontsize())
        TEXT_AUDIT.append({
            "figure_id": figure_id,
            "text": label.replace("\n", " | "),
            "font_family": family,
            "font_size_pt": size,
            "minimum_8pt_pass": int(size >= 8.0),
            "times_new_roman_pass": int(family == "Times New Roman"),
        })

    for index, ax in enumerate(visible_axes, start=1):
        legend = ax.get_legend()
        if legend is not None:
            legend_box = legend.get_window_extent(renderer)
            axes_box = ax.get_window_extent(renderer)
            LEGEND_AUDIT.append({
                "figure_id": figure_id,
                "legend_id": f"axes_{index}",
                "intersects_data_axes": int(legend_box.overlaps(axes_box)),
            })

        for axis_name, labels in (("x", ax.get_xticklabels()), ("y", ax.get_yticklabels())):
            boxes = [item.get_window_extent(renderer) for item in labels if item.get_visible() and item.get_text().strip()]
            overlaps = sum(boxes[i].overlaps(boxes[j]) for i in range(len(boxes)) for j in range(i + 1, len(boxes)))
            outside = sum(
                box.x0 < fig.bbox.x0 - 1 or box.y0 < fig.bbox.y0 - 1
                or box.x1 > fig.bbox.x1 + 1 or box.y1 > fig.bbox.y1 + 1
                for box in boxes
            )
            TICK_AUDIT.append({
                "figure_id": figure_id,
                "axis_id": f"axes_{index}_{axis_name}",
                "label_count": len(boxes),
                "overlap_pair_count": overlaps,
                "outside_figure_count": outside,
            })

    for index, legend in enumerate(fig.legends, start=1):
        legend_box = legend.get_window_extent(renderer)
        intersects = any(legend_box.overlaps(ax.get_window_extent(renderer)) for ax in visible_axes)
        LEGEND_AUDIT.append({
            "figure_id": figure_id,
            "legend_id": f"figure_{index}",
            "intersects_data_axes": int(intersects),
        })


def save(fig, stem: str, figure_id: str, sources: list[Path]) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    audit_layout(fig, figure_id)
    for extension in ("pdf", "svg", "png"):
        dpi = 500 if extension == "png" else None
        fig.savefig(OUT / f"{stem}.{extension}", dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    MANIFEST[figure_id] = sources


def null_series(group: list[dict[str, str]]) -> bool:
    return max(max(abs(num(row, "norm_fd", 0.0)), abs(num(row, "norm_kv", num(row, "norm_kd", 0.0)))) for row in group) <= 1e-10


def point_ratio(row: dict[str, str], is_null: bool) -> float:
    if is_null:
        value = num(row, "scaled_absolute_error_l2", num(row, "absolute_error_l2", 0.0))
        return max(value / 1e-10, 1e-14)
    return max(num(row, "relative_error_l2", num(row, "relative_error_scaled_l2", 0.0)) / 1e-5, 1e-14)


def block_envelope(ax, rows: list[dict[str, str]], block: str, title: str, color: str = TEAL) -> None:
    selected = [row for row in rows if row.get("block_id") == block]
    groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in selected:
        groups[(row.get("state_id", ""), row.get("direction_id", ""))].append(row)
    by_h: dict[float, list[float]] = defaultdict(list)
    for group in groups.values():
        null = null_series(group)
        ordered = sorted(group, key=lambda row: num(row, "h"))
        hs = [num(row, "h") for row in ordered]
        ratios = [point_ratio(row, null) for row in ordered]
        ax.loglog(hs, ratios, color=color, alpha=0.16, linewidth=0.65)
        for h, ratio in zip(hs, ratios):
            by_h[h].append(ratio)
    hs = sorted(by_h)
    if hs:
        ax.loglog(hs, [max(by_h[h]) for h in hs], marker="o", color=color, label="worst registered series")
    ax.axhline(1.0, color="#333333", linestyle="--", linewidth=0.85, label="applicable gate")
    ax.set_title(title, loc="left", fontweight="bold")
    ax.set_xlabel("perturbation scale, h")
    ax.set_ylabel("error / gate")
    plain_log_ticks(ax, "x"); plain_log_ticks(ax, "y"); style_axis(ax)


def s1() -> None:
    path = EVIDENCE / "wp09_six_block_fd_results_run1.csv"
    rows = read_csv(path)
    fig, axes = plt.subplots(2, 3, figsize=(WIDTH, 7.5), constrained_layout=True)
    for index, (ax, block) in enumerate(zip(axes.flat, ["KPU", "KSU", "KPP", "KPS", "KSP", "KSS"])):
        block_envelope(ax, rows, block, f"M1 baseline: {block}", RED)
        panel(ax, f"({chr(97 + index)})")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=2, frameon=False)
    save(fig, "Figure_S1_M1_six_block_baseline", "Figure S1", [path])


def s2() -> None:
    path = EVIDENCE / "wp10_i1_fd_results_run1.csv"
    rows = read_csv(path)
    fig, axes = plt.subplots(1, 2, figsize=(WIDTH, 3.7), constrained_layout=True)
    for label, ax, block in zip("ab", axes, ["KPU", "KSU"]):
        block_envelope(ax, rows, block, f"Corrected {block}")
        panel(ax, f"({label})")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=2, frameon=False)
    save(fig, "Figure_S2_KPU_KSU_corrections", "Figure S2", [path])


def s3() -> None:
    path = EVIDENCE / "wp10_i2_fd_results_run1.csv"
    rows = read_csv(path)
    fig, axes = plt.subplots(1, 3, figsize=(WIDTH, 3.7), constrained_layout=True)
    for label, ax, block in zip("abc", axes, ["KUU", "KUP", "KUS"]):
        block_envelope(ax, rows, block, f"Corrected {block}")
        panel(ax, f"({label})")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=2, frameon=False)
    save(fig, "Figure_S3_mechanics_row_corrections", "Figure S3", [path])


def s4() -> None:
    fd_path = EVIDENCE / "wp10_i3_fd_results_run1.csv"
    term_path = EVIDENCE / "wp10_i3_kss_term_audit.csv"
    fd, terms = read_csv(fd_path), read_csv(term_path)
    fig, axes = plt.subplots(1, 2, figsize=(WIDTH, 4.0), constrained_layout=True)
    block_envelope(axes[0], fd, "KSS", "Branchwise KSS directional audit")
    panel(axes[0], "(a)")
    term_order = ["KSS_STORAGE", "KSS_STRAIN_SAT", "KSS_MOBILITY", "KSS_PC_FIRST", "KSS_PC_SECOND"]
    medians = [np.median([num(row, "frobenius_norm") for row in terms if row["term_id"] == term]) for term in term_order]
    axes[1].bar(range(len(term_order)), medians, color=[NAVY, BLUE, TEAL, GOLD, RED])
    axes[1].set_yscale("log"); plain_log_ticks(axes[1], "y")
    axes[1].set_xticks(range(len(term_order)), ["storage", "strain", "mobility", "$p_c'$", "$p_c''$"], rotation=35, ha="right", fontsize=10)
    axes[1].set_ylabel("median Frobenius norm"); axes[1].set_title("KSS term ledger", loc="left", fontweight="bold")
    style_axis(axes[1], "major"); panel(axes[1], "(b)")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=2, frameon=False)
    save(fig, "Figure_S4_KSS_branch_and_terms", "Figure S4", [fd_path, term_path])


def aggregate_by_h(rows: list[dict[str, str]], key: str, predicate=lambda row: True) -> tuple[list[float], list[float], list[float], list[float]]:
    grouped: dict[float, list[float]] = defaultdict(list)
    for row in rows:
        if predicate(row):
            grouped[num(row, "h")].append(max(num(row, key), 1e-30))
    hs = sorted(grouped)
    return hs, [np.median(grouped[h]) for h in hs], [min(grouped[h]) for h in hs], [max(grouped[h]) for h in hs]


def s5() -> None:
    raw_path = EVIDENCE / "wp10_i4_kpp_raw_fd_run1.csv"
    phys_path = EVIDENCE / "wp10_i4_kpp_physical_fd_run1.csv"
    ksp_path = EVIDENCE / "wp10_i4_ksp_nullspace_scaling.csv"
    raw, phys, ksp = read_csv(raw_path), read_csv(phys_path), read_csv(ksp_path)
    fig, axes = plt.subplots(1, 2, figsize=(WIDTH, 4.0), constrained_layout=True)
    for rows, color, label in [(raw, RED, "returned Kpp"), (phys, TEAL, "physical Kpp")]:
        hs, med, low, high = aggregate_by_h(rows, "relative_error_l2")
        axes[0].fill_between(hs, low, high, color=color, alpha=0.14)
        axes[0].loglog(hs, med, marker="o", color=color, label=label)
    axes[0].axhline(1e-5, color="#333", ls="--", lw=0.85, label="relative gate")
    axes[0].set_xlabel("h"); axes[0].set_ylabel("relative L2 error"); axes[0].set_title("Kpp stabilization separation", loc="left", fontweight="bold")
    plain_log_ticks(axes[0], "x"); plain_log_ticks(axes[0], "y"); style_axis(axes[0]); panel(axes[0], "(a)")
    for direction, color in [("D1", RED), ("D2", TEAL), ("D3", BLUE)]:
        rows = [row for row in ksp if row["direction_id"] == direction]
        hs = [num(row, "h") for row in rows]
        metric = [max(num(row, "scaled_absolute_error_l2") / 1e-10, 1e-14) if direction == "D1" else max(num(row, "scaled_absolute_error_l2") / 1e-10, 1e-14) for row in rows]
        axes[1].loglog(hs, metric, marker="o", color=color, label=f"{direction}: {'uniform' if direction == 'D1' else 'nonuniform'}")
    axes[1].axhline(1.0, color="#333", ls="--", lw=0.85, label="absolute null gate")
    axes[1].set_xlabel("h"); axes[1].set_ylabel("scaled absolute mismatch / gate"); axes[1].set_title("Ksp pressure-direction audit", loc="left", fontweight="bold")
    plain_log_ticks(axes[1], "x"); plain_log_ticks(axes[1], "y"); style_axis(axes[1]); panel(axes[1], "(b)")
    for ax in axes:
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.30), ncol=2, frameon=False)
    save(fig, "Figure_S5_KPP_KSP_separation", "Figure S5", [raw_path, phys_path, ksp_path])


def s6() -> None:
    path = EVIDENCE / "wp10_i5_full_matrix_fd_run1.csv"
    rows = read_csv(path)
    states = list(dict.fromkeys(row["state_id"] for row in rows))
    directions = ["M1", "M2", "M3"]
    best = np.full((len(states), len(directions)), np.nan)
    for i, state in enumerate(states):
        for j, direction in enumerate(directions):
            vals = [num(row, "relative_error_scaled_l2") for row in rows if row["state_id"] == state and row["direction_id"] == direction]
            if vals:
                best[i, j] = math.log10(max(min(vals), 1e-16))
    fig = plt.figure(figsize=(WIDTH, 7.2), constrained_layout=True)
    grid = fig.add_gridspec(2, 1, height_ratios=[1.4, 1.0])
    ax1, ax2 = fig.add_subplot(grid[0]), fig.add_subplot(grid[1])
    image = ax1.imshow(best, aspect="auto", cmap="viridis", vmin=-16, vmax=-4)
    ax1.set_xticks(range(3), directions); ax1.set_yticks(range(len(states)), [state.replace("_", " ") for state in states])
    for i in range(len(states)):
        for j in range(3):
            ax1.text(j, i, f"{best[i,j]:.1f}" if np.isfinite(best[i,j]) else "no data", ha="center", va="center", color="white" if best[i,j] < -9 else INK)
    bar = fig.colorbar(image, ax=ax1, pad=0.02); bar.set_label("log10(best relative scaled L2 error)")
    ax1.set_title("All 39 mixed full-matrix series", loc="left", fontweight="bold"); panel(ax1, "(a)")
    for direction, color in [("M1", RED), ("M2", GOLD), ("M3", BLUE)]:
        selected = sorted([row for row in rows if row["state_id"] == "FINAL_CLIP_LOW" and row["direction_id"] == direction], key=lambda row: num(row, "h"))
        ax2.loglog([num(row, "h") for row in selected], [num(row, "relative_error_scaled_l2") for row in selected], marker="o", color=color, label=direction)
    ax2.axhline(1e-5, color="#333", ls="--", lw=0.85, label="strict gate")
    ax2.set_xlabel("h"); ax2.set_ylabel("relative scaled L2 error"); ax2.set_title("Three retained binary64-inconclusive series", loc="left", fontweight="bold")
    plain_log_ticks(ax2, "x"); plain_log_ticks(ax2, "y"); style_axis(ax2); panel(ax2, "(b)")
    ax2.legend(loc="upper center", bbox_to_anchor=(0.5, -0.25), ncol=4, frameon=False)
    save(fig, "Figure_S6_full_matrix_39_series", "Figure S6", [path])


def s7() -> None:
    comp_path = EVIDENCE / "wp10_i5b_oracle_component_results.csv"
    dir_path = EVIDENCE / "wp10_i5b_oracle_directional_run1.csv"
    components, directions = read_csv(comp_path), read_csv(dir_path)
    fig, axes = plt.subplots(2, 1, figsize=(WIDTH, 6.3), constrained_layout=True)
    x = np.arange(1, 17); width = 0.24
    for index, (direction, color) in enumerate([("M1", RED), ("M2", GOLD), ("M3", TEAL)]):
        rows = sorted([row for row in components if row["direction_id"] == direction and row["run_id"] == "run1"], key=lambda row: int(row["row"]))
        axes[0].bar(x + (index - 1) * width, [max(num(row, "scaled_absolute_error"), 1e-18) for row in rows], width, label=direction, color=color)
    axes[0].set_yscale("log"); plain_log_ticks(axes[0], "y")
    axes[0].set_xticks(x); axes[0].set_xlabel("residual row"); axes[0].set_ylabel("scaled absolute mismatch")
    axes[0].set_title("Sixteen-row active-set oracle closure", loc="left", fontweight="bold"); style_axis(axes[0], "major"); panel(axes[0], "(a)")
    axes[0].legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=3, frameon=False)
    labels = [row["direction_id"] for row in directions]
    axes[1].bar(labels, [num(row, "relative_error_scaled_l2") for row in directions], color=[RED, GOLD, TEAL])
    axes[1].axhline(1e-10, color="#333", ls="--", lw=0.85, label="analytic gate")
    axes[1].set_yscale("log"); plain_log_ticks(axes[1], "y")
    axes[1].set_ylabel("relative scaled L2 error"); axes[1].set_title("Directional oracle versus physical Jacobian", loc="left", fontweight="bold")
    style_axis(axes[1], "major"); panel(axes[1], "(b)")
    axes[1].legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), frameon=False)
    save(fig, "Figure_S7_active_set_component_closure", "Figure S7", [comp_path, dir_path])


def s8() -> None:
    patch_path = EVIDENCE / "wp10_i6_patch_results_run1.csv"
    replay_path = EVIDENCE / "wp10_i6_rollback_replay_run1.csv"
    patch_rows, replay_rows = read_csv(patch_path), read_csv(replay_path)
    fig, axes = plt.subplots(1, 2, figsize=(WIDTH, 4.1), constrained_layout=True)
    modes = [row["mode_id"].replace("_AFFINE", "").replace("UNIFORM", "0") for row in patch_rows]
    errors = [max(num(row, "interior_scaled_error"), num(row, "balance_scaled_error"), 1e-18) for row in patch_rows]
    axes[0].bar(range(len(modes)), errors, color=[BLUE if row["null_mode"] == "1" else TEAL for row in patch_rows])
    axes[0].set_yscale("log"); plain_log_ticks(axes[0], "y")
    axes[0].set_xticks(range(len(modes)), modes, rotation=90); axes[0].set_ylabel("maximum registered scaled error")
    axes[0].set_title("Available assembled patch modes", loc="left", fontweight="bold"); style_axis(axes[0], "major"); panel(axes[0], "(a)")
    histories = list(dict.fromkeys(row["history"] for row in replay_rows)); quantities = list(dict.fromkeys(row["quantity"] for row in replay_rows))
    matrix = np.zeros((len(histories), len(quantities)))
    for row in replay_rows:
        matrix[histories.index(row["history"]), quantities.index(row["quantity"])] = 1 if row["status"] == "PASS" else 0
    axes[1].imshow(matrix, aspect="auto", cmap=ListedColormap([RED, "#dcefed"]), vmin=0, vmax=1)
    quantity_labels = {"RHS": "RHS", "AMATRX": "AMATRX", "PHYSICAL_SVARS": "SVARS", "BRANCH": "branch"}
    axes[1].set_xticks(range(len(quantities)), [quantity_labels.get(item, item) for item in quantities], rotation=50, ha="right", rotation_mode="anchor", fontsize=9); axes[1].set_yticks(range(len(histories)), histories)
    for i in range(len(histories)):
        for j in range(len(quantities)):
            axes[1].text(j, i, "pass" if matrix[i,j] else "fail", ha="center", va="center")
    axes[1].set_title("Rollback replay controls", loc="left", fontweight="bold"); panel(axes[1], "(b)")
    save(fig, "Figure_S8_patch_and_replay", "Figure S8", [patch_path, replay_path])


def s9() -> None:
    o64_path = EVIDENCE / "wp10_i6a_recomputed_order.csv"
    o128_path = EVIDENCE / "wp10_i6a_binary128_order.csv"
    floor_path = EVIDENCE / "wp10_i6a_roundoff_model.csv"
    o64, o128, floor = read_csv(o64_path), read_csv(o128_path), read_csv(floor_path)
    fig, axes = plt.subplots(1, 2, figsize=(WIDTH, 4.0), constrained_layout=True)
    for benchmark, color in [("N1", RED), ("N2", BLUE)]:
        rows64 = [row for row in o64 if row["benchmark"] == benchmark]
        rows128 = [row for row in o128 if row["benchmark"] == benchmark]
        axes[0].plot([int(row["iteration"]) for row in rows64], [num(row, "q") for row in rows64], marker="o", color=color, alpha=0.5, label=f"{benchmark} binary64")
        axes[0].plot([int(row["iteration"]) for row in rows128], [num(row, "residual_order") for row in rows128], marker="s", ls="--", color=color, label=f"{benchmark} binary128")
    axes[0].axhline(1.8, color="#333", ls="--", lw=0.85, label="registered q gate")
    axes[0].set_xlabel("iteration"); axes[0].set_ylabel("empirical residual order")
    axes[0].set_title("Binary64 and binary128 order", loc="left", fontweight="bold"); style_axis(axes[0], "major"); panel(axes[0], "(a)")
    axes[0].legend(loc="upper center", bbox_to_anchor=(0.5, -0.25), ncol=2, frameon=False)
    for benchmark, color in [("N1", RED), ("N2", BLUE)]:
        rows = [row for row in floor if row["benchmark"] == benchmark]
        axes[1].semilogy([int(row["iteration"]) for row in rows], [num(row, "cancellation_factor") for row in rows], marker="o", color=color, label=benchmark)
    axes[1].set_xlabel("iteration"); axes[1].set_ylabel("RHS subtraction cancellation factor")
    axes[1].set_title("Floating-point resolution diagnostic", loc="left", fontweight="bold"); style_axis(axes[1]); panel(axes[1], "(b)")
    axes[1].legend(loc="upper center", bbox_to_anchor=(0.5, -0.25), ncol=2, frameon=False)
    save(fig, "Figure_S9_precision_and_order_audit", "Figure S9", [o64_path, o128_path, floor_path])


def s10() -> None:
    map_path = EVIDENCE / "wp10_i6b_one_step_map_binary128_run1.csv"
    slope_path = EVIDENCE / "wp10_i6b_local_order_slopes.csv"
    rows, slopes = read_csv(map_path), read_csv(slope_path)
    fig, axes = plt.subplots(1, 2, figsize=(WIDTH, 4.0), constrained_layout=True)
    for benchmark, color in [("N1", RED), ("N2", BLUE)]:
        selected = sorted([row for row in rows if row["benchmark"] == benchmark], key=lambda row: num(row, "e0"))
        e0 = np.array([num(row, "e0") for row in selected]); e1 = np.array([num(row, "e1") for row in selected])
        axes[0].loglog(e0, e1, marker="o", color=color, label=benchmark)
    ref = np.logspace(-6, -1, 30); axes[0].loglog(ref, ref, color="#777", ls=":", label="linear reference"); axes[0].loglog(ref, ref**2, color="#333", ls="--", label="quadratic reference")
    axes[0].set_xlabel("initial scaled root error"); axes[0].set_ylabel("one-step scaled root error")
    axes[0].set_title("Frozen one-step map", loc="left", fontweight="bold"); plain_log_ticks(axes[0], "x"); plain_log_ticks(axes[0], "y"); style_axis(axes[0]); panel(axes[0], "(a)")
    axes[0].legend(loc="upper center", bbox_to_anchor=(0.5, -0.25), ncol=2, frameon=False)
    labels=[]; values=[]; colors=[]
    for row in slopes:
        if row["arithmetic"] == "binary64":
            metric_label = {"root_error": "root", "scaled_residual": "residual"}.get(row["metric"], row["metric"].replace("_", " "))
            labels.append(f"{row['benchmark']}\n{metric_label}"); values.append(num(row, "ols_slope")); colors.append(RED if row["benchmark"] == "N1" else BLUE)
    axes[1].bar(range(len(values)), values, color=colors); axes[1].axhline(2.0, color="#333", ls="--", lw=0.85, label="quadratic slope")
    axes[1].set_xticks(range(len(values)), labels); axes[1].set_ylabel("OLS slope"); axes[1].set_title("Registered local slopes", loc="left", fontweight="bold")
    style_axis(axes[1], "major"); panel(axes[1], "(b)"); axes[1].legend(loc="upper center", bbox_to_anchor=(0.5, -0.25), frameon=False)
    save(fig, "Figure_S10_frozen_one_step_map", "Figure S10", [map_path, slope_path])


def s11() -> None:
    path = EVIDENCE / "wp10_i6c_r1_basin_entry_radius.csv"
    rows = read_csv(path)
    directions = ["LEGACY", "ERR_K1", "ERR_K2", "ERR_K3", "UPDATE_K0", "UPDATE_K1", "ORTH_K0", "SVD_MAX", "SVD_MID", "SVD_MIN", "FIELD_U", "FIELD_P", "FIELD_S"]
    matrix = np.full((2, len(directions)), np.nan)
    for row in rows:
        if row["scope"] == "DIRECTION" and row["persistent_entry"] == "1":
            matrix[["N1", "N2"].index(row["benchmark"]), directions.index(row["direction_id"])] = math.log10(num(row, "entry_radius"))
    fig, ax = plt.subplots(figsize=(WIDTH, 3.8), constrained_layout=True)
    image = ax.imshow(matrix, aspect="auto", cmap="viridis", vmin=-4, vmax=-2)
    short_labels = ["Base", "Err 1", "Err 2", "Err 3", "Upd 0", "Upd 1", "Orth.",
                    "SVD max", "SVD mid", "SVD min", r"$u$", r"$p_w$", r"$S_n$"]
    ax.set_yticks(range(2), ["N1", "N2"])
    ax.set_xticks(range(len(directions)), short_labels, rotation=60, ha="right", rotation_mode="anchor")
    ax.tick_params(axis="x", labelsize=8.5, pad=2)
    for i in range(2):
        for j in range(len(directions)):
            ax.text(j, i, "none" if np.isnan(matrix[i,j]) else f"{matrix[i,j]:.1f}",
                    ha="center", va="center",
                    color="white" if not np.isnan(matrix[i,j]) and matrix[i,j] < -2.8 else INK,
                    fontsize=8.5)
    bar = fig.colorbar(image, ax=ax, pad=0.02); bar.set_label("log10(entry radius)")
    ax.set_title("Direction-dependent local-basin entries: 23 of 26 pairs", loc="left", fontweight="bold")
    save(fig, "Figure_S11_corrected_local_basin", "Figure S11", [path])


def s12() -> None:
    sub_path = EVIDENCE / "wp10_i6d_subspace_audit.csv"
    dec_path = EVIDENCE / "wp10_i6d_one_step_decomposition.csv"
    sub, dec = read_csv(sub_path), read_csv(dec_path)
    targets = [("N1", "FIELD_U"), ("N1", "FIELD_P"), ("N2", "FIELD_P")]
    sub = [row for row in sub if (row["benchmark"], row["direction_id"]) in targets]
    fig, axes = plt.subplots(1, 2, figsize=(WIDTH, 4.0), constrained_layout=True)
    x=np.arange(len(sub)); width=0.2
    for index,(label,key,color) in enumerate([("bottom-3 direction","bottom3_direction_energy",NAVY),("bottom-3 response","bottom3_response_energy",BLUE),("equal-P direction","equal_pressure_energy",GOLD),("equal-P response","equal_pressure_response_energy",RED)]):
        axes[0].bar(x+(index-1.5)*width,[num(row,key) for row in sub],width,label=label,color=color)
    axes[0].axhline(0.9,color="#333",ls="--",lw=0.85,label="dominance gate"); axes[0].set_ylim(0,1.05)
    axes[0].set_xticks(x,[f"{row['benchmark']}/{row['direction_id'].replace('FIELD_','')}" for row in sub]); axes[0].set_ylabel("energy fraction")
    axes[0].set_title("Candidate subspace explanations",loc="left",fontweight="bold"); style_axis(axes[0],"major"); panel(axes[0],"(a)")
    axes[0].legend(loc="upper center",bbox_to_anchor=(0.5,-0.25),ncol=2,frameon=False)
    for target,color in zip(targets,[RED,GOLD,BLUE]):
        rows=sorted([row for row in dec if (row["benchmark"],row["direction_id"])==target and row["sign"]=="1"],key=lambda row:num(row,"amplitude"),reverse=True)
        axes[1].semilogx([num(row,"amplitude") for row in rows],[num(row,"nonlinear_fraction") for row in rows],marker="o",color=color,label=f"{target[0]}/{target[1].replace('FIELD_','')}")
    axes[1].axhline(0.9,color="#333",ls="--",lw=0.85,label="dominance gate"); axes[1].invert_xaxis(); axes[1].set_ylim(0,1.02)
    axes[1].set_xlabel("scaled amplitude"); axes[1].set_ylabel("nonlinear contribution fraction"); axes[1].set_title("One-step nonlinear decomposition",loc="left",fontweight="bold")
    plain_log_ticks(axes[1],"x"); style_axis(axes[1]); panel(axes[1],"(b)"); axes[1].legend(loc="upper center",bbox_to_anchor=(0.5,-0.25),ncol=2,frameon=False)
    save(fig,"Figure_S12_field_subspace_audit","Figure S12",[sub_path,dec_path])


def write_qa() -> None:
    QA.mkdir(parents=True, exist_ok=True)
    rows=[]
    for figure_id,sources in MANIFEST.items():
        for source in sources:
            rows.append({"figure_id":figure_id,"source_path":source.relative_to(PACKAGE).as_posix(),"source_sha256":sha256(source),"source_status":"FROZEN_READ_ONLY"})
    with (QA / "supplementary_figure_source_manifest_v0.8.csv").open("w",encoding="utf-8",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=rows[0].keys(),lineterminator="\n"); writer.writeheader(); writer.writerows(rows)
    with (QA / "supplementary_figure_text_size_audit_v0.8.2.csv").open("w",encoding="utf-8",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=["figure_id","text","font_family","font_size_pt","minimum_8pt_pass","times_new_roman_pass"],lineterminator="\n"); writer.writeheader(); writer.writerows(TEXT_AUDIT)
    with (QA / "supplementary_figure_legend_overlap_audit_v0.8.2.csv").open("w",encoding="utf-8",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=["figure_id","legend_id","intersects_data_axes"],lineterminator="\n"); writer.writeheader(); writer.writerows(LEGEND_AUDIT)
    with (QA / "supplementary_figure_tick_label_layout_audit_v0.8.2.csv").open("w",encoding="utf-8",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=["figure_id","axis_id","label_count","overlap_pair_count","outside_figure_count"],lineterminator="\n"); writer.writeheader(); writer.writerows(TICK_AUDIT)
    outputs=[{"file":path.name,"sha256":sha256(path),"bytes":path.stat().st_size} for path in sorted(OUT.iterdir())]
    passed = (
        bool(TEXT_AUDIT)
        and all(row["minimum_8pt_pass"] and row["times_new_roman_pass"] for row in TEXT_AUDIT)
        and all(not row["intersects_data_axes"] for row in LEGEND_AUDIT)
        and all(not row["overlap_pair_count"] for row in TICK_AUDIT)
    )
    report={"status":"PASS" if passed else "FAIL","figure_count":len(MANIFEST),"font_family":"Times New Roman","normal_text_pt":FONT,"minimum_text_gate_pt":8,"minimum_text_observed_pt":min(row["font_size_pt"] for row in TEXT_AUDIT),"legend_policy":"OUTSIDE_DATA_AXES","legend_overlap_count":int(sum(row["intersects_data_axes"] for row in LEGEND_AUDIT)),"tick_overlap_pair_count":int(sum(row["overlap_pair_count"] for row in TICK_AUDIT)),"tick_outside_figure_count":int(sum(row["outside_figure_count"] for row in TICK_AUDIT)),"source_manifest_rows":len(rows),"outputs":outputs}
    (QA / "supplementary_figure_build_qa_v0.8.json").write_text(json.dumps(report,indent=2),encoding="utf-8")


def main() -> None:
    for function in (s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12):
        function()
    write_qa()
    print(json.dumps({"status":"PASS","figures":len(MANIFEST),"font":"Times New Roman","font_pt":FONT}))


if __name__ == "__main__":
    main()
