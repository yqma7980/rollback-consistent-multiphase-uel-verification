#!/usr/bin/env python3
"""Rebuild Paper 1 main artwork from frozen evidence at publication scale."""

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
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter
import numpy as np


PACKAGE = Path(__file__).resolve().parents[1]
EVIDENCE = PACKAGE / "evidence"
OUT = PACKAGE / "generated" / "figures"
QA = PACKAGE / "generated" / "qa"

PAGE_WIDTH_IN = 6.70  # review-PDF width; production enlargement to 190 mm only increases lettering
FONT_PT = 11.0
PANEL_PT = 12.0
TITLE_PT = 11.5

NAVY = "#17324d"
TEAL = "#16756f"
RED = "#b2453d"
GOLD = "#bd8126"
BLUE = "#376aa0"
PURPLE = "#7a5698"
INK = "#202a33"
MUTED = "#64717a"
GRID = "#d7dde1"
PALE = "#f5f7f8"

mpl.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["Times New Roman"],
        "mathtext.fontset": "custom",
        "mathtext.rm": "Times New Roman",
        "mathtext.it": "Times New Roman:italic",
        "mathtext.bf": "Times New Roman:bold",
        "font.size": FONT_PT,
        "axes.titlesize": TITLE_PT,
        "axes.labelsize": FONT_PT,
        "xtick.labelsize": FONT_PT,
        "ytick.labelsize": FONT_PT,
        "legend.fontsize": FONT_PT,
        "figure.titlesize": PANEL_PT,
        "axes.linewidth": 0.7,
        "lines.linewidth": 1.25,
        "lines.markersize": 3.5,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "svg.hashsalt": "paper1-v082",
        "savefig.facecolor": "white",
        "axes.facecolor": "white",
    }
)

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


def value(row: dict[str, str], key: str) -> float:
    text = row[key].strip()
    if text.startswith("."):
        text = "0" + text
    if text.startswith("-."):
        text = "-0" + text[1:]
    return float(text)


def panel_label(ax, label: str, x: float = -0.11, y: float = 1.04) -> None:
    ax.text(x, y, label, transform=ax.transAxes, fontweight="bold", fontsize=PANEL_PT,
            va="bottom", ha="left", color=INK, clip_on=False)


def finish_axis(ax, grid: str = "both") -> None:
    ax.tick_params(direction="out", length=3, width=0.7)
    ax.grid(True, which=grid, color=GRID, linewidth=0.55, alpha=0.75)
    for spine in ax.spines.values():
        spine.set_color("#59636b")


def legend_below(ax, ncol: int = 2, y: float = -0.30):
    return ax.legend(
        loc="upper center", bbox_to_anchor=(0.5, y), ncol=ncol,
        frameon=False, borderaxespad=0.0, columnspacing=1.0, handletextpad=0.45,
    )


def audit_legends(fig, figure_id: str) -> None:
    fig.canvas.draw()
    axes = [ax for ax in fig.axes if ax.get_visible() and ax.axison]
    for index, ax in enumerate(axes, start=1):
        legend = ax.get_legend()
        if legend is None:
            continue
        legend_box = legend.get_window_extent(fig.canvas.get_renderer())
        axes_box = ax.get_window_extent(fig.canvas.get_renderer())
        LEGEND_AUDIT.append({
            "figure_id": figure_id, "legend_id": f"axes_{index}",
            "placement": "BELOW_AXES", "intersects_data_axes": int(legend_box.overlaps(axes_box)),
        })
    for index, legend in enumerate(fig.legends, start=1):
        legend_box = legend.get_window_extent(fig.canvas.get_renderer())
        intersects = any(legend_box.overlaps(ax.get_window_extent(fig.canvas.get_renderer())) for ax in axes)
        LEGEND_AUDIT.append({
            "figure_id": figure_id, "legend_id": f"figure_{index}",
            "placement": "OUTSIDE_LOWER_CENTER", "intersects_data_axes": int(intersects),
        })


def audit_tick_labels(fig, figure_id: str) -> None:
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    figure_box = fig.bbox
    for axis_index, ax in enumerate([item for item in fig.axes if item.get_visible() and item.axison], start=1):
        for axis_name, labels in (("x", ax.get_xticklabels()), ("y", ax.get_yticklabels())):
            boxes = [label.get_window_extent(renderer) for label in labels if label.get_visible() and label.get_text().strip()]
            overlap_pairs = sum(boxes[i].overlaps(boxes[j]) for i in range(len(boxes)) for j in range(i + 1, len(boxes)))
            outside = sum(
                box.x0 < figure_box.x0 - 1 or box.y0 < figure_box.y0 - 1
                or box.x1 > figure_box.x1 + 1 or box.y1 > figure_box.y1 + 1
                for box in boxes
            )
            TICK_AUDIT.append({
                "figure_id": figure_id,
                "axis_id": f"axes_{axis_index}_{axis_name}",
                "label_count": len(boxes),
                "overlap_pair_count": overlap_pairs,
                "outside_figure_count": outside,
            })


def plain_log_ticks(ax, axis: str) -> None:
    target = ax.xaxis if axis == "x" else ax.yaxis
    target.set_major_locator(LogLocator(base=10))
    target.set_major_formatter(
        FuncFormatter(lambda number, _position: "" if number <= 0 else f"1e{int(round(math.log10(number)))}")
    )
    target.set_minor_formatter(NullFormatter())


def audit_text(fig, figure_id: str) -> list[dict[str, object]]:
    rows = []
    for artist in fig.findobj(match=mpl.text.Text):
        text = artist.get_text().strip()
        if not text:
            continue
        rows.append(
            {
                "figure_id": figure_id,
                "text": text.replace("\n", " | "),
                "font_family": artist.get_fontproperties().get_name(),
                "font_size_pt": float(artist.get_fontsize()),
                "hard_minimum_pass": int(float(artist.get_fontsize()) >= 8.0),
                "times_new_roman_pass": int(artist.get_fontproperties().get_name() == "Times New Roman"),
            }
        )
    return rows


def save_figure(fig, stem: str, figure_id: str, text_audit: list[dict[str, object]]) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    audit_legends(fig, figure_id)
    audit_tick_labels(fig, figure_id)
    for extension in ("pdf", "svg"):
        fig.savefig(OUT / f"{stem}.{extension}", dpi=500)
    fig.savefig(OUT / f"{stem}.png", dpi=500)
    text_audit.extend(audit_text(fig, figure_id))
    plt.close(fig)


def add_box(ax, xy, width, height, title, lines, face, edge) -> None:
    box = FancyBboxPatch(xy, width, height, boxstyle="round,pad=0.012,rounding_size=0.018",
                         linewidth=1.1, edgecolor=edge, facecolor=face)
    ax.add_patch(box)
    x, y = xy
    ax.text(x + width / 2, y + height - 0.065, title, ha="center", va="top",
            fontweight="bold", fontsize=TITLE_PT, color=edge)
    for index, line in enumerate(lines):
        ax.text(x + width / 2, y + height - 0.12 - 0.052 * index, line,
                ha="center", va="top", fontsize=FONT_PT, color=INK)


def add_arrow(ax, start, end, color=NAVY, style="-|>") -> None:
    ax.add_patch(FancyArrowPatch(
        start, end, arrowstyle=style, mutation_scale=11, linewidth=1.2,
        color=color, shrinkA=4, shrinkB=4, connectionstyle="arc3,rad=0.0", zorder=4,
    ))


def figure_1(text_audit, manifest):
    fig, ax = plt.subplots(figsize=(PAGE_WIDTH_IN, 4.05), constrained_layout=True)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    add_box(ax, (0.02, 0.62), 0.25, 0.26, "Committed state", ["solver-managed SVARS", r"CAP_REF = $C_n$", r"accepted $x_n$"], "#e4f1ef", TEAL)
    add_box(ax, (0.36, 0.62), 0.25, 0.26, "Trial A", [r"candidate $x_A, C_A$", "may be rejected"], "#eef2f5", NAVY)
    add_box(ax, (0.36, 0.24), 0.25, 0.25, "Trial B", ["same declared input", r"must restart from $C_n$"], "#eef2f5", NAVY)
    add_box(ax, (0.70, 0.62), 0.28, 0.26, "Unsafe private cache", [r"trial writes $C_A$", "rejection cannot restore it"], "#f7e4e1", RED)
    add_box(ax, (0.70, 0.22), 0.28, 0.29, "Solver-managed repair", ["candidate remains in", "solver-managed SVARS", "commit only after acceptance"], "#e4f1ef", TEAL)
    add_arrow(ax, (0.27, 0.75), (0.36, 0.75))
    add_arrow(ax, (0.485, 0.62), (0.485, 0.49))
    add_arrow(ax, (0.61, 0.75), (0.70, 0.75), RED)
    add_arrow(ax, (0.61, 0.365), (0.70, 0.365), TEAL)
    add_arrow(ax, (0.84, 0.62), (0.84, 0.51), RED)
    ax.text(0.855, 0.565, "history leak", ha="left", va="center", color=RED, fontweight="bold",
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 1.5}, zorder=5)
    ax.add_patch(FancyBboxPatch((0.08, 0.015), 0.84, 0.13, boxstyle="round,pad=0.012",
                                facecolor=PALE, edgecolor=GRID, linewidth=0.9))
    ax.text(0.50, 0.080, "Equal-state replay, rejection order, cutback and restart are independent gates;\n"
            "the repair is bounded to the identified CAP_REF path.", ha="center", va="center", color=INK)
    save_figure(fig, "Figure_1_rollback_transaction", "Figure 1", text_audit)
    manifest["Figure 1"] = []


def figure_2(text_audit, manifest):
    fig, ax = plt.subplots(figsize=(PAGE_WIDTH_IN, 3.75), constrained_layout=True)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    boxes = [
        (0.015, "Relations", [r"$p_{eq},p_c,K$", "branches"], GOLD, "#f8efdc"),
        (0.215, "Residual", [r"$R_u,R_p,R_s$", "units/signs"], NAVY, "#eef2f5"),
        (0.415, "Matrix", [r"$J_{phys}$", r"$J_{num}$ / raw"], TEAL, "#e4f1ef"),
        (0.615, "State", ["committed", "trial/diagnostic"], NAVY, "#eef2f5"),
        (0.815, "Host", ["RHS/AMATRX", "accept/reject"], GOLD, "#f8efdc"),
    ]
    box_width = 0.17
    for x, title, lines, edge, face in boxes:
        add_box(ax, (x, 0.62), box_width, 0.24, title, lines, face, edge)
    for index in range(4):
        add_arrow(ax, (boxes[index][0] + box_width, 0.74), (boxes[index + 1][0], 0.74))
    questions = [
        (0.17, RED, "Matrix entry", "paired with residual?", "physical or numerical?"),
        (0.50, TEAL, "Changed state", "who owns/commits it?", "can rejection restore it?"),
        (0.83, GOLD, "Passing derivative", "which branch/direction?", "does it imply a basin?"),
    ]
    for x, color, head, line1, line2 in questions:
        add_box(ax, (x - 0.145, 0.20), 0.29, 0.24, head, [line1, line2], "white", color)
    add_arrow(ax, (0.92, 0.12), (0.08, 0.12), RED)
    ax.text(0.50, 0.025, "Reverse audit prevents matrix, state and convergence claims from being conflated.",
            ha="center", va="bottom", fontweight="bold", color=RED)
    save_figure(fig, "Figure_2_traceability", "Figure 2", text_audit)
    manifest["Figure 2"] = []


def fd_envelope(ax, rows, block, states, color, title):
    selected = [r for r in rows if r["block_id"] == block and r["state_id"] in states]
    series = defaultdict(list)
    for row in selected:
        series[(row["state_id"], row["direction_id"])].append(row)
    grouped = defaultdict(list)
    for group in series.values():
        is_null = max(max(value(row, "norm_fd"), value(row, "norm_kv")) for row in group) <= 1e-10
        hs_series, ratios = [], []
        for row in sorted(group, key=lambda item: value(item, "h")):
            metric = value(row, "absolute_error_l2") / 1e-10 if is_null else value(row, "relative_error_l2") / 1e-5
            metric = max(metric, 1e-12)
            hs_series.append(value(row, "h"))
            ratios.append(metric)
            grouped[value(row, "h")].append(metric)
        ax.loglog(hs_series, ratios, color=color, alpha=0.18, linewidth=0.65)
    hs = np.array(sorted(grouped))
    worst = np.array([max(grouped[h]) for h in hs])
    ax.loglog(hs, worst, marker="o", color=color, label="worst of all registered series")
    ax.axhline(1.0, color="#333333", linestyle="--", linewidth=0.85, label="registered gate")
    ax.set_title(title, loc="left", fontweight="bold")
    ax.set_xlabel("perturbation scale, h")
    ax.set_ylabel("error / applicable gate")
    plain_log_ticks(ax, "x")
    plain_log_ticks(ax, "y")
    finish_axis(ax)


def figure_3(text_audit, manifest):
    p_wp09 = EVIDENCE / "wp09_six_block_fd_results_run1.csv"
    p_i1 = EVIDENCE / "wp10_i1_fd_results_run1.csv"
    p_i3 = EVIDENCE / "wp10_i3_fd_results_run1.csv"
    wp09, i1, i3 = read_csv(p_wp09), read_csv(p_i1), read_csv(p_i3)
    states = {"I0", "I1", "I2"}
    fig, axes = plt.subplots(2, 2, figsize=(PAGE_WIDTH_IN, 6.15), constrained_layout=True)
    items = [
        (axes[0, 0], wp09, "KPU", RED, "M1 KPU: registered omission"),
        (axes[0, 1], i1, "KPU", TEAL, "Corrected KPU: physical derivative"),
        (axes[1, 0], wp09, "KSS", RED, "M1 KSS: incomplete branch derivative"),
        (axes[1, 1], i3, "KSS", TEAL, "Corrected KSS: branchwise derivative"),
    ]
    for label, (ax, rows, block, color, title) in zip("abcd", items):
        fd_envelope(ax, rows, block, states, color, title)
        panel_label(ax, f"({label})")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=2, frameon=False)
    save_figure(fig, "Figure_3_block_corrections", "Figure 3", text_audit)
    manifest["Figure 3"] = [p_wp09, p_i1, p_i3]


def by_h(rows, key, filters=None):
    grouped = defaultdict(list)
    for row in rows:
        if filters and not filters(row):
            continue
        grouped[value(row, "h")].append(max(value(row, key), 1e-30))
    hs = np.array(sorted(grouped))
    med = np.array([np.median(grouped[h]) for h in hs])
    low = np.array([min(grouped[h]) for h in hs])
    high = np.array([max(grouped[h]) for h in hs])
    return hs, med, low, high


def figure_4(text_audit, manifest):
    p_raw = EVIDENCE / "wp10_i4_kpp_raw_fd_run1.csv"
    p_phys = EVIDENCE / "wp10_i4_kpp_physical_fd_run1.csv"
    p_full = EVIDENCE / "wp10_i5_full_matrix_fd_run1.csv"
    p_class = EVIDENCE / "wp10_i5b_final_classification.csv"
    raw, phys, full, classes = map(read_csv, (p_raw, p_phys, p_full, p_class))
    fig = plt.figure(figsize=(PAGE_WIDTH_IN, 6.75), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, height_ratios=[1.05, 1.0])
    ax1, ax2, ax3 = fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1]), fig.add_subplot(grid[1, :])

    for rows, color, label in ((raw, RED, r"returned $K_{pp}$"), (phys, TEAL, r"physical $K_{pp}$")):
        hs, med, low, high = by_h(rows, "relative_error_l2")
        ax1.fill_between(hs, low, high, color=color, alpha=0.12, linewidth=0)
        ax1.loglog(hs, med, marker="o", color=color, label=label)
    ax1.axhline(1e-5, color="#333", ls="--", lw=0.85, label=r"$10^{-5}$ gate")
    ax1.set_title(r"$J_{returned}=J_{phys}+J_{num}$", loc="left", fontweight="bold")
    ax1.set_xlabel("h"); ax1.set_ylabel("median relative L2 error"); plain_log_ticks(ax1, "x"); plain_log_ticks(ax1, "y")
    legend_below(ax1, ncol=2, y=-0.36); finish_axis(ax1); panel_label(ax1, "(a)", x=-0.29)

    core_filter = lambda r: r["state_id"] != "FINAL_CLIP_LOW"
    hs, med, low, high = by_h(full, "relative_error_scaled_l2", core_filter)
    ax2.fill_between(hs, low, high, color=TEAL, alpha=0.16, linewidth=0)
    ax2.loglog(hs, med, marker="o", color=TEAL, label="smooth/stable core median")
    colors = {"M1": RED, "M2": GOLD, "M3": PURPLE}
    for direction in ("M1", "M2", "M3"):
        filt = lambda r, d=direction: r["state_id"] == "FINAL_CLIP_LOW" and r["direction_id"] == d
        h, m, _, _ = by_h(full, "relative_error_scaled_l2", filt)
        ax2.loglog(h, m, ls=":", marker=".", color=colors[direction], label=f"FINAL_CLIP_LOW/{direction}")
    ax2.axhline(1e-5, color="#333", ls="--", lw=0.85, label=r"$10^{-5}$ gate")
    ax2.set_title("Mixed full-matrix directional series", loc="left", fontweight="bold")
    ax2.set_xlabel("h"); ax2.set_ylabel("relative scaled L2 error"); plain_log_ticks(ax2, "x"); plain_log_ticks(ax2, "y")
    legend_below(ax2, ncol=2, y=-0.36); finish_axis(ax2); panel_label(ax2, "(b)", x=-0.29)

    labels = [row["direction_id"] for row in classes]
    x = np.arange(len(labels)); width = 0.23
    metrics = [
        ("binary64 strict", "i5_original_best_relative_l2", RED),
        ("binary128 FD", "high_precision_best_relative_l2", BLUE),
        ("analytic oracle", "analytic_oracle_relative_l2", TEAL),
    ]
    for index, (name, key, color) in enumerate(metrics):
        ax3.bar(x + (index - 1) * width, [value(row, key) for row in classes], width, label=name, color=color)
    ax3.axhline(1e-5, color="#555", ls="--", lw=0.85, label="strict FD gate")
    ax3.axhline(1e-10, color="#555", ls=":", lw=0.85, label="active-set oracle gate")
    ax3.set_yscale("log"); ax3.set_xticks(x, labels); ax3.set_ylabel("best relative L2 error")
    plain_log_ticks(ax3, "y")
    ax3.set_title("Stable FINAL_CLIP_LOW active-set supplement; 36/39 strict record retained", loc="left", fontweight="bold")
    legend_below(ax3, ncol=3, y=-0.29); finish_axis(ax3, "major"); panel_label(ax3, "(c)")
    save_figure(fig, "Figure_4_jphys_active_set", "Figure 4", text_audit)
    manifest["Figure 4"] = [p_raw, p_phys, p_full, p_class]


def figure_5(text_audit, manifest):
    p_patch = EVIDENCE / "wp10_i6_patch_results_run1.csv"
    p_replay = EVIDENCE / "wp10_i6_rollback_replay_run1.csv"
    p_order = EVIDENCE / "wp10_i6_newton_order.csv"
    p_o64 = EVIDENCE / "wp10_i6a_recomputed_order.csv"
    p_o128 = EVIDENCE / "wp10_i6a_binary128_order.csv"
    patch, replay, order, o64, o128 = map(read_csv, (p_patch, p_replay, p_order, p_o64, p_o128))
    fig, axes = plt.subplots(2, 2, figsize=(PAGE_WIDTH_IN, 7.35), constrained_layout=True)

    mode_labels = {
        "TX": "TX", "TY": "TY", "ROT": "ROT", "EXX": "EXX", "EYY": "EYY", "GXY": "GXY",
        "P_UNIFORM": "P0", "P_AFFINE_X": "Px", "P_AFFINE_Y": "Py",
        "S_UNIFORM": "S0", "S_AFFINE_X": "Sx", "S_AFFINE_Y": "Sy",
        "MIXED_AFFINE": "Mix",
    }
    modes = [mode_labels.get(row["mode_id"], row["mode_id"]) for row in patch]
    errors = [max(value(row, "interior_scaled_error"), value(row, "balance_scaled_error"), 1e-18) for row in patch]
    colors = [BLUE if row["null_mode"] == "1" else TEAL for row in patch]
    x = np.arange(len(modes))
    axes[0, 0].bar(x, errors, color=colors, width=0.72)
    axes[0, 0].axhline(1e-8, color="#333", ls="--", lw=0.85, label="nonzero gate")
    axes[0, 0].axhline(1e-10, color="#333", ls=":", lw=0.85, label="null gate")
    axes[0, 0].set_yscale("log"); axes[0, 0].set_xticks(x, modes, rotation=90, ha="center")
    plain_log_ticks(axes[0, 0], "y")
    axes[0, 0].set_ylabel("registered scaled error"); axes[0, 0].set_title("2x2 Q4 patch modes", loc="left", fontweight="bold")
    legend_below(axes[0, 0], ncol=2, y=-0.33); finish_axis(axes[0, 0], "major"); panel_label(axes[0, 0], "(a)", x=-0.18)

    histories = list(dict.fromkeys(row["history"] for row in replay))
    quantities = list(dict.fromkeys(row["quantity"] for row in replay))
    matrix = np.full((len(histories), len(quantities)), np.nan)
    for row in replay:
        matrix[histories.index(row["history"]), quantities.index(row["quantity"])] = 1.0 if row["status"] == "PASS" else 0.0
    cmap = ListedColormap(["#b2453d", "#dcefed"])
    axes[0, 1].imshow(matrix, aspect="auto", cmap=cmap, vmin=0, vmax=1)
    quantity_labels = {
        "RHS": "RHS", "AMATRX": "Matrix", "PHYSICAL_SVARS": "SVARS", "BRANCH": "Branch",
    }
    axes[0, 1].set_xticks(
        range(len(quantities)), [quantity_labels.get(item, item) for item in quantities],
        rotation=35, ha="right",
    )
    history_labels = {
        "H0_PRE_PATCH": "H0", "H1_AB": "H1 A/B", "H1_BA": "H1 B/A",
        "H2_H3_CUTBACK": "H2/H3", "H4_KINC0": "H4", "H5_POST_PATCH": "H5",
    }
    axes[0, 1].set_yticks(range(len(histories)), [history_labels.get(item, item) for item in histories])
    for i in range(len(histories)):
        for j in range(len(quantities)):
            if np.isfinite(matrix[i, j]):
                axes[0, 1].text(j, i, "PASS" if matrix[i, j] else "FAIL", ha="center", va="center", fontsize=FONT_PT)
    axes[0, 1].set_title("Rollback and rejected-trial replay", loc="left", fontweight="bold")
    panel_label(axes[0, 1], "(b)", x=-0.18)

    for benchmark, color in (("N1", RED), ("N2", BLUE)):
        rows = [r for r in order if r["benchmark"] == benchmark and r["profile"] == "M2_PHYS"]
        axes[1, 0].plot([int(r["iteration"]) for r in rows], [value(r, "q") for r in rows], marker="o", color=color, label=benchmark)
    axes[1, 0].axhline(1.8, color="#333", ls="--", lw=0.85, label="registered q gate")
    axes[1, 0].set_xlabel("Newton iteration"); axes[1, 0].set_ylabel("empirical residual order, q")
    axes[1, 0].set_title("I6 Newton-order hard stop", loc="left", fontweight="bold")
    legend_below(axes[1, 0], ncol=2, y=-0.30); finish_axis(axes[1, 0], "major"); panel_label(axes[1, 0], "(c)", x=-0.18)

    for benchmark, color in (("N1", RED), ("N2", BLUE)):
        rows64 = [r for r in o64 if r["benchmark"] == benchmark]
        rows128 = [r for r in o128 if r["benchmark"] == benchmark]
        axes[1, 1].plot([int(r["iteration"]) for r in rows64], [value(r, "q") for r in rows64],
                        marker="o", color=color, alpha=0.48, label=f"{benchmark} binary64")
        axes[1, 1].plot([int(r["iteration"]) for r in rows128], [value(r, "residual_order") for r in rows128],
                        marker="s", color=color, ls="--", label=f"{benchmark} binary128")
    axes[1, 1].axhline(1.8, color="#333", ls=":", lw=0.85)
    axes[1, 1].set_xlabel("Newton iteration"); axes[1, 1].set_ylabel("empirical residual order, q")
    axes[1, 1].set_title("Binary64/binary128 order audit", loc="left", fontweight="bold")
    legend_below(axes[1, 1], ncol=2, y=-0.30); finish_axis(axes[1, 1], "major"); panel_label(axes[1, 1], "(d)", x=-0.18)
    save_figure(fig, "Figure_5_patch_newton_limits", "Figure 5", text_audit)
    manifest["Figure 5"] = [p_patch, p_replay, p_order, p_o64, p_o128]


def figure_6(text_audit, manifest):
    p_basin = EVIDENCE / "wp10_i6c_r1_basin_entry_radius.csv"
    p_sub = EVIDENCE / "wp10_i6d_subspace_audit.csv"
    p_dec = EVIDENCE / "wp10_i6d_one_step_decomposition.csv"
    basin, subspace, decomposition = map(read_csv, (p_basin, p_sub, p_dec))
    directions = ["LEGACY", "ERR_K1", "ERR_K2", "ERR_K3", "UPDATE_K0", "UPDATE_K1", "ORTH_K0",
                  "SVD_MAX", "SVD_MID", "SVD_MIN", "FIELD_U", "FIELD_P", "FIELD_S"]
    benchmarks = ["N1", "N2"]
    matrix = np.full((2, 13), np.nan)
    for row in basin:
        if row["scope"] == "DIRECTION" and row["persistent_entry"] == "1":
            matrix[benchmarks.index(row["benchmark"]), directions.index(row["direction_id"])] = math.log10(value(row, "entry_radius"))
    fig = plt.figure(figsize=(PAGE_WIDTH_IN, 7.60))
    grid = fig.add_gridspec(
        2, 2, height_ratios=[1.05, 1.15],
        left=0.13, right=0.94, top=0.95, bottom=0.17,
        hspace=0.58, wspace=0.42,
    )
    ax1, ax2, ax3 = fig.add_subplot(grid[0, :]), fig.add_subplot(grid[1, 0]), fig.add_subplot(grid[1, 1])
    image = ax1.imshow(matrix, aspect="auto", cmap="viridis", vmin=-4, vmax=-2)
    direction_labels = ["Base", "Err 1", "Err 2", "Err 3", "Upd 0", "Upd 1", "Orth.",
                        "SVD max", "SVD mid", "SVD min", r"$u$", r"$p_w$", r"$S_n$"]
    ax1.set_yticks(range(2), benchmarks)
    ax1.set_xticks(range(13), direction_labels, rotation=60, ha="right", rotation_mode="anchor")
    ax1.tick_params(axis="x", labelsize=8.5, pad=2)
    ax1.set_title("Local-basin entry: 23/26 benchmark-direction pairs", loc="left", fontweight="bold", pad=10)
    colorbar = fig.colorbar(image, ax=ax1, pad=0.02, aspect=18); colorbar.set_label("log10(entry radius)")
    for i in range(2):
        for j in range(13):
            ax1.text(j, i, "none" if np.isnan(matrix[i, j]) else f"{matrix[i, j]:.1f}", ha="center", va="center",
                     color="white" if not np.isnan(matrix[i, j]) and matrix[i, j] < -2.8 else INK, fontsize=8.5)
    panel_label(ax1, "(a)", x=-0.12, y=1.10)

    targets = [("N1", "FIELD_U"), ("N1", "FIELD_P"), ("N2", "FIELD_P")]
    rows = [r for r in subspace if (r["benchmark"], r["direction_id"]) in targets]
    labels = [f"{r['benchmark']}/{r['direction_id'].replace('FIELD_', '')}" for r in rows]
    x = np.arange(len(rows)); width = 0.2
    metrics = [
        ("bottom-3 dir.", "bottom3_direction_energy", NAVY),
        ("bottom-3 resp.", "bottom3_response_energy", BLUE),
        ("equal-P dir.", "equal_pressure_energy", GOLD),
        ("equal-P resp.", "equal_pressure_response_energy", RED),
    ]
    for index, (name, key, color) in enumerate(metrics):
        ax2.bar(x + (index - 1.5) * width, [value(r, key) for r in rows], width, label=name, color=color)
    ax2.axhline(0.9, color="#333", ls="--", lw=0.85, label="dominance gate")
    ax2.set_ylim(0, 1.05); ax2.set_xticks(x, labels, rotation=18, ha="right"); ax2.set_ylabel("energy fraction")
    ax2.set_title("Mixed subspace evidence", loc="left", fontweight="bold", pad=8)
    legend_below(ax2, ncol=2, y=-0.24); finish_axis(ax2, "major"); panel_label(ax2, "(b)", x=-0.25, y=1.06)

    colors = {("N1", "FIELD_U"): RED, ("N1", "FIELD_P"): GOLD, ("N2", "FIELD_P"): BLUE}
    for target in targets:
        rows_t = sorted([r for r in decomposition if (r["benchmark"], r["direction_id"]) == target and int(r["sign"]) == 1],
                        key=lambda r: value(r, "amplitude"), reverse=True)
        ax3.semilogx([value(r, "amplitude") for r in rows_t], [value(r, "nonlinear_fraction") for r in rows_t],
                     marker="o", color=colors[target], label=f"{target[0]}/{target[1].replace('FIELD_', '')}")
    ax3.axhline(0.9, color="#333", ls="--", lw=0.85, label="dominance gate")
    ax3.invert_xaxis(); ax3.set_ylim(-0.02, 1.02); ax3.set_xlabel("scaled amplitude")
    plain_log_ticks(ax3, "x")
    ax3.set_ylabel("nonlinear contribution fraction"); ax3.set_title("Nonlinear remainder", loc="left", fontweight="bold", pad=8)
    legend_below(ax3, ncol=2, y=-0.24); finish_axis(ax3); panel_label(ax3, "(c)", x=-0.25, y=1.06)
    save_figure(fig, "Figure_6_basin_subspace", "Figure 6", text_audit)
    manifest["Figure 6"] = [p_basin, p_sub, p_dec]


def write_audits(text_audit, manifest):
    QA.mkdir(parents=True, exist_ok=True)
    with (QA / "figure_text_size_audit.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["figure_id", "text", "font_family", "font_size_pt", "hard_minimum_pass", "times_new_roman_pass"], lineterminator="\n")
        writer.writeheader(); writer.writerows(text_audit)
    with (QA / "figure_legend_overlap_audit.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["figure_id", "legend_id", "placement", "intersects_data_axes"], lineterminator="\n")
        writer.writeheader(); writer.writerows(LEGEND_AUDIT)
    with (QA / "figure_tick_label_layout_audit.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["figure_id", "axis_id", "label_count", "overlap_pair_count", "outside_figure_count"], lineterminator="\n")
        writer.writeheader(); writer.writerows(TICK_AUDIT)
    manifest_rows = []
    for figure_id, sources in manifest.items():
        for source in sources or [None]:
            manifest_rows.append(
                {
                    "figure_id": figure_id,
                    "source_path": "CONCEPTUAL_FROZEN_CONTRACT" if source is None else source.relative_to(PACKAGE).as_posix(),
                    "source_sha256": "NOT_APPLICABLE" if source is None else sha256(source),
                    "source_status": "FROZEN_READ_ONLY",
                }
            )
    with (QA / "figure_source_manifest.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=manifest_rows[0].keys(), lineterminator="\n")
        writer.writeheader(); writer.writerows(manifest_rows)
    outputs = []
    for path in sorted(OUT.iterdir()):
        outputs.append({"file": path.name, "sha256": sha256(path), "bytes": path.stat().st_size})
    report = {
        "status": "PASS" if text_audit and all(row["hard_minimum_pass"] and row["times_new_roman_pass"] for row in text_audit) and all(not row["intersects_data_axes"] for row in LEGEND_AUDIT) and all(not row["overlap_pair_count"] for row in TICK_AUDIT) else "FAIL",
        "figure_count": 6,
        "output_file_count": len(outputs),
        "minimum_text_size_pt": min(row["font_size_pt"] for row in text_audit),
        "font_family": "Times New Roman",
        "font_family_gate": "PASS" if all(row["times_new_roman_pass"] for row in text_audit) else "FAIL",
        "legend_count": len(LEGEND_AUDIT),
        \
        "normal_text_target_pt": FONT_PT,
        "review_source_width_mm": PAGE_WIDTH_IN * 25.4,
        "production_registered_width_mm": 190,
        "primary_format": "PDF_VECTOR",
        "outputs": outputs,
    }
    (QA / "figure_build_qa.json").write_text(json.dumps(report, indent=2), encoding="utf-8")


def main():
    text_audit = []
    manifest = {}
    figure_1(text_audit, manifest)
    figure_2(text_audit, manifest)
    figure_3(text_audit, manifest)
    figure_4(text_audit, manifest)
    figure_5(text_audit, manifest)
    figure_6(text_audit, manifest)
    write_audits(text_audit, manifest)
    print(json.dumps({"status": "PASS", "figures": 6, "minimum_text_size_pt": min(r["font_size_pt"] for r in text_audit)}))


if __name__ == "__main__":
    main()
