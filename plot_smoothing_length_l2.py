"""Journal-style plots of regional SPH potential-flow velocity-gradient errors.

Place next to SmoothingLengthL2.csv and run:
    python plot_potential_flow_regions_journal.py

Creates six figures (linear, semilogx, loglog for h/dp and mean h/D),
each saved as PDF, PNG, and SVG in the same folder as this script.
"""

import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, LogLocator, NullFormatter, ScalarFormatter

FOLDER = Path(__file__).resolve().parent
CSV_FILE = FOLDER / "SmoothingLengthL2.csv"

# Disjoint regions first, then overlapping aggregates.
REGIONS = {
    "wall_0.5D_1D": r"$0.5D\leq r<D$",
    "near_1D_2D": r"$D\leq r<2D$",
    "outer_2D_10D": r"$2D\leq r\leq10D$",
    "excluding_wall_1D_10D": r"$D\leq r\leq10D$",
    "full_0.5D_10D": r"$0.5D\leq r\leq10D$",
}

# Five curves remain identifiable in monochrome printing.
STYLES = {
    "wall_0.5D_1D": dict(color="0.05", linestyle="-", marker="o"),
    "near_1D_2D": dict(color="0.30", linestyle="--", marker="s"),
    "outer_2D_10D": dict(color="0.48", linestyle="-.", marker="^"),
    "excluding_wall_1D_10D": dict(color="0.25", linestyle=":", marker="D"),
    "full_0.5D_10D": dict(color="0.60", linestyle="-", marker="v"),
}


def load_results():
    if not CSV_FILE.is_file():
        raise FileNotFoundError(
            f"Missing CSV: {CSV_FILE}\n"
            "Put this script next to the newly generated SmoothingLengthL2.csv."
        )

    data = {region: {} for region in REGIONS}
    with CSV_FILE.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        required = {"ch", "field", "method", "region", "h_over_D", "L2_area"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"CSV is missing {sorted(missing)}.\n"
                "The C++ summary header must include the 'region' column."
            )

        for line, row in enumerate(reader, start=2):
            region = (row.get("region") or "").strip()
            if region not in data:
                continue
            if (row.get("field") or "").strip().lower() != "potential_flow":
                continue
            if (row.get("method") or "").strip().lower() != "standard":
                continue
            try:
                ch = float(row["ch"])
                hd = float(row["h_over_D"])
                err = float(row["L2_area"])
            except (ValueError, TypeError) as exc:
                raise ValueError(f"Invalid numerical value on CSV line {line}") from exc
            if not all(math.isfinite(v) for v in (ch, hd, err)) or ch <= 0 or hd <= 0 or err < 0:
                raise ValueError(f"Invalid data on CSV line {line}")
            # When the CSV contains repeated runs, retain the last entry for each ch.
            data[region][ch] = {"h_over_dp": ch, "h_over_D": hd, "error": err}

    if not any(data.values()):
        raise ValueError("No potential_flow / standard region results found in the CSV.")
    return data


def set_style():
    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 9,
        "mathtext.fontset": "stix",
        "axes.labelsize": 10,
        "axes.linewidth": 0.8,
        "axes.spines.top": True,
        "axes.spines.right": True,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.major.size": 4,
        "ytick.major.size": 4,
        "xtick.minor.size": 2,
        "ytick.minor.size": 2,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.minor.width": 0.6,
        "ytick.minor.width": 0.6,
        "legend.fontsize": 8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "savefig.dpi": 400,
    })


def format_axes(ax, scale, x_label):
    ax.set_xlabel(x_label)
    ax.set_ylabel(r"$L_{2,\nabla\mathbf{u}}\; (\mathrm{s}^{-1})$")
    ax.set_xscale("log" if scale in ("semilogx", "loglog") else "linear")
    ax.set_yscale("log" if scale == "loglog" else "linear")

    if scale in ("semilogx", "loglog"):
        ax.xaxis.set_minor_locator(LogLocator(base=10, subs=range(2, 10)))
        ax.xaxis.set_minor_formatter(NullFormatter())
    else:
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))

    if scale == "loglog":
        ax.yaxis.set_minor_locator(LogLocator(base=10, subs=range(2, 10)))
        ax.yaxis.set_minor_formatter(NullFormatter())
    else:
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        sf = ScalarFormatter(useOffset=False, useMathText=True)
        sf.set_powerlimits((-3, 3))
        ax.yaxis.set_major_formatter(sf)
        ax.set_ylim(bottom=0)

    # Subtle horizontal guides; no dense full grid.
    ax.grid(axis="y", which="major", color="0.88", linewidth=0.45)
    ax.set_axisbelow(True)
    ax.margins(x=0.035)


def plot_one(data, x_key, x_label, scale):
    fig, ax = plt.subplots(figsize=(6.6, 3.65))
    any_curve = False

    for region, label in REGIONS.items():
        rows = sorted(data[region].values(), key=lambda row: row[x_key])
        if scale == "loglog":
            rows = [row for row in rows if row["error"] > 0]
        if not rows:
            continue
        ax.plot(
            [row[x_key] for row in rows],
            [row["error"] for row in rows],
            label=label,
            linewidth=1.15,
            markersize=3.4,
            markerfacecolor="white",
            markeredgewidth=0.85,
            markevery=max(1, len(rows) // 12),
            **STYLES[region],
        )
        any_curve = True

    if not any_curve:
        plt.close(fig)
        print(f"Skipped {x_key} {scale}: no valid values")
        return

    format_axes(ax, scale, x_label)
    # A right-side legend keeps overlapping results visible.
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5),
              frameon=False, handlelength=2.5, labelspacing=0.85,
              borderaxespad=0)
    fig.subplots_adjust(left=0.13, right=0.68, bottom=0.17, top=0.96)

    basename = f"Nr200_L2_regions_vs_{x_key}_{scale}"
    for ext in ("png", "pdf", "svg"):
        output = FOLDER / f"{basename}.{ext}"
        fig.savefig(output, bbox_inches="tight", facecolor="white")
        print("Saved:", output)
    plt.close(fig)


def main():
    print("Reading:", CSV_FILE)
    data = load_results()
    for region in REGIONS:
        print(f"  {region}: {len(data[region])} smoothing lengths")
    set_style()
    for x_key, x_label in (
        ("h_over_dp", r"$h/\Delta p$"),
        ("h_over_D", r"$\overline{h}/D$"),
    ):
        for scale in ("linear", "semilogx", "loglog"):
            plot_one(data, x_key, x_label, scale)
    print("\nAll plots saved in:", FOLDER)


if __name__ == "__main__":
    main()
