
"""
2D SPH potential-flow velocity-gradient L2 error plots.

Standard uncorrected SPH gradient:
    grad(A_i) = sum_j V_j * A_j * grad(W_ij)

Analytical potential flow around a circular cylinder.
No time-dependent simulation.

Produces both:
    1. Semi-log plots (logarithmic x-axis, linear y-axis)
    2. Log-log plots (logarithmic x-axis and y-axis)

Run:
    python plot_potential_flow_gradient.py SmoothingLengthL2.csv

Optional:
    python plot_potential_flow_gradient.py SmoothingLengthL2.csv --error L2

    python plot_potential_flow_gradient.py SmoothingLengthL2.csv --error L2_area
"""

import argparse
import csv
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

from matplotlib.ticker import (
    LogLocator,
    LogFormatterMathtext,
    FuncFormatter,
    MaxNLocator,
    ScalarFormatter,
)


# ====================================================
# 1. Logarithmic x-axis tick labels
# ====================================================

def tick_label(value, position):

    if value <= 0 or not math.isfinite(value):
        return ''

    exponent = math.floor(math.log10(value))
    mantissa = value / 10**exponent

    if not any(
        abs(mantissa - target) < 1e-6
        for target in (1, 2, 5)
    ):
        return ''

    if 0.001 <= value < 10000:
        return f'{value:g}'

    return rf'${value:.0e}$'


# ====================================================
# 2. Load potential-flow L2 results
# ====================================================

def load_results(path):

    results = []

    with path.open(
        newline='',
        encoding='utf-8-sig'
    ) as handle:

        reader = csv.DictReader(handle)

        required = {
            'kernel',
            'ch',
            'field',
            'method',
            'h_over_D',
            'delta_over_h',
            'L2',
            'L2_area'
        }

        missing = required - set(reader.fieldnames or [])

        if missing:
            raise ValueError(
                f'Missing CSV columns: {sorted(missing)}'
            )

        for row in reader:

            if row['field'].strip().lower() != 'potential_flow':
                continue

            if row['method'].strip().lower() != 'standard':
                continue

            try:

                ch = float(row['ch'])
                h_over_D = float(row['h_over_D'])
                dp_over_h = float(row['delta_over_h'])
                L2 = float(row['L2'])
                L2_area = float(row['L2_area'])

            except (ValueError, TypeError):
                continue

            values = (
                ch,
                h_over_D,
                dp_over_h,
                L2,
                L2_area
            )

            if not all(math.isfinite(v) for v in values):
                continue

            if ch <= 0 or h_over_D <= 0 or dp_over_h <= 0:
                continue

            if L2 < 0 or L2_area < 0:
                continue

            # results.append({
            #     'kernel': row['kernel'].strip().lower(),
            #     'ch': ch,
            #     'h_over_D': h_over_D,
            #     'dp_over_h': dp_over_h,
            #     'L2': L2,
            #     'L2_area': L2_area
            # })
            
            results.append({
                'kernel': row['kernel'].strip().lower(),
                'ch': ch,
                'h_over_D': h_over_D,
                'h_over_dp': 1.0 / dp_over_h,
                'L2': L2,
                'L2_area': L2_area
            })

    if not results:
        raise ValueError(
            'No potential-flow standard SPH results found. '
            'Check the CSV field and method columns.'
        )

    return results


# ====================================================
# 3. Research-paper plotting style
# ====================================================

def paper_style():

    plt.rcParams.update({

        'font.family': 'serif',
        'font.size': 10,

        'axes.labelsize': 11,
        'axes.titlesize': 11,
        'axes.linewidth': 0.9,

        'xtick.labelsize': 9,
        'ytick.labelsize': 9,

        'xtick.direction': 'in',
        'ytick.direction': 'in',

        'xtick.top': True,
        'ytick.right': True,

        'xtick.major.size': 5,
        'ytick.major.size': 5,

        'xtick.minor.size': 2.5,
        'ytick.minor.size': 2.5,

        'legend.fontsize': 9,

        'pdf.fonttype': 42,
        'savefig.dpi': 350,

    })


# ====================================================
# 4. Format axes
# ====================================================

def format_axis(ax, xlabel, errors, log_y=False):

    # ====================================================
    # X-axis: logarithmic for both plot types
    # ====================================================

    ax.set_xscale('log')

    ax.xaxis.set_major_locator(
        LogLocator(
            base=10,
            subs=(1.0, 2.0, 5.0),
            numticks=40
        )
    )

    ax.xaxis.set_major_formatter(
        FuncFormatter(tick_label)
    )

    ax.xaxis.set_minor_locator(
        LogLocator(
            base=10,
            subs=(3, 4, 6, 7, 8, 9),
            numticks=100
        )
    )

    # ====================================================
    # Y-axis: linear or logarithmic
    # ====================================================

    if log_y:

        positive_errors = [
            error
            for error in errors
            if error > 0
        ]

        if not positive_errors:
            raise ValueError(
                'No positive L2 errors available '
                'for logarithmic plotting.'
            )

        ymin = min(positive_errors)
        ymax = max(positive_errors)

        ax.set_yscale('log')

        ax.set_ylim(
            ymin / 1.2,
            ymax * 1.2
        )

        ax.yaxis.set_major_locator(
            LogLocator(
                base=10,
                numticks=12
            )
        )

        ax.yaxis.set_major_formatter(
            LogFormatterMathtext(base=10)
        )

        ax.yaxis.set_minor_locator(
            LogLocator(
                base=10,
                subs=tuple(range(2, 10)),
                numticks=100
            )
        )

    else:

        ymin = min(errors)
        ymax = max(errors)

        ax.set_yscale('linear')

        if ymax > ymin:

            padding = 0.08 * (ymax - ymin)

            ax.set_ylim(
                max(0.0, ymin - padding),
                ymax + padding
            )

        else:

            padding = max(
                0.05 * abs(ymin),
                1e-6
            )

            ax.set_ylim(
                max(0.0, ymin - padding),
                ymax + padding
            )

        ax.yaxis.set_major_locator(
            MaxNLocator(nbins=8)
        )

        ax.yaxis.set_major_formatter(
            ScalarFormatter(useOffset=False)
        )

        ax.ticklabel_format(
            axis='y',
            style='plain',
            useOffset=False
        )

    # ====================================================
    # Labels and grid
    # ====================================================

    ax.set_xlabel(xlabel)

    ax.set_ylabel(
        r'$L_2$ velocity-gradient error (s$^{-1}$)'
    )

    ax.grid(
        True,
        which='major',
        color='0.85',
        linewidth=0.5,
        alpha=0.6
    )

    ax.grid(
        False,
        which='minor'
    )


# ====================================================
# 5. Plot potential-flow gradient L2 error
# ====================================================

def plot_results(
    results,
    folder,
    prefix,
    x_key,
    xlabel,
    error_key,
    log_y=False
):

    fig, ax = plt.subplots(
        figsize=(5.5, 4.0)
    )

    kernels = sorted({
        row['kernel']
        for row in results
    })

    colors = [
        'black',
        '0.4',
        '0.65',
        '0.25'
    ]

    styles = [
        '-',
        '--',
        '-.',
        ':'
    ]

    all_errors = []

    for index, kernel in enumerate(kernels):

        selected = [
            row
            for row in results
            if row['kernel'] == kernel
        ]

        # Retain the final result for duplicate x-values.
        by_x = {
            row[x_key]: row[error_key]
            for row in selected
        }

        xs = sorted(by_x)

        if log_y:

            # Logarithmic axes cannot display zero or negative errors.
            xs = [
                x for x in xs
                if by_x[x] > 0
            ]

        ys = [
            by_x[x]
            for x in xs
        ]

        if not xs:
            continue

        all_errors.extend(ys)

        ax.plot(
            xs,
            ys,

            color=colors[index % len(colors)],
            linestyle=styles[index % len(styles)],

            linewidth=1.3,

            marker='o',
            markersize=3.5,

            markerfacecolor='white',
            markeredgewidth=0.8,

            label=kernel.capitalize()
        )

    if not all_errors:
        plt.close(fig)
        print(
            f'Skipped {prefix}: no valid data '
            f'for this axis scale.'
        )
        return

    format_axis(
        ax,
        xlabel,
        all_errors,
        log_y=log_y
    )

    ax.set_title(
        'Potential-flow velocity-gradient error',
        pad=8
    )

    ax.legend(
        loc='best',
        frameon=False,
        handlelength=2.5
    )

    fig.tight_layout()

    # ====================================================
    # Save PNG and PDF
    # ====================================================

    scale_suffix = (
        '_loglog'
        if log_y
        else '_semilogx'
    )

    for ext in ('.png', '.pdf'):

        target = (
            folder /
            f'{prefix}{scale_suffix}{ext}'
        )

        fig.savefig(
            target,
            bbox_inches='tight',
            facecolor='white'
        )

        print('Saved', target)

    plt.close(fig)


# ====================================================
# 6. Main
# ====================================================

def main():

    parser = argparse.ArgumentParser(
        description=__doc__
    )

    parser.add_argument(
        'csv',
        nargs='?',
        default='SmoothingLengthL2.csv',
        help='CSV output from the C++ potential-flow gradient test'
    )

    parser.add_argument(
        '--output',
        help='Output folder; default: same as CSV'
    )

    parser.add_argument(
        '--error',
        choices=[
            'L2',
            'L2_area'
        ],
        default='L2_area',
        help='Error metric to plot'
    )

    args = parser.parse_args()

    csv_path = Path(args.csv)

    results = load_results(csv_path)

    paper_style()

    folder = (
        Path(args.output)
        if args.output
        else csv_path.parent
    )

    folder.mkdir(
        parents=True,
        exist_ok=True
    )

    # ====================================================
    # Generate both semi-log and log-log figures
    # ====================================================

    for log_y in (False, True):

        # L2 versus h/D

        plot_results(
            results,
            folder,
            'L2_vs_h_over_D',
            'h_over_D',
            r'$h/D$',
            args.error,
            log_y=log_y
        )
        
        # L2 versus dp/h

        # plot_results(
        #     results,
        #     folder,
        #     'L2_vs_dp_over_h',
        #     'dp_over_h',
        #     r'$dp/h$',
        #     args.error,
        #     log_y=log_y
        # )
        
        plot_results(
            results,
            folder,
            'L2_vs_h_over_dp',
            'h_over_dp',
            r'$h/dp$',
            args.error,
            log_y=log_y
        )

    print(
        'Finished: potential-flow gradient L2 plots generated.'
    )


if __name__ == '__main__':
    main()
