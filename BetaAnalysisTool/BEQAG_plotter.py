import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse

def plot_with_uncertainties(ax, df, x_col, y_col, yerr_col, group_col='Channel'):
    groups = df.groupby(group_col)
    markers = ['^', 'v']
    colours = ['dodgerblue', 'limegreen']
    for i, (name, group) in enumerate(groups):
        marker = markers[i % len(markers)]
        colour = colours[i % len(colours)]
        if (group[x_col] == 0).all(): continue
        ax.errorbar(
            group[x_col], group[y_col], yerr=group[yerr_col],
            fmt=marker,
            markersize=9,
            elinewidth=1.5,
            capsize=4,
            markerfacecolor=colour,
            markeredgewidth=1.2,
            markeredgecolor='black',
            ecolor=colour,
            label=name
        )
    if x_col == 'Bias':
        x_col = 'Bias / V'
    if x_col == 'Gain by Log Fluence':
        x_col = r'G$\ast$ln($\Phi$)'
    ax.set_xlabel(x_col, fontsize = 14)
    ax.set_ylabel(y_col, fontsize = 14)
    ax.grid(True)
    ax.legend(fontsize=14, loc='upper right')

def plot_without_uncertainties(ax, df, x_col, y_col, group_col='Channel'):
    groups = df.groupby(group_col)
    markers = ['^', 'v']
    colours = ['dodgerblue', 'limegreen']
    for i, (name, group) in enumerate(groups):
        marker = markers[i % len(markers)]
        colour = colours[i % len(colours)]
        ax.scatter(
            group[x_col], group[y_col],
            marker=marker,
            s=81,
            facecolors=colour,
            edgecolors='black',
            linewidths=1.2,
            label=name
        )
    if x_col == 'Bias':
        x_col = 'Bias / V'
    ax.set_xlabel(x_col, fontsize = 14)
    ax.set_ylabel(y_col, fontsize = 14)
    ax.grid(True)
    ax.legend(fontsize=14, loc='upper right')

def annotate_plot(ax, filename, subtitle):
    base_text = filename.split('W')[0] if 'W' in filename else filename
    ax.text(0.05, 0.95, base_text, transform=ax.transAxes,
            fontsize=20, fontweight='bold', va='top')
    ax.text(0.05, 0.87, subtitle, transform=ax.transAxes,
            fontsize=14, style='italic', va='top')

def compute_gain_by_log_fluence(row, fluence_map):
    fluence = fluence_map.get(row['Channel'], 0)
    return row['Gain'] * np.log(fluence) if fluence > 0 else 0

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_file", help="Path to the CSV file")
    parser.add_argument("subtitle", help="Subtitle for the plots")
    parser.add_argument("--fluences", help="Comma-separated list of fluences", default=None)
    args = parser.parse_args()

    df = pd.read_csv(args.csv_file)
    base_name = os.path.splitext(os.path.basename(args.csv_file))[0]

    if args.fluences:
        fluence_values = [float(f) for f in args.fluences.split(',')]
        channels = df['Channel'].unique()
        assert len(fluence_values) == len(channels), "Number of fluences must match number of unique channels"
        fluence_map = dict(zip(channels, fluence_values))
        df['Gain by Log Fluence'] = df.apply(lambda row: compute_gain_by_log_fluence(row, fluence_map), axis=1)
    else:
        df['Gain by Log Fluence'] = np.nan

    columns_to_plot = [
        [('Bias', 'Amplitude / mV', None), ('Bias', 'Rise time / ps', 'Rise time Unc / ps'), ('Bias', 'RMS Noise / mV', 'RMS Unc / mV')],
        [('E field / V/cm', 'Amplitude / mV', None), ('E field / V/cm', 'Rise time / ps', 'Rise time Unc / ps'), ('E field / V/cm', 'RMS Noise / mV', 'RMS Unc / mV')],
        [('Gain', 'Amplitude / mV', None), ('Gain', 'Rise time / ps', 'Rise time Unc / ps'), ('Gain', 'RMS Noise / mV', 'RMS Unc / mV')],
        [('Bias', 'Charge / fC', None), ('Bias', 'Frac Charge >1.5xMPV', None), ('Bias', 'Frac Charge >1.5xQmax', None)],
        [('E field / V/cm', 'Charge / fC', None), ('E field / V/cm', 'Frac Charge >1.5xMPV', None), ('E field / V/cm', 'Frac Charge >1.5xQmax', None)],
        [('Gain', 'Charge / fC', None), ('Gain', 'Frac Charge >1.5xMPV', None), ('Gain', 'Frac Charge >1.5xQmax', None)],
        [('Bias', 'Gain', None), ('Bias', 'TR @ 30% / ps', 'TR Unc @ 30% / ps'), ('Bias', 'TR @ 50% / ps', 'TR Unc @ 50% / ps')],
        [('E field / V/cm', 'Gain', None), ('E field / V/cm', 'TR @ 30% / ps', 'TR Unc @ 30% / ps'), ('E field / V/cm', 'TR @ 50% / ps', 'TR Unc @ 50% / ps')],
        [('Gain by Log Fluence', 'RMS Noise / mV', 'RMS Unc / mV'), ('Charge / fC', 'TR @ 30% / ps', 'TR Unc @ 30% / ps'), ('Charge / fC', 'TR @ 50% / ps', 'TR Unc @ 50% / ps')],
        [('Amplitude / mV', 'Gain', None), ('Amplitude / mV', 'TR @ 30% / ps', 'TR Unc @ 30% / ps'), ('Amplitude / mV', 'TR @ 50% / ps', 'TR Unc @ 50% / ps')],
        [('Bias', 'Approx Jitter / ps', 'Approx Jitter Unc / ps'), ('Bias', 'Jitter / ps', 'Jitter Unc / ps'), ('Bias', 'Jitter[20%:80%] / ps', 'Jitter[20%:80%] Unc / ps')],
        [('E field / V/cm', 'Approx Jitter / ps', 'Approx Jitter Unc / ps'), ('E field / V/cm', 'Jitter / ps', 'Jitter Unc / ps'), ('E field / V/cm', 'Jitter[20%:80%] / ps', 'Jitter[20%:80%] Unc / ps')],
        [('Charge / fC', 'Approx Jitter / ps', 'Approx Jitter Unc / ps'), ('Charge / fC', 'Jitter / ps', 'Jitter Unc / ps'), ('Charge / fC', 'Jitter[20%:80%] / ps', 'Jitter[20%:80%] Unc / ps')],
        [('Amplitude / mV', 'Approx Jitter / ps', 'Approx Jitter Unc / ps'), ('Amplitude / mV', 'Jitter / ps', 'Jitter Unc / ps'), ('Amplitude / mV', 'Jitter[20%:80%] / ps', 'Jitter[20%:80%] Unc / ps')],
        [('Bias', 'Landau TR Cpt / ps', 'Landau TR Unc / ps'), ('Bias', 'WF6 Param / ps/um', 'WF6 Param Unc / ps/um'), ('Bias', 'PMAX low / mV', None)],
        [('E field / V/cm', 'Landau TR Cpt / ps', 'Landau TR Unc / ps'), ('E field / V/cm', 'WF6 Param / ps/um', 'WF6 Param Unc / ps/um'), ('E field / V/cm', 'PMAX low / mV', None)],
        [('Charge / fC', 'Landau TR Cpt / ps', 'Landau TR Unc / ps'), ('Charge / fC', 'WF6 Param / ps/um', 'WF6 Param Unc / ps/um'), ('Charge / fC', 'PMAX low / mV', None)],
        [('Amplitude / mV', 'Landau TR Cpt / ps', 'Landau TR Unc / ps'), ('Amplitude / mV', 'WF6 Param / ps/um', 'WF6 Param Unc / ps/um'), ('Amplitude / mV', 'PMAX low / mV', None)],
    ]



    output_dir = f"{base_name}_plots"
    os.makedirs(output_dir, exist_ok=True)

    for i, plot_group in enumerate(columns_to_plot, start=1):
        fig, axs = plt.subplots(1, 3, figsize=(18, 5), sharey=False)

        for ax, (x, y, yerr) in zip(axs, plot_group):
            if yerr is not None and yerr in df.columns:
                plot_with_uncertainties(ax, df, x, y, yerr)
            else:
                plot_without_uncertainties(ax, df, x, y)
            annotate_plot(ax, args.csv_file, args.subtitle)

        plt.tight_layout()
        save_name = os.path.join(output_dir, f"{base_name}_plotset_{i}.png")
        fig.savefig(save_name, facecolor='white')
        plt.close(fig)

if __name__ == "__main__":
    main()
