#!/usr/bin/env python3

#COUNTS - python Automated_Plotting.py count --data_dir C:\Users\lucca\Desktop\GeneratedOutput\SimOutput\test\CA_PA\repeat0
#SNAPSHOTS - python Automated_Plotting.py snapshots --parent_dir C:\Users\lucca\Desktop\GeneratedOutput\SimOutput\test\CA_PA

import os
import re
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse


#importing defined modules/functions
from RodShapedBacteria import RodShapedBacterium
import utilities as ut
import fastCellPlotting as fcp
import drag_functions as dfunc

ut.setMPL()

DEFAULT_OUTPUT_DIR = "C:\\Users\\lucca\\Desktop\\GeneratedOutput"

def plot_counts_over_time(data_dir, save_path=None, output_dir=DEFAULT_OUTPUT_DIR):
    file_pattern = os.path.join(data_dir, "biofilm_*.dat")
    files = sorted(glob.glob(file_pattern), key=lambda x: int(re.search(r'(\d+)', x).group(1)))

    time_steps = []
    zeta1_counts = []
    zetanot1_counts = []

    is_double= 'double' in data_dir  # check if path indicates 2 different zeta values

    plt.figure(figsize=(5, 3.5))
    for file_path in files:
        match = re.search(r'(\d+)', os.path.basename(file_path))
        if not match:
            continue
        time_step = int(match.group(1))
        df = pd.read_csv(file_path, sep="\t")
        zetas = df["non_dimzeta"].unique()

        time_steps.append(time_step * 0.1)

        for zeta in zetas:
            if zeta == 1.0:
                zeta1_count = df[df["non_dimzeta"] == 1.0].shape[0]
                zeta1_counts.append(zeta1_count)
            else:
                zetanot1_count = df[df["non_dimzeta"] != 1.0].shape[0]
                zetanot1_counts.append(zetanot1_count)

    for zeta in zetas:
        if zeta == 1.0 and len(zeta1_counts) == len(time_steps):
            plt.plot(time_steps, zeta1_counts, 'o', color='cyan', label="Zeta = 1")

        if zeta != 1.0 and len(zetanot1_counts) == len(time_steps):
            plt.plot(time_steps, zetanot1_counts, 'o', color='#9e003a', label="Zeta != 1")
            
            
    plt.xlabel("Time (h)")
    plt.ylabel("Cell/Segment Count")
    #plt.yscale("log")
    plt.legend()
    plt.tight_layout()


    if save_path is None:
        parts = os.path.normpath(data_dir).split(os.sep)
        tag = "_".join(parts[-2:])
        os.makedirs(output_dir, exist_ok=True)
        save_path = os.path.join(output_dir, f"counts_{tag}.pdf")

    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved cell count plot to: {save_path}")
    plt.show()

def plot_average_counts_over_repeats(parent_dir, output_dir=DEFAULT_OUTPUT_DIR):
    repeat_dirs = sorted(glob.glob(os.path.join(parent_dir, "repeat*")))
    all_zeta1 = {}
    all_zetanot1 = {}
    all_times = set()

    is_double = 'double' in parent_dir  # check if path indicates mutant mode


    for repeat_dir in repeat_dirs:
        file_pattern = os.path.join(repeat_dir, "biofilm_*.dat")
        files = sorted(glob.glob(file_pattern), key=lambda x: int(re.search(r'(\d+)', x).group(1)))
        for file_path in files:
            match = re.search(r'(\d+)', os.path.basename(file_path))
            if not match:
                continue
            time_step = int(match.group(1))
            time_h = time_step * 0.1
            df = pd.read_csv(file_path, sep="\t")

            zeta1 = df[df["non_dimzeta"] == 1.0].shape[0]
            all_zeta1.setdefault(time_h, []).append(zeta1)

            if is_double:
                cells = ut.getCells(file_path)
                zetanot1 = 0
                for cell in cells:
                    if cell.non_dimzeta != 1.0:
                        zetanot1 += 1
                all_zetanot1.setdefault(time_h, []).append(zetanot1)
            all_times.add(time_h)
        print(all_zeta1,"A")
        print("")
        print(all_zetanot1,"B")
        print("")

    sorted_times = sorted(all_times)
    avg_zeta1, sem_zeta1 = [], []
    avg_zetanot1, sem_zetanot1 = [], []

    for t in sorted_times:
        zeta1_vals = all_zeta1.get(t, [])
        while len(zeta1_vals) < len(repeat_dirs):
            zeta1_vals.append(0)
        avg_zeta1.append(np.mean(zeta1_vals))
        sem_zeta1.append(np.std(zeta1_vals, ddof=1) / np.sqrt(len(zeta1_vals)))

        if is_double:
            zetanot1_vals = all_zetanot1.get(t, [])
            while len(zetanot1_vals) < len(repeat_dirs):
                zetanot1_vals.append(0)
            avg_zetanot1.append(np.mean(zetanot1_vals))
            sem_zetanot1.append(np.std(zetanot1_vals, ddof=1) / np.sqrt(len(zetanot1_vals)))

    if is_double:
        plt.figure(figsize=(5.5, 4))
        plt.errorbar(sorted_times, avg_zeta1, yerr=sem_zeta1, fmt='o', color='cyan', label="Zeta = 1")
        plt.errorbar(sorted_times, avg_zetanot1, yerr=sem_zetanot1, fmt='o', color='#9e003a', label="Zeta != 1")
    else:
        plt.figure(figsize=(5.5, 4))
        plt.errorbar(sorted_times, avg_zeta1, yerr=sem_zeta1, fmt='o', color='cyan', label="Zeta = 1")

    plt.xlabel("Time (h)")
    plt.ylabel("Average Cell/Segment Count")
    plt.legend()
    plt.tight_layout()

    parts = os.path.normpath(parent_dir).split(os.sep)
    tag = "_".join(parts[-2:])
    os.makedirs(output_dir, exist_ok=True)
    save_path = os.path.join(output_dir, f"avg_counts_{tag}.pdf")
    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved average cell count plot to: {save_path}")
    plt.show()

def plot_cells_grid(parent_dir, output_path=None, num_snapshots=7, output_dir=DEFAULT_OUTPUT_DIR):
    repeat_dirs = sorted(glob.glob(os.path.join(parent_dir, "repeat*")))
    num_repeats = len(repeat_dirs)

    fig, axes = plt.subplots(num_repeats, num_snapshots, figsize=(15, 3 * num_repeats),
                             constrained_layout=True, facecolor='w')

    def plotCells(ax, file):
        dat = pd.read_csv(file, sep='\t')
        cells = ut.getCells(file)
        x_center,y_center = 0, 0

        maxx, minx = dat['pos_x'].max() + 5 - x_center, dat['pos_x'].min() - 5 - x_center
        maxy, miny = dat['pos_y'].max() + 5 - y_center, dat['pos_y'].min() - 5 - y_center
        X, Y = maxx - minx, maxy - miny
        X_c, Y_c = 0.5 * (maxx + minx), 0.5 * (maxy + miny)
        if X >= Y:
            Y = X
            miny, maxy = Y_c - 0.5 * Y, Y_c + 0.5 * Y
        else:
            X = Y
            minx, maxx = X_c - 0.5 * X, X_c + 0.5 * X

        print(cells[0].pos_x, cells[0].pos_y)
        cells = dfunc.centerCells(cells)
        print(cells[0].pos_x, cells[0].pos_y)

        fcp.addAllCellsToPlot(cells, ax, ax_rng=maxx - minx, show_id=False, ec='w')


        #Plots walls
        scale = 1
        wall_color = 'k'
        ax.plot([minx/scale, maxx/scale], [miny/scale, miny/scale], color=wall_color, alpha=0.6)
        ax.plot([minx/scale, maxx/scale], [maxy/scale, maxy/scale], color=wall_color, alpha=0.6)
        ax.plot([minx/scale, minx/scale], [miny/scale, maxy/scale], color=wall_color, alpha=0.6)
        ax.plot([maxx/scale, maxx/scale], [miny/scale, maxy/scale], color=wall_color, alpha=0.6)


        ax.set_xlim([-80, 80])
        ax.set_ylim([-5, 160])
        ax.axis('scaled')
        ax.axis('off')

    for r, repeat_dir in enumerate(repeat_dirs):
        file_pattern = os.path.join(repeat_dir, "biofilm_*.dat")
        files = sorted(glob.glob(file_pattern))
        if len(files) < num_snapshots:
            selected_files = files
        else:
            selected_indices = np.linspace(5, len(files) - 1, num_snapshots, dtype=int)
            selected_files = [files[i] for i in selected_indices]


            # # Let's take 4 snapshots spaced by len(files) // 4
            # selected_indices = [5]  # start from the beginning (or you can use 1 if you want to skip t=0)

            # # Add 3 more points evenly spaced
            # quarter = len(files) // 4
            # selected_indices += [quarter, 2 * quarter, 3 * quarter]

            # Get corresponding files
            selected_files = [files[i] for i in selected_indices]
        for c, (ax, file) in enumerate(zip(axes[r], selected_files)):
            plotCells(ax, file)
            match = re.search(r'biofilm_(\d+)\.dat$', os.path.basename(file))
            if match:
                frame_number = int(match.group(1))
                time_hours = frame_number * 0.1
                ax.set_title(f"{time_hours:.1f} h", color='k', fontsize=12)
        axes[r, 0].set_ylabel(f"Repeat {r+1}", fontsize=12, color='k')

    if output_path is None:
        parts = os.path.normpath(parent_dir).split(os.sep)
        tag = "_".join(parts[-2:])
        os.makedirs(output_dir, exist_ok=True)
        output_path_pdf = os.path.join(output_dir, f"snapshots_{tag}.pdf")
        output_path_png = os.path.join(output_dir, f"snapshots_{tag}.png")

    fig.savefig(output_path_pdf, format="pdf", dpi=600, bbox_inches="tight")
    fig.savefig(output_path_png, dpi=600, bbox_inches="tight")
    print(f"Saved snapshot grid to: {output_path_pdf}")
    plt.show()

def plot_Rg_over_time(data_dir, save_path=None, output_dir=DEFAULT_OUTPUT_DIR):
    file_pattern = os.path.join(data_dir, "biofilm_*.dat")
    files = sorted(glob.glob(file_pattern), key=lambda x: int(re.search(r'(\d+)', x).group(1)))

    time_steps = []
    Rg_tot = []
    Rg_zeta1 = []
    Rg_zetanot1 = []

    is_double= 'double' in data_dir  # check if path indicates 2 different zeta values


    plt.figure(figsize=(5, 3.5))
    for file_path in files:
        match = re.search(r'(\d+)', os.path.basename(file_path))
        if not match:
            continue

        time_step = int(match.group(1))
        cells = ut.getCells(file_path)
        zetas = set(cell.non_dimzeta for cell in cells)

        Rg_tot.append(dfunc.radiusGyration(cells))
        time_steps.append(time_step * 0.1)

        if len(zetas) != 1: #If only one zeta, only Rg_tot is required, so skip the rest
            for zeta in zetas:
                if zeta == 1.0:
                    Rg_zeta1.append(dfunc.RgZeta(cells, zeta))
                else:
                    Rg_zetanot1.append(dfunc.RgZeta(cells, zeta)) 
    
    if len(Rg_tot) == len(time_steps):
        plt.plot(time_steps, Rg_tot, 'o', color='black', label="Total Rg")

    if is_double:
        if len(Rg_zetanot1) == len(time_steps):
            plt.plot(time_steps, Rg_zetanot1, 'o', color='#9e003a', label="Zeta != 1")    

        if len(Rg_zeta1) == len(time_steps):
            plt.plot(time_steps, Rg_zeta1, 'o', color='cyan', label="Zeta = 1")
            
    plt.xlabel("Time (h)")
    plt.ylabel("Rg (microm)")
    plt.legend()
    plt.tight_layout()


    if save_path is None:
        parts = os.path.normpath(data_dir).split(os.sep)
        tag = "_".join(parts[-2:])
        os.makedirs(output_dir, exist_ok=True)
        save_path = os.path.join(output_dir, f"Rg_{tag}.pdf")

    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved Rg plot to: {save_path}")
    plt.show()


def main():
    parser = argparse.ArgumentParser(description="Automated plotting script for biofilm simulation data")
    subparsers = parser.add_subparsers(dest="command")

    count_parser = subparsers.add_parser("count", help="Plot cell counts over time")
    count_parser.add_argument("--data_dir", required=True)
    count_parser.add_argument("--save_path", required=False)
    count_parser.add_argument("--output_dir", required=False)

    snap_parser = subparsers.add_parser("snapshots", help="Plot snapshot grid")
    snap_parser.add_argument("--parent_dir", required=True)
    snap_parser.add_argument("--output_path", required=False)
    snap_parser.add_argument("--num_snapshots", type=int, default=5)
    snap_parser.add_argument("--output_dir", required=False)

    avg_parser = subparsers.add_parser("avg_counts", help="Plot averaged cell counts with SEM")
    avg_parser.add_argument("--parent_dir", required=True)
    avg_parser.add_argument("--output_dir", required=False)

    rg_parser = subparsers.add_parser("rg", help="Plot radius of gyration over time")
    rg_parser.add_argument("--data_dir", required=True)
    rg_parser.add_argument("--save_path", required=False)
    rg_parser.add_argument("--output_dir", required=False)

    all_parser = subparsers.add_parser("all", help="Plot all: counts, growth, snapshots, averages")
    all_parser.add_argument("--data_dir", required=True)
    all_parser.add_argument("--parent_dir", required=True)
    all_parser.add_argument("--num_snapshots", type=int, default=7)
    all_parser.add_argument("--output_dir", required=False)

    args = parser.parse_args()
    output_dir = DEFAULT_OUTPUT_DIR

    if args.command == "count":
        plot_counts_over_time(args.data_dir, args.save_path, output_dir)
    elif args.command == "snapshots":
        plot_cells_grid(args.parent_dir, args.output_path, args.num_snapshots, output_dir)
    elif args.command == "avg_counts":
        plot_average_counts_over_repeats(args.parent_dir, output_dir)
    elif args.command == "rg":
        plot_Rg_over_time(args.data_dir, args.save_path, output_dir)

    elif args.command == "all":
        plot_average_counts_over_repeats(args.parent_dir, output_dir)
        plot_cells_grid(args.parent_dir, None, args.num_snapshots, output_dir)  
    else:
        parser.print_help()

if __name__ == "__main__":
    main()

