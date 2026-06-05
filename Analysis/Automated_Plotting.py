#!/usr/bin/env python3

#COUNTS - python Automated_Plotting.py count --data_dir C:\Users\lucca\Desktop\GeneratedOutput\SimOutput\test\simple_FreeGrowth_double\repeat8
#SNAPSHOTS - python Automated_Plotting.py snapshots --parent_dir C:\Users\lucca\Desktop\GeneratedOutput\SimOutput\test\simple_FreeGrowth_double\ 
# AVERAGE COUNTS - python Automated_Plotting.py avg_counts --parent_dir C:\Users\lucca\Desktop\GeneratedOutput\SimOutput\test\simple_FreeGrowth_double\   
#rg - python Automated_Plotting.py rg --data_dir C:\Users\lucca\Desktop\GeneratedOutput\SimOutput\test\double_FreeGrowth\GR5\repeat1\

#photo (snapshot) - python Automated_Plotting.py single --file_path C:\Users\lucca\Desktop\GeneratedOutput\SimOutput\test\simple_FreeGrowth_double\repeat8\biofilm_100.dat

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
    time_steps, Lambda1_counts, Lambdanot1_counts, Lambdas = dfunc.counts(data_dir)
    
    plt.figure(figsize=(5, 3.5))
    for Lambda in Lambdas:
        if Lambda == 1.0 and len(Lambda1_counts) == len(time_steps):
            plt.plot(time_steps, Lambda1_counts, 'o', color=dfunc.colour_Lambda(Lambda), label="Lambda = 1")

        if Lambda != 1.0 and len(Lambdanot1_counts) == len(time_steps):
            plt.plot(time_steps, Lambdanot1_counts, 'o', color=dfunc.colour_Lambda(Lambda), label="Lambda = "+str(Lambda))

            
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

def plot_cells_grid(parent_dir, output_path=None, num_snapshots=5, output_dir=DEFAULT_OUTPUT_DIR):
    repeat_dirs = sorted(glob.glob(os.path.join(parent_dir, "repeat*")))
    num_repeats = len(repeat_dirs)

    fig, axes = plt.subplots(num_repeats, num_snapshots, figsize=(10, 2* num_repeats),
                             constrained_layout=True, facecolor='w')

    for r, repeat_dir in enumerate(repeat_dirs):
        file_pattern = os.path.join(repeat_dir, "biofilm_*.dat")
        files = sorted(glob.glob(file_pattern))
        if len(files) < num_snapshots:
            selected_files = files
        else:
            selected_indices = np.linspace(0, len(files) - 1, num_snapshots, dtype=int)
            selected_files = [files[i] for i in selected_indices]


            # # Let's take 4 snapshots spaced by len(files) // 4
            # selected_indices = [5]  # start from the beginning (or you can use 1 if you want to skip t=0)

            # # Add 3 more points evenly spaced
            # quarter = len(files) // 4
            # selected_indices += [quarter, 2 * quarter, 3 * quarter]

            # Get corresponding files
            selected_files = [files[i] for i in selected_indices]
        if len(axes.shape) == 1:  # If only one repeat, axes will be 1D
            for c, (ax, file) in enumerate(zip(axes, selected_files)):
                dfunc.plotCells(ax, file)
                match = re.search(r'biofilm_(\d+)\.dat$', os.path.basename(file))
                if match:
                    frame_number = int(match.group(1))
                    time_hours = frame_number * 0.1
                    ax.set_title(f"{time_hours:.1f} h", color='k', fontsize=12)
        else:
            for c, (ax, file) in enumerate(zip(axes[r], selected_files)):
                dfunc.plotCells(ax, file)
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
    time_steps, Rg_tot, Lambdas, Rg_Lambda1, Rg_Lambdanot1 = dfunc.RgLambda_time(data_dir)

    plt.figure(figsize=(5, 3.5))
    if len(Rg_tot) == len(time_steps):
        plt.plot(time_steps, Rg_tot, 'o', color='black', label="Total Rg")

    if len(Lambdas) != 1:
        if len(Rg_Lambdanot1) == len(time_steps):
            plt.plot(time_steps, Rg_Lambdanot1, 'o', color=dfunc.colour_Lambda(Lambdas[Lambdas != 1.0][0]), label="Lambda =" + str(Lambdas[Lambdas != 1.0][0]))    

        if len(Rg_Lambda1) == len(time_steps):
            plt.plot(time_steps, Rg_Lambda1, 'o', color=dfunc.colour_Lambda(1.0), label="Lambda = 1")
            
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

def plot_single_snapshot(file_path, output_dir=DEFAULT_OUTPUT_DIR):
    fig, axes = plt.subplots(1, 1, figsize=(15, 3),
                             constrained_layout=True, facecolor='w')
    
    dfunc.plotCells(axes, file_path)

    plt.tight_layout()
    

    parts = os.path.normpath(file_path).split(os.sep)
    tag = "_".join(parts[-2:])
    os.makedirs(output_dir, exist_ok=True)
    save_path = os.path.join(output_dir, f"avg_counts_{tag}.pdf")
    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved snapshot grid to: {save_path}")
    plt.show()

def plot_growth_rate(data_dir, save_path=None, output_dir=DEFAULT_OUTPUT_DIR):
    time_steps, Lambda1_counts, Lambdanot1_counts, Lambdas = dfunc.counts(data_dir)
    growth_rate_Lambda1,err_Lambda1 = dfunc.estimate_growth_rate(Lambda1_counts, time_step=0.1)
    growth_rate_Lambdanot1, err_Lambdanot1 = dfunc.estimate_growth_rate(Lambdanot1_counts, time_step=0.1)
    print(f"Estimated doubling time for Lambda=1: {growth_rate_Lambda1}+-{err_Lambda1}, Lambda={Lambdas[Lambdas != 1.0][0]}: {growth_rate_Lambdanot1}+-{err_Lambdanot1}")

    plt.figure(figsize=(5, 3.5))
    if len(Lambda1_counts) == len(time_steps):
        plt.plot(time_steps, Lambda1_counts, 'o', color=dfunc.colour_Lambda(1.0), label=" Counts Lambda = 1")
        plt.plot(time_steps, 2**(np.array(time_steps)/growth_rate_Lambda1), '-', color=dfunc.colour_Lambda(1.0), label="Fit Lambda1: $N(t)= 2^{(t/"+str(round(growth_rate_Lambda1, 2)) +")}$")

    if len(Lambdanot1_counts) == len(time_steps):
        plt.plot(time_steps, Lambdanot1_counts, 'o', color=dfunc.colour_Lambda(Lambdas[Lambdas != 1.0][0]), label=" Counts Lambda =" + str(Lambdas[Lambdas != 1.0][0]))
        plt.plot(time_steps, 2**(np.array(time_steps)/growth_rate_Lambdanot1), '-', color=dfunc.colour_Lambda(Lambdas[Lambdas != 1.0][0]), label="Fit Lambda=" + str(Lambdas[Lambdas != 1.0][0]) + ": $N(t)= 2^{(t/"+str(round(growth_rate_Lambdanot1, 2)) +")}$")

    plt.xlabel("Time (h)")
    plt.ylabel("Count N(t)")
    plt.legend()
    plt.tight_layout()

    if save_path is None:
        parts = os.path.normpath(data_dir).split(os.sep)
        tag = "_".join(parts[-2:])
        os.makedirs(output_dir, exist_ok=True)
        save_path = os.path.join(output_dir, f"growth_rate_{tag}.pdf")

    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved growth rate plot to: {save_path}")
    plt.show()


def plot_average_counts_over_repeats(parent_dir, output_dir=DEFAULT_OUTPUT_DIR):
    repeat_dirs = sorted(glob.glob(os.path.join(parent_dir, "repeat*")))
    all_Lambda1 = {}
    all_Lambdanot1 = {}
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

            Lambda1 = dfunc.find_Lambda_cells(df, 1.0).shape[0]
            all_Lambda1.setdefault(time_h, []).append(Lambda1)

            if is_double:
                cells = ut.getCells(file_path)
                Lambdanot1 = dfunc.find_Lambda_cells(cells, 1.0).shape[0]
                all_Lambdanot1.setdefault(time_h, []).append(Lambdanot1)
            all_times.add(time_h)

    sorted_times = sorted(all_times)
    avg_Lambda1, sem_Lambda1 = [], []
    avg_Lambdanot1, sem_Lambdanot1 = [], []

    for t in sorted_times:
        Lambda1_vals = all_Lambda1.get(t, [])
        while len(Lambda1_vals) < len(repeat_dirs):
            Lambda1_vals.append(0)
        avg_Lambda1.append(np.mean(Lambda1_vals))
        sem_Lambda1.append(np.std(Lambda1_vals, ddof=1) / np.sqrt(len(Lambda1_vals)))

        if is_double:
            Lambdanot1_vals = all_Lambdanot1.get(t, [])
            while len(Lambdanot1_vals) < len(repeat_dirs):
                Lambdanot1_vals.append(0)
            avg_Lambdanot1.append(np.mean(Lambdanot1_vals))
            sem_Lambdanot1.append(np.std(Lambdanot1_vals, ddof=1) / np.sqrt(len(Lambdanot1_vals)))

    if is_double:
        plt.figure(figsize=(5.5, 4))
        plt.errorbar(sorted_times, avg_Lambda1, yerr=sem_Lambda1, fmt='o', color=dfunc.colour_Lambda(1.0), label="Lambda = 1")
        plt.errorbar(sorted_times, avg_Lambdanot1, yerr=sem_Lambdanot1, fmt='o', color=dfunc.colour_Lambda(Lambdas[Lambdas != 1.0][0]), label="Lambda = " + str(Lambdas[Lambdas != 1.0][0]))
    else:
        plt.figure(figsize=(5.5, 4))
        plt.errorbar(sorted_times, avg_Lambda1, yerr=sem_Lambda1, fmt='o', color=dfunc.colour_Lambda(1.0), label="Lambda = 1")

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
#FIX THIS ONE!!!!


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

    single_parser = subparsers.add_parser("single", help="Plot single snapshot")
    single_parser.add_argument("--file_path", required=True)
    single_parser.add_argument("--output_dir", required=False)

    growth_parser = subparsers.add_parser("growth_rate", help="Plot growth rate over time")
    growth_parser.add_argument("--data_dir", required=True)
    growth_parser.add_argument("--save_path", required=False)
    growth_parser.add_argument("--output_dir", required=False)

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
    elif args.command == "single":
        plot_single_snapshot(args.file_path, output_dir)
    elif args.command == "growth_rate":
        plot_growth_rate(args.data_dir, args.save_path, output_dir)


    elif args.command == "all":
        plot_average_counts_over_repeats(args.parent_dir, output_dir)
        plot_cells_grid(args.parent_dir, None, args.num_snapshots, output_dir)  
    else:
        parser.print_help()

if __name__ == "__main__":
    main()

