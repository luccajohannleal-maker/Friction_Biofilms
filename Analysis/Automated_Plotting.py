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
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches


#importing defined modules/functions
import utilities as ut
import drag_functions as dfunc
import plotting_functions as pfunc

ut.setMPL()

DEFAULT_OUTPUT_DIR = "C:\\Users\\lucca\\Desktop\\GeneratedOutput"

def plot_counts_over_time(data_dirs, output_dir=DEFAULT_OUTPUT_DIR):
    plt.figure(figsize=(5, 3.5))

    if type(data_dirs) == str:
        data_dirs = [data_dirs]
    else:
        data_dirs = list(data_dirs)

    Lambdas = []

    for data_dir in data_dirs:
        Lambdas += list(pfunc.plot_count(data_dir))
    
    Lambdas = set(Lambdas)
    legend_elements = []
    for Lambda in Lambdas:
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label='$\Lambda =$'+str(Lambda)) ]
    plt.legend(handles=legend_elements)

    plt.xlabel("Time (h)")
    plt.ylabel("Cell/Segment Count")
    plt.yscale("log")
    plt.tight_layout()


    os.makedirs(output_dir, exist_ok=True)
    save_path = os.path.join(output_dir, f"counts.pdf")

    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved cell count plot to: {save_path}")
    plt.show()

def plot_cells_grid(data_dirs, num_snapshots=5, output_dir=DEFAULT_OUTPUT_DIR):
    #repeat_dirs = sorted(glob.glob(os.path.join(parent_dir, "repeat*")))
    #num_repeats = len(repeat_dirs)
    if type(data_dirs) == str:
        num_repeats = 1
        data_dirs = [data_dirs]
    
    else:
        data_dirs = list(data_dirs)
        num_repeats = len(data_dirs)

    fig, axes = plt.subplots(num_repeats, num_snapshots, figsize=(10, 2* num_repeats),
                             constrained_layout=True, facecolor='w')

    for r, data_dir in enumerate(data_dirs):
        file_pattern = os.path.join(data_dir, "biofilm_*.dat")
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


    os.makedirs(output_dir, exist_ok=True)
    output_path_pdf = os.path.join(output_dir, f"snapshots.pdf")

    fig.savefig(output_path_pdf, format="pdf", dpi=600, bbox_inches="tight")
    print(f"Saved snapshot grid to: {output_path_pdf}")
    plt.show()

def plot_channels(data_dirs, num_snapshots=5,width=120, output_dir=DEFAULT_OUTPUT_DIR):
    #repeat_dirs = sorted(glob.glob(os.path.join(parent_dir, "repeat*")))
    #num_repeats = len(repeat_dirs)
    if type(data_dirs) == str:
        num_repeats = 1
        data_dirs = [data_dirs]
    
    else:
        data_dirs = list(data_dirs)
        num_repeats = len(data_dirs)

    fig, axes = plt.subplots(num_repeats, num_snapshots, figsize=(10, 2* num_repeats),
                             constrained_layout=True, facecolor='w')

    for r, data_dir in enumerate(data_dirs):
        file_pattern = os.path.join(data_dir, "biofilm_*.dat")
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
                dfunc.plotCells_channel(ax, file, width)
                match = re.search(r'biofilm_(\d+)\.dat$', os.path.basename(file))
                if match:
                    frame_number = int(match.group(1))
                    time_hours = frame_number * 0.1
                    ax.set_title(f"{time_hours:.1f} h", color='k', fontsize=12)
        else:
            for c, (ax, file) in enumerate(zip(axes[r], selected_files)):
                dfunc.plotCells_channel(ax, file, width)
                match = re.search(r'biofilm_(\d+)\.dat$', os.path.basename(file))
                if match:
                    frame_number = int(match.group(1))
                    time_hours = frame_number * 0.1
                    ax.set_title(f"{time_hours:.1f} h", color='k', fontsize=12)
            axes[r, 0].set_ylabel(f"Repeat {r+1}", fontsize=12, color='k')


    os.makedirs(output_dir, exist_ok=True)
    output_path_pdf = os.path.join(output_dir, f"channel_snapshot.pdf")

    fig.savefig(output_path_pdf, format="pdf", dpi=600, bbox_inches="tight")
    print(f"Saved channel snapshot grid to: {output_path_pdf}")
    plt.show()


def plot_Rg_over_time(data_dirs, save_path=None, output_dir=DEFAULT_OUTPUT_DIR):
    plt.figure(figsize=(5, 3.5))

    if type(data_dirs) == str:
        data_dirs = [data_dirs]
    else:
        data_dirs = list(data_dirs)

    Lambdas = []

    for data_dir in data_dirs:
        Lambdas += list(pfunc.plot_rg_linear(data_dir))

    Lambdas = set(Lambdas)
    legend_elements = []
    for Lambda in Lambdas:
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label='$\Lambda =$'+str(Lambda)) ]
    plt.legend(handles=legend_elements)
            
    plt.xlabel("Time (h)")
    plt.ylabel("$log_2 [R_g] $ (microm)")
    plt.tight_layout()


    os.makedirs(output_dir, exist_ok=True)
    save_path = os.path.join(output_dir, f"Rg.pdf")

    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved Rg plot to: {save_path}")
    plt.show()

def plot_single_snapshot(file_path, output_dir=DEFAULT_OUTPUT_DIR):
    fig, axes = plt.subplots(1, 1, figsize=(15, 3),
                             constrained_layout=True, facecolor='w')
    
    Lambdas = list(dfunc.plotCells(axes, file_path))
    legend_elements = []
    for Lambda in Lambdas:
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label='$\Lambda =$'+str(Lambda)) ]
    axes.legend(handles=legend_elements)

    plt.tight_layout()
    

    parts = os.path.normpath(file_path).split(os.sep)
    tag = "_".join(parts[-2:])
    os.makedirs(output_dir, exist_ok=True)
    save_path = os.path.join(output_dir, f"avg_counts_{tag}.pdf")
    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved snapshot grid to: {save_path}")
    plt.show()

def plot_growth_rate(data_dirs, output_dir=DEFAULT_OUTPUT_DIR):
    plt.figure(figsize=(5, 3.5))
    if type(data_dirs) == str:
        data_dirs = [data_dirs]
    else:
        data_dirs = list(data_dirs)

    Lambdas = []

    for data_dir in data_dirs:
        Lambdas= Lambdas + list(pfunc.plot_GR_linear(data_dir)) #can replace plot_GR_exp with plot_GR_linear

    Lambdas = set(Lambdas)
    legend_elements = []
    for Lambda in Lambdas:
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label='$\Lambda =$'+str(Lambda)) ]
    plt.legend(handles=legend_elements)

    plt.xlabel("Time (h)")
    plt.tight_layout()

    save_path = os.path.join(output_dir, f"growth_rate_GR.pdf")

    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved growth rate plot to: {save_path}")
    plt.show()

def plot_shape_asphericity(data_dirs, output_dir=DEFAULT_OUTPUT_DIR):
    plt.figure(figsize=(5, 3.5))

    if type(data_dirs) == str:
        data_dirs = [data_dirs]
    else:
        data_dirs = list(data_dirs)

    Lambdas =[]

    for data_dir in data_dirs:
        Lambdas= Lambdas + list(pfunc.plot_shape_asphericity_time(data_dir))

    Lambdas = set(Lambdas)
    legend_elements = []
    for Lambda in Lambdas:
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label='$\Lambda =$'+str(Lambda)) ]
    plt.legend(handles=legend_elements)

    plt.xlabel("Time (h)")
    plt.ylabel("Shape asphericity $A$")
    plt.tight_layout()

    save_path = os.path.join(output_dir, f"asphericity_colony.pdf")

    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved asphericity plot to: {save_path}")
    plt.show()

def plot_dAsph(data_dirs, output_dir=DEFAULT_OUTPUT_DIR):
    plt.figure(figsize=(5, 3.5))

    if type(data_dirs) == str:
        data_dirs = [data_dirs]
    else:
        data_dirs = list(data_dirs)

    Lambdas =[]

    for data_dir in data_dirs:
        Lambdas= Lambdas + pfunc.plot_delta_asph_time(data_dir)

    Lambdas = set(Lambdas)
    legend_elements = []
    for Lambda in Lambdas:
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label='$\Lambda$ ='+str(Lambda)) ]
    plt.legend(handles=legend_elements)

    plt.xlabel("Time (h)")
    plt.ylabel("Difference in $A$")
    plt.tight_layout()

    save_path = os.path.join(output_dir, f"Dasphericity_colony.pdf")

    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved difference in asphericity plot to: {save_path}")
    plt.show()

def plot_stress_t(data_dirs, output_dir=DEFAULT_OUTPUT_DIR):
    plt.figure(figsize=(5, 3.5))

    if type(data_dirs) == str:
        data_dirs = [data_dirs]
    else:
        data_dirs = list(data_dirs)

    Lambdas = []
    for data_dir in data_dirs:
        Lambdas = Lambdas + list(pfunc.plot_stress_time(data_dir))
    
    Lambdas = set(Lambdas)
    legend_elements = [Line2D([0], [0], color='k', marker="v",label= r"$\tau_2$", markersize=15),
                    Line2D([0], [0], color='k', marker="*",label= r"$\tau_1$", markersize=15),
                    Line2D([0], [0], color='k',marker="o",label= "$\sigma_{\parallel}$", markersize=15),
                    Line2D([0], [0], color='k',marker="x",label= "$\sigma_{\perp}$", markersize=15),
                    ]
    for Lambda in Lambdas:
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label='$\Lambda =$'+str(Lambda)) ]

    plt.legend(handles=legend_elements)
    plt.xlabel("Time (h)")
    plt.ylabel(r"$\sigma_{\parallel}$, $\sigma_{\perp}$ and $\tau$")
    plt.tight_layout()

    save_path = os.path.join(output_dir, f"stress_colony_time.pdf")

    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved stress over time plot to: {save_path}")
    plt.show()

def plot_pressure_t(data_dirs, output_dir=DEFAULT_OUTPUT_DIR):
    plt.figure(figsize=(5, 3.5))

    if type(data_dirs) == str:
        data_dirs = [data_dirs]
    else:
        data_dirs = list(data_dirs)

    Lambdas = []
    for data_dir in data_dirs:
        Lambdas = Lambdas + list(pfunc.plot_pressure_time(data_dir))
    
    Lambdas = set(Lambdas)
    legend_elements = [
                    Line2D([0], [0], color='k',marker="o",label= r"$\alpha$", markersize=15),
                    Line2D([0], [0], color='k',marker="x",label= "$p$", markersize=15),
                    ]
    for Lambda in Lambdas:
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label='$\Lambda =$'+str(Lambda)) ]

    plt.legend(handles=legend_elements)
    plt.xlabel("Time (h)")
    plt.ylabel(r'<pressure> and <$\alpha$> (Pa m)')
    plt.tight_layout()

    save_path = os.path.join(output_dir, f"pressure_time.pdf")

    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved pressure plot to: {save_path}")
    plt.show()



def plot_stress_dist(data_dirs, output_dir=DEFAULT_OUTPUT_DIR):
    plt.figure(figsize=(5, 3.5))

    if type(data_dirs) == str:
        data_dirs = [data_dirs]
    else:
        data_dirs = list(data_dirs)

    Lambdas = []
    for data_dir in data_dirs:
        files = np.asarray(pfunc.get_file_paths(data_dir))
        filepath = files[-1]
        print(filepath)
        Lambdas = Lambdas + list(pfunc.plot_stress_distance(filepath))

    Lambdas = set(Lambdas)
    legend_elements = [Line2D([0], [0], color='k', marker="v",label= r"$\tau_2$", markersize=15),
                    Line2D([0], [0], color='k', marker="*",label= r"$\tau_1$", markersize=15),
                    Line2D([0], [0], color='k',marker="o",label= "$\sigma_{\parallel}$", markersize=15),
                    Line2D([0], [0], color='k',marker="x",label= "$\sigma_{\perp}$", markersize=15),
                    ]
    for Lambda in Lambdas:
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label='$\Lambda =$'+str(Lambda)) ]

    plt.legend(handles=legend_elements)
    plt.xlabel("Distance from centre (microns)")
    plt.ylabel(r"$\sigma_{\parallel}$, $\sigma_{\perp}$ and $\tau$")
    plt.tight_layout()

    save_path = os.path.join(output_dir, f"stress_colony_distance.pdf")

    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved stress over distance plot to: {save_path}")
    plt.show()

def plot_pressure_dist(data_dirs, output_dir=DEFAULT_OUTPUT_DIR):
    plt.figure(figsize=(5, 3.5))

    if type(data_dirs) == str:
        data_dirs = [data_dirs]
    else:
        data_dirs = list(data_dirs)

    Lambdas = []
    for data_dir in data_dirs:
        files = np.asarray(pfunc.get_file_paths(data_dir))
        filepath = files[-1]
        print(filepath)
        Lambdas = Lambdas + list(pfunc.plot_pressure_distance(filepath))

    Lambdas = set(Lambdas)
    legend_elements = [
                    Line2D([0], [0], color='k',marker="o",label= r"$\alpha$", markersize=15),
                    Line2D([0], [0], color='k',marker="x",label= "$p$", markersize=15),
                    ]
    for Lambda in Lambdas:
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label='$\Lambda =$'+str(Lambda)) ]

    plt.legend(handles=legend_elements)
    plt.xlabel("Distance from centre (microns)")
    plt.ylabel(r'pressure and $\alpha$ (Pa m)')
    plt.tight_layout()

    save_path = os.path.join(output_dir, f"pressure_distance.pdf")

    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved pressure vs distance plot to: {save_path}")
    plt.show()


def average_growth_time(dirs_averaged,output_dir=DEFAULT_OUTPUT_DIR):
    """
    This function should have as input a LIST containing LISTS of the
    directories to be average and plotted at the same time
    e.g. [[file1, file2], [file3, file4]]
    1 and 2 will be averaged and plloted and 3 and 4 will be averaged and
    plotted.
    """
    plt.figure(figsize=(5, 3.5))
    Lambdas = []

    for data_dirs in dirs_averaged:
        Lambdas = Lambdas + list(pfunc.plot_average_growth(data_dirs)) 

    Lambdas = set(Lambdas)
    legend_elements = []
    for Lambda in Lambdas:
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label='$\Lambda =$'+str(Lambda)) ]
    plt.legend(handles=legend_elements)

    plt.xlabel("Time (h)")
    plt.ylabel("$log_2 <N(t)>$")
    plt.tight_layout()

    save_path = os.path.join(output_dir, f"average_GR.pdf")

    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved average growth rate plot to: {save_path}")
    plt.show()

def average_asphericity_time(dirs_averaged,output_dir=DEFAULT_OUTPUT_DIR):
    """
    This function should have as input a LIST containing LISTS of the
    directories to be average and plotted at the same time
    e.g. [[file1, file2], [file3, file4]]
    1 and 2 will be averaged and plloted and 3 and 4 will be averaged and
    plotted.
    """
    plt.figure(figsize=(5, 3.5))
    Lambdas = []

    for data_dirs in dirs_averaged:
        Lambdas = Lambdas + list(pfunc.plot_average_asphericity(data_dirs)) 

    Lambdas = set(Lambdas)
    legend_elements = []
    for Lambda in Lambdas:
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label='$\Lambda =$'+str(Lambda)) ]
    plt.legend(handles=legend_elements)

    plt.xlabel("Time (h)")
    plt.ylabel("$<A>$")
    plt.tight_layout()

    save_path = os.path.join(output_dir, f"average_asphericity.pdf")

    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved average asphericity plot to: {save_path}")
    plt.show()

def average_dasph_time(dirs_averaged,output_dir=DEFAULT_OUTPUT_DIR):
    """
    This function should have as input a LIST containing LISTS of the
    directories to be average and plotted at the same time
    e.g. [[file1, file2], [file3, file4]]
    1 and 2 will be averaged and plloted and 3 and 4 will be averaged and
    plotted.
    """
    plt.figure(figsize=(5, 3.5))
    Lambdas = []

    for data_dirs in dirs_averaged:
        Lambdas = Lambdas + list(pfunc.plot_average_dasph(data_dirs)) 

    Lambdas = set(Lambdas)
    legend_elements = []
    for Lambda in Lambdas:
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label='$\Lambda =$'+str(Lambda)) ]
    plt.legend(handles=legend_elements)

    plt.xlabel("Time (h)")
    plt.ylabel("$<\Delta A>$")
    plt.tight_layout()

    save_path = os.path.join(output_dir, f"average_dasphericity.pdf")

    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved average delta asphericity plot to: {save_path}")
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



def plot_t_doub_Lambda_av(dirs_averaged, output_dir=DEFAULT_OUTPUT_DIR, lines=False):
    """
    This function should have as input a LIST containing LISTS of the
    directories to be average and plotted at the same time
    e.g. [[file1, file2], [file3, file4]]
    1 and 2 will be averaged and plloted and 3 and 4 will be averaged and
    plotted.
    """
    plt.figure(figsize=(5, 3.5))
    Lambdas = []

    for data_dirs in dirs_averaged:
        Lambdas = Lambdas + list(pfunc.t_doub_Lambda(data_dirs)) 

    Lambdas = set(Lambdas)
    legend_elements = []
    for Lambda in Lambdas:
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label='$\Lambda =$'+str(Lambda)) ]
    plt.legend(handles=legend_elements)


    #plots growth rate lines
    if lines:
        plt.plot([0.5,10.5],[0.75,0.75],"--", color="k")
        plt.plot([0.5,10.5],[1.25,1.25],":", color="k")


    plt.xlabel("$\Lambda$")
    plt.ylabel('$t_{doubling}$ (h)')
    plt.xlim([0.5,10.5])

    plt.tight_layout()

    save_path = os.path.join(output_dir, f"tdoub_lambda.pdf")

    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved tdoub vs lambda plot to: {save_path}")
    plt.show()






def main():
    parser = argparse.ArgumentParser(description="Automated plotting script for biofilm simulation data")
    subparsers = parser.add_subparsers(dest="command")

    count_parser = subparsers.add_parser("count", help="Plot cell counts over time")
    count_parser.add_argument("--data_dirs", required=True)
    count_parser.add_argument("--output_dir", required=False)

    snap_parser = subparsers.add_parser("snapshots", help="Plot snapshot grid")
    snap_parser.add_argument("--data_dirs", required=True)
    snap_parser.add_argument("--num_snapshots", type=int, default=5)
    snap_parser.add_argument("--output_dir", required=False)

    avg_parser = subparsers.add_parser("avg_counts", help="Plot averaged cell counts with SEM")
    avg_parser.add_argument("--parent_dir", required=True)
    avg_parser.add_argument("--output_dir", required=False)

    rg_parser = subparsers.add_parser("rg", help="Plot radius of gyration over time")
    rg_parser.add_argument("--data_dirs", required=True)
    rg_parser.add_argument("--output_dir", required=False)

    single_parser = subparsers.add_parser("single", help="Plot single snapshot")
    single_parser.add_argument("--file_path", required=True)
    single_parser.add_argument("--output_dir", required=False)

    growth_parser = subparsers.add_parser("growth_rate", help="Plot growth rate over time")
    growth_parser.add_argument("--data_dirs", required=True)
    growth_parser.add_argument("--output_dir", required=False)

    args = parser.parse_args()
    output_dir = DEFAULT_OUTPUT_DIR

    if args.command == "count":
        plot_counts_over_time(args.data_dirs, output_dir)
    elif args.command == "snapshots":
        plot_cells_grid(args.data_dirs, args.num_snapshots, output_dir)
    elif args.command == "avg_counts":
        plot_average_counts_over_repeats(args.parent_dir, output_dir)
    elif args.command == "rg":
        plot_Rg_over_time(args.data_dirs, output_dir)
    elif args.command == "single":
        plot_single_snapshot(args.file_path, output_dir)
    elif args.command == "growth_rate":
        plot_growth_rate(args.data_dirs, output_dir)

    else:
        parser.print_help()

if __name__ == "__main__":
    main()

