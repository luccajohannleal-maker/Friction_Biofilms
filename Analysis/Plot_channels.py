import os
import re
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import argparse
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches


#importing defined modules/functions
import utilities as ut
import drag_functions as dfunc
import plotting_functions as pfunc

ut.setMPL()

DEFAULT_OUTPUT_DIR = "C:\\Users\\lucca\\Desktop\\GeneratedOutput"




def plot_channels(data_dirs, num_snapshots=5,width=60, output_dir=DEFAULT_OUTPUT_DIR):
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

def plot_exit_y(data_dirs,width=60,plot_tot=False, output_dir=DEFAULT_OUTPUT_DIR):
    if type(data_dirs) == str:
        num_repeats = 1
        data_dirs = [data_dirs]
    
    else:
        data_dirs = list(data_dirs)
        num_repeats = len(data_dirs)

    fig, axes = plt.subplots(1,num_repeats, figsize=(num_repeats*2.5, 2.5),
                             constrained_layout=True, facecolor='w')

    if num_repeats == 1:
        pfunc.plot_leaving_cells(axes, data_dirs[0], plot_tot, width)

    else:
        for (ax, file) in zip(axes, data_dirs):
            pfunc.plot_leaving_cells(ax, file, plot_tot, width)

    plt.tight_layout()

    os.makedirs(output_dir, exist_ok=True)
    output_path_pdf = os.path.join(output_dir, f"channel_y_exit.pdf")

    fig.savefig(output_path_pdf, format="pdf", dpi=600, bbox_inches="tight")
    print(f"Saved channel y exit to: {output_path_pdf}")
    plt.show()

def plot_frac_time(data_dirs,width=60, output_dir=DEFAULT_OUTPUT_DIR):
    plt.figure(figsize=(5, 3.5))

    if type(data_dirs) == str:
        data_dirs = [data_dirs]
    else:
        data_dirs = list(data_dirs)

    ratio = []

    for data_dir in data_dirs:
        #ratio += [pfunc.plot_fraction_time_ratio(data_dir,width)]
        ratio += pfunc.plot_fraction_time(data_dir,width)
    
    ratio = set(ratio)
    legend_elements = []
    for r in ratio:
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(r),label=r'Lambda = '+str(r)) ]
    plt.legend(handles=legend_elements)

    plt.xlabel("Time (h)")
    plt.ylabel("Fraction of cells")
    plt.tight_layout()

    os.makedirs(output_dir, exist_ok=True)
    save_path = os.path.join(output_dir, f"counts.pdf")
    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved cell count plot to: {save_path}")
    plt.show()

def average_frac(dirs_averaged,width=60,output_dir=DEFAULT_OUTPUT_DIR):
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
        Lambdas = Lambdas + list(pfunc.plot_average_fraction(data_dirs,width)) 

    Lambdas = set(Lambdas)
    legend_elements = []
    for Lambda in Lambdas:
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label=r'$\Lambda =$'+str(Lambda)) ]
    plt.legend(handles=legend_elements)

    plt.xlabel("Time (h)")
    plt.ylabel("Fraction of cells")
    plt.tight_layout()

    save_path = os.path.join(output_dir, f"average_fraction.pdf")

    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved average fraction plot to: {save_path}")
    plt.show()

def plot_initial_final_fraction(data_dirs,width=40, output_dir=DEFAULT_OUTPUT_DIR):
    plt.figure(figsize=(5, 3.5))

    if type(data_dirs) == str:
        data_dirs = [data_dirs]
    else:
        data_dirs = list(data_dirs)

    Lambdas = []

    for data_dir in data_dirs:
        Lambdas += list(pfunc.plot_initial_final(data_dir,width))
    Lambdas = set(Lambdas)
    legend_elements = [Line2D([0], [0], color='k',marker=".",label= r"$\Lambda$ lower", markersize=15),
                    Line2D([0], [0], color='k',marker="x",label= r"$\Lambda$ higher", markersize=15)]
    for Lambda in Lambdas:
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label=r'$\Lambda = 1$ and '+str(Lambda)) ]
    plt.legend(handles=legend_elements)

    plt.ylabel("Initial fraction")
    plt.xlabel("Final fraction")
    plt.tight_layout()

    os.makedirs(output_dir, exist_ok=True)
    save_path = os.path.join(output_dir, f"initialxfinal_fractions.pdf")
    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved initial x final fractions plot to: {save_path}")
    plt.show()

def plot_orientation_time(data_dirs,width=40, output_dir=DEFAULT_OUTPUT_DIR):
    plt.figure(figsize=(5, 3.5))

    if type(data_dirs) == str:
        data_dirs = [data_dirs]
    else:
        data_dirs = list(data_dirs)

    Lambdas = []

    for data_dir in data_dirs:
        Lambdas += list(pfunc.orientation_time(data_dir,width))
    Lambdas = set(Lambdas)
    legend_elements = []
    for Lambda in Lambdas:
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label=r'$\Lambda = 1$ and '+str(Lambda)) ]
    plt.legend(handles=legend_elements)

    plt.ylabel("Average angle (degrees)")
    plt.xlabel("Time (h)")
    plt.tight_layout()

    os.makedirs(output_dir, exist_ok=True)
    save_path = os.path.join(output_dir, f"initialxfinal_fractions.pdf")
    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved initial x final fractions plot to: {save_path}")
    plt.show()

def average_orientation(dirs_averaged,width=40,output_dir=DEFAULT_OUTPUT_DIR):
    """
    This function should have as input a LIST containing LISTS of the
    directories to be average and plotted at the same time
    e.g. [[file1, file2], [file3, file4]]
    1 and 2 will be averaged and plloted and 3 and 4 will be averaged and
    plotted.
    """
    plt.figure(figsize=(5, 3.5))

    ratio = []
    for data_dirs in dirs_averaged:
        ratio += pfunc.plot_average_orientation(data_dirs,width=width)

    ratio = sorted(list(set(ratio)))
    legend_elements = []
    for r in ratio:
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_ratio(r),label=r'ratio = '+str(r)) ]
    plt.legend(handles=legend_elements)

    plt.xlabel("Time (h)")
    plt.ylabel(r"$<\theta>$ (degrees)")
    plt.tight_layout()

    save_path = os.path.join(output_dir, f"average_orientation.pdf")

    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved average fraction plot to: {save_path}")
    plt.show()

def whisker_diagram(dirs_averaged,width=40, output_dir=DEFAULT_OUTPUT_DIR):
    plt.figure(figsize=(5, 3.5))

    Lambdas = []

    for data_dirs in dirs_averaged:
        Lambdas = Lambdas + list(pfunc.plot_whisker_initial_final(data_dirs,width)) 
    Lambdas = set(Lambdas)

    plt.xlabel("Number of lower friction cells at simulation start (total = 30)")
    plt.ylabel("Final fraction")
    plt.tight_layout()

    os.makedirs(output_dir, exist_ok=True)
    save_path = os.path.join(output_dir, f"initialxfinal_fractions.pdf")
    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved initial x final fractions plot to: {save_path}")
    plt.show()


def plot_x_middle(data_dirs,width=40, output_dir=DEFAULT_OUTPUT_DIR):
    plt.figure(figsize=(5, 3.5))

    if type(data_dirs) == str:
        data_dirs = [data_dirs]
    else:
        data_dirs = list(data_dirs)

    ratio = []

    for data_dir in data_dirs:
        ratio += [pfunc.fraction_center(data_dir,width)]
    
    ratio = set(ratio)
    legend_elements = []
    for r in ratio:
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_ratio(r),label=r'ratio = '+str(r)) ]
    plt.legend(handles=legend_elements)

    plt.xlabel("Time (h)")
    plt.ylabel("Fraction of cells in centre of channel")
    plt.tight_layout()

    os.makedirs(output_dir, exist_ok=True)
    save_path = os.path.join(output_dir, f"counts.pdf")
    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved cell count plot to: {save_path}")
    plt.show()

def plot_exit_frac(data_dirs, output_dir=DEFAULT_OUTPUT_DIR):
    plt.figure(figsize=(5, 3.5))

    if type(data_dirs) == str:
        data_dirs = [data_dirs]
    else:
        data_dirs = list(data_dirs)

    ratio = []

    for data_dir in data_dirs:
        r = pfunc.plot_leaving_cells(data_dir)
        if r == 0:
            continue
        ratio += r
    
    ratio = set(ratio)
    legend_elements = []
    for r in ratio:
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(r),label=r'ratio = '+str(r)) ]
    plt.legend(handles=legend_elements)

    plt.xlabel("Time (h)")
    plt.ylabel(r"$df_2^{left}/dt$")
    plt.tight_layout()

    os.makedirs(output_dir, exist_ok=True)
    save_path = os.path.join(output_dir, f"counts.pdf")
    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved cell count plot to: {save_path}")
    plt.show()

def plot_frac_IC_normalised(data_dirs,width=60, output_dir=DEFAULT_OUTPUT_DIR,norm=False):
    plt.figure(figsize=(5, 3.5))

    if type(data_dirs) == str:
        data_dirs = [data_dirs]
    else:
        data_dirs = list(data_dirs)

    ratio = []
    
    t0=0
    tf=0
    for data_dir in data_dirs:
        r,t = pfunc.cells_change_fraction(data_dir,width,norm)
        if t[0]<t0:
            t0 = t[0]
        if t[-1]>tf:
            tf = t[-1]
        ratio += [r]

    t0 -=t0*0.1
    tf -=tf*0.1
    
    plt.axline((0, 0), (0,1), linewidth=0.5, color='k')
    plt.axline((t0, 0), (tf,0), linewidth=0.5, color='k')
    
    ratio = set(ratio)
    legend_elements = [Line2D([0], [0], color='k', marker="v",label= r"$f_0 = 0.1$", markersize=10),
                    Line2D([0], [0], color='k', marker="*",label= r"$f_0 = 0.25$", markersize=10),
                    Line2D([0], [0], color='k',marker="o",label= r"$f_0 = 0.33$", markersize=10),
                    Line2D([0], [0], color='k',marker="x",label= r"$f_0 = 0.5$", markersize=10)]
    for r in ratio:
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(r),label=r'ratio = '+str(r)) ]
    plt.legend(handles=legend_elements)

    plt.xlabel(r"$(t-t_0)/\tau^*$")
    plt.ylabel(r"$(f_2(t) - f_0)/(1-f_0)$")
    plt.tight_layout()

    os.makedirs(output_dir, exist_ok=True)
    save_path = os.path.join(output_dir, f"counts.pdf")
    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved cell count plot to: {save_path}")
    #plt.show()

def fraction_IC_average(dirs_averaged,width=40, output_dir=DEFAULT_OUTPUT_DIR):
    plt.figure(figsize=(5, 3.5))

    Lambdas = []

    for data_dirs in dirs_averaged:
        Lambdas = Lambdas + list(pfunc.plot_whisker_initial_final(data_dirs,width)) 
    Lambdas = set(Lambdas)

    plt.xlabel("Number of lower friction cells at simulation start (total = 30)")
    plt.ylabel("Final fraction")
    plt.tight_layout()

    os.makedirs(output_dir, exist_ok=True)
    save_path = os.path.join(output_dir, f"initialxfinal_fractions.pdf")
    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved initial x final fractions plot to: {save_path}")
    plt.show()