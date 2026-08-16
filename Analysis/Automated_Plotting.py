#!/usr/bin/env python3

#COUNTS - python Automated_Plotting.py count --data_dir C:\Users\lucca\Desktop\GeneratedOutput\SimOutput\test\simple_FreeGrowth_double\repeat8
#SNAPSHOTS - python Automated_Plotting.py snapshots --parent_dir C:\Users\lucca\Desktop\GeneratedOutput\SimOutput\test\simple_FreeGrowth_double\ 
# AVERAGE COUNTS - python Automated_Plotting.py avg_counts --parent_dir C:\Users\lucca\Desktop\GeneratedOutput\SimOutput\test\simple_FreeGrowth_double\   
#rg - python Automated_Plotting.py rg --data_dir C:\Users\lucca\Desktop\GeneratedOutput\SimOutput\test\double_FreeGrowth\GR5\repeat1\

#photo (snapshot) - python Automated_Plotting.py single --file_path C:\Users\lucca\Desktop\GeneratedOutput\SimOutput\test\simple_FreeGrowth_double\repeat8\biofilm_100.dat

from logging import raiseExceptions
import os
import re
import glob
from turtle import color
from wsgiref.validate import ErrorWrapper
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
from DistributionFunctions import computeColonyContour
from shapely.geometry import Polygon
from generalPlotting import addNematicDirector
import fastCellPlotting as fcp
from generalPlotting import addNematicDirector
from matplotlib.animation import FuncAnimation, PillowWriter, FFMpegWriter
import matplotlib.animation as animation


#importing defined modules/functions
import utilities as ut
import drag_functions as dfunc
import plotting_functions as pfunc

ut.setMPL()

DEFAULT_OUTPUT_DIR = "C:\\Users\\lucca\\Desktop\\GeneratedOutput"

def plot_cells_grid(data_dirs, num_snapshots=5,director=False, defects=False, output_dir=DEFAULT_OUTPUT_DIR,label=False):
    if type(data_dirs) == str:
        num_repeats = 1
        data_dirs = [data_dirs]
    
    else:
        data_dirs = list(data_dirs)
        num_repeats = len(data_dirs)

    fig, axes = plt.subplots(num_repeats, num_snapshots, figsize=(num_snapshots*2, 2* num_repeats),
                             constrained_layout=True, facecolor='w')
    
    initial_time= 20
    Lambdas = []

    for r, data_dir in enumerate(data_dirs):
        file_pattern = os.path.join(data_dir, "biofilm_*.dat")
        files = sorted(glob.glob(file_pattern))
        if len(files) < num_snapshots:
            selected_files = files
        else:
            selected_indices = np.linspace(initial_time, len(files) - 1, num_snapshots, dtype=int)
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
                Lambdas = dfunc.plotCells(ax, file, director,defects)
                Lambdas = sorted(set(Lambdas))
                legend_elements = []
                for Lambda in Lambdas:
                    legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label=str(Lambda)) ]
                fig.legend(handles=legend_elements,title=r'$\Lambda_i$')
                match = re.search(r'biofilm_(\d+)\.dat$', os.path.basename(file))
                if match:
                    frame_number = int(match.group(1))
                    time_hours = frame_number * 0.1
                    fig.set_title(f"{time_hours:.1f} h", color='k', fontsize=12)
        else:
            for c, (ax, file) in enumerate(zip(axes[r], selected_files)):
                Lambdas += list(dfunc.plotCells(ax, file, director,defects))
                match = re.search(r'biofilm_(\d+)\.dat$', os.path.basename(file))
                if match:
                    frame_number = int(match.group(1))
                    time_hours = frame_number * 0.1
                    ax.set_title(f"{time_hours:.1f} h", color='k', fontsize=12)
            axes[r, 0].set_ylabel(f"Repeat {r+1}", fontsize=12, color='k')

    Lambdas = sorted(set(Lambdas))
    legend_elements = []
    for Lambda in Lambdas:
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label=str(round(Lambda,1))) ]
    fig.legend(handles=legend_elements,title=r'$\Lambda_i$',loc='outside right center',fontsize = 20)
    


    os.makedirs(output_dir, exist_ok=True)
    output_path_pdf = os.path.join(output_dir, f"snapshots.pdf")

    fig.savefig(output_path_pdf, format="pdf", dpi=600, bbox_inches="tight")
    print(f"Saved snapshot grid to: {output_path_pdf}")
    plt.style.use("dark_background")
    plt.show()

def plot_cells_com(data_dirs, num_snapshots=5, output_dir=DEFAULT_OUTPUT_DIR):
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
                dfunc.plotCellsCOM(ax, file)
                match = re.search(r'biofilm_(\d+)\.dat$', os.path.basename(file))
                if match:
                    frame_number = int(match.group(1))
                    time_hours = frame_number * 0.1
                    ax.set_title(f"{time_hours:.1f} h", color='k', fontsize=12)
        else:
            for c, (ax, file) in enumerate(zip(axes[r], selected_files)):
                dfunc.plotCellsCOM(ax, file)
                match = re.search(r'biofilm_(\d+)\.dat$', os.path.basename(file))
                if match:
                    frame_number = int(match.group(1))
                    time_hours = frame_number * 0.1
                    ax.set_title(f"{time_hours:.1f} h", color='k', fontsize=12)
            axes[r, 0].set_ylabel(f"Repeat {r+1}", fontsize=12, color='k')


    os.makedirs(output_dir, exist_ok=True)
    output_path_pdf = os.path.join(output_dir, f"snapshots.pdf")

    fig.savefig(output_path_pdf, format="pdf", dpi=600, bbox_inches="tight")
    print(f"Saved snapshot grid + com to: {output_path_pdf}")
    plt.show()

def make_animation(data_dir,filename="simulation", output_dir=DEFAULT_OUTPUT_DIR,width=0):
    fig, ax = plt.subplots(figsize=(8, 8))
    print(FFMpegWriter.isAvailable())
    print(animation.writers.list())

    def draw_frame(file):
        ax.clear()
        director = False
        defects = False

        streamplot = True and director

        dat = pd.read_csv(file, sep='\t')
        cells = ut.getCells(file)
        cells = dfunc.centerCells(cells)

        x_center, y_center = 0, 0

        maxx = dat['pos_x'].max() + 5 - x_center
        minx = dat['pos_x'].min() - 5 - x_center
        maxy = dat['pos_y'].max() + 5 - y_center
        miny = dat['pos_y'].min() - 5 - y_center

        X, Y = maxx-minx, maxy-miny
        X_c, Y_c = 0.5*(maxx+minx), 0.5*(maxy+miny)

        if X >= Y:
            Y = X
            miny, maxy = Y_c-0.5*Y, Y_c+0.5*Y
        else:
            X = Y
            minx, maxx = X_c-0.5*X, X_c+0.5*X

        fcp.addAllCellsToPlot(cells, ax, ax_rng=maxx-minx,
                            show_id=False, ec='w')

        if director:
            q_name = f"{file[:-4].replace('/','_')}_Q.npy"
            addNematicDirector(ax, cells, q_name,
                            streamplot=streamplot, dr=5)

        if defects:
            x_pos = np.asarray(dat["pos_x"])
            y_pos = np.asarray(dat["pos_y"])
            x_or = np.asarray(dat["ori_x"])
            y_or = np.asarray(dat["ori_y"])

            list_defects, s, phi = dfunc.locate_nematic_defects_vector(
                x_pos, y_pos, x_or, y_or)

            for defect in list_defects:
                pos = defect["pos"]
                color = "r" if float(defect["charge"]) == 0.5 else "b"
                ax.scatter(pos[0], pos[1], c=color)

        # walls
        wall_color = 'k'
        ax.set_xlim([minx,maxx])
        ax.set_ylim([miny,maxy])

        if width !=0:
            y_top= width/2
            y_bottom = -width/2
            scale = 1
            wall_color = 'k'
            ax.plot([y_bottom*1.5/scale, y_top*1.5/scale], [y_top/scale, y_top/scale], color=wall_color, alpha=0.6)
            ax.plot([y_bottom*1.5/scale, y_top*1.5/scale], [y_bottom/scale, y_bottom/scale], color=wall_color, alpha=0.6)
            ax.plot([y_bottom*1.5/scale, y_bottom*1.5/scale], [y_top/scale, y_bottom/scale],"--", color=wall_color, alpha=0.6)
            ax.plot([y_top*1.5/scale, y_top*1.5/scale], [y_top/scale, y_bottom/scale],"--", color=wall_color, alpha=0.6)

            ax.set_xlim([y_bottom*1.75/scale, y_top*1.75/scale])
            ax.set_ylim([y_bottom*1.2/scale,y_top*1.2/scale])
        ax.set_aspect('equal')
        ax.axis('off')

    files = pfunc.get_file_paths(data_dir)

    anim = FuncAnimation(
    fig,
    draw_frame,
    frames=files,
    interval=50,      # milliseconds between frames
    repeat=False
)
    writer = PillowWriter(fps=10)
    save_path = os.path.join(output_dir, str(filename) + f"GIF.gif")
    anim.save(save_path, writer=writer, dpi=200)

    writer = FFMpegWriter(fps=10)
    save_path = os.path.join(output_dir, str(filename) + f"MP4.mp4")
    anim.save(filename=save_path, writer=writer)


def plot_single_snapshot(file_path, output_dir=DEFAULT_OUTPUT_DIR,director=False,defects=False):
    fig, axes = plt.subplots(1, 1, figsize=(15, 15),
                             constrained_layout=True, facecolor='w')
    
    Lambdas = list(dfunc.plotCells(axes, file_path,director=director,defects=defects))
    c=False
    if c:
        cells = ut.getCells(file_path)
        cellshigh = dfunc.find_Lambda_cells(cells,Lambda=Lambdas[1])
        contour_total = pd.DataFrame(computeColonyContour(cells),columns=["x","y"])
        contour = np.asarray(computeColonyContour(cellshigh))

        Lambda_point =[]
        sensibility = 0.5

        for point in contour:
            d = np.sqrt((contour_total["x"]-point[0])**2 + (contour_total["y"]-point[1])**2)

            if d.min() < sensibility: #at least one point of the total contour is close to the high contour
                Lambda_point.append(Lambdas[1])
            else:
                Lambda_point.append(Lambdas[0])

        cx = contour[:,0]
        cy = contour[:,1]
        for Lambda in Lambdas:
            mask = np.where(np.asarray(Lambda_point)==Lambda)
            if Lambda == 1.0:
                axes.scatter(cx[mask],cy[mask],color="b",s=30)
            else:
                axes.scatter(cx[mask],cy[mask],color="g",s=30)

    """ legend_elements = []
    for Lambda in Lambdas:
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label=r'$\Lambda =$'+str(Lambda)) ]
    axes.legend(handles=legend_elements)"""

    plt.tight_layout()
    

    parts = os.path.normpath(file_path).split(os.sep)
    tag = "_".join(parts[-2:])
    os.makedirs(output_dir, exist_ok=True)
    save_path = os.path.join(output_dir, f"avg_counts_{tag}.pdf")
    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved snapshot grid to: {save_path}")
    plt.show()

def snapshot_director(data_dirs, output_dir=DEFAULT_OUTPUT_DIR):
    n_dirs = len(data_dirs)
    fig, axes = plt.subplots(n_dirs, 1, figsize=(3,3*n_dirs),
                             constrained_layout=True, facecolor='w')
    i=0
    for data_dir in data_dirs:
        files = np.asarray(pfunc.get_file_paths(data_dir))
        filepath = files[-1]
        print(filepath)
        Lambdas = list(dfunc.plotCells(axes[i], filepath, director=True))
        i +=1

    legend_elements = []
    for Lambda in Lambdas:
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label=r'$\Lambda =$'+str(Lambda)) ]
    #axes.legend(handles=legend_elements)

    plt.tight_layout()
    

    parts = os.path.normpath(data_dirs[0]).split(os.sep)
    tag = "_".join(parts[-2:])
    os.makedirs(output_dir, exist_ok=True)
    save_path = os.path.join(output_dir, f"director_snapshots_{tag}.pdf")
    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved snapshot with directors to: {save_path}")
    plt.show()

def snapshot_defects(data_dirs, output_dir=DEFAULT_OUTPUT_DIR):
    n_dirs = len(data_dirs)
    fig, axes = plt.subplots(n_dirs, 1, figsize=(3,3*n_dirs),
                             constrained_layout=True, facecolor='w')
    i=0
    for data_dir in data_dirs:
        files = np.asarray(pfunc.get_file_paths(data_dir))
        filepath = files[-1]
        print(filepath)
        Lambdas = list(dfunc.plotCells(axes[i], filepath, defects=True))
        i +=1

    legend_elements = []
    for Lambda in Lambdas:
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label=r'$\Lambda =$'+str(Lambda)) ]
    #axes.legend(handles=legend_elements)

    plt.tight_layout()
    

    parts = os.path.normpath(data_dirs[0]).split(os.sep)
    tag = "_".join(parts[-2:])
    os.makedirs(output_dir, exist_ok=True)
    save_path = os.path.join(output_dir, f"director_snapshots_{tag}.pdf")
    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved snapshot with directors to: {save_path}")
    plt.show()

def plot_aspect_ratio_time(data_dirs, output_dir=DEFAULT_OUTPUT_DIR):
    plt.figure(figsize=(5, 3.5))

    if type(data_dirs) == str:
        data_dirs = [data_dirs]
    else:
        data_dirs = list(data_dirs)

    for data_dir in data_dirs:
            time_steps, aspect_ratio = dfunc.calc_aspect_ratio(data_dir)
            plt.plot(time_steps, aspect_ratio, color="k")

    plt.xlabel("Time (h)")
    plt.ylabel("Aspect ratio")
    plt.tight_layout()


    os.makedirs(output_dir, exist_ok=True)
    save_path = os.path.join(output_dir, f"asp_ratio.pdf")

    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved aspect ratio plot to: {save_path}")
    plt.show()

def plot_counts_over_time(data_dirs, output_dir=DEFAULT_OUTPUT_DIR,channels=False,width=60,total=False):
    plt.figure()
    fig, axes = plt.subplots(1,1, figsize=(5, 3.5),)

    if type(data_dirs) == str:
        data_dirs = [data_dirs]
    else:
        data_dirs = list(data_dirs)

    if total:
        Lambdas = []
        for data_dir in data_dirs:
            Lambdas += list(pfunc.plot_count_tot(data_dir,channels,width))
        plt.ylabel(r"Total Cell Count, $N_\textrm{total}(t)$")
            
    else:
        Lambdas = []

        for data_dir in data_dirs:
            Lambdas += list(pfunc.plot_count(data_dir,channels,width))
        
        
        plt.ylabel(r"Cell Count, $N(t)$")
    plt.yscale("log")
    Lambdas = sorted(set(Lambdas))
    
    legend_elements = []
    for Lambda in Lambdas:
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label=str(Lambda)) ]
    first = axes.legend(handles=legend_elements,title=r'$\Lambda_i \neq 1$')
    axes.add_artist(first)

    second, = axes.plot([0,12],[480,480],"k--", label=r"Expected max $N_\textrm{total}$")
    axes.legend(handles=[second])
    plt.xlabel("Time (h)")
    plt.tight_layout()


    os.makedirs(output_dir, exist_ok=True)
    save_path = os.path.join(output_dir, f"counts.pdf")

    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved cell count plot to: {save_path}")
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

    Lambdas = sorted(set(Lambdas))
    legend_elements = []
    for Lambda in Lambdas:
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label=str(Lambda)) ]
    plt.legend(handles=legend_elements,title=r'$\Lambda_i$')
            
    plt.xlabel("Time (h)")
    plt.ylabel("$log_2 [R_g]$ (microm)")
    plt.tight_layout()


    os.makedirs(output_dir, exist_ok=True)
    save_path = os.path.join(output_dir, f"Rg.pdf")

    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved Rg plot to: {save_path}")
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
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label=r'$\Lambda =$'+str(Lambda)) ]
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

    Lambdas = sorted(set(Lambdas))
    legend_elements = []
    for Lambda in Lambdas:
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label=str(Lambda)) ]
    plt.legend(handles=legend_elements,title=r'$\Lambda_i$')

    plt.xlabel("Time (h)")
    plt.ylabel(r"Shape asphericity, $\Delta$")
    plt.tight_layout()

    save_path = os.path.join(output_dir, f"asphericity_colony.pdf")

    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved asphericity plot to: {save_path}")
    plt.show()

def plot_avl(data_dirs, output_dir=DEFAULT_OUTPUT_DIR):
    plt.figure(figsize=(5, 3.5))

    if type(data_dirs) == str:
        data_dirs = [data_dirs]
    else:
        data_dirs = list(data_dirs)

    Lambdas =[]

    for data_dir in data_dirs:
        Lambdas += pfunc.plot_av_length_time(data_dir)

    Lambdas = sorted(set(Lambdas))
    legend_elements = []
    for Lambda in Lambdas:
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label=str(Lambda)) ]
    plt.legend(handles=legend_elements,title=r'$\Lambda_i$')

    plt.xlabel("Time (h)")
    plt.ylabel(r"Average length, $\bar \ell$")
    plt.tight_layout()

    save_path = os.path.join(output_dir, f"average_length.pdf")

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
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label=r'$\Lambda$ ='+str(Lambda)) ]
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
                    Line2D([0], [0], color='k',marker="o",label= r"$\sigma_{\parallel}$", markersize=15),
                    Line2D([0], [0], color='k',marker="x",label= r"$\sigma_{\perp}$", markersize=15)
                    ]
    for Lambda in Lambdas:
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label=r'$\Lambda =$'+str(Lambda)) ]

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
                    Line2D([0], [0], color='k',marker="x",label= r"$p$", markersize=15),
                    ]
    for Lambda in Lambdas:
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label=r'$\Lambda =$'+str(Lambda)) ]

    plt.legend(handles=legend_elements)
    plt.xlabel("Time (h)")
    plt.ylabel(r'<pressure> and <$\alpha$> (Pa m)')
    plt.tight_layout()

    save_path = os.path.join(output_dir, f"pressure_time.pdf")

    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved pressure plot to: {save_path}")
    plt.show()

def plot_COM_time(data_dirs, output_dir=DEFAULT_OUTPUT_DIR,log=False):
    plt.figure(figsize=(5, 3.5))

    if type(data_dirs) == str:
        data_dirs = [data_dirs]
    else:
        data_dirs = list(data_dirs)

    Lambdas = []

    for data_dir in data_dirs:
        Lambdas += list(pfunc.plot_COM_interacting(data_dir))
    
    Lambdas = sorted(set(Lambdas))
    legend_elements = []
    if len(Lambdas)==1:
        legend_elements = [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambdas[0]),label=r'$\Lambda =$'+str(Lambdas[0])) ]
    else:
        #legend_elements.append(Line2D([0], [0], color='k',linestyle="--",label='Total'))
        for Lambda in Lambdas:
            legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label=Lambda) ]
    plt.legend(handles=legend_elements,title=r'$\Lambda_i$')

    plt.xlabel(r'Time after collision, $t - t_{collision}$ (h)')
    plt.ylabel(r"$\vec r_{COM} - \vec r_{0}$")
    if log:
        plt.yscale("log")
        plt.xscale("log")
    plt.tight_layout()


    os.makedirs(output_dir, exist_ok=True)
    save_path = os.path.join(output_dir, f"com_movement.pdf")

    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved COM movement plot to: {save_path}")
    plt.show()

def plot_IQ(data_dirs, output_dir=DEFAULT_OUTPUT_DIR):
    plt.figure(figsize=(5, 3.5))

    if type(data_dirs) == str:
        data_dirs = [data_dirs]
    else:
        data_dirs = list(data_dirs)

    Lambdas = []

    for data_dir in data_dirs:
        Lambdas += list(pfunc.plot_IQ_time(data_dir))
    
    Lambdas = set(Lambdas)
    legend_elements = []
    if len(Lambdas)==1:
        legend_elements = [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambdas[0]),label=r'$\Lambda =$'+str(Lambdas[0])) ]
    else:
        legend_elements.append(Line2D([0], [0], color='k',linestyle="--",label='Total'))
        for Lambda in Lambdas:
            legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label=r'$\Lambda =$'+str(Lambda)) ]
    plt.legend(handles=legend_elements)

    plt.xlabel("Time (h)")
    plt.ylabel(r"$IQ$")
    plt.tight_layout()


    os.makedirs(output_dir, exist_ok=True)
    save_path = os.path.join(output_dir, f"IQ.pdf")

    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved IQ plot to: {save_path}")
    plt.show()

def plot_surface_fraction_time(data_dirs, output_dir=DEFAULT_OUTPUT_DIR):
    plt.figure(figsize=(5, 3.5))

    if type(data_dirs) == str:
        data_dirs = [data_dirs]
    else:
        data_dirs = list(data_dirs)

    Lambdas = []

    for data_dir in data_dirs:
        print(data_dir)
        Lambdas += list(pfunc.surface_fraction_time(data_dir))
    
    Lambdas = set(Lambdas)
    legend_elements = []
    for Lambda in Lambdas:
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label=r'$\Lambda =$'+str(Lambda)) ]
    plt.legend(handles=legend_elements)

    plt.xlabel("Time after collision (h)")
    plt.ylabel("fraction of cells occupying surface")
    plt.tight_layout()

    os.makedirs(output_dir, exist_ok=True)
    save_path = os.path.join(output_dir, f"com_movement.pdf")
    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved surface fraction plot to: {save_path}")
    plt.show()

def plot_interface_fraction_time(data_dirs, output_dir=DEFAULT_OUTPUT_DIR):
    plt.figure(figsize=(5, 3.5))

    if type(data_dirs) == str:
        data_dirs = [data_dirs]
    else:
        data_dirs = list(data_dirs)

    Lambdas = []

    for data_dir in data_dirs:
        print(data_dir)
        Lambdas += list(pfunc.interfacexsurface_time(data_dir))
    
    Lambdas = set(Lambdas)
    legend_elements = []
    for Lambda in Lambdas:
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label=r'$\Lambda =$'+str(Lambda)) ]
    plt.legend(handles=legend_elements)

    plt.xlabel("Time after collision (h)")
    plt.ylabel("fraction of the contour occupying interface")
    plt.tight_layout()

    os.makedirs(output_dir, exist_ok=True)
    save_path = os.path.join(output_dir, f"inner_contour_fraction.pdf")
    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved surface fraction plot to: {save_path}")
    plt.show()

def average_interface_fraction_ratio(dirs_averaged,output_dir=DEFAULT_OUTPUT_DIR):
    """
    This function should have as input a LIST containing LISTS of the
    directories to be average and plotted at the same time
    e.g. [[file1, file2], [file3, file4]]
    1 and 2 will be averaged and plloted and 3 and 4 will be averaged and
    plotted.
    """
    plt.figure(figsize=(5, 3.5))
    for data_dirs in dirs_averaged:
        pfunc.plot_average_interfacexsurface_ratio(data_dirs)

    legend_elements = [mpatches.Patch(facecolor=dfunc.colour_Lambda(1.0),label=r'$\Lambda_1/\Lambda_2 =$'+str(1.0)),
                       mpatches.Patch(facecolor=dfunc.colour_Lambda(1.5),label=r'$\Lambda_1/\Lambda_2 =$'+str(1.5)),
                       mpatches.Patch(facecolor=dfunc.colour_Lambda(2.0),label=r'$\Lambda_1/\Lambda_2 =$'+str(2.0)),
                       mpatches.Patch(facecolor=dfunc.colour_Lambda(5.0),label=r'$\Lambda_1/\Lambda_2 =$'+str(5.0)),
                       mpatches.Patch(facecolor=dfunc.colour_Lambda(10.0),label=r'$\Lambda_1/\Lambda_2 =$'+str(10.0))]
    legend_elements.append(Line2D([0], [0], color='k',linestyle="-",label=r'$\Lambda_1$'))
    legend_elements.append(Line2D([0], [0], color='k',linestyle="--",label=r'$\Lambda_2$'))

    
    plt.legend(handles=legend_elements)

    plt.xlabel(r'Time after collision, $t - t_{collision}$ (h)')
    plt.ylabel(r"<fraction of the contour occupying interface>")

    plt.tight_layout()

    save_path = os.path.join(output_dir, f"average_interface_fraction_ratio.pdf")

    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved average xCOM plot to: {save_path}")
    plt.show()

def average_interface_fraction_time(dirs_averaged,output_dir=DEFAULT_OUTPUT_DIR,log=False):
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
        Lambdas = Lambdas + list(pfunc.plot_average_interfacexsurface(data_dirs)) 

    Lambdas = sorted(set(Lambdas))
    legend_elements = []
    for Lambda in Lambdas:
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label=str(round(Lambda,1))) ]
    plt.legend(handles=legend_elements,title=r'$\Lambda$')

    plt.xlabel(r'Time after collision, $t - t_{collision}$ (h)')
    plt.ylabel(r"Fraction of interface, $ \phi$")

    plt.tight_layout()

    save_path = os.path.join(output_dir, f"surfacefrac.pdf")

    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved average xCOM plot to: {save_path}")
    plt.show()



def COM_Ncells_long(data_dirs, output_dir=DEFAULT_OUTPUT_DIR,log=False,Nlong1=30,Nlong2=30):
    plt.figure(figsize=(5, 3.5))

    if type(data_dirs) == str:
        data_dirs = [data_dirs]
    else:
        data_dirs = list(data_dirs)

    Lambdas = []

    Nlong11 = [130,200,150,150,90,160,270]
    Nlong1 =[170,200,250,100,90,200,270]

    Nlong11_5= [90,110,100,150,160]
    Nlong1_5 = [110,70,80,150,170]
    
    Nlong15= [160,180,50,100,150,150,100,300]
    Nlong5 = [100,100,160,100,250,140,100,150]

    Nlong110= [120,170,150,80,190,100,190,350]
    Nlong10 = [150,150,110,130,100,100,200,120]
    #Nlong10 = [200,200,200,200,200,200,200,200]
    
    Nlong1_combine = Nlong11+Nlong11_5+Nlong15+Nlong110
    Nlong2_combine= Nlong1+Nlong1_5+Nlong5+Nlong10

    for i,data_dir in enumerate(data_dirs):
        #Lambda, k,x = pfunc.COM_N(data_dir,Nlong1_combine[i],Nlong2_combine[i])
        Lambda, k,x = pfunc.COM_N(data_dir,300,300)
        if Lambda == [0]:
            continue
        Lambdas += list(Lambda)
    
    Lambdas = set(Lambdas)
    legend_elements = []
    for Lambda in Lambdas:
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label=r'$\Lambda =$'+str(Lambda)) ]
    plt.legend(handles=legend_elements)

    if log:
        plt.yscale("log")
        plt.xscale("log")

    plt.xlabel(r"$N_{cells}$")
    plt.ylabel(r"$(\vec R_{COM} - \vec r_{0})/A$")
    plt.tight_layout()

    os.makedirs(output_dir, exist_ok=True)
    save_path = os.path.join(output_dir, f"com_movement.pdf")
    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved surface fraction plot to: {save_path}")
    plt.show()

def COM_Ncells_parameter_dist(data_dirs, output_dir=DEFAULT_OUTPUT_DIR):
    plt.figure(figsize=(5, 3.5))

    if type(data_dirs) == str:
        data_dirs = [data_dirs]
    else:
        data_dirs = list(data_dirs)

    Lambdas = []

    Nlong11 = [130,200,150,150,90,160,270]
    Nlong1 =[170,200,250,100,90,200,270]

    Nlong11_5= [90,110,100,150,160]
    Nlong1_5 = [110,70,80,150,170]
    
    Nlong15= [160,180,50,100,150,150,100,300]
    Nlong5 = [100,100,160,100,250,140,100,150]

    Nlong110= [120,170,150,80,190,100,190,350]
    Nlong10 = [150,150,110,130,100,100,200,120]
    #Nlong10 = [200,200,200,200,200,200,200,200]
    
    Nlong1_combine = Nlong11+Nlong11_5+Nlong15+Nlong110
    Nlong2_combine= Nlong1+Nlong1_5+Nlong5+Nlong10

    xpar=[]
    kpar=[]
    Lambdas = []

    for i,data_dir in enumerate(data_dirs):
        #Lambda, k,x = pfunc.COM_N(data_dir,Nlong1_combine[i],Nlong2_combine[i],plot=False)
        Lambda, k,x = pfunc.COM_N(data_dir,200,200,plot=False)
        if Lambda == [0]:
            continue
        Lambdas +=Lambda
        xpar += list(x)
        kpar += list(k)
    data = np.asarray([np.asarray(Lambdas),np.asarray(kpar),np.asarray(xpar)])
    header = ["Lambda","k","x"]
    df = pd.DataFrame(np.transpose(data),columns=header)

    Lambdas = df["Lambda"].unique()
    df_x = df["x"]
    df_k = df["k"]

    N_data = df.shape[0]
    legend_elements=[]

    for Lambda in Lambdas:
        x_av,x_err = dfunc.average_std(df_x[df['Lambda']==Lambda])
        k_av,k_err = dfunc.average_std(df_k[df['Lambda']==Lambda])
        print(f"Scaling Lambda={Lambda}: Rcom(N) = ({round(k_av,4)}+-{round(k_err,4)})*N^({round(x_av,4)}+-{round(x_err,4)})")
        plt.errorbar(Lambda,x_av,yerr=x_err,fmt="x",color=dfunc.colour_Lambda(Lambda))
        plt.scatter(Lambda*np.ones(df_x[df['Lambda']==Lambda].shape[0]),df_x[df['Lambda']==Lambda],color=dfunc.colour_Lambda(Lambda),s=5,alpha=0.3)

        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label=r'$\Lambda =$'+str(Lambda)) ]

    plt.legend(handles=legend_elements)

    #print(f"Scaling total: Rcom(N) = ({df_k.mean()}+-{df_k.std()/np.sqrt(N_data)})*N^({df_x.mean()}+-{df_x.std()/np.sqrt(N_data)})")
    
    plt.xlabel(r"$\Lambda$")
    plt.ylabel(r"Scaling exponent $x$ for $R_{com} = kN^x$")
    plt.show()

    os.makedirs(output_dir, exist_ok=True)
    save_path = os.path.join(output_dir, f"com_movement.pdf")
    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved surface fraction plot to: {save_path}")
    plt.show()

def plot_total_perim_area(data_dirs):
    plt.figure(figsize=(5, 3.5))

    if type(data_dirs) == str:
        data_dirs = [data_dirs]
    else:
        data_dirs = list(data_dirs)

    for data_dir in data_dirs:
        pfunc.plot_perimeter_total_colony(data_dir)

    plt.xlabel("Time (h)")
    plt.ylabel(r"Total perimeter")
    plt.tight_layout()

    plt.show()
    
def average_surfacefrac_ratio(dirs_averaged,output_dir=DEFAULT_OUTPUT_DIR):
    """
    This function should have as input a LIST containing LISTS of the
    directories to be average and plotted at the same time
    e.g. [[file1, file2], [file3, file4]]
    1 and 2 will be averaged and plloted and 3 and 4 will be averaged and
    plotted.
    """
    plt.figure(figsize=(5, 3.5))
    for data_dirs in dirs_averaged:
        pfunc.plot_average_surface_fraction_ratio(data_dirs)

    legend_elements = [mpatches.Patch(facecolor=dfunc.colour_Lambda(1.001),label=str(1.0)),
                       mpatches.Patch(facecolor=dfunc.colour_Lambda(1.5),label=str(1.5)),
                       mpatches.Patch(facecolor=dfunc.colour_Lambda(5.0),label=str(5.0)),
                       mpatches.Patch(facecolor=dfunc.colour_Lambda(10.0),label=str(10.0))]

    
    plt.legend(handles=legend_elements,title=r"$\Lambda$")

    plt.xlabel(r'Time after collision, $t - t_{collision}$ (h)')
    plt.ylabel(r"$Fraction of surface occupied by $\Lambda_2$, $S_2$")

    plt.tight_layout()

    save_path = os.path.join(output_dir, f"average_surfacefrac.pdf")

    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved average xCOM plot to: {save_path}")
    plt.show()

def average_surfacefrac_time(dirs_averaged,output_dir=DEFAULT_OUTPUT_DIR,log=False):
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
        Lambdas = Lambdas + list(pfunc.plot_average_surface_fraction(data_dirs)) 

    legend_elements = [mpatches.Patch(facecolor=dfunc.colour_Lambda(1.001),label=str(1.0)),
                           mpatches.Patch(facecolor=dfunc.colour_Lambda(1.5),label=str(1.5)),
                           mpatches.Patch(facecolor=dfunc.colour_Lambda(5.0),label=str(5.0)),
                           mpatches.Patch(facecolor=dfunc.colour_Lambda(10.0),label=str(10.0))]
    
        
    plt.legend(handles=legend_elements,title=r"$\Lambda$")
    plt.xlabel(r'Time after collision, $t - t_{collision}$ (h)')
    plt.ylabel(r"Fraction of surface occupied by $\Lambda_2, \ S_2$")

    plt.tight_layout()

    save_path = os.path.join(output_dir, f"surfacefrac.pdf")

    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved average xCOM plot to: {save_path}")
    plt.show()



def PerimArea_Ncells(data_dirs, output_dir=DEFAULT_OUTPUT_DIR):
    fig, axes = plt.subplots(1, 2, figsize=(10,5),
                             constrained_layout=True, facecolor='w')

    if type(data_dirs) == str:
        data_dirs = [data_dirs]
    else:
        data_dirs = list(data_dirs)

    Lambdas = []
    xPerim = []
    kPerim = []
    kArea = []
    xArea =[]

    for data_dir in data_dirs:
        Lambda,kP,xP,kA,xA = pfunc.plot_perimeter_area_N(data_dir,ax=axes)
        Lambdas += list(Lambda)
        xPerim += list(xP)
        kPerim += list(kP)
        kArea += list(kA)
        xArea += list(xA)

    data = np.asarray([np.asarray(Lambdas),np.asarray(kPerim),np.asarray(xPerim),np.asarray(kArea),np.asarray(xArea)])
    header = "Lambda\tk_perim\tx_perim\tk_area\tx_area"
    #np.savetxt("Area_Perimeter_exponents", np.transpose(data), delimiter="\t", fmt="%g",header=header)

    Lambdas = set(Lambdas)
    legend_elements = []
    for Lambda in Lambdas:
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label=r'$\Lambda =$'+str(Lambda)) ]

    for i,ax in enumerate(axes):
        ax.legend(handles=legend_elements)
        ax.set_yscale("log")
        ax.set_xscale("log")

        ax.set_xlabel("$N_{cells}$")
        if i == 0:
            ax.set_ylabel(r"P(N) $\mu m$")
        else:
            ax.set_ylabel(r"A(N) $\mu m^2$")
    plt.tight_layout()

    os.makedirs(output_dir, exist_ok=True)
    save_path = os.path.join(output_dir, f"com_movement.pdf")
    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved surface fraction plot to: {save_path}")
    plt.show()

def estimate_PerimArea_growth_mean_exp(filename="Area_Perimeter_exponents"):
    plt.figure(figsize=(5, 3.5))
    df = pd.read_csv('C:\\Users\\lucca\\Desktop\\Friction_Biofilms\\Analysis\\'+str(filename), sep='\t')
    Lambdas = df["Lambda"].unique()
    df_xA = df["x_area"]
    df_kA = df["k_area"]

    df_xP = df["x_perim"]
    df_kP = df["k_perim"]

    N_data = df.shape[0]
    legend_elements=[Line2D([0], [0], color='k',marker="o",label= r"$x_{Perim}$", markersize=15),
                    Line2D([0], [0], color='k',marker="x",label= r"$x_{Area}$", markersize=15)]

    for Lambda in Lambdas:
        xA_av,xA_err = dfunc.average_std(df_xA[df['Lambda']==Lambda])
        kA_av,kA_err = dfunc.average_std(df_kA[df['Lambda']==Lambda])
        print(f"Scaling Lambda={Lambda}: A(N) = ({round(kA_av,4)}+-{round(kA_err,4)})*N^({round(xA_av,4)}+-{round(xA_err,4)})")
        plt.errorbar(Lambda,xA_av,yerr=xA_err,fmt="x",color=dfunc.colour_Lambda(Lambda))

        xP_av,xP_err = dfunc.average_std(df_xP[df['Lambda']==Lambda])
        kP_av,kP_err = dfunc.average_std(df_kP[df['Lambda']==Lambda])
        print(f"Scaling Lambda={Lambda}: P(N) = ({round(kP_av,4)}+-{round(kP_err,4)})*N^({round(xP_av,4)}+-{round(xP_err,4)})")
        plt.errorbar(Lambda,xP_av,yerr=xP_err,fmt="x",color=dfunc.colour_Lambda(Lambda))

        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label=r'$\Lambda =$'+str(Lambda)) ]

    plt.legend(handles=legend_elements)

    print(f"Scaling total: A(N) = ({df_kA.mean()}+-{df_kA.std()/np.sqrt(N_data)})*N^({df_xA.mean()}+-{df_xA.std()/np.sqrt(N_data)})")
    print(f"Scaling total: P(N) = ({df_kP.mean()}+-{df_kP.std()/np.sqrt(N_data)})*N^({df_xP.mean()}+-{df_xP.std()/np.sqrt(N_data)})")
    
    plt.xlabel(r"$\Lambda$")
    plt.ylabel(r"Area and perimeter scaling coefficient")
    plt.show()

def PerimArea_Comparison(data_dirs, output_dir=DEFAULT_OUTPUT_DIR,log=True):
    fig, axes = plt.subplots(1, 4, figsize=(12,3),
                             constrained_layout=True, facecolor='w')

    if type(data_dirs) == str:
        data_dirs = [data_dirs]
    else:
        data_dirs = list(data_dirs)
        
    Lambdas = []
    for data_dir in data_dirs:
        Lambdas = Lambdas + list(pfunc.perim_area_compare_FreeXInteracting(data_dir,axes))
    
    Lambdas = set(Lambdas)
    legend_elements = []
    for Lambda in Lambdas:
        legend_elements = legend_elements+ [mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label=r'$\Lambda =$'+str(Lambda))]

    for i,ax in enumerate(axes):
        ax.legend(handles=legend_elements)
        if log:
            ax.set_yscale("log")
            ax.set_xscale("log")

        ax.set_xlabel(r"$N_{cells}$")
        if i == 0:
            ax.set_ylabel(r"$P/P_{free}$")
        elif i == 1:
            ax.set_ylabel(r"$A/A_{free}$")
        elif i == 2:
            ax.set_ylabel(r"$IQ/IQ_{free}$")
            ax.plot([0,4000],[1,1],"k") #IQ for a circle
        else:
            ax.set_ylabel(r"$P^2/A / P^2_{free}/A_{free}$")
    plt.tight_layout()

    os.makedirs(output_dir, exist_ok=True)
    save_path = os.path.join(output_dir, f"perimeter_area(N).pdf")
    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved Perimeter and Area over N plot to: {save_path}")
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
                    Line2D([0], [0], color='k',marker="o",label= r"$\sigma_{\parallel}$", markersize=15),
                    Line2D([0], [0], color='k',marker="x",label= r"$\sigma_{\perp}$", markersize=15),
                    ]
    for Lambda in Lambdas:
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label=r'$\Lambda =$'+str(Lambda)) ]

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
                    Line2D([0], [0], color='k',marker="x",label= r"$p$", markersize=15),
                    ]
    for Lambda in Lambdas:
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label=r'$\Lambda =$'+str(Lambda)) ]

    plt.legend(handles=legend_elements)
    plt.xlabel("Distance from centre (microns)")
    plt.ylabel(r'pressure and $\alpha$ (Pa m)')
    plt.tight_layout()

    save_path = os.path.join(output_dir, f"pressure_distance.pdf")

    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved pressure vs distance plot to: {save_path}")
    plt.show()

def plot_frac_dist(data_dirs, output_dir=DEFAULT_OUTPUT_DIR):
    plt.figure(figsize=(5, 3.5))

    if type(data_dirs) == str:
        data_dirs = [data_dirs]
    else:
        data_dirs = list(data_dirs)

    ratio = []
    for data_dir in data_dirs:
        ratio = ratio + [pfunc.plot_fraction_distance(data_dir)]

    ratio = set(ratio)
    legend_elements = []
    for r in ratio:
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(r),label=r'ratio = '+str(r)) ]

    plt.legend(handles=legend_elements)
    plt.xlabel(r"$R_{COM}/R_{max}$")
    plt.ylabel(r'Fraction of higher cells')
    plt.tight_layout()

    save_path = os.path.join(output_dir, f"pressure_distance.pdf")

    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved pressure vs distance plot to: {save_path}")
    plt.show()

def av_frac_dist(dirs_averaged, output_dir=DEFAULT_OUTPUT_DIR,com="total"):
    plt.figure(figsize=(5, 3.5))

    ratio = []
    for data_dirs in dirs_averaged:
        ratio = ratio + [pfunc.average_fraction_distance(data_dirs,com=com)]

    ratio = sorted(set(ratio))
    legend_elements = []
    for r in ratio:
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(r),label=str(r)) ]

    plt.legend(handles=legend_elements,title=r"ratio, $ \Lambda$")
    plt.xlabel(r"Radial separation, $\bar{\rho}$")
    if com=="higher":
        plt.ylabel(r'Fraction of lower friction cells, $f_1$')
    elif com=="total":
        plt.ylabel(r'Fraction of higher friction cells, $f_2$')
    plt.xlim((0,1))
    plt.tight_layout()

    save_path = os.path.join(output_dir, f"radial_distribution_cells_from_{com}_com.pdf")

    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved plot to: {save_path}")
    plt.show()

def interface_width(dirs_averaged, output_dir=DEFAULT_OUTPUT_DIR):
    fig, axes = plt.subplots(1, 2, figsize=(10, 5),
                                 constrained_layout=True, facecolor='w')
    w = []
    d_star = []
    ratio = []
    errorw = []
    errord =[]
 
    for data_dirs in dirs_averaged:
        try:
            params,err,r = pfunc.average_fraction_distance(data_dirs,ax=None,com="higher",width=True)
        except:
            continue
        
        ratio += [r]
        w += [params[0]]
        d_star += [params[1]]
        errorw +=[err[0]]
        errord +=[err[1]]

    axes[1].errorbar(ratio,d_star,yerr=errord,marker="x",c="k", linestyle='')
    axes[0].errorbar(ratio,w,yerr=errorw,marker="x",c="k", linestyle='')
    
    axes[1].set_xlabel(r"Ratio, $\Lambda$",fontsize=15)
    axes[0].set_xlabel(r"Ratio, $\Lambda$",fontsize=15)
    axes[0].set_ylabel(r"Width of interface, $w$",fontsize=15)
    axes[1].set_ylabel(r"Multiplying constant, $d^*$",fontsize=15)

    os.makedirs(output_dir, exist_ok=True)
    save_path = os.path.join(output_dir, f"distribution_fraction_radial_parameters.pdf")
    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved plot to: {save_path}")



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

    plt.legend(title=r"$\Lambda_i: \tau_{doub}$")

    plt.xlabel("Time (h)")
    plt.ylabel(r"Number of cells, $\log_2 N(t)$")
    plt.tight_layout()

    save_path = os.path.join(output_dir, f"average_GR.pdf")

    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved average growth rate plot to: {save_path}")
    plt.show()

def average_Rg_time(dirs_averaged,output_dir=DEFAULT_OUTPUT_DIR):
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
        Lambdas = Lambdas + list(pfunc.plot_average_rg(data_dirs)) 

    plt.legend(title=r"$\Lambda_i: \tau_{doub}$")

    plt.xlabel("Time, $t$ (h)")
    plt.ylabel(r"Radius of gyration, $\log_2 R_\mathrm{g}(t)$")
    plt.tight_layout()

    save_path = os.path.join(output_dir, f"average_Rg.pdf")

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

    Lambdas = sorted(set(Lambdas))
    legend_elements = []
    for Lambda in Lambdas:
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label=str(Lambda)) ]
    plt.legend(handles=legend_elements,title=r'$\Lambda_i$')
            

    plt.xlabel("Time (h)")
    plt.ylabel(r"Shape asphericity, $\Delta$")
    plt.tight_layout()
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
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label=r'$\Lambda =$'+str(Lambda)) ]
    plt.legend(handles=legend_elements)

    plt.xlabel("Time (h)")
    plt.ylabel(r"$<\Delta A>$")
    plt.tight_layout()

    save_path = os.path.join(output_dir, f"average_dasphericity.pdf")

    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved average delta asphericity plot to: {save_path}")
    plt.show()

def average_IQ_time(dirs_averaged,output_dir=DEFAULT_OUTPUT_DIR):
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
        Lambdas = Lambdas + list(pfunc.plot_average_COM(data_dirs)) 

    Lambdas = set(Lambdas)
    legend_elements = []
    for Lambda in Lambdas:
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label=r'$\Lambda =$'+str(Lambda)) ]
    plt.legend(handles=legend_elements)

    plt.xlabel("Time (h)")
    plt.ylabel(r"$<x_{COM}-x_{0}>$ (microns)")
    plt.tight_layout()

    save_path = os.path.join(output_dir, f"average_xCOM.pdf")

    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved average xCOM plot to: {save_path}")
    plt.show()



def average_COM_time(dirs_averaged,output_dir=DEFAULT_OUTPUT_DIR,log=False):
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
        Lambdas = Lambdas + list(pfunc.plot_average_COM(data_dirs)) 

    Lambdas = set(Lambdas)
    legend_elements = []
    for Lambda in Lambdas:
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label=r'$\Lambda =$'+str(Lambda)) ]
    plt.legend(handles=legend_elements)

    if log:
        plt.xlabel(r"Time after collision, $t - t_{collision}$ (h)")
        plt.ylabel(r"$<\vec r_{COM} - \vec r_{0}>$")
        plt.yscale("log")
        plt.xscale("log")
    else:
        plt.xlabel(r"Time after collision, $t - t_{collision}$ (h)")
        plt.ylabel(r"$<\vec r_{COM} - \vec r_{0}>$")

    plt.tight_layout()

    save_path = os.path.join(output_dir, f"average_xCOM.pdf")

    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved average xCOM plot to: {save_path}")
    plt.show()

def average_COM_Ncells(dirs_averaged,output_dir=DEFAULT_OUTPUT_DIR,log=False):
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
        Lambdas = Lambdas + list(pfunc.average_COM_N(data_dirs)) 

    Lambdas = set(Lambdas)
    legend_elements = []
    for Lambda in Lambdas:
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label=r'$\Lambda =$'+str(Lambda)) ]
    plt.legend(handles=legend_elements)

    if log:
        plt.yscale("log")
        plt.xscale("log")

    plt.xlabel(r"$N_{cells}$")
    plt.ylabel(r"$\vec r_{COM} - \vec r_{0}$")

    plt.tight_layout()

    save_path = os.path.join(output_dir, f"average_xCOM.pdf")

    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved average xCOM plot to: {save_path}")
    plt.show()


def average_COM_ratio(dirs_averaged,output_dir=DEFAULT_OUTPUT_DIR):
    """
    This function should have as input a LIST containing LISTS of the
    directories to be average and plotted at the same time
    e.g. [[file1, file2], [file3, file4]]
    1 and 2 will be averaged and plloted and 3 and 4 will be averaged and
    plotted.
    """
    plt.figure(figsize=(5, 3.5))
    for data_dirs in dirs_averaged:
        pfunc.plot_average_COM_ratio(data_dirs)

    legend_elements = [mpatches.Patch(facecolor=dfunc.colour_Lambda(1.0),label=r'$\Lambda_1/\Lambda_2 =$'+str(1.0)),
                       mpatches.Patch(facecolor=dfunc.colour_Lambda(1.5),label=r'$\Lambda_1/\Lambda_2 =$'+str(1.5)),
                       mpatches.Patch(facecolor=dfunc.colour_Lambda(2.0),label=r'$\Lambda_1/\Lambda_2 =$'+str(2.0)),
                       mpatches.Patch(facecolor=dfunc.colour_Lambda(5.0),label=r'$\Lambda_1/\Lambda_2 =$'+str(5.0)),
                       mpatches.Patch(facecolor=dfunc.colour_Lambda(10.0),label=r'$\Lambda_1/\Lambda_2 =$'+str(10.0))]
    legend_elements.append(Line2D([0], [0], color='k',linestyle="-",label=r'$\Lambda_1$'))
    legend_elements.append(Line2D([0], [0], color='k',linestyle="--",label=r'$\Lambda_2$'))

    
    plt.legend(handles=legend_elements)

    plt.xlabel(r'Time after collision, $t - t_{collision}$ (h)')
    plt.ylabel(r"$<\vec r_{COM} - \vec r_{0}>$")
    plt.yscale("log")
    plt.xscale("log")
    plt.tight_layout()

    save_path = os.path.join(output_dir, f"average_xCOM.pdf")

    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved average xCOM plot to: {save_path}")
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

    Lambdas = sorted(set(Lambdas))
    legend_elements = [Line2D([0], [0], color='k',linestyle="--",label=r'Prediction')]
    


    #plots growth rate lines
    if lines:
        #plt.plot([0.5,10.5],[0.7,0.7],"--", color="k")
        plt.plot([0.5,10.5],[3.5/3,3.5/3],"--", color="k")

    plt.legend(handles=legend_elements)


    plt.xlabel(r"Friction coefficient, $\Lambda_i$")
    plt.ylabel(r'Doubling time, $\tau_{doub}$ (h)')
    plt.xlim([0.5,5.5])
    plt.title(r"$\mu = 3$ microns/h")

    plt.tight_layout()

    save_path = os.path.join(output_dir, f"tdoub_lambda.pdf")

    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved tdoub vs lambda plot to: {save_path}")
    plt.show()

def plot_fraction_y(data_dirs, width=0, output_dir=DEFAULT_OUTPUT_DIR):
    plt.figure(figsize=(5, 3.5))

    if type(data_dirs) == str:
        data_dirs = [data_dirs]
    else:
        data_dirs = list(data_dirs)

    Lambdas = []
    timestep = 110 #timestep to consider!

    for data_dir in data_dirs:
        files = pfunc.get_file_paths(data_dir)
        filepath = files[-1]
        print(filepath)
        Lambdas = Lambdas + list(pfunc.plot_yfraction(filepath,width))
        #for filepath in files:
            #if str(timestep) in filepath:
                #Lambdas = Lambdas + list(pfunc.plot_yfraction(filepath))
                #print(filepath)

    Lambdas = set(Lambdas)
    legend_elements = []
    for Lambda in Lambdas:
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label=R'$\Lambda =$'+str(Lambda)) ]

    plt.legend(handles=legend_elements)
    plt.xlabel("y position (microns)")
    plt.ylabel(r"fraction of cells")
    plt.tight_layout()

    save_path = os.path.join(output_dir, f"fraction_y.pdf")

    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved fraction over colony plot to: {save_path}")
    plt.show()

def plot_yfraction_repeats(dirs_averaged, width=0, output_dir=DEFAULT_OUTPUT_DIR):
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
        Lambdas = Lambdas + list(pfunc.yfraction_repeats(data_dirs,width)) 

    Lambdas = set(Lambdas)
    legend_elements = [Line2D([0], [0], color='k',linestyle="--",label=r'$\Lambda_1$')]
    for Lambda in Lambdas:
        if Lambda == 1.0:
            continue
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label=r'$\Lambda =$'+str(Lambda)) ]
    plt.legend(handles=legend_elements)

    plt.xlabel("y position (microns)")
    plt.ylabel("fraction of cells")
    plt.tight_layout()

    save_path = os.path.join(output_dir, f"frac_ytotal.pdf")

    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved fraction vs position plot to: {save_path}")
    plt.show()

def plot_fraction_x_evolution(data_dir, width=0, output_dir=DEFAULT_OUTPUT_DIR):
    plt.figure(figsize=(5, 3.5))
    Lambdas = []

    files = pfunc.get_file_paths(data_dir)
    filepaths_consider = [1,15,30,45,80]
    markers = [".","x","v","*","p"]
    legend_elements = []
    for i,timestep in enumerate(filepaths_consider):
        try:
            print(files[timestep])
            Lambdas = Lambdas + list(pfunc.plot_xfraction(files[timestep],width,markers[i]))
            legend_elements +=Line2D([0], [0], color='k',marker=markers[i],label= r"t="+str(round(timestep*0.1,1))+"h", markersize=15),
        except:
            print("timestep: "+str(timestep*0.1)+"(h) not found")

    Lambdas = set(Lambdas)
    
    for Lambda in Lambdas:
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label=R'$\Lambda =$'+str(Lambda)) ]

    plt.legend(handles=legend_elements)
    plt.xlabel("x position (microns)")
    plt.ylabel(r"fraction of cells")
    plt.tight_layout()

    save_path = os.path.join(output_dir, f"fraction_y.pdf")

    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved fraction over colony plot to: {save_path}")
    plt.show()

from analytical_model_functions import evolve_time
def plot_wrapping_analytical_varying(dirs_averaged,output_dir=DEFAULT_OUTPUT_DIR):
    """
        This function should have as input a LIST containing LISTS of the
        directories to be average and plotted at the same time
        e.g. [[file1, file2], [file3, file4]]
        1 and 2 will be averaged and plloted and 3 and 4 will be averaged and
        plotted.
    """
    plt.figure(figsize=(5, 3.5))
    colours = ["k",dfunc.colour_Lambda(10.0),"purple","b"]
    
    for i, data_dirs in enumerate(dirs_averaged):
        pfunc.surface_frac_analytical_param(data_dirs,colours[i])


    dt = 0.0001
    tend = 10
    time = np.arange(0,tend,dt)

    frac_s0 = 0.1
    
    av_l=[5,5,5,6.5]
    t0 = [1,3,5,3]
    tau=[3,1.8,1.35,2.15]

    x = [3,5,7,5]
    fraction = np.zeros((len(time),len(x)))
    for i in range(0,len(x)):
        r0 = np.sqrt(2)*2**(t0[i]/tau[i])
        sgreen,sblue = evolve_time(time,frac_s0=frac_s0,dt=dt,mu=x[i],av_l=av_l[i],R0=r0,tau=tau[i])

        fraction[:,i] = sblue/sgreen
        plt.plot(time,fraction[:,i],c=colours[i])

    legend_elements = []
    for i,c in enumerate(colours):
        legend_elements = legend_elements+ [ mpatches.Patch(facecolor=c,label=f"{x[i]}, {av_l[i]}") ]

    legend_elements = legend_elements + [Line2D([0], [0], color='k',linestyle="-",label=r'Model'),Line2D([0], [0], color='k',marker="x",label= r"Simulation", markersize=15)]

    plt.legend(handles=legend_elements,title=r'$\mu, \ \bar{\ell} $')
    
    plt.xlabel(r'Time after collision, $t - t_{collision}$ (h)')
    plt.ylabel(r"Fraction of interface, $\phi$")
    
    plt.tight_layout()
    
    save_path = os.path.join(output_dir, f"surfacefrac_paramschange.pdf")
    plt.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved average xCOM plot to: {save_path}")
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

"""if __name__ == "__main__":
    main()"""