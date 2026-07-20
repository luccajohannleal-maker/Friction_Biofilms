import os
import re
import glob
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
import analytical_model_functions as analytical
ut.setMPL()


global data_loc,Free0_8,Free1,Free1_5,Free5,Free10

data_loc = "C:\\Users\\lucca\\Desktop\\GeneratedOutput\\SimOutput\\data\\"

Free0_8 = [data_loc+"FreeGrow\\Lambda0_8\\repeat0"]
Free1 = pfunc.repeat_files(data_loc+"FreeGrow\\Lambda1",5)+ pfunc.repeat_files(data_loc+"FreeGrow\\stress\\Lambda1",3)
Free1_5 = [data_loc+"FreeGrow\\stress\\Lambda1_5\\repeat0",data_loc+"FreeGrow\\Lambda1_5\\repeat0"]
Free5 = pfunc.repeat_files(data_loc+"FreeGrow\\Lambda5",5)+ pfunc.repeat_files(data_loc+"FreeGrow\\stress\\Lambda5",4)
Free10 = pfunc.repeat_files(data_loc+"FreeGrow\\Lambda10",1)+pfunc.repeat_files(data_loc+"FreeGrow\\stress\\Lambda10",4)




def figure1():
    initial_time= 0
    files = pfunc.get_file_paths(data_loc+"FreeGrow\\stress\\Lambda1_5\\repeat0")
    if len(files) < 6:
        selected_files = files
    else:
        selected_indices = np.linspace(initial_time, len(files) - 1, 6, dtype=int)
        selected_files = [files[i] for i in selected_indices]
    for i,file in enumerate(selected_files):
        ax = plt.subplot2grid((3,6),(0,i))
        dfunc.plotCells(ax, file)
        match = re.search(r'biofilm_(\d+)\.dat$', os.path.basename(file))
        if match:
            frame_number = int(match.group(1))
            time_hours = frame_number * 0.1
            ax.set_title(f"{time_hours:.1f} h", color='k', fontsize=12)




    b = plt.subplot2grid((3,6),(1,0),rowspan=2,colspan=2)
    dirs_averaged = [Free1,Free1_5,Free5,Free10]
    Lambdas = []
    tdoub = []
    terr = []

    for data_dirs in dirs_averaged:
        Lambda,t,err = pfunc.plot_average_growth(data_dirs,ax=b)
        Lambdas += list(Lambda) 
        tdoub += [round(t,2)]
        terr += [round(err,2)]

    b.plot([0,8],dfunc.doubling_linear_growth(np.array([0,8]),0.7,0),"r--")

    legend_elements = []
    for i,Lambda in enumerate(Lambdas):
        legend_elements += [mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label=r'$\Lambda =$'+str(Lambda)+r", $t_{d}=$"+str(tdoub[i])+r"$\pm$"+str(terr[i]))]

    
    b.legend(handles=legend_elements)

    b.set_xlabel("Time (h)")
    b.set_ylabel(r"$log_2 <N(t)>$")



    c = plt.subplot2grid((3,6),(1,2),rowspan=2,colspan=2)
    Lambdas = []
    tdoub = []
    terr = []

    for data_dirs in dirs_averaged:
        Lambda,t,err = pfunc.plot_average_growth(data_dirs,ax=c)
        Lambdas += list(Lambda) 
        tdoub += [round(t,2)]
        terr += [round(err,2)]

    c.plot([0,8],dfunc.doubling_linear_growth(np.array([0,8]),1.4,0),"r--")

    legend_elements = []
    for i,Lambda in enumerate(Lambdas):
        legend_elements += [mpatches.Patch(facecolor=dfunc.colour_Lambda(Lambda),label=r'$\Lambda =$'+str(Lambda)+r", $2t_{d}=$"+str(tdoub[i])+r"$\pm$"+str(terr[i]))]
    c.legend(handles=legend_elements)





    d = plt.subplot2grid((3,6),(1,4),rowspan=2,colspan=2)

    plt.show()





def figureX():
    fig, axes = plt.subplots(1,2, figsize=(8,3), constrained_layout=True, facecolor='w')
    interacting10 = pfunc.repeat_files(data_loc+"Interacting_colonies\\stress\\Lambda1AND10",5)+pfunc.repeat_files(data_loc+"Interacting_colonies\\sameIC\\Lambda1AND10",2)
    interacting5 = pfunc.repeat_files(data_loc+"Interacting_colonies\\stress\\Lambda1AND5",4)+pfunc.repeat_files(data_loc+"Interacting_colonies\\sameIC\\Lambda1AND5",2) + pfunc.repeat_files(data_loc+"Interacting_colonies\\long",1)

    t0_10=[]
    t0_5=[]
    tend10 = 0
    tend5 = 0 

    for data_dir in interacting10:
        files = pfunc.get_file_paths(data_dir)
        t_collision = dfunc.colonies_collided(files)
        t0_10.append(t_collision)

        final_path = files[-1]
        tfinal= final_path[-9:-4]
        if tend10 < round(float(tfinal)*0.1-t_collision,1):
            tend10 = round(float(tfinal)*0.1-t_collision,1)
        
    for data_dir in interacting5:
        files = pfunc.get_file_paths(data_dir)
        t_collision = dfunc.colonies_collided(files)
        t0_5.append(t_collision)

        final_path = files[-1]
        tfinal= final_path[-9:-4]
        if tend5 < round(float(tfinal)*0.1-t_collision,1):
            tend5 = round(float(tfinal)*0.1-t_collision,1)

    analytical.plot_average_inner_fraction(t0_10,tend10,10,axes[0])
    analytical.plot_average_inner_fraction(t0_5,tend5,5,axes[1])

    axes[0].set_title(r"$\Lambda=10$")
    axes[1].set_title(r"$\Lambda=5$")
    axes[0].legend()
    axes[1].legend()
    for ax in axes.flat:
        ax.set(xlabel=r'$t - t_{collision}$', ylabel=r'<interface fraction>')

    # Hide x labels and tick labels for top plots and y ticks for right plots.
    for ax in axes.flat:
        ax.label_outer()

    plt.show()



figureX()