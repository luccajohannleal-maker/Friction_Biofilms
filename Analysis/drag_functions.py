import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import fastCellPlotting as fcp
import utilities as ut
import re
import glob
import os

from scipy.optimize import curve_fit


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

        cells = centerCells(cells)

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

def centerBiofilm(cells):
    """
        Parameters:
            cells: list of ChainingRodShapedBacterium
                cells to be centred

        Returns:
            position of center of cell
    """
    
    x_centre = np.mean([cell.pos_x for cell in cells])
    y_centre = np.mean([cell.pos_y for cell in cells])
    z_centre = np.mean([cell.pos_y for cell in cells])

    return x_centre, y_centre, z_centre

def centerCells(cells):
    """
        Parameters:
            cells: list of ChainingRodShapedBacterium
                cells to be centred

        Returns:
            cells with positions shifted to be centred around (0,0)
    """
    x_centre, y_centre, z_centre = centerBiofilm(cells)

    for cell in cells:
        cell.pos_x -= x_centre
        cell.pos_y -= y_centre
        cell.pos_z -= z_centre
    return cells

def distance_from_origin(cells):

    cells = centerCells(cells)
    dist=np.zeros(len(cells))
    for i in range(0,len(cells)):
        dist[i] = np.sqrt(cells[i].pos_x**2+cells[i].pos_y**2+cells[i].pos_z**2)

    return dist

def max_radius(data_dir):
    file_pattern = os.path.join(data_dir, "biofilm_*.dat")
    files = sorted(glob.glob(file_pattern), key=lambda x: int(re.search(r'(\d+)', x).group(1)))

    max_distances=[]
    time_steps = []

    for file_path in files:
        match = re.search(r'(\d+)', os.path.basename(file_path))
        if not match:
            continue
        
        time_step = int(match.group(1))
        cells = ut.getCells(file_path)

        time_steps.append(time_step * 0.1)
        max_distances.append(max(distance_from_origin(cells)))

    return time_steps, max_distances

def Gyration_values(cells):
    """
        Parameters:
            cells: list of ChainingRodShapedBacterium
                cells to be analysed

        Returns:
            radius of gyration of the total cell population
    """
    Tensor = Gyration_tensor_2D(cells)

    Rg=np.sqrt(np.trace(Tensor))
    eig_vals,eig_vect = np.linalg.eigvals(Tensor)

    return Rg,eig_vals


def Gyration_tensor_2D(cells):
    Gyr_tensor = np.zeros((2,2))
    cells = centerCells(cells)
    x,y = position_cells(cells)
    Gyr_tensor[0,0]= np.mean(np.dot(x,x))
    Gyr_tensor[1,1]= np.mean(np.dot(y,y))
    Gyr_tensor[0,1]= np.mean(np.dot(x,y))
    Gyr_tensor[1,0]= np.mean(np.dot(x,y))
    return Gyr_tensor


def position_cells(cells):
    """
        Parameters:
            cells: list of ChainingRodShapedBacterium
                cells to be analysed
            Lambda: float
                non-dimensionalised drag coefficient

        Returns:
            radius of gyration of the cell population for a given Lambda
    """
    x_cell,y_cell = [],[]

    for cell in cells:
            x_cell.append(cell.pos_x)
            y_cell.append(cell.pos_y)

    x_cell = np.array(x_cell)
    y_cell = np.array(y_cell)
    return x_cell,y_cell

def RgLambda_time(data_dir):
    file_pattern = os.path.join(data_dir, "biofilm_*.dat")
    files = sorted(glob.glob(file_pattern), key=lambda x: int(re.search(r'(\d+)', x).group(1)))
    

    Rg_Lambda1 = []
    Rg_Lambdanot1 = []
    time_steps = []
    Rg_tot = []

    for file_path in files:
        match = re.search(r'(\d+)', os.path.basename(file_path))
        if not match:
            continue

        time_step = int(match.group(1))
        cells = ut.getCells(file_path)
        Lambdas = find_Lambdas(cells)

        Rg_tot.append(Gyration_values(cells)[0])
        time_steps.append(time_step * 0.1)

        if len(Lambdas) != 1: #If only one Lambda, only Rg_tot is required, so skip the rest
            for Lambda in Lambdas:
                if Lambda == 1.0:
                    cells1 = find_Lambda_cells(cells,1)
                    Rg_Lambda1.append(Gyration_values(cells1)[0])
                else:
                    cellsnot1 = find_Lambda_cells(cells,Lambda)
                    Rg_Lambdanot1.append(Gyration_values(cellsnot1)[0]) 
    
    return time_steps, Rg_tot, Lambdas, Rg_Lambda1, Rg_Lambdanot1

def find_Lambda_cells(cells, Lambda=1.0):
    """
        Parameters:
            cells: list of ChainingRodShapedBacterium
                cells to be analysed
            Lambda: float
                cell drag coefficient/ base drag coefficient
        Returns:
            list of cells with given Lambda
    """

    if type(cells) == pd.DataFrame:
        filtered_Lambda_cells = cells[cells["Lambda"] == Lambda]
    
    elif type(cells) == list:
        filtered_Lambda_cells = [cell for cell in cells if cell.Lambda == Lambda]
        
    return filtered_Lambda_cells

def find_Lambdas(cells):
    """
        Parameters:
            cells: list of ChainingRodShapedBacterium
                cells to be analysed or a DataFrame containing cell data

        Returns:
            list of unique Lambda values in the cell population
    """

    if type(cells) == pd.DataFrame:
        Lambdas = cells["Lambda"].unique()
        return np.array(Lambdas)
    

    if type(cells) == list:
        Lambdas = set()
        for cell in cells:
            Lambdas.add(cell.Lambda)
        return np.array(list(Lambdas))

def colour_Lambda(Lambda):
    """
        Parameters:
            Lambda: float
                cell drag coefficient/ base drag coefficient

        Returns:
            colour for plotting cell with given Lambda

            If Lambda is 1.0, return light blue. 
            Otherwise, return a shade of red based on the Lambda value. 
            Darker red for higher Lambda, lighter red for lower Lambda.

    """
    if Lambda == 1.0:
        return "#00ffff"
    elif Lambda < 1.0:
        return (255*Lambda/255, 0, 0, 1)
    elif Lambda > 1.0:
        return (0, 255/(Lambda*255), 0, 1)

def counts(data_dir):
    file_pattern = os.path.join(data_dir, "biofilm_*.dat")
    files = sorted(glob.glob(file_pattern), key=lambda x: int(re.search(r'(\d+)', x).group(1)))

    time_steps = []
    Lambda1_counts = []
    Lambdanot1_counts = []

    for file_path in files:
        match = re.search(r'(\d+)', os.path.basename(file_path))
        if not match:
            continue
        time_step = int(match.group(1))
        df = pd.read_csv(file_path, sep="\t")
        Lambdas = find_Lambdas(df)

        time_steps.append(time_step * 0.1)

        for Lambda in Lambdas:
            if Lambda == 1.0:
                Lambda1_count = find_Lambda_cells(df, 1.0).shape[0]
                Lambda1_counts.append(Lambda1_count)
            else:
                Lambdanot1_count =find_Lambda_cells(df, Lambda).shape[0]
                Lambdanot1_counts.append(Lambdanot1_count)

    return time_steps, Lambda1_counts, Lambdanot1_counts, Lambdas

def estimate_growth_rate(counts, time_step=0.1):
    """
        Parameters:
            counts: list of int
                cell counts over time
            time_step: int
                time step of the simulation

        Returns:
            doubling time parameter popt[0] and its error np.sqrt(np.diag(pcov))[0]
    """
    def doubling_growth(t, t_doub):
        return 2**(t/t_doub)
    
    t = np.arange(0, len(counts) * time_step, time_step)
    popt, pcov= curve_fit(doubling_growth, t, counts,p0=0.175)

    return popt[0],np.sqrt(np.diag(pcov))[0]


"""files = ["C:\\Users\\lucca\\Desktop\\GeneratedOutput\\SimOutput\\Free_Growing\\comparison\\repeat0","C:\\Users\\lucca\\Desktop\\GeneratedOutput\\SimOutput\\Free_Growing\\comparison\\repeat1","C:\\Users\\lucca\\Desktop\\GeneratedOutput\\SimOutput\\Free_Growing\\comparison\\repeat2"]

t08,d08 = max_radius(files[0])
t1,d1 = max_radius(files[1])
t15,d15 = max_radius(files[2])
plt.plot(t08, d08, 'o', color=colour_Lambda(0.8), label="Lambda =" + str(0.8))
plt.plot(t1, d1, 'o', color=colour_Lambda(1), label="Lambda =" + str(1))
plt.plot(t15, d15, 'o', color=colour_Lambda(1.5), label="Lambda =" + str(1.5))
plt.ylabel("max distance from origin")


time_steps08, Lambda08_counts, Lambdanot08_counts, Lambda08 = counts(files[0]) 
time_steps1, Lambda1_counts, Lambdanot1_counts, Lambda1 = counts(files[1])
time_steps15, Lambda15_counts, Lambdanot15_counts, Lambda15 = counts(files[2])

growth_rate_Lambda1,err_Lambda1 = estimate_growth_rate(Lambda1_counts, 0.1)
growth_rate_Lambda08, err_Lambda08 = estimate_growth_rate(Lambdanot08_counts, 0.1)
growth_rate_Lambda15, err_Lambda15 = estimate_growth_rate(Lambdanot15_counts, 0.1)

print(f"Estimated growth rate for Lambda=1: {growth_rate_Lambda1}+-{err_Lambda1}")
print(f"Estimated growth rate for Lambda=0.8: {growth_rate_Lambda08}+-{err_Lambda08}")
print(f"Estimated growth rate for Lambda=1.5: {growth_rate_Lambda15}+-{err_Lambda15}")

plt.figure(figsize=(5, 3.5))
plt.plot(time_steps08, Lambdanot08_counts, 'o', color=colour_Lambda(0.8), label=" Counts Lambda =" + str(0.8))
plt.plot(time_steps08, np.exp(growth_rate_Lambda08 * np.array(time_steps08)), '-', color=colour_Lambda(0.8), label="Fit Lambda=" + str(0.8) + ": N(t)= 2^(t/"+str(round(growth_rate_Lambda08, 2)) +")")
plt.plot(time_steps1, Lambda1_counts, 'o', color=colour_Lambda(1), label=" Counts Lambda =" + str(1))
plt.plot(time_steps1, np.exp(growth_rate_Lambda1 * np.array(time_steps1)), '-', color=colour_Lambda(1), label="Fit Lambda=" + str(1) + ": N(t)= 2^(t/"+str(round(growth_rate_Lambda1, 2)) +")")
plt.plot(time_steps15, Lambdanot15_counts, 'o', color=colour_Lambda(1.5), label=" Counts Lambda =" + str(1.5))
plt.plot(time_steps15, np.exp(growth_rate_Lambda15 * np.array(time_steps15)), '-', color=colour_Lambda(1.5), label="Fit Lambda=" + str(1.5) + ": N(t)= 2^(t/"+str(round(growth_rate_Lambda15, 2)) +")")

plt.ylabel("Count N(t)")
plt.yscale("log")


time_steps08, Rg_08, Lambdas, Rg_Lambda1, Rg_Lambdanot1=RgLambda_time(files[0])
time_steps1, Rg_1, Lambdas, Rg_Lambda1, Rg_Lambdanot1=RgLambda_time(files[1])
time_steps15, Rg_15, Lambdas, Rg_Lambda1, Rg_Lambdanot1=RgLambda_time(files[2])

plt.plot(time_steps08, Rg_08, 'o', color=colour_Lambda(0.8), label="Lambda =" + str(0.8))
plt.plot(time_steps1, Rg_1, 'o', color=colour_Lambda(1), label="Lambda =" + str(1))
plt.plot(time_steps15, Rg_15, 'o', color=colour_Lambda(1.5), label="Lambda =" + str(1.5))
plt.ylabel("Rg")

plt.xlabel("Time (h)")
plt.legend()
plt.tight_layout()
plt.show()
"""


"""files = ["C:\\Users\\lucca\\Desktop\\GeneratedOutput\\SimOutput\\Free_Growing\\Lambda1\\repeat1","C:\\Users\\lucca\\Desktop\\GeneratedOutput\\SimOutput\\Free_Growing\\Lambda5\\repeat0"]

plt.figure(figsize=(5, 3.5))
t1,d1 = max_radius(files[0])
t15,d15 = max_radius(files[1])
plt.plot(t1, d1, 'o', color=colour_Lambda(1), label="Lambda =" + str(1))
plt.plot(t15, d15, 'o', color=colour_Lambda(5), label="Lambda =" + str(5))
plt.ylabel("max distance from origin")
plt.xlabel("Time (h)")
plt.legend()
plt.tight_layout()
plt.show()


time_steps1, Lambda1_counts, Lambdanot1_counts, Lambda1 = counts(files[0])
time_steps15, Lambda15_counts, Lambdanot15_counts, Lambda15 = counts(files[1])

growth_rate_Lambda1,err_Lambda1 = estimate_growth_rate(Lambda1_counts, 0.1)
growth_rate_Lambda15, err_Lambda15 = estimate_growth_rate(Lambdanot15_counts, 0.1)

print(f"Estimated growth rate for Lambda=1: {growth_rate_Lambda1}+-{err_Lambda1}")
print(f"Estimated growth rate for Lambda=5: {growth_rate_Lambda15}+-{err_Lambda15}")

plt.figure(figsize=(5, 3.5))
plt.plot(time_steps1, Lambda1_counts, 'o', color=colour_Lambda(1), label=" Counts Lambda =" + str(1))
plt.plot(time_steps1, np.exp(growth_rate_Lambda1 * np.array(time_steps1)), '-', color=colour_Lambda(1), label="Fit Lambda=" + str(1) + ": N(t)= exp("+str(round(growth_rate_Lambda1, 2)) +"*t)")
plt.plot(time_steps15, Lambdanot15_counts, 'o', color=colour_Lambda(5), label=" Counts Lambda =" + str(5))
plt.plot(time_steps15, np.exp(growth_rate_Lambda15 * np.array(time_steps15)), '-', color=colour_Lambda(5), label="Fit Lambda=" + str(5) + ": N(t)= exp("+str(round(growth_rate_Lambda15, 2)) +"*t)")

plt.ylabel("Count N(t)")
plt.yscale("log")
plt.xlabel("Time (h)")
plt.legend()
plt.tight_layout()
plt.show()



plt.figure(figsize=(5, 3.5))
time_steps1, Rg_1, Lambdas, Rg_Lambda1, Rg_Lambdanot1=RgLambda_time(files[0])
time_steps15, Rg_15, Lambdas, Rg_Lambda1, Rg_Lambdanot1=RgLambda_time(files[1])

plt.plot(time_steps1, Rg_1, 'o', color=colour_Lambda(1), label="Lambda =" + str(1))
plt.plot(time_steps15, Rg_15, 'o', color=colour_Lambda(5), label="Lambda =" + str(5))
plt.ylabel("Rg")





plt.xlabel("Time (h)")
plt.legend()
plt.tight_layout()
plt.show()"""