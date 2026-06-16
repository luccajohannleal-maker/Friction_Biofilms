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


def plotCells_channel(ax, file, width=120):
        dat = pd.read_csv(file, sep='\t')
        cells = ut.getCells(file)
        x_center,y_center = 0, 0

        y_top= width/2
        y_bottom = -width/2

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

        fcp.addAllCellsToPlot(cells, ax, ax_rng=maxx - minx, show_id=False, ec='w')


        #Plots walls
        scale = 1
        wall_color = 'k'
        ax.plot([y_bottom*1.5/scale, y_top*1.5/scale], [y_top/scale, y_top/scale], color=wall_color, alpha=0.6)
        ax.plot([y_bottom*1.5/scale, y_top*1.5/scale], [y_bottom/scale, y_bottom/scale], color=wall_color, alpha=0.6)
        ax.plot([y_bottom*1.5/scale, y_bottom*1.5/scale], [y_top/scale, y_bottom/scale],"--", color=wall_color, alpha=0.6)
        ax.plot([y_top*1.5/scale, y_top*1.5/scale], [y_top/scale, y_bottom/scale],"--", color=wall_color, alpha=0.6)


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
    dist=[]
    for cell in cells:
        dist.append(np.sqrt(cell.pos_x**2 + cell.pos_y**2 + cell.pos_z**2))
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
    eig_val1,eig_val2 = np.linalg.eigvals(Tensor)

    return Rg,[eig_val1,eig_val2]

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

    Rg_1 = []
    Rg_2 = []
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

        Lambdas = sorted(Lambdas)
        if len(Lambdas) != 1: #If only one Lambda, only Rg_tot is required, so skip the rest
            cells1 = find_Lambda_cells(cells,Lambdas[0])
            Rg_1.append(Gyration_values(cells1)[0])
        
            cells2 = find_Lambda_cells(cells,Lambdas[1])
            Rg_2.append(Gyration_values(cells2)[0])

    return time_steps, Rg_tot, Lambdas, Rg_1, Rg_2

def RgLambda_time_est(Rg, time_step=0.1):
    Rg = np.asarray(Rg)
    Rg[Rg<1] = 1
    linear_rg = np.log2(Rg)

    t = np.arange(0, (len(Rg) - 0.5)* time_step, time_step) #0.5 INCLUDED TO AVOID ANY SMALL ERRORS IN FLOAT POINTS
    popt, pcov= curve_fit(Rg_Lambda_growth, t, linear_rg, p0=(0,2.5)) # tau_rg = 2*tdouble - see notes
    return popt,np.sqrt(np.diag(pcov))

def Rg_Lambda_growth(t,R0,tau):
    return R0 + t/tau

def Gyr_eig_time(data_dir):
    file_pattern = os.path.join(data_dir, "biofilm_*.dat")
    files = sorted(glob.glob(file_pattern), key=lambda x: int(re.search(r'(\d+)', x).group(1)))
    
    Rg1_eig = []
    Rg2_eig = []
    time_steps = []
    Rg_tot_eig = []

    for file_path in files:
        match = re.search(r'(\d+)', os.path.basename(file_path))
        if not match:
            continue

        time_step = int(match.group(1))
        cells = ut.getCells(file_path)
        Lambdas = find_Lambdas(cells)
        Lambdas = sorted(Lambdas)

        Rg_tot_eig.append(Gyration_values(cells)[1])
        time_steps.append(time_step * 0.1)

        if len(Lambdas) != 1: #If only one Lambda, only Rg_tot is required, so skip the rest
            cells1 = find_Lambda_cells(cells,Lambdas[0])
            Rg1_eig.append(Gyration_values(cells1)[1])

            cells2 = find_Lambda_cells(cells,Lambdas[1])
            Rg2_eig.append(Gyration_values(cells2)[1]) 
    
    return time_steps, np.asarray(Rg_tot_eig), Lambdas, np.asarray(Rg1_eig), np.asarray(Rg2_eig)

def calc_asphericity(data_dir):
    file_pattern = os.path.join(data_dir, "biofilm_*.dat")
    files = sorted(glob.glob(file_pattern), key=lambda x: int(re.search(r'(\d+)', x).group(1)))
    
    asph1 = []
    asph2 = []
    time_steps = []


    for file_path in files:
        match = re.search(r'(\d+)', os.path.basename(file_path))
        if not match:
            continue

        time_step = int(match.group(1))
        cells = ut.getCells(file_path)
        Lambdas = find_Lambdas(cells)
        Lambdas = sorted(Lambdas)
        time_steps.append(time_step * 0.1)

        if len(Lambdas) == 1: #If only one Lambda, only Rg_tot is required, so skip the rest
            cells1 = find_Lambda_cells(cells,Lambdas[0])
            rg,eigenv = Gyration_values(cells1)
            if eigenv[0] and eigenv[1] == 0: #one cell
                asph1.append(1) 
            else:
                asph1.append((eigenv[0]-eigenv[1])**2 /((eigenv[0]+eigenv[1])**2))
            
        else:
            cells1 = find_Lambda_cells(cells,Lambdas[0])
            rg1,eigenv1 = Gyration_values(cells1)

            cells2 = find_Lambda_cells(cells,Lambdas[1])
            rg2,eigenv2 = Gyration_values(cells2)

            if eigenv1[0] and eigenv1[1] == 0: #one cell
                asph1.append(1) 
            else:
                asph1.append((eigenv1[0]-eigenv1[1])**2 /((eigenv1[0]+eigenv1[1])**2))
        
            if eigenv2[0] and eigenv2[1] == 0: #one cell
                asph2.append(1) 
            else:
                asph2.append((eigenv2[0]-eigenv2[1])**2 /((eigenv2[0]+eigenv2[1])**2))

    return time_steps, Lambdas, asph1, asph2


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
    elif Lambda == 1.001:
        return "#1e00ff"
    elif Lambda == 10.001:
        return "#4f7a55"
    elif Lambda < 1.0:
        return (Lambda, 0, 0, 1)
    elif Lambda > 1.0:
        return (0,(Lambda-1)/(10-1), 0, 1)

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

def doubling_linear_growth(t, t_doub, N0):
        return N0 + t/t_doub
def doubling_exp_growth(t, t_doub):
        return  2**(t/t_doub)


def estimate_exp_growth_rate(counts, time_step=0.1):
    """
        Parameters:
            counts: list of int
                cell counts over time
            time_step: int
                time step of the simulation

        Returns:
            doubling time parameter popt[0] and its error np.sqrt(np.diag(pcov))[0]
    """
    t = np.arange(0, (len(counts) - 0.5)* time_step, time_step) #0.5 INCLUDED TO AVOID ANY SMALL ERRORS IN FLOAT POINTS
    popt, pcov= curve_fit(doubling_exp_growth, t, counts, p0=(0.7))

    return popt,np.sqrt(np.diag(pcov))

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
    linear_counts = np.log2(counts)

    t = np.arange(0, (len(counts) - 0.5)* time_step, time_step) #0.5 INCLUDED TO AVOID ANY SMALL ERRORS IN FLOAT POINTS
    popt, pcov= curve_fit(doubling_linear_growth, t, linear_counts, p0=(1.25, 0))

    return popt,np.sqrt(np.diag(pcov))


def find_stress(cells_data): #divide by two as expected from paper "Growing microdomains"
    par_stress = abs(cells_data["st_par"])/2
    perp_stress = abs(cells_data["st_perp"])/2
    shear1 = abs(cells_data["st_shear1"])/2
    shear2 = abs(cells_data["st_shear2"])/2

    return par_stress, perp_stress, shear1, shear2

def find_pressure(cells_data):
   par,perp,shear1,shear2 = find_stress(cells_data)

   par = np.asarray(par)
   perp = np.asarray(perp)

   p = abs(par+perp)/2 #pressure
   alpha = abs(par - perp) #deviatoric active stress
   
   return p, alpha


