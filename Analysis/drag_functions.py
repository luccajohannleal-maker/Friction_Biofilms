import numpy as np
import pandas as pd
import fastCellPlotting as fcp
import utilities as ut
import re
import glob
import os
import matplotlib.pyplot as plt
from shapely.geometry import Polygon


from DistributionFunctions import computeColonyContour
from scipy.optimize import curve_fit
from generalPlotting import addNematicDirector


def plotCells(ax, file,director=False,defects=False):
        streamplot=(True & director)
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
        Lambdas = find_Lambdas(cells)

        fcp.addAllCellsToPlot(cells, ax, ax_rng=maxx - minx, show_id=False, ec='w')
        #if len(cells) > 200:
            #surface_fraction(cells,ax)

        if director:
            q_name=f"{file[:-4].replace('/','_')}_Q.npy"
            addNematicDirector(ax,cells,q_name,streamplot=streamplot,dr=5)
        
        if defects:
            x_pos = np.asarray(dat["pos_x"])
            y_pos = np.asarray(dat["pos_y"])
            x_or = np.asarray(dat["ori_x"])
            y_or = np.asarray(dat["ori_y"])
            list_defects,s,phi = locate_nematic_defects_vector(x_pos,y_pos,x_or,y_or)
            print("n defects =", len(list_defects))
            for defect in list_defects:
                pos = defect["pos"]
                if float(defect["charge"]) == 0.5:
                    ax.scatter(pos[0],pos[1],c="r")
                else:
                    ax.scatter(pos[0],pos[1],c="b")

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
        return Lambdas

def plotCellsCOM(ax, file):
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


        Lambdas = find_Lambdas(cells)
        for Lambda in Lambdas:
            cells1 = find_Lambda_cells(cells,Lambda=Lambda)
            xcom,ycom,zcom = centerBiofilm(cells1)
            ax.scatter(xcom,ycom,s=10, c="r")

        #Plots walls
        scale = 1
        wall_color = 'k'
        ax.plot([minx/scale, maxx/scale], [miny/scale, miny/scale], color=wall_color, alpha=0.6)
        ax.plot([minx/scale, maxx/scale], [maxy/scale, maxy/scale], color=wall_color, alpha=0.6)
        ax.plot([minx/scale, minx/scale], [miny/scale, maxy/scale], color=wall_color, alpha=0.6)
        ax.plot([maxx/scale, maxx/scale], [miny/scale, maxy/scale], color=wall_color, alpha=0.6)
        ax.plot([0,0],[miny/scale, maxy/scale],"--", color=wall_color, alpha=0.6)


        ax.set_xlim([-80, 80])
        ax.set_ylim([-5, 160])
        ax.axis('scaled')
        ax.axis('off')

def plotCells_channel(ax, file, width=40, MM=False):
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
        
        cells = remove_out_channel(cells,width)

        fcp.addAllCellsToPlot(cells, ax, ax_rng=maxx - minx, show_id=False, ec='w')


        #Plots walls
        if MM:
            scale = 1
            wall_color = 'k'
            ax.plot([-60, 60], [y_top/scale, y_top/scale],"--", color=wall_color, alpha=0.6)
            ax.plot([-60, 60], [y_bottom/scale, y_bottom/scale],"--", color=wall_color, alpha=0.6)

        else:
            scale = 0.3
            wall_color = 'k'
            ax.plot([y_bottom*1.5/scale, y_top*1.5/scale], [y_top, y_top], color=wall_color, alpha=0.6)
            ax.plot([y_bottom*1.5/scale, y_top*1.5/scale], [y_bottom, y_bottom], color=wall_color, alpha=0.6)
            #ax.plot([y_bottom*1.5/scale, y_bottom*1.5/scale], [y_top, y_bottom],"--", color=wall_color, alpha=0.6)
            #ax.plot([y_top*1.5/scale, y_top*1.5/scale], [y_top, y_bottom],"--", color=wall_color, alpha=0.6)
            


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
    if type(cells) == pd.DataFrame:
        x_centre = cells["pos_x"].mean()
        y_centre = cells["pos_y"].mean()
        z_centre = cells["pos_z"].mean()
    
    elif type(cells) == list:
        x_centre = np.mean([cell.pos_x for cell in cells])
        y_centre = np.mean([cell.pos_y for cell in cells])
        z_centre = np.mean([cell.pos_z for cell in cells])

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

    if type(cells) == pd.DataFrame:
        cells["pos_x"] =cells["pos_x"] - x_centre
        cells["pos_y"] = cells["pos_y"] - y_centre
        cells["pos_z"] = cells["pos_z"] - z_centre
    
    elif type(cells) == list:
        for cell in cells:
            cell.pos_x = cell.pos_x- x_centre
            cell.pos_y = cell.pos_y- y_centre
            cell.pos_z = cell.pos_z- z_centre
    return cells
    
def colonies_collided_COM(files):
    """
    This function performs 2 checks to see if two separate colonies have collided.

    We first analyise the COM of the colony with lower friction. Once a collision has happened,
    it will start to mvoe. Hence, we first look for this colony of COM related to the starting position and when
    it moves by more than 2 microns (max movement for a  free growing colony is 1.5), we consider that the 
    colonies have interacted.

    However, to avoid any large fluctuations that may happen due to division events, a second test is done. This
    second test may be looked at as more accurate than the first one but it is much more computationally costly. 
    This step involves looking at all cells from colony 1 and computing their distance to every cell in colony 2.
    The colonies will be interacting if ANY cell from 1 is within 2 microns of colony 2 (av length of a bacteria is
    around 4-5 micoons).

    
    """
    cells = ut.getCells(files[0])
    Lambdas = sorted(find_Lambdas(cells))
    xinitial,yinitial,zinitial = centerBiofilm(find_Lambda_cells(cells,Lambdas[0]))

    t = 0
    for file in files[1:]:
        t+=0.1
        cells = ut.getCells(file)
        cells1 = find_Lambda_cells(cells,Lambdas[0])
        x,y,z = centerBiofilm(cells1)
        d = np.sqrt((x-xinitial)**2 + (y-yinitial)**2)
        if  d > 1.2: # see if center of mass of smaller colony has moved (check 1)
            dat = pd.read_csv(file, sep='\t')
            cells2 = find_Lambda_cells(dat,Lambdas[1])
            ypos2 = cells2["pos_y"]
            xpos2 = cells2["pos_x"]
            for cell in cells1:
                d_point = np.sqrt((cell.pos_x-xpos2)**2 + (cell.pos_y-ypos2)**2)
                if np.sum(d_point < 5) != 0:
                    return t

def colonies_collided(files):
    """
    This function performs 2 checks to see if two separate colonies have collided.

    We first 

    However, to avoid any large fluctuations that may happen due to division events, a second test is done. This
    second test may be looked at as more accurate than the first one but it is much more computationally costly. 
    This step involves looking at all cells from colony 1 and computing their distance to every cell in colony 2.
    The colonies will be interacting if ANY cell from 1 is within 2 microns of colony 2 (av length of a bacteria is
    around 4-5 micoons).

    
    """
    growth_rate = 5 #growth rate in microns/h

    cdat = pd.read_csv(files[0], sep='\t')
    Lambdas = sorted(find_Lambdas(cdat))
    d0 = np.sqrt(cdat["pos_x"].diff()**2 +cdat["pos_y"].diff()**2)

    t_est = round(d0/growth_rate,1) # estimated collision time
    t = t_est[1] 

    """if t < 3: #impose a minimum of t_est > 3, prior there is too much variation
        t = 3
    print(t)"""

    collision_detection ="contour" #"position" or "perimeter"

    for file in files[int(t*10):]:
        t+=0.1

        cells = ut.getCells(file) #easier to loop through
        cells = centerCells(cells)

        cells1 = find_Lambda_cells(cells,Lambdas[0])
        dat = pd.read_csv(file, sep='\t') # easier to calculate position differences
        cells2 = find_Lambda_cells(dat,Lambdas[1])
        ypos2 = cells2["pos_y"]
        xpos2 = cells2["pos_x"]

        if collision_detection == "cells":
            for cell in cells1:
                d_point = np.sqrt((cell.pos_x-xpos2)**2 + (cell.pos_y-ypos2)**2)
                if np.sum(d_point < 2) > 1:
                    return t
                
        elif collision_detection =="contour": #when colonies interact, perimeter increases a lot
            contour = computeColonyContour(cells1)
            for position in contour:
                d_point = np.sqrt((position[0]-xpos2)**2 + (position[1]-ypos2)**2)
                if np.sum(d_point < 1.75) !=0:
                    #print(t)
                    return t
        else:
            print("please set the algorythm to either cells or contour")

def remove_out_channel(cells,width=60):
    if type(cells) == pd.DataFrame:
        cells = cells[abs(cells["pos_y"]) < (width/2)]
    elif type(cells) == list:
    
        cells1 = cells
        for cell in cells1:
            if abs(cell.pos_y) > (width/2):
                cells.remove(cell)
    return cells

def average_length(cells):
    if type(cells) == pd.DataFrame:
        l_av = cells["length"].mean()
    
    elif type(cells) == list:
        l_av = np.mean([cell.length for cell in cells])
    return l_av



def distance_from_origin(cells):
    if type(cells) == pd.DataFrame:
        cells = centerCells(cells)
        x = cells["pos_x"]
        y = cells["pos_y"]
        z = cells["pos_z"]
        distance = np.sqrt(x**2 + y**2 + z**2)
    
    elif type(cells) == list:
        cells = centerCells(cells)
        distance=[]
        for cell in cells:
            distance.append(np.sqrt(cell.pos_x**2 + cell.pos_y**2 + cell.pos_z**2))
    return distance

def distance_from_point(cells,point=[0,0,0]):
    if type(cells) == pd.DataFrame:
        cells = centerCells(cells)
        x = cells["pos_x"] - point[0]
        y = cells["pos_y"] - point[1]
        z = cells["pos_z"] - point[2]
        distance = np.sqrt(x**2 + y**2 + z**2)
    
    elif type(cells) == list:
        cells = centerCells(cells)
        distance=[]
        for cell in cells:
            distance.append(np.sqrt((cell.pos_x - point[0])**2 + (cell.pos_y - point[1])**2 + (cell.pos_z - point[2])**2))
    return distance

def calc_perimeter_area(cells):
    positions = computeColonyContour(cells)
    np.asarray(positions)

    polygon =  Polygon(positions)

    area = polygon.area
    perimeter = polygon.length
    """i = 0
    #plt.scatter(positions[:,0],positions[:,1]) Plot that shows it works

    for i in range(0,len(positions)-1):
        perimeter += np.sqrt((positions[i+1,0]-positions[i,0])**2 + (positions[i+1,1]-positions[i,1])**2)

    perimeter += np.sqrt((positions[0,0]-positions[-1,0])**2 + (positions[0,1]-positions[-1,1])**2) #final point"""
    return perimeter, area

def perimeter_area_time(data_dir):
    print(data_dir)
    file_pattern = os.path.join(data_dir, "biofilm_*.dat")
    files = sorted(glob.glob(file_pattern), key=lambda x: int(re.search(r'(\d+)', x).group(1)))
    time_steps = []

    ptot = []
    perimeter1 = []
    area1 = []
    perimeter2 = []
    area2 = []
    atot = []

    for file_path in files:
        match = re.search(r'(\d+)', os.path.basename(file_path))
        if not match:
            continue
        
        time_step = int(match.group(1))
        cells = ut.getCells(file_path)
        p,a = calc_perimeter_area(cells)
        ptot.append(p)
        atot.append(a)


        Lambdas = find_Lambdas(cells)
        cells1 =find_Lambda_cells(cells,Lambdas[0])

        time_steps.append(time_step * 0.1)
        p,a = calc_perimeter_area(cells1)
        perimeter1.append(p)
        area1.append(a)
    
        if len(Lambdas) != 1:
            cells2 =find_Lambda_cells(cells,Lambdas[1])
            p,a = calc_perimeter_area(cells2)
            perimeter2.append(p)
            area2.append(a)

    return time_steps, perimeter1, perimeter2, ptot, area1, area2,atot, Lambdas

def surface_fraction(data_dir):
    file_pattern = os.path.join(data_dir, "biofilm_*.dat")
    files = sorted(glob.glob(file_pattern), key=lambda x: int(re.search(r'(\d+)', x).group(1)))
    time_steps = []
    
    t_collision = colonies_collided(files) #finds time when colonies collide

    frac1 = []
    frac2 = []
    
    i = 0
    for file_path in files:
        match = re.search(r'(\d+)', os.path.basename(file_path))
        if not match:
            continue
        time_step = int(match.group(1))

        if time_step*0.1 < t_collision:
            continue
        time_steps.append(i * 0.1) #only appends after colonies have collided
        
        cells = ut.getCells(file_path)
        dat = pd.read_csv(file_path, sep='\t')
        Lambdas = sorted(find_Lambdas(cells))

        contour = np.asarray(computeColonyContour(cells))

        Lambda_point =[]

        for point in contour:
            d = np.sqrt((dat["pos_x"]-point[0])**2 + (dat["pos_y"]-point[1])**2)
            idxmin = d.idxmin()
            Lambda_point.append(dat.loc[idxmin, 'Lambda'])

                
        count1 = Lambda_point.count(Lambdas[0])
        count2 = Lambda_point.count(Lambdas[1])

        frac1.append(count1/(count1+count2))
        frac2.append(count2/(count1+count2))

        """if ax != None:
            cx = contour[:,0]
            cy = contour[:,1]
            for Lambda in Lambdas:
                mask = np.where(np.asarray(Lambda_point)==Lambda)
                ax.scatter(cx[mask],cy[mask],color=colour_Lambda(Lambda))"""
        i+=1
    return time_steps, frac1, frac2, Lambdas

def surface_fraction_high(data_dir):
    """
    Finds the fraction of the the contour that is part of the interface of the colony for the higher Lambda.

    It calculates the total contour of the colonies combined and then checks the fraction of points of the higher Lambda 
    contour that belongs to the total contour. This gives us the points of the HIGHER LAMBDA contour that are exposed
    to the exterior of the colony, and the points that are part of the interface.

    From this we calculate the total fraction of the points to find the fraction of the contour that is part of the
    interface.

    returns:
        frac_interface: 
        frac_external:

    
    """
    file_pattern = os.path.join(data_dir, "biofilm_*.dat")
    files = sorted(glob.glob(file_pattern), key=lambda x: int(re.search(r'(\d+)', x).group(1)))
    time_steps = []
    
    t_collision = colonies_collided(files) #finds time when colonies collide

    frac_interface = []
    frac_external = []
    
    i = 0
    for file_path in files:
        match = re.search(r'(\d+)', os.path.basename(file_path))
        if not match:
            continue
        time_step = int(match.group(1))

        if time_step*0.1 < t_collision:
            continue
        time_steps.append(i * 0.1) #only appends after colonies have collided
        
        cells = ut.getCells(file_path)
        Lambdas = sorted(find_Lambdas(cells))
        cellshigh = find_Lambda_cells(cells,Lambda=Lambdas[1])
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

                
        count1 = Lambda_point.count(Lambdas[0])
        count2 = Lambda_point.count(Lambdas[1])

        frac_interface.append(count1/(count1+count2))
        frac_external.append(count2/(count1+count2))
        i+=1
    return time_steps, frac_interface, frac_external, Lambdas

def s0_interfaceXsurface(file_path):
    cells = ut.getCells(file_path)
    Lambdas = sorted(find_Lambdas(cells))
    cellshigh = find_Lambda_cells(cells,Lambda=Lambdas[1])
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

                
    count1 = Lambda_point.count(Lambdas[0])
    count2 = Lambda_point.count(Lambdas[1])

    frac_interface = count1/(count1+count2)
    frac_external = count2/(count1+count2)
    return frac_interface

def is_in_circle(x,y,cx,cy,r):
    return (x-cx)**2 +(y-cy)**2 <= r**2

def find_interface(ax,cells):
    Lambdas = find_Lambdas(cells)
    cells1 = find_Lambda_cells(cells,Lambdas[0])
    cells2 = find_Lambda_cells(cells,Lambdas[1])

    c1_pos = []
    for cell in cells1:
        c1_pos.append([cell.pos_x,cell.pos_y])
    c1_pos = np.asarray(c1_pos)

    c2_pos = []
    for cell in cells2:
        c2_pos.append([cell.pos_x,cell.pos_y])
    c2_pos = np.asarray(c2_pos)

    contour1 = computeColonyContour(cells1)
    contour2 = computeColonyContour(cells2)

    c1_interface = []
    c1_rest = []

    radius = 10

    for point in contour1:
        if is_in_circle(c1_pos[:,0],c1_pos[:,1],point[0],point[1],radius).any() and is_in_circle(c2_pos[:,0],c2_pos[:,1],point[0],point[1],radius).any():
            c1_interface.append(point)
        else:
            c1_rest.append(point)
    c1_interface = np.asarray(c1_interface)
    c1_rest = np.asarray(c1_rest)

    c2_interface = []
    c2_rest = []

    for point in contour2:
        if is_in_circle(c1_pos[:,0],c1_pos[:,1],point[0],point[1],radius).any() and is_in_circle(c2_pos[:,0],c2_pos[:,1],point[0],point[1],radius).any():
            c2_interface.append(point)
        else:
            c2_rest.append(point)
    c2_interface = np.asarray(c2_interface)
    c2_rest = np.asarray(c2_rest)

    c1_interface = remove_spatial_outliers(c1_interface)
    c1_rest = remove_spatial_outliers(c1_rest)
    c2_interface = remove_spatial_outliers(c2_interface)
    c2_rest = remove_spatial_outliers(c2_rest)


    ax.plot(c1_interface[:,0],c1_interface[:,1],"r.")
    ax.plot(c1_rest[:,0],c1_rest[:,1],"b.")

    ax.plot(c2_interface[:,0],c2_interface[:,1],"g.")
    ax.plot(c2_rest[:,0],c2_rest[:,1],"k.")



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
    eig_val1,eig_val2 = np.linalg.eigvals(Tensor)
    Rg=np.sqrt(eig_val1+eig_val2) #Rg = sqrt(lambda1 + lambda2)

    return Rg,[eig_val1,eig_val2]

def Gyration_tensor_2D(cells):
    Gyr_tensor = np.zeros((2,2))
    cells = centerCells(cells)
    x,y = position_cells(cells)
    Gyr_tensor[0,0]= np.mean(x*x)
    Gyr_tensor[1,1]= np.mean(y*y)
    Gyr_tensor[0,1]= np.mean(x*y)
    Gyr_tensor[1,0]= np.mean(y*x)
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

def calc_aspect_ratio(data_dir):
    file_pattern = os.path.join(data_dir, "biofilm_*.dat")
    files = sorted(glob.glob(file_pattern), key=lambda x: int(re.search(r'(\d+)', x).group(1)))
    
    asp_ratio = []
    time_steps = []


    for file_path in files:
        match = re.search(r'(\d+)', os.path.basename(file_path))
        if not match:
            continue

        time_step = int(match.group(1))
        cells = ut.getCells(file_path)
        time_steps.append(time_step * 0.1)

        rg,eigenv = Gyration_values(cells)
        
        eigenv=sorted(eigenv)
        if (eigenv[0] < 0.1) or (eigenv[1] < 0.1): #one cell
            asp_ratio.append(1) 

        else:
            asp_ratio.append(eigenv[1]/eigenv[0])
            print(eigenv[1]/eigenv[0])
        print(eigenv)
            
    return time_steps, asp_ratio

def get_orientation_cells(cells):
    if type(cells) == pd.DataFrame:
        theta = np.atan2(cells["ori_y"],cells["ori_x"])    
    elif type(cells) == list:
        theta = []
        for cell in cells:
            theta.append(np.atan2(cell.ori_y,cell.ori_x))
    return np.asarray(theta)*180/np.pi



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

def RgLambda_time_est(Rg, time_step=0.1,total=False):
    Rg = np.asarray(Rg)
    Rg[Rg<1] = 1
    linear_rg = np.log2(Rg)
    t = np.arange(0, (len(Rg) - 0.5)* time_step, time_step) #0.5 INCLUDED TO AVOID ANY SMALL ERRORS IN FLOAT POINTS
    if total: #total starts higher so becomes circular later
        t = np.arange(4, (len(Rg) - 0.5)* time_step, time_step)
        linear_rg = linear_rg[40:]
    popt, pcov= curve_fit(Rg_Lambda_growth, t, linear_rg, p0=(2,2.5)) # tau_rg = 2*tdouble - see notes
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
            if (eigenv[0]== 0) and (eigenv[1] == 0): #one cell
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
    elif Lambda == 5.001:
        return "#4c564d"
    elif Lambda == 2.0:
        return "#ff0000"
    elif Lambda < 1.0:
        return (0.1/Lambda, 0, 0, 1)
    elif Lambda > 1.0:
        return (0,(Lambda-1)/(10-1), 0, 1)

def colour_ratio(ratio):
    if (abs(ratio)- 1.0)<0.1:
        return "#00ffff"
        """elif ratio == 2.0:
            return "#ff0000"""
    elif ratio > 1.0:
        return (0,(ratio-1)/(10-1), 0, 1)


def counts(data_dir,channels=False,width=60):
    file_pattern = os.path.join(data_dir, "biofilm_*.dat")
    files = sorted(glob.glob(file_pattern), key=lambda x: int(re.search(r'(\d+)', x).group(1)))
    
    time_steps = []
    Lambda1_counts = []
    Lambda2_counts = []

    Lambdas = sorted(find_Lambdas(pd.read_csv(files[0], sep="\t")))

    for file_path in files:
        match = re.search(r'(\d+)', os.path.basename(file_path))
        if not match:
            continue
        time_step = int(match.group(1))
        df = pd.read_csv(file_path, sep="\t")
        if channels:
            df = remove_out_channel(df,width)

        time_steps.append(time_step * 0.1)

        Lambda1_count = find_Lambda_cells(df, Lambdas[0]).shape[0]
        Lambda1_counts.append(Lambda1_count)
        if len(Lambdas) > 1:
            Lambda2_count = find_Lambda_cells(df, Lambdas[1]).shape[0]
            Lambda2_counts.append(Lambda2_count)

    return time_steps, Lambda1_counts, Lambda2_counts, Lambdas

def doubling_linear_growth(t, t_doub, N0):
        return N0 + t/t_doub
def doubling_exp_growth(t, t_doub):
        return  2**(t/t_doub)
def linear_exp(t,A,x,m):
    return m+A*(t)**x
def tanhx_simple(t,t_star):
    return np.tanh(t/t_star)
def tanhx(t,t_star,A):
    return A*np.tanh(t/t_star)
def tanhx_fraction(rho,w,d_star):
    return d_star*(0.5+ np.tanh(rho/w))

def find_tanh_1param(xdata,ydata):
    """
        Parameters:
            counts: list of int
                cell counts over time
            time_step: int
                time step of the simulation

        Returns:
            doubling time parameter popt[0] and its error np.sqrt(np.diag(pcov))[0]
    """
    popt, pcov= curve_fit(tanhx_simple, xdata, ydata,maxfev = 1000000,p0=(1))

    return popt,np.sqrt(np.diag(pcov))

def find_tanh_2param(xdata,ydata,t_guess=None):
    """
        Parameters:
            counts: list of int
                cell counts over time
            time_step: int
                time step of the simulation

        Returns:
            doubling time parameter popt[0] and its error np.sqrt(np.diag(pcov))[0]
    """
    if t_guess is None:
        popt, pcov= curve_fit(tanhx, xdata, ydata,maxfev = 1000000,p0=(1,1))
    else:
        popt, pcov= curve_fit(tanhx, xdata, ydata,maxfev = 1000000,p0=(t_guess,1),bounds=([0, 0], [np.inf, 5]))

    return popt,np.sqrt(np.diag(pcov))

def find_tanh_1param(xdata,ydata):
    """
        Parameters:
            counts: list of int
                cell counts over time
            time_step: int
                time step of the simulation

        Returns:
            doubling time parameter popt[0] and its error np.sqrt(np.diag(pcov))[0]
    """
    popt, pcov= curve_fit(tanhx_simple, xdata, ydata,maxfev = 1000000,p0=(1))

    return popt,np.sqrt(np.diag(pcov))

def find_tanh_fraction_dist(xdata,ydata):
    """
        Parameters:
            counts: list of int
                cell counts over time
            time_step: int
                time step of the simulation

        Returns:
            doubling time parameter popt[0] and its error np.sqrt(np.diag(pcov))[0]
    """
    popt, pcov= curve_fit(tanhx_fraction, xdata, ydata,maxfev = 1000000,p0=(1,1))

    return popt,np.sqrt(np.diag(pcov))


def find_scaling_law(xdata,ydata):
    """
        Parameters:
            counts: list of int
                cell counts over time
            time_step: int
                time step of the simulation

        Returns:
            doubling time parameter popt[0] and its error np.sqrt(np.diag(pcov))[0]
    """
    popt, pcov= curve_fit(linear_exp, xdata, ydata,maxfev = 1000000,p0=(1,0.5,0))

    return popt,np.sqrt(np.diag(pcov))

def average_std(dataset):
    if type(dataset) == pd.DataFrame or type(dataset) == pd.Series:
        av = dataset.mean()
        error = dataset.std()/np.sqrt(dataset.shape[0])
        return av,error

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
    popt, pcov= curve_fit(doubling_exp_growth, t, counts, p0=(0,1.8))

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



#def calc_microdomain_size():


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

def calc_fraction(n_L1,n_L2):
    frac1 = n_L1/(n_L1+n_L2)
    frac2 = n_L2/(n_L1+n_L2)
    return frac1, frac2


def find_fraction_distance(filepath,n_points=10,com="total"):
    df = pd.read_csv(filepath, sep="\t")
    df = centerCells(df)
    Lambdas = sorted(find_Lambdas(df))
    if com=="total":
        d_cells = np.asarray(distance_from_origin(df))
    if com=="higher":
        centre = centerBiofilm(find_Lambda_cells(df,Lambdas[-1]))
        d_cells = np.asarray(distance_from_point(df,centre))
    
    dmax = max(d_cells)

    d_cells_norm = d_cells/dmax

    counts = np.zeros((2,n_points))

    ds = 1/n_points
    distance = np.linspace(0, 1, n_points+1)
    frac_higher = []
    frac_lower = []
 
    for i in range(0,n_points):
        bin_mask = (d_cells_norm >= distance[i]) & (d_cells_norm < distance[i+1])

        counts[0, i] = np.count_nonzero(bin_mask & (df["Lambda"] == Lambdas[0]))
        counts[1, i] = np.count_nonzero(bin_mask & (df["Lambda"] == Lambdas[1]))
        if (counts[0, i] != 0) or (counts[1, i]) != 0:
            frac_higher.append(counts[1,i]/(np.add(counts[1,i],counts[0,i])))
            frac_lower.append(counts[0,i]/(np.add(counts[1,i],counts[0,i])))
        else: 
            frac_higher.append(0)
            frac_lower.append(0)
    
    distance +=ds/2
    ratio = round(Lambdas[1]/Lambdas[0],1)

    return distance[:-1],np.asarray(frac_higher),np.asarray(frac_lower),ratio

def find_invasion_time(data_dir):
    time_steps, Lambda1_counts, Lambda2_counts, Lambdas = counts(data_dir)
    Lambdas = sorted(Lambdas)
    frac1,frac2 = calc_fraction(np.asarray(Lambda1_counts), np.asarray(Lambda2_counts))
    
    t_fill = 3.2
    frac_target = 0.25
    if sum(frac1 <= frac_target) != 0:
        index = next((i for i, val in enumerate(frac1[int(t_fill*10):]) if val < frac_target ), -1)
        t_invasion = round(index*0.1,1)
        return t_invasion,round(Lambdas[1]/Lambdas[0],1)
    else:
        print(f"{round(Lambdas[1]/Lambdas[0],1)} no invasion tf = {time_steps[-1]} h")
    



from scipy.ndimage import minimum_filter

def locate_nematic_defects_vector(x, y, ux, uy, r_cut=8):
    """
    Finds defect positions and charges from particle vector data.
    x, y: 1D arrays of bacteria positions
    ux, uy: 1D arrays of bacteria orientation unit vector components
    """
    # 1. Define Spatial Grid boundaries
    grid_size = int((max(x)-min(x))/4)
    x_grid = np.linspace(x.min(), x.max(), grid_size)
    y_grid = np.linspace(y.min(), y.max(), grid_size)
    X, Y = np.meshgrid(x_grid, y_grid)
    
    Qxx = np.zeros_like(X)
    Qxy = np.zeros_like(Y)

    Density = np.zeros_like(X)
    
    # 2. Compute Coarse-Grained Q-Tensor Field directly from vector components
    for i in range(grid_size):
        for j in range(grid_size):
            # Find particles inside local radius
            dist = np.sqrt((x - X[i,j])**2 + (y - Y[i,j])**2)
            mask = dist < r_cut
            num_particles = np.sum(mask)
            Density[i,j] = num_particles # Track how many bacteria are here
    
            if np.sum(mask) > 0:
                Qxx[i,j] = np.mean(ux[mask]**2 - 0.5)
                Qxy[i,j] = np.mean(ux[mask] * uy[mask])

    # 3. Compute Scalar Order Parameter S and Director Angle Phi
    S = 2.0 * np.sqrt(Qxx**2 + Qxy**2)
    Phi = 0.5 * np.arctan2(Qxy, Qxx) # Needed purely for loop integration
    
    # 4. Find Local Minima of S (Defect Cores)
    valid_colony_zone = (Density >= 2)
    local_min = (S == minimum_filter(S, footprint=np.ones((3,3))))
    
    # Combine masks: must be a local minimum, below threshold, AND inside the cluster
    defect_mask = local_min & (S < 0.25) & valid_colony_zone
    defect_indices = np.argwhere(defect_mask)
    defect_indices = np.argwhere(defect_mask)
    
    defects = []
    
    # 5. Calculate Winding Number (k) for each candidate core
    for idx in defect_indices:
        r, c = idx[0], idx[1]
        if r == 0 or r == grid_size-1 or c == 0 or c == grid_size-1:
            continue # Skip edge boundaries
            
        # Extract 4 corners of a small loop around the core
        loop_phi = [Phi[r, c-1], Phi[r+1, c], Phi[r, c+1], Phi[r-1, c]]
        
        # Calculate differences along the loop, wrapping to [-pi/2, pi/2]
        dphi = 0
        for m in range(4):
            diff = loop_phi[(m+1)%4] - loop_phi[m]
            diff = (diff + np.pi/2) % np.pi - np.pi/2 # Nematic periodic wrapping
            dphi += diff
            
        charge = dphi / (2 * np.pi)
        
        # Classify as +1/2 or -1/2 defect
        if np.isclose(abs(charge), 0.5, atol=0.1):
            defects.append({
                'pos': (X[r,c], Y[r,c]),
                'charge': np.sign(charge) * 0.5
            })
            
    return defects, S, Phi



test_exit = "C:\\Users\\lucca\\Desktop\\GeneratedOutput\\SimOutput\\data\\InfiniteChannel\\\MM\\width2_15\\Lambda0_1\\repeat20"
#print(find_invasion_time(test_exit))