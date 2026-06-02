import numpy as np
import matplotlib.pyplot as plt


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

    return x_centre, y_centre

def centerCells(cells):
    """
        Parameters:
            cells: list of ChainingRodShapedBacterium
                cells to be centred

        Returns:
            cells with positions shifted to be centred around (0,0)
    """
    x_centre, y_centre = centerBiofilm(cells)

    for cell in cells:
        cell.pos_x -= x_centre
        cell.pos_y -= y_centre
    return cells

def radiusGyration(cells):
    """
        Parameters:
            cells: list of ChainingRodShapedBacterium
                cells to be analysed

        Returns:
            radius of gyration of the total cell population
    """
    x_centre, y_centre = centerBiofilm(cells)

    return np.sqrt(np.mean([(cell.pos_x - x_centre)**2 + (cell.pos_y - y_centre)**2 for cell in cells]))


def RgZeta(cells,zeta):
    """
        Parameters:
            cells: list of ChainingRodShapedBacterium
                cells to be analysed
            zeta: float
                non-dimensionalised drag coefficient

        Returns:
            radius of gyration of the cell population for a given zeta
    """
    x_cell,y_cell = [],[]

    for cell in cells:
            if cell.non_dimzeta == zeta:
                x_cell.append(cell.pos_x)
                y_cell.append(cell.pos_y)

    x_cell = np.array(x_cell)
    y_cell = np.array(y_cell)

    x_centre = np.mean(x_cell)
    y_centre = np.mean(y_cell)

    return np.sqrt(np.mean((x_cell - x_centre)**2 + (y_cell - y_centre)**2 ))
