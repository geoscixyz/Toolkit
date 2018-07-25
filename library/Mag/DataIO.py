import numpy as np
import scipy as sp
from . import Simulator
from scipy.spatial import cKDTree
from SimPEG.Utils import mkvc
import geosoft.gxpy.grid as gxgrd
import geosoft.gxpy.gx as gx
import matplotlib.pyplot as plt


def loadGRDFile(fileName, plotIt=True):
    """
        Load a data matrix in Geosoft *.grd format and return
        a dictionary with gridded data
    """
    gxc = gx.GXpy()
    data = dataGrid()

    with gxgrd.Grid(fileName) as grid:

        lim = grid.extent_2d()
        data.limits = np.r_[lim[0], lim[2], lim[1], lim[3]]
        # coordinate_system = grid.coordinate_system
        data.values = grid.xyzv()[:, :, 3]
        data.nx, data.ny = grid.nx, grid.ny
        data.dx, data.dy = grid.dx, grid.dy
        data.x0, data.y0 = grid.x0, grid.y0

    if plotIt:
        xLoc = np.asarray(range(data.nx))*data.dx+data.x0
        yLoc = np.asarray(range(data.ny))*data.dy+data.y0

        fig, axs = plt.figure(figsize=(8, 8)), plt.subplot()
        Simulator.plotData2D(
            xLoc, yLoc, data.values, marker=False, fig=fig, ax=axs
        )
        axs.grid(True)

        plt.show()
        fig.savefig('./images/SearchQuestII.png', bbox_inches='tight')

    return data


class dataGrid(object):
    """
        Grid data object
    """

    x0, y0 = 0., 0.
    nx, ny = 1, 1
    dx, dy = 1., 1.
    limits = np.r_[0, 1, 0, 1]
    values = None

    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)

        return


