import numpy as np
import scipy as sp
from . import Simulator
from scipy.spatial import cKDTree
from SimPEG.Utils import mkvc
import geosoft.gxpy.grid as gxgrd
import geosoft.gxpy.gx as gx
import matplotlib.pyplot as plt
import gdal
import osr
import ogr
import os
from shapely.geometry import mapping, LineString
import fiona
from fiona.crs import from_epsg


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
        fig, im, cbar = Simulator.plotData2D(
            xLoc, yLoc, data.values, marker=False, fig=fig, ax=axs,
            colorbar=True
        )
        axs.grid(True)
        cbar.set_label('TMI (nT)')
        plt.show()
        fig.savefig('./images/SearchQuestII.png', bbox_inches='tight')

    return data


def loadGeoTiffFile(fileName, plotIt=True):
    """
        Load a data matrix in Geosoft *.grd format and return
        a dictionary with gridded data
    """
    data = dataGrid()

    rasterObject = gdal.Open(fileName)
    band = rasterObject.GetRasterBand(1)
    data.values = band.ReadAsArray()
    data.nx, data.ny = data.values.shape[1], data.values.shape[0]
    data.x0, y0 = rasterObject.GetGeoTransform()[0], rasterObject.GetGeoTransform()[3]
    data.dx, data.dy = np.round(rasterObject.GetGeoTransform()[1]), np.round(np.abs(rasterObject.GetGeoTransform()[5]))
    data.y0 = y0 - data.ny*data.dy
    data.limits = np.r_[data.x0, data.x0+data.nx*data.dx, data.y0, y0]

    if plotIt:
        xLoc = np.asarray(range(data.nx))*data.dx+data.x0
        yLoc = np.asarray(range(data.ny))*data.dy+data.y0

        fig, axs = plt.figure(figsize=(8, 8)), plt.subplot()
        fig, im, cbar = Simulator.plotData2D(
            xLoc, yLoc, data.values, marker=False, fig=fig, ax=axs,
            colorbar=True
        )
        axs.grid(True)
        cbar.set_label('TMI (nT)')
        plt.show()
        fig.savefig('./images/GeoTIFFSynthetic.png', bbox_inches='tight')

    return data


class dataGrid(object):
    """
        Grid data object
    """

    x0, y0 = 0., 0.
    nx, ny = 1, 1
    dx, dy = 1., 1.
    values = None

    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)

        return

    @property
    def hx(self):

        if getattr(self, '_hx', None) is None:
            self._hx = np.asarray(range(self.nx)) * self.dx + self.x0

        return self._hx

    @property
    def hy(self):

        if getattr(self, '_hy', None) is None:
            self._hy = np.asarray(range(self.ny)) * self.dy + self.y0

        return self._hy


def arrayToRaster(array, fileName, EPSGCode, xMin, xMax, yMin, yMax, numBands, dataType='image'):
    """
        Source:

            Cameron Cooke: http://cgcooke.github.io/GDAL/

    """

    xPixels = array.shape[1]  # number of pixels in x
    yPixels = array.shape[0]  # number of pixels in y
    pixelXSize = (xMax-xMin)/(xPixels-1)  # size of the pixel in X direction
    pixelYSize = -(yMax-yMin)/(yPixels-1)  # size of the pixel in Y direction

    driver = gdal.GetDriverByName('GTiff')

    # Chose type
    if dataType=='image':
        encodeType = gdal.GDT_Byte
    else:
        encodeType = gdal.GDT_Float32

    dataset = driver.Create(
        fileName, xPixels, yPixels, numBands,
        encodeType,
    )

 #options=['PHOTOMETRIC=RGB']
    dataset.SetGeoTransform((xMin, pixelXSize, 0, yMax, 0, pixelYSize))

    datasetSRS = osr.SpatialReference()
    datasetSRS.ImportFromEPSG(EPSGCode)
    dataset.SetProjection(datasetSRS.ExportToWkt())

    if numBands == 1:
        dataset.GetRasterBand(1).WriteArray(array)
    else:
        for i in range(0, numBands):
            dataset.GetRasterBand(i+1).WriteArray(array[:, :, i])

    dataset.FlushCache()  # Write to disk.
    print('Image saved as: ' + fileName + ' Click box again to continue...')


def exportShapefile(polylines, attributes, EPSGCode=26909, fileName='MyShape', label='AverageDepth', attType='int'):
    """
        Function to export polylines to shape file with attribute

    """
    crs = from_epsg(EPSGCode)

    # Define a polygon feature geometry with one attribute
    schema = {
        'geometry': 'LineString',
        'properties': {label: attType},

    }

    with fiona.open(fileName + '.shp', 'w', 'ESRI Shapefile', schema, crs=crs) as c:
        ## If there are multiple geometries, put the "for" loop here
        for poly, att in zip(polylines, attributes):

            if np.all([np.any(att), len(poly)>1]):
                pline = LineString(list(tuple(map(tuple, poly))))

                res = {}
                res['properties'] = {}
                res['properties'][label] = np.mean(att)
                # geometry of of the original polygon shapefile
                res['geometry'] = mapping(pline)
                c.write(res)
