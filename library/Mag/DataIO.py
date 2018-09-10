import numpy as np
import scipy as sp
from . import (Simulator, MathUtils)
from scipy.spatial import cKDTree
from SimPEG.Utils import mkvc
import matplotlib.pyplot as plt
import gdal
import osr
import ogr
import os
from shapely.geometry import mapping, LineString
import fiona
from fiona.crs import from_epsg


class dataGrid(object):
    """
        Grid data object
    """

    x0, y0 = 0., 0.
    nx, ny = 1, 1
    dx, dy = 1., 1.
    values = None
    valuesFilled = None
    valuesFilledUC = None
    inc = 90.
    dec = 90.
    indVal = None
    indNan = None

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

    @property
    def npadx(self):

        if getattr(self, '_npadx', None) is None:
            self._npadx = int(np.floor(self.values.shape[1]))

        return self._npadx

    @property
    def gridCC(self):

        if getattr(self, '_gridCC', None) is None:

            X, Y = np.meshgrid(self.hx, self.hy)

            self._gridCC = np.c_[mkvc(X), mkvc(Y)]

        return self._gridCC

    @property
    def npady(self):

        if getattr(self, '_npady', None) is None:
            self._npady = int(np.floor(self.values.shape[0]))

        return self._npady

    @property
    def gridPadded(self):

        if getattr(self, '_gridPadded', None) is None:
            if getattr(self, '_valuesFilledUC', None) is not None:
                # if ~np.array_equal(
                #             self.valuesFilled[~np.isnan(self.valuesFilled)],
                #             self.values[~np.isnan(self.indNan)]
                #         ):

                #             self._valuesFilled = None
                grid = self.valuesFilledUC
            else:

                if np.any(np.isnan(self.values)):
                    self.indNan = np.isnan(mkvc(self.values))
                    grid = self.valuesFilled
                else:
                    grid = self.values

            # Add paddings
            dpad = np.c_[
                np.fliplr(grid[:, 0:self.npadx]),
                grid,
                np.fliplr(grid[:, -self.npadx:])
                ]

            dpad = np.r_[
                np.flipud(dpad[0:self.npady, :]),
                dpad,
                np.flipud(dpad[-self.npady:, :])
                ]

            # Tapper the paddings
            rampx = -np.cos(np.pi*np.asarray(range(self.npadx))/self.npadx)
            rampx = np.r_[rampx, np.ones(grid.shape[1]), -rampx]/2. + 0.5

            rampy = -np.cos(np.pi*np.asarray(range(self.npady))/self.npady)
            rampy = np.r_[rampy, np.ones(grid.shape[0]), -rampy]/2. + 0.5
            tapperx, tappery = np.meshgrid(rampx, rampy)

            self._gridPadded = tapperx*tappery*dpad

        return self._gridPadded

    @property
    def valuesFilled(self):

        if getattr(self, '_valuesFilled', None) is None:
            values = mkvc(self.values)
            indVal = np.where(~self.indNan)[0]

            tree = cKDTree(self.gridCC[indVal, :])
            dists, indexes = tree.query(self.gridCC[self.indNan, :])

            uInd = np.unique(indVal[indexes])

            xyz = self.gridCC[uInd, :]

            _, values[self.indNan] = MathUtils.minCurvatureInterp(
                            xyz, values[uInd], xyzOut=self.gridCC[self.indNan, :])

            self._valuesFilled = values.reshape(self.values.shape, order='F')

        return self._valuesFilled

    @property
    def valuesFilledUC(self):
        if getattr(self, '_valuesFilledUC', None) is None:
            self._valuesFilledUC = self.upwardContinuation()

        return self._valuesFilledUC

    @property
    def Kx(self):

        if getattr(self, '_Kx', None) is None:

            dx = (self.dx)/(self.gridPadded.shape[1]-1)
            dy = (self.dy)/(self.gridPadded.shape[0]-1)

            kx = np.fft.fftfreq(self.gridPadded.shape[1], dx)
            ky = np.fft.fftfreq(self.gridPadded.shape[0], dy)

            Ky, Kx = np.meshgrid(ky, kx)
            self._Ky, self._Kx = Ky.T, Kx.T

        return self._Kx

    @property
    def Ky(self):

        if getattr(self, '_Ky', None) is None:

            dx = (self.dx)/(self.gridPadded.shape[1]-1)
            dy = (self.dy)/(self.gridPadded.shape[0]-1)

            kx = np.fft.fftfreq(self.gridPadded.shape[1], dx)
            ky = np.fft.fftfreq(self.gridPadded.shape[0], dy)

            Ky, Kx = np.meshgrid(ky, kx)
            self._Ky, self._Kx = Ky.T, Kx.T

        return self._Ky

    @property
    def gridFFT(self):

        if getattr(self, '_gridFFT', None) is None:

            self._gridFFT = np.fft.fft2(self.gridPadded)

        return self._gridFFT

    @property
    def derivativeX(self):

        if getattr(self, '_derivativeX', None) is None:

            FHxD = (self.Kx*1j)*self.gridFFT
            fhxd_pad = np.fft.ifft2(FHxD)
            derivX = np.real(
                fhxd_pad[self.npady:-self.npady, self.npadx:-self.npadx])

            if self.indNan is not None:
                derivX = mkvc(derivX)

                derivX[self.indNan] = np.nan
                derivX = derivX.reshape(self.values.shape, order='F')

            self._derivativeX = derivX

        return self._derivativeX

    @property
    def derivativeY(self):

        if getattr(self, '_derivativeY', None) is None:

            FHyD = (self.Ky*1j)*self.gridFFT

            fhyd_pad = np.fft.ifft2(FHyD)
            derivY = np.real(
                fhyd_pad[self.npady:-self.npady, self.npadx:-self.npadx])

            if self.indNan is not None:
                derivY = mkvc(derivY)

                derivY[self.indNan] = np.nan
                derivY = derivY.reshape(self.values.shape, order='F')

            self._derivativeY = derivY

        return self._derivativeY

    @property
    def firstVertical(self):

        if getattr(self, '_firstVertical', None) is None:

            FHzD = self.gridFFT*np.sqrt(self.Kx**2. + self.Ky**2.)
            fhzd_pad = np.fft.ifft2(FHzD)
            firstVD = np.real(
                fhzd_pad[self.npady:-self.npady, self.npadx:-self.npadx])

            if self.indNan is not None:
                firstVD = mkvc(firstVD)

                firstVD[self.indNan] = np.nan
                firstVD = firstVD.reshape(self.values.shape, order='F')

            self._firstVertical = firstVD
        return self._firstVertical

    @property
    def totalHorizontal(self):

        if getattr(self, '_totalHorizontal', None) is None:

            self._totalHorizontal = np.sqrt(
                self.derivativeX**2. + self.derivativeY**2. + 1e-8
            )

        return self._totalHorizontal

    @property
    def tiltAngle(self):

        if getattr(self, '_tiltAngle', None) is None:

            self._tiltAngle = np.arctan2(
                self.firstVertical, self.totalHorizontal
            )

        return self._tiltAngle

    @property
    def analyticSignal(self):

        if getattr(self, '_analyticSignal', None) is None:

            self._analyticSignal = np.sqrt(
                self.derivativeX**2. +
                self.derivativeY**2. +
                self.firstVertical**2. + 1e-8
            )

        return self._analyticSignal

    @property
    def RTP(self):

        if getattr(self, '_RTP', None) is None:

            h0_xyz = MathUtils.dipazm_2_xyz(self.inc, self.dec)
            Frtp = self.gridFFT / (
                (h0_xyz[2] + 1j*(self.Kx*h0_xyz[0] + self.Ky*h0_xyz[1]))**2.
            )
            rtp_pad = np.fft.ifft2(Frtp)
            rtp = np.real(
                rtp_pad[self.npady:-self.npady, self.npadx:-self.npadx])

            if self.indNan is not None:
                rtp = mkvc(rtp)

                rtp[self.indNan] = np.nan
                rtp = rtp.reshape(self.values.shape, order='F')

            self._RTP = rtp

        return self._RTP

    def upwardContinuation(self, z=0):
        """
            Function to calculate upward continued data
        """

        upFact = -(
            np.sqrt(self.Kx**2. + self.Ky**2.) *
            z /
            np.sqrt(self.dx**2. + self.dy**2.)
        )

        FzUpw = self.gridFFT*np.exp(upFact)
        zUpw_pad = np.fft.ifft2(FzUpw)
        zUpw = np.real(
            zUpw_pad[self.npady:-self.npady, self.npadx:-self.npadx])

        self._valuesFilledUC = zUpw.copy()

        if self.indNan is not None:
            zUpw = mkvc(zUpw)

            zUpw[self.indNan] = np.nan
            zUpw = zUpw.reshape(self.values.shape, order='F')

        return zUpw


def loadGRDFile(fileName, plotIt=True):
    """
        Load a data matrix in Geosoft *.grd format and return
        a dictionary with gridded data
    """

    try:
        import geosoft.gxpy.grid as gxgrd
        import geosoft.gxpy.gx as gx
    except ImportError:
        print("geosoft module not installed. loadGRDFile ")
        return

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
    data.nx = data.values.shape[1]
    data.ny = data.values.shape[0]
    data.x0 = rasterObject.GetGeoTransform()[0]
    y0 = rasterObject.GetGeoTransform()[3]
    data.dx = np.round(rasterObject.GetGeoTransform()[1])
    data.dy = np.round(np.abs(rasterObject.GetGeoTransform()[5]))
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


def arrayToRaster(
    array, fileName, EPSGCode,
    xMin, xMax, yMin, yMax,
    numBands, dataType='image'
):
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
    if dataType == 'image':
        encodeType = gdal.GDT_Byte
    else:
        encodeType = gdal.GDT_Float32

    dataset = driver.Create(
        fileName, xPixels, yPixels, numBands,
        encodeType,
    )

    print(xMin, pixelXSize, 0, yMax, 0, pixelYSize)
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


def exportShapefile(
    polylines, attributes, EPSGCode=26909,
    saveAs='MyShape', label='AvgDepth', attType='int', directory="Output"
):

    """
        Function to export polylines to shape file with attribute

    """

    if not os.path.isdir(directory):
        os.makedirs(directory)

    crs = from_epsg(EPSGCode)

    # Define a polygon feature geometry with one attribute
    schema = {
        'geometry': 'LineString',
        'properties': {label: attType},

    }

    with fiona.open(
        directory + os.path.sep + saveAs + '.shp', 'w', 'ESRI Shapefile', schema, crs=crs
    ) as c:

        # If there are multiple geometries, put the "for" loop here
        for poly, att in zip(polylines, attributes):

            if np.all([np.any(att), len(poly) > 1]):
                pline = LineString(list(tuple(map(tuple, poly))))

                res = {}
                res['properties'] = {}
                res['properties'][label] = np.mean(att)
                # geometry of of the original polygon shapefile
                res['geometry'] = mapping(pline)
                c.write(res)
