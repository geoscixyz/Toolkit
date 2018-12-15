import numpy as np
import scipy as sp
from . import (Simulator, MathUtils)
from scipy.spatial import cKDTree
# from SimPEG.Utils import mkvc
import matplotlib.pyplot as plt
import gdal
import osr
import ogr
import os
import re
from shapely.geometry import mapping, LineString
import fiona
from fiona.crs import from_epsg
from download import download
import ipywidgets as widgets
import shapefile
import zipfile


gridProps = [
    'valuesFilledUC', 'valuesFilled',
    'derivativeX', 'derivativeY', 'firstVertical',
    'totalHorizontal', 'tiltAngle', 'analyticSignal',
    'gridFFT', 'gridPadded',
]


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
    heightUC = None
    inc = np.nan
    dec = np.nan
    indVal = None
    indNan = None
    EPSGcode = None

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

            self._gridCC = np.c_[X.flatten(order='F'), Y.flatten(order='F')]

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
                    self.indNan = np.isnan(self.values.flatten(order='F'))
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
    def values(self):
        """
            Data values
        """

        if getattr(self, '_values', None) is None:

            print("Values on the grid are not set")

            return

        if getattr(self, '_RTP', None) is not None:

            return self._RTP
        else:
            return self._values

    @property
    def valuesFilled(self):

        if getattr(self, '_valuesFilled', None) is None:
            values = self.values.flatten(order='F')

            # Do a minimum curvature extrapolation
            isNan = np.isnan(values)

            # Get the real values on the outer edge
            indVal = np.where(~isNan)[0]
            tree = cKDTree(self.gridCC[indVal, :])
            dists, indexes = tree.query(self.gridCC[isNan, :])

            uInd = np.unique(indVal[indexes])

            xyz = self.gridCC[uInd, :]

            _, values[isNan] = MathUtils.minCurvatureInterp(
                            xyz, values[uInd],
                            xyzOut=self.gridCC[isNan, :])

            # If there are still NDV, just do a neareest neighbour

            while np.isnan(values).sum() > 0:
                isNan = np.isnan(values)
                # Get the real values on the outer edge
                indVal = np.where(~isNan)[0]
                tree = cKDTree(self.gridCC[indVal, :])
                dists, indexes = tree.query(self.gridCC[isNan, :])
                values[isNan] = values[indVal[indexes]]

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
                derivX = derivX.flatten(order='F')

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
                derivY = derivY.flatten(order='F')

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
                firstVD = firstVD.flatten(order='F')

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

    def setRTP(self, isRTP):

        if isRTP:

            if np.isnan(self.inc):
                print("Attibute 'inc' needs to be set")
            if np.isnan(self.dec):
                print("Attibute 'dec' needs to be set")

            h0_xyz = MathUtils.dipazm_2_xyz(self.inc, self.dec)
            Frtp = self.gridFFT / (
                (h0_xyz[2] + 1j*(self.Kx*h0_xyz[0] + self.Ky*h0_xyz[1]))**2.
            )
            rtp_pad = np.fft.ifft2(Frtp)
            rtp = np.real(
                rtp_pad[self.npady:-self.npady, self.npadx:-self.npadx])

            if self.indNan is not None:
                rtp = rtp.flatten(order='F')

                rtp[self.indNan] = np.nan
                rtp = rtp.reshape(self.values.shape, order='F')

            self._RTP = rtp

            for prop in gridProps:
                setattr(self, '_{}'.format(prop), None)

        else:
            self._RTP = None
            for prop in gridProps:
                setattr(self, '_{}'.format(prop), None)


    def upwardContinuation(self, z=0):
        """
            Function to calculate upward continued data
        """

        self.heightUC = z
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
            zUpw = zUpw.flatten(order='F')

            zUpw[self.indNan] = np.nan
            zUpw = zUpw.reshape(self.values.shape, order='F')

        return zUpw


def loadGRDFile(fileName, plotIt=False):
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
        temp = grid.xyzv()[:, :, 3]
        temp[temp == -99999] = np.nan
        data._values = temp
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


def loadGeoTiffFile(fileName, plotIt=False):
    """
        Load a data matrix in Geosoft *.grd format and return
        a dictionary with gridded data
    """
    data = dataGrid()

    rasterObject = gdal.Open(fileName)
    band = rasterObject.GetRasterBand(1)
    temp = np.flipud(band.ReadAsArray())
    temp[temp == -99999] = np.nan
    data._values = temp
    data.nx = data.values.shape[1]
    data.ny = data.values.shape[0]
    data.x0 = rasterObject.GetGeoTransform()[0]
    y0 = rasterObject.GetGeoTransform()[3]
    data.dx = rasterObject.GetGeoTransform()[1]
    data.dy = np.abs(rasterObject.GetGeoTransform()[5])
    data.y0 = y0 - data.ny*data.dy
    data.limits = np.r_[data.x0, data.x0+data.nx*data.dx, data.y0, y0]

    proj = osr.SpatialReference(wkt=rasterObject.GetProjection())
    data.EPSGcode = proj.GetAttrValue('AUTHORITY', 1)
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


def writeGeotiff(
    array, fileName, EPSGcode,
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

    dataset.SetGeoTransform((xMin, pixelXSize, 0, yMax, 0, pixelYSize))

    datasetSRS = osr.SpatialReference()

    datasetSRS.ImportFromEPSG(int(EPSGcode))

    dataset.SetProjection(datasetSRS.ExportToWkt())

    if numBands == 1:
        dataset.GetRasterBand(1).WriteArray(np.flipud(array))
    else:
        for i in range(0, numBands):
            dataset.GetRasterBand(i+1).WriteArray(np.flipud(array[:, :, i]))

    dataset.FlushCache()  # Write to disk.


def readShapefile(fileName):

    world = shapefile.Reader(fileName)
    # Extract lines from shape file
    X, Y = [], []
    for shape in world.shapeRecords():

        for ii, part in enumerate(shape.shape.parts):

            if ii != len(shape.shape.parts)-1:
                x = [i[0] for i in shape.shape.points[shape.shape.parts[ii]:shape.shape.parts[ii+1]]]
                y = [i[1] for i in shape.shape.points[shape.shape.parts[ii]:shape.shape.parts[ii+1]]]

            else:
                x = [i[0] for i in shape.shape.points[shape.shape.parts[ii]:]]
                y = [i[1] for i in shape.shape.points[shape.shape.parts[ii]:]]

            if len(x) > 10:
                X.append(np.vstack(x))
                Y.append(np.vstack(y))

    return X, Y


def exportShapefile(
    polylines, attributes, EPSGcode=3156,
    saveAs='./Output/MyShape', label='AvgDepth', attType='int'
):

    """
        Function to export polylines to shape file with attribute

    """

    # if not os.path.isdir(directory):
    #     os.makedirs(directory)

    crs = from_epsg(int(EPSGcode))

    # Define a polygon feature geometry with one attribute
    schema = {
        'geometry': 'LineString',
        'properties': {label: attType},

    }

    with fiona.open(
        saveAs + '.shp', 'w', 'ESRI Shapefile', schema, crs=crs
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


def fetchData(
    path="./assets/Search/", checkDir=False, file="",
    localCloud='Local', dtype='CSV', loadDir=False, loadFile=False
):

    def fileLoader(
        localCloud, path, loadDir, checkDir, files,
        dtype, loadFile
    ):

        if loadDir:

            if localCloud == 'Cloud':
                print("Downloading... wait for it...")

                if "?dl=0" in path:
                    fileName = re.split('[/?]', path)[-2]
                else:
                    fileName = re.split('[/]', path)[-1]
                out = download(path, './Output/' + fileName, replace=True)
                path = './Output/'

                if '.zip' in fileName:
                    with zipfile.ZipFile(path+fileName) as zf:
                        zf.extractall(path)

                    if os.path.isdir(path+fileName):

                        path = path+fileName + os.path.sep

                    zf.close()

        if loadFile:
            print('Loading file:' + files)
            print('Please standby ...')
            if dtype == 'CSV':
                data = np.loadtxt(files)
                print('CSV file loaded. You will need to grid your data')

            elif dtype == 'GeoTiff':
                data = loadGeoTiffFile(files)
                print('Load complete')
            elif dtype == 'GRD':

                assert os.name == 'nt', "GRD file reader only available for Windows users. Sorry, you can complain to Geosoft"
                data = loadGRDFile(files)
                print('Load complete')
            return data, dtype

    def changeFileList(_):
        loadDir.value = False
        checkDir.value = False
        if localCloud.value == 'Cloud':
            lookinto = "./Output/"

        else:
            lookinto = path.value

        fileList = []
        for pathWalk, subdirs, fileNames in os.walk(lookinto):
            for name in fileNames:
               fileList += [os.path.join(pathWalk, name)]

        files.options = fileList

    if not isinstance(file, list):
        print(file)
        file = [file]

    def loadIt(_):

        if loadFile.value:
            loadFile.value = False

    localCloud = widgets.RadioButtons(
        options=['Local', 'Cloud'],
        description='File Type:',
        value=localCloud,
        disabled=False
    )
    path = widgets.Text(
        value=path,
        description='Path:',
        disabled=False
    )
    loadDir = widgets.ToggleButton(
        value=loadDir,
        description='Download',
        disabled=False,
        button_style='',
        tooltip='Fetch file on Cloud or Locally',
        icon='check'
    )
    checkDir = widgets.ToggleButton(
        value=checkDir,
        description='Check folder',
        disabled=False,
        button_style='',
        tooltip='Fetch files in Local folder',
        icon='check'
    )
    checkDir.observe(changeFileList)
    # files = os.listdir(path.value)

    loadFile = widgets.ToggleButton(
        value=loadFile,
        description='Load File',
        disabled=False,
        button_style='',
        tooltip='Load data to memory',
        icon='check'
    )
    # loadFile.observe(loadIt)
    files = widgets.Dropdown(
        options=file,
        description='Files:',
        disabled=False,
    )
    dtype = widgets.RadioButtons(
        options=['CSV', 'GeoTiff', 'GRD'],
        value=dtype,
        description='File Type:',
        disabled=False
    )
    out = widgets.interactive(fileLoader,
                              localCloud=localCloud,
                              path=path,
                              loadDir=loadDir,
                              checkDir=checkDir,
                              files=files,
                              dtype=dtype,
                              loadFile=loadFile
    )

    return out


def gdalWarp(fileNameOut, fileNameIn, EPSGcode):
    """
        Load a file and reproject
    """

    grid = gdal.Open(fileNameIn)

    gdal.Warp(fileNameOut, grid, dstSRS='EPSG:' + str(int(EPSGcode)))
