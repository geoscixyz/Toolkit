import numpy as np
import scipy as sp
from . import (Simulator, MathUtils)
from scipy.spatial import cKDTree
from scipy import ndimage
from matplotlib.contour import QuadContourSet
import matplotlib.pyplot as plt
from osgeo import gdal, osr
import os
import re
from shapely.geometry import mapping, LineString
from fiona.crs import from_epsg
from download import download
import ipywidgets as widgets
import fiona
import zipfile

gridProps = [
    'valuesFilledUC', 'valuesFilled',
    'derivativeX', 'derivativeY', 'firstVertical',
    'totalHorizontal', 'tiltAngle', 'analyticSignal',
    'TDXderivative', 'gridFFT', 'gridPadded', 'RTP'
]

try:
    from scipy.fft import fftfreq

except ModuleNotFoundError:
    from scipy.fftpack import fftfreq


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
    fourier_gaussian = 0

    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)

        return

    def set_gaussian_filter(self, sigma):

        for prop in gridProps:
            setattr(self, '_{}'.format(prop), None)

        self.fourier_gaussian = sigma / self.dx

        return

    @property
    def hx(self):
        """
            Cell center location along x-axis
        """
        if getattr(self, '_hx', None) is None:
            self._hx = (
                np.asarray(range(self.nx)) *
                self.dx + self.x0 + self.dx/2
            )

        return self._hx

    @property
    def hy(self):
        """
            Cell center location along y-axis
        """

        if getattr(self, '_hy', None) is None:
            self._hy = (
                np.asarray(range(self.ny)) *
                self.dy + self.y0 + self.dy/2
            )

        return self._hy

    @property
    def gridCC(self):
        """
            Cell center location [nC x 2]
        """
        if getattr(self, '_gridCC', None) is None:

            X, Y = np.meshgrid(self.hx, self.hy)

            self._gridCC = np.c_[X.flatten(order='F'), Y.flatten(order='F')]

        return self._gridCC

    @property
    def npadx(self):
        """
            Number of padding cells  along x-axis for FFT
        """
        if getattr(self, '_npadx', None) is None:
            self._npadx = int(np.floor(self.values.shape[1]))

        return self._npadx

    @property
    def npady(self):
        """
            Number of padding cells  along y-axis for FFT
        """
        if getattr(self, '_npady', None) is None:
            self._npady = int(np.floor(self.values.shape[0]))

        return self._npady

    @property
    def gridPadded(self):
        """
            Grid values padded with mirror image cosin tapper
        """
        if getattr(self, '_gridPadded', None) is None:
            if getattr(self, '_valuesFilledUC', None) is not None:
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

            if self.fourier_gaussian > 0:

                return ndimage.gaussian_filter(
                    self._values, sigma=self.fourier_gaussian
                    )

            else:

                return self._values

    @property
    def valuesFilled(self):
        """
            Data values with no-data-values interpolated by minimum curvature
        """
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
                            xyzOut=self.gridCC[isNan, :],
                            overlap=np.max([self.dx, self.dy]))

            # If there are still NDV, just do a nearest neighbour

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
        """
            Data values with no-data-values interpolated by minimum curvature
            and upward continued.
        """
        if getattr(self, '_valuesFilledUC', None) is None:
            self._valuesFilledUC = self.upwardContinuation()

        return self._valuesFilledUC

    @property
    def Kx(self):
        """
            Wavenumber along the x-axis.
        """
        if getattr(self, '_Kx', None) is None:

            nx = self.gridPadded.shape[1]
            kx = 2. * np.pi * fftfreq(nx, self.dx)

            ny = self.gridPadded.shape[0]
            ky = 2. * np.pi * fftfreq(ny, self.dy)

            Ky, Kx = np.meshgrid(ky, kx)
            self._Ky, self._Kx = Ky.T, Kx.T

        return self._Kx

    @property
    def Ky(self):
        """
            Wavenumber along the y-axis.
        """
        if getattr(self, '_Ky', None) is None:

            nx = self.gridPadded.shape[1]
            kx = 2. * np.pi * fftfreq(nx, self.dx)

            ny = self.gridPadded.shape[0]
            ky = 2. * np.pi * fftfreq(ny, self.dy)

            Ky, Kx = np.meshgrid(ky, kx)
            self._Ky, self._Kx = Ky.T, Kx.T

        return self._Ky

    @property
    def gridFFT(self):
        """
            Padded grid values in the Fourier domain
        """
        if getattr(self, '_gridFFT', None) is None:

            self._gridFFT = np.fft.fft2(self.gridPadded)

        return self._gridFFT

    @property
    def derivativeX(self):
        """
            Derivative of grid along the x-axis
        """
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
        """
            Derivative of grid along the y-axis
        """
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
        """
            First vertical derivative of grid
        """
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
        """
            Total horizontal derivative of grid
        """
        if getattr(self, '_totalHorizontal', None) is None:

            self._totalHorizontal = np.sqrt(
                self.derivativeX**2. + self.derivativeY**2. + 1e-8
            )

        return self._totalHorizontal

    @property
    def tiltAngle(self):
        """
            Tilt angle of grid
        """
        if getattr(self, '_tiltAngle', None) is None:

            self._tiltAngle = np.arctan2(
                self.firstVertical, self.totalHorizontal
            )

        return self._tiltAngle

    @property
    def analyticSignal(self):
        """
            Analytic signal of grid
        """
        if getattr(self, '_analyticSignal', None) is None:

            self._analyticSignal = np.sqrt(
                self.derivativeX**2. +
                self.derivativeY**2. +
                self.firstVertical**2. + 1e-8
            )

        return self._analyticSignal

    @property
    def TDXderivative(self):
        """
            First vertical derivative of grid
        """
        if getattr(self, '_TDXderivative', None) is None:

            self._TDXderivative = np.arctan2(
                self.totalHorizontal, self.firstVertical
            )

        return self._TDXderivative

    def setRTP(self, isRTP):
        """
            Reduce to pole the grid
        """
        if isRTP:

            if np.isnan(self.inc):
                print("Attibute 'inc' needs to be set")
            if np.isnan(self.dec):
                print("Attibute 'dec' needs to be set")

            h_xyz = MathUtils.dipazm_2_xyz(-self.inc, self.dec)

            a = 2*h_xyz[0]*h_xyz[2]
            b = 2*h_xyz[1]*h_xyz[2]

            R = (self.Kx**2. + self.Ky**2.) / (
                (h_xyz[2]**2 - h_xyz[0]**2) * self.Kx**2. +
                (h_xyz[2]**2 - h_xyz[1]**2) * self.Ky**2. -
                (2*h_xyz[0]*h_xyz[1]) * self.Kx * self.Ky +
                1j * (self.Kx**2. + self.Ky**2.)**0.5 *
                (a*self.Kx + b*self.Ky) + 1e-8
            )

            R[0, 0] = 0
            Frtp = R * self.gridFFT
            rtp_pad = np.fft.ifft2(Frtp)
            rtp = np.real(
                rtp_pad[self.npady:-self.npady, self.npadx:-self.npadx])

            if self.indNan is not None:
                rtp = rtp.flatten(order='F')

                rtp[self.indNan] = np.nan
                rtp = rtp.reshape(self.values.shape, order='F')

            for prop in gridProps:
                setattr(self, '_{}'.format(prop), None)

            self._RTP = rtp

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
            z/2
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
        data.x0, data.y0 = grid.x0-grid.dx/2, grid.y0-grid.dy/2

    if plotIt:
        xLoc = np.asarray(range(data.nx))*data.dx+data.x0
        yLoc = np.asarray(range(data.ny))*data.dy+data.y0

        fig, axs = plt.figure(figsize=(8, 8)), plt.subplot()
        fig, im, cbar = Simulator.plotData2D(
            xLoc, yLoc, data.values, marker=False, fig=fig, ax=axs,
            colorbar=True
        )
        axs.grid(True)
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
    temp = band.ReadAsArray()

    data.nx = temp.shape[1]
    data.ny = temp.shape[0]

    data.x0 = rasterObject.GetGeoTransform()[0]

    data.dx = rasterObject.GetGeoTransform()[1]
    data.dy = np.abs(rasterObject.GetGeoTransform()[5])

    # temp[temp == -99999] = np.nan
    if np.sign(rasterObject.GetGeoTransform()[5]) == -1:
        ymax = rasterObject.GetGeoTransform()[3]
        data.y0 = ymax - data.ny*data.dy
        data._values = np.flipud(temp)
    else:
        data.y0 = rasterObject.GetGeoTransform()[3]
        ymax = data.y0 + data.ny*data.dy
        data._values = temp

    data.limits = np.r_[data.x0, data.x0+data.nx*data.dx, data.y0, ymax]

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
    pixelXSize = (xMax-xMin)/(xPixels)  # size of the pixel in X direction
    pixelYSize = -(yMax-yMin)/(yPixels)  # size of the pixel in Y direction

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

    shape = fiona.open(fileName)
    # Extract lines from shape file
    X, Y = [], []
    for item in shape.items():
        xy = item[1]['geometry']['coordinates']

        if len(xy) > 10:
            xy = np.vstack(xy)
            X.append(xy[:, 0])
            Y.append(xy[:, 1])

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

    if isinstance(polylines, QuadContourSet):

        temp, attributes = [], []
        for segs, level in zip(polylines.allsegs, polylines.levels):

            for poly in segs:

                temp += [poly]
                attributes += [level]

        polylines = temp

    with fiona.open(
        saveAs + '.shp', 'w', driver='ESRI Shapefile', schema=schema, crs=crs
    ) as c:

        # If there are multiple geometries, put the "for" loop here
        for poly, att in zip(polylines, attributes):

            if np.all([len(poly) > 0, len(poly) > 1]):
                pline = LineString(list(tuple(map(tuple, poly))))

                res = {}
                res['properties'] = {}
                res['properties'][label] = np.mean(att)

                # geometry of of the original polygon shapefile
                res['geometry'] = mapping(pline)
                c.write(res)

    # function to generate .prj file information using spatialreference.org
    def getWKT_PRJ(epsg_code):
        import urllib
        # access projection information
        wkt = urllib.request.urlopen("http://spatialreference.org/ref/epsg/{0}/prettywkt/".format(epsg_code))
        # remove spaces between charachters
        remove_spaces = wkt.read().replace(b" ", b"")
        # place all the text on one line
        output = remove_spaces.replace(b"\n", b"")
        return output

    # create the .prj file
    prj = open(saveAs + ".prj", "w")
    # call the function and supply the epsg code
    epsg = getWKT_PRJ(str(int(EPSGcode)))
    prj.write(epsg.decode("utf-8"))
    prj.close()

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
