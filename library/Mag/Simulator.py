from . import Mag
from . import MathUtils
from . import DataIO
from . import ProblemSetter
import SimPEG.PF as PF
import shapefile
from SimPEG.Utils import mkvc
from scipy.constants import mu_0
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import ipywidgets as widgets
from ipywidgets import Layout
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.interpolate import griddata, interp1d, RegularGridInterpolator
import scipy.sparse as sp
from scipy.spatial import cKDTree
from scipy.interpolate.interpnd import _ndim_coords_from_arrays
from matplotlib.colors import LightSource, Normalize
from library.graphics import graphics
from matplotlib.ticker import FormatStrFormatter
import matplotlib as mpl
from skimage import exposure
from matplotlib.patches import Rectangle
import webbrowser
from osgeo import ogr, osr
import PIL

def PFSimulator(prism, survey):

    def PFInteract(update, susc, comp, irt, Q, RemInc, RemDec,
                   Profile_npt, Profile_azm, Profile_len,
                   Profile_ctx, Profile_cty):

        # # Get the line extent from the 2D survey for now
        prob = Mag.Problem()
        prob.prism = prism
        prob.survey = survey

        return PlotFwrSim(prob, susc, comp, irt, Q, RemInc, RemDec,
                          Profile_azm, Profile_len, Profile_npt,
                          Profile_ctx, Profile_cty)

    locs = survey.srcField.rxList[0].locs
    xlim = np.asarray([locs[:, 0].min(), locs[:, 0].max()])
    ylim = np.asarray([locs[:, 1].min(), locs[:, 1].max()])

    Lx = xlim[1] - xlim[0]
    Ly = ylim[1] - ylim[0]
    diag = (Lx**2. + Ly**2.)**0.5/2.

    ctx = np.mean(xlim)
    cty = np.mean(ylim)

    out = widgets.interactive(PFInteract,
                              update=widgets.ToggleButton(description='Refresh', value=False),
                              susc=widgets.FloatSlider(min=0, max=2, step=0.001, value=0.1, continuous_update=False),
                              comp=widgets.ToggleButtons(options=['tf', 'bx', 'by', 'bz']),
                              irt=widgets.ToggleButtons(options=['induced', 'remanent',  'total']),
                              Q=widgets.FloatSlider(min=0., max=10, step=1, value=0, continuous_update=False),
                              RemInc=widgets.FloatSlider(min=-90., max=90, step=5, value=0, continuous_update=False),
                              RemDec=widgets.FloatSlider(min=-90., max=90, step=5, value=0, continuous_update=False),
                              Profile_npt=widgets.BoundedFloatText(min=10, max=100, step=1, value=20, continuous_update=False),
                              Profile_azm=widgets.FloatSlider(min=-90, max=90, step=5, value=45., continuous_update=False),
                              Profile_len=widgets.FloatSlider(min=10, max=diag, step=10, value= Ly, continuous_update=False),
                              Profile_ctx=widgets.FloatSlider(value=ctx, min=xlim[0], max=xlim[1], step=0.1, continuous_update=False, color='black'),
                              Profile_cty=widgets.FloatSlider(value=cty, min=ylim[0], max=ylim[1], step=0.1, continuous_update=False, color='black'), )
    return out


def cmaps():
    """
      Return some pre-selected colormaps from matplotlib
    """

    return [
        'viridis', 'plasma', 'magma', 'RdBu_r',
        'Greys_r', 'jet', 'rainbow', 'pink',
        'bone', 'hsv', 'nipy_spectral'
        ]


def units():
    """
        Returns a dictionary of units for all filters
    """

    return {
        'derivativeX': '[nT/m]', 'derivativeY': '[nT/m]',
        'firstVertical': '[nT/m]', 'totalHorizontal': '[nT/m]',
        'tiltAngle': '[Degree]', 'analyticSignal': '[nT/m]',
        'RTP': '[nT]', 'TMI': '[nT]'
        }


def PlotFwrSim(prob, susc, comp, irt, Q, rinc, rdec,
               Profile_azm, Profile_len, Profile_npt,
               Profile_ctx, Profile_cty):

    def MagSurvey2D(x, y, data, Profile_ctx, Profile_cty, Profile_azm,
                    Profile_len, Profile_npt,
                    fig=None, ax=None,
                    vmin=None, vmax=None, pred=None):

        # Get the line extent from the 2D survey for now
        Profile_azm /= 180./np.pi
        Profile_len /= 2.*0.98

        dx = np.cos(-Profile_azm)*Profile_len
        dy = np.sin(-Profile_azm)*Profile_len

        a = [Profile_ctx - dx, Profile_cty - dy]
        b = [Profile_ctx + dx, Profile_cty + dy]

        return plotMagSurvey2D(x, y,
                               data, a, b, Profile_npt,
                               fig=fig, ax=ax,
                               vmin=vmin, vmax=vmax, pred=pred)

    def MagSurveyProfile(survey, Profile_ctx, Profile_cty, Profile_azm,
                         Profile_len, Profile_npt,
                         data=None, fig=None, ax=None):

        # Get the line extent from the 2D survey for now
        Profile_azm /= 180./np.pi
        Profile_len /= 2.*0.98

        dx = np.cos(-Profile_azm)*Profile_len
        dy = np.sin(-Profile_azm)*Profile_len

        a = [Profile_ctx - dx, Profile_cty - dy]
        b = [Profile_ctx + dx, Profile_cty + dy]

        xyz = survey.srcField.rxList[0].locs
        dobs = survey.dobs

        return plotProfile2D(xyz[:, 0], xyz[:, 1], [dobs, data], a, b, Profile_npt,
                             fig=fig, ax=ax, ylabel='nT')

    survey = prob.survey
    rxLoc = survey.srcField.rxList[0].locs
    prob.Q, prob.rinc, prob.rdec = Q, rinc, rdec
    prob.uType, prob.mType = comp, irt
    prob.susc = susc

    # Compute fields from prism
    fields = prob.fields()

    dpred = np.zeros_like(fields[0])
    for b in fields:
        dpred += b

    vmin = survey.dobs.min()
    vmax = survey.dobs.max()
    rxLoc = survey.srcField.rxList[0].locs
    x, y = rxLoc[:, 0], rxLoc[:, 1]

    f = plt.figure(figsize=(8, 8))

    ax0 = plt.subplot(1, 2, 1)
    MagSurvey2D(x, y, survey.dobs, Profile_ctx, Profile_cty, Profile_azm,
                Profile_len, Profile_npt, fig=f, ax=ax0, pred=dpred,
                vmin=survey.dobs.min(), vmax=survey.dobs.max())

    ax0 = plt.subplot(1, 2, 2)
    MagSurvey2D(x, y, dpred, Profile_ctx, Profile_cty, Profile_azm,
                Profile_len, Profile_npt, fig=f, ax=ax0, pred=dpred,
                vmin=survey.dobs.min(), vmax=survey.dobs.max())

    f = plt.figure(figsize=(12, 5))
    ax2 = plt.subplot()
    MagSurveyProfile(survey, Profile_ctx, Profile_cty, Profile_azm,
                     Profile_len, Profile_npt, data=dpred, fig=f, ax=ax2)

    plt.show()


def ViewMagSurveyWidget(survey, shapeFile=None):

    def MagSurvey2D(East, North, Azimuth, Length, Sampling, ):

        # Calculate the original map extents
        if isinstance(survey, DataIO.dataGrid):
            xLoc = survey.hx
            yLoc = survey.hy
            data = survey.values
        else:
            xLoc = survey.srcField.rxList[0].locs[:, 0]
            yLoc = survey.srcField.rxList[0].locs[:, 1]
            data = survey.dobs
        # Get the line extent from the 2D survey for now
        ColorMap = "RdBu_r"
        Azimuth = np.deg2rad((450 - Azimuth) % 360)
        Length /= 2.*0.98
        a = [East - np.cos(Azimuth)*Length, North - np.sin(Azimuth)*Length]
        b = [East + np.cos(Azimuth)*Length, North + np.sin(Azimuth)*Length]

        fig = plt.figure(figsize=(10, 6))
        ax1 = plt.subplot(1, 2, 1)

        plotMagSurvey2D(
         xLoc, yLoc, data, a, b, Sampling,
         fig=fig, ax=ax1, cmap=ColorMap, marker=False, shapeFile=shapeFile
        )

        # if Profile:
        ax2 = plt.subplot(1, 2, 2)
        plotProfile2D(xLoc, yLoc, data, a, b, Sampling,
                      fig=fig, ax=ax2, ylabel='nT')

        # ax2.set_aspect(0.5)
        pos = ax2.get_position()
        ax2.set_position([pos.x0, pos.y0+0.25, pos.width*2.0, pos.height*0.5])
        ax2.set_title('A', loc='left', fontsize=14)
        ax2.set_title("A'", loc='right', fontsize=14)

        plt.show()
        return survey

    # Calculate the original map extents
    if isinstance(survey, DataIO.dataGrid):
        xLoc = survey.hx
        yLoc = survey.hy
        data = survey.values
    else:
        xLoc = survey.srcField.rxList[0].locs[:, 0]
        yLoc = survey.srcField.rxList[0].locs[:, 1]
        data = survey.dobs

    Lx = xLoc.max() - xLoc.min()
    Ly = yLoc.max() - yLoc.min()
    diag = (Lx**2. + Ly**2.)**0.5

    cntr = [np.mean(xLoc), np.mean(yLoc)]

    out = widgets.interactive(
        MagSurvey2D,
        East=widgets.FloatSlider(min=cntr[0]-Lx, max=cntr[0]+Lx, step=10, value=cntr[0],continuous_update=False),
        North=widgets.FloatSlider(min=cntr[1]-Ly, max=cntr[1]+Ly, step=10, value=cntr[1],continuous_update=False),
        Azimuth=widgets.FloatSlider(min=0, max=180, step=5, value=90, continuous_update=False),
        Length=widgets.FloatSlider(min=20, max=diag, step=20, value=diag/2., continuous_update=False),
        Sampling=widgets.BoundedFloatText(min=10, max=1000, step=5, value=100, continuous_update=False)
        # ColorMap=widgets.Dropdown(
        #           options=cmaps(),
        #           value='RdBu_r',
        #           description='ColorMap',
        #           disabled=False,
        #         )
    )

    return out


def plotMagSurvey2D(x, y, data, a, b, npts, pred=None, marker=True,
                    fig=None, ax=None, vmin=None, vmax=None, shapeFile=None,
                    cmap='RdBu_r', equalizeHist='HistEqualized'):
    """
    Plot the data and line profile inside the spcified limits
    """

    if fig is None:
        fig = plt.figure()

    if ax is None:
        ax = plt.subplot()

    xLine, yLine = linefun(a[0], b[0], a[1], b[1], npts)

    # Use SimPEG.PF ploting function
    fig, im, cbar = plotData2D(
        x, y, data, fig=fig,  ax=ax,
        vmin=vmin, vmax=vmax,
        marker=marker, cmap=cmap,
        colorbar=False, equalizeHist=equalizeHist,
        shapeFile=shapeFile
    )

    ax.plot(xLine, yLine, 'k.', ms=5)
    cbar = plt.colorbar(im, orientation='horizontal')
    cbar.set_label('TMI (nT)')
    ax.text(xLine[0], yLine[0], 'A', fontsize=16, color='k', ha='left')
    ax.text(xLine[-1], yLine[-1], "A'", fontsize=16,
            color='k', ha='right')
    ax.grid(True)
    ax.set_xlim([x.min(), x.max()])
    ax.set_ylim([y.min(), y.max()])
    return


def linefun(x1, x2, y1, y2, nx, tol=1e-3):
    dx = x2-x1
    dy = y2-y1

    if np.abs(dx) < tol:
        y = np.linspace(y1, y2, nx)
        x = np.ones_like(y)*x1
    elif np.abs(dy) < tol:
        x = np.linspace(x1, x2, nx)
        y = np.ones_like(x)*y1
    else:
        x = np.linspace(x1, x2, nx)
        slope = (y2-y1)/(x2-x1)
        y = slope*(x-x1)+y1
    return x, y


def ViewPrism(survey):

    def Prism(update, dx, dy, dz, x0, y0, elev, prism_inc, prism_dec, View_dip, View_azm, View_lim):

        prism = definePrism()
        prism.dx, prism.dy, prism.dz, prism.z0 = dx, dy, dz, elev
        prism.x0, prism.y0 = x0, y0
        prism.pinc, prism.pdec = prism_inc, prism_dec

        # Display the prism and survey points
        plotObj3D([prism], survey, View_dip, View_azm, View_lim)

        return prism

    rxLoc = survey.srcField.rxList[0].locs
    cntr = np.mean(rxLoc[:, :2], axis=0)

    xlim = rxLoc[:, 0].max() - rxLoc[:, 0].min()
    ylim = rxLoc[:, 1].max() - rxLoc[:, 1].min()

    lim = np.max([xlim, ylim])/2.

    out = widgets.interactive(Prism,
                              update=widgets.ToggleButton(description='Refresh', value=False),
                              dx=widgets.FloatSlider(min=.01, max=1000., step=.01, value=lim/4, continuous_update=False),
                              dy=widgets.FloatSlider(min=.01, max=1000., step=.01, value=lim/4, continuous_update=False),
                              dz=widgets.FloatSlider(min=.01, max=1000., step=.01, value=lim/4, continuous_update=False),
                              x0=widgets.FloatSlider(min=cntr[0]-1000, max=cntr[0]+1000, step=1., value=cntr[0], continuous_update=False),
                              y0=widgets.FloatSlider(min=cntr[1]-1000, max=cntr[1]+1000, step=1., value=cntr[1], continuous_update=False),
                              elev=widgets.FloatSlider(min=-1000., max=1000., step=1., value=0., continuous_update=False),
                              prism_inc=(-90., 90., 5.),
                              prism_dec=(-90., 90., 5.),
                              View_dip=widgets.FloatSlider(min=0, max=90, step=1, value=30, continuous_update=False),
                              View_azm=widgets.FloatSlider(min=0, max=360, step=1, value=0, continuous_update=False),
                              View_lim=widgets.FloatSlider(min=1, max=2*lim, step=1, value=lim, continuous_update=False),
                              )

    return out


def plotObj3D(prisms, survey, View_dip, View_azm, View_lim, fig=None, axs=None, title=None, colors=None):

    """
    Plot the prism in 3D
    """

    rxLoc = survey.srcField.rxList[0].locs

    if fig is None:
        fig = plt.figure(figsize=(7, 7))

    if axs is None:
        axs = fig.add_subplot(111, projection='3d')

    if title is not None:
        axs.set_title(title)

    plt.rcParams.update({'font.size': 13})

    cntr = np.mean(rxLoc[:, :2], axis=0)

    axs.set_xlim3d(-View_lim + cntr[0], View_lim + cntr[0])
    axs.set_ylim3d(-View_lim + cntr[1], View_lim + cntr[1])
#     axs.set_zlim3d(depth+np.array(surveyArea[:2]))
    axs.set_zlim3d(rxLoc[:, 2].max()*1.1-View_lim*2, rxLoc[:, 2].max()*1.1)

    if colors is None:
        colors = ['w']*len(prisms)

    for prism, color in zip(prisms, colors):

        x1, x2 = prism.xn[0], prism.xn[1]
        y1, y2 = prism.yn[0], prism.yn[1]
        z1, z2 = prism.zn[0], prism.zn[1]
        pinc, pdec = prism.pinc, prism.pdec

        # Create a rectangular prism, rotate and plot
        block_xyz = np.asarray([[x1, x1, x2, x2, x1, x1, x2, x2],
                               [y1, y2, y2, y1, y1, y2, y2, y1],
                               [z1, z1, z1, z1, z2, z2, z2, z2]])

        xyz = MathUtils.rotate(block_xyz.T, np.r_[prism.xc, prism.yc, prism.zc], pinc, pdec)
        # R = MagUtils.rotationMatrix(pinc, pdec)

        # xyz = R.dot(block_xyz).T

        # Offset the prism to true coordinate
        # offx = prism.xc
        # offy = prism.yc
        # offz = prism.zc

        #print xyz
        # Face 1
        axs.add_collection3d(Poly3DCollection([list(zip(xyz[:4, 0],
                                                   xyz[:4, 1],
                                                   xyz[:4, 2]))]))

        # Face 2
        axs.add_collection3d(Poly3DCollection([list(zip(xyz[4:, 0],
                                                   xyz[4:, 1],
                                                   xyz[4:, 2]))], facecolors=color))

        # Face 3
        axs.add_collection3d(Poly3DCollection([list(zip(xyz[[0, 1, 5, 4], 0],
                                                   xyz[[0, 1, 5, 4], 1],
                                                   xyz[[0, 1, 5, 4], 2]))]))

        # Face 4
        axs.add_collection3d(Poly3DCollection([list(zip(xyz[[3, 2, 6, 7], 0],
                                                   xyz[[3, 2, 6, 7], 1],
                                                   xyz[[3, 2, 6, 7], 2]))]))

       # Face 5
        axs.add_collection3d(Poly3DCollection([list(zip(xyz[[0, 4, 7, 3], 0],
                                                   xyz[[0, 4, 7, 3], 1],
                                                   xyz[[0, 4, 7, 3], 2]))]))

       # Face 6
        axs.add_collection3d(Poly3DCollection([list(zip(xyz[[1, 5, 6, 2], 0],
                                                   xyz[[1, 5, 6, 2], 1],
                                                   xyz[[1, 5, 6, 2], 2]))]))


    axs.set_xlabel('Easting (X; m)')
    axs.set_ylabel('Northing (Y; m)')
    axs.set_zlabel('Depth (Z; m)')

    if survey.dobs is not None:
        color = survey.dobs
    else:
        color = 'k'
    axs.scatter(rxLoc[:, 0], rxLoc[:, 1], zs=rxLoc[:, 2], c=color, s=20, cmap='RdBu_r', zorder=100)

    # Convert from geographic
    azmDeg = (450 - View_azm) % 360 + 180

    axs.view_init(View_dip, azmDeg)
    plt.show()

    return True


class definePrism(object):
    """
        Define a prism and its attributes

        Prism geometry:
            - dx, dy, dz: width, length and height of prism
            - depth : depth to top of prism
            - susc : susceptibility of prism
            - x0, y0 : center of prism in horizontal plane
            - pinc, pdec : inclination and declination of prism
    """

    x0, y0, z0, dx, dy, dz = 0., 0., 0., 1., 1., 1.
    pinc, pdec = 0., 0.


    # Define the nodes of the prism
    @property
    def xn(self):
        xn = np.asarray([-self.dx/2. + self.x0, self.dx/2. + self.x0])

        return xn

    @property
    def yn(self):
        yn = np.asarray([-self.dy/2. + self.y0, self.dy/2. + self.y0])

        return yn

    @property
    def zn(self):
        zn = np.asarray([-self.dz + self.z0, self.z0])

        return zn

    @property
    def xc(self):
        xc = (self.xn[0] + self.xn[1]) / 2.

        return xc

    @property
    def yc(self):
        yc = (self.yn[0] + self.yn[1]) / 2.

        return yc

    @property
    def zc(self):
        zc = (self.zn[0] + self.zn[1]) / 2.

        return zc


def fitline(prism, survey):

    def profiledata(Binc, Bdec, Bigrf, depth,
                    susc, comp, irt, Q, rinc, rdec, update):

        # Get the line extent from the 2D survey for now
        prob = Mag.problem()
        prob.prism = prism.result

        xyzLoc = survey.srcField.rxList[0].locs.copy()
        xyzLoc[:, 2] += depth

        rxLoc = PF.BaseMag.RxObs(xyzLoc)
        srcField = PF.BaseMag.SrcField([rxLoc], param=[Bigrf, Binc, Bdec])
        survey2D = PF.BaseMag.LinearSurvey(srcField)
        survey2D.dobs = survey.dobs
        prob.survey = survey2D

        prob.Q, prob.rinc, prob.rdec = Q, rinc, rdec
        prob.uType, prob.mType = comp, irt
        prob.susc = susc

        # Compute fields from prism
        fields = prob.fields()

        dpred = np.zeros_like(fields[0])
        for b in fields:
            dpred += (b + Bigrf)

        a = np.r_[xyzLoc[:, 0].min(), 0]
        b = np.r_[xyzLoc[:, 0].max(), 0]
        return plotProfile2D(
                  xyzLoc, [survey2D.dobs, dpred], a, b, 10,
                  dType='2D', ylabel='nT',
                )

    Q = widgets.interactive(
        profiledata, Binc=widgets.FloatSlider(min=-90., max=90, step=5, value=90, continuous_update=False),
        Bdec=widgets.FloatSlider(min=-90., max=90, step=5, value=0, continuous_update=False),
        Bigrf=widgets.FloatSlider(min=54000., max=55000, step=10, value=54500, continuous_update=False),
        depth=widgets.FloatSlider(min=0., max=2., step=0.05, value=0.5),
        susc=widgets.FloatSlider(min=0.,  max=800., step=5.,  value=1.),
        comp=widgets.ToggleButtons(options=['tf', 'bx', 'by', 'bz']),
        irt=widgets.ToggleButtons(options=['induced', 'remanent', 'total']),
        Q=widgets.FloatSlider(min=0.,  max=10., step=0.1,  value=0.),
        rinc=widgets.FloatSlider(min=-180.,  max=180., step=1.,  value=0.),
        rdec=widgets.FloatSlider(min=-180.,  max=180., step=1.,  value=0.),
        update=widgets.ToggleButton(description='Refresh', value=False)
    )
    return Q


class MidPointNorm(Normalize):
    """
      Color range normalization based on a mid-point
      Provided from:
      https://stackoverflow.com/questions/7404116/defining-the-midpoint-of-a-colormap-in-matplotlib
    """
    def __init__(self, midpoint=None, vmin=None, vmax=None, clip=False):
        Normalize.__init__(self, vmin, vmax, clip)
        self.midpoint = midpoint

    def __call__(self, value, clip=None):
        if clip is None:
            clip = self.clip

        result, is_scalar = self.process_value(value)

        self.autoscale_None(result)

        if self.midpoint is None:
            self.midpoint = (self.vmin + self.vmax)/2.
        vmin, vmax, midpoint = self.vmin, self.vmax, self.midpoint

        if not (vmin < midpoint < vmax):
            raise ValueError("midpoint must be between maxvalue and minvalue.")
        elif vmin == vmax:
            result.fill(0) # Or should it be all masked? Or 0.5?
        elif vmin > vmax:
            raise ValueError("maxvalue must be bigger than minvalue")
        else:
            vmin = float(vmin)
            vmax = float(vmax)
            if clip:
                mask = np.ma.getmask(result)
                result = np.ma.array(np.clip(result.filled(vmax), vmin, vmax),
                                  mask=mask)

            # ma division is very slow; we can take a shortcut
            resdat = result.data

            # First scale to -1 to 1 range, than to from 0 to 1.
            resdat -= midpoint
            resdat[resdat > 0] /= abs(vmax - midpoint)
            resdat[resdat < 0] /= abs(vmin - midpoint)

            resdat /= 2.
            resdat += 0.5
            result = np.ma.array(resdat, mask=result.mask, copy=False)

        if is_scalar:
            result = result[0]
        return result

    def inverse(self, value):
        if not self.scaled():
            raise ValueError("Not invertible until scaled")
        vmin, vmax, midpoint = self.vmin, self.vmax, self.midpoint

        if cbook.iterable(value):
            val = ma.asarray(value)
            val = 2 * (val-0.5)
            val[val > 0] *= abs(vmax - midpoint)
            val[val < 0] *= abs(vmin - midpoint)
            val += midpoint
            return val
        else:
            val = 2 * (val - 0.5)
            if val < 0:
                return val*abs(vmin-midpoint) + midpoint
            else:
                return val*abs(vmax-midpoint) + midpoint


def plotDataHillside(x, y, z, axs=None, fill=True, contours=25,
                     vmin=None, vmax=None, levels=None, resolution=25,
                     clabel=True, cmap='RdBu_r', ve=1., alpha=0.5, alphaHS=0.5,
                     distMax=1000, midpoint=None, azdeg=315, altdeg=45,
                     equalizeHist='HistEqualized', minCurvature=True,
                     scatterData=None, shapeFile=None):

    ls = LightSource(azdeg=azdeg, altdeg=altdeg)

    if z.ndim == 1:

        if minCurvature:
            gridCC, d_grid = MathUtils.minCurvatureInterp(
                np.c_[x, y], z,
                vectorX=None, vectorY=None, vectorZ=None, gridSize=resolution,
                tol=1e-5, iterMax=None, method='spline',
            )
            X = gridCC[:, 0].reshape(d_grid.shape, order='F')
            Y = gridCC[:, 1].reshape(d_grid.shape, order='F')

        else:
            npts_x = int((x.max() - x.min())/Resolution)
            npts_y = int((y.max() - y.min())/Resolution)
            # Create grid of points
            vectorX = np.linspace(x.min(), x.max(), npts_x)
            vectorY = np.linspace(y.min(), y.max(), npts_y)

            Y, X = np.meshgrid(vectorY, vectorX)

            d_grid = griddata(np.c_[x, y], z, (X, Y), method='cubic')

        # Remove points beyond treshold
        tree = cKDTree(np.c_[x, y])
        xi = _ndim_coords_from_arrays((X, Y), ndim=2)
        dists, indexes = tree.query(xi)

        # Copy original result but mask missing values with NaNs
        d_grid[dists > distMax] = np.nan

    else:

        X, Y, d_grid = x, y, z

    im, CS = [], []
    if axs is None:
        axs = plt.subplot()

    if fill:

        if vmin is None:
            vmin = np.floor(d_grid[~np.isnan(d_grid)].min())

        if vmax is None:
            vmax = np.ceil(d_grid[~np.isnan(d_grid)].max())

        if equalizeHist == 'HistEqualized':

            subGrid = d_grid[~np.isnan(d_grid)]
            cdf, bins = exposure.cumulative_distribution(
                    subGrid[
                       (subGrid < vmax) *
                       (subGrid > vmin)
                       ].flatten(), nbins=256
            )
            my_cmap = graphics.equalizeColormap(cmap, bins, cdf)
        else:
            my_cmap = cmap

        extent = x.min(), x.max(), y.min(), y.max()
        im = axs.imshow(d_grid, vmin=vmin, vmax=vmax,
                       cmap=my_cmap, clim=[vmin, vmax],
                       alpha=alpha,
                       extent=extent, origin='lower')
        if np.all([alpha != 1, alphaHS != 0]):
            axs.imshow(ls.hillshade(d_grid, vert_exag=ve,
                       dx=resolution, dy=resolution),
                       cmap='gray_r', alpha=alphaHS,
                       extent=extent, origin='lower')

        # clevels = np.linspace(vmin, vmax, contours)
        # im = axs.contourf(
        #     X, Y, d_grid, contours, levels=clevels,
        #     cmap=my_cmap, alpha=alpha
        # )


    if levels is not None:
        CS = axs.contour(
            X, Y, d_grid, levels.shape[0],
            levels=levels, colors='k', linewidths=0.5
        )

        if clabel:
            plt.clabel(CS, inline=1, fontsize=10, fmt='%i')

    if scatterData is not None:
        plt.scatter(
          scatterData['x'], scatterData['y'],
          scatterData['size'], c=scatterData['c'],
          cmap=scatterData['cmap'],
          vmin=scatterData['clim'][0],
          vmax=scatterData['clim'][1]
        )

    if shapeFile is not None:
        plotShapeFile(shapeFile, ax=axs)

    axs.set_xlim([extent[0], extent[1]])
    axs.set_ylim([extent[2], extent[3]])
    axs.set_aspect('auto')
    return X, Y, d_grid, im, CS


def plotData2D(x, y, d, title=None,
               vmin=None, vmax=None, contours=None, fig=None, ax=None,
               colorbar=True, marker=True, cmap="RdBu_r",
               equalizeHist='HistEqualized', shapeFile=None):
    """ Function plot_obs(rxLoc,d)
    Generate a 2d interpolated plot from scatter points of data

    INPUT
    rxLoc       : Observation locations [x,y,z]
    d           : Data vector

    OUTPUT
    figure()

    Created on Dec, 27th 2015

    @author: dominiquef

    """

    from scipy.interpolate import griddata
    import pylab as plt

    # Plot result
    if fig is None:
        fig = plt.figure()

    if ax is None:
        ax = plt.subplot()

    if d.ndim == 1:
        assert x.shape[0] == d.shape[0], "Data and x locations must be consistant"

        assert y.shape[0] == d.shape[0], "Data and y locations must be consistant"


    plt.sca(ax)
    if marker:
        plt.scatter(x, y, c='k', s=10)

    if d is not None:

        ndv = np.isnan(d) == False
        if (vmin is None):
            vmin = d[ndv].min()

        if (vmax is None):
            vmax = d[ndv].max()

        # If data not a grid, create points evenly sampled
        if d.ndim == 1:
            # Create grid of points
            xGrid = np.linspace(x.min(), x.max(), 100)
            yGrid = np.linspace(y.min(), y.max(), 100)

            X, Y = np.meshgrid(xGrid, yGrid)
            d_grid = griddata(np.c_[x, y], d, (X, Y), method='linear')

        # Already a grid
        else:

            # if x.ndim == 1:
            #     X, Y = np.meshgrid(x, y)
            #     d_grid = d

            # else:
            X, Y, d_grid = x, y, d

        if equalizeHist == 'HistEqualized':
            cdf, bins = exposure.cumulative_distribution(
                d_grid[~np.isnan(d_grid)].flatten(), nbins=256
            )
            my_cmap = graphics.equalizeColormap(cmap, bins, cdf)
        else:

            my_cmap = cmap

        im = plt.imshow(
                d_grid, extent=[x.min(), x.max(), y.min(), y.max()],
                origin='lower', vmin=vmin, vmax=vmax, cmap=my_cmap
              )

        cbar = []
        if colorbar:
            cbar = plt.colorbar(fraction=0.02)

        if contours is not None:
            plt.contour(X, Y, d_grid, levels=contours, colors='k',
                        vmin=vmin, vmax=vmax)

    if title is not None:
        plt.title(title)

    if shapeFile is not None:
        plotShapeFile(shapeFile, ax=ax)

    plt.yticks(rotation='vertical')
    ylabel = np.round(np.linspace(y.min(), y.max(), 5) * 1e-3) * 1e+3
    ax.set_yticklabels(ylabel[1:4], size=12, rotation=90, va='center')
    ax.set_yticks(ylabel[1:4])
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    xlabel = np.round(np.linspace(x.min(), x.max(), 5) * 1e-3) * 1e+3
    ax.set_xticklabels(xlabel[1:4], size=12, va='center')
    ax.set_xticks(xlabel[1:4])
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax.set_xlabel('Easting')
    ax.set_ylabel('Northing')
    ax.grid(True)
    ax.set_aspect('equal')
    # plt.gca().set_aspect('equal', adjustable='box')

    return fig, im, cbar


def plotProfile2D(x, y, data, a, b, npts,
                  fig=None, ax=None, plotStr=['b', 'r'],
                  coordinate_system='local',
                  ylabel='Data'):
    """
    Plot the data and line profile inside the spcified limits
    """
    def linefun(x1, x2, y1, y2, nx, tol=1e-3):
        dx = x2-x1
        dy = y2-y1

        if np.abs(dx) <= tol:
            y = np.linspace(y1, y2, nx)
            x = np.ones_like(y)*x1
        elif np.abs(dy) <= tol:
            x = np.linspace(x1, x2, nx)
            y = np.ones_like(x)*y1
        else:
            x = np.linspace(x1, x2, nx)
            slope = (y2-y1)/(x2-x1)
            y = slope*(x-x1)+y1
        return x, y

    if fig is None:
        fig = plt.figure(figsize=(6, 9))

        plt.rcParams.update({'font.size': 14})

    if ax is None:
        ax = plt.subplot()

    xLine, yLine = linefun(a[0], b[0], a[1], b[1], npts)

    ind = (xLine > x.min()) * (xLine < x.max()) * (yLine > y.min()) * (yLine < y.max())

    xLine = xLine[ind]
    yLine = yLine[ind]

    distance = np.sqrt((xLine-a[0])**2.+(yLine-a[1])**2.)
    if coordinate_system == 'xProfile':
        distance += a[0]
    elif coordinate_system == 'yProfile':
        distance += a[1]

    if not isinstance(data, list):
        data = [data]

    for ii, d in enumerate(data):
        if d.ndim == 1:
            dline = griddata(np.c_[x, y], d, (xLine, yLine), method='linear')

        else:
            F = RegularGridInterpolator((x, y), d.T)
            dline = F(np.c_[xLine, yLine])

        # Check for nan
        ind = np.isnan(dline)==False

        if plotStr[ii]:
            ax.plot(distance[ind], dline[ind], plotStr[ii])
        else:
            ax.plot(distance[ind], dline[ind])

    # ax.set_xlim(distance.min(), distance.max())
    ax.set_ylabel(ylabel)
    ax.set_aspect('auto')
    ax.grid(True)
    return ax


def dataHillsideWidget(
    survey, EPSGcode=26909, HSTransp=0.5, SunAzimuth=270,
    saveAs='DataHillshade', dpi=300,
    scatterData=None, shapeFile=None,
  ):

    def plotWidget(
            SunAzimuth, SunAngle,
            ColorTransp, HSTransp, vScale,
            Contours, ColorMap, VminVmax, Equalize,
            SaveGeoTiff, saveAs
         ):

        # Calculate the original map extents
        if isinstance(survey, DataIO.dataGrid):
            xLoc = survey.hx
            yLoc = survey.hy
            data = survey.values

        else:
            xLoc = survey.srcField.rxList[0].locs[:, 0]
            yLoc = survey.srcField.rxList[0].locs[:, 1]
            data = survey.dobs

        if SaveGeoTiff:
            fig = plt.figure()
            fig.set_size_inches(7, 7)
            axs = plt.Axes(fig, [0., 0., 1., 1.])
            axs.set_axis_off()
            fig.add_axes(axs)

        else:
            fig = plt.figure(figsize=(7, 7))
            axs = plt.subplot()

        # Add shading
        X, Y, d_grid, im, CS = plotDataHillside(
          xLoc, yLoc, data,
          axs=axs, cmap=ColorMap,
          clabel=False, contours=Contours,
          vmax=VminVmax[1], vmin=VminVmax[0],
          alpha=ColorTransp, alphaHS=HSTransp,
          ve=vScale, azdeg=SunAzimuth, altdeg=SunAngle,
          equalizeHist=Equalize, scatterData=scatterData,
          shapeFile=shapeFile)

        # Add points at the survey locations
        # plt.scatter(xLoc, yLoc, s=2, c='k')
        if SaveGeoTiff:
            plt.savefig("Output/" + saveAs + '.png', dpi=dpi)
            plt.close()

            img = np.asarray(PIL.Image.open("Output/" + saveAs + '.png'))

            if survey.EPSGcode is not None:
                EPSGcode = survey.EPSGcode

            DataIO.arrayToRaster(
                img, "Output/" + saveAs + '.tiff',
                EPSGcode, np.min(X), np.max(X), np.min(Y), np.max(Y), 3
            )

        else:
            axs.set_aspect('equal')
            cbar = plt.colorbar(im, fraction=0.02)
            cbar.set_label('TMI (nT)')
            plt.yticks(rotation='vertical')
            ylabel = np.round(np.linspace(Y.min(), Y.max(), 5) * 1e-3) * 1e+3
            axs.set_yticklabels(ylabel[1:4], size=12, rotation=90, va='center')
            axs.set_yticks(ylabel[1:4])
            axs.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
            xlabel = np.round(np.linspace(X.min(), X.max(), 5) * 1e-3) * 1e+3
            axs.set_xticklabels(xlabel[1:4], size=12, va='center')
            axs.set_xticks(xlabel[1:4])
            axs.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
            axs.set_xlabel("Easting (m)", size=14)
            axs.set_ylabel("Northing (m)", size=14)
            axs.grid('on', color='k', linestyle='--')

            if scatterData is not None:
                pos = axs.get_position()
                cbarax = fig.add_axes([pos.x0+0.875, pos.y0+0.225,  pos.width*.025, pos.height*0.4])
                norm = mpl.colors.Normalize(vmin=scatterData['clim'][0], vmax=scatterData['clim'][1])
                cb = mpl.colorbar.ColorbarBase(
                  cbarax, cmap=scatterData['cmap'],
                  norm=norm,
                  orientation="vertical")
                cb.set_label("Depth (m)", size=12)

            plt.show()

    # Calculate the original map extents
    if isinstance(survey, DataIO.dataGrid):
        xLoc = survey.hx
        yLoc = survey.hy
        data = survey.values

    else:
        xLoc = survey.srcField.rxList[0].locs[:, 0]
        yLoc = survey.srcField.rxList[0].locs[:, 1]
        data = survey.dobs
    SunAzimuth = widgets.FloatSlider(
        min=0, max=360, step=5, value=SunAzimuth, continuous_update=False,
        description='SunAzimuth'
        )
    SunAngle = widgets.FloatSlider(
        min=0, max=90, step=5, value=15, continuous_update=False,
        description='SunAngle'
        )
    ColorTransp = widgets.FloatSlider(
        min=0, max=1, step=0.05, value=0.9, continuous_update=False,
        description='ColorTransp'
        )
    HSTransp = widgets.FloatSlider(
        min=0, max=1, step=0.05, value=HSTransp, continuous_update=False,
        description='HSTransp'
        )

    vScale = widgets.FloatSlider(
        min=1, max=10, step=1., value=5.0, continuous_update=False,
        description='vScale'
        )

    Contours = widgets.IntSlider(
        min=10, max=100, step=10, value=50, continuous_update=False,
        description='Contours'
        )
    ColorMap = widgets.Dropdown(
        options=cmaps(),
        value='RdBu_r',
        description='ColorMap',
        disabled=False,
    )

    VminVmax = widgets.FloatRangeSlider(
                                    value=[data.min(), data.max()],
                                    min=data.min(),
                                    max=data.max(),
                                    step=1.0,
                                    description='Color Range',
                                    disabled=False,
                                    continuous_update=False,
                                    orientation='horizontal',
                                    readout=True,
                                    readout_format='.1f',
                                )

    Equalize = widgets.Dropdown(
                                  options=['Linear', 'HistEqualized'],
                                  value='HistEqualized',
                                  description='Color Normalization',
                                  disabled=False,
                                )

    SaveGeoTiff = widgets.ToggleButton(
                                  value=False,
                                  description='Export geoTiff',
                                  disabled=False,
                                  button_style='',
                                  tooltip='Description',
                                  icon='check'
                                )

    saveAs = widgets.Text(
        value="MyGeoTiff",
        description='GeoTiff name:',
        disabled=False
        )
    out = widgets.interactive_output(plotWidget,
                              {
                                    'SunAzimuth': SunAzimuth,
                                    'SunAngle': SunAngle,
                                    'ColorTransp': ColorTransp,
                                    'HSTransp': HSTransp,
                                    'vScale': vScale,
                                    'Contours': Contours,
                                    'ColorMap': ColorMap,
                                    'VminVmax': VminVmax,
                                    'Equalize': Equalize,
                                    'saveAs': saveAs,
                                    'SaveGeoTiff': SaveGeoTiff
                                }
                              )

    left = widgets.VBox(
            [SunAzimuth, SunAngle, ColorTransp, HSTransp, vScale,
             Contours, ColorMap, VminVmax, Equalize, saveAs, SaveGeoTiff],
            layout=Layout(
                width='35%', height='400px', margin='60px 0px 0px 0px'
            )
        )

    image = widgets.VBox(
              [out],
              layout=Layout(
                  width='65%', height='400px', margin='0px 0px 0px 0px'
              )
            )

    return widgets.HBox([left, out])


def gridFiltersWidget(
    survey, gridFilter='derivativeX',
    saveAs=None, ColorTransp=0.9, HSTransp=0.5,
    EPSGcode=np.nan, dpi=300, scatterData=None,
    inc=np.nan, dec=np.nan,
    SunAzimuth=90, SunAngle=15, vScale=5.,
    ColorMap='RdBu_r', shapeFile=None
):

    gridProps = [
        'valuesFilledUC', 'valuesFilled',
        'derivativeX', 'derivativeY', 'firstVertical',
        'totalHorizontal', 'tiltAngle', 'analyticSignal',
        'gridFFT', 'gridPadded',
      ]

    def plotWidget(
            SunAzimuth, SunAngle,
            ColorTransp, HSTransp, vScale,
            ColorMap, Filters, UpDist, SaveGrid, saveAs
         ):

        # If changed upward distance, reset the FFT
        if UpDist != survey.heightUC:
            for prop in gridProps:
                    setattr(survey, '_{}'.format(prop), None)

            data = survey.upwardContinuation(z=UpDist)
            survey._gridPadded = None
            survey._gridFFT = None

        if Filters == 'TMI':
            data = survey.upwardContinuation(z=UpDist)
        else:
            data = getattr(survey, '{}'.format(Filters))

        ind = ~np.isnan(data)
        vmin, vmax = np.percentile(data[ind], 5), np.percentile(data[ind], 95)

        vScale *= np.abs(survey.values[ind].max() - survey.values[ind].min()) * np.abs(data[ind].max() - data[ind].min())

        plotIt(
            data, SunAzimuth, SunAngle,
            ColorTransp, HSTransp, vScale,
            ColorMap, Filters, vmin, vmax, 'HistEqualized', SaveGrid, saveAs
        )

    def plotIt(
            data, SunAzimuth, SunAngle,
            ColorTransp, HSTransp, vScale,
            ColorMap, Filters, vmin, vmax, equalizeHist, SaveGrid, saveAs
         ):

        if SaveGrid:
            fig = plt.figure()
            fig.set_size_inches(7, 7)
            axs = plt.Axes(fig, [0., 0., 1., 1.])
            axs.set_axis_off()
            fig.add_axes(axs)

        else:

            fig = plt.figure(figsize=(7, 7))
            axs = plt.subplot()

        # Add shading
        X, Y, data, im, CS = plotDataHillside(
            survey.hx, survey.hy, data,
            axs=axs, cmap=ColorMap,
            clabel=False, resolution=10,
            vmin=vmin, vmax=vmax, contours=50,
            alpha=ColorTransp, alphaHS=HSTransp,
            ve=vScale, azdeg=SunAzimuth, altdeg=SunAngle,
            equalizeHist=equalizeHist, scatterData=scatterData,
            shapeFile=shapeFile
        )

        if SaveGrid:

            if saveAs is None:
                saveAs = Filters

            plt.savefig("Output/" + saveAs + '.png', dpi=dpi)
            plt.close()

            img = np.asarray(PIL.Image.open("Output/" + saveAs + '.png'))

            if survey.EPSGcode is not None:
                EPSGcode = survey.EPSGcode

            DataIO.arrayToRaster(
                img, "Output/" + saveAs + '.tiff',
                EPSGcode, np.min(X), np.max(X), np.min(Y), np.max(Y), 3
            )

        else:
            # Add points at the survey locations
            # plt.scatter(xLoc, yLoc, s=2, c='k')
            axs.set_aspect('equal')
            cbar = plt.colorbar(im, fraction=0.02)
            cbar.set_label(Filters + " " +units()[Filters])

            axs.set_xlabel("Easting (m)", size=14)
            axs.set_ylabel("Northing (m)", size=14)
            axs.grid('on', color='k', linestyle='--')

            plt.show()

    assert isinstance(survey, DataIO.dataGrid), "Only implemented for objects of class DataIO.dataGrid"
    SunAzimuth = widgets.FloatSlider(
        min=0, max=360, step=5, value=SunAzimuth,
        continuous_update=False,
        description='SunAzimuth'
        )
    SunAngle = widgets.FloatSlider(
        min=0, max=90, step=5, value=SunAngle,
        description='SunAngle', continuous_update=False
        )
    ColorTransp = widgets.FloatSlider(
        min=0, max=1, step=0.05, value=ColorTransp,
        description='ColorTransp', continuous_update=False
        )
    HSTransp = widgets.FloatSlider(
        min=0, max=1, step=0.05, value=HSTransp,
        description='HSTransp', continuous_update=False
        )
    vScale = widgets.FloatSlider(
        min=1, max=10, step=1., value=vScale,
        description='vScale', continuous_update=False
        )
    ColorMap = widgets.Dropdown(
        options=cmaps(),
        value=ColorMap,
        description='ColorMap',
        disabled=False,
        )
    Filters = widgets.Dropdown(
        options=[
            'TMI',
            'derivativeX', 'derivativeY', 'firstVertical',
            'totalHorizontal', 'tiltAngle', 'analyticSignal',
            'RTP'],
        value=gridFilter,
        description='Grid Filters',
        disabled=False,
        )
    UpDist = widgets.FloatSlider(
        min=0, max=200, step=10, value=0,
        continuous_update=False, description='UpC Height'
        )
    SaveGrid = widgets.ToggleButton(
        value=False,
        description='Export Grid',
        disabled=False,
        button_style='',
        tooltip='Description',
        icon='check'
        )
    saveAs = widgets.Text(
        value="MyGeoTiff",
        description='GeoTiff name:',
        disabled=False
        )

    out = widgets.interactive_output(plotWidget,
                            {
                              'SunAzimuth': SunAzimuth,
                              'SunAngle': SunAngle,
                              'ColorTransp': ColorTransp,
                              'HSTransp': HSTransp,
                              'vScale': vScale,
                              'ColorMap': ColorMap,
                              'Filters': Filters,
                              'UpDist': UpDist,
                              'saveAs': saveAs,
                              'SaveGrid': SaveGrid,
                            }
                        )

    left = widgets.VBox(
            [SunAzimuth, SunAngle, ColorTransp, HSTransp, vScale,
             ColorMap, Filters, UpDist, saveAs, SaveGrid],
            layout=Layout(
                width='35%', height='400px', margin='60px 0px 0px 0px'
            )
        )

    image = widgets.VBox(
              [out],
              layout=Layout(
                  width='65%', height='400px', margin='0px 0px 0px 0px'
              )
            )

    return widgets.HBox([left, out])

    return out


def gridTilt2Depth(
    survey, gridFilter='tiltAngle',
    GridFileName=None, ColorTransp=0.9, HSTransp=0.5,
    ShapeFileName=None,
    EPSGcode=26909, dpi=300, scatterData=None,
    SunAzimuth=90, SunAngle=15, vScale=5., shapeFile=None,
    ColorMap='RdBu_r', ColorDepth='viridis_r', depthRange=[0, 500],
    markerSize=1
):

    gridProps = [
        'valuesFilledUC', 'valuesFilled',
        'derivativeX', 'derivativeY', 'firstVertical',
        'totalHorizontal', 'tiltAngle', 'analyticSignal',
        'RTP', 'gridFFT', 'gridPadded',
      ]

    def plotWidget(
            SunAzimuth, SunAngle,
            ColorTransp, HSTransp, vScale,
            ColorMap, Filters, UpDist,
            ContourColor, ContourSize,
            GridFileName, SaveGrid,
            ShapeFileName, SaveShape,

         ):

        # If changed upward distance, reset the FFT
        if UpDist != survey.heightUC:
            for prop in gridProps:
                    setattr(survey, '_{}'.format(prop), None)

            data = survey.upwardContinuation(z=UpDist)
            survey._gridPadded = None
            survey._gridFFT = None

        if Filters == 'TMI':
            data = survey.upwardContinuation(z=UpDist)
        else:
            data = getattr(survey, '{}'.format(Filters))

        ind = ~np.isnan(data)
        vmin, vmax = np.percentile(data[ind], 5), np.percentile(data[ind], 95)

        # Compute estimated depth
        polylines, attributes = MathUtils.estimateDepth(survey)

        if SaveShape:
            # Export to shapefile
            DataIO.exportShapefile(polylines, attributes, EPSGcode=EPSGcode, saveAs=ShapeFileName, label='AvgDepth')

        scatterData = {}
        scatterData['x'] = np.vstack(polylines)[:, 0]
        scatterData['y'] = np.vstack(polylines)[:, 1]
        scatterData['size'] = ContourSize
        scatterData['c'] = np.concatenate(attributes)-UpDist
        scatterData['cmap'] = ContourColor
        scatterData['clim'] = [np.percentile(scatterData['c'], 25), np.percentile(scatterData['c'], 75)]

        vScale *= np.abs(survey.values[ind].max() - survey.values[ind].min()) * np.abs(data[ind].max() - data[ind].min())

        plotIt(
            data, SunAzimuth, SunAngle,
            ColorTransp, HSTransp, vScale,
            ColorMap, Filters, vmin, vmax, 'HistEqualized',
            SaveGrid, GridFileName,
            scatterData, shapeFile
        )

    def plotIt(
            data, SunAzimuth, SunAngle,
            ColorTransp, HSTransp, vScale,
            ColorMap, Filters, vmin, vmax, equalizeHist, SaveGrid, saveAs,
            scatterData, shapeFile
         ):

        if SaveGrid:
            fig = plt.figure()
            fig.set_size_inches(7, 7)
            axs = plt.Axes(fig, [0., 0., 1., 1.])
            axs.set_axis_off()
            fig.add_axes(axs)

        else:

            fig = plt.figure(figsize=(7, 7))
            axs = plt.subplot()

        # Add shading
        X, Y, data, im, CS = plotDataHillside(
            survey.hx, survey.hy, data,
            axs=axs, cmap=ColorMap,
            clabel=False, resolution=10,
            vmin=vmin, vmax=vmax, contours=50,
            alpha=ColorTransp, alphaHS=HSTransp,
            ve=vScale, azdeg=SunAzimuth, altdeg=SunAngle,
            equalizeHist=equalizeHist, scatterData=scatterData,
            shapeFile=shapeFile
        )

        if SaveGrid:

            if saveAs is None:
                saveAs = Filters

            plt.savefig("Output/" + saveAs + '.png', dpi=dpi)
            plt.close()

            img = np.asarray(PIL.Image.open("Output/" + saveAs + '.png'))

            if survey.EPSGcode is not None:
                EPSGcode = survey.EPSGcode

            DataIO.arrayToRaster(
                img, "Output/" + saveAs + '.tiff',
                EPSGcode, np.min(X), np.max(X), np.min(Y), np.max(Y), 3
            )

        else:
            # Add points at the survey locations
            # plt.scatter(xLoc, yLoc, s=2, c='k')
            axs.set_aspect('equal')
            cbar = plt.colorbar(im, fraction=0.02)
            # cbar.set_label(Filters + " " +units()[Filters])
            plt.yticks(rotation='vertical')
            ylabel = np.round(np.linspace(Y.min(), Y.max(), 5) * 1e-3) * 1e+3
            axs.set_yticklabels(ylabel[1:4], size=12, rotation=90, va='center')
            axs.set_yticks(ylabel[1:4])
            axs.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
            xlabel = np.round(np.linspace(X.min(), X.max(), 5) * 1e-3) * 1e+3
            axs.set_xticklabels(xlabel[1:4], size=12, va='center')
            axs.set_xticks(xlabel[1:4])
            axs.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
            axs.set_xlabel("Easting (m)", size=14)
            axs.set_ylabel("Northing (m)", size=14)
            axs.grid('on', color='k', linestyle='--')

            pos = axs.get_position()
            cbarax = fig.add_axes([pos.x0+0.875, pos.y0+0.225,  pos.width*.025, pos.height*0.4])
            norm = mpl.colors.Normalize(vmin=scatterData['clim'][0], vmax=scatterData['clim'][1])
            cb = mpl.colorbar.ColorbarBase(
              cbarax, cmap=scatterData['cmap'],
              norm=norm,
              orientation="vertical")
            cb.set_label("Depth (m)", size=12)

    assert isinstance(survey, DataIO.dataGrid), "Only implemented for objects of class DataIO.dataGrid"

    SunAzimuth = widgets.FloatSlider(
        min=0, max=360, step=5, value=SunAzimuth,
        continuous_update=False,
        description='SunAzimuth'
        )
    SunAngle = widgets.FloatSlider(
        min=0, max=90, step=5, value=SunAngle,
        description='SunAngle', continuous_update=False
        )
    ColorTransp = widgets.FloatSlider(
        min=0, max=1, step=0.05, value=ColorTransp,
        description='ColorTransp', continuous_update=False
        )
    HSTransp = widgets.FloatSlider(
        min=0, max=1, step=0.05, value=HSTransp,
        description='HSTransp', continuous_update=False
        )
    vScale = widgets.FloatSlider(
        min=1, max=10, step=1., value=vScale,
        description='vScale', continuous_update=False
        )
    ColorMap = widgets.Dropdown(
        options=cmaps(),
        value=ColorMap,
        description='ColorMap',
        disabled=False,
        )
    Filters = widgets.Dropdown(
        options=[
            'TMI',
            'tiltAngle'],
        value=gridFilter,
        description='Grid Filters',
        disabled=False,
        )
    UpDist = widgets.FloatSlider(
        min=0, max=200, step=10, value=0,
        continuous_update=False, description='UpC Height'
        )
    SaveGrid = widgets.ToggleButton(
        value=False,
        description='Export GeoTiff',
        disabled=False,
        button_style='',
        tooltip='Description',
        icon='check'
        )
    SaveShape = widgets.ToggleButton(
        value=False,
        description='Export Shapefile',
        disabled=False,
        button_style='',
        tooltip='Description',
        icon='check'
        )
    GridFileName = widgets.Text(
        value="MyGeoTiff",
        description='GeoTiff name:',
        disabled=False
        )
    ShapeFileName = widgets.Text(
        value="EstimatedDepth",
        description='Shapefile name:',
        disabled=False
        )
    ContourColor = widgets.Dropdown(
        options=cmaps(),
        value='viridis',
        description='ContourColor',
        disabled=False,
        )
    ContourSize = widgets.FloatSlider(
        min=1, max=10, step=1., value=1,
        description='ContourSize', continuous_update=False
        )

    out = widgets.interactive_output(plotWidget,
                            {
                              'SunAzimuth': SunAzimuth,
                              'SunAngle': SunAngle,
                              'ColorTransp': ColorTransp,
                              'HSTransp': HSTransp,
                              'vScale': vScale,
                              'ColorMap': ColorMap,
                              'Filters': Filters,
                              'UpDist': UpDist,
                              'ContourColor': ContourColor,
                              'ContourSize': ContourSize,
                              'GridFileName': GridFileName,
                              'SaveGrid': SaveGrid,
                              'ShapeFileName': ShapeFileName,
                              'SaveShape': SaveShape
                            }
                        )

    left = widgets.VBox(
            [SunAzimuth, SunAngle, ColorTransp, HSTransp, vScale,
             ColorMap, Filters, UpDist, GridFileName,
             SaveGrid,
             ContourColor, ContourSize, ShapeFileName, SaveShape],
            layout=Layout(
                width='35%', margin='0px 0px 0px 0px'
            )
        )

    # right = widgets.VBox(
    #         [ContourColor, ContourSize, ShapeFileName, SaveShape],
    #         layout=Layout(
    #             width='35%', margin='0px 0px 0px 0px'
    #         )
    #     )

    image = widgets.VBox(
              [out],
              layout=Layout(
                  width='65%', height='400px', margin='0px 0px 0px 0px'
              )
            )

    return widgets.HBox([left, image])

    return out


def plotShapeFile(shapeFile, ax=None, fill=True, linewidth=1):

    X, Y = DataIO.readShapefile(shapeFile)

    if ax is None:
        ax = plt.subplot()

    for x, y in zip(X, Y):
        ax.plot(x, y, 'k', linewidth=linewidth)

    return ax


def worldViewerWidget(worldFile, data, grid, z=0, shapeFile=None):

    def plotLocs(placeID):

        selection = int(np.r_[[ii for ii, s in enumerate(list(data.keys())) if placeID in s]])
        dataVals = list(data.values())[selection]

        Xloc, Yloc = np.meshgrid(grid.hx[::5], grid.hy[::5])
        Zloc = np.ones_like(Xloc)*z

        locs = np.c_[mkvc(Xloc), mkvc(Yloc), mkvc(Zloc)]
        survey, _, _ = ProblemSetter.setSyntheticProblem(locs, EarthField=dataVals[-3:])

        xyz = survey.srcField.rxList[0].locs
        plt.figure(figsize=(10, 8))
        ax1 = plt.subplot(1, 2, 1)
        fig, im, cbar = plotData2D(
          xyz[:, 0], xyz[:, 1], survey.dobs,
          ax=ax1, cmap='RdBu_r', marker=False, colorbar=False,
          shapeFile=shapeFile
        )

        ax1.set_xticks([0])
        ax1.set_xticklabels([MathUtils.decimalDegrees2DMS(dataVals[1], "Longitude")])
        ax1.set_xlabel('Longitude')
        ax1.set_yticks([0])
        ax1.set_yticklabels([MathUtils.decimalDegrees2DMS(dataVals[0], "Latitude")])
        ax1.set_ylabel('Latitude')
        ax1.grid(True)
        ax1.set_xlim([Xloc.min(), Xloc.max()])
        ax1.set_ylim([Yloc.min(), Yloc.max()])
        cbar = plt.colorbar(im, orientation='horizontal')
        cbar.set_label('TMI (nT)')

        axs = plt.subplot(1, 2, 2)
        axs = plotShapeFile(worldFile, ax=axs, fill=False)
        # axs.set_axis_off()
        axs.set_aspect('equal')
        pos = axs.get_position()
        axs.set_position([pos.x0, pos.y0,  pos.width*1.5, pos.height*1.5])
        axs.patch.set_alpha(0.0)
        # xydata = np.loadtxt("./assets/country-capitals.csv", delimiter=",")
        for key, entry in zip(list(data.keys()), list(data.values())):
            axs.scatter(entry[1], entry[0], c='k')
        axs.scatter(dataVals[1], dataVals[0], s=50, c='r', marker='s', )

        axs.set_title(
          "Earth's Field: " + str(int(dataVals[-3])) + "nT, "
          "Inc: " + str(int(dataVals[-2])) + "$^\circ$, "
          "Dec: " + str(int(dataVals[-1])) + "$^\circ$"
        )

        # Add axes with rotating arrow
        pos = axs.get_position()
        arrowAxs = fig.add_axes(
          [7, 7,  pos.width*.5, pos.height*0.5], projection='3d'
        )
        block_xyz = np.asarray([
                        [-.2, -.2, .2, .2, 0],
                        [-.25, -.25, -.25, -.25, 0.5],
                        [-.2, .2, .2, -.2, 0]
                    ])

        # rot = Utils.mkvc(Utils.dipazm_2_xyz(pinc, pdec))

        # xyz = Utils.rotatePointsFromNormals(block_xyz.T, np.r_[0., 1., 0.], rot,
        #                                     np.r_[p.xc, p.yc, p.zc])

        R = MathUtils.rotationMatrix(dataVals[-2], dataVals[-1])

        xyz = R.dot(block_xyz).T

        #print xyz
        # Face 1
        arrowAxs.add_collection3d(Poly3DCollection([list(zip(xyz[[1, 2, 4], 0],
                                                   xyz[[1, 2, 4], 1],
                                                   xyz[[1, 2, 4], 2]))],
                                                   facecolors='w')
                                  )

        arrowAxs.add_collection3d(Poly3DCollection([list(zip(xyz[[0, 1, 4], 0],
                                                   xyz[[0, 1, 4], 1],
                                                   xyz[[0, 1, 4], 2]))],
                                                   facecolors='k')
                                  )

        arrowAxs.add_collection3d(Poly3DCollection([list(zip(xyz[[2, 3, 4], 0],
                                                   xyz[[2, 3, 4], 1],
                                                   xyz[[2, 3, 4], 2]))],
                                                   facecolors='w')
                                  )

        arrowAxs.add_collection3d(Poly3DCollection([list(zip(xyz[[0, 3, 4], 0],
                                                   xyz[[0, 3, 4], 1],
                                                   xyz[[0, 3, 4], 2]))],
                                                   facecolors='k')
                                  )

        arrowAxs.add_collection3d(Poly3DCollection([list(zip(xyz[:4, 0],
                                                   xyz[:4, 1],
                                                   xyz[:4, 2]))],
                                                   facecolors='r')
                                  )

        arrowAxs.view_init(30, -90)
        arrowAxs.set_xlim([-0.5, 0.5])
        arrowAxs.set_ylim([-0.5, 0.5])
        arrowAxs.set_zlim([-0.5, 0.5])
        arrowAxs.set_xticks([])
        arrowAxs.set_yticks([])
        arrowAxs.set_zticks([])
        # arrowAxs.set_aspect('equal')

        plt.show()

        return axs

    out = widgets.interactive(plotLocs,
                              placeID=widgets.Dropdown(
                                options=list(data.keys()),
                                value=list(data.keys())[0],
                                description='Location:',
                                disabled=False,
                                )
                              )

    return out


def dataGriddingWidget(
    survey, EPSGcode=np.nan, saveAs="MyGeoTiff",
    shapeFile=None, inc=np.nan, dec=np.nan,
    Method='linear'
):

    def plotWidget(
            Resolution, Method,
            ColorMap,
            EPSGcode, inc, dec,
            GetIncDec, saveAs, SaveGrid
         ):

        if Method == 'minimumCurvature':
            gridCC, d_grid = MathUtils.minCurvatureInterp(
                np.c_[xLoc, yLoc], data,
                gridSize=Resolution, method='spline'
                )
            X = gridCC[:, 0].reshape(d_grid.shape, order='F')
            Y = gridCC[:, 1].reshape(d_grid.shape, order='F')

        else:
            npts_x = int((xLoc.max() - xLoc.min())/Resolution)
            npts_y = int((yLoc.max() - yLoc.min())/Resolution)

            # Create grid of points
            vectorX = np.linspace(xLoc.min(), xLoc.max(), npts_x)
            vectorY = np.linspace(yLoc.min(), yLoc.max(), npts_y)

            X, Y = np.meshgrid(vectorX, vectorY)

            d_grid = griddata(np.c_[xLoc, yLoc], data, (X, Y), method=Method)

        gridOut = DataIO.dataGrid()

        gridOut._values = d_grid
        gridOut.nx, gridOut.ny = gridOut.values.shape[1], gridOut.values.shape[0]
        gridOut.x0, gridOut.y0 = X.min(), Y.min()
        gridOut.dx = (X.max() - X.min()) / gridOut.values.shape[1]
        gridOut.dy = (Y.max() - Y.min()) / gridOut.values.shape[0]
        gridOut.limits = np.r_[gridOut.x0, gridOut.x0+gridOut.nx*gridOut.dx, gridOut.y0, gridOut.y0+gridOut.ny*gridOut.dy]

        if not np.isnan(EPSGcode):
            gridOut.EPSGcode = int(EPSGcode)
        gridOut.inc, gridOut.dec = inc, dec

        if SaveGrid:
            if np.isnan(EPSGcode):
                print("Need to assign an EPSGcode before exporting")
                return
            DataIO.arrayToRaster(
                d_grid, saveAs + '.tiff',
                gridOut.EPSGcode, X.min(), X.max(),
                Y.min(), Y.max(), 1,
                dataType='grid')

        else:
            fig = plt.figure(figsize=(7, 7))
            axs = plt.subplot()
            # Add shading
            X, Y, d_grid, im, CS = plotDataHillside(
                X, Y, d_grid, alpha=1., contours=Contours,
                axs=axs, cmap=ColorMap, clabel=False, shapeFile=shapeFile)

            # Add points at the survey locations
            plt.scatter(xLoc, yLoc, s=2, c='k')
            axs.set_aspect('auto')
            cbar = plt.colorbar(im, fraction=0.02)
            cbar.set_label('TMI (nT)')
            plt.yticks(rotation='vertical')
            ylabel = np.round(np.linspace(Y.min(), Y.max(), 5) * 1e-3) * 1e+3
            axs.set_yticklabels(ylabel[1:4], size=12, rotation=90, va='center')
            axs.set_yticks(ylabel[1:4])
            axs.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
            xlabel = np.round(np.linspace(X.min(), X.max(), 5) * 1e-3) * 1e+3
            axs.set_xticklabels(xlabel[1:4], size=12, va='center')
            axs.set_xticks(xlabel[1:4])
            axs.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
            axs.set_xlabel("Easting (m)", size=14)
            axs.set_ylabel("Northing (m)", size=14)
            axs.grid('on', color='k', linestyle='--')
            plt.show()

        # Create grid object
        return gridOut

    # Calculate the original map extents
    xLoc = survey[:, 0]
    yLoc = survey[:, 1]
    data = survey[:, -1]

    def fetchURL(_):
        if GetIncDec.value:
            GetIncDec.value = False
            if np.isnan(EPSGcode.value):
                print("Enter EPSGcode first")
                return

            x, y, z = np.mean(xLoc), np.mean(yLoc), 0.
            # input SpatialReference
            inSpatialRef = osr.SpatialReference()
            inSpatialRef.ImportFromEPSG(int(EPSGcode.value))

            # output SpatialReference
            outSpatialRef = osr.SpatialReference()
            outSpatialRef.ImportFromEPSG(4326)

            # create the CoordinateTransformation
            coordTrans = osr.CoordinateTransformation(inSpatialRef, outSpatialRef)

            latLon = coordTrans.TransformPoint(x, y, z)
            url = (
                "https://www.ngdc.noaa.gov/geomag-web/calculators/calculateIgrfwmm?lat1=" +
                str(latLon[1]) + "&lon1=" + str(latLon[0]) +
                "&model=IGRF&startYear=2000&endYear=2019&resultFormat=html"
                )
            out = webbrowser.open(url)

    Resolution = widgets.FloatText(
        value=25,
        description='Grid (m):',
        disabled=False
        )
    Method = widgets.Dropdown(
        options=[
          'nearest', 'linear', 'cubic',
          'minimumCurvature'
          ],
        value=Method,
        description='Method',
        disabled=False,
        )
    Contours = widgets.IntSlider(
        min=10, max=100, step=10,
        value=50, continuous_update=False
        )
    ColorMap = widgets.Dropdown(
        options=cmaps(),
        value='RdBu_r',
        description='ColorMap',
        disabled=False,
        )
    EPSGcode = widgets.FloatText(
        value=EPSGcode,
        description='EPSG code:',
        disabled=False
    )
    inc = widgets.FloatText(
        value=inc,
        description='Inclination angle positive downward from horizontal:',
        disabled=False
        )
    dec = widgets.FloatText(
        value=dec,
        description='Declination angle positve clockwise from North:',
        disabled=False
        )
    GetIncDec = widgets.ToggleButton(
        value=False,
        description='Fetch Inc/Dec',
        disabled=False,
        button_style='',
        tooltip='Connect to NOAA API',
        icon='check'
        )

    GetIncDec.observe(fetchURL)
    saveAs = widgets.Text(
        value=saveAs,
        description='GeoTiff name:',
        disabled=False
        )
    SaveGrid = widgets.ToggleButton(
        value=False,
        description='Export Grid',
        disabled=False,
        button_style='',
        tooltip='Write file',
        icon='check'
        )

    out = widgets.interactive(plotWidget,
                              Resolution=Resolution,
                              Method=Method,
                              ColorMap=ColorMap,
                              EPSGcode=EPSGcode,
                              inc=inc, dec=dec,
                              GetIncDec=GetIncDec,
                              saveAs=saveAs,
                              SaveGrid=SaveGrid
                              )

    return out


def dataGridGeoref(
    survey, EPSGcode=np.nan, saveAs="MyGeoTiff",
    shapeFile=None, inc=np.nan, dec=np.nan, applyRTP=False,
):

    def plotWidget(
            ColorMap,
            EPSGcode, inc, dec,
            GetIncDec, applyRTP, saveAs, SaveGrid
         ):

        if not np.isnan(EPSGcode):
            survey.EPSGcode = int(EPSGcode)

        survey.inc, survey.dec = inc, dec

        survey.setRTP(applyRTP)

        if SaveGrid:
            if np.isnan(EPSGcode):
                print("Need to assign an EPSGcode before exporting")
                return
            DataIO.arrayToRaster(
                survey.values, saveAs + '.tiff',
                survey.EPSGcode, X.min(), X.max(),
                Y.min(), Y.max(), 1,
                dataType='grid')

        else:
            fig = plt.figure(figsize=(7, 7))
            axs = plt.subplot()
            # Add shading
            X, Y, d_grid, im, CS = plotDataHillside(
                survey.hx, survey.hy, survey.values, alpha=1.,
                axs=axs, cmap=ColorMap, clabel=False, shapeFile=shapeFile)

            # Add points at the survey locations
            # plt.scatter(xLoc, yLoc, s=2, c='k')
            axs.set_aspect('auto')
            cbar = plt.colorbar(im, fraction=0.02)
            cbar.set_label('TMI (nT)')
            plt.yticks(rotation='vertical')
            ylabel = np.round(np.linspace(Y.min(), Y.max(), 5) * 1e-3) * 1e+3
            axs.set_yticklabels(ylabel[1:4], size=12, rotation=90, va='center')
            axs.set_yticks(ylabel[1:4])
            axs.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
            xlabel = np.round(np.linspace(X.min(), X.max(), 5) * 1e-3) * 1e+3
            axs.set_xticklabels(xlabel[1:4], size=12, va='center')
            axs.set_xticks(xlabel[1:4])
            axs.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
            axs.set_xlabel("Easting (m)", size=14)
            axs.set_ylabel("Northing (m)", size=14)
            axs.grid('on', color='k', linestyle='--')
            plt.show()

        # Create grid object
        return survey

    assert isinstance(survey, DataIO.dataGrid), "Only implemented for objects of class DataIO.dataGrid"

    def fetchURL(_):
        if GetIncDec.value:
            GetIncDec.value = False
            if np.isnan(EPSGcode.value):
                print("Enter EPSGcode first")
                return

            x, y, z = np.mean(survey.hx), np.mean(survey.hy), 0.
            # input SpatialReference
            inSpatialRef = osr.SpatialReference()
            inSpatialRef.ImportFromEPSG(int(EPSGcode.value))

            # output SpatialReference
            outSpatialRef = osr.SpatialReference()
            outSpatialRef.ImportFromEPSG(4326)

            # create the CoordinateTransformation
            coordTrans = osr.CoordinateTransformation(inSpatialRef, outSpatialRef)

            latLon = coordTrans.TransformPoint(x, y, z)
            url = (
                "https://www.ngdc.noaa.gov/geomag-web/calculators/calculateIgrfwmm?lat1=" +
                str(latLon[1]) + "&lon1=" + str(latLon[0]) +
                "&model=IGRF&startYear=2000&endYear=2019&resultFormat=html"
                )
            out = webbrowser.open(url)

    ColorMap = widgets.Dropdown(
        options=cmaps(),
        value='RdBu_r',
        description='ColorMap',
        disabled=False,
        )

    if survey.EPSGcode is not None:
        EPSGcode = survey.EPSGcode

    EPSGcode = widgets.FloatText(
        value=EPSGcode,
        description='EPSG code:',
        disabled=False
    )
    inc = widgets.FloatText(
        value=inc,
        description='Inclination angle positive downward from horizontal:',
        disabled=False
        )
    dec = widgets.FloatText(
        value=dec,
        description='Declination angle positve clockwise from North:',
        disabled=False
        )
    GetIncDec = widgets.ToggleButton(
        value=False,
        description='Fetch Inc/Dec',
        disabled=False,
        button_style='',
        tooltip='Connect to NOAA API',
        icon='check'
        )

    GetIncDec.observe(fetchURL)

    applyRTP = widgets.ToggleButton(
        value=applyRTP,
        description='Reduce to pole',
        disabled=False,
        button_style='',
        tooltip='Transform to RTP data',
        icon='check'
        )

    saveAs = widgets.Text(
        value=saveAs,
        description='GeoTiff name:',
        disabled=False
        )
    SaveGrid = widgets.ToggleButton(
        value=False,
        description='Export Grid',
        disabled=False,
        button_style='',
        tooltip='Write file',
        icon='check'
        )

    out = widgets.interactive(plotWidget,
                              ColorMap=ColorMap,
                              EPSGcode=EPSGcode,
                              inc=inc, dec=dec,
                              GetIncDec=GetIncDec,
                              applyRTP=applyRTP,
                              saveAs=saveAs,
                              SaveGrid=SaveGrid
                              )

    return out

def setDataExtentWidget(survey, East=None, North=None):
    """
        Small application to carve out a subset of a larger data set
    """

    def dataSelector(East, North, Width, Height):

        lims = np.r_[
            East-Width/2, East+Width/2,
            North-Height/2, North+Height/2
        ]

        fig, axs = plt.figure(figsize=(12, 6)), plt.subplot(1, 2, 1)
        plotData2D(
            xLoc, yLoc, survey.values, marker=False, fig=fig, ax=axs,
            colorbar=False
        )

        # axs.scatter(East, North, 20, marker='+')
        axs.add_patch(
            Rectangle(
                (East-Width/2, North-Height/2),
                Width,
                Height,
                facecolor='none', edgecolor='k',
                zorder=3
                )
            )

        # Extract data within window and plot
        indx = np.logical_and(xLoc > lims[0], xLoc < lims[1])

        indy = np.logical_and(yLoc > lims[2], yLoc < lims[3])

        nx, ny = np.count_nonzero(indx), np.count_nonzero(indy)

        # Create new dataGrid object
        dataSub = DataIO.dataGrid()
        dataSub.limits = lims
        # coordinate_system = grid.coordinate_system
        temp = survey.values[:, indx]
        temp = temp[indy, :]

        dataSub._values = temp
        dataSub.nx, dataSub.ny = nx, ny
        dataSub.dx, dataSub.dy = survey.dx, survey.dy
        dataSub.x0, dataSub.y0 = East-Width/2, North-Height/2

        # fig,
        axs = plt.subplot(1, 2, 2)
        fig, im, cbar = plotData2D(
            xLoc[indx], yLoc[indy], dataSub.values, marker=False, fig=fig, ax=axs
        )
        cbar.set_label('TMI (nT)')
        return dataSub

    if isinstance(survey, DataIO.dataGrid):

        xLoc = np.asarray(range(survey.nx))*survey.dx+survey.x0
        yLoc = np.asarray(range(survey.ny))*survey.dy+survey.y0
        xlim = survey.limits[:2]
        ylim = survey.limits[2:]

        if East is None:
            East = np.mean(xLoc)

        if North is None:
            North = np.mean(yLoc)

    else:
        print("Only implemented for class 'dataGrid'")
        # xLoc = survey.srcField.rxList[0].locs[:, 0]
        # yLoc = survey.srcField.rxList[0].locs[:, 1]
        # xlim = np.asarray([xLoc.min(), xLoc.max()])
        # ylim = np.asarray([yLoc.min(), yLoc.max()])
        # data = survey.dobs

    out = widgets.interactive(
            dataSelector,
            East=widgets.FloatSlider(min=xlim[0], max=xlim[1], step=500, value=East, continuous_update=False),
            North=widgets.FloatSlider(min=ylim[0], max=ylim[1], step=10, value=North, continuous_update=False),
            Width=widgets.FloatSlider(min=1000, max=100000, step=1000, value=30000, continuous_update=False),
            Height=widgets.FloatSlider(min=1000, max=100000, step=1000, value=30000, continuous_update=False)
            )

    return out
