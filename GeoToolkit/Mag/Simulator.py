from . import Mag
from . import MathUtils
from . import DataIO
from . import ProblemSetter
import re
from matplotlib import pyplot as plt
import numpy as np
import ipywidgets as widgets
from ipywidgets import Layout
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.interpolate import griddata, RegularGridInterpolator
from scipy.spatial import cKDTree
from scipy.interpolate.interpnd import _ndim_coords_from_arrays
from matplotlib.colors import LightSource, Normalize
from GeoToolkit.graphics import graphics
from matplotlib.ticker import FormatStrFormatter
import matplotlib as mpl
from skimage import exposure
from matplotlib.patches import Rectangle
import webbrowser
from osgeo import osr
import os
import PIL

np.seterr(divide='ignore', invalid='ignore')

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

    locs = survey.rxLoc
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
        'viridis', 'plasma', 'magma', 'Spectral_r',
        'Greys_r', 'jet', 'rainbow', 'pink', 'RdBu_r',
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
        'TDXderivative': '[Degree]','RTP': '[nT]', 'TMI': '[nT]'
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

        xyz = survey.rxLoc
        dobs = survey.dobs

        return plotProfile2D(xyz[:, 0], xyz[:, 1], [dobs, data], a, b, Profile_npt,
                             fig=fig, ax=ax, ylabel='nT')

    survey = prob.survey
    rxLoc = survey.rxLoc
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
    rxLoc = survey.rxLoc
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
            xLoc = survey.rxLoc[:, 0]
            yLoc = survey.rxLoc[:, 1]
            data = survey.dobs
        # Get the line extent from the 2D survey for now
        ColorMap = "Spectral_r"
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
        xLoc = survey.rxLoc[:, 0]
        yLoc = survey.rxLoc[:, 1]
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
        #           value='Spectral_r',
        #           description='ColorMap',
        #           style={'description_width': 'initial'},
        #         )
    )

    return out


def plotMagSurvey2D(x, y, data, a, b, npts, pred=None, marker=True,
                    fig=None, ax=None, vmin=None, vmax=None, shapeFile=None,
                    cmap='Spectral_r', equalizeHist='HistEqualized'):
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
    ax.text(xLine[0], yLine[0], 'A', fontsize=16, color='k', ha='right')
    ax.text(xLine[-1], yLine[-1], "A'", fontsize=16,
            color='k', ha='left')
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

    rxLoc = survey.rxLoc
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

    rxLoc = survey.rxLoc

    if fig is None:
        fig = plt.figure(figsize=(9, 9))

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
    axs.scatter(rxLoc[:, 0], rxLoc[:, 1], zs=rxLoc[:, 2], c=color, s=20, cmap='Spectral_r', zorder=100)

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

        xyzLoc = survey.rxLoc.copy()
        xyzLoc[:, 2] += depth

        # rxLoc = PF.BaseMag.RxObs(xyzLoc)
        # srcField = PF.BaseMag.SrcField([rxLoc], param=[Bigrf, Binc, Bdec])
        # survey2D = PF.BaseMag.LinearSurvey(srcField)

        survey2D = Mag.Survey(np.c_[Bigrf, Binc, Bdec])
        survey2D._rxLoc = xyz

        survey2D._dobs = survey.dobs
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


def plotDataHillside(x, y, z, axs=None, fill=True, contours=None,
                     vmin=None, vmax=None, resolution=25,
                     clabel=True, cmap='Spectral_r', ve=1., alpha=0.5, alphaHS=0.5,
                     distMax=1000, midpoint=None, azdeg=315, altdeg=45,
                     equalizeHist='HistEqualized', minCurvature=True,
                     scatterData=None, shapeFile=None):

    ls = LightSource(azdeg=azdeg, altdeg=altdeg)

    if z.ndim == 1:

        if minCurvature:

            gridCC, d_grid = MathUtils.minCurvatureInterp(
                np.c_[x, y], z,
                vectorX=None, vectorY=None, vectorZ=None,
                gridSize=resolution,
                tol=1e-5, iterMax=None, method='spline',
            )
            X = gridCC[:, 0].reshape(d_grid.shape, order='F')
            Y = gridCC[:, 1].reshape(d_grid.shape, order='F')

        else:
            npts_x = int((x.max() - x.min())/resolution)
            npts_y = int((y.max() - y.min())/resolution)
            # Create grid of points
            vectorX = np.linspace(x.min(), x.max(), npts_x)
            vectorY = np.linspace(y.min(), y.max(), npts_y)

            X, Y = np.meshgrid(vectorX, vectorY)

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

    extent = x.min(), x.max(), y.min(), y.max()

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


        im = axs.imshow(d_grid, vmin=vmin, vmax=vmax,
                       cmap=my_cmap, clim=[vmin, vmax],
                       alpha=alpha,
                       extent=extent, origin='lower')

        if np.all([alpha != 1, alphaHS != 0]):

            axs.imshow(ls.hillshade(d_grid, vert_exag=ve,
                       dx=resolution, dy=resolution),
                       cmap='gray_r', alpha=alphaHS,
                       extent=extent, origin='lower')

    if contours is not None:
        # clevels = np.round(np.linspace(vmin, vmax, contours) * 1e-1) * 1e+1

        # if np.all(clevels == 0):
        #     clevels = np.linspace(vmin, vmax, contours)

        # clevels = np.unique(clevels)
        # # Insert zero contour
        # if ~np.any(clevels == 0):
        #     clevels = np.sort(np.r_[clevels, 0])
        CS = axs.contour(
            X, Y, d_grid, len(contours), levels=contours,
            colors='k', linewidths=0.5
        )

            # plt.clabel(CS, inline=1, fontsize=5, fmt='%i')

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
               vmin=None, vmax=None, contours=0, fig=None, ax=None,
               colorbar=True, marker=True, cmap="Spectral_r",
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

        if contours > 0:
            clevels = np.round(np.linspace(vmin, vmax, contours) * 1e-1) * 1e+1

            if np.all(clevels == 0):
                clevels = np.linspace(vmin, vmax, contours)

            clevels = np.unique(clevels)
            # Insert zero contour
            if ~np.any(clevels == 0):
                clevels = np.sort(np.r_[clevels, 0])
            CS = axs.contour(
                X, Y, d_grid, contours, levels=clevels,
                colors='k', linewidths=0.5
            )
            # plt.clabel(CS, inline=1, fontsize=5, fmt='%i')

    if title is not None:
        plt.title(title)

    if shapeFile is not None:
        plotShapeFile(shapeFile, ax=ax)

    plt.yticks(rotation='vertical')
    roundFact = 10**(np.floor(np.log10(np.abs(y.max() - y.min()))) - 2)
    ylabel = np.round(np.linspace(y.min(), y.max(), 5) / roundFact) * roundFact
    ax.set_yticklabels(ylabel[1:4], size=12, rotation=90, va='center')
    ax.set_yticks(ylabel[1:4])
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    roundFact = 10**(np.floor(np.log10(np.abs(x.max() - x.min()))) - 2)
    xlabel = np.round(np.linspace(x.min(), x.max(), 5) / roundFact) * roundFact
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
                  ylabel='Data', linewidth=0.5):
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
            ax.plot(distance[ind], dline[ind], plotStr[ii], linewidth=linewidth)
        else:
            ax.plot(distance[ind], dline[ind], linewidth=linewidth)

    # ax.set_xlim(distance.min(), distance.max())
    ax.set_ylabel(ylabel)
    ax.set_aspect('auto')
    ax.grid(True)
    return ax


def dataHillsideWidget(
    gridObject, EPSGcode=None, HSTransp=0.5, SunAzimuth=270,
    saveAs='./Output/DataHillshade',
    ShapeFileName="./Output/Contours",
    dpi=300, Contours=None,
    scatterData=None, shapeFile=None, omit=[], units='tiltAngle'
  ):

    def plotWidget(
                ColorMap,
                VminVmax,
                Equalize,
                ColorTransp,
                SunAzimuth,
                SunAngle,
                HSTransp,
                vScale,
                saveAs,
                SaveGrid,
                Contours,
                ShapeFileName,
                SaveShape,
                EPSGcode,
         ):

        if SaveGrid:
            lims = gridObject.limits
            DataIO.writeGeotiff(
                gridObject.values, saveAs + '_GRID.tiff',
                gridObject.EPSGcode, lims[0], lims[1],
                lims[2], lims[3], 1,
                dataType='grid')

            if gridObject.EPSGcode != EPSGcode:

                print(
                    "Output EPSG code differ from input grid definition."
                    "The geotiff will be reprojected"
                    )
                DataIO.gdalWarp(
                    saveAs + '_EPSG' + str(int(EPSGcode)) + '_GRID.tiff',
                    saveAs + '_GRID.tiff', int(EPSGcode)
                )
                print(
                    "New file written:" +
                    saveAs + '_EPSG' + str(int(EPSGcode)) + '_GRID.tiff'
                    )

        # Parse contour values
        if Contours is not "":
            vals = re.split(',', Contours)
            cntrs = []
            for val in vals:
                if ":" in val:
                    param = np.asarray(re.split(":", val), dtype='int')
                    cntrs += [np.arange(param[0], param[2], param[1])]

                else:
                    cntrs += [np.float(val)]
            Contours = np.unique(np.sort(np.hstack(cntrs)))
        else:
            Contours = None

        X, Y, data, im, CS = plotSave(
            gridObject, gridObject.values, scatterData, shapeFile,
            SunAzimuth, SunAngle, ColorTransp, HSTransp, vScale, Contours,
            ColorMap, units, VminVmax[0], VminVmax[1], Equalize,
            saveAs, EPSGcode, SaveGrid, dpi=dpi
        )

        if Contours is not "":

            if SaveShape:

                # Export to shapefile
                DataIO.exportShapefile(
                        CS, [],
                        EPSGcode=EPSGcode,
                        saveAs=ShapeFileName,
                        label=units,
                        attType='float'
                )
    # Calculate the original map extents
    assert isinstance(gridObject, DataIO.dataGrid), "Only implemented for objects of class DataIO.dataGrid"
    xLoc = gridObject.hx
    yLoc = gridObject.hy
    data = gridObject.values

    # Trigger the save and uncheck button
    def saveIt(_):

        if SaveGrid.value:
            SaveGrid.value = False
            print('Image saved as: ' + saveAs.value)

    def saveShape(_):

        if SaveShape.value:
            SaveShape.value = False
            print('Shapefile saved as: ' + ShapeFileName.value)

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

    # Contours = widgets.IntSlider(
    #     min=0, max=100, step=2, value=contours, continuous_update=False,
    #     description='Contours'
    #     )

    Contours = widgets.Text(
        value=Contours,
        description='Contours',
        style={'description_width': 'initial'}, continuous_update=False)

    ColorMap = widgets.Dropdown(
        options=cmaps(),
        value='Spectral_r',
        description='ColorMap',
        style={'description_width': 'initial'},
    )

    vmin = data[~np.isnan(data)].min()
    vmax = data[~np.isnan(data)].max()
    VminVmax = widgets.FloatRangeSlider(
                                    value=[vmin, vmax],
                                    min=vmin,
                                    max=vmax,
                                    step=1.0,
                                    description='Color Range',
                                    style={'description_width': 'initial'},
                                    continuous_update=False,
                                    orientation='horizontal',
                                    readout=True,
                                    readout_format='.1f',
                                )

    Equalize = widgets.Dropdown(
                                  options=['Linear', 'HistEqualized'],
                                  value='HistEqualized',
                                  description='Color Normalization',
                                  style={'description_width': 'initial'},
                                )

    saveAs = widgets.Text(
        value=saveAs,
        description='File Name:',
        style={'description_width': 'initial'}
        )
    SaveGrid = widgets.ToggleButton(
        value=False,
        description='Export GeoTiff',
        style={'description_width': 'initial'},
        button_style='',
        tooltip='Description',
        icon='check'
        )

    SaveGrid.observe(saveIt)
    EPSGcode = widgets.FloatText(
        value=gridObject.EPSGcode,
        description='EPSG code:',
        style={'description_width': 'initial'}
    )

    SaveShape = widgets.ToggleButton(
        value=False,
        description='Export Shapefile',
        style={'description_width': 'initial'},
        button_style='',
        tooltip='Description',
        icon='check'
        )

    SaveShape.observe(saveShape)

    ShapeFileName = widgets.Text(
        value=ShapeFileName,
        description='Save as:',
        style={'description_width': 'initial'}
        )

    keys = {
            'ColorMap': ColorMap,
            'VminVmax': VminVmax,
            'Equalize': Equalize,
            'ColorTransp': ColorTransp,
            'SunAzimuth': SunAzimuth,
            'SunAngle': SunAngle,
            'HSTransp': HSTransp,
            'vScale': vScale,
            'saveAs': saveAs,
            'SaveGrid': SaveGrid,
            'Contours': Contours,
            'ShapeFileName': ShapeFileName,
            'SaveShape': SaveShape,
            'EPSGcode': EPSGcode,
        }

    widgList = []
    for key in list(keys.keys()):

        if key not in omit:
            widgList += [keys[key]]

    out = widgets.interactive_output(plotWidget, keys)

    left = widgets.VBox(
            widgList,
            layout=Layout(
                width='35%', height='600px', margin='60px 0px 0px 0px'
            )
        )

    return widgets.HBox([left, out])


def gridFiltersWidget(
    gridObject, gridFilter='derivativeX',
    ColorTransp=0.9, HSTransp=0.5,
    EPSGcode=None, dpi=300, scatterData=None,
    inc=np.nan, dec=np.nan, Contours=None,
    SunAzimuth=270, SunAngle=15, vScale=5.,
    ColorMap='Spectral_r', shapeFile=None,
    saveAs="./Output/MyGeoTiff_" + 'derivativeX',
    ShapeFileName="./Output/Contours_" + 'derivativeX',
    omit=[]
):

    gridProps = [
        'valuesFilledUC', 'valuesFilled',
        'derivativeX', 'derivativeY', 'firstVertical',
        'totalHorizontal', 'tiltAngle', 'analyticSignal',
        'TDXderivative', 'gridFFT', 'gridPadded',
      ]

    def plotWidget(
            Filters, UpDist,
            ColorMap, ColorTransp,
            SunAzimuth, SunAngle,
            HSTransp, vScale,
            saveAs, SaveGrid,
            Contours, ShapeFileName,
            SaveShape, EPSGcode,
         ):

        # If changed upward distance, reset the FFT
        if UpDist != gridObject.heightUC:
            for prop in gridProps:
                    setattr(gridObject, '_{}'.format(prop), None)

            data = gridObject.upwardContinuation(z=UpDist)
            gridObject._gridPadded = None
            gridObject._gridFFT = None

        if Filters == 'TMI':
            data = gridObject.upwardContinuation(z=UpDist)

        else:
            data = getattr(gridObject, '{}'.format(Filters))

        ind = ~np.isnan(data)

        vmin, vmax = np.percentile(data[ind], 5), np.percentile(data[ind], 95)

        # vScale *= (
        #     np.abs(gridObject.values[ind].max() - gridObject.values[ind].min()) *
        #     np.abs(data[ind].max() - data[ind].min())
        # )

        if SaveGrid:
            lims = gridObject.limits
            DataIO.writeGeotiff(
                data,
                saveAs + '_GRID.tiff',
                gridObject.EPSGcode, lims[0], lims[1],
                lims[2], lims[3], 1,
                dataType='grid')

            if gridObject.EPSGcode != EPSGcode:

                print(
                    "Output EPSG code differ from input grid definition."
                    "The geotiff will be reprojected"
                    )
                DataIO.gdalWarp(
                    saveAs + '_EPSG' + str(int(EPSGcode)) + '_GRID.tiff',
                    saveAs + '_GRID.tiff', int(EPSGcode)
                )
                print(
                    "New file written:" +
                    saveAs + '_EPSG' + str(int(EPSGcode)) + '_GRID.tiff'
                    )

                # Parse contour values
        if Contours is not "":
            vals = re.split(',', Contours)
            cntrs = []
            for val in vals:
                if ":" in val:
                    param = np.asarray(re.split(":", val), dtype='int')
                    cntrs += [np.arange(param[0], param[2], param[1])]

                else:
                    cntrs += [np.float(val)]

            Contours = np.unique(np.sort(np.hstack(cntrs)))
        else:
            Contours = None

        X, Y, data, im, CS = plotSave(
            gridObject, data, scatterData, shapeFile,
            SunAzimuth, SunAngle, ColorTransp, HSTransp, vScale, Contours,
            ColorMap, Filters, vmin, vmax, 'HistEqualized',
            saveAs, EPSGcode, SaveGrid, dpi=dpi
        )

        if Contours is not "":

            if SaveShape:

                # Export to shapefile
                DataIO.exportShapefile(
                        CS, [],
                        EPSGcode=EPSGcode,
                        saveAs=ShapeFileName,
                        label=Filters,
                        attType='float'
                )

    assert isinstance(gridObject, DataIO.dataGrid), "Only implemented for objects of class DataIO.dataGrid"

    def saveIt(_):
        if SaveGrid.value:
            SaveGrid.value = False
            print('Image saved as: ' + saveAs.value)

    def saveShape(_):

        if SaveShape.value:
            SaveShape.value = False
            print('Shapefile saved as: ' + ShapeFileName.value)

    def labelUpdate(_):

        saveAs.value = "./Output/MyGeoTiff_" + Filters.value
        ShapeFileName.value = "./Output/Contours_" + Filters.value

    SunAzimuth = widgets.FloatSlider(
        min=0, max=360, step=5, value=SunAzimuth,
        continuous_update=False,
        description='Sun Azimuth'
        )
    SunAngle = widgets.FloatSlider(
        min=0, max=90, step=5, value=SunAngle,
        description='Sun Angle', continuous_update=False
        )
    ColorTransp = widgets.FloatSlider(
        min=0, max=1, step=0.05, value=ColorTransp,
        description='ColorTransp', continuous_update=False
        )
    HSTransp = widgets.FloatSlider(
        min=0, max=1, step=0.05, value=HSTransp,
        description='Sun Transp', continuous_update=False
        )
    vScale = widgets.FloatSlider(
        min=1, max=10, step=2., value=vScale,
        description='V scale', continuous_update=False
        )
    ColorMap = widgets.Dropdown(
        options=cmaps(),
        value=ColorMap,
        description='ColorMap',
        style={'description_width': 'initial'},
        )
    Contours = widgets.Text(
        value=Contours,
        description='Contours',
        style={'description_width': 'initial'}, continuous_update=False)
    Filters = widgets.Dropdown(
        options=[
            'TMI',
            'derivativeX', 'derivativeY', 'firstVertical',
            'totalHorizontal', 'tiltAngle', 'analyticSignal',
        'TDXderivative'],
        value=gridFilter,
        description='Grid Filters',
        style={'description_width': 'initial'},
        )

    Filters.observe(labelUpdate)
    UpDist = widgets.FloatSlider(
        min=0, max=5000, step=10, value=0,
        continuous_update=False, description='UpC Height'
        )
    SaveGrid = widgets.ToggleButton(
        value=False,
        description='Export GeoTiff',
        style={'description_width': 'initial'},
        button_style='',
        tooltip='Description',
        icon='check'
        )

    SaveGrid.observe(saveIt)

    EPSGcode = widgets.FloatText(
        value=gridObject.EPSGcode,
        description='EPSG code:',
        style={'description_width': 'initial'}
    )

    saveAs = widgets.Text(
        value=saveAs,
        description='Save as:',
        style={'description_width': 'initial'}
        )

    SaveShape = widgets.ToggleButton(
        value=False,
        description='Export Shapefile',
        style={'description_width': 'initial'},
        button_style='',
        tooltip='Description',
        icon='check'
        )

    SaveShape.observe(saveShape)

    ShapeFileName = widgets.Text(
        value=ShapeFileName,
        description='Save as:',
        style={'description_width': 'initial'}
        )

    keys = {
        'Filters': Filters,
        'UpDist': UpDist,
        'ColorMap': ColorMap,
        'ColorTransp': ColorTransp,
        'SunAzimuth': SunAzimuth,
        'SunAngle': SunAngle,
        'HSTransp': HSTransp,
        'vScale': vScale,
        'saveAs': saveAs,
        'SaveGrid': SaveGrid,
        'Contours': Contours,
        'ShapeFileName': ShapeFileName,
        'SaveShape': SaveShape,
        'EPSGcode': EPSGcode,
        }

    widgList = []
    for key in list(keys.keys()):

        if key not in omit:
            widgList += [keys[key]]

    out = widgets.interactive_output(plotWidget, keys)

    left = widgets.VBox(
            widgList,
            layout=Layout(
                width='35%', height='600px', margin='60px 0px 0px 0px'
            )
        )

    image = widgets.VBox(
              [out],
              layout=Layout(
                  width='65%', height='600px', margin='0px 0px 0px 0px'
              )
            )

    return widgets.HBox([left, out])

    # return out


def gridTilt2Depth(
    gridObject, gridFilter='tiltAngle',
    ColorTransp=0.9, HSTransp=0.5,
    EPSGcode=None, dpi=300, scatterData=None,
    SunAzimuth=270, SunAngle=15, vScale=5., shapeFile=None,
    ColorMap='Spectral_r', ColorDepth='viridis_r', depthRange=[0, 500],
    markerSize=1, omit=[],
    ShapeFileName="./Output/EstimatedDepth",
    GridFileName="./Output/MyGeoTiff",
    CSVFileName="./Output/EstimatedDepth.csv",
):

    gridProps = [
        'valuesFilledUC', 'valuesFilled',
        'derivativeX', 'derivativeY', 'firstVertical',
        'totalHorizontal', 'tiltAngle', 'analyticSignal',
        'TDXderivative','RTP', 'gridFFT', 'gridPadded',
      ]

    def plotWidget(
            SunAzimuth, SunAngle,
            ColorTransp, HSTransp, vScale,
            ColorMap, Filters, UpDist,
            ContourColor, ContourSize,
            GridFileName, EPSGcode, SaveGrid,
            ShapeFileName, SaveShape,
            CSVFileName, SaveCSV

         ):

        # If changed upward distance, reset the FFT
        if UpDist != gridObject.heightUC:
            for prop in gridProps:
                    setattr(gridObject, '_{}'.format(prop), None)

            data = gridObject.upwardContinuation(z=UpDist)
            gridObject._gridPadded = None
            gridObject._gridFFT = None

        if Filters == 'TMI':
            data = gridObject.upwardContinuation(z=UpDist)
        else:
            data = getattr(gridObject, '{}'.format(Filters))

        ind = ~np.isnan(data)
        vmin, vmax = np.percentile(data[ind], 5), np.percentile(data[ind], 95)

        # Compute estimated depth
        polylines, attributes = MathUtils.estimateDepth(gridObject)

        if SaveShape:

            if len(polylines) > 0:
                # Export to shapefile
                DataIO.exportShapefile(
                    polylines, attributes,
                    EPSGcode=EPSGcode, saveAs=ShapeFileName,
                    label='AvgDepth')
                print('Shapefile saved as: ' + ShapeFileName)


            else:
                print("No [-45, 45] contour found")

        if SaveCSV:
            if len(polylines) > 0:
                # Export to CSV file
                np.savetxt(CSVFileName, np.c_[np.vstack(polylines)[:, 0:2], np.concatenate(attributes)], fmt="%.2f", delimiter=",")
                print('CSV saved as: ' + CSVFileName)
            else:
                print("No [-45, 45] contour found")

        if len(polylines) > 0:
            scatterData = {}
            scatterData['x'] = np.vstack(polylines)[:, 0]
            scatterData['y'] = np.vstack(polylines)[:, 1]
            scatterData['size'] = ContourSize
            scatterData['c'] = np.concatenate(attributes)
            scatterData['cmap'] = ContourColor
            scatterData['clim'] = [
                np.percentile(scatterData['c'], 25),
                np.percentile(scatterData['c'], 75)
            ]
            scatterData['colorbar'] = True
        else:
            scatterData = None

        vScale *= (
            np.abs(gridObject.values[ind].max() - gridObject.values[ind].min()) *
            np.abs(data[ind].max() - data[ind].min())
        )
        plotSave(
            gridObject, data, scatterData, shapeFile,
            SunAzimuth, SunAngle,
            ColorTransp, HSTransp, vScale, None,
            ColorMap, Filters, vmin, vmax, 'HistEqualized',
            GridFileName, EPSGcode, SaveGrid,
            dpi=dpi
        )


    assert isinstance(gridObject, DataIO.dataGrid), "Only implemented for objects of class DataIO.dataGrid"

    def saveIt(_):

        if SaveGrid.value:
            SaveGrid.value = False
            print('Image saved as: ' + GridFileName.value)

    def saveCSV(_):

        if SaveCSV.value:
            SaveCSV.value = False


    def saveShape(_):

        if SaveShape.value:
            SaveShape.value = False


    SunAzimuth = widgets.FloatSlider(
        min=0, max=360, step=5, value=SunAzimuth,
        continuous_update=False,
        description='Sun Azimuth'
        )
    SunAngle = widgets.FloatSlider(
        min=0, max=90, step=5, value=SunAngle,
        description='Sun Angle', continuous_update=False
        )
    ColorTransp = widgets.FloatSlider(
        min=0, max=1, step=0.05, value=ColorTransp,
        description='ColorTransp', continuous_update=False
        )
    HSTransp = widgets.FloatSlider(
        min=0, max=1, step=0.05, value=HSTransp,
        description='Sun Transp', continuous_update=False
        )
    vScale = widgets.FloatSlider(
        min=1, max=200, step=5., value=vScale,
        description='V scale', continuous_update=False
        )
    ColorMap = widgets.Dropdown(
        options=cmaps(),
        value=ColorMap,
        description='ColorMap',
        style={'description_width': 'initial'},
        )
    Filters = widgets.Dropdown(
        options=[
            'TMI',
            'tiltAngle'],
        value=gridFilter,
        description='Grid Filters',
        style={'description_width': 'initial'},
        )
    UpDist = widgets.FloatSlider(
        min=0, max=2000, step=10, value=0,
        continuous_update=False, description='UpC Height'
        )
    SaveGrid = widgets.ToggleButton(
        value=False,
        description='Export GeoTiff',
        style={'description_width': 'initial'},
        button_style='',
        tooltip='Export GeoTiff image',
        icon='check'
        )

    SaveGrid.observe(saveIt)
    SaveShape = widgets.ToggleButton(
        value=False,
        description='Export Shapefile',
        style={'description_width': 'initial'},
        button_style='',
        tooltip='Export Contours',
        icon='check'
        )

    SaveShape.observe(saveShape)
    SaveCSV = widgets.ToggleButton(
        value=False,
        description='Export CSV',
        style={'description_width': 'initial'},
        button_style='',
        tooltip='Export estimated depth points',
        icon='check'
        )

    SaveCSV.observe(saveCSV)

    GridFileName = widgets.Text(
        value=GridFileName,
        description='Save as:',
        style={'description_width': 'initial'}
        )
    ShapeFileName = widgets.Text(
        value=ShapeFileName,
        description='Shapefile name:',
        style={'description_width': 'initial'}
        )

    CSVFileName = widgets.Text(
        value=CSVFileName,
        description='CSV file name:',
        style={'description_width': 'initial'}
        )
    ContourColor = widgets.Dropdown(
        options=cmaps(),
        value='viridis',
        description='Depth Color',
        style={'description_width': 'initial'},
        )
    ContourSize = widgets.FloatSlider(
        min=0.1, max=10, step=0.1, value=1,
        description='Marker Size', continuous_update=False
        )
    EPSGcode = widgets.FloatText(
        value=gridObject.EPSGcode,
        description='EPSG code:',
        style={'description_width': 'initial'}
    )
    keys = {
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
        "EPSGcode": EPSGcode,
        'SaveGrid': SaveGrid,
        'ShapeFileName': ShapeFileName,
        'SaveShape': SaveShape,
        'CSVFileName': CSVFileName,
        'SaveCSV': SaveCSV
        }

    widgList = []
    for key in list(keys.keys()):

        if key not in omit:
            widgList += [keys[key]]

    out = widgets.interactive_output(plotWidget, keys)

    left = widgets.VBox(
            widgList,
            layout=Layout(
                width='35%', height='600px', margin='60px 0px 0px 0px'
            )
        )

    return widgets.HBox([left, out])


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

        locs = np.c_[Xloc.flatten(order='F'), Yloc.flatten(order='F'), Zloc.flatten(order='F')]
        survey, _, _ = ProblemSetter.setSyntheticProblem(locs, EarthField=dataVals[-3:])

        xyz = survey.rxLoc
        plt.figure(figsize=(10, 8))
        ax1 = plt.subplot(1, 2, 1)
        fig, im, cbar = plotData2D(
          xyz[:, 0], xyz[:, 1], survey.dobs,
          ax=ax1, cmap='Spectral_r', marker=False, colorbar=False,
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
                                style={'description_width': 'initial'},
                                )
                              )

    return out


def dataGriddingWidget(
    survey, EPSGcode=np.nan, saveAs="Output/MyGeoTiff", marker=True,
    shapeFile=None, inc=np.nan, dec=np.nan, dataColumn=-1, overlap=0,
    Method='minimumCurvature', Contours=None, omit=[], units="TMI",
    dpi=200, resolution=25, maxDistance=200, nPoints=5,
):

    def plotWidget(
            Resolution, MaxDistance, Method,
            ColorMap,
            EPSGcode,
            GetIncDec, saveAs, SaveGrid
         ):

        if Method == 'minimumCurvature':
            gridCC, d_grid = MathUtils.minCurvatureInterp(
                np.c_[xLoc, yLoc], data, maxDistance=MaxDistance,
                gridSize=Resolution, method='spline', overlap=overlap,
                nPoints=nPoints,
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

            tree = cKDTree(np.c_[xLoc, yLoc])
            # xi = _ndim_coords_from_arrays((gridCC[:,0], gridCC[:,1]), ndim=2)
            dists, _ = tree.query(
                np.c_[X.flatten(order='F'), Y.flatten(order='F')]
            )

            # Copy original result but mask missing values with NaNs
            d_grid[(dists > MaxDistance).reshape(d_grid.shape, order='F')] = np.nan

        gridObject = DataIO.dataGrid()

        gridObject._values = d_grid
        gridObject.nx, gridObject.ny = gridObject.values.shape[1], gridObject.values.shape[0]

        gridObject.dx = (X.max() - X.min()) / (gridObject.values.shape[1] - 1)
        gridObject.dy = (Y.max() - Y.min()) / (gridObject.values.shape[0] - 1)

        gridObject.x0, gridObject.y0 = X.min()-gridObject.dx/2., Y.min()-gridObject.dy/2.

        gridObject.limits = np.r_[gridObject.x0, gridObject.x0+gridObject.nx*gridObject.dx, gridObject.y0, gridObject.y0+gridObject.ny*gridObject.dy]

        if not np.isnan(EPSGcode):
            gridObject.EPSGcode = int(EPSGcode)

        if SaveGrid:
            if np.isnan(EPSGcode):
                print("Need to assign an EPSGcode before exporting")
                return
            DataIO.writeGeotiff(
                d_grid, saveAs + '.tiff',
                gridObject.EPSGcode,
                gridObject.limits[0], gridObject.limits[1],
                gridObject.limits[2], gridObject.limits[3], 1,
                dataType='grid')

        if marker:
            scatterData = {}
            scatterData['x'] = xLoc
            scatterData['y'] = yLoc
            scatterData['size'] = 10
            scatterData['c'] = 'k'
            scatterData['cmap'] = 'Greys'
            scatterData['clim'] = [None, None]
            scatterData['colorbar'] = False
        else:
            scatterData = None

        plotSave(
            gridObject, d_grid, scatterData, None,
            90, 15, 1, 0, 0, None,
            ColorMap, units, None, None, 'HistEqualized',
            saveAs, EPSGcode, SaveGrid, dpi=dpi
        )
        # Create grid object
        return gridObject

    # Calculate the original map extents
    xLoc = survey[:, 0]
    yLoc = survey[:, 1]
    data = survey[:, int(dataColumn-1)]

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

    def saveIt(_):

        if SaveGrid.value:
            SaveGrid.value = False
            print('Image saved as: ' + saveAs.value)

    Resolution = widgets.FloatText(
        value=resolution,
        description='Grid (m):',
        style={'description_width': 'initial'}
        )

    MaxDistance = widgets.FloatText(
        value=maxDistance,
        description='Dist Max (m):',
        style={'description_width': 'initial'},
        )

    Method = widgets.Dropdown(
        options=[
          'nearest', 'linear', 'cubic',
          'minimumCurvature'
          ],
        value=Method,
        description='Method',
        style={'description_width': 'initial'},
        )
    ColorMap = widgets.Dropdown(
        options=cmaps(),
        value='Spectral_r',
        description='ColorMap',
        style={'description_width': 'initial'},
        )
    EPSGcode = widgets.FloatText(
        value=EPSGcode,
        description='EPSG code:',
        style={'description_width': 'initial'}
    )

    GetIncDec = widgets.ToggleButton(
        value=False,
        description='Fetch Inc/Dec',
        style={'description_width': 'initial'},
        button_style='',
        tooltip='Connect to NOAA API',
        icon='check'
        )

    GetIncDec.observe(fetchURL)
    saveAs = widgets.Text(
        value=saveAs,
        description='Save as:',
        style={'description_width': 'initial'}
        )
    SaveGrid = widgets.ToggleButton(
        value=False,
        description='Export Grid',
        style={'description_width': 'initial'},
        button_style='',
        tooltip='Write file',
        icon='check'
        )

    SaveGrid.observe(saveIt)

    for key in omit:
        locals()[key].disabled=True

    out = widgets.interactive(plotWidget,
                              Resolution=Resolution,
                              MaxDistance=MaxDistance,
                              Method=Method,
                              ColorMap=ColorMap,
                              EPSGcode=EPSGcode,
                              GetIncDec=GetIncDec,
                              saveAs=saveAs,
                              SaveGrid=SaveGrid
                              )

    return out


def dataGridGeoref(
    gridObject, EPSGcode=np.nan, saveAs="./Output/MyGeoTiff",
    shapeFile=None, inc=np.nan, dec=np.nan, ndv=np.nan, applyRTP=False,
    omit=[]
):

    values = gridObject.values.copy()
    def plotWidget(
            ColorMap,
            EPSGcode, inc, dec, ndv,
            GetIncDec, applyRTP, units="TMI"
         ):

        gridObject._values = values

        if ~np.isnan(ndv):
            gridObject._values[values==ndv] = np.nan

        if not np.isnan(EPSGcode):
            gridObject.EPSGcode = int(EPSGcode)

        gridObject.inc, gridObject.dec = inc, dec

        gridObject.setRTP(applyRTP)

        plotSave(
            gridObject, gridObject.values, None, None,
            90, 15, 1, 0, 1, None,
            "Spectral_r", units, None, None, 'HistEqualized', "", EPSGcode,
            False, dpi=200
        )

        # Create grid object
        return gridObject

    assert isinstance(gridObject, DataIO.dataGrid), "Only implemented for objects of class DataIO.dataGrid"

    def fetchURL(_):
        if GetIncDec.value:
            GetIncDec.value = False
            if np.isnan(EPSGcode.value):
                print("Enter EPSGcode first")
                return

            x, y, z = np.mean(gridObject.hx), np.mean(gridObject.hy), 0.
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

    def saveIt(_):

        if SaveGrid.value:
            SaveGrid.value = False
            print('Image saved as: ' + saveAs.value)

    ColorMap = widgets.Dropdown(
        options=cmaps(),
        value='Spectral_r',
        description='ColorMap',
        style={'description_width': 'initial'},
        )

    if gridObject.EPSGcode is not None:
        EPSGcode = gridObject.EPSGcode

    EPSGcode = widgets.FloatText(
        value=EPSGcode,
        description='EPSG code:',
        style={'description_width': 'initial'}
    )
    inc = widgets.FloatText(
        value=inc,
        description='Inclination angle',
        tooltip='Positive downward from horizontal:',
        style={'description_width': 'initial'}
        )
    dec = widgets.FloatText(
        value=dec,
        description='Declination angle',
        tooltip='Positve clockwise from North:',
        style={'description_width': 'initial'}
        )

    if np.any(gridObject.values > 1e+10):

        ndv = np.unique(gridObject.values[gridObject.values > 1e+10])

    elif np.any(gridObject.values == -99999):

        ndv = -99999

    elif np.any(np.isnan(gridObject.values)):

        ndv = np.nan

    else:

        ndv = np.nan

    ndv = widgets.FloatText(
        value=ndv,
        description='No-data-Value',
        style={'description_width': 'initial'}
        )

    GetIncDec = widgets.ToggleButton(
        value=False,
        description='Fetch Inc/Dec',
        style={'description_width': 'initial'},
        button_style='',
        tooltip='Connect to NOAA API',
        icon='check'
        )

    GetIncDec.observe(fetchURL)

    applyRTP = widgets.ToggleButton(
        value=applyRTP,
        description='Reduce to pole',
        style={'description_width': 'initial'},
        button_style='',
        tooltip='Transform to RTP data',
        icon='check'
        )

    saveAs = widgets.Text(
        value=saveAs,
        description='Save as:',
        style={'description_width': 'initial'}
        )
    SaveGrid = widgets.ToggleButton(
        value=False,
        description='Export GeoTiff',
        style={'description_width': 'initial'},
        button_style='',
        tooltip='Write file',
        icon='check'
        )
    SaveGrid.observe(saveIt)

    for key in omit:
        locals()[key].disabled=True
    out = widgets.interactive(plotWidget,
                              ColorMap=ColorMap,
                              EPSGcode=EPSGcode,
                              inc=inc, dec=dec,
                              ndv=ndv,
                              GetIncDec=GetIncDec,
                              applyRTP=applyRTP,
                              saveAs=saveAs,
                              SaveGrid=SaveGrid
                              )

    return out


def setDataExtentWidget(
    survey, East=None, North=None, nCx=100, nCy=100,
    EPSGcode=None, saveAs="./Output/MyGeoTiff", omit=[]
):
    """
        Small application to carve out a subset of a larger data set
    """

    def dataSelector(East, North, Width, Height, saveAs, EPSGcode, SaveGrid):

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
                linewidth=3,
                linestyle='--',
                zorder=3
                )
            )

        # Extract data within window and plot
        indx = np.logical_and(xLoc > lims[0], xLoc < lims[1])

        indy = np.logical_and(yLoc > lims[2], yLoc < lims[3])

        nx, ny = np.count_nonzero(indx), np.count_nonzero(indy)

        # Create new dataGrid object
        dataSub = DataIO.dataGrid()

        # coordinate_system = grid.coordinate_system
        temp = survey.values[:, indx]
        temp = temp[indy, :]


        dataSub._values = temp
        dataSub.nx, dataSub.ny = nx, ny
        dataSub.dx, dataSub.dy = survey.dx, survey.dy
        dataSub.x0, dataSub.y0 = xLoc[indx].min()-survey.dx/2, yLoc[indy].min()-survey.dy/2
        dataSub.EPSGcode = survey.EPSGcode

        dataSub.limits = np.r_[
            dataSub.x0, dataSub.x0+dataSub.nx*dataSub.dx,
            dataSub.y0, dataSub.y0+dataSub.ny*dataSub.dy
        ]

        if SaveGrid:

            if np.isnan(EPSGcode) and getattr(dataSub, 'EPSGcode', None) is None:
                print("Need to assign an EPSGcode before exporting")
                return

            elif getattr(dataSub, 'EPSGcode', None) is None:
                dataSub.EPSGcode = int(EPSGcode)

            DataIO.writeGeotiff(
                dataSub.values, saveAs + '.tiff',
                dataSub.EPSGcode,
                dataSub.limits[0], dataSub.limits[1],
                dataSub.limits[2], dataSub.limits[3], 1,
                dataType='grid')

            if dataSub.EPSGcode != EPSGcode:

                print(
                    "Output EPSG code differ from input grid definition."
                    "The geotiff will be reprojected"
                    )
                DataIO.gdalWarp(
                    saveAs + 'EPSG' + str(EPSGcode) + '.tiff',
                    saveAs + '.tiff', int(EPSGcode)
                )
                print(
                    "New file written:" +
                    saveAs + 'EPSG' + str(int(EPSGcode)) + '.tiff'
                    )


        # fig,
        axs = plt.subplot(1, 2, 2)
        fig, im, cbar = plotData2D(
            xLoc[indx], yLoc[indy], dataSub.values, marker=False, fig=fig, ax=axs
        )
        return dataSub

    if isinstance(survey, DataIO.dataGrid):

        xLoc = survey.hx
        yLoc = survey.hy
        xlim = survey.limits[:2]
        ylim = survey.limits[2:]
        dx = survey.dx
        dy = survey.dy
        nx = survey.nx
        ny = survey.ny


        if East is None:
            East = np.mean(xLoc)

        if North is None:
            North = np.mean(yLoc)

    else:
        print("Only implemented for class 'dataGrid'")

    def saveIt(_):

        if SaveGrid.value:
            SaveGrid.value = False
            print('Image saved as: ' + saveAs.value)

    East = widgets.FloatSlider(
        min=xlim[0], max=xlim[1], step=dx, value=East, continuous_update=False
        )
    North = widgets.FloatSlider(
        min=ylim[0], max=ylim[1], step=dy, value=North, continuous_update=False
        )
    Width = widgets.FloatSlider(
        min=2*dx, max=np.abs(xlim[1] - xlim[0])+2*dx, step=dx, value=nCx*dx,
        continuous_update=False
        )
    Height = widgets.FloatSlider(
        min=2*dy, max=np.abs(ylim[1] - ylim[0])+2*dy, step=dy, value=nCy*dy,
        continuous_update=False
        )
    saveAs = widgets.Text(
        value=saveAs,
        description='Save as:',
        style={'description_width': 'initial'}
        )

    if survey.EPSGcode is not None:
        EPSGcode = survey.EPSGcode

    EPSGcode = widgets.FloatText(
        value=EPSGcode,
        description='EPSG code:',
        style={'description_width': 'initial'}
    )
    SaveGrid = widgets.ToggleButton(
        value=False,
        description='Export Grid',
        button_style='',
        tooltip='Write file',
        icon='check'
        )
    SaveGrid.observe(saveIt)

    for key in omit:
        locals()[key].disabled=True
    out = widgets.interactive(
            dataSelector,
            East=East,
            North=North,
            Width=Width,
            Height=Height,
            saveAs=saveAs,
            EPSGcode=EPSGcode,
            SaveGrid=SaveGrid
            )

    return out


def plotSave(
        gridObject, data, scatterData, shapeFile,
        SunAzimuth, SunAngle, ColorTransp, HSTransp, vScale, Contours,
        ColorMap, Filters, vmin, vmax, equalizeHist, saveAs, EPSGcode,
        SaveGrid, dpi=200
     ):

    if SaveGrid:
        fig = plt.figure()
        fig.set_size_inches(9, 9)
        axs = plt.Axes(fig, [0., 0., 1., 1.])
        axs.set_axis_off()
        fig.add_axes(axs)

    else:

        fig = plt.figure(figsize=(9, 9))
        axs = plt.subplot()

    # Add shading
    X, Y, data, im, CS = plotDataHillside(
        gridObject.hx, gridObject.hy, data,
        axs=axs, cmap=ColorMap,
        clabel=False, resolution=10,
        vmin=vmin, vmax=vmax, contours=Contours,
        alpha=ColorTransp, alphaHS=HSTransp,
        ve=vScale, azdeg=SunAzimuth, altdeg=SunAngle,
        equalizeHist=equalizeHist, scatterData=scatterData,
        shapeFile=shapeFile
    )

    if SaveGrid:

        if saveAs is None:
            saveAs = Filters

        plt.savefig(saveAs + '.png', dpi=dpi)
        plt.close()

        img = np.asarray(PIL.Image.open(saveAs + '.png'))

        if (EPSGcode is None) and (getattr(gridObject, 'EPSGcode', None) is None):
            print("Need to assign an EPSGcode before exporting")
            return

        elif getattr(gridObject, 'EPSGcode', None) is None:
            gridObject.EPSGcode = int(EPSGcode)

        DataIO.writeGeotiff(
            np.flipud(img), saveAs + '.tiff',
            gridObject.EPSGcode,
            gridObject.limits[0], gridObject.limits[1],
            gridObject.limits[2], gridObject.limits[3], 3
        )

        if gridObject.EPSGcode != EPSGcode:

            print(
                "Output EPSG code differ from input grid definition."
                "The geotiff will be reprojected"
                )
            DataIO.gdalWarp(
                saveAs + '_EPSG' + str(int(EPSGcode)) + '.tiff',
                saveAs + '.tiff', int(EPSGcode)
            )
            print(
                "New file written:" +
                saveAs + '_EPSG' + str(int(EPSGcode)) + '.tiff'
                )

        os.remove(saveAs + '.png')

        fig, ax = plt.figure(), plt.subplot()
        plt.gca().set_visible(False)
        cbar = plt.colorbar(im, fraction=0.02)
        plt.savefig(saveAs + 'Colorbar.png', dpi=dpi, bbox_inches='tight')

    else:
        # Add points at the gridObject locations
        # plt.scatter(xLoc, yLoc, s=2, c='k')



        axs.set_aspect('equal')
        cbar = plt.colorbar(im, fraction=0.02)

        axs.set_xlabel("Easting (m)", size=14)
        axs.set_ylabel("Northing (m)", size=14)
        axs.grid('on', color='k', linestyle='--')

        if (scatterData is not None) and scatterData['colorbar']:
            pos = axs.get_position()
            cbarax = fig.add_axes([pos.x0+0.875, pos.y0+0.225,  pos.width*.025, pos.height*0.4])
            norm = mpl.colors.Normalize(vmin=scatterData['clim'][0], vmax=scatterData['clim'][1])
            cmap = mpl.cm.get_cmap(name=scatterData['cmap'])
            cb = mpl.colorbar.ColorbarBase(
              cbarax, cmap=cmap,
              norm=norm,
              orientation="vertical")
            cb.set_label("Depth (m)", size=12)

    return X, Y, data, im, CS
