from . import Mag
from . import MagUtils
import SimPEG.PF as PF
from SimPEG.Utils import mkvc
from scipy.constants import mu_0
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import ipywidgets as widgets
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.interpolate import griddata, interp1d
import scipy.sparse as sp
from scipy.spatial import cKDTree
from scipy.interpolate.interpnd import _ndim_coords_from_arrays
from matplotlib.colors import LightSource, Normalize
from scipy.sparse.linalg import bicgstab

def PFSimulator(prism, survey):

    def PFInteract(update, susc, comp, irt, Q, RemInc, RemDec,
                   Profile_npt, Profile_azm, Profile_len,
                   Profile_ctx, Profile_cty):

        # # Get the line extent from the 2D survey for now
        prob = Mag.problem()
        prob.prism = prism.result
        prob.survey = survey.result

        return PlotFwrSim(prob, susc, comp, irt, Q, RemInc, RemDec,
                          Profile_azm, Profile_len, Profile_npt,
                          Profile_ctx, Profile_cty)

    locs = survey.result.srcField.rxList[0].locs
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

    # Create problem

def PlotFwrSim(prob, susc, comp, irt, Q, rinc, rdec,
               Profile_azm, Profile_len, Profile_npt,
               Profile_ctx, Profile_cty):

    def MagSurvey2D(survey, Profile_ctx, Profile_cty, Profile_azm,
                    Profile_len, Profile_npt,
                    data=None, fig=None, ax=None,
                    vmin=None, vmax=None, pred=None):

        # Get the line extent from the 2D survey for now
        Profile_azm /= 180./np.pi
        Profile_len /= 2.*0.98

        dx = np.cos(-Profile_azm)*Profile_lenf
        dy = np.sin(-Profile_azm)*Profile_len

        a = [Profile_ctx - dx, Profile_cty - dy]
        b = [Profile_ctx + dx, Profile_cty + dy]

        return plotMagSurvey2D(survey, a, b, Profile_npt,
                               data=data, fig=fig, ax=ax,
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

        return plotProfile(xyz, dobs, a, b, Profile_npt,
                           data=data, fig=fig, ax=ax)

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

    f = plt.figure(figsize=(8, 8))

    ax0 = plt.subplot(1, 2, 1)
    MagSurvey2D(survey, Profile_ctx, Profile_cty, Profile_azm,
                Profile_len, Profile_npt, fig=f, ax=ax0, pred=dpred,
                vmin=survey.dobs.min(), vmax=survey.dobs.max())

    f = plt.figure(figsize=(12, 5))
    ax2 = plt.subplot()
    MagSurveyProfile(survey, Profile_ctx, Profile_cty, Profile_azm,
                     Profile_len, Profile_npt, data=dpred, fig=f, ax=ax2)

    plt.show()


def ViewMagSurvey2D(survey):

    def MagSurvey2D(East, North, Width, Height, Azimuth, Length, Npts, Profile):

        # Get the line extent from the 2D survey for now
        Azimuth /= 180./np.pi
        Length /= 2.*0.98
        a = [East - np.cos(-Azimuth)*Length, North - np.sin(-Azimuth)*Length]
        b = [East + np.cos(-Azimuth)*Length, North + np.sin(-Azimuth)*Length]

        xlim = East + np.asarray([-Width/2., Width/2.])
        ylim = North + np.asarray([-Height/2., Height/2.])

        # Re-sample the survey within the region
        rxLoc = survey.srcField.rxList[0].locs

        ind = np.all([rxLoc[:, 0] > xlim[0], rxLoc[:, 0] < xlim[1],
                      rxLoc[:, 1] > ylim[0], rxLoc[:, 1] < ylim[1]], axis=0)

        rxLoc = PF.BaseMag.RxObs(rxLoc[ind, :])
        srcField = PF.BaseMag.SrcField([rxLoc], param=survey.srcField.param)
        surveySim = PF.BaseMag.LinearSurvey(srcField)
        surveySim.dobs = survey.dobs[ind]

        fig = plt.figure(figsize=(6, 9))
        ax1 = plt.subplot(2, 1, 1)

        plotMagSurvey2D(surveySim, a, b, Npts, fig=fig, ax=ax1)

        if Profile:

            ax2 = plt.subplot(2, 1, 2)

            xyz = surveySim.srcField.rxList[0].locs
            dobs = surveySim.dobs
            plotProfile(xyz, dobs, a, b, Npts, data=None,
                        fig=fig, ax=ax2)

        return surveySim

    # Calculate the original map extents
    locs = survey.srcField.rxList[0].locs
    xlim = np.asarray([locs[:, 0].min(), locs[:, 0].max()])
    ylim = np.asarray([locs[:, 1].min(), locs[:, 1].max()])

    Lx = xlim[1] - xlim[0]
    Ly = ylim[1] - ylim[0]
    diag = (Lx**2. + Ly**2.)**0.5 /2.

    East = np.mean(xlim)
    North = np.mean(ylim)
    cntr = [East, North]

    out = widgets.interactive(MagSurvey2D,
                    East=widgets.FloatSlider(min=cntr[0]-Lx, max=cntr[0]+Lx, step=10, value=cntr[0],continuous_update=False),
                    North=widgets.FloatSlider(min=cntr[1]-Ly, max=cntr[1]+Ly, step=10, value=cntr[1],continuous_update=False),
                    Width=widgets.FloatSlider(min=10, max=Lx*1.05, step=10, value=Lx*1.05, continuous_update=False),
                    Height=widgets.FloatSlider(min=10, max=Ly*1.05, step=10, value=Ly*1.05, continuous_update=False),
                    Azimuth=widgets.FloatSlider(min=-90, max=90, step=5, value=0, continuous_update=False),
                    Length=widgets.FloatSlider(min=10, max=diag, step=10, value= Ly, continuous_update=False),
                    Npts=widgets.BoundedFloatText(min=10, max=100, step=1, value=20, continuous_update=False),
                    Profile=widgets.ToggleButton(description='Profile', value=False))

    return out


def plotMagSurvey2D(survey, a, b, npts, data=None, pred=None,
                    fig=None, ax=None, vmin=None, vmax=None):
    """
    Plot the data and line profile inside the spcified limits
    """

    if fig is None:
        fig = plt.figure()

    if ax is None:
        ax = plt.subplot(1, 2, 1)

    x, y = linefun(a[0], b[0], a[1], b[1], npts)
    rxLoc = survey.srcField.rxList[0].locs

    if data is None:
        data = survey.dobs

    # Use SimPEG.PF ploting function
    PF.Magnetics.plot_obs_2D(rxLoc, d=data, fig=fig,  ax=ax,
                             vmin=vmin, vmax=vmax,
                             marker=False, cmap='RdBu_r')

    ax.plot(x, y, 'w.', ms=10)
    ax.text(x[0], y[0], 'A', fontsize=16, color='w', ha='left')
    ax.text(x[-1], y[-1], 'B', fontsize=16,
            color='w', ha='right')
    ax.grid(True)

    if pred is not None:
        ax2 = plt.subplot(1, 2, 2)

        if pred.min() != pred.max():
            PF.Magnetics.plot_obs_2D(rxLoc, d=pred, fig=fig,  ax=ax2,
                                     vmin=vmin, vmax=vmax,
                                     marker=False, cmap='RdBu_r')

        else:
            PF.Magnetics.plot_obs_2D(rxLoc, d=pred, fig=fig,  ax=ax2,
                                     vmin=pred.min(), vmax=pred.max(),
                                     marker=False, cmap='RdBu_r')
        ax2.plot(x, y, 'w.', ms=10)
        ax2.text(x[0], y[0], 'A', fontsize=16, color='w',
                ha='left')
        ax2.text(x[-1], y[-1], 'B', fontsize=16,
                color='w', ha='right')
        ax2.set_yticks([])
        ax2.set_yticklabels("")
        ax2.grid(True)

    plt.show()
    return


def plotProfile(xyz, dobs, a, b, npts, data=None,
                fig=None, ax=None, dType='3D'):
    """
    Plot the data and line profile inside the spcified limits
    """

    if fig is None:
        fig = plt.figure(figsize=(6, 9))

        plt.rcParams.update({'font.size': 14})

    if ax is None:
        ax = plt.subplot()

    rxLoc = xyz

    x, y = linefun(a[0], b[0], a[1], b[1], npts)

    distance = np.sqrt((x-a[0])**2.+(y-a[1])**2.)

    if dType == '2D':
        distance = rxLoc[:, 0]
        dline = dobs

    else:
        dline = griddata(rxLoc[:, :2], dobs, (x, y), method='linear')

    ax.plot(distance, dline, 'b.-')

    if data is not None:

        if dType == '2D':
            distance = rxLoc[:, 0]
            dline = data

        else:
            dline = griddata(rxLoc[:, :2], data, (x, y), method='linear')

        ax.plot(distance, dline, 'r.-')

    ax.set_xlim(distance.min(), distance.max())

    ax.set_xlabel("Distance (m)")
    ax.set_ylabel("Magnetic field (nT)")

    #ax.text(distance.min(), dline.max()*0.8, 'A', fontsize = 16)
    # ax.text(distance.max()*0.97, out_linei.max()*0.8, 'B', fontsize = 16)
    ax.legend(("survey", "simulated"), bbox_to_anchor=(0.5, -0.3))
    ax.grid(True)
    plt.show()

    return True


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
                              View_azm=widgets.FloatSlider(min=0, max=360, step=1, value=220, continuous_update=False),
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
        depth = prism.z0
        x1, x2 = prism.xn[0]-prism.xc, prism.xn[1]-prism.xc
        y1, y2 = prism.yn[0]-prism.yc, prism.yn[1]-prism.yc
        z1, z2 = prism.zn[0]-prism.zc, prism.zn[1]-prism.zc
        pinc, pdec = prism.pinc, prism.pdec

        # Create a rectangular prism, rotate and plot
        block_xyz = np.asarray([[x1, x1, x2, x2, x1, x1, x2, x2],
                               [y1, y2, y2, y1, y1, y2, y2, y1],
                               [z1, z1, z1, z1, z2, z2, z2, z2]])

        R = MagUtils.rotationMatrix(pinc, pdec)

        xyz = R.dot(block_xyz).T

        # Offset the prism to true coordinate
        offx = prism.xc
        offy = prism.yc
        offz = prism.zc

        #print xyz
        # Face 1
        axs.add_collection3d(Poly3DCollection([list(zip(xyz[:4, 0] + offx,
                                                   xyz[:4, 1] + offy,
                                                   xyz[:4, 2] + offz))]))

        # Face 2
        axs.add_collection3d(Poly3DCollection([list(zip(xyz[4:, 0] + offx,
                                                   xyz[4:, 1] + offy,
                                                   xyz[4:, 2] + offz))], facecolors=color))

        # Face 3
        axs.add_collection3d(Poly3DCollection([list(zip(xyz[[0, 1, 5, 4], 0] + offx,
                                                   xyz[[0, 1, 5, 4], 1] + offy,
                                                   xyz[[0, 1, 5, 4], 2] + offz))]))

        # Face 4
        axs.add_collection3d(Poly3DCollection([list(zip(xyz[[3, 2, 6, 7], 0] + offx,
                                                   xyz[[3, 2, 6, 7], 1] + offy,
                                                   xyz[[3, 2, 6, 7], 2] + offz))]))

       # Face 5
        axs.add_collection3d(Poly3DCollection([list(zip(xyz[[0, 4, 7, 3], 0] + offx,
                                                   xyz[[0, 4, 7, 3], 1] + offy,
                                                   xyz[[0, 4, 7, 3], 2] + offz))]))

       # Face 6
        axs.add_collection3d(Poly3DCollection([list(zip(xyz[[1, 5, 6, 2], 0] + offx,
                                                   xyz[[1, 5, 6, 2], 1] + offy,
                                                   xyz[[1, 5, 6, 2], 2] + offz))]))


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
    print(azmDeg)
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
        return plotProfile(xyzLoc, survey2D.dobs, a, b, 10, data=dpred, dType='2D')

    Q = widgets.interactive(profiledata, Binc=widgets.FloatSlider(min=-90., max=90, step=5, value=90, continuous_update=False),
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

def minCurvatureInterp(
    locs, data,
    vectorX=None, vectorY=None, vectorZ=None, gridSize=10,
    tol=1e-5, iterMax=None, method='spline'
):
    """
    Interpolate properties with a minimum curvature interpolation
    :param locs:  numpy.array of size n-by-3 of point locations
    :param data: numpy.array of size n-by-m of values to be interpolated
    :param vectorX: numpy.ndarray Gridded locations along x-axis [Default:None]
    :param vectorY: numpy.ndarray Gridded locations along y-axis [Default:None]
    :param vectorZ: numpy.ndarray Gridded locations along z-axis [Default:None]
    :param gridSize: numpy float Grid point seperation in meters [DEFAULT:10]
    :param method: 'relaxation' || 'spline' [Default]
    :param tol: float tol=1e-5 [Default] Convergence criteria
    :param iterMax: int iterMax=None [Default] Maximum number of iterations

    :return: numpy.array of size nC-by-m of interpolated values

    """

    def av_extrap(n):
        """Define 1D averaging operator from cell-centers to nodes."""
        Av = (
            sp.spdiags(
                (0.5 * np.ones((n, 1)) * [1, 1]).T,
                [-1, 0],
                n + 1, n,
                format="csr"
            )
        )
        Av[0, 1], Av[-1, -2] = 0.5, 0.5
        return Av

    def aveCC2F(grid):
        "Construct the averaging operator on cell cell centers to faces."
        if grid.ndim == 1:
            aveCC2F = av_extrap(grid.shape[0])
        elif grid.ndim == 2:
            aveCC2F = sp.vstack((
                sp.kron(speye(grid.shape[1]), av_extrap(grid.shape[0])),
                sp.kron(av_extrap(grid.shape[1]), speye(grid.shape[0]))
            ), format="csr")
        elif grid.ndim == 3:
            aveCC2F = sp.vstack((
                kron3(
                    speye(grid.shape[2]), speye(grid.shape[1]), av_extrap(grid.shape[0])
                ),
                kron3(
                    speye(grid.shape[2]), av_extrap(grid.shape[1]), speye(grid.shape[0])
                ),
                kron3(
                    av_extrap(grid.shape[2]), speye(grid.shape[1]), speye(grid.shape[0])
                )
            ), format="csr")
        return aveCC2F

    assert locs.shape[0] == data.shape[0], ("Number of interpolated locs " +
                                            "must match number of data")

    if vectorY is not None:
        assert locs.shape[1] >= 2, (
                "Found vectorY as an input." +
                " Point locations must contain X and Y coordinates."
            )

    if vectorZ is not None:
        assert locs.shape[1] == 3, (
                "Found vectorZ as an input." +
                " Point locations must contain X, Y and Z coordinates."
            )

    ndim = locs.shape[1]

    # Define a new grid based on data extent
    if vectorX is None:
        xmin, xmax = locs[:, 0].min(), locs[:, 0].max()
        nCx = int((xmax-xmin)/gridSize)
        vectorX = xmin+np.cumsum(np.ones(nCx) * gridSize)

    if vectorY is None and ndim >= 2:
        ymin, ymax = locs[:, 1].min(), locs[:, 1].max()
        nCy = int((ymax-ymin)/gridSize)
        vectorY = ymin+np.cumsum(np.ones(nCy) * gridSize)

    if vectorZ is None and ndim == 3:
        zmin, zmax = locs[:, 2].min(), locs[:, 2].max()
        nCz = int((zmax-zmin)/gridSize)
        vectorZ = zmin+np.cumsum(np.ones(nCz) * gridSize)

    if ndim == 3:
        gridCy, gridCx, gridCz = np.meshgrid(vectorY, vectorX, vectorZ)
        gridCC = np.c_[mkvc(gridCx), mkvc(gridCy), mkvc(gridCz)]
    elif ndim == 2:
        gridCy, gridCx = np.meshgrid(vectorY, vectorX)
        gridCC = np.c_[mkvc(gridCx), mkvc(gridCy)]
    else:
        gridCC = vectorX

    # Build the cKDTree for distance lookup
    tree = cKDTree(locs)
    # Get the grid location
    d, ind = tree.query(gridCC, k=1)

    if method == 'relaxation':

        Ave = aveCC2F(gridCx)

        count = 0
        residual = 1.

        m = np.zeros((gridCC.shape[0], data.shape[1]))

        # Begin with neighrest primers
        for ii in range(m.shape[1]):
            # F = NearestNDInterpolator(mesh.gridCC[ijk], data[:, ii])
            m[:, ii] = data[ind, ii]

        while np.all([count < iterMax, residual > tol]):
            for ii in range(m.shape[1]):
                # F = NearestNDInterpolator(mesh.gridCC[ijk], data[:, ii])
                m[:, ii] = data[ind, ii]
            mtemp = m
            m = Ave.T * (Ave * m)
            residual = np.linalg.norm(m-mtemp)/np.linalg.norm(mtemp)
            count += 1

        print(count)
        return gridCC, m

    elif method == 'spline':

        ndat = locs.shape[0]
        # nC = int(nCx*nCy)

        A = np.zeros((ndat, ndat))
        for i in range(ndat):

            r = (locs[i, 0] - locs[:, 0])**2. + (locs[i, 1] - locs[:, 1])**2.
            A[i, :] = r.T * (np.log((r.T + 1e-8)**0.5) - 1.)

        # Solve system for the weights
        w = bicgstab(A, data, tol=1e-6)

        # Compute new solution
        # Reformat the line data locations but skip every n points for test
        nC = gridCC.shape[0]
        m = np.zeros(nC)

        # We can parallelize this part later
        for i in range(nC):

            r = (gridCC[i, 0] - locs[:, 0])**2. + (gridCC[i, 1] - locs[:, 1])**2.
            m[i] = np.sum(w[0] * r.T * (np.log((r.T + 1e-8)**0.5) - 1.))

        return gridCC, m.reshape(gridCx.shape, order='F')

    else:

        NotImplementedError("Only methods 'relaxation' || 'spline' are available" )


def plotDataHillside(x, y, z, axs=None, fill=True, contour=0,
                     vmin=None, vmax=None,
                     clabel=True, cmap='RdBu_r', ve=1., alpha=1., alphaHS=1.,
                     distMax=1000, midpoint=None, azdeg=315, altdeg=45):

    ls = LightSource(azdeg=azdeg, altdeg=altdeg)

    if x.ndim == 1:
        # Create grid of points
        vectorX = np.linspace(x.min(), x.max(), 1000)
        vectorY = np.linspace(y.min(), y.max(), 1000)

        X, Y = np.meshgrid(vectorX, vectorY)

        # Interpolate
        d_grid = griddata(np.c_[x, y], z, (X, Y), method='cubic')

        # Remove points beyond treshold
        tree = cKDTree(np.c_[x, y])
        xi = _ndim_coords_from_arrays((X, Y), ndim=2)
        dists, indexes = tree.query(xi)

        # Copy original result but mask missing values with NaNs
        d_grid[dists > distMax] = np.nan

    else:

        X, Y, d_grid = x, y, z

    class MidPointNorm(Normalize):
        def __init__(self, midpoint=None, vmin=None, vmax=None, clip=False):
            Normalize.__init__(self, vmin, vmax, clip)
            self.midpoint = midpoint

        def __call__(self, value, clip=None):
            if clip is None:
                clip = self.clip

            result, is_scalar = self.process_value(value)

            self.autoscale_None(result)

            if self.midpoint is None:
                self.midpoint = np.mean(value)
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

    im, CS = [], []
    if axs is None:
        axs = plt.subplot()

    if fill:
        extent = x.min(), x.max(), y.min(), y.max()
        im = axs.contourf(
            X, Y, d_grid, 50, vmin=vmin, vmax=vmax,
            cmap=cmap, norm=MidPointNorm(midpoint=midpoint), alpha=alpha
        )

        axs.imshow(ls.hillshade(d_grid.T, vert_exag=ve, dx=5., dy=5.),
                   cmap='gray_r', alpha=alphaHS,
                   extent=extent, origin='lower')

    if contour > 0:
        CS = axs.contour(
            X, Y, d_grid, int(contour), colors='k', linewidths=0.5
        )

        if clabel:
            plt.clabel(CS, inline=1, fontsize=10, fmt='%i')
    return im, CS
