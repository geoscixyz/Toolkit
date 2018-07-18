
import numpy as np
from . import Mag
from . import MathUtils
from . import Simulator
from SimPEG import PF, Utils, Mesh, Maps

import matplotlib.pyplot as plt


def blockModel():
    """
        List of paramters defining the synthetic block model
        [xCenter, yCenter, zCenter, Width, Height, Depth, rotationX, rotationZ]
    """

    parameters = [
            [2000, 500, -100, 5000,  4000,  1000, 60, 0],
            [-500, 0, -100, 300, 300, 300, -30, 0],
            [200, 100, -100, 4000, 100, 1000, 55, 10],
           ]

    susceptibility = [0.075, 0.1, -0.05,  0.005]

    return parameters, susceptibility


def setSyntheticProblem(
        rxLocs, EarthField=[50000, 90, 0],
        topo=None, discretize=False, plotSections=False
     ):
    """
        Set the synthetic problem with multiple blocks.
        Output the figure used in the doc
    """

    # Import the data
    # topo = np.genfromtxt('TKCtopoDwnS.dat', skip_header=1)
    # fileName = './DIGHEM_Mag_floor10nt_25m.obs'

    # Shift everything centered at origin
    cntr = np.mean(rxLocs, axis=0)
    rxLocs -= np.kron(np.ones((rxLocs.shape[0], 1)), cntr)

    if topo is not None:
        topo -= np.kron(np.ones((topo.shape[0], 1)), cntr)

    # Create survey
    survey = Mag.createMagSurvey(rxLocs, EarthField=EarthField)
    cntr = np.mean(rxLocs, axis=0)

    if discretize:
        # Mesh discretization for plotting
        hx = [(10, 320)]
        hy = [(10, 320)]
        hz = [(10, 80)]
        x0 = np.min(rxLocs, axis=0)
        x0[2] -= 600

        # Create a mesh
        mesh = Mesh.TensorMesh([hx, hy, hz], x0=x0)
        model = np.zeros(mesh.nC)

        if topo is not None:
            actv = Utils.modelutils.surface2ind_topo(mesh, topo)
        else:
            actv = np.ones(mesh.nC, dtype='bool')
        actvMap = Maps.InjectActiveCells(mesh, actv, np.nan)

    # Cycle through the parameters, create blocks for forward and
    # discretize on to the mesh
    prisms = []

    # User defined parameters for the blocks
    params, suscs = blockModel()

    # Create the synthetic blocks model and place
    # it at the center of the survey
    for param, susc in zip(params, suscs):

        prism = Simulator.definePrism()
        prism.x0, prism.y0, prism.z0 = cntr[0]+param[0], cntr[1]+param[1], rxLocs[:, 2].min() +param[2]
        prism.dx, prism.dy, prism.dz = param[3], param[4], param[5]
        prism.pdec, prism.pinc = param[6], param[7]
        prisms.append(prism)

        # Forward model data
        prob = Mag.Problem(prism=prism, survey=survey)
        prob.susc = susc
        survey.dobs += prob.fields()[0]

        if discretize:
            # Discretize onto mesh
            X, Y, Z = np.meshgrid(prism.xn, prism.yn, prism.zn)
            pts = np.c_[Utils.mkvc(X), Utils.mkvc(Y), Utils.mkvc(Z)]

            xyz = MathUtils.rotate(
                pts, np.r_[prism.xc, prism.yc, prism.zc],
                prism.pinc, prism.pdec
            )
            ind = Utils.ModelBuilder.PolygonInd(mesh, xyz)

            model[ind] += susc

    if np.all([discretize, plotSections]):
        model = model[actv]

        # All ploting functions
        fig = plt.figure(figsize=(10, 6))
        axs = plt.subplot(1, 2, 1)
        indy = int(mesh.vnC[1]/2)-18
        indz = -40

        # Plot horizontal section
        im = mesh.plotSlice(
            actvMap*model, normal='Z', ax=axs,
            ind=indz, clim=[0.0, 0.1], pcolorOpts={'cmap': 'jet'}
        )

        a = np.r_[rxLocs[:, 0].min(), mesh.vectorCCy[indy]]
        b = np.r_[rxLocs[:, 0].max(), mesh.vectorCCy[indy]]

        plt.scatter(rxLocs[:, 0], rxLocs[:, 1], 10, c='k', marker='.')
        plt.plot(np.r_[a[0], b[0]], np.r_[a[1], b[1]], 'r--')

        axs.set_title(
            'Plan view'
        )
        axs.set_xlabel('Easting (m)')
        axs.set_ylabel('Northing (m)')
        axs.set_aspect('equal')
        axs.set_xlim(rxLocs[:, 0].min()-100, rxLocs[:, 0].max()+100)
        axs.set_ylim(rxLocs[:, 1].min()-100, rxLocs[:, 1].max()+100)

        # Plot vertical section
        axs = plt.subplot(1, 2, 2)
        indy = int(mesh.vnC[1]/2)-18
        im = mesh.plotSlice(
            actvMap*model, normal='Y', ax=axs,
            ind=indy, clim=[0.0, 0.1], pcolorOpts={'cmap': 'jet'}
        )

        cbar = plt.colorbar(im[0], orientation='horizontal')
        cbar.set_label('SI')

        Simulator.plotProfile2D(
            rxLocs, rxLocs[:, -1], a, b, 10, ax=axs,
            coordinate_system='xProfile', ylabel='k:'
        )

        if topo is not None:
            Simulator.plotProfile2D(
                topo, topo[:, -1], a, b, 10, ax=axs,
                plotStr=['k-'],
                coordinate_system='xProfile', ylabel=''
            )

        axs.set_title(
            'EW Section'
        )
        axs.set_ylim(-1500, 100)
        axs.set_aspect('equal')
        axs.set_xlabel('Easting (m)')
        axs.set_ylabel('Depth (m)')
        axs.yaxis.set_label_position("right")
        fig.savefig('./images/SyntheticModel.png', bbox_inches='tight')

        plt.close()
    return survey
