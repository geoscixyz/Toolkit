
import numpy as np
from . import Mag
from . import MagUtils
from . import Simulator
from SimPEG import PF, Utils, Mesh, Maps

import matplotlib.pyplot as plt


def setSyntheticProblem():
    """
        Set the synthetic problem with multiple blocks.
        Output the figure used in the doc
    """

    # Import the data
    topo = np.genfromtxt('TKCtopoDwnS.dat', skip_header=1)
    fileName = './DIGHEM_Mag_floor10nt_25m.obs'

    # Create the survey
    xyzd = np.genfromtxt(fileName, skip_header=3)
    xyzd[:, -2] = 0
    B = np.r_[60308, 83.8, 25.4]
    survey = Mag.createMagSurvey(xyzd, B)
    cntr = np.mean(xyzd[:, :2], axis=0)
    rxLocs = survey.srcField.rxList[0].locs

    # User defined parameters for the blocks

    params = [[2000, 500, -100, 5000,  4000,  500, 60, 0],
              [-500, 0, -100, 300, 300, 300, -30, 0],
              [200, 100, -100, 4000, 100, 500, 55, 10],
              ]
    suscs = [0.075, 0.1, -0.05,  0.005]

    # Mesh discretization for plotting
    hx = [(10, 320)]
    hy = [(10, 320)]
    hz = [(10, 80)]
    x0 = np.min(xyzd[:, :3], axis=0)
    x0[2] -= 600

    # Create a mesh
    mesh = Mesh.TensorMesh([hx, hy, hz], x0=x0)
    actv = Utils.modelutils.surface2ind_topo(mesh, topo)
    actvMap = Maps.InjectActiveCells(mesh, actv, np.nan)

    model = np.zeros(mesh.nC)

    # Cycle through the parameters, create blocks for forward and
    # discretize on to the mesh
    prisms = []

    # Create the synthetic blocks model and place
    # it at the center of the survey
    for param, susc in zip(params, suscs):

        prism = Simulator.definePrism()
        prism.x0, prism.y0, prism.z0 = cntr[0]+param[0], cntr[1]+param[1], xyzd[:, 2].min() +param[2]
        prism.dx, prism.dy, prism.dz = param[3], param[4], param[5]
        prism.pdec, prism.pinc = param[6], param[7]

        prisms.append(prism)

        # Forward model data
        prob = Mag.problem(prism=prism, survey=survey)
        prob.susc = susc
        survey.dobs += prob.fields()[0]

        # Discretize onto mesh
        X, Y, Z = np.meshgrid(prism.xn, prism.yn, prism.zn)
        pts = np.c_[Utils.mkvc(X), Utils.mkvc(Y), Utils.mkvc(Z)]

        xyz = MagUtils.rotate(pts, np.r_[prism.xc, prism.yc, prism.zc], prism.pinc, prism.pdec)
        ind = Utils.ModelBuilder.PolygonInd(mesh, xyz)

        model[ind] = susc

    model = model[actv]

    # All ploting functions
    fig = plt.figure(figsize=(10, 6))
    axs = plt.subplot(1, 2, 1)
    indy = int(mesh.vnC[1]/2)-18
    indz = -40

    # Plot horizontal section
    im = mesh.plotSlice(actvMap*model, normal='Z', ax=axs, ind=indz,clim=[-0.1, 0.1], pcolorOpts={'cmap':'jet'})

    a, b = np.r_[rxLocs[:, 0].min(), mesh.vectorCCy[indy]], np.r_[rxLocs[:, 0].max(), mesh.vectorCCy[indy]]

    plt.scatter(rxLocs[:, 0], rxLocs[:, 1], 10, c='k', marker='.')
    plt.plot(np.r_[a[0], b[0]], np.r_[a[1], b[1]], 'r--')    # Simulator.plotProfile2D(topo, np.r_[rxLocs[:, 0].min(), mesh.vectorCCy[indy]], np.r_[rxLocs[:, 0].max(), mesh.vectorCCy[indy]], 10, ax=axs, plotStr='k-', coordinate_system = 'xProfile')

    axs.set_title('Plan view: ' + str(int(mesh.vectorCCz[indz]))+ " m depth")
    axs.set_aspect('equal')
    axs.set_xlim(rxLocs[:, 0].min()-100, rxLocs[:, 0].max()+100)
    axs.set_ylim(rxLocs[:, 1].min()-100, rxLocs[:, 1].max()+100)

    # Plot vertical section
    axs = plt.subplot(1, 2, 2)
    indy = int(mesh.vnC[1]/2)-18
    im = mesh.plotSlice(actvMap*model, normal='Y', ax=axs, ind=indy,clim=[-0.1, 0.1], pcolorOpts={'cmap':'jet'})
    plt.colorbar(im[0], orientation='horizontal')
    rxLocs = survey.srcField.rxList[0].locs
    Simulator.plotProfile2D(rxLocs, a, b, 10, ax=axs, coordinate_system = 'xProfile')
    # Simulator.plotProfile2D(topo, np.r_[rxLocs[:, 0].min(), mesh.vectorCCy[indy]], np.r_[rxLocs[:, 0].max(), mesh.vectorCCy[indy]], 10, ax=axs, plotStr='k-', coordinate_system = 'xProfile')

    axs.set_title('EW Section: ' + str(int(mesh.vectorCCy[indy]))+ " N")
    axs.set_aspect('equal')

    # plt.show()

    fig.savefig('./images/SyntheticModel.png', bbox_inches='tight')

    return survey
