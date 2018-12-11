
import numpy as np
from GeoToolkit.Mag import Simulator, DataIO, MathUtils, Mag
# from SimPEG import PF, Utils, Mesh, Maps
import ipywidgets as widgets
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
from scipy.interpolate import griddata, interp1d
from scipy.interpolate import NearestNDInterpolator, LinearNDInterpolator
# from SimPEG.Utils import mkvc, ndgrid, uniqueRows


def blockModel():
    """
        List of paramters defining the synthetic block model
        [xCenter, yCenter, zCenter, Width, Height, Depth, rotationX, rotationZ]
    """

    parameters = [
            [2000, 500, -100, 5000,  4000,  1000, 60, 0],
            [-500, 0, -100, 300, 300, 300, -30, 0],
            [400, 100, -100, 4000, 100, 1000, 55, 10],
           ]

    susceptibility = [0.075, 0.1, -0.05,  0.005]

    return parameters, susceptibility


def setSyntheticProblem(
        rxLocs, EarthField=[50000, 90, 0], discretize=False
     ):
    """
        Set the synthetic problem with multiple blocks.
        Output the figure used in the doc
    """

    if discretize:
        # Mesh discretization for plotting
        hx = [(10, 320)]
        hy = [(10, 320)]
        hz = [(10, 120)]
        x0 = np.min(rxLocs, axis=0)
        x0[2] -= 1000

        # Create a mesh
        mesh = Mesh.TensorMesh([hx, hy, hz], x0=x0)
        model = np.zeros(mesh.nC)
    else:
        mesh = []
        model = []

    # Shift everything centered at origin
    cntr = np.mean(rxLocs, axis=0)
    rxLocs -= np.kron(np.ones((rxLocs.shape[0], 1)), cntr)

    # Create survey
    survey = Mag.createMagSurvey(rxLocs, EarthField=EarthField)
    cntr = np.mean(rxLocs, axis=0)

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
        survey._dobs = survey.dobs + prob.fields()[0]

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

    return survey, mesh, model


def meshBuilder(xyz, h, padDist, meshGlobal=None,
                expFact=1.3,
                meshType='TENSOR',
                verticalAlignment='top'):
    """
        Function to quickly generate a Tensor mesh
        given a cloud of xyz points, finest core cell size
        and padding distance.
        If a meshGlobal is provided, the core cells will be centered
        on the underlaying mesh to reduce interpolation errors.

        :param numpy.ndarray xyz: n x 3 array of locations [x, y, z]
        :param numpy.ndarray h: 1 x 3 cell size for the core mesh
        :param numpy.ndarray padDist: 2 x 3 padding distances [W,E,S,N,Down,Up]
        [OPTIONAL]
        :param numpy.ndarray padCore: Number of core cells around the xyz locs
        :object SimPEG.Mesh: Base mesh used to shift the new mesh for overlap
        :param float expFact: Expension factor for padding cells [1.3]
        :param string meshType: Specify output mesh type: "TensorMesh"

        RETURNS:
        :object SimPEG.Mesh: Mesh object

    """

    assert meshType in ['TENSOR', 'TREE'], ('Revise meshType. Only ' +
                                            ' TENSOR | TREE mesh ' +
                                            'are implemented')

    # Get extent of points
    limx = np.r_[xyz[:, 0].max(), xyz[:, 0].min()]
    limy = np.r_[xyz[:, 1].max(), xyz[:, 1].min()]
    limz = np.r_[xyz[:, 2].max(), xyz[:, 2].min()]

    # Get center of the mesh
    midX = np.mean(limx)
    midY = np.mean(limy)
    midZ = np.mean(limz)

    nCx = int(limx[0]-limx[1]) / h[0]
    nCy = int(limy[0]-limy[1]) / h[1]
    nCz = int(limz[0]-limz[1]+int(np.min(np.r_[nCx, nCy])/3)) / h[2]

    if meshType == 'TENSOR':
        # Make sure the core has odd number of cells for centereing
        # on global mesh
        if meshGlobal is not None:
            nCx += 1 - int(nCx % 2)
            nCy += 1 - int(nCy % 2)
            nCz += 1 - int(nCz % 2)

        # Figure out paddings
        def expand(dx, pad):
            L = 0
            nC = 0
            while L < pad:
                nC += 1
                L = np.sum(dx * expFact**(np.asarray(range(nC))+1))

            return nC

        # Figure number of padding cells required to fill the space
        npadEast = expand(h[0], padDist[0, 0])
        npadWest = expand(h[0], padDist[0, 1])
        npadSouth = expand(h[1], padDist[1, 0])
        npadNorth = expand(h[1], padDist[1, 1])
        npadDown = expand(h[2], padDist[2, 0])
        npadUp = expand(h[2], padDist[2, 1])

        # Create discretization
        hx = [(h[0], npadWest, -expFact),
              (h[0], nCx),
              (h[0], npadEast, expFact)]
        hy = [(h[1], npadSouth, -expFact),
              (h[1], nCy), (h[1],
              npadNorth, expFact)]
        hz = [(h[2], npadDown, -expFact),
              (h[2], nCz),
              (h[2], npadUp, expFact)]

        # Create mesh
        mesh = Mesh.TensorMesh([hx, hy, hz], 'CC0')

        # Re-set the mesh at the center of input locations
        # Set origin
        if verticalAlignment == 'center':
            mesh.x0 = [midX-np.sum(mesh.hx)/2., midY-np.sum(mesh.hy)/2., midZ-np.sum(mesh.hz)/2.]
        elif verticalAlignment == 'top':
            mesh.x0 = [midX-np.sum(mesh.hx)/2., midY-np.sum(mesh.hy)/2., limz[0]-np.sum(mesh.hz)]
        else:
            assert NotImplementedError("verticalAlignment must be 'center' | 'top'")

    elif meshType == 'TREE':

        # Figure out full extent required from input
        extent = np.max(np.r_[nCx * h[0] + padDist[0, :].sum(),
                              nCy * h[1] + padDist[1, :].sum(),
                              nCz * h[2] + padDist[2, :].sum()])

        maxLevel = int(np.log2(extent/h[0]))+1

        # Number of cells at the small octree level
        # For now equal in 3D

        nCx, nCy, nCz = 2**(maxLevel), 2**(maxLevel), 2**(maxLevel)
        # nCy = 2**(int(np.log2(extent/h[1]))+1)
        # nCz = 2**(int(np.log2(extent/h[2]))+1)

        # Define the mesh and origin
        # For now cubic cells
        mesh = Mesh.TreeMesh([np.ones(nCx)*h[0],
                              np.ones(nCx)*h[1],
                              np.ones(nCx)*h[2]])

        # Set origin
        if verticalAlignment == 'center':
            mesh.x0 = np.r_[-nCx*h[0]/2.+midX, -nCy*h[1]/2.+midY, -nCz*h[2]/2.+midZ]
        elif verticalAlignment == 'top':
            mesh.x0 = np.r_[-nCx*h[0]/2.+midX, -nCy*h[1]/2.+midY, -(nCz-1)*h[2] + limz.max()]
        else:
            assert NotImplementedError("verticalAlignment must be 'center' | 'top'")

    return mesh


def refineTree(mesh, xyz, finalize=False, dtype="point", nCpad=[1, 1, 1]):

    maxLevel = int(np.log2(mesh.hx.shape[0]))

    if dtype == "point":

        mesh.insert_cells(xyz, np.ones(xyz.shape[0])*maxLevel, finalize=False)

        stencil = np.r_[
                np.ones(nCpad[0]),
                np.ones(nCpad[1])*2,
                np.ones(nCpad[2])*3
            ]

        # Reflect in the opposite direction
        vec = np.r_[stencil[::-1], 1, stencil]
        vecX, vecY, vecZ = np.meshgrid(vec, vec, vec)
        gridLevel = np.maximum(np.maximum(vecX,
                               vecY), vecZ)
        gridLevel = np.kron(np.ones((xyz.shape[0], 1)), gridLevel)

        # Grid the coordinates
        vec = np.r_[-stencil[::-1], 0, stencil]
        vecX, vecY, vecZ = np.meshgrid(vec, vec, vec)
        offset = np.c_[
            mkvc(np.sign(vecX)*2**np.abs(vecX) * mesh.hx.min()),
            mkvc(np.sign(vecY)*2**np.abs(vecY) * mesh.hx.min()),
            mkvc(np.sign(vecZ)*2**np.abs(vecZ) * mesh.hx.min())
        ]

        # Replicate the point locations in each offseted grid points
        newLoc = (
            np.kron(xyz, np.ones((offset.shape[0], 1))) +
            np.kron(np.ones((xyz.shape[0], 1)), offset)
        )

        mesh.insert_cells(
            newLoc, maxLevel-mkvc(gridLevel)+1, finalize=finalize
        )

    elif dtype == 'surface':

        # Get extent of points
        limx = np.r_[xyz[:, 0].max(), xyz[:, 0].min()]
        limy = np.r_[xyz[:, 1].max(), xyz[:, 1].min()]

        F = NearestNDInterpolator(xyz[:, :2], xyz[:, 2])
        zOffset = 0
        # Cycle through the first 3 octree levels
        for ii in range(3):

            dx = mesh.hx.min()*2**ii

            # Horizontal offset
            xyOff = dx * 2

            nCx = int(limx[0]-limx[1] + 2 * xyOff) / dx
            nCy = int(limy[0]-limy[1] + 2 * xyOff) / dx

            # Create a grid at the octree level in xy
            CCx, CCy = np.meshgrid(
                np.linspace(limx[1]-xyOff, limx[0]+xyOff, nCx),
                np.linspace(limy[1]-xyOff, limy[0]+xyOff, nCy)
            )

            z = F(mkvc(CCx), mkvc(CCy))

            for level in range(int(nCpad[ii])):

                mesh.insert_cells(
                    np.c_[mkvc(CCx), mkvc(CCy), z-zOffset], np.ones_like(z)*maxLevel-ii,
                    finalize=False
                )

                zOffset += dx

        if finalize:
            mesh.finalize()

    else:
        NotImplementedError("Only dtype='points' has been implemented")

    return mesh


if __name__ == '__main__':

    # Load data and topo and build default model
    assetDir = './../../docs/Notebooks/assets/TKC/'
    survey = Mag.readMagneticsObservations('DIGHEM_Mag_floor10nt_25m.obs')
    topo = np.genfromtxt('TKCtopoDwnS.dat', skip_header=1)
    locs = survey.rxLoc

    # Build the problem
    survey, mesh, model = setSyntheticProblem(locs, topo=topo, discretize=True)

    if topo is not None:
        topo -= np.kron(np.ones((topo.shape[0], 1)), cntr)

    if topo is not None:
        actv = Utils.modelutils.surface2ind_topo(mesh, topo)
    else:
        actv = np.ones(mesh.nC, dtype='bool')

    model = model[actv]
    actvMap = Maps.InjectActiveCells(mesh, actv, np.nan)

    # All ploting functions
    fig = plt.figure(figsize=(10, 6))
    axs = plt.subplot(1, 2, 1)
    indy = int(mesh.vnC[1]/2)-18
    indz = -32

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
        rxLocs[:, 0], rxLocs[:, 1], rxLocs[:, -1], a, b, 10, ax=axs,
        coordinate_system='xProfile', ylabel='k:'
    )

    if topo is not None:
        Simulator.plotProfile2D(
            topo[:, 0], topo[:, 1], topo[:, -1], a, b, 10, ax=axs,
            plotStr=['k-'],
            coordinate_system='xProfile', ylabel=''
        )

    axs.set_title(
        'EW Section'
    )
    axs.set_ylim(-1000, 100)
    axs.set_aspect('equal')
    axs.set_xlabel('Easting (m)')
    axs.set_ylabel('Depth (m)')
    axs.yaxis.set_label_position("right")
    fig.savefig('./images/SyntheticModel.png', bbox_inches='tight')

    plt.close()

    # Save contours to shapefile
    zSlice = (actvMap*model).reshape(mesh.vnC, order='F')[:,:,indz]
    contours = plt.contour(mesh.vectorCCx,mesh.vectorCCy, zSlice.T, [0.03,0.075,0.09]).allsegs[0]

    DataIO.exportShapefile(
            contours, [1]*len(contours),
            saveAs='Synthetic_Zcontours',
            directory="./assets/Synthetic")
