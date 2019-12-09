import numpy as np
import scipy as sp
from scipy.spatial import cKDTree
# from SimPEG.Utils import mkvc, speye
from scipy.sparse.linalg import bicgstab
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
from scipy.interpolate import LinearNDInterpolator
from tqdm import tqdm


def rotationMatrix(inc, dec, normal=True):
    """
        Take an inclination and declination angle and return a rotation matrix

    """

    phi = -np.deg2rad(np.asarray(inc))
    theta = -np.deg2rad(np.asarray(dec))

    Rx = np.asarray([[1, 0, 0],
                    [0, np.cos(phi), -np.sin(phi)],
                    [0, np.sin(phi), np.cos(phi)]])

    Rz = np.asarray([[np.cos(theta), -np.sin(theta), 0],
                    [np.sin(theta), np.cos(theta), 0],
                    [0, 0, 1]])

    if normal:
        R = Rz.dot(Rx)
    else:
        R = Rx.dot(Rz)

    return R


def dipazm_2_xyz(dip, azm_N):
    """
    dipazm_2_xyz(dip,azm_N)

    Function converting degree angles for dip and azimuth from north to a
    3-components in cartesian coordinates.

    INPUT
    dip     : Value or vector of dip from horizontal in DEGREE
    azm_N   : Value or vector of azimuth from north in DEGREE

    OUTPUT
    M       : [n-by-3] Array of xyz components of a unit vector in cartesian

    Created on Dec, 20th 2015

    @author: dominiquef
    """

    # Modify azimuth from North to Cartesian-X
    azm_X = (450. - azm_N) % 360.

    M = np.zeros(3)

    inc = -np.deg2rad(np.asarray(dip))
    dec = np.deg2rad(azm_X)

    M[0] = np.cos(inc) * np.cos(dec)
    M[1] = np.cos(inc) * np.sin(dec)
    M[2] = np.sin(inc)

    return M


def rotate(xyz, center, theta, phi):
    """
      Rotate scatter points in column format around a center location

      INPUT
      :param: xyz nDx3 matrix
      :param: center xyz location of rotation
      :param: theta angle rotation around x-axis
      :param: phi angle rotation around z-axis

    """
    xyz -= np.kron(np.ones((xyz.shape[0], 1)), np.r_[center])

    R = rotationMatrix(-theta, phi)

    xyzRot = R.dot(xyz.T).T + np.kron(np.ones((xyz.shape[0], 1)), np.r_[center])

    return xyzRot


def minCurvatureInterp(
    xyz, data, xyzOut=None,
    vectorX=None, vectorY=None, vectorZ=None, gridSize=10,
    tol=1e-5, iterMax=None, method='spline', maxDistance=None,
    nPoints=5, overlap=0
):
    """
    Interpolate properties with a minimum curvature interpolation
    :param xyz:  numpy.array of size n-by-3 of point locations
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

    # def av_extrap(n):
    #     """Define 1D averaging operator from cell-centers to nodes."""
    #     Av = (
    #         sp.sparse.spdiags(
    #             (0.5 * np.ones((n, 1)) * [1, 1]).T,
    #             [-1, 0],
    #             n + 1, n,
    #             format="csr"
    #         )
    #     )
    #     Av[0, 1], Av[-1, -2] = 0.5, 0.5
    #     return Av

    # def aveCC2F(grid):
    #     "Construct the averaging operator on cell cell centers to faces."
    #     if grid.ndim == 1:
    #         aveCC2F = av_extrap(grid.shape[0])
    #     elif grid.ndim == 2:
    #         aveCC2F = sp.vstack((
    #             sp.kron(speye(grid.shape[1]), av_extrap(grid.shape[0])),
    #             sp.kron(av_extrap(grid.shape[1]), speye(grid.shape[0]))
    #         ), format="csr")
    #     elif grid.ndim == 3:
    #         aveCC2F = sp.vstack((
    #             kron3(
    #                 speye(grid.shape[2]), speye(grid.shape[1]), av_extrap(grid.shape[0])
    #             ),
    #             kron3(
    #                 speye(grid.shape[2]), av_extrap(grid.shape[1]), speye(grid.shape[0])
    #             ),
    #             kron3(
    #                 av_extrap(grid.shape[2]), speye(grid.shape[1]), speye(grid.shape[0])
    #             )
    #         ), format="csr")
    #     return aveCC2F

    # assert xyz.shape[0] == data.shape[0], ("Number of interpolated xyz " +
    #                                         "must match number of data")

    if vectorY is not None:
        assert xyz.shape[1] >= 2, (
                "Found vectorY as an input." +
                " Point locations must contain X and Y coordinates."
            )

    if vectorZ is not None:
        assert xyz.shape[1] == 3, (
                "Found vectorZ as an input." +
                " Point locations must contain X, Y and Z coordinates."
            )

    ndim = xyz.shape[1]

    if xyzOut is None:
        # Define a new grid based on data extent
        if vectorX is None:
            xmin, xmax = xyz[:, 0].min(), xyz[:, 0].max()
            nCx = int((xmax-xmin)/gridSize)
            vectorX = xmin+np.cumsum(np.ones(nCx) * gridSize)

        if vectorY is None and ndim >= 2:
            ymin, ymax = xyz[:, 1].min(), xyz[:, 1].max()
            nCy = int((ymax-ymin)/gridSize)
            vectorY = ymin+np.cumsum(np.ones(nCy) * gridSize)

        if vectorZ is None and ndim == 3:
            zmin, zmax = xyz[:, 2].min(), xyz[:, 2].max()
            nCz = int((zmax-zmin)/gridSize)
            vectorZ = zmin+np.cumsum(np.ones(nCz) * gridSize)

        if ndim == 3:
            gridCx, gridCy, gridCz = np.meshgrid(vectorX, vectorY, vectorZ)
            gridCC = np.c_[gridCx.flatten(order='F'), gridCy.flatten(order='F'), gridCz.flatten(order='F')]
        elif ndim == 2:
            gridCx, gridCy = np.meshgrid(vectorX, vectorY)
            gridCC = np.c_[gridCx.flatten(order='F'), gridCy.flatten(order='F')]
        else:
            gridCC = vectorX
    else:
        # Use the locations given by user
        gridCC = xyzOut

    if method == 'relaxation':

        # Ave = aveCC2F(gridCx)

        # count = 0
        # residual = 1.

        # m = np.zeros((gridCC.shape[0], data.shape[1]))

        # # Begin with neighrest primers
        # for ii in range(m.shape[1]):
        #     # F = NearestNDInterpolator(mesh.gridCC[ijk], data[:, ii])
        #     m[:, ii] = data[ind, ii]

        # while np.all([count < iterMax, residual > tol]):
        #     for ii in range(m.shape[1]):
        #         # F = NearestNDInterpolator(mesh.gridCC[ijk], data[:, ii])
        #         m[:, ii] = data[ind, ii]
        #     mtemp = m
        #     m = Ave.T * (Ave * m)
        #     residual = np.linalg.norm(m-mtemp)/np.linalg.norm(mtemp)
        #     count += 1

        # return gridCC, m
        NotImplementedError("method='relaxation' not yet implemented")

    elif method == 'spline':

        tree = cKDTree(xyz[:, :2])

        # First tile the problem to avoid large memory issues
        tiles = tileSurveyPoints(xyz, 1000, overlap=[overlap, overlap])

        # Get the XY limits of each tile
        X1, Y1 = tiles[0][:, 0], tiles[0][:, 1]
        X2, Y2 = tiles[1][:, 0], tiles[1][:, 1]

        # If max distance given, cut out points
        if maxDistance is not None:

            tree = cKDTree(xyz[:, :2])
            # xi = _ndim_coords_from_arrays((gridCC[:,0], gridCC[:,1]), ndim=2)
            dists, _ = tree.query(gridCC)

            # Copy original result but mask missing values with NaNs
            inRadius = dists < maxDistance

        else:
            inRadius = np.ones(gridCC.shape[0], dtype='bool')

        m = np.zeros(gridCC.shape[0])
        maskOut = np.zeros(gridCC.shape[0], dtype='bool')
        weights = np.zeros(gridCC.shape[0])
        # Loop over the tiles and compute interpolation
        for tt in range(X1.shape[0]):

            # Grab the interpolated points within a tile
            lims = np.r_[X1[tt], X2[tt], Y1[tt], Y2[tt]]

            inTile = np.all(
                [
                    gridCC[:, 0] >= lims[0], gridCC[:, 0] <= lims[1],
                    gridCC[:, 1] >= lims[2], gridCC[:, 1] <= lims[3]
                ], axis=0
            )

            inAll = (inTile * inRadius) == 1
            maskOut[inAll] = True

            # Tapper the grid
            r = (
                (gridCC[inAll, 0] - np.mean(gridCC[inAll, 0]))**2. +
                (gridCC[inAll, 1] - np.mean(gridCC[inAll, 1]))**2.
            )**0.5

            tapper = (1.01 - r/r.max())

            rQuery, indexes = tree.query(gridCC[inAll, :], k=nPoints)

            # Get closest querry points
            indx = np.unique(indexes)

            if tt == 0:
                baseLine = len(indx)

            ndat = xyz.shape[0]

            X, Y = np.meshgrid(xyz[indx, 0], xyz[indx, 1])
            # for i in range(ndat):

            r = (X.T - X)**2. + (Y.T - Y)**2. + 1e-8
            A = r.T * (np.log((r.T)**0.5) - 1.)

            # # Solve system for the green parameters
            # Scale the tolerance on the number data points
            g = sp.sparse.linalg.bicgstab(A, data[indx], tol=tol*baseLine/len(indx))

            # We can parallelize this part later
            mInterp = np.zeros(inAll.sum())
            for ii in range(indx.shape[0]):

                r = (gridCC[inAll, 0] - xyz[indx[ii], 0])**2. + (gridCC[inAll, 1] - xyz[indx[ii], 1])**2. + 1e-8
                mInterp += g[0][ii] * (r).T * (np.log((r.T)**0.5) - 1.)

            m[inAll] += (mInterp * tapper)
            weights[inAll] += tapper

        m[maskOut] /= weights[maskOut]
        m[maskOut == 0] = np.nan

        if xyzOut is None:
            m = m.reshape(gridCx.shape, order='F')

        return gridCC, m

    else:

        NotImplementedError("Only methods 'relaxation' || 'spline' are available" )


def decimalDegrees2DMS(value, type):
    """
        Converts a Decimal Degree Value into
        Degrees Minute Seconds Notation.

        Pass value as double
        type = {Latitude or Longitude} as string

        returns a string as D:M:S:Direction
        created by: anothergisblog.blogspot.com
    """
    degrees = int(value)
    submin = abs( (value - int(value) ) * 60)
    minutes = int(submin)
    subseconds = abs((submin-int(submin)) * 60)
    direction = ""
    if type == "Longitude":
        if degrees < 0:
            direction = "W"
        elif degrees > 0:
            direction = "E"
        else:
            direction = ""
    elif type == "Latitude":
        if degrees < 0:
            direction = "S"
        elif degrees > 0:
            direction = "N"
        else:
            direction = ""
    notation = str(degrees) + ":" + str(minutes) + ":" +\
               str(subseconds)[0:5] + "" + direction
    return notation


def estimateDepth(grid, method="tiltAngleDerivative"):
    """
        Function to estimate depth of anomalies
        from tilt angle. Distance between the 0 and
        pi/4 contour.

        INPUT
        :object: grid Grid object from

        OUTPUT
        :xyDepth:
    """

    assert method.lower() in ["tiltAngle".lower(), "tiltAngleDerivative".lower()], "'method' should be 'tiltAngle' or 'tiltAngleDerivative'"

    upward_height = 0
    print(grid.heightUC)
    if getattr(grid, 'heightUC', None) is not None:
        upward_height += grid.heightUC

    tilt = grid.tiltAngle
    X, Y = np.meshgrid(grid.hx, grid.hy)
    C_0 = plt.contour(X, Y, tilt, levels=[0], colors='k')
    plt.close()
    xy = []
    depth = []
    if method == 'tiltAngle':

        C_45 = plt.contour(X, Y, tilt, levels=[np.pi/4.], colors='r')
        plt.close()

        # Get zero contour nodes
        # xy0 = np.vstack(C_0.allsegs[0])

        # Get 45 contour nodes
        xy45 = np.vstack(C_45.allsegs[0])

        # Create ckDtree for shortest distance
        tree = cKDTree(xy45)

        # Compute shortest distance between pair of points
        for contour in C_0.allsegs[0]:

            if contour.shape[0] == 1:
                continue

            # Query two closest points to each nodes of zero contour
            d, indx = tree.query(contour, k=2)

            length = (
                (xy45[indx[:, 1], 0] - xy45[indx[:, 0], 0])**2. +
                (xy45[indx[:, 1], 1] - xy45[indx[:, 0], 1])**2.)

            indL = length > 0

            depth += [upward_height + np.abs(
                (xy45[indx[indL, 1], 1] - xy45[indx[indL, 0], 1])*contour[indL, 0] -
                (xy45[indx[indL, 1], 0] - xy45[indx[indL, 0], 0])*contour[indL, 1] +
                xy45[indx[indL, 1], 0] * xy45[indx[indL, 0], 1] -
                xy45[indx[indL, 1], 1] * xy45[indx[indL, 0], 0]) / length[indL]**0.5]

            xy += [contour[indL, :]]
    else:

        # Compute the total derivative of the tilt angle
        # Compute tilt angle derivative
        grad_tilt = np.gradient(grid.tiltAngle, grid.dx, grid.dy)

        THD_tilt = (grad_tilt[0]**2. + grad_tilt[1]**2. + 1e-8)**0.5

        # Interpolate value on 0 contour
        tri2D = Delaunay(grid.gridCC[:, :2])
        F = LinearNDInterpolator(tri2D, THD_tilt.flatten(order='F'))

        for contour in C_0.allsegs[0]:

            estimate_z = 1./F(contour[:, 0], contour[:, 1])
            depth += [upward_height + estimate_z[~np.isnan(estimate_z)]]
            xy += [contour[~np.isnan(estimate_z), :]]

        # if getattr(self, 'heightUC', None) is None
    return xy, depth


def tileSurveyPoints(xyLocs, maxNpoints, overlap=[0, 0]):
    """
        Function to tile an survey points into smaller square subsets of points

        :param numpy.ndarray xyLocs: n x 2 array of locations [x,y]
        :param integer maxNpoints: maximum number of points in each tile

        RETURNS:
        :param numpy.ndarray: Return a list of arrays  for the SW and NE
                            limits of each tiles

    """

    # Initialize variables
    nNx = 2
    nNy = 1
    nObs = 1e+8
    countx = 0
    county = 0
    xlim = [xyLocs[:, 0].min(), xyLocs[:, 0].max()]
    ylim = [xyLocs[:, 1].min(), xyLocs[:, 1].max()]

    # Refine the brake recursively
    while nObs > maxNpoints:

        nObs = 0

        if countx > county:
            nNx += 1
        else:
            nNy += 1

        countx = 0
        county = 0
        xtiles = np.linspace(xlim[0], xlim[1], nNx)
        ytiles = np.linspace(ylim[0], ylim[1], nNy)

        # Remove tiles without points in
        filt = np.ones((nNx-1)*(nNy-1), dtype='bool')

        for ii in range(xtiles.shape[0]-1):
            for jj in range(ytiles.shape[0]-1):
                # Mask along x axis
                maskx = np.all([xyLocs[:, 0] >= xtiles[ii],
                               xyLocs[:, 0] <= xtiles[int(ii+1)]], axis=0)

                # Mask along y axis
                masky = np.all([xyLocs[:, 1] >= ytiles[jj],
                               xyLocs[:, 1] <= ytiles[int(jj+1)]], axis=0)

                # Remember which axis has the most points (for next split)
                countx = np.max([np.sum(maskx), countx])
                county = np.max([np.sum(masky), county])

                # Full mask
                mask = np.all([maskx, masky], axis=0)
                nObs = np.max([nObs, np.sum(mask)])

                # Check if at least one point is picked
                if np.sum(mask) == 0:
                    filt[jj + ii*(nNy-1)] = False

    x1, x2 = xtiles[:-1], xtiles[1:]
    y1, y2 = ytiles[:-1], ytiles[1:]

    X1, Y1 = np.meshgrid(x1, y1)
    xy1 = np.c_[
        X1.flatten(order="F")[filt] - overlap[0],
        Y1.flatten(order="F")[filt] - overlap[1]
    ]

    X2, Y2 = np.meshgrid(x2, y2)
    xy2 = np.c_[
        X2.flatten(order="F")[filt] + overlap[0],
        Y2.flatten(order="F")[filt] + overlap[1]
    ]

    return [xy1, xy2]


def downsample_xy(locations, radius):
    """
    downsample_xy(locations)

    Function to downsample a cloud of points in 2D based on
    distance between neighbours

    Parameter
    ---------

    locations: numpy.ndarray
        Point locations [nx2]

    radius: float
        Minimum radial distance between points


    Return
    ------

    index: bool
        Array of bool of shape n for points to stay

    """

    tree = cKDTree(locations[:, :2])

    nstn = locations.shape[0]
    # Initialize the filter
    index = np.ones(nstn, dtype='bool')

    count = -1
    print("Begin filtering for radius= " + str(radius))

    for ii in tqdm(range(nstn)):

        if index[ii]:

            ind = tree.query_ball_point(locations[ii, :2], radius)

            index[ind] = False
            index[ii] = True

        # count = progress(ii, count, nstn)

    return index


def progress(iter, prog, final):
    """
    progress(iter,prog,final)

    Function measuring the progress of a process and print to screen the %.
    Useful to estimate the remaining runtime of a large problem.

    Parameters
    ----------

    iter: int
        Current interation

    prog: float
        Fractional progress

    final: int
        Final interation count


    """
    arg = np.floor(float(iter)/float(final)*10.)

    if arg > prog:

        print("Done " + str(arg*10) + " %")
        prog = arg

    return prog
