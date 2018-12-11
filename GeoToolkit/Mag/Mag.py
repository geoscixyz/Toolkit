from . import MathUtils
from scipy.constants import mu_0
import re
import numpy as np
# from SimPEG import Utils, PF
# from SimPEG.PF import BaseMag


class Problem(object):
    """
            Earth's field:
            - Binc, Bdec : inclination and declination of Earth's mag field
            - Bigrf : amplitude of earth's field in units of nT

        Remnance:
            - Q : Koenigsberger ratio
            - Rinc, Rdec : inclination and declination of remnance in block

    """
    #Bdec, Binc, Bigrf = 90., 0., 50000.
    Q, rinc, rdec = 0., 0., 0.
    uType, mType = 'tf', 'induced'
    susc = 1.
    prism = None
    survey = None

    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)

        return

    @property
    def Mind(self):
        # Define magnetization direction as sum of induced and remanence
        mind = MathUtils.dipazm_2_xyz(self.Hinc, self.Hdec)
        R = MathUtils.rotationMatrix(-self.prism.pinc, -self.prism.pdec, normal=False)
        Mind = self.susc*self.Higrf*R.dot(mind.T)
        # Mind = self.susc*self.Higrf*PF.Magnetics.dipazm_2_xyz(self.Binc - self.prism.pinc,
        #                                                self.Bdec - self.prism.pdec)
        return Mind

    @property
    def Mrem(self):

        mrem = MathUtils.dipazm_2_xyz(self.rinc, self.rdec)
        R = MathUtils.rotationMatrix(-self.prism.pinc, -self.prism.pdec, normal=False)
        Mrem = self.Q*self.susc*self.Higrf * R.dot(mrem.T)

        return Mrem

    @property
    def Higrf(self):

        if getattr(self, '_Higrf', None) is None:
            self._Higrf = self.survey.srcFieldParam[0]

        return self._Higrf * 1e-9 / mu_0

    @property
    def Hinc(self):

        if getattr(self, '_Hinc', None) is None:
            self._Hinc = self.survey.srcFieldParam[1]

        return self._Hinc

    @property
    def Hdec(self):

        if getattr(self, '_Hdec', None) is None:
            self._Hdec = self.survey.srcFieldParam[2]

        return self._Hdec

    @property
    def G(self):

        if getattr(self, '_G', None) is None:

            rxLoc = self.survey.rxLoc

            xLoc = rxLoc[:, 0] - self.prism.xc
            yLoc = rxLoc[:, 1] - self.prism.yc
            zLoc = rxLoc[:, 2] - self.prism.zc

            R = MathUtils.rotationMatrix(-self.prism.pinc, -self.prism.pdec, normal=False)

            rxLoc = R.dot(np.c_[xLoc, yLoc, zLoc].T).T

            rxLoc = np.c_[rxLoc[:, 0] + self.prism.xc, rxLoc[:, 1] + self.prism.yc, rxLoc[:, 2] + self.prism.zc]

            # Create the linear forward system
            self._G = Intrgl_Fwr_Op(self.prism.xn, self.prism.yn, self.prism.zn, rxLoc)

        return self._G

    def fields(self):

        if (self.mType == 'induced') or (self.mType == 'total'):

            b = self.G.dot(self.Mind)
            self.fieldi = self.extractFields(b)

        if (self.mType == 'remanent') or (self.mType == 'total'):

            b = self.G.dot(self.Mrem)

            self.fieldr = self.extractFields(b)

        if self.mType == 'induced':
            return [self.fieldi]
        elif self.mType == 'remanent':
            return [self.fieldr]
        elif self.mType == 'total':
            return [self.fieldi, self.fieldr]

    def extractFields(self, bvec):

        nD = int(bvec.shape[0]/3)
        bvec = np.reshape(bvec, (3, nD))

        R = MathUtils.rotationMatrix(self.prism.pinc, self.prism.pdec)
        bvec = R.dot(bvec)

        if self.uType == 'bx':
            u = bvec[0, :].flatten()

        if self.uType == 'by':
            u = bvec[1, :].flatten()

        if self.uType == 'bz':
            u = bvec[2, :].flatten()

        if self.uType == 'tf':
            # Projection matrix
            Ptmi = MathUtils.dipazm_2_xyz(self.Hinc,
                                         self.Hdec)

            u = Ptmi.dot(bvec).flatten()

        return u


def Intrgl_Fwr_Op(xn, yn, zn, rxLoc):

    """

    Magnetic forward operator in integral form

    flag        = 'ind' | 'full'

      1- ind : Magnetization fixed by user

      3- full: Full tensor matrix stored with shape([3*ndata, 3*nc])

    Return
    _G = Linear forward modeling operation

     """

    yn2, xn2, zn2 = np.meshgrid(yn[1:], xn[1:], zn[1:])
    yn1, xn1, zn1 = np.meshgrid(yn[0:-1], xn[0:-1], zn[0:-1])

    Yn = np.c_[yn1.flatten(), yn2.flatten()]
    Xn = np.c_[xn1.flatten(), xn2.flatten()]
    Zn = np.c_[zn1.flatten(), zn2.flatten()]

    ndata = rxLoc.shape[0]

    # Pre-allocate forward matrix
    G = np.zeros((int(3*ndata), 3))

    for ii in range(ndata):

        tx, ty, tz = calcRow(Xn, Yn, Zn, rxLoc[ii, :])

        G[ii, :] = tx / 1e-9 * mu_0
        G[ii+ndata, :] = ty / 1e-9 * mu_0
        G[ii+2*ndata, :] = tz / 1e-9 * mu_0

    return G


def createMagSurvey(xyz, EarthField=np.r_[50000, 90, 0], data=None):
    """
        Create SimPEG magnetic survey pbject

        INPUT
        :param array: xyz, n-by-4 array of observation locations
        :param array: EarthField [Default 50000,90,0], 1-by-3 array of inducing field param [|B|, Inc, Dec]

        OPTIONAL
        :param array: data, n-by-4 array of data
    """

    # rxLoc = BaseMag.RxObs(xyz)
    # srcField = BaseMag.SrcField([rxLoc], param=EarthField)


    survey = Survey(EarthField)
    survey._rxLoc = xyz
    if data is not None:
        survey._dobs = data
    else:
        survey._dobs = np.zeros(xyz.shape[0])

    return survey


def readMagneticsObservations(obs_file):
        """
            Read and write UBC mag file format

            INPUT:
            :param fileName, path to the UBC obs mag file

            OUTPUT:
            :param survey
            :param M, magnetization orentiaton (MI, MD)
        """

        fid = open(obs_file, 'r')

        # First line has the inclination,declination and amplitude of B0
        line = fid.readline()
        B = np.array(line.split(), dtype=float)

        # Second line has the magnetization orientation and a flag
        line = fid.readline()
        M = np.array(line.split(), dtype=float)

        # Third line has the number of rows
        line = fid.readline()
        ndat = int(line.strip())

        # Pre-allocate space for obsx, obsy, obsz, data, uncert
        line = fid.readline()
        temp = np.array(line.split(), dtype=float)

        d = np.zeros(ndat, dtype=float)
        wd = np.zeros(ndat, dtype=float)
        locXYZ = np.zeros((ndat, 3), dtype=float)

        for ii in range(ndat):

            temp = np.array(line.split(), dtype=float)
            locXYZ[ii, :] = temp[:3]

            if len(temp) > 3:
                d[ii] = temp[3]

                if len(temp) == 5:
                    wd[ii] = temp[4]

            line = fid.readline()

        # rxLoc = BaseMag.RxObs(locXYZ)
        # srcField = BaseMag.SrcField([rxLoc], param=(B[2], B[0], B[1]))
        survey = Survey(np.r_[B[2], B[0], B[1]])
        survey._rxLoc = locXYZ
        survey._dobs = d
        survey._std = wd
        return survey


def calcRow(Xn, Yn, Zn, rxLoc):
    """
    Load in the active nodes of a tensor mesh and computes the magnetic tensor
    for a given observation location rxLoc[obsx, obsy, obsz]

    INPUT:
    Xn, Yn, Zn: Node location matrix for the lower and upper most corners of
                all cells in the mesh shape[nC,2]
    M
    OUTPUT:
    Tx = [Txx Txy Txz]
    Ty = [Tyx Tyy Tyz]
    Tz = [Tzx Tzy Tzz]

    where each elements have dimension 1-by-nC.
    Only the upper half 5 elements have to be computed since symetric.
    Currently done as for-loops but will eventually be changed to vector
    indexing, once the topography has been figured out.

    Created on Oct, 20th 2015

    @author: dominiquef

     """

    eps = 1e-8  # add a small value to the locations to avoid /0

    nC = Xn.shape[0]

    # Pre-allocate space for 1D array
    Tx = np.zeros((1, 3*nC))
    Ty = np.zeros((1, 3*nC))
    Tz = np.zeros((1, 3*nC))

    dz2 = Zn[:, 1] - rxLoc[2] + eps
    dz1 = Zn[:, 0] - rxLoc[2] + eps

    dy2 = Yn[:, 1] - rxLoc[1] + eps
    dy1 = Yn[:, 0] - rxLoc[1] + eps

    dx2 = Xn[:, 1] - rxLoc[0] + eps
    dx1 = Xn[:, 0] - rxLoc[0] + eps

    dx2dx2 = dx2**2.
    dx1dx1 = dx1**2.

    dy2dy2 = dy2**2.
    dy1dy1 = dy1**2.

    dz2dz2 = dz2**2.
    dz1dz1 = dz1**2.

    R1 = (dy2dy2 + dx2dx2)
    R2 = (dy2dy2 + dx1dx1)
    R3 = (dy1dy1 + dx2dx2)
    R4 = (dy1dy1 + dx1dx1)

    arg1 = np.sqrt(dz2dz2 + R2)
    arg2 = np.sqrt(dz2dz2 + R1)
    arg3 = np.sqrt(dz1dz1 + R1)
    arg4 = np.sqrt(dz1dz1 + R2)
    arg5 = np.sqrt(dz2dz2 + R3)
    arg6 = np.sqrt(dz2dz2 + R4)
    arg7 = np.sqrt(dz1dz1 + R4)
    arg8 = np.sqrt(dz1dz1 + R3)

    Tx[0, 0:nC] = (
        np.arctan2(dy1 * dz2, (dx2 * arg5 + eps)) -
        np.arctan2(dy2 * dz2, (dx2 * arg2 + eps)) +
        np.arctan2(dy2 * dz1, (dx2 * arg3 + eps)) -
        np.arctan2(dy1 * dz1, (dx2 * arg8 + eps)) +
        np.arctan2(dy2 * dz2, (dx1 * arg1 + eps)) -
        np.arctan2(dy1 * dz2, (dx1 * arg6 + eps)) +
        np.arctan2(dy1 * dz1, (dx1 * arg7 + eps)) -
        np.arctan2(dy2 * dz1, (dx1 * arg4 + eps))
    )

    Ty[0, 0:nC] = (
        np.log((dz2 + arg2 + eps) / (dz1 + arg3 + eps)) -
        np.log((dz2 + arg1 + eps) / (dz1 + arg4 + eps)) +
        np.log((dz2 + arg6 + eps) / (dz1 + arg7 + eps)) -
        np.log((dz2 + arg5 + eps) / (dz1 + arg8 + eps))
    )

    Ty[0, nC:2*nC] = (
        np.arctan2(dx1 * dz2, (dy2 * arg1 + eps)) -
        np.arctan2(dx2 * dz2, (dy2 * arg2 + eps)) +
        np.arctan2(dx2 * dz1, (dy2 * arg3 + eps)) -
        np.arctan2(dx1 * dz1, (dy2 * arg4 + eps)) +
        np.arctan2(dx2 * dz2, (dy1 * arg5 + eps)) -
        np.arctan2(dx1 * dz2, (dy1 * arg6 + eps)) +
        np.arctan2(dx1 * dz1, (dy1 * arg7 + eps)) -
        np.arctan2(dx2 * dz1, (dy1 * arg8 + eps))
    )

    R1 = (dy2dy2 + dz1dz1)
    R2 = (dy2dy2 + dz2dz2)
    R3 = (dy1dy1 + dz1dz1)
    R4 = (dy1dy1 + dz2dz2)

    Ty[0, 2*nC:] = (
        np.log((dx1 + np.sqrt(dx1dx1 + R1) + eps) /
               (dx2 + np.sqrt(dx2dx2 + R1) + eps)) -
        np.log((dx1 + np.sqrt(dx1dx1 + R2) + eps) /
               (dx2 + np.sqrt(dx2dx2 + R2) + eps)) +
        np.log((dx1 + np.sqrt(dx1dx1 + R4) + eps) /
               (dx2 + np.sqrt(dx2dx2 + R4) + eps)) -
        np.log((dx1 + np.sqrt(dx1dx1 + R3) + eps) /
               (dx2 + np.sqrt(dx2dx2 + R3) + eps))
    )

    R1 = (dx2dx2 + dz1dz1)
    R2 = (dx2dx2 + dz2dz2)
    R3 = (dx1dx1 + dz1dz1)
    R4 = (dx1dx1 + dz2dz2)

    Tx[0, 2*nC:] = (
        np.log((dy1 + np.sqrt(dy1dy1 + R1) + eps) /
               (dy2 + np.sqrt(dy2dy2 + R1) + eps)) -
        np.log((dy1 + np.sqrt(dy1dy1 + R2) + eps) /
               (dy2 + np.sqrt(dy2dy2 + R2) + eps)) +
        np.log((dy1 + np.sqrt(dy1dy1 + R4) + eps) /
               (dy2 + np.sqrt(dy2dy2 + R4) + eps)) -
        np.log((dy1 + np.sqrt(dy1dy1 + R3) + eps) /
               (dy2 + np.sqrt(dy2dy2 + R3) + eps))
    )

    Tz[0, 2*nC:] = -(Ty[0, nC:2*nC] + Tx[0, 0:nC])
    Tz[0, nC:2*nC] = Ty[0, 2*nC:]
    Tx[0, nC:2*nC] = Ty[0, 0:nC]
    Tz[0, 0:nC] = Tx[0, 2*nC:]

    Tx = Tx/(4*np.pi)
    Ty = Ty/(4*np.pi)
    Tz = Tz/(4*np.pi)

    return Tx, Ty, Tz


class Survey():
    """Base Magnetics Survey"""

    rxType = 'tmi'  #: receiver type

    def __init__(self, srcField, **kwargs):
        self._srcFieldParam = srcField
        # Survey.BaseSurvey.__init__(self, **kwargs)

    def eval(self, u):
        return u

    @property
    def nD(self):
        if getattr(self, '_rxLoc', None) is not None:
            return self.rxLoc.shape[0]
        else:
            return None

    @property
    def rxLoc(self):

        if getattr(self, '_rxLoc', None) is not None:
            return self._rxLoc
        else:
            return None

    @property
    def srcFieldParam(self):
        if getattr(self, '_srcFieldParam', None) is not None:
            return self._srcFieldParam
        else:
            return None

    @property
    def dobs(self):
        if getattr(self, '_dobs', None) is not None:
            return self._dobs
        else:
            return None

    @property
    def std(self):
        if getattr(self, '_std', None) is not None:
            return self._std
        else:
            return None
