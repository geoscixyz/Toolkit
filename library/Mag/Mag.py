from . import MathUtils
from scipy.constants import mu_0
import re
import numpy as np
from SimPEG import Utils, PF
from SimPEG.PF import BaseMag


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
            self._Higrf = self.survey.srcField.param[0]

        return self._Higrf * 1e-9 / mu_0

    @property
    def Hinc(self):

        if getattr(self, '_Hinc', None) is None:
            self._Hinc = self.survey.srcField.param[1]

        return self._Hinc

    @property
    def Hdec(self):

        if getattr(self, '_Hdec', None) is None:
            self._Hdec = self.survey.srcField.param[2]

        return self._Hdec

    @property
    def G(self):

        if getattr(self, '_G', None) is None:

            rxLoc = self.survey.srcField.rxList[0].locs

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
            u = Utils.mkvc(bvec[0, :])

        if self.uType == 'by':
            u = Utils.mkvc(bvec[1, :])

        if self.uType == 'bz':
            u = Utils.mkvc(bvec[2, :])

        if self.uType == 'tf':
            # Projection matrix
            Ptmi = MathUtils.dipazm_2_xyz(self.Hinc,
                                         self.Hdec)

            u = Utils.mkvc(Ptmi.dot(bvec))

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

    Yn = np.c_[Utils.mkvc(yn1), Utils.mkvc(yn2)]
    Xn = np.c_[Utils.mkvc(xn1), Utils.mkvc(xn2)]
    Zn = np.c_[Utils.mkvc(zn1), Utils.mkvc(zn2)]

    ndata = rxLoc.shape[0]

    # Pre-allocate forward matrix
    G = np.zeros((int(3*ndata), 3))

    for ii in range(ndata):

        tx, ty, tz = PF.Magnetics.get_T_mat(Xn, Yn, Zn, rxLoc[ii, :])

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

    rxLoc = BaseMag.RxObs(xyz)
    srcField = BaseMag.SrcField([rxLoc], param=EarthField)
    survey = BaseMag.LinearSurvey(srcField)

    if data is not None:
        survey.dobs = data
    else:
        survey.dobs = np.zeros(xyz.shape[0])

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

        rxLoc = BaseMag.RxObs(locXYZ)
        srcField = BaseMag.SrcField([rxLoc], param=(B[2], B[0], B[1]))
        survey = BaseMag.LinearSurvey(srcField)
        survey.dobs = d
        survey.std = wd
        return survey
