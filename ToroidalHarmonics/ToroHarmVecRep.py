# Written by Vassili Savinov on 05/02/19
# class for containing the representation of a 3d vector
# field in terms of toroidal harmonics
# see Toroidal_Harmonics_Def.pdf

import numpy as np
import numpy.linalg as npla
from DTORH import DTORH
from ToroCoords import ToroCoords

class ToroHarmVecRep:

    # build a rank 6 tensor lnms
    # l - type (P,Q)
    # n - legendre n
    # m - legendre m
    # s - all the spatial values flattened in an array
    def __prep_contraction_tens(self):

        # prepare new nm tensors


        # truncate pCube and qCube to right size
        #self.

        firstKind_Psi = np.sqrt(np.cosh(self.etaTens)-np.cos(self.thetaTens)) * self.pCube * np.exp()

        # vector for all the quantities
        contrPsiTens = np.zeros( (2, self.nCount, self.mCount, len(self.etaVec)), dtype=np.complex)
        #!!!!!!!!!!!!!!!

        contrPsiTens[0,:,:,:]=self.pCube[None,:,:,:]

    def __prep_psi_tens(self):
        # get the necessary legendre functions
        with DTORH() as dtorh: (pCube, qCube) = \
            dtorh.GetCubeNMZ(self.nCount + 2, self.mCount, np.cosh(self.etaVec))


        !!!!!!!!!!!!!!!! need to add exponentials etc!!!!!!!!!!

        # basic
        self.psiTens = np.zeros( (2, self.nCount, self.mCount, len(self.etaVec)), dtype=np.complex)
        self.psiTens[0, :, :, :] = pCube[0:self.nCount, :, :]#  loose tail
        self.psiTens[1, :, :, :] = qCube[0:self.nCount, :, :]  # loose tail

        # shifted n -> n+1
        self.psiTens_n_plus_1 = np.zeros( (2, self.nCount, self.mCount, len(self.etaVec)), dtype=np.complex)
        self.psiTens_n_plus_1[0, :, :, :] = pCube[1:(self.nCount+1), :, :]#  loose tail
        self.psiTens_n_plus_1[1, :, :, :] = qCube[1:(self.nCount+1), :, :]  # loose tail

        # shifted n -> n+2
        self.psiTens_n_plus_2 = np.zeros((2, self.nCount, self.mCount, len(self.etaVec)), dtype=np.complex)
        self.psiTens_n_plus_2[0, :, :, :] = pCube[2:(self.nCount + 2), :, :]  # loose tail
        self.psiTens_n_plus_2[1, :, :, :] = qCube[2:(self.nCount + 2), :, :]  # loose tail

        # need to define the basic tensors
        self.etaTens = np.tile(self.etaVec, (2, self.nCount, self.mCount, 1))
        self.thetaTens = np.tile(self.thetaVec, (2, self.nCount, self.mCount, 1))
        self.phiTens = np.tile(self.phiVec, (2, self.nCount, self.mCount, 1))
        self.jacTens = np.tile(self.jacobian, (2, self.nCount, self.mCount, 1))

        # nm tensors
        nVec = np.arange(self.nCount)
        self.nTens = np.tile(nVec[None,:,None,None], (2, 1, self.mCount, len(self.etaVec)))
        #
        mVec = np.arange(self.mCount)
        self.mTens = np.tile(mVec[None, None, :, None], (2, self.nCount, 1, len(self.etaVec)))

        # now we can work, the basic psi is defined, but let as also define the derivatives
        # first order
        self.D_psi_D_eta_tens = \
            ( np.cosh(self.etaTens)*((2*self.nTens+1)*np.cos(self.thetaTens)-2*self.nTens*np.cosh(self.etaTens)) - 1 )*self.psiTens/\
            (2 * np.sinh(self.etaTens) * (np.cosh(self.etaTens)-np.cos(self.thetaTens)) ) + \
            (self.nTens+0.5-mTens)*np.exp(-1j*self.thetaTens)*self.psiTens_n_plus_1/np.sinh(self.etaTens)




    # remove spatial points that are bad in toroidal coordinates and flatten arrays
    def __extract_clean_points(self, xTens, yTens, zTens, divFTens, tiny):
        # first of all get the toroidal coordinates
        torCrd = ToroCoords()
        etaTens, thetaTens, phiTens = torCrd.Cart_to_Toro(xTens, yTens, zTens)

        # prepare the legendre polynomials
        temp_coshEtaVec = (np.cosh(etaTens)).flatten()

        # some values are bad for numerical calculations, e.g. coshEtaVec=1.0 filter them
        # add all the filters here
        # there is no problem with dropping some values, that's why we have many!
        goodLogical = np.isfinite(temp_coshEtaVec)
        # further checks are impossible non-finite numbers, so remove them before any further logical checks
        temp_coshEtaVec[np.where(goodLogical != True)] = 0.0
        # good to contninue
        goodLogical = np.logical_and((temp_coshEtaVec - 1.0) >= tiny, goodLogical)
        self.goodInds = np.where(goodLogical)  # find indices , keep this in case any further cleanup is needed

        # store good flattened values
        self.etaVec = etaTens.flatten()[self.goodInds]
        self.thetaVec = thetaTens.flatten()[self.goodInds]
        self.phiVec = phiTens.flatten()[self.goodInds]
        self.divFVec = divFTens.flatten()[self.goodInds]

        # also store the sqrt_metric for further computations
        self.jacobian = torCrd.sqrt_metric(self.etaVec, self.thetaVec, self.phiVec)

    # to initialize itwe need is
    # three rank-3 tensors giving the XYZ coordinates of the points in 3d space
    # then we need 3 rank three tensors giving the div(F), r.F and L.F, where F is the field of interest
    # nCount and mCount is the maxiumum orders of toroidal harmonics to consider
    # tiny = a finite small number taken to be equivalent to zero
    def __init__(self, raw_xTens, raw_yTens, raw_zTens, raw_divFTens, nCount=2, mCount=3, tiny=1e-12):
        # !!!!!!!! for now only work with divergence  , rFTens, LFTens, later

        self.nCount = nCount
        self.mCount = mCount

        # convert tensors for spatial points into flattened vectors
        # and remove all the points that have bad cooridinates in toroidal coordinate system
        self.__extract_clean_points(raw_xTens, raw_yTens, raw_zTens, raw_divFTens, tiny)

        # prepare psi tensors
        self.__prep_psi_tens()

        # get the contraction tensor
        #self.__prep_contraction_tens()


