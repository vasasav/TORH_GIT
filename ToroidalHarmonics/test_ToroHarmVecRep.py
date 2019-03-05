# Written by Vassili Savinov on 04/03/2019
# step by step test all the parts of the ToroHarmVecRep

import unittest
from ToroHarmVecRep import ToroHarmVecRep
import numpy as np
import numpy.random as npr
import pandas as pd

class test_ToroHarmVecRep(unittest.TestCase):

    # convert toroidal coordinates to cartesian coordinates
    # this is implemented already in ToroCoords class
    # but I want to keep the testing implementation independent
    @staticmethod
    def toro_to_cart(ETA, THETA, PHI, a_val=1):
        # generate the Cartesian coordinates : Moon & Spencer "Field theory handbook" p 112
        X = a_val * (np.sinh(ETA) * np.cos(PHI)) / (np.cosh(ETA) - np.cos(THETA))
        Y = a_val * (np.sinh(ETA) * np.sin(PHI)) / (np.cosh(ETA) - np.cos(THETA))
        Z = a_val * (np.sin(THETA)) / (np.cosh(ETA) - np.cos(THETA))

        return X, Y, Z

    # the constructor
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)# initialize parent class
        # prepare testing variables

    """
    test the psi tensors generated by the program
    this include the basic psi functions as well as their derivatives up to second order
    the idea is to use Mathematica as an alternative to computing these harmonics
    """
    def test_psi_tensors(self, try_count = 50, dim_size = 3, eta_max = 13, a_val = 1.0, rtol = 1e-6):
        # load the datat computed by mathematica
        col_names = ['eta', 'theta', 'phi', 'n', 'm', 'psi_first_re', 'psi_first_im', 'psi_second_re', 'psi_second_im',
                     'dpsi_deta_first_re', 'dpsi_deta_first_im', 'dpsi_deta_second_re', 'dpsi_deta_second_im',
                     'dpsi_dtheta_first_re', 'dpsi_dtheta_first_im', 'dpsi_dtheta_second_re', 'dpsi_dtheta_second_im',
                     'd2psi_deta2_first_re', 'd2psi_deta2_first_im', 'd2psi_deta2_second_re', 'd2psi_deta2_second_im',
                     'd2psi_dtheta2_first_re', 'd2psi_dtheta2_first_im', 'd2psi_dtheta2_second_re', 'd2psi_dtheta2_second_im',
                     'd2psi_detatheta_first_re', 'd2psi_detatheta_first_im', 'd2psi_detatheta_second_re', 'd2psi_detatheta_second_im']
        #
        validData = pd.read_csv('ToroidalHarmonicsDefinition_Convention\\torHarm.csv', names=col_names)

        # load data
        etaVec = np.array(validData['eta'])
        thetaVec = np.array(validData['theta'])
        phiVec = np.array(validData['phi'])

        nVec = np.array(validData['n'], dtype=np.uint)
        mVec = np.array(validData['m'], dtype=np.uint)

        # get cartesian versions
        xVec, yVec, zVec = test_ToroHarmVecRep.toro_to_cart(etaVec, thetaVec, phiVec, a_val=a_val)

        # prepare the tensors
        toro_vec_rep = ToroHarmVecRep(xVec, yVec, zVec, 0 * xVec, 0 * xVec, 0 * xVec, nCount=np.max(nVec)+1, mCount=np.max(mVec)+1)

        # now compare the psi_tensor
        # computed by Python
        comp_psi_first = np.array(
            [toro_vec_rep.psiTens[0, nVec[iCell], mVec[iCell], iCell] for iCell in range(len(xVec))])
        comp_psi_second = np.array(
            [toro_vec_rep.psiTens[1, nVec[iCell], mVec[iCell], iCell] for iCell in range(len(xVec))])
        # tagret
        tgt_psi_first_re = np.array(validData['psi_first_re'])
        tgt_psi_first_im = np.array(validData['psi_first_im'])
        tgt_psi_second_re = np.array(validData['psi_second_re'])
        tgt_psi_second_im = np.array(validData['psi_second_im'])
        #
        self.assertTrue( np.all(np.isclose(np.real(comp_psi_first), tgt_psi_first_re, rtol=rtol)) )
        self.assertTrue( np.all(np.isclose(np.imag(comp_psi_first), tgt_psi_first_im, rtol=rtol)) )
        self.assertTrue( np.all(np.isclose(np.real(comp_psi_second), tgt_psi_second_re, rtol=rtol)) )
        self.assertTrue( np.all(np.isclose(np.imag(comp_psi_second), tgt_psi_second_im, rtol=rtol)) )

        # next the derivatives of the psi tensor dpsi_deta
        comp_dpsi_deta_first = np.array(
            [toro_vec_rep.D_psi_D_eta_Tens[0, nVec[iCell], mVec[iCell], iCell] for iCell in range(len(xVec))])
        comp_dpsi_deta_second = np.array(
            [toro_vec_rep.D_psi_D_eta_Tens[1, nVec[iCell], mVec[iCell], iCell] for iCell in range(len(xVec))])
        # tagret
        tgt_dpsi_deta_first_re = np.array(validData['dpsi_deta_first_re'])
        tgt_dpsi_deta_first_im = np.array(validData['dpsi_deta_first_im'])
        tgt_dpsi_deta_second_re = np.array(validData['dpsi_deta_second_re'])
        tgt_dpsi_deta_second_im = np.array(validData['dpsi_deta_second_im'])
        #
        self.assertTrue(np.all(np.isclose(np.real(comp_dpsi_deta_first), tgt_dpsi_deta_first_re, rtol=rtol)))
        self.assertTrue(np.all(np.isclose(np.imag(comp_dpsi_deta_first), tgt_dpsi_deta_first_im, rtol=rtol)))
        self.assertTrue(np.all(np.isclose(np.real(comp_dpsi_deta_second), tgt_dpsi_deta_second_re, rtol=rtol)))
        self.assertTrue(np.all(np.isclose(np.imag(comp_dpsi_deta_second), tgt_dpsi_deta_second_im, rtol=rtol)))

        # next the derivatives of the psi tensor dpsi_dtheta
        comp_dpsi_dtheta_first = np.array(
            [toro_vec_rep.D_psi_D_theta_Tens[0, nVec[iCell], mVec[iCell], iCell] for iCell in range(len(xVec))])
        comp_dpsi_dtheta_second = np.array(
            [toro_vec_rep.D_psi_D_theta_Tens[1, nVec[iCell], mVec[iCell], iCell] for iCell in range(len(xVec))])
        # tagret
        tgt_dpsi_dtheta_first_re = np.array(validData['dpsi_dtheta_first_re'])
        tgt_dpsi_dtheta_first_im = np.array(validData['dpsi_dtheta_first_im'])
        tgt_dpsi_dtheta_second_re = np.array(validData['dpsi_dtheta_second_re'])
        tgt_dpsi_dtheta_second_im = np.array(validData['dpsi_dtheta_second_im'])
        #
        self.assertTrue(np.all(np.isclose(np.real(comp_dpsi_dtheta_first), tgt_dpsi_dtheta_first_re, rtol=rtol)))
        self.assertTrue(np.all(np.isclose(np.imag(comp_dpsi_dtheta_first), tgt_dpsi_dtheta_first_im, rtol=rtol)))
        self.assertTrue(np.all(np.isclose(np.real(comp_dpsi_dtheta_second), tgt_dpsi_dtheta_second_re, rtol=rtol)))
        self.assertTrue(np.all(np.isclose(np.imag(comp_dpsi_dtheta_second), tgt_dpsi_dtheta_second_im, rtol=rtol)))

        # next the derivatives of the psi tensor d2psi_deta2
        comp_d2psi_deta2_first = np.array(
            [toro_vec_rep.D2_psi_D_eta2_Tens[0, nVec[iCell], mVec[iCell], iCell] for iCell in range(len(xVec))])
        comp_d2psi_deta2_second = np.array(
            [toro_vec_rep.D2_psi_D_eta2_Tens[1, nVec[iCell], mVec[iCell], iCell] for iCell in range(len(xVec))])
        # tagret
        tgt_d2psi_deta2_first_re = np.array(validData['d2psi_deta2_first_re'])
        tgt_d2psi_deta2_first_im = np.array(validData['d2psi_deta2_first_im'])
        tgt_d2psi_deta2_second_re = np.array(validData['d2psi_deta2_second_re'])
        tgt_d2psi_deta2_second_im = np.array(validData['d2psi_deta2_second_im'])
        #
        self.assertTrue(np.all(np.isclose(np.real(comp_d2psi_deta2_first), tgt_d2psi_deta2_first_re, rtol=rtol)))
        self.assertTrue(np.all(np.isclose(np.imag(comp_d2psi_deta2_first), tgt_d2psi_deta2_first_im, rtol=rtol)))
        self.assertTrue(np.all(np.isclose(np.real(comp_d2psi_deta2_second), tgt_d2psi_deta2_second_re, rtol=rtol)))
        self.assertTrue(np.all(np.isclose(np.imag(comp_d2psi_deta2_second), tgt_d2psi_deta2_second_im, rtol=rtol)))

        # next the derivatives of the psi tensor d2psi_dtheta2
        comp_d2psi_dtheta2_first = np.array(
            [toro_vec_rep.D2_psi_D_theta2_Tens[0, nVec[iCell], mVec[iCell], iCell] for iCell in range(len(xVec))])
        comp_d2psi_dtheta2_second = np.array(
            [toro_vec_rep.D2_psi_D_theta2_Tens[1, nVec[iCell], mVec[iCell], iCell] for iCell in range(len(xVec))])
        # tagret
        tgt_d2psi_dtheta2_first_re = np.array(validData['d2psi_dtheta2_first_re'])
        tgt_d2psi_dtheta2_first_im = np.array(validData['d2psi_dtheta2_first_im'])
        tgt_d2psi_dtheta2_second_re = np.array(validData['d2psi_dtheta2_second_re'])
        tgt_d2psi_dtheta2_second_im = np.array(validData['d2psi_dtheta2_second_im'])
        #
        self.assertTrue(np.all(np.isclose(np.real(comp_d2psi_dtheta2_first), tgt_d2psi_dtheta2_first_re, rtol=rtol)))
        self.assertTrue(np.all(np.isclose(np.imag(comp_d2psi_dtheta2_first), tgt_d2psi_dtheta2_first_im, rtol=rtol)))
        self.assertTrue(np.all(np.isclose(np.real(comp_d2psi_dtheta2_second), tgt_d2psi_dtheta2_second_re, rtol=rtol)))
        self.assertTrue(np.all(np.isclose(np.imag(comp_d2psi_dtheta2_second), tgt_d2psi_dtheta2_second_im, rtol=rtol)))

        # next the derivatives of the psi tensor d2psi_detatheta
        comp_d2psi_detatheta_first = np.array(
            [toro_vec_rep.D2_psi_D_eta_theta_Tens[0, nVec[iCell], mVec[iCell], iCell] for iCell in range(len(xVec))])
        comp_d2psi_detatheta_second = np.array(
            [toro_vec_rep.D2_psi_D_eta_theta_Tens[1, nVec[iCell], mVec[iCell], iCell] for iCell in range(len(xVec))])
        # tagret
        tgt_d2psi_detatheta_first_re = np.array(validData['d2psi_detatheta_first_re'])
        tgt_d2psi_detatheta_first_im = np.array(validData['d2psi_detatheta_first_im'])
        tgt_d2psi_detatheta_second_re = np.array(validData['d2psi_detatheta_second_re'])
        tgt_d2psi_detatheta_second_im = np.array(validData['d2psi_detatheta_second_im'])
        #
        self.assertTrue(np.all(np.isclose(np.real(comp_d2psi_detatheta_first), tgt_d2psi_detatheta_first_re, rtol=rtol)))
        self.assertTrue(np.all(np.isclose(np.imag(comp_d2psi_detatheta_first), tgt_d2psi_detatheta_first_im, rtol=rtol)))
        self.assertTrue(np.all(np.isclose(np.real(comp_d2psi_detatheta_second), tgt_d2psi_detatheta_second_re, rtol=rtol)))
        self.assertTrue(np.all(np.isclose(np.imag(comp_d2psi_detatheta_second), tgt_d2psi_detatheta_second_im, rtol=rtol)))

    """
    One of the first things to be done in the test_ToroHarmVecRep
    is to convert from Cartesian coordinates to toroidal coordinates
    
    test that this works as it should
    
    try_count = number of time I will try to test
    dim_size = number of the single dimension for the coordinate cubes
    eta_max = maximum size for eta
    a_val the scale for toroidal coordinates
    rtol the relative tolrance to use in comparisons
    """
    def test_toro_coord(self, try_count = 5, dim_size = 2, eta_max = 13, a_val = 1.0, rtol = 1e-6):

        for iTry in range(try_count):
            # generate the positions in toroidal coordinates
            ETA = npr.rand(dim_size, dim_size, dim_size) * eta_max + 1e-9 # 0...eta_max
            THETA = (npr.rand(dim_size, dim_size, dim_size) * 2 * np.pi) - np.pi # -pi...pi
            PHI = npr.rand(dim_size, dim_size, dim_size) * 2 * np.pi # 0 ... 2pi

            X, Y, Z = test_ToroHarmVecRep.toro_to_cart(ETA, THETA, PHI, a_val=a_val)

            # now get the decomposition
            toro_vec_rep = ToroHarmVecRep(X, Y, Z, 0 * X, 0 * X, 0 * X)

            # compare the new values, use cos/sin to avoid phase unwrapping on angles
            self.assertTrue(np.all(np.isclose(toro_vec_rep.etaVec, ETA.flatten(), rtol=rtol)))
            # use of sin and cos allows to avoid unwrapping issues
            self.assertTrue(np.all(np.isclose( np.sin(toro_vec_rep.thetaVec), np.sin(THETA.flatten()), rtol=rtol)))
            self.assertTrue(np.all(np.isclose( np.cos(toro_vec_rep.thetaVec), np.cos(THETA.flatten()), rtol=rtol)))
            # use of sin and cos allows to avoid unwrapping issues
            self.assertTrue(np.all(np.isclose( np.sin(toro_vec_rep.phiVec), np.sin(PHI.flatten()), rtol=rtol)))
            self.assertTrue(np.all(np.isclose( np.cos(toro_vec_rep.phiVec), np.cos(PHI.flatten()), rtol=rtol)))

            # also within this time we compute the r2 within the initialization
            self.assertTrue(np.all(np.isclose(toro_vec_rep.r2Vec, (X**2+Y**2+Z**2).flatten(), rtol=rtol)))

            # also computing the Jacobian for the toroidal coordinates, i.e.
            # \frac{\partial(xyz)}{\partial(\eta\theta\phi)}
            # see the corresponding Mathematica file
            jacobian = a_val**3 * np.sinh(ETA)/(np.cosh(ETA)-np.cos(THETA))**3
            self.assertTrue(np.all(np.isclose(toro_vec_rep.jacobian, jacobian.flatten(), rtol=rtol)))

            print('test_toro_coord run = %d' % iTry)


