# Written by Vassili Savinov on 04/03/2019
# step by step test all the parts of the ToroHarmVecRep

import unittest
from ToroHarmVecRep import ToroHarmVecRep
import numpy as np
import numpy.random as npr

class test_ToroHarmVecRep(unittest.TestCase):

    # the constructor
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)# initialize parent class
        # prepare testing variables

    """
    One of the first things to be done in the test_ToroHarmVecRep
    is to convert from Cartesian coordinates to toroidal coordinates
    
    test that this works as it should
    
    try_count = number of time I will try to test
    dim_size = number of the single dimension for the coordinate cubes
    eta_max = maximum size for eta
    a_val the scale for toroidal coordinates
    """
    def test_toro_coord(self, try_count = 50, dim_size = 3, eta_max = 13, a_val = 1.0):

        for iTry in range(try_count):
            # generate the positions in toroidal coordinates
            ETA = npr.rand(dim_size, dim_size, dim_size) * eta_max + 1e-9 # 0...eta_max
            THETA = (npr.rand(dim_size, dim_size, dim_size) * 2 * np.pi) - np.pi # -pi...pi
            PHI = npr.rand(dim_size, dim_size, dim_size) * 2 * np.pi # 0 ... 2pi

            # generate the Cartesian coordinates : Moon & Spencer "Field theory handbook" p 112
            X = a_val * ( np.sinh(ETA) * np.cos(PHI) ) / ( np.cosh(ETA) - np.cos(THETA) )
            Y = a_val * (np.sinh(ETA) * np.sin(PHI) ) / (np.cosh(ETA) - np.cos(THETA) )
            Z = a_val * ( np.sin(THETA) ) / (np.cosh(ETA) - np.cos(THETA) )

            # now get the decomposition
            toro_vec_rep = ToroHarmVecRep(X, Y, Z, 0 * X, 0 * X, 0 * X)

            # compare the new values, use cos/sin to avoid phase unwrapping on angles
            self.assertTrue(
                np.all(
                np.isclose(toro_vec_rep.etaVec, ETA.flatten())                      &
                np.isclose(np.cos(toro_vec_rep.thetaVec), np.cos(THETA.flatten()))  &
                np.isclose(np.sin(toro_vec_rep.thetaVec), np.sin(THETA.flatten()))  &
                np.isclose(np.cos(toro_vec_rep.phiVec), np.cos(PHI.flatten()))      &
                np.isclose(np.sin(toro_vec_rep.phiVec), np.sin(PHI.flatten()))
                )
            )

            print('test_toro_coord run = %d' % iTry)

