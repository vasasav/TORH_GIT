# Written by Vassili Savinov on 05/02/2019
# this class will be used to handle toroidal coordinates
# note that toroidal coordinates contain a single length parameter that scales the coordinates

# using Mathematica and Moon & Spencer "Field theory handbook"

import numpy as np

class ToroCoords:

    # initialize
    def __init__(self, aParam = 1.0):

        self.aParam = aParam # the scale for toroidal coordinates

    # convert toroidal coordinates to cartesian
    def Toro_to_Cart(self, etaVals, thetaVals, phiVals):
        ########## check sanity of the inputs
        if not (isinstance(etaVals, np.ndarray) and np.isreal(etaVals).all() and
                isinstance(thetaVals, np.ndarray) and np.isreal(thetaVals).all() and
                isinstance(phiVals, np.ndarray) and np.isreal(phiVals).all()):
            raise ValueError('ToroCoords.Toro_to_Cart: xVals, yVals, zVals have to be real numpy arrays')

        if (etaVals.shape != thetaVals.shape) or (etaVals.shape != phiVals.shape):
            raise ValueError('ToroCoords.Toro_to_Cart: xVals, yVals, zVals have to be of the same shape')

        preFac = self.aParam / (np.cosh(etaVals) - np.cos(thetaVals))

        xVals = preFac * np.cos(phiVals) * np.sinh(etaVals)
        yVals = preFac * np.sin(phiVals) * np.sinh(etaVals)
        zVals = preFac * np.sin(thetaVals)

        return (xVals, yVals, zVals)

    # convert cartesian coordinates to toroidal coordinates
    # accept xVals, yVals, zVals of any shape as long as it is the same shape
    def Cart_to_Toro(self, xVals, yVals, zVals):
        ########## check sanity of the inputs
        if not (    isinstance(xVals, np.ndarray) and np.isreal(xVals).all() and
                    isinstance(yVals, np.ndarray) and np.isreal(yVals).all() and
                    isinstance(zVals, np.ndarray) and np.isreal(zVals).all()):
            raise ValueError('ToroCoords.Cart_to_Toro: xVals, yVals, zVals have to be real numpy arrays')

        if (xVals.shape != yVals.shape) or (xVals.shape != zVals.shape):
            raise ValueError('ToroCoords.Cart_to_Toro: xVals, yVals, zVals have to be of the same shape')

        # now compute the coordintes
        # Moon & Spences p 112
        #{\[CapitalEta], \[CapitalTheta], \[CapitalPhi]} = FullSimplify[CoordinateTransform[ "Cartesian" -> "Toroidal", {X, Y, Z}]](**)
        # and
        etaVals = np.arctanh( (2 * self.aParam * np.sqrt(xVals**2 + yVals**2))/(xVals**2 + yVals**2 + zVals**2 + self.aParam**2) )
        #
        phiVals = np.arctan2(yVals, xVals)

        # now work out theta
        thetaVals = np.arctan2(2.0 * self.aParam * zVals, (xVals**2 + yVals**2 + zVals**2-self.aParam**2))

        return (etaVals, thetaVals, phiVals)