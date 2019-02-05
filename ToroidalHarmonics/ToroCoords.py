# Written by Vassili Savinov on 05/02/2019
# this class will be used to handle toroidal coordinates
# note that toroidal coordinates contain a single length parameter that scales the coordinates

# using Mathematica and Moon & Spencer "Field theory handbook"

import numpy as np

class ToroCoords:

    # check wether the three arrays are real and of the same shape
    # could be toroidal or cartesian, hence the naming convention
    @staticmethod
    def __check_valid_coords(arr1, arr2, arr3):
        if not (isinstance(arr1, np.ndarray) and np.isreal(arr1).all() and
                isinstance(arr2, np.ndarray) and np.isreal(arr2).all() and
                isinstance(arr3, np.ndarray) and np.isreal(arr3).all()):
            raise ValueError('ToroCoords: coordinates have to be real numpy arrays')

        if (arr1.shape != arr2.shape) or (arr1.shape != arr3.shape):
            raise ValueError('ToroCoords: coordinates have to be of the same shape')

    # initialize
    def __init__(self, aParam = 1.0):

        self.aParam = aParam # the scale for toroidal coordinates

    # convert toroidal coordinates to cartesian
    def Toro_to_Cart(self, etaVals, thetaVals, phiVals):
        ########## check sanity of the inputs
        ToroCoords.__check_valid_coords(etaVals, thetaVals, phiVals)

        preFac = self.aParam / (np.cosh(etaVals) - np.cos(thetaVals))

        xVals = preFac * np.cos(phiVals) * np.sinh(etaVals)
        yVals = preFac * np.sin(phiVals) * np.sinh(etaVals)
        zVals = preFac * np.sin(thetaVals)

        return (xVals, yVals, zVals)

    # convert cartesian coordinates to toroidal coordinates
    # accept xVals, yVals, zVals of any shape as long as it is the same shape
    def Cart_to_Toro(self, xVals, yVals, zVals):
        ########## check sanity of the inputs
        ToroCoords.__check_valid_coords(xVals, yVals, zVals)

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

    # get the radius squared given the toroidal tuple (eta, theta, phi)
    # do check is a switch to turn on/off the checking of input sanity
    # helps if there are nested calculations, where input only needs to be checked once
    def Rad2(self, etaVals, thetaVals, phiVals, doCheck=True):
        if doCheck: ToroCoords.__check_valid_coords(etaVals, thetaVals, phiVals)

        return ( self.aParam**2 * (np.cosh(etaVals)+np.cos(thetaVals))/(np.cosh(etaVals)-np.cos(thetaVals)) )

    # derivative of radius squared with respect to eta
    def Rad2_eta(self, etaVals, thetaVals, phiVals, doCheck=True):
        if doCheck: ToroCoords.__check_valid_coords(etaVals, thetaVals, phiVals)

        return (
                    -2*self.aParam**2 *
                    ( np.sinh(etaVals)*np.cos(thetaVals) ) /
                                        ( np.cosh(etaVals)-np.cos(thetaVals) )**2
               )

    # derivative of radius squared with respect to theta
    def Rad2_theta(self, etaVals, thetaVals, phiVals, doCheck=True):
        if doCheck: ToroCoords.__check_valid_coords(etaVals, thetaVals, phiVals)

        return (
                -2 * self.aParam ** 2 *
                (np.cosh(etaVals) * np.sin(thetaVals)) /
                (np.cosh(etaVals) - np.cos(thetaVals)) ** 2
        )

    # second order derivative of radius squared with respect to eta, eta
    def Rad2_eta_eta(self, etaVals, thetaVals, phiVals, doCheck=True):
        if doCheck: ToroCoords.__check_valid_coords(etaVals, thetaVals, phiVals)

        return (
                self.aParam ** 2 *
                ( np.cosh(2*etaVals) + 2*np.cosh(etaVals)*np.cos(thetaVals) - 3 )*np.cos(thetaVals) /
                (np.cosh(etaVals) - np.cos(thetaVals)) ** 3
        )

    # second order derivative of radius squared with respect to theta, theta
    def Rad2_theta_theta(self, etaVals, thetaVals, phiVals, doCheck=True):
        if doCheck: ToroCoords.__check_valid_coords(etaVals, thetaVals, phiVals)

        return (
                2*(self.aParam ** 2) *
                np.cosh(etaVals)*( np.sin(thetaVals)**2 + 1 - np.cosh(etaVals)*np.cos(thetaVals) ) /
                (np.cosh(etaVals) - np.cos(thetaVals)) ** 3
        )

    # second order derivative of radius squared with respect to eta, theta
    def Rad2_eta_theta(self, etaVals, thetaVals, phiVals, doCheck=True):
        if doCheck: ToroCoords.__check_valid_coords(etaVals, thetaVals, phiVals)

        return (
                2*(self.aParam ** 2) *
                np.sinh(etaVals)*( np.cosh(etaVals) + np.cos(thetaVals) )*np.sin(thetaVals) /
                (np.cosh(etaVals) - np.cos(thetaVals)) ** 3
        )