# Written by Vassili Savinov on 05/02/19
# test ToroCoords

import numpy as np
import numpy.random as npr
from ToroCoords import ToroCoords

########### Testing coordinate changes
# prepare data
xVals = (2*npr.rand(100, 35, 78)-1)*100#np.array([0, 1.5], dtype=np.double)
yVals = (2*npr.rand(100, 35, 78)-1)*100#np.array([0, 1.5], dtype=np.double)
zVals = (2*npr.rand(100, 35, 78)-1)*100#np.array([0, 1.5], dtype=np.double)
# prep coord system
tC = ToroCoords()

# forward
etaVals, thetaVals, phiVals = tC.Cart_to_Toro(xVals, yVals, zVals)

# back
xVals1, yVals1, zVals1 = tC.Toro_to_Cart(etaVals, thetaVals, phiVals)

# compare
diffX = np.max( np.abs( ( xVals-xVals1)/(xVals+1e-10) ) )
diffY = np.max( np.abs( ( yVals-yVals1)/(yVals+1e-10) ) )
diffZ = np.max( np.abs( ( zVals-zVals1)/(zVals+1e-10) ) )

print( 'Max differences: %.5e, %.5e, %.5e' % (diffX, diffY, diffZ) )