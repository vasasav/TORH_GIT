# Written by Vassili Savinov on 05/02/19
# test ToroCoords

import numpy as np
import numpy.random as npr
from ToroCoords import ToroCoords

########### Testing coordinate changes
# prepare data
xVals = np.array([0, 1.7, -8], dtype=np.double)#(2*npr.rand(100, 35, 78)-1)*100#np.array([0, 1.5], dtype=np.double)
yVals = np.array([0, -3.8, 0.23], dtype=np.double)#(2*npr.rand(100, 35, 78)-1)*100#np.array([0, 1.5], dtype=np.double)
zVals = np.array([0, 0.2, +7.0], dtype=np.double)#(2*npr.rand(100, 35, 78)-1)*100#np.array([0, 1.5], dtype=np.double)
# prep coord system
tC = ToroCoords(np.pi)

# forward
torVals = tC.Cart_to_Toro(xVals, yVals, zVals)

# back
# NB! operator * unpacks tuples, nice way to pass parameters around
xVals1, yVals1, zVals1 = tC.Toro_to_Cart( *torVals )

# compare
diffX = np.max( np.abs( ( xVals-xVals1)/(xVals+1e-10) ) )
diffY = np.max( np.abs( ( yVals-yVals1)/(yVals+1e-10) ) )
diffZ = np.max( np.abs( ( zVals-zVals1)/(zVals+1e-10) ) )

################ aux methods
r2 = tC.Rad2(*torVals)
r2_eta = tC.Rad2_eta(*torVals)
r2_theta = tC.Rad2_theta(*torVals)
r2_eta_eta = tC.Rad2_eta_eta(*torVals)
r2_theta_theta = tC.Rad2_theta_theta(*torVals)
r2_eta_theta = tC.Rad2_eta_theta(*torVals)

print(*torVals)
print( 'Max differences: %.5e, %.5e, %.5e' % (diffX, diffY, diffZ) )
print( 'to check' )
print( r2_eta_theta )
