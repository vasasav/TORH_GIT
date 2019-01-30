# Written by Vassili Savinov on 30/01/2019
# It seems that Mathematica does support the right kind of associated Legendre functions
# P^m_{n-1/2 }(z) = LegendreP[n-1/2, m, 3, z], 3 here is the kind of function,
# note that Segura and Gil use notation P^m_n-1/2 (same as Mathematica)
# we want the one defined for z>=1

# I just wanted to have this script as a simple dump to load the fortran library

from DTORH import DTORH
import numpy as np
import pylab as pl

# prep z
zVec=np.array([200], dtype=np.double)

# get the data
mode=0
m=4
n=2

with DTORH() as dtorh: (plVec, qlVec) = dtorh.FixedM(zVec, m, n, mode)

print('m=%d n=%d\t plVec[0] = %.5e' % (m, n, plVec[0]))
print('m=%d n=%d\t qlVec[0] = %.5e' % (m, n, qlVec[0]))

# Checked few values, there seems to be agreement. Good.