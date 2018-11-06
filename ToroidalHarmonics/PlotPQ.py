from DTORH import DTORH
import numpy as np
import pylab as pl

# prep z
zVec=np.linspace(1.1, 10, 500, dtype=np.double)

# get the data
mode=0
m=0
n=0
with DTORH() as dtorh: (plVec, qlVec) = dtorh.FixedM(zVec, m, n, mode)

# plot it
pl.figure(1)
pl.plot(zVec, plVec, label=('P m=%d, n=%d' % (m,n)))
pl.plot(zVec, qlVec, label=('Q m=%d, n=%d' % (m,n)))

m=3
n=5
with DTORH() as dtorh: (plVec, qlVec) = dtorh.FixedM(zVec, m, n, mode)
pl.plot(zVec, plVec, label=('P m=%d, n=%d' % (m,n)))
pl.plot(zVec, qlVec, label=('Q m=%d, n=%d' % (m,n)))

pl.legend()
pl.xlabel('z')
pl.ylabel('P,Q value')
pl.show()
