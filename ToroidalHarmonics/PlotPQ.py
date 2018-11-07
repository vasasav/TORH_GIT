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

m=0
n=1
with DTORH() as dtorh: (plVec, qlVec) = dtorh.FixedM(zVec, m, n, mode)
pl.plot(zVec, plVec, label=('P m=%d, n=%d' % (m,n)))
pl.plot(zVec, qlVec, label=('Q m=%d, n=%d' % (m,n)))

pl.legend()
pl.xlabel('z')
pl.ylabel('P,Q value')


### test matrix

nCount=3
mCount=4
zRange=np.array((1.1, 1.5, 2))


m=1
n=2
iZ=2
with DTORH() as dtorh: (pCube, qCube) = dtorh.GetCubeNMZ(nCount, mCount, zRange, mode)
with DTORH() as dtorh: (plVec, qlVec) = dtorh.FixedM(zRange, m, n, mode)

print('test  \t p[%d,%d]: %.6e' % (n, m, pCube[n,m,iZ]))
print('sanity\t p[%d,%d]: %.6e' % (n, m, plVec[iZ]))
print('diff=%.6f' % ((pCube[n,m,iZ]-plVec[iZ])/plVec[iZ]))

print('test  \t q[%d,%d]: %.6e' % (n, m, qCube[n,m,iZ]))
print('sanity\t q[%d,%d]: %.6e' % (n, m, qlVec[iZ]))
print('diff=%.6f' % ((qCube[n,m,iZ]-qlVec[iZ])/qlVec[iZ]))

#pl.show()