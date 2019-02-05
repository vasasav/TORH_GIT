# Written by Vassili Savinov on 05/02/2019
# getting familiar with tensorsolve

import numpy as np
import numpy.linalg as npla
import numpy.random as npr

# I will be primarily interested in rank three tensors, and rank 6 mixing tensors

inTens = 2*npr.rand(2, 20, 35)-1

# mixing tensor
mixTens = 2*npr.rand(*(*inTens.shape, *inTens.shape))-1

# constract
outTens = np.einsum('ijkstu,stu', mixTens, inTens)

print(outTens)

# now invert
solTens = npla.tensorsolve(mixTens, outTens)

# compare
diff=np.max( np.abs((inTens-solTens)/(inTens+1e-12)) )

print('max difference is %.3e %%' % (100*diff))