# Written by Vassili Savinov on 05/02/2019
# test ToroHarmVecRep

import numpy as np
from ToroHarmVecRep import ToroHarmVecRep

# we need a vector field of some sort
# lets well-behaved but non-central field \vec{F}=\vec{r}*x*exp(-(x^2+(y-2)^2+(z-1)^2))
# so \vec{\nabla}.\vec{F}=2x*(2+z-r^2)*exp(-(x^2+(y-2)^2+(z-1)^2))

def divF_func(x, y, z): return 2*x*(2+z-(x**2+y**2+z**2))*np.exp(-(x**2+(y-2)**2+(z-1)**2))


xVec = np.linspace(-2, 2, 13)
yVec = np.linspace(-2, 2, 17)
zVec = np.linspace(-2, 2, 19)

xTens, yTens, zTens = np.meshgrid(xVec, yVec, zVec, indexing='ij')
divFTens=divF_func(xTens, yTens, zTens)
rDotFTens = (xTens**2+yTens**2+zTens**2)*xTens*np.exp(-(xTens**2+(yTens-2)**2+(zTens-1)**2))
LDotFTens = 0*xTens

toroHarmVec = ToroHarmVecRep(xTens, yTens, zTens, divFTens, rDotFTens, LDotFTens)

pass
