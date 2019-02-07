# Written by Vassili Savinov on 05/02/2019
# test ToroHarmVecRep

import pylab as pl
import numpy as np
from ToroHarmVecRep import ToroHarmVecRep
import matplotlib.colors as mpc

# we need a vector field of some sort
# use a well-behaved vector field F={y, x, z}*Exp[-(x^2 + (y - 1)^2 + z^2)]
# div(F) = exp(-x^2 - (-1 + y)^2 - z^2)* (1 + x (2 - 4 y) - 2 z^2)
# r.F = exp(-x^2 - (-1 + y)^2 - z^2) * (2 x y + z^2)
# L.F = -i2*exp(-x^2 - (-1 + y)^2 - z^2) * (x - y) z

xVec = np.linspace(-2.5, 2.5, 35)
yVec = np.linspace(-2.5, 2.5, 35)
zVec = np.linspace(-2.5, 2.5, 35)

xTens, yTens, zTens = np.meshgrid(xVec, yVec, zVec, indexing='ij')
#
divFTens = np.exp(- (xTens**2 + (yTens-1)**2 + zTens**2) ) * (1 + xTens*(2 - 4*yTens) - 2*zTens**2)
#
rDotFTens = np.exp(- (xTens**2 + (yTens-1)**2 + zTens**2) ) * (2*xTens*yTens + zTens**2)
#
LDotFTens = -1j * 2 * np.exp(- (xTens**2 + (yTens-1)**2 + zTens**2) )  * (xTens - yTens) * zTens
#

toroHarmVec = ToroHarmVecRep(xTens, yTens, zTens, divFTens, rDotFTens, LDotFTens)

############# now plot

aCoeff_max_val=10*np.log10(np.max( np.abs(toroHarmVec.aCoeff_Tens) ))
aCoeff_col_abs_norm = mpc.Normalize(vmin=aCoeff_max_val-60, vmax=aCoeff_max_val)
#
bCoeff_max_val=10*np.log10(np.max( np.abs(toroHarmVec.bCoeff_Tens) ))
bCoeff_col_abs_norm = mpc.Normalize(vmin=bCoeff_max_val-60, vmax=bCoeff_max_val)
#
cCoeff_max_val=10*np.log10(np.max( np.abs(toroHarmVec.cCoeff_Tens) ))
cCoeff_col_abs_norm = mpc.Normalize(vmin=cCoeff_max_val-60, vmax=cCoeff_max_val)
#
col_pha_norm = mpc.Normalize(vmin=-180, vmax=180)

pl.style.use('dark_background')

tickStep = 5#dB

pl.figure(1)
###
pl.subplot(321)
pl.imshow(10*np.log10( np.abs(np.squeeze( toroHarmVec.aCoeff_Tens[0,:,:] )) ), norm=aCoeff_col_abs_norm, cmap=pl.cm.hot)
pl.colorbar().set_label(' (dB)')
pl.xlabel('m - order')
pl.ylabel('n - order')
pl.yticks(np.arange(0, toroHarmVec.aCoeff_Tens.shape[2], tickStep))
pl.xticks(np.arange(0, toroHarmVec.aCoeff_Tens.shape[1], tickStep))
#pl.title('$\\left|a^{(1)}\\right|$')
###
pl.subplot(322)
pl.imshow(10*np.log10( np.abs(np.squeeze( toroHarmVec.aCoeff_Tens[1,:,:] )) ), norm=aCoeff_col_abs_norm, cmap=pl.cm.hot)
pl.colorbar().set_label(' (dB)')
pl.xlabel('m - order')
pl.ylabel('n - order')
pl.yticks(np.arange(0, toroHarmVec.aCoeff_Tens.shape[2], tickStep))
pl.xticks(np.arange(0, toroHarmVec.aCoeff_Tens.shape[1], tickStep))
#pl.title('$\\left|a^{(2)}\\right|$')
###
pl.subplot(323)
pl.imshow(10*np.log10( np.abs(np.squeeze( toroHarmVec.bCoeff_Tens[0,:,:] )) ), norm=bCoeff_col_abs_norm, cmap=pl.cm.hot)
pl.colorbar().set_label(' (dB)')
pl.xlabel('m - order')
pl.ylabel('n - order')
pl.yticks(np.arange(0, toroHarmVec.aCoeff_Tens.shape[2], tickStep))
pl.xticks(np.arange(0, toroHarmVec.aCoeff_Tens.shape[1], tickStep))
#pl.title('$\\left|b^{(1)}\\right|$')
###
pl.subplot(324)
pl.imshow(10*np.log10( np.abs(np.squeeze( toroHarmVec.bCoeff_Tens[1,:,:] )) ), norm=bCoeff_col_abs_norm, cmap=pl.cm.hot)
pl.colorbar().set_label(' (dB)')
pl.xlabel('m - order')
pl.ylabel('n - order')
pl.yticks(np.arange(0, toroHarmVec.aCoeff_Tens.shape[2], tickStep))
pl.xticks(np.arange(0, toroHarmVec.aCoeff_Tens.shape[1], tickStep))
#pl.title('$\\left|b^{(2)}\\right|$')
###
pl.subplot(325)
pl.imshow(10*np.log10( np.abs(np.squeeze( toroHarmVec.cCoeff_Tens[0,:,:] )) ), norm=cCoeff_col_abs_norm, cmap=pl.cm.hot)
pl.colorbar().set_label(' (dB)')
pl.xlabel('m - order')
pl.ylabel('n - order')
pl.yticks(np.arange(0, toroHarmVec.aCoeff_Tens.shape[2], tickStep))
pl.xticks(np.arange(0, toroHarmVec.aCoeff_Tens.shape[1], tickStep))
#pl.title('$\\left|c^{(1)}\\right|$')
###
pl.subplot(326)
pl.imshow(10*np.log10( np.abs(np.squeeze( toroHarmVec.cCoeff_Tens[1,:,:] )) ), norm=cCoeff_col_abs_norm, cmap=pl.cm.hot)
pl.colorbar().set_label(' (dB)')
pl.xlabel('m - order')
pl.ylabel('n - order')
pl.yticks(np.arange(0, toroHarmVec.aCoeff_Tens.shape[2], tickStep))
pl.xticks(np.arange(0, toroHarmVec.aCoeff_Tens.shape[1], tickStep))
#pl.title('$\\left|c^{(2)}\\right|$')

########## phase

pl.figure(2)

pl.subplot(321)
pl.imshow(np.angle(np.squeeze(toroHarmVec.aCoeff_Tens[0,:,:]), deg=True), norm=col_pha_norm, cmap=pl.cm.jet)
pl.colorbar().set_label('(deg)')
pl.xlabel('m - order')
pl.ylabel('n - order')
pl.yticks(np.arange(0, toroHarmVec.aCoeff_Tens.shape[2], tickStep))
pl.xticks(np.arange(0, toroHarmVec.aCoeff_Tens.shape[1], tickStep))
#pl.title('$\\left|c^{(2)}\\right|$')

pl.subplot(322)
pl.imshow(np.angle(np.squeeze(toroHarmVec.aCoeff_Tens[1,:,:]), deg=True), norm=col_pha_norm, cmap=pl.cm.jet)
pl.colorbar().set_label('(deg)')
pl.xlabel('m - order')
pl.ylabel('n - order')
pl.yticks(np.arange(0, toroHarmVec.aCoeff_Tens.shape[2], tickStep))
pl.xticks(np.arange(0, toroHarmVec.aCoeff_Tens.shape[1], tickStep))
#pl.title('$\\left|c^{(2)}\\right|$')

pl.subplot(323)
pl.imshow(np.angle(np.squeeze(toroHarmVec.bCoeff_Tens[0,:,:]), deg=True), norm=col_pha_norm, cmap=pl.cm.jet)
pl.colorbar().set_label('(deg)')
pl.xlabel('m - order')
pl.ylabel('n - order')
pl.yticks(np.arange(0, toroHarmVec.aCoeff_Tens.shape[2], tickStep))
pl.xticks(np.arange(0, toroHarmVec.aCoeff_Tens.shape[1], tickStep))
#pl.title('$\\left|c^{(2)}\\right|$')

pl.subplot(324)
pl.imshow(np.angle(np.squeeze(toroHarmVec.bCoeff_Tens[1,:,:]), deg=True), norm=col_pha_norm, cmap=pl.cm.jet)
pl.colorbar().set_label('(deg)')
pl.xlabel('m - order')
pl.ylabel('n - order')
pl.yticks(np.arange(0, toroHarmVec.aCoeff_Tens.shape[2], tickStep))
pl.xticks(np.arange(0, toroHarmVec.aCoeff_Tens.shape[1], tickStep))
#pl.title('$\\left|c^{(2)}\\right|$')

pl.subplot(325)
pl.imshow(np.angle(np.squeeze(toroHarmVec.cCoeff_Tens[0,:,:]), deg=True), norm=col_pha_norm, cmap=pl.cm.jet)
pl.colorbar().set_label('(deg)')
pl.xlabel('m - order')
pl.ylabel('n - order')
pl.yticks(np.arange(0, toroHarmVec.aCoeff_Tens.shape[2], tickStep))
pl.xticks(np.arange(0, toroHarmVec.aCoeff_Tens.shape[1], tickStep))
#pl.title('$\\left|c^{(2)}\\right|$')

pl.subplot(326)
pl.imshow(np.angle(np.squeeze(toroHarmVec.cCoeff_Tens[1,:,:]), deg=True), norm=col_pha_norm, cmap=pl.cm.jet)
pl.colorbar().set_label('(deg)')
pl.xlabel('m - order')
pl.ylabel('n - order')
pl.yticks(np.arange(0, toroHarmVec.aCoeff_Tens.shape[2], tickStep))
pl.xticks(np.arange(0, toroHarmVec.aCoeff_Tens.shape[1], tickStep))
#pl.title('$\\left|c^{(2)}\\right|$')

pl.show()