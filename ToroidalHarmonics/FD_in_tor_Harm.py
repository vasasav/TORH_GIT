# Written by Vassili Savinov on 11/03/19
# decompose the field of the flying doughtnut in toroidal harmonics
# The relevant mathematical operations have been obtained by the Mathematica field

import unittest
from ToroHarmVecRep import ToroHarmVecRep
import numpy as np
import pylab as pl
import pandas as pd
import matplotlib.colors as mpc

# load the datat computed by mathematica
col_names = ['X', 'Y', 'Z', 'q1', 'q2',
             'r_dot_E', 'L_dot_E_im', 'r_dot_B', 'L_dot_B_im']
fdData = pd.read_csv('FD_numData.csv', names=col_names)

# prep data
X = np.array(fdData['X'])
Y = np.array(fdData['Y'])
Z = np.array(fdData['Z'])
#
r_dot_E = np.array(fdData['r_dot_E'], dtype=np.complex)
L_dot_E = 1j*np.array(fdData['L_dot_E_im'], dtype=np.complex)
#
r_dot_B = np.array(fdData['r_dot_B'], dtype=np.complex)
L_dot_B = 1j*np.array(fdData['L_dot_B_im'], dtype=np.complex)

##### do the decomposition
fdToroHarmVecRep = ToroHarmVecRep(X, Y, Z,
                                    raw_divFTens=0*X,
                                    raw_rDotFTens=r_dot_E,
                                    raw_LDotFTens=L_dot_E,
                                    nCount=15, mCount=19)

########### display
aCoeff_max_val=10*np.log10(np.max( np.abs(fdToroHarmVecRep.aCoeff_Tens) ))
aCoeff_col_abs_norm = mpc.Normalize(vmin=aCoeff_max_val-60, vmax=aCoeff_max_val)
#
bCoeff_max_val=10*np.log10(np.max( np.abs(fdToroHarmVecRep.bCoeff_Tens) ))
bCoeff_col_abs_norm = mpc.Normalize(vmin=bCoeff_max_val-60, vmax=bCoeff_max_val)
#
cCoeff_max_val=10*np.log10(np.max( np.abs(fdToroHarmVecRep.cCoeff_Tens) ))
cCoeff_col_abs_norm = mpc.Normalize(vmin=cCoeff_max_val-60, vmax=cCoeff_max_val)
#
col_pha_norm = mpc.Normalize(vmin=-180, vmax=180)

pl.style.use('dark_background')

tickStep = 5#dB

pl.figure(1)
###
pl.subplot(121)
pl.imshow(10*np.log10( np.abs(np.squeeze( fdToroHarmVecRep.bCoeff_Tens[0,:,:] )) ),cmap=pl.cm.hot)
pl.colorbar().set_label(' (dB)')
pl.xlabel('m - order')
pl.ylabel('n - order')
pl.yticks(np.arange(0, fdToroHarmVecRep.bCoeff_Tens.shape[2], tickStep))
pl.xticks(np.arange(0, fdToroHarmVecRep.bCoeff_Tens.shape[1], tickStep))

pl.subplot(122)
pl.imshow(10*np.log10( np.abs(np.squeeze( fdToroHarmVecRep.bCoeff_Tens[1,:,:] )) ),  cmap=pl.cm.hot)
pl.colorbar().set_label(' (dB)')
pl.xlabel('m - order')
pl.ylabel('n - order')
pl.yticks(np.arange(0, fdToroHarmVecRep.bCoeff_Tens.shape[2], tickStep))
pl.xticks(np.arange(0, fdToroHarmVecRep.bCoeff_Tens.shape[1], tickStep))

pl.figure(2)

pl.subplot(121)
pl.imshow(np.angle(np.squeeze(fdToroHarmVecRep.bCoeff_Tens[0,:,:]), deg=True), norm=col_pha_norm, cmap=pl.cm.jet)
pl.colorbar().set_label('(deg)')
pl.xlabel('m - order')
pl.ylabel('n - order')
pl.yticks(np.arange(0, fdToroHarmVecRep.aCoeff_Tens.shape[2], tickStep))
pl.xticks(np.arange(0, fdToroHarmVecRep.aCoeff_Tens.shape[1], tickStep))

pl.subplot(122)
pl.imshow(np.angle(np.squeeze(fdToroHarmVecRep.bCoeff_Tens[1,:,:]), deg=True), norm=col_pha_norm, cmap=pl.cm.jet)
pl.colorbar().set_label('(deg)')
pl.xlabel('m - order')
pl.ylabel('n - order')
pl.yticks(np.arange(0, fdToroHarmVecRep.aCoeff_Tens.shape[2], tickStep))
pl.xticks(np.arange(0, fdToroHarmVecRep.aCoeff_Tens.shape[1], tickStep))


pl.show()
