# Written by Vassili Savinov on 17/09/2018
# a class wrapper to load and drop the DLL with the DTORTH functions

import ctypes as ct
import numpy as np

class DTORH:
    class NewNTooSmall(Exception): pass# new N is too small in dtorh1
    class Internal(Exception): pass  # new N is too small in dtorh1

    dllTorHarm = None

    # keep it simple for now, the DLL must be in the same folder
    #gets called by __enter__ amongst other things
    def __init__(self):
        self.dllTorHarm = ct.CDLL('wrapDTORH64.dll')

    def __enter__(self):# with statement
        return self

    def __exit__(self, exc_type, exc_value, traceback):# with statement
        dllHandle=self.dllTorHarm._handle
        del self.dllTorHarm#release object
        ct.windll.kernel32.FreeLibrary(dllHandle)#release lib

    # handle the errors raised by the DTORH functions
    def __HandleErrorCode(self, errCode):
        if errCode==0: return# no error

        if errCode   == 1 : raise DTORH.Internal('IPRE MUST BE 1 OR 2')
        elif errCode == 2 : raise DTORH.Internal('YOU MUST CHOOSE MODE=2')
        elif errCode == 3 : raise DTORH.Internal('M IS TOO LARGE FOR MODE=0. BETTER TRY MODE=1')
        elif errCode == 4 : raise DTORH.Internal('M IS TOO LARGE FOR MODE=0. BETTER TRY MODE=1')
        elif errCode == 5 : raise DTORH.Internal('IMPROPER ARGUMENT. Z MUST BE GREATER THAN 1')
        elif errCode == 6 : raise DTORH.Internal('M IS TOO LARGE FOR MODE=1,2')
        elif errCode == 10: raise DTORH.Internal('Invalid mode (need 0,1,2)!')
        elif errCode == 11: raise DTORH.Internal('Continued fraction convergence failed!')
        else:               raise DTORH.Internal('Uknown error')

    # mode 0 means evalute as equired
    # mode 1 normalize results by gamma(m+1/2)
    # mode 2, similar to mode 1, but without restriction on large z
    def FixedM(self, zVec, mVal, nVal, mode=0):
        plVec=np.zeros(len(zVec), dtype=np.double)
        qlVec = np.zeros(len(zVec), dtype=np.double)

        # allow to go one higher and to start at zero
        # +100 just for extra memory, it sometimes failes without it
        qlMem = np.zeros(nVal + 2+100, dtype=np.double)  # scaled by 1/gamma(m+1/2)
        plMem = np.zeros(nVal + 2+100, dtype=np.double)

        # diagonostics
        c_errCode=ct.c_int(0)
        c_newN=ct.c_int(0)

        iZ=0

        for zVal in zVec:
            self.dllTorHarm.wrapDTORH1((ct.c_double)(zVal),
                                  (ct.c_int)(mVal),
                                  (ct.c_int)(nVal),
                                  plMem.ctypes.data_as(ct.POINTER(ct.c_double)),
                                  qlMem.ctypes.data_as(ct.POINTER(ct.c_double)),
                                  ct.pointer(c_newN),
                                  ct.pointer(c_errCode),
                                  ct.c_int(mode))

            if c_newN.value<nVal:
                raise DTORH.NewNTooSmall('FixedM: New N is smaller than the target!')

            self.__HandleErrorCode(c_errCode.value)

            plVec[iZ] = plMem[nVal]
            qlVec[iZ] = qlMem[nVal]

            iZ=iZ+1

        return (plVec, qlVec)

#    !!!!!!!!!!!!! need code to extract NM matrices for P, and Q