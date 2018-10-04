/*
Written by Vassili Savinov on 10/09/2018
source file for the c-wrappers around the DORTH functions by

Computer Physics Communications 124 (2000) 104–122
Evaluation of toroidal harmonics, J. Segura, A. Gil 
*/

#include <stdio.h>
#include "wrapDTORH.h"

/*
 SUBROUTINE DTORH2(Z,MDIM,NDIM,MMAX,NMAX,PL,QL,NEWM,NEWN)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  INPUT :                                                            C      
C    Z       ARGUMENT OF THE FUNCTIONS                                C                                                                         
C    MDIM    M-DIMENSION OF THE ARRAYS: MDIM MUST BE GREATER THAN     C
C            MMAX                                                     C                       
C    NDIM    N-DIMENSION OF THE ARRAYS: NDIM MUST BE GREATER THAN     C
C            NMAX                                                     C
C    MMAX    MAXIMUM ORDER OF THE FUNCTIONS :                         C   
C            WE CALCULATE FUNCTIONS OF ALL ORDERS BELOW               C   
C            MIN(NEWM,MMAX). NEWM IS DEFINED BELOW.                   C   
C    NMAX    MAXIMUM DEGREE OF THE FUNCTIONS :                        C   
C            WE GET  FUNCTIONS OF ALL THE DEGREES BELOW               C   
C            MIN(NEWN(M),NMAX). NEWN(M) IS DEFINED BELOW .            C                                                                           
C  OUTPUT :                                                           C    
C   *IF MODE IS EQUAL TO 0:                                           C                                                                            
C    PL(M,N)                                                          C    
C            THESE VALUES ARE KEPT IN AN ARRAY                        C   
C    QL(M,N)                                                          C    
C            THESE VALUES ARE KEPT IN AN ARRAY                        C   
C                                                                     C    
C    NEWM    MAXIMUM  ORDER  OF FUNCTIONS CALCULATED WHEN             C   
C            QL (MMAX+1,0)   IS LARGER THAN 1/TINY                    C   
C            (OVERFLOW LIMIT = 1/TINY, TINY IS DEFINED BELOW)         C   
C    NEWN(M) MAXIMUM  DEGREE  OF FUNCTIONS CALCULATED FOR A           C   
C            GIVEN ORDER M WHEN PL (M,NMAX+1) IS LARGER THAN          C
C            1/TINY (OVERFLOW LIMIT = 1/TINY, TINY IS DEFINED BELOW)  C          
C    NOTE1:  FOR A PRECISION OF 10**(-12), IF Z>5 AND (Z/M)>0.22 THE  C
C            CODE USES A SERIES EXPANSION FOR PL(M,0).                C
C            WHEN Z<20 AND (Z/M)<0.22 A CONTINUED FRACTION            C
C            IS APPLIED.                                              C
C    NOTE2:  FOR A PRECISION OF 10**(-8), IF Z>5 AND (Z/M)>0.12       C
C            THE CODE USES A SERIES EXPANSION FOR PL(M,0).WHEN Z<20   C
C            AND (Z/M)<0.12 A CONTINUED FRACTION IS APPLIED.          C    
*/

// SUBROUTINE DTORH2(Z,MDIM,NDIM,MMAX,NMAX,PL,QL,NEWM,NEWN)
/*
C ERR is set to inform of any errors during execution
C ERR is set to inform of any errors during execution
C ERR=0 means OK
C ERR=1 WRITE(6,*)'IPRE MUST BE 1 OR 2'
C ERR=2 WRITE(*,*)'YOU MUST CHOOSE MODE=2'
C ERR=5 WRITE(6,*)'IMPROPER ARGUMENT. Z MUST BE GREATER THAN 1'
C ERR=7 WRITE(6,*)'MMAX IS TOO LARGE FOR MODE=0' WRITE(6,*)'BETTER TRY MODE=1'
*/
extern void dtorh2m0_(double*, int*, int*, int*, int*, double**, double**, int*, int*, int*);
extern void dtorh2m1_(double*, int*, int*, int*, int*, double**, double**, int*, int*, int*);
extern void dtorh2m2_(double*, int*, int*, int*, int*, double**, double**, int*, int*, int*);


// could link directly to fortran, but this is easier to debug
void __cdecl wrapDTORH2(double zVal, int mDim, int nDim, int mMax, int nMax, 
									double* PLVec, double* QLVec, int* newM, int* newNM, int* err, int mode)
{
	// pass into FORTRAN
	// I will convert from 1d to fortran 2d separatelly - too much headache
	
	
	switch(mode)
	{
		// pass into FORTRAN
		case 0: dtorh2m0_(&zVal, &mDim, &nDim, &mMax, &nMax, (double**)PLVec, (double**)QLVec, newM, newNM, err);break;
		case 1: dtorh2m1_(&zVal, &mDim, &nDim, &mMax, &nMax, (double**)PLVec, (double**)QLVec, newM, newNM, err);break;
		case 2: dtorh2m2_(&zVal, &mDim, &nDim, &mMax, &nMax, (double**)PLVec, (double**)QLVec, newM, newNM, err);break;
		default: *err=10; break;
	};	
	
	return;
} 

