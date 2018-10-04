/*
Written by Vassili Savinov on 10/09/2018
source file for the c-wrappers around the DORTH functions by

Computer Physics Communications 124 (2000) 104–122
Evaluation of toroidal harmonics, J. Segura, A. Gil 
*/

#include <stdio.h>
#include "wrapDTORH.h"

/*
          SUBROUTINE DTORH3(Z,MDIM,NDIM,MMAX,NMAX,PL,QL,NEWM,NEWN)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  INPUT :                                                         C         
C    Z        ARGUMENT OF THE FUNCTIONS                            C       
C    MDIM     M-DIMENSION OF THE ARRAYS. MDIM MUST BE GREATER      C
C             THAN MMAX                                            C                    
C    NDIM     N-DIMENSION OF THE ARRAYS. NDIM MUST BE GREATER      C
C             THAN NMAX                                            C                                
C    MMAX     MAXIMUM ORDER  OF THE FUNCTIONS :                    C                  
C              WE CALCULATE FUNCTIONS OF ALL ORDERS BELOW MMAX.    C                   
C    NMAX     MAXIMUM DEGREE  OF THE FUNCTIONS :                   C                  
C              WE CALCULATE FUNCTIONS OF ALL DEGREES BELOW         C       
C              MIN(NMAX,NEWN).                                     C                                                                                C
C  OUTPUT :                                                        C       
C   *IF MODE IS EQUAL TO 0:                                        C                                                                               
C    PL(M,N)                                                       C       
C             THESE VALUES ARE KEPT IN AN ARRAY                    C       
C    QL(M,N)                                                       C       
C             THESE VALUES ARE KEPT IN AN ARRAY                    C       
C                                                                  C       
C    NEWM     MAXIMUM  ORDER OF FUNCTIONS CALCULATED WHEN          C       
C             QL (MMAX,0)   IS LARGER THAN 1/TINY                  C       
C              (OVERFLOW LIMIT = 1/TINY, TINY IS DEFINED BELOW).   C       
C    NEWN     MAXIMUM  DEGREE OF FUNCTIONS CALCULATED WHEN         C                 
C             PL (M,NMAX)   IS LARGER THAN 1/TINY  FOR SOME        C       
C             M=0,...,NEWM                                         C       
C              (OVERFLOW LIMIT = 1/TINY, TINY IS DEFINED BELOW).   C       
C    NOTE1: FOR A PRECISION OF 10**(-12), IF Z>5 AND (Z/M)>0.22    C
C           THE CODE USES A SERIES EXPANSION FOR PL(M,0). WHEN     C
C           Z<20 AND (Z/M)<0.22 A CONTINUED FRACTION IS APPLIED.   C                           
C    NOTE2: FOR A PRECISION OF 10**(-8), IF Z>5 AND (Z/M)>0.12     C
C           THE CODE USES A SERIES EXPANSION FOR PL(M,0). WHEN     C
C           Z<20 AND (Z/M)<0.12 A CONTINUED FRACTION IS APPLIED.   C                                                                                                    
C   *IF MODE IS EQUAL TO 1:                                        C       
C      THE SET OF FUNCTIONS EVALUATED IS:                          C       
C                 PL(M,N)/GAMMA(M+1/2),QL(M,N)/GAMMA(M+1/2),       C       
C      WHICH ARE RESPECTIVELY STORED IN THE ARRAYS PL(M,N),QL(M,N) C       
C      NEWM AND NEWN REFER TO THIS NEW SET OF FUNCTIONS            C       
C      NOTE1 AND NOTE2 ALSO APPLY IN THIS CASE                     C       
C   *IF MODE IS EQUAL TO 2:                                        C       
C      THE CODE PERFORMS AS FOR MODE 1, BUT THE RESTRICTION Z<20   C
C      FOR THE EVALUATION OF THE CONTINUED FRACTION IS NOT         C
C      CONSIDERED                                                  C
C      WARNING: USE ONLY IF HIGH M'S FOR Z>20 ARE REQUIRED. THE    C      
C      EVALUATION OF THE CF MAY FAIL TO CONVERGE FOR TOO HIGH Z'S  C      
C  PARAMETERS:                                                     C       
C   MODE: ENABLES THE ENLARGEMENT OF THE RANGE OF ORDERS AND       C
C         DEGREES THAT CAN BE EVALUATED.                           C                    
C   EPS:  CONTROLS THE ACCURACY OF THE CONTINUED FRACTIONS         C
C         AND SERIES.                                              C
C   IPRE: REQUIRED PRECISION IN THE EVALUATION OF TOROIDAL         C
C         HARMONICS.                                               C
C           *IF IPRE=1, PRECISION=10**(-12) (TAKING EPS<10**(-12)) C                                                                                C
C           *IF IPRE=2, PRECISION=10**(-8) (TAKING EPS<10**(-8))   C       
C   TINY: SMALL PARAMETER NEAR THE UNDERFLOW LIMIT.                C                   
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   
*/

// SUBROUTINE DTORH3(Z,MDIM,NDIM,MMAX,NMAX,PL,QL,NEWM,NEWN)
/*
C ERR is set to inform of any errors during execution
C ERR=0 means OK
C ERR=1 WRITE(6,*)'IPRE MUST BE 1 OR 2'
C ERR=2 WRITE(*,*)'YOU MUST CHOOSE MODE=2'
C ERR=5 WRITE(6,*)'IMPROPER ARGUMENT. Z MUST BE GREATER THAN 1'
C ERR=7 WRITE(6,*)'MMAX IS TOO LARGE FOR MODE=0' WRITE(6,*)'BETTER TRY MODE=1'
*/
extern void dtorh3m0_(double*, int*, int*, int*, int*, double**, double**, int*, int*, int*);
extern void dtorh3m1_(double*, int*, int*, int*, int*, double**, double**, int*, int*, int*);
extern void dtorh3m2_(double*, int*, int*, int*, int*, double**, double**, int*, int*, int*);


// could link directly to fortran, but this is easier to debug
void __cdecl wrapDTORH3(double zVal, int mDim, int nDim, int mMax, int nMax, 
									double* PLVec, double* QLVec, int* newM, int* newN, int* err, int mode)
{
	// pass into FORTRAN
	// I will convert from 1d to fortran 2d separatelly - too much headache
	;
	
	switch(mode)
	{
		// pass into FORTRAN
		case 0: dtorh3m0_(&zVal, &mDim, &nDim, &mMax, &nMax, (double**)PLVec, (double**)QLVec, newM, newN, err);break;
		case 1: dtorh3m1_(&zVal, &mDim, &nDim, &mMax, &nMax, (double**)PLVec, (double**)QLVec, newM, newN, err);break;
		case 2: dtorh3m2_(&zVal, &mDim, &nDim, &mMax, &nMax, (double**)PLVec, (double**)QLVec, newM, newN, err);break;
		default: *err=10; break;
	};	
	
	return;
} 

