/*
Written by Vassili Savinov on 10/09/2018
header file for the c-wrappers around the DORTH functions by

Computer Physics Communications 124 (2000) 104–122
Evaluation of toroidal harmonics, J. Segura, A. Gil 
*/

#ifndef WRAPDTORH_H___
#define WRAPDTORH_H___

//DTORH1(Z,M,NMAX,PL,QL,NEWN,ERR)   
/*
INTEGER M,NMAX,NMAXP,MODE,IPRE,MP,NP,N,ICAL,NEWN,I
       DOUBLE PRECISION Z,PI,EPS,TINY,OVER,TINYSQ,QZ,PISQ,DPPI,
     *   FL,CC,AR,GAMMA,FC,QDC1,QARGU,ARGU1,DFACQS,FCP,DD,
     *   GAMMAH,ELLIP1,ELLIP2,D1,QM0,DFAC3,DFAC4,PL0,FACTCO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     
C THE DIMENSION OF QLMM (INTERNAL ARRAY) MUST BE GREATER THAN M  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
       DOUBLE PRECISION  QLMM(0:1001),PL(0:NMAX+1),QL(0:NMAX+1),PR(2)
       PARAMETER(PI=3.14159265358979323D0,EPS=1.D-14,TINY=1.D-290,
     *           MODE=1,IPRE=1)
       OVER=1.D0/TINY
       TINYSQ=DSQRT(TINY)
*/
/*C ERR is set to inform of any errors during execution
C ERR=0 means OK
C ERR=1 WRITE(6,*)'IPRE MUST BE 1 OR 2'
C ERR=2 WRITE(*,*)'YOU MUST CHOOSE MODE=2'
C ERR=3 WRITE(*,*)'M IS TOO LARGE FOR MODE=0' WRITE(*,*)'BETTER TRY MODE=1'
C ERR=4 WRITE(*,*)'M IS TOO LARGE FOR MODE=0' WRITE(*,*)'BETTER TRY MODE=1'
C ERR=5 WRITE(6,*)'IMPROPER ARGUMENT. Z MUST BE GREATER THAN 1'
C ERR=6 WRITE(*,*)'M IS TOO LARGE FOR MODE=1,2'
C ERR=10 invalid mode, has to be 0,1,2
*/
void __cdecl wrapDTORH1(double zVal, int mVal, int nMax, double* PLVec, double* QLVec, int* newN, int* err, int mode);

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
//  SUBROUTINE DTORH2(Z,MDIM,NDIM,MMAX,NMAX,PL,QL,NEWM,NEWN)
/*C ERR is set to inform of any errors during execution
C ERR is set to inform of any errors during execution
C ERR=0 means OK
C ERR=1 WRITE(6,*)'IPRE MUST BE 1 OR 2'
C ERR=2 WRITE(*,*)'YOU MUST CHOOSE MODE=2'
C ERR=5 WRITE(6,*)'IMPROPER ARGUMENT. Z MUST BE GREATER THAN 1'
C ERR=7 WRITE(6,*)'MMAX IS TOO LARGE FOR MODE=0' WRITE(6,*)'BETTER TRY MODE=1'
C ERR=10 invalid mode, has to be 0,1,2
*/
void __cdecl wrapDTORH2(double zVal, int mDim, int nDim, int mMax, 
								int nMax, double* PLVec, double* QLVec, int* newM, int* newNM, int* err, int mode);
								
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
C ERR=10 invalid mode, has to be 0,1,2
*/
void __cdecl wrapDTORH3(double zVal, int mDim, int nDim, int mMax, int nMax, 
									double* PLVec, double* QLVec, int* newM, int* newNM, int* err, int mode);

#endif