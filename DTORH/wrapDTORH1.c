/*
Written by Vassili Savinov on 10/09/2018
source file for the c-wrappers around the DORTH functions by

Computer Physics Communications 124 (2000) 104–122
Evaluation of toroidal harmonics, J. Segura, A. Gil 
*/

#include <stdio.h>
#include "wrapDTORH.h"

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
// fortran functions are atomatically underlined
/*C ERR is set to inform of any errors during execution
C ERR=0 means OK
C ERR=1 WRITE(6,*)'IPRE MUST BE 1 OR 2'
C ERR=2 WRITE(*,*)'YOU MUST CHOOSE MODE=2'
C ERR=3 WRITE(*,*)'M IS TOO LARGE FOR MODE=0' WRITE(*,*)'BETTER TRY MODE=1'
C ERR=4 WRITE(*,*)'M IS TOO LARGE FOR MODE=0' WRITE(*,*)'BETTER TRY MODE=1'
C ERR=5 WRITE(6,*)'IMPROPER ARGUMENT. Z MUST BE GREATER THAN 1'
C ERR=6 WRITE(*,*)'M IS TOO LARGE FOR MODE=1,2'
C ERR=10 invalid mode, has to be 0,1,2*/
extern void dtorh1m0_(double*, int*, int*, double*, double*, int*, int*);// hard-coded modes mode 0
extern void dtorh1m1_(double*, int*, int*, double*, double*, int*, int*);
extern void dtorh1m2_(double*, int*, int*, double*, double*, int*, int*);


//!!!!!!!!!!!!!!!!!!!! account for mode selection (different fortran functions + one wrapper)

// could link directly to fortran, but this is easier to debug
void __cdecl wrapDTORH1(double zVal, int mVal, int nMax, double* PLVec, double* QLVec, int* newN, int* err, int mode)
{
	switch(mode)
	{
		// pass into FORTRAN
		case 0: dtorh1m0_(&zVal, &mVal, &nMax, PLVec, QLVec, newN, err);break;
		case 1: dtorh1m1_(&zVal, &mVal, &nMax, PLVec, QLVec, newN, err);break;
		case 2: dtorh1m2_(&zVal, &mVal, &nMax, PLVec, QLVec, newN, err);break;
		default: *err=10; break;
	};	
	
	return;
} 


