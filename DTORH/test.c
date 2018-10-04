#include <stdio.h>
#include "wrapDTORH.h"

#define NMAX 300

// assuming that in fortran one can use 0..dim as indices
#define FORT2D_TOC1D(M,N,MDIM) (N*(MDIM+1)+M)

void HandleErrorDTORH1(int err)
{
		/*C ERR is set to inform of any errors during execution
C ERR=0 means OK
C ERR=1 WRITE(6,*)'IPRE MUST BE 1 OR 2'
C ERR=2 WRITE(*,*)'YOU MUST CHOOSE MODE=2'
C ERR=3 WRITE(*,*)'M IS TOO LARGE FOR MODE=0' WRITE(*,*)'BETTER TRY MODE=1'
C ERR=4 WRITE(*,*)'M IS TOO LARGE FOR MODE=0' WRITE(*,*)'BETTER TRY MODE=1'
C ERR=5 WRITE(6,*)'IMPROPER ARGUMENT. Z MUST BE GREATER THAN 1'
C ERR=6 WRITE(*,*)'M IS TOO LARGE FOR MODE=1,2'*/

	switch(err)
	{
		case 0:		break;// all good
		case 1:  	printf("dtorh1: IPRE MUST BE 1 OR 2\n"); break;
		case 2:  	printf("dtorh1: YOU MUST CHOOSE MODE=2"); break;
		case 3:  	
		case 4:  	printf("dtorh1: M IS TOO LARGE FOR MODE=0'. WRITE(*,*)'BETTER TRY MODE=1\n"); break;
		case 5:  	printf("dtorh1: IMPROPER ARGUMENT. Z MUST BE GREATER THAN 1\n"); break;
		case 6:  	printf("dtorh1: M IS TOO LARGE FOR MODE=1,2\n"); break;
		case 10:	printf("dtorh1: Invalid mode (need 0,1,2)!\n"); break;
		default:	printf("dtorh1: unknown error!\n"); break;
	};
	
	return;
}

void HandleErrorDTORH2(int err)
{
/*
/*C ERR is set to inform of any errors during execution
C ERR=0 means OK
C ERR=1 WRITE(6,*)'IPRE MUST BE 1 OR 2'
C ERR=2 WRITE(*,*)'YOU MUST CHOOSE MODE=2'
C ERR=5 WRITE(6,*)'IMPROPER ARGUMENT. Z MUST BE GREATER THAN 1'
C ERR=6 WRITE(6,*)'MMAX IS TOO LARGE FOR MODE=0' WRITE(6,*)'BETTER TRY MODE=1'
*/
	switch(err)
	{
		case 0:		break;
		case 1:  	printf("dtorh2: IPRE MUST BE 1 OR 2\n"); break;
		case 2:  	printf("dtorh2: YOU MUST CHOOSE MODE=2"); break;
		case 5:  	printf("dtorh2: IMPROPER ARGUMENT. Z MUST BE GREATER THAN 1\n"); break;
		case 7:  	printf("dtorh2: 'MMAX IS TOO LARGE FOR MODE=0' WRITE(6,*)'BETTER TRY MODE=1'\n"); break;
		case 10:	printf("dtorh1: Invalid mode (need 0,1,2)!\n"); break;
		default:	printf("dtorh1: unknown error!\n"); break;
	};
	
	return;
}

void HandleErrorDTORH3(int err)
{
/*
C ERR is set to inform of any errors during execution
C ERR=0 means OK
C ERR=1 WRITE(6,*)'IPRE MUST BE 1 OR 2'
C ERR=2 WRITE(*,*)'YOU MUST CHOOSE MODE=2'
C ERR=5 WRITE(6,*)'IMPROPER ARGUMENT. Z MUST BE GREATER THAN 1'
C ERR=7 WRITE(6,*)'MMAX IS TOO LARGE FOR MODE=0' WRITE(6,*)'BETTER TRY MODE=1'
*/
	switch(err)
	{
		case 0:		break;
		case 1:  	printf("dtorh3: IPRE MUST BE 1 OR 2\n"); break;
		case 2:  	printf("dtorh3: YOU MUST CHOOSE MODE=2"); break;
		case 5:  	printf("dtorh3: IMPROPER ARGUMENT. Z MUST BE GREATER THAN 1\n"); break;
		case 7:  	printf("dtorh3: 'MMAX IS TOO LARGE FOR MODE=0' WRITE(6,*)'BETTER TRY MODE=1'\n"); break;
		case 10:	printf("dtorh1: Invalid mode (need 0,1,2)!\n"); break;
		default:	printf("dtorh1: unknown error!\n"); break;
	};
	
	return;
}

int main()
{
	double zVal=1.5;
	int mVal=0;
	int nMax=NMAX;
	int newN=0;
	
	int err=5;
	int mode=0;
	
	// some weird 'fails to converge' comes on if the 
	 // PLVec, QLVec is set to be what it is supposed to be
	 // increasing it by even 2 seems to fix the error. But i want to be on the safe side
	 // tested it all by running the author's program in gdb-debugger and 
	 // mine, and ensuring the equality between going-out given the same going-in params
	 double PLVec[NMAX*2+10]={0.0};
	 double QLVec[NMAX*2+10]={0.0};
	
	////// test DTORH1
	wrapDTORH1(zVal, mVal, nMax, PLVec, QLVec, &newN, &err, mode);
	HandleErrorDTORH1(err);
	
	printf("z=%1.5f\n", zVal);
	printf("PL[0]=%.5e,\tQL[0]=%.5e\n", PLVec[0], QLVec[0]);
	printf("PL[1]=%.5e,\tQL[1]=%.5e\n", PLVec[1], QLVec[1]);
	printf("PL[last]=%.5e,\tQL[last]=%.5e\n", PLVec[newN+1], QLVec[newN+1]);
	
	wrapDTORH1(0.5, mVal, nMax, PLVec, QLVec, &newN, &err, mode);
	HandleErrorDTORH1(err);
	
	////// test DTORH2
	// here dim means that fortran will use indices 0...dim, so there are dim+1 elements
	err=5;
	
	const int d2_mDim=51;
	const int d2_nDim=301;
	const int d2_mMax=50;
	const int d2_nMax=300;
	//
	int newM=0;
	int newNM[d2_mDim+1];
	// columns and rows are back-to-front in fortran
	double PL1Vec[(d2_mDim+1)*(d2_nDim+1)];
	double QL1Vec[(d2_mDim+1)*(d2_nDim+1)];
	
	wrapDTORH2(zVal, d2_mDim, d2_nDim, d2_mMax, d2_nMax, PL1Vec, QL1Vec, &newM, newNM, &err, mode);
	HandleErrorDTORH2(err);
	
	printf("\n-----------------\nDTORH2\n");
	printf("Z=%.3f\t newM=%d\t NM=%d\t PL(M,N)=%.5e\t QL[M,N]=%.5e PL[0,0]=%.5e\n", 
				zVal, newM, newNM[newM], 
				PL1Vec[FORT2D_TOC1D(newM,newNM[newM],d2_mDim)], 
				QL1Vec[FORT2D_TOC1D(newM,newNM[newM],d2_mDim)], 
				PL1Vec[FORT2D_TOC1D(0,0,d2_mDim)]);
				
	wrapDTORH2(-zVal, d2_mDim, d2_nDim, d2_mMax, d2_nMax, PL1Vec, QL1Vec, &newM, newNM, &err, mode);
	HandleErrorDTORH2(err);
				
	// test DTORH3
	err=5;
	wrapDTORH3(zVal, d2_mDim, d2_nDim, d2_mMax, d2_nMax, PL1Vec, QL1Vec, &newM, &newN, &err, mode);
	HandleErrorDTORH3(err);
	
	printf("\n-----------------\nDTORH3\n");
	printf("Z=%.3f\t newM=%d\t NM=%d\t PL(M,N)=%.5e\t QL[M,N]=%.5e PL[0,0]=%.5e\n", 
				zVal, newM, newNM[newM], 
				PL1Vec[FORT2D_TOC1D(newM,newN,d2_mDim)], 
				QL1Vec[FORT2D_TOC1D(newM,newN,d2_mDim)], 
				PL1Vec[FORT2D_TOC1D(0,0,d2_mDim)]);
	
	wrapDTORH3(-zVal, d2_mDim, d2_nDim, d2_mMax, d2_nMax, PL1Vec, QL1Vec, &newM, &newN, &err, mode);
	HandleErrorDTORH3(err);
	
	return 0;
}
