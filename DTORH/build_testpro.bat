:: clean up
del testpro.exe
:: compile and run
gfortran -g dtorh1.f dtorh2.f dtorh3.f rout.f testpro.f -o testpro -lblas -llapack
testpro.exe
pause