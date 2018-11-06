:: clean up
del wrapDTORH.dll
del wrapDTORH64.dll
del test.exe
:: compile wrapper
gfortran -c -g -O0 dtorh1_mod0.f &&^
gfortran -c -g -O0 dtorh1_mod1.f &&^
gfortran -c -g -O0 dtorh1_mod2.f &&^
::
gfortran -c -g -O0 dtorh2_mod0.f &&^
gfortran -c -g -O0 dtorh2_mod1.f &&^
gfortran -c -g -O0 dtorh2_mod2.f &&^
::
gfortran -c -g -O0 dtorh3_mod0.f &&^
gfortran -c -g -O0 dtorh3_mod1.f &&^
gfortran -c -g -O0 dtorh3_mod2.f &&^
::
gfortran -c -g -O0 rout_mod.f &&^
gcc -c -g wrapDTORH1.c &&^
gcc -c -g wrapDTORH2.c &&^
gcc -c -g wrapDTORH3.c &&^
gcc dtorh1_mod0.o dtorh1_mod1.o dtorh1_mod2.o ^
	dtorh2_mod0.o dtorh2_mod1.o dtorh2_mod2.o ^
	dtorh3_mod0.o dtorh3_mod1.o dtorh3_mod2.o ^
	rout_mod.o ^
	wrapDTORH1.o wrapDTORH2.o wrapDTORH3.o ^
	-g -llapack -lblas -lgfortran -shared -o wrapDTORH64.dll
:: compile c tester
gcc test.c wrapDTORH64.dll -g -o test.exe
:: launch c tester
test.exe
pause