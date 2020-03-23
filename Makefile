SRC = $(wildcard modCor.f90 modCommon.f90 jsub.f90 test.f90 )

all : 
	gfortran -O3 $(SRC) -o jsp #-Werror 
clean :
	rm *.mod jsp
