SRC = $(wildcard modCor.f90 template.f90 )

all : 
	gfortran -O3 $(SRC) -o jsp #-Werror 
clean :
	rm *.mod jsp
