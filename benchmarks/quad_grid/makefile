### suffix rule
.SUFFIXES:
.SUFFIXES: .f90 .o .mod
%.o:%.f90
	$(GF) $(COMPFLAGS) -c $^ -o $@

### compiler
GF = /usr/bin/gfortran
VPATH = ./:../../src/:
COMMONFLAGS = -lblas -llapack
COMPFLAGS = -Wall -c -fbounds-check $(COMMONFLAGS)
LINKFLAGS = $(COMMONFLAGS)

### objects
OBJ = globals.o nurbs.o spline.o fekete.o ncc_triangle.o fem_reordering.o grid_opt.o quad_gen.o grid_vis.o tester.o

### compile and link
quadgrid: $(OBJ)
	$(GF) -o $@ $(LINKFLAGS) $(OBJ)
clean: $(OBJ)
	rm *.o *.mod
