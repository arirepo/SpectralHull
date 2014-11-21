### suffix rule
.SUFFIXES:
.SUFFIXES: .f90 .o .mod
%.o:%.f90
	$(GF) $(COMPFLAGS) -c $^ -o $@

### compiler
GF = gfortran
VPATH = ./:../../src/:
COMMONFLAGS = -lblas -llapack
COMPFLAGS = -c -fbounds-check $(COMMONFLAGS)
LINKFLAGS = $(COMMONFLAGS)

### objects
OBJ = globals.o spline.o nurbs.o grid_wrapper.o quad_gen.o grid_vis.o tester.o

### compile and link
quadgrid: $(OBJ)
	$(GF) -o $@ $(LINKFLAGS) $(OBJ)
clean: $(OBJ)
	rm *.o *.mod