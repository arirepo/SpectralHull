### suffix rule
.SUFFIXES:
.SUFFIXES: .f90 .o .mod
%.o:%.f90
	$(GF) $(COMPFLAGS) -c $^ -o $@
%.o:%.f
	$(GF) -c $^ -o $@

### compiler
##GF = /usr/bin/gfortran
GF = ~/gcc-trunk/bin/gfortran
VPATH = ./:../../src/:
COMMONFLAGS = -lblas -llapack
COMPFLAGS = -Wall -c -fbounds-check $(COMMONFLAGS)
LINKFLAGS = $(COMMONFLAGS)

### objects
OBJ = splib.o globals.o condition.o dunavant.o nurbs.o spline.o fekete.o ncc_triangle.o fem_reordering.o quadri_elem.o approx_fekete.o grid_opt.o quad_gen.o grid_vis.o gen_basis.o element_opt.o bn_integral.o bcs.o lag_basis.o fem_utils.o spem_2d.o tester.o

### compile and link
quadgrid: $(OBJ)
	$(GF) -o $@ $(LINKFLAGS) $(OBJ)
clean: $(OBJ)
	rm *.o *.mod
