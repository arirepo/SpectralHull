#makefile for "Terrible" Spectral Element Solver  
# by Arash Ghasemi and Mary Barker
### NOTE!!! Please leave this section unaltered
###############################################
VPATH = ./:../../src:
%.o: %.f90
	$(F90) $(COMPFLAGS) $^ -o $@
%.o: %.f
	$(F77) -c  $(OMP_FLAGS) $^ -o $@
%.o: %.c
	$(CC) $(CCOMPILEFLAGS) $^ -o $@ 
###############################################
###############################################
## NOTE!!! Now select your compiler from below 
### compiler
# F90 = /usr/bin/gfortran
#F90 = /opt/gcc/bin/gfortran
F90 = /home/aghasemi/gcc-trunk/install_dir/bin/gfortran
# F90 = /usr/local/gcc-5-20150306/bin/gfortran
F77 = $(F90)
CC = cc
# CC = /usr/local/gcc-5-20150306/bin/gcc
# CC = icc
##F90 = gfortran
# F90 = ifort
# COMMONFLAGS  =  -fbounds-check -Wall -Wconversion -Wunderflow -L $(LIBS) -I $(INCS)
# COMMONFLAGS  =  -O3 -L $(LIBS) -I $(INCS)
# COMMONFLAGS  = -O3
# COMMONFLAGS  = 
# CCOMPILEFLAGS  =  -c 
COMMONFLAGS  = $(OMP_FLAGS) -O3 
CCOMPILEFLAGS  =  -c $(OMP_FLAGS) 
##COMMONFLAGS  = 
##COMMONFLAGS  =  -O3 -check all 
# COMMONFLAGS  =  -O3 -parallel  
# COMMONFLAGS  =  $(OMP_FLAGS) -O3 -xHost -ipo
COMPFLAGS    =  -c  $(COMMONFLAGS)
LIBS    = 
INCS    = 
LINKFLAGS    =      $(COMMONFLAGS) 
LINK_SUFFIX  =     -llapack -lblas
#LINK_SUFFIX  =   -L/home/aghasemi/blueprint_lapack/gfort_lapack/lapack-3.5.0 -llapack -lblas

OMP_FLAGS = -fopenmp
# OMP_FLAGS = 

## NOTE!!! Now add your object files that you need for 
#          your test case problem. 

### objects
OBJ = renka_trimesh_lib.o renka_trimesh.o polygauss.o quads.o euler2d_eqs.o triangle.o custom_Tri_wrapper.o condition.o fekete.o ncc_triangle.o globals.o nurbs.o spline.o splib.o ortho_poly.o fem_reordering.o quadri_elem.o approx_fekete.o grid_opt.o grd2hull.o lag_basis.o fem_utils.o dunavant.o gen_basis.o element_opt.o element_opt_dg2d.o bn_integral.o bcs.o spem_2d.o interp.o ps.o trimesher.o triangulation_quad.o dg_workspace.o plot_curved_elems.o plot_curved_elems_dg.o tester.o   
###

### compile and link
EXE = exe.run

all: $(EXE)
exe.run:  $(OBJ)
	$(F90)  -o $@ $(LINKFLAGS)  $(OBJ) $(LINK_SUFFIX)
clean:
	rm *.run *~ *# $(OBJ) *.phy *.mod *.o
