#makefile for "Terrible" Spectral Element Solver  
# by Arash Ghasemi and Mary Barker
### NOTE!!! Please leave this section unaltered
###############################################
VPATH = ../../src:
%.o: %.f90
	$(F90) $(COMPFLAGS) $^ -o $@
%.o: %.f
	$(F90) -c $^ -o $@
%.o: %.c
	$(CC) $(CCOMPILEFLAGS) $^ -o $@ 
###############################################
###############################################
## NOTE!!! Now select your compiler from below 
### compiler
# F90 = /usr/bin/gfortran
# F90 = /opt/gcc/bin/gfortran
F90 = /home/aghasemi/gcc-trunk/bin/gfortran
CC = cc
# CC = icc
##F90 = gfortran
# F90 = ifort
COMMONFLAGS  =  -O3 -fbounds-check -Wall -Wconversion -Wunderflow -L $(LIBS) -I $(INCS)
# COMMONFLAGS  =  -O1 -L $(LIBS) -I $(INCS)
# COMMONFLAGS  = -O3
CCOMPILEFLAGS  =  -c 
##COMMONFLAGS  = 
##COMMONFLAGS  =  -O3 -check all 
#COMMONFLAGS  =  -O3 -parallel  
COMPFLAGS    =  -c  $(COMMONFLAGS)
LIBS    = 
INCS    = 
LINKFLAGS    =      $(COMMONFLAGS) 
LINK_SUFFIX  =     -llapack -lblas

## NOTE!!! Now add your object files that you need for 
#          your test case problem. 

### objects
OBJ = triangle.o custom_Tri_wrapper.o condition.o rev_cuthill_mckee.o fekete.o ncc_triangle.o globals.o nurbs.o spline.o splib.o ortho_poly.o fem_reordering.o quadri_elem.o approx_fekete.o grid_opt.o lag_basis.o fem_utils.o dunavant.o gen_basis.o element_opt.o bn_integral.o bcs.o normal_rand.o sparse_opt.o precond_opt.o reorder_opt.o sparse_assembly.o mfree.o mfree_st.o spem_2d.o interp.o ps.o trimesher.o plot_curved_elems.o  benchmark_geom.o   
###

### compile and link
EXE = exe.run

all: $(EXE)
exe.run:  $(OBJ)
	$(F90)  -o $@ $(LINKFLAGS)  $(OBJ) $(LINK_SUFFIX)
clean:
	rm *.run *~ *# $(OBJ) *.phy *.mod *.o
