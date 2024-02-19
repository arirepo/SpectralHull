program elast2d
  use globals
  use grid_opt
  use trimesher
  use ps
  use spem_2d
  use fem_utils
  use plot_curved_elems
  use precond_opt
  use solidMat

  implicit none

  type(fem_struct) :: fem
  integer, dimension(:), allocatable :: pin, eltypein
  real*8 :: tolerance
  type(precond) :: tprecond
  type(isoMat2D)  :: alum

  call read_segment_file('../../geom/elas_plt_hole_tri_coarse_seg_mjm.dat', 'parabolic', fem%grd)
  call trigen('pnq34.0jY', fem%grd)

  ! Setting the p-refinement variables.
  allocate(pin(fem%grd%ntri), eltypein(fem%grd%ntri))
  pin = 6     !6th order element
  eltypein = 1     ! 1:Chebyshev (Fekete), 0:Lagrange
  tolerance = 1.0d-13 ! Spatial tolerance for two neighboring points.
  call add_more_points(fem%grd, pin, eltypein, tolerance)
  deallocate(pin, eltypein)

  !Setting elastic material properties
  alum%E = 10000000.0d0
  alum%nu = .33d0

  call init_mat2D(alum, state = 0) !state=0 --> plane stress

  ! Initialize the FEM object: calculate and assemble the K matrices
  ! and imposes boundary conditions.
  call fem_init(fem, neqs = 2, lin_solve_method = 'LAPACK_LU', tag = 1, fem_type = 1, alum)

  ! solver the Ax = b problem
  call fem_solve(fem, .True.)

  call write_u_tecplot('elast_p6.dat',fem%grd,fem%u)


  print *, 'Solution completed!'
end program elast2d
