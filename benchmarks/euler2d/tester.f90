program benchmark_geom
  use globals
  use grid_opt
  use trimesher
  use spem_2d
  use plot_curved_elems
  use element_opt_dg2d
  implicit none

  ! local vars
  type(fem_struct) :: fem
  integer, dimension(:), allocatable :: pin, eltypein
  real*8 :: tolerance
  type(element_dg2d) :: telemdg

  ! grid generation
  call read_segment_file('../../geom/simp_rect_offset.dat', 'parabolic', fem%grd)
  call trigen('pnq34.0jY', fem%grd)
  allocate( pin(fem%grd%ntri), eltypein(fem%grd%ntri) )
  pin = 2
  eltypein = 0
  tolerance = 1.0d-13
  call add_more_points(fem%grd, pin, eltypein, tolerance, galtype = 'DG')

  ! fem initialization
  call fem_init(fem, neqs = 1, lin_solve_method = 'LAPACK_LU', tag = 1) 
  ! call vis_curved_grid(fem, 'euler_grd.dat')

  ! test dg elements
  call telemdg%init(telemdg%element, 1, fem%grd, 4, 1.4d0)

  ! celan ups
  deallocate(pin, eltypein)

  ! print
  print *, 'OK' 

end program benchmark_geom
