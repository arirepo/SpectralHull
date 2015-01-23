program tester
  use quad_gen
  use grid_vis
  use grid_opt
  use spem_2d
  use plot_curved_elems

  implicit none

  integer, parameter :: rk = selected_real_kind(12)
  integer, dimension(:), allocatable :: elem, elemidx, etype
  real(rk), dimension(:), allocatable :: x, y
  real*8, dimension(:,:), allocatable :: u
  integer :: nq, nt
  ! type(grid) :: grd
  real(rk) :: dx, dy, idx, idy
  integer, dimension(:), allocatable :: pin, eltypein
  type(fem_struct) :: fem

  ! dx = 3.0_rk
  dx = 3.0_rk
  idx = dx
  dy = dx
  idy = dy

  ! ! call quad_grid_gen('../../geom/segments_tmp.dat', elem, x, y, nq, nt)
  ! ! call quad_grid_gen('../../geom/naca_edges.dat', elem, x, y, nq, nt)
  ! call quad_grid_gen('../../geom/coarse_cylinder_tri.dat', elem, elemidx, etype, x, y, nq, nt)

  ! call printout(elem, x, y, nq, nt)
  ! deallocate(elem, x, y, elemidx, etype)

  ! testing type(grid) convertor
! call quadgen('../../geom/coarse_cylinder_tri.dat', grd, 1)
  call quadgen('../../geom/unit_quadri_for_Cp.dat', fem%grd, 2, dx, dy, idx, idy)
  ! call quadgen('../../geom/circles.dat', fem%grd, 1, dx, dy, idx, idy, 1)
  ! call quadgen('../../geom/simp_rect.dat', fem%grd, 1, dx, dy, idx, idy, 1)

  call print_grid_props(fem%grd)
  call visual_curved_bonds('./bnd.dat', fem%grd, 10)
  allocate(u(1, fem%grd%nnodesg))
  u = 1.0d0
  call write_u_tecplot('./grd_P1.dat', fem%grd, u)
  deallocate(u)


  ! subroutine add_more_points(grd, pin, eltypein, tolerance, galtype, n1, n2)
  allocate(pin(fem%grd%ncellsg), eltypein(fem%grd%ncellsg))
  pin = 8
  eltypein = 1
  call add_more_points(fem%grd, pin, eltypein, 1.0d-12, 'PG')

  call fem_init(fem, neqs = 1, lin_solve_method = 'LAPACK_LU', tag = 1) 

  call fem_solve(fem, echo = .true.)

  call write_u_tecplot('./grd.dat', fem%grd, fem%u)

  call vis_curved_grid(fem, 'mixed.dat')

  ! done here
  print *, 'OK'

end program tester
