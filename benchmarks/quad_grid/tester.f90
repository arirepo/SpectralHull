program tester
  use quad_gen
  use grid_vis
  use grid_opt

  implicit none

  integer, parameter :: rk = selected_real_kind(12)
  integer, dimension(:), allocatable :: elem, elemidx, etype
  real(rk), dimension(:), allocatable :: x, y
  real*8, dimension(:,:), allocatable :: u
  integer :: nq, nt
  type(grid) :: grd
  real(rk) :: dx, dy, idx, idy
  integer, dimension(:), allocatable :: pin, eltypein

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
  call quadgen('../../geom/circles.dat', grd, 1, dx, dy, idx, idy)
  call print_grid_props(grd)
  call visual_curved_bonds('./bnd.dat', grd, 10)


  ! subroutine add_more_points(grd, pin, eltypein, tolerance, galtype, n1, n2)
  allocate(pin(grd%ncellsg), eltypein(grd%ncellsg))
  pin = 6
  eltypein = 1
  call add_more_points(grd, pin, eltypein, 1.0d-12, 'PG')

  allocate(u(1, grd%nnodesg))
  u = 1.0d0
  call write_u_tecplot('./grd.dat', grd, u)

  ! done here
  print *, 'OK'

end program tester
