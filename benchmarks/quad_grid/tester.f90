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
  ! ! call quad_grid_gen('../../geom/segments_tmp.dat', elem, x, y, nq, nt)
  ! ! call quad_grid_gen('../../geom/naca_edges.dat', elem, x, y, nq, nt)
  ! call quad_grid_gen('../../geom/coarse_cylinder_tri.dat', elem, elemidx, etype, x, y, nq, nt)

  ! call printout(elem, x, y, nq, nt)
  ! deallocate(elem, x, y, elemidx, etype)

  ! testing type(grid) convertor
  call quadgen('../../geom/coarse_cylinder_tri.dat', grd, 1)
  call print_grid_props(grd)

  allocate(u(1, grd%nnodesg))
  u = 1.0d0
  call write_u_tecplot('./grd.dat', grd, u)

  call visual_curved_bonds('./bnd.dat', grd, 10)

  ! done here
  print *, 'OK'

end program tester
