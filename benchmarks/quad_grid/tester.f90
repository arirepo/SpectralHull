program tester
  use quad_gen
  use grid_vis

  implicit none

  integer, parameter :: rk = selected_real_kind(12)
  integer, dimension(:), allocatable :: elem, elemidx, etype
  real(rk), dimension(:), allocatable :: x, y
  integer :: nq, nt

! call quad_grid_gen('../../geom/segments_tmp.dat', elem, x, y, nq, nt)
! call quad_grid_gen('../../geom/naca_edges.dat', elem, x, y, nq, nt)
  call quad_grid_gen('../../geom/coarse_cylinder_tri.dat', elem, elemidx, etype, x, y, nq, nt)

  call printout(elem, x, y, nq, nt)
  deallocate(elem, x, y, elemidx, etype)
end program tester
