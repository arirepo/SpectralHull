program tester
  use quad_gen
  use grid_vis
  use grid_opt
  use spem_2d
  use plot_curved_elems
  use mfree

  implicit none

  integer, parameter :: rk = selected_real_kind(12)
  real*8, dimension(:,:), allocatable :: u
  ! type(grid) :: grd
  real(rk) :: dx, dy, idx, idy
  integer, dimension(:), allocatable :: pin, eltypein
  type(fem_struct) :: fem
  type(mf_struct) :: tmf
  real*8, dimension(:), allocatable :: x2, b, N, res
  real*8 :: epsil0
  integer :: num, nrst

  ! dx = 3.0_rk
  dx = 3.0_rk
  idx = dx
  dy = dx
  idy = dy

  ! generate quad-tri grid ...
  call quadgen('../../geom/unit_quadri_for_Cp.dat', fem%grd, 2, dx, dy, idx, idy, 1)

  ! show the generated grid properties and plot grid and boundaries
  call print_grid_props(fem%grd)
  call visual_curved_bonds('./bnd.dat', fem%grd, 10)
  allocate(u(1, fem%grd%nnodesg))
  u = 1.0d0
  call write_u_tecplot('./grd_P1.dat', fem%grd, u)
  deallocate(u)


  ! add high-order points ...
  !
  ! subroutine add_more_points(grd, pin, eltypein, tolerance, galtype, n1, n2)
  allocate(pin(fem%grd%ncellsg), eltypein(fem%grd%ncellsg))
  pin = 6
  eltypein = 0
  call add_more_points(fem%grd, pin, eltypein, 1.0d-12, 'PG')

  ! initalize the Finite Element generic object
  call fem_init(fem, neqs = 1, lin_solve_method = 'LAPACK_LU', tag = 1) 

  ! solve using full matrix assembled approach and write the solution 
  call fem_solve(fem, echo = .true.)
  call write_u_tecplot('./full_grd.dat', fem%grd, fem%u)
  call vis_curved_grid(fem, 'full_matrix.dat')

  ! test matrix-free implementation for Laplace equation
  call tmf%init(fem%grd, fem%elems, (/1, 2/) &
       , (/1.0d0, 2.0d0/) ) ! tags 1,2 are fixed with Dirichlet, the rest are Nuimann

  ! prepare initial condition
  allocate(x2(fem%grd%nnodesg))
  x2 = 0.0d0
  call init_cond(fem%grd, tmf%tags, tmf%values &
       , 'none', x2, bnval)

  ! parameters for GMRES
  epsil0 = 1.0d-15
  num = 100 !num. inner loops
  nrst = 20 !n-restarts
  allocate(b(size(x2)), N(size(tmf%free)))
  b = 0.0d0 ! rhs
  call init_cond(fem%grd, tmf%tags, tmf%values &
       , 'none', b, bnval)
  N = 1.0d0 ! preconditioner
 
  ! subroutine gmres_orig(mf, b, N, epsil, num, nrst, sol, res)
  call gmres_orig(tmf, b, N, epsil0, num, nrst, x2, res)
  print *, 'GMRES_res = ', res

  ! quickly find the error between full matrix approach and matrix free
  print *, 'Linf(error(Full matrix - Matrix free))' &
       , maxval(abs(fem%u - reshape(x2, (/1, fem%grd%nnodesg /))))

  ! plot the matrix free results
  fem%u = reshape(x2, (/1, fem%grd%nnodesg /) )
  call vis_curved_grid(fem, 'u_mfree.dat')

  ! clean ups
  if(allocated(u)) deallocate(u)
  if(allocated(pin)) deallocate(pin)
  if(allocated(eltypein)) deallocate(eltypein)
  if(allocated(x2)) deallocate(x2)
  if(allocated(b)) deallocate(b)
  if(allocated(N)) deallocate(N)
  if(allocated(res)) deallocate(res)

  ! done here
  print *, 'OK'

contains

function bnval(grd, pt)
implicit none
type(grid), intent(in) :: grd
integer, intent(in) :: pt
real*8 :: bnval

bnval = -1.0d0 * grd%x(pt)

! done here
end function bnval


end program tester
