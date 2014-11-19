program benchmark_geom
  use globals
  use grid_opt
  use trimesher
  use ps
  use spem_2d
  use fem_utils
  use plot_curved_elems
  use mfree
  use mfree_st
  use reorder_opt
  use sparse_opt
  use sparse_assembly
  use precond_opt
  implicit none

  type(fem_struct) :: fem
  type(grid) :: grd2

  integer, dimension(:), allocatable :: pin, eltypein
  real*8 :: tolerance
  real*8, dimension(:,:), allocatable :: utmp
  integer :: ii, jj
  real*8 :: xx, yy, err
  type(mf_struct) :: tmf
  real*8, dimension(:), allocatable :: x, x2, x3, b, N, res
  real*8 :: epsil0
  integer :: num, nrst

  type(mf_struct_st) :: tmfst
  real*8, dimension(:), allocatable :: u0, du0
  real*8, dimension(:, :), allocatable :: rhs, xst, xst0
  real*8 :: dtt
  real*8 :: r0, theta0
  type(reorder) :: treorder
  type(smat) :: tsmat, tsmat2
  real*8, dimension(:, :), allocatable :: full_mat
  type(precond) :: tprecond
  integer :: flag, iter
  real*8 :: relres



  call read_segment_file('../../geom/coarse_cylinder_tri.dat', 'parabolic', fem%grd)
  call trigen('pnq34.0jY', fem%grd)
  allocate( pin(fem%grd%ntri), eltypein(fem%grd%ntri) )
  pin = 6
  eltypein = 1
  ! eltypein(1) = 0
  tolerance = 1.0d-13
  call add_more_points(fem%grd, pin, eltypein, tolerance)
  deallocate(pin, eltypein)

  call fem_init(fem, neqs = 1, lin_solve_method = 'LAPACK_LU', tag = 1) 

  ! sparsity pattern stuff and preconditioning
  call treorder%init(fem%grd)
  call treorder%spy('spy.dat')
  call treorder%reorder(fem%grd)
  call treorder%spy('spy_after.dat')

  call tsmat%init(treorder%rows, 1)
  call tsmat2%init(treorder%rows, 1)
  call sparse_assem(fem%grd, fem%elems, tsmat)
  call sparse_assem_M(fem%grd, fem%elems, tsmat2)
  call tsmat%sp2full(full_mat)

  print *, ' maxval(abs(K_full - K_sparse)) = ', maxval(abs(fem%KK_mat - full_mat))


  call tmf%init(fem%grd, fem%elems, (/1, 2, 3, 4 /) &
       , (/1.0d0, 1.0d0, 2.0d0, 2.0d0 /) )


  ! space-time
  !
  !
  dtt = 0.12d0
  num = 20 !size(x)
  nrst = 30 !1


  allocate(u0(fem%neqs*fem%grd%nnodesg) , du0(fem%neqs*fem%grd%nnodesg) )
  u0 = 0.0d0; du0 = 0.0d0

  do ii = 1, fem%grd%nnodesg
     r0 = sqrt(fem%grd%x(ii)**2.0d0 + fem%grd%y(ii)**2.0d0)
     theta0 = acos(fem%grd%x(ii) / r0 )
     u0(ii) = exp(-4.0d0 * (r0-2.0d0)**2.0d0) * exp(-4.0d0 * (theta0-PI)**2.0d0)

     ! u0(ii) = exp(-25.0d0 * ((fem%grd%x(ii) - 0.5d0)**2.0d0 &
     !      + (fem%grd%y(ii)-0.5d0)**2.0d0))
  end do
  call init_cond(fem%grd, (/1, 2/), (/0.0d0, 0.0d0/) &
       , 'none', u0)

  call tmfst%init_st(fem%grd, fem%elems, (/1, 2/) &
       , (/0.0d0, 0.0d0 /) &
       , 5, -0.5d0, -0.5d0, 0.0d0, dtt, u0, du0)

  allocate(rhs(fem%neqs*fem%grd%nnodesg, tmfst%nS) &
       , xst(fem%neqs*fem%grd%nnodesg, tmfst%nS) &
       , xst0(fem%neqs*fem%grd%nnodesg, tmfst%nS))
  if(allocated(N)) deallocate(N)
  allocate(N(size(tmfst%free) * tmfst%nS))
  N = 1.0d0

  call tprecond%init(tsmat2, tmfst%zeta, tsmat, tmfst%free, 4, 4)
  ! tprecond%enabled = .false. 
  epsil0 = 1.0d-12

  do ii = 1, 100
     ! test call rhs
     call tmfst%rhs(rhs)
     do jj = 1, tmfst%nS
        xst(:, jj) = u0
     end do
     if(allocated(res)) deallocate(res)
     ! call tmfst%gmres(rhs, tprecond, epsil0, num, nrst, xst, res)
     ! call tmfst%stbicg(rhs, tprecond, epsil0, 2000, xst, res)
     xst0 = xst
     call tmfst%idrs(rhs, 10, epsil0, 1000, tprecond, xst0 &
          , 1, xst, flag, res, iter, relres)

     call tmfst%window(xst)

     ! ! print the residuals
     print *, 'ii = ', ii, maxval(xst(:, tmfst%nS)) , minval(xst(:, tmfst%nS)), minval(res)
  end do

   
  fem%u = reshape(xst(:, tmfst%nS), (/ 1, tmfst%grd%nnodesg /) )
  call vis_curved_grid(fem, 'ahah.dat')
  call scatter_tecplot('reordered_scattered.dat', .false., 'sample' &
       , fem%grd%x, fem%grd%y, fem%u)
  ! call adaptive_vis_curved_grid(fem, 'gg.dat', 1.0d-14, 1.0d-4)
        call adaptive_vis_curved_grid(fem, 'gg.dat' &
             , 1.0d-13, 1.0d-5, (/0, 100000/), 10.0d0, 1.0d-6, '3pt')

  ! print
  print *, 'OK' 


end program benchmark_geom
