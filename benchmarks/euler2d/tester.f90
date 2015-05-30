program benchmark_geom
  use globals
  use grid_opt
  use trimesher
  use dg_workspace
  use euler2d_eqs
  use spem_2d
  use plot_curved_elems
  use plot_curved_elems_dg
  use element_opt_dg2d
  implicit none

  ! local vars
  integer :: i
  type(dg_wspace) :: wspace
  integer, dimension(:), allocatable :: pin, eltypein
  real*8 :: tolerance, gamma0, rho0, u0, v0, P0
  real*8, dimension(4, 8) :: bc_vals0
  character(len = 128), dimension(8) :: bc_names0
  real*8, dimension(:, :), allocatable :: utmp
  type(fem_struct) :: tfem
  type(phys_props) :: tprops

  ! grid generation
  ! call read_segment_file('../../geom/naca0012_euler.dat', 'parabolic', wspace%grd)
  ! call read_segment_file('../../geom/cylinder_euler_tetrex.dat', 'parabolic', wspace%grd)
  !call read_segment_file('../../geom/naca0012_euler_tetrex.dat', 'parabolic', wspace%grd)
  ! call trigen('pnq32.0jY', wspace%grd)
  ! call trigen_based_TETREX('../../geom/cylinder_euler_tetrex.grd', wspace%grd)
  !call trigen_based_TETREX('../../geom/naca0012_euler_tetrex.grd', wspace%grd)
  ! call read_segment_file('../../geom/coarse_cylinder_tri2.dat', 'parabolic', wspace%grd)
!  call read_segment_file('../../geom/box_box.dat', 'parabolic', wspace%grd)
  call read_segment_file('../../geom/1box.dat', 'parabolic', wspace%grd)

  ! call read_segment_file('../../geom/small_cylinder.dat', 'parabolic', wspace%grd)
  ! call read_segment_file('../../geom/mms_square.dat', 'parabolic', wspace%grd)

  ! call read_segment_file('../../geom/chambered_airfoil.dat', 'parabolic', wspace%grd)
  call trigen('pnq36.0jY', wspace%grd)
  ! call trigen('pnq36.1jY', wspace%grd)

! print grid properties
! call print_grid_props(wspace%grd)
  ! write the skeleton grid ...
  allocate(utmp(1, wspace%grd%nnodesg))
  utmp = 0.0d0 
  call write_u_tecplot('p1grd.dat', wspace%grd, utmp)
  deallocate(utmp)

! stop
  ! proceed to add higher order points
  allocate( pin(wspace%grd%ntri), eltypein(wspace%grd%ntri) )
  pin = 6
  eltypein = 1
  tolerance = 1.0d-13
!    pin(4) = 3
! eltypein(4) = 1
! pin((/12, 22, 15, 2, 19, 16/)) = 3
! eltypein((/12, 22, 15, 2, 19, 16/)) = 1
! eltypein(26) = 0
! pin((/257, 15, 61, 5, 32/)) = 5

  call add_more_points(wspace%grd, pin, eltypein, tolerance, galtype = 'DG')

  ! print grid properties
  call print_grid_props(wspace%grd)

  ! write higher-order grid and data to view
  allocate(utmp(1, wspace%grd%nnodesg))
  utmp = 0.0d0 
  call write_u_tecplot('high_order_grd.dat', wspace%grd, utmp)
  deallocate(utmp)

  ! DG initialization
  gamma0 = 1.4d0
  rho0 = 0.1d0
  u0 = 0.2d0
  v0 = 0.0d0
  P0 = 0.8d0
  do i = 1, 4
     call u2U(rho0, u0, v0, P0, gamma0, bc_vals0(:, i))
  end do
  bc_names0(1) = 'outflow'
  bc_names0(2) = 'outflow'
  bc_names0(3) = 'outflow'
  bc_names0(4) = 'outflow'
  bc_names0(5) = 'outflow'
  bc_names0(6) = 'outflow'
  bc_names0(7) = 'outflow'
  bc_names0(8) = 'outflow'

  ! bc_names0(1) = 'wall'
  ! bc_names0(2) = 'wall'
  ! bc_names0(3) = 'outflow'
  ! bc_names0(4) = 'outflow'

  ! bc_names0(1) = 'outflow'
  ! bc_names0(2) = 'outflow'
  ! bc_names0(3) = 'outflow'
  ! bc_names0(4) = 'outflow'

  ! init workspace
  tprops%is_viscous = .true.
  tprops%mu = 1.d0
  tprops%lambda = -2.0d0/ 3.0d0
  tprops%Pr = 0.71d0
  tprops%adia = .true.

  call wspace%init(neqs = 4, gamma = gamma0 &
       , bc_names = bc_names0, bc_vals = bc_vals0, tol = tolerance, tprops = tprops)

  ! init the field
  call wspace%init_field(rho0, u0, v0, P0)

  ! march only one time step
  call wspace%march_field(dt = 1.0d-5, itrs = 1)
  ! call wspace%march_euler_implicit(dt = 4.0d-3, itrs = 100, inewtons = 2, num = 20, nrst = 1, epsil = 1.d-14)
! print *, 'heyhey'
  ! call wspace%march_field(dt = 1.0d-4, itrs = 30000)

  ! ! high-order visulization for continious curved elements
  ! ! NEEDS TO BE REPLACED with the ones for DG elements!
  ! tfem%grd = wspace%grd
  ! tfem%neqs = 4
  ! tfem%tag = 1
  ! ! allocate solution vars
  ! allocate(tfem%u(4,wspace%grd%nnodesg) &
  !        , tfem%dudx(4,wspace%grd%nnodesg) &
  !        , tfem%dudy(4,wspace%grd%nnodesg) )

  ! call wspace%udg2upg(tfem%u)
  ! call vis_curved_grid(tfem, 'curved.dat')

  ! write higher-order grid and data to view
  ! allocate(utmp(4, wspace%grd%nnodesg))
  ! call wspace%udg2upg(utmp)
  ! call write_u_tecplot('bad_high_order_grd.dat', wspace%grd, utmp)
  ! deallocate(utmp)

  call wspace%scatter_tecplot_dg('scatter.dat', .false., 'hey')
  call vis_curved_grid_dg(wspace, 'dg_vis.dat')
 

  ! celan ups
  deallocate(pin, eltypein)

  ! print
  print *, 'OK' 

end program benchmark_geom
