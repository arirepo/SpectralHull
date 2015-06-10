program benchmark_geom
  use globals
  use grid_opt
  use trimesher
  use dg_workspace
  use euler2d_eqs
  use spem_2d
  use plot_curved_elems
  use plot_curved_elems_dg
  implicit none

  ! local vars
  integer :: i, j, tpt
  type(dg_wspace) :: wspace
  integer, dimension(:), allocatable :: pin, eltypein
  real*8 :: tolerance, gamma0, rho0, u0, v0, P0
  real*8, dimension(4, 7) :: bc_vals0
  character(len = 128), dimension(7) :: bc_names0
  real*8, dimension(:, :), allocatable :: utmp
  type(fem_struct) :: tfem
  real*8 :: xx, yy
  logical :: adp_flag
  ! for animation
  character(len = 100) :: anim_name
  integer :: anim_itr

  ! grid generation
  ! call read_segment_file('../../geom/naca0012_euler.dat', 'parabolic', wspace%grd)
  ! call read_segment_file('../../geom/cylinder_euler_tetrex.dat', 'parabolic', wspace%grd)
  ! call read_segment_file('../../geom/naca0012_euler_tetrex.dat', 'parabolic', wspace%grd)
  ! call read_segment_file('../../geom/triangle_TETREX.dat', 'parabolic', wspace%grd)
  call read_segment_file('../../geom/triangle_TETREX2.dat', 'parabolic', wspace%grd)
  ! call trigen('pnq32.0jY', wspace%grd)
  ! call trigen_based_TETREX('../../geom/cylinder_euler_tetrex.grd', wspace%grd)
  !call trigen_based_TETREX('../../geom/naca0012_euler_tetrex.grd', wspace%grd)
  ! call trigen_based_TETREX('../../geom/triangle_TETREX.grd', wspace%grd)
  call trigen_based_TETREX('../../geom/triangle_TETREX2.grd', wspace%grd)
  ! call read_segment_file('../../geom/coarse_cylinder_tri2.dat', 'parabolic', wspace%grd)
  ! call read_segment_file('../../geom/chambered_airfoil.dat', 'parabolic', wspace%grd)
  ! call read_segment_file('../../geom/triangle.dat', 'parabolic', wspace%grd)

  ! call trigen('pnq35.0jY', wspace%grd)
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
  pin = 3
  eltypein = 0


  ! adpat <p> regionally
  do i = 1, size(wspace%grd%icon, 1) !ncellsg

     adp_flag = .true.
     do j = 1, size(wspace%grd%icon, 2) !npe

        tpt = wspace%grd%icon(i, j)
        xx = wspace%grd%x(tpt)
        yy = wspace%grd%y(tpt)
        adp_flag = adp_flag .and. (xx >= 0.3d0) .and. (xx <= 2.0d0) .and. (yy >= -0.4d0)  .and. (yy <= 0.4d0)
     end do

     if (adp_flag ) then
        pin(i) = 7
        eltypein(i) = 0
     end if

  end do


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
  ! rho0 = 0.1d0
  rho0 = 1.0d0
  u0 = 0.2d0
  v0 = 0.0d0
  ! P0 = 0.8d0
  P0 = 1.0d0 / gamma0

  do i = 1, 4
     call u2U(rho0, u0, v0, P0, gamma0, bc_vals0(:, i))
  end do
  bc_names0(1) = 'outflow'
  bc_names0(2) = 'outflow'
  bc_names0(3) = 'outflow'
  bc_names0(4) = 'outflow'
  bc_names0(5) = 'wall'
  bc_names0(6) = 'wall'
  bc_names0(7) = 'wall'

  ! bc_names0(1) = 'wall'
  ! bc_names0(2) = 'wall'
  ! bc_names0(3) = 'outflow'
  ! bc_names0(4) = 'outflow'

  ! bc_names0(1) = 'outflow'
  ! bc_names0(2) = 'outflow'
  ! bc_names0(3) = 'outflow'
  ! bc_names0(4) = 'outflow'

  ! init workspace
  call wspace%init(neqs = 4, gamma = gamma0 &
       , bc_names = bc_names0, bc_vals = bc_vals0, tol = tolerance)

  ! init the field
  call wspace%init_field(rho0, u0, v0, P0)

  ! march only one time step
  ! call wspace%march_field(dt = 1.0d-4, itrs = 300)
  do anim_itr = 1, 400

     ! compute the field evolution
     call wspace%tvd_rk(dt = 1.0d-4, itrs = 500)

     ! make animation
     write (anim_name, "(A5,I0.3)") "dgvis", anim_itr
     call vis_curved_grid_dg(wspace, anim_name)

  end do
  
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
