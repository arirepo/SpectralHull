module dg_workspace
  use grid_opt
  use euler2d_eqs
  use spline
  use element_opt_dg2d
  implicit none

  private

  type bc

     character(len = 128) :: name
     ! one value per tag at this version
     ! val(1:neqs)
     real*8, dimension(:), allocatable :: val 

  end type bc

  type dg_wspace
     private
     real*8 :: dt
     type(grid), public :: grd
     type(element_dg2d), dimension(:), allocatable, public :: elems
     type(bc), dimension(:), allocatable :: bcs ! zero-based

   contains
     procedure, public :: init => init_wspace 
     procedure, public :: init_edg_quadrat => init_elem_edg_quadratures
     procedure, public :: comp_Fstar
     procedure, public :: comp_elem_rhs
     procedure, public :: init_field
     procedure, public :: march_field
     procedure, public :: udg2upg
     procedure, public :: scatter_tecplot_dg
     procedure, public :: comp_dFpm
     procedure, public :: comp_jac_bnd_integral
     procedure, public :: comp_full_Ax
     procedure, public :: march_euler_implicit
  
  end type dg_wspace


  public :: dg_wspace

contains

  ! initilizes the workspace
  ! HINT : the grid in wspace%grd should already
  ! been allocated and initialized properly
  subroutine init_wspace(wspace, neqs, gamma, bc_names, bc_vals, tol, tprops)
    implicit none
    class(dg_wspace), intent(inout) :: wspace
    integer, intent(in) :: neqs
    real*8, intent(in) :: gamma
    character(len = 128), dimension(:), intent(in) :: bc_names
    real*8, dimension(:, :), intent(in) :: bc_vals
    real*8, intent(in) :: tol
    class(phys_props), intent(in) :: tprops

    ! local vars
    integer :: i, nbcs, j

    ! bullet proofing ...
    ! check see if the grid wrapper is initialized/allocated
    ! appropriately 
    if ( (wspace%grd%nnodesg <= 0) .or. &
         (.not. allocated(wspace%grd%icon)) ) then
       print *, '<wspace%grd> not allocated/initialized! stop'
       stop
    end if

    ! allocate all elements in the workspace
    allocate(wspace%elems(wspace%grd%ncellsg))

    wspace%elems(:)%local_props%is_viscous = tprops%is_viscous
    wspace%elems(:)%local_props%mu = tprops%mu
    wspace%elems(:)%local_props%lambda = tprops%lambda
    wspace%elems(:)%local_props%Pr = tprops%Pr
    wspace%elems(:)%local_props%adia = tprops%adia

    ! proceed to initialize the elements one by one
    do i = 1, size(wspace%elems)
       call wspace%elems(i)%init(wspace%elems(i) &
            , i, wspace%grd, neqs, gamma)

    end do

    ! initializing boundary conditions
    nbcs = size(bc_names)
    ! this is a zero-based array and the zeroth entry is just
    ! showing that it belongs to the interior edge. 
    allocate(wspace%bcs(0:nbcs))
    wspace%bcs(0)%name = 'interior'
    do i = 1, nbcs
       allocate(wspace%bcs(i)%val(neqs))
       wspace%bcs(i)%val = bc_vals(:, i)
       wspace%bcs(i)%name = trim(bc_names(i))
    end do

    ! initializing the edges and tneigh(s) on the edges
    do i = 1, size(wspace%elems)
       do j = 1, wspace%elems(i)%nedgs
          call wspace%init_edg_quadrat(wspace%elems(i), j, tol)
       end do
    end do

    ! done here
  end subroutine init_wspace

  subroutine init_elem_edg_quadratures(wspace, elem, iedg, tol)
    implicit none
    class(dg_wspace), intent(inout), target :: wspace ! work space
    class(element_dg2d), intent(inout), target :: elem
    integer, intent(in) :: iedg !edge number
    real*8, intent(in) :: tol

    ! local vars
    integer :: i, j, max_npedg, r, tag
    type(edg_dg), pointer :: tedg => null()
    type(neigh_dg), pointer :: tneigh => null()
    real*8, parameter :: alpha = 0.0d0, beta  = 0.0d0
    logical :: do_snap = .false.
    type(curve), pointer :: tcurve=>null()
    real*8, dimension(:), pointer :: t=>null(), xs=>null(), ys=>null()
    real*8 :: tt, x1, y1, x2, y2, xq, yq, xdot, ydot
    type(grid), pointer :: grd => null()
    real*8 :: t1, t2
    real*8, dimension(:), allocatable :: derjac

    ! HARD reset
    do_snap = .false.

    grd => wspace%grd
    tedg => elem%edgs(iedg)
    tag = tedg%tag

    if ( tag .ne. 0 ) then
       ! is a boundary element 
       do_snap = .true.
       tcurve => grd%bn_curves(tag)
       t => tcurve%t
       xs => tcurve%x
       ys => tcurve%y
    end if

    do j = 1, size(tedg%neighs) ! loop over all neighbors to that edge (shared segments)

       tneigh => tedg%neighs(j) ! select this neighbor 

       ! find the maximum number of 1d interpolation (Lagrange) points per this shared 
       ! segment with the j^th neighbor of this edge.
       !
       ! HINT : use npedg info to see how many points per edge we have in this elem
       ! and this neighbor. and take the maximum as 
       ! the highest possible degree of a 1d polynomial that can be defined on that
       ! part of the shared segment 
       ! 
       if (tneigh%elnum .eq. -1)  then !wall
          max_npedg = elem%npedg ! count on element itself
       else
          max_npedg = max(elem%npedg, wspace%elems(tneigh%elnum)%npedg)
       end if

       ! degree of exactness relation for 1d Gauss legendre: 2r - 1 = max_npedg -1
       r = 2 * elem%p + 1 ! 3.5 is safety factor :)
       tneigh%ngpseg = r

       ! allocate neigh specific arrays
       allocate(tneigh%xi(r), tneigh%W(r))
       allocate(tneigh%xloc_in(2, r), tneigh%xloc_out(2, r))
       allocate(tneigh%x(2, r), tneigh%dx(2, r))
       allocate(tneigh%n(2, r), tneigh%s(r))
       allocate(tneigh%psi_in(elem%npe, r))
       ! alloc/init Fstar(1:neqs, 1:ngpseg) 
       allocate(tneigh%Fstar(elem%neqs, r))
       ! alloc/init dFpm(1:neqs, 1:neqs, 2, 1:ngpseg) 
       allocate(tneigh%dFpm(elem%neqs, elem%neqs, 2, r))

       ! alloc/init Hs(1:neqs, 1:ngpseg) 
       allocate(tneigh%Hs(elem%neqs, r))

       tneigh%xi = 0.0d0; tneigh%W = 0.0d0
       tneigh%xloc_in = 0.0d0; tneigh%xloc_out = 0.0d0
       tneigh%x = 0.0d0; tneigh%dx = 0.0d0
       tneigh%n = 0.0d0; tneigh%s = 0.0d0
       tneigh%psi_in = 0.0d0
       tneigh%Fstar = 0.0d0
       tneigh%dFpm = 0.0d0
       tneigh%Hs = 0.0d0

       ! computing Legendre-Gauss-Jacobi points for integration
       ! and corresponding weight functions
       allocate(derjac(r))
       call ZEJAGA(r, alpha, beta, tneigh%xi, derjac)
       call WEJAGA(r, alpha, beta, tneigh%xi, derjac, tneigh%W)
       deallocate(derjac)

       ! 1 is the start and 2 is the end of this edge segment; just for convention
       x1 = tneigh%xs(1); x2 = tneigh%xe(1) 
       y1 = tneigh%xs(2); y2 = tneigh%xe(2) 

       ! find the constant derivatives
       xdot = 0.5d0 * (x2 - x1)
       ydot = 0.5d0 * (y2 - y1)

       do i = 1, r

          ! Mapping of coordinates of Gauss-Legendre quadrature
          ! for straight edges to physical space
          xq = x1 + (tneigh%xi(i) + 1.0d0) / 2.0d0 * (x2 - x1) 
          yq = y1 + (tneigh%xi(i) + 1.0d0) / 2.0d0 * (y2 - y1) 

          ! if boundary edge then also snapp <x> and <dx>
          if ( do_snap ) then
             call find_t(grd, tag, x1, y1, tol, t1)
             if ( t1 .eq. -1.0d0) then
                print *, 'error : the parameter t can not be computed' &
                     , ' for point (', x1, ',', y1,') on ' &
                     , 'edg # ', iedg,'. stop.'
                stop
             end if
             call find_t(grd, tag, x2, y2, tol, t2)
             if ( t2 .eq. -1.0d0) then
                print *, 'error : the parameter t can not be computed' &
                     , ' for point (', x2, ',', y2,') on ' &
                     , 'edg # ', iedg,'. stop.'
                stop
             end if

             call find_t(grd, tag, xq, yq, tol, tt)
             if ( tt .eq. -1.0d0) then
                print *, 'error : the parameter t can not be computed' &
                     , ' for point (', xq, ',', yq,') on ' &
                     , 'edg # ', iedg,'. stop.'
                stop
             end if
             ! find new snapped xq and yq on the curve
             xq = 0.0d0; yq = 0.0d0
             call spline_nurbs_eval2(tt, xs, ys, tcurve%a, tcurve%b, tcurve%c, &
                  & tcurve%Mx, tcurve%My, t, xq, yq, 'interp', tcurve%btype)

             ! computing derivatives xdot and ydot ...
             xdot = 0.0d0; ydot = 0.0d0
             call spline_nurbs_eval2(tt, xs, ys, tcurve%a, tcurve%b, tcurve%c, &
                  & tcurve%Mx, tcurve%My, t, xdot, ydot, 'diff1', tcurve%btype)

             ! add stuff coming from chain rule for differentiation
             xdot = xdot * 0.5d0 * (t2 - t1)
             ydot = ydot * 0.5d0 * (t2 - t1)

          end if

          ! store physical coords and derivatives
          tneigh%x(1,i) = xq; tneigh%x(2,i) = yq
          tneigh%dx(1,i) = xdot; tneigh%dx(2,i) = ydot
          tneigh%s(i) = sqrt(xdot**2 + ydot**2)
          tneigh%n(1, i) = ydot / tneigh%s(i)
          tneigh%n(2, i) = -xdot / tneigh%s(i) 

          ! convert physical (xq, yq) to local coords
          ! (r,s) in this element store in this element
          !
          ! xy2rs_dg(elem, x, y, maxitr, tolrs, r, s)
          call elem%xy2rs(xq , yq, 40, tol &
               , tneigh%xloc_in(1,i), tneigh%xloc_in(2,i))
          ! (r, s) in the neighbor element stored in this element
          if (tneigh%elnum .ne. -1)  then
             call wspace%elems(tneigh%elnum)%xy2rs(xq , yq &
                  , 40, tol, tneigh%xloc_out(1,i), tneigh%xloc_out(2,i))
          end if

          ! store the full span basis functions of this element 
          ! evaluated at the current local coordinates of this element
          call elem%tbasis%eval(tneigh%xloc_in(1, i), tneigh%xloc_in(2, i) &
               , 0, tneigh%psi_in(:, i))

       end do ! next quadr. point per the current shared segment

    end do ! next neighbor per the current edge

    ! done here
  end subroutine init_elem_edg_quadratures

  ! computes the on edge Rankine–Hugoniot value approximated 
  ! by either a Rieman solver or flux splitting algorithm
  ! and store it in tneigh%Fstar
  !
  subroutine comp_Fstar(wspace, elem, tedg, tneigh)
    implicit none
    class(dg_wspace), target :: wspace ! work space
    class(element_dg2d), target :: elem
    type(edg_dg) :: tedg
    type(neigh_dg) :: tneigh

    ! local vars
    integer :: k
    real*8 :: r, s, nx, ny
    real*8, dimension(elem%neqs) :: UL, UR, FP, FM
    logical, dimension(4) :: f_select
    real*8, dimension(elem%neqs) :: fvl_p, fvl_m
    real*8, dimension(elem%neqs, elem%neqs) :: d_fvl_p, d_fvl_m
    real*8, dimension(elem%neqs, 2) :: Wk, FvL, FvR
    logical :: edge_adiabatic
    class(element_dg2d), pointer :: telem => null()

    ! init
    edge_adiabatic = .false.

    ! loop over Gauss points per this neighboring segment
    do k = 1, tneigh%ngpseg

       ! HARD reset
       UL = 0.0d0; UR = 0.0d0; FP = 0.0d0; FM = 0.0d0

       ! find outgoing unit normals at that gauss points
       nx = tneigh%n(1, k); ny = tneigh%n(2, k)

       !
       !          <<< compute UL and UR procedure >>>
       ! 
       ! evaluate UL at the current Gauss point using interior (r, s) 
       ! coordinates and the basis function of this element
       r = tneigh%xloc_in(1, k); s = tneigh%xloc_in(2, k) 
       call elem%comp_u(r, s, UL)

       ! now decide on UR ... 
! print *, 'tedg%tag = ', tedg%tag
! print *, 'elem%number = ', elem%number
! print *, 'wspace%bcs(tedg%tag)%name = ', wspace%bcs(tedg%tag)%name
       select case (wspace%bcs(tedg%tag)%name)

       case ('interior') ! then UR is in the other neighboring element

          ! evaluate UR at the current Gauss point using 
          ! the neighboring element's local coordinates (r, s)
          ! and the basis functions of the neighboring element
          r = tneigh%xloc_out(1, k); s = tneigh%xloc_out(2, k)  
          call wspace%elems(tneigh%elnum)%comp_u(r, s, UR)
! print *, '6664441221360', UR - UL

       case ('inflow', 'outflow') ! boundary edge; then UR is preset in bcs value

          UR = wspace%bcs(tedg%tag)%val
! print *, 'UR = ', UR
! stop
UR = UL
       case ('wall')

          ! ! impose adiabatic condition
          ! if ( elem%local_props%adia ) then
          !    edge_adiabatic = .true.
          ! end if
UR = UL
          ! do nothing for now!

       case default

          print *, 'could not decide on UR in comp_Fstar! stop'
          stop

       end select

       !        <<< compute Flux procedure >>>
       ! 
       !      special situation : Wall treatment
       if (wspace%bcs(tedg%tag)%name .eq. 'wall') then

          call calc_wall_flux(UL, elem%neqs, elem%gamma, nx, ny, fvl_p, d_fvl_p)

          ! store F+
          FP = fvl_p
          FM = 0.0d0

       else

          ! compute F+ using UL
          f_select = (/ .true. , .false., .false., .false. /)
          call calc_van_leer(UL, elem%neqs, elem%gamma, nx, ny, f_select &
               , fvl_p, fvl_m, d_fvl_p, d_fvl_m)

          ! store F+
          FP = fvl_p

          ! compute F- using UR
          f_select = (/ .false. , .true., .false., .false. /)
          call calc_van_leer(UR, elem%neqs, elem%gamma, nx, ny, f_select &
               , fvl_p, fvl_m, d_fvl_p, d_fvl_m)

          ! store F-
          FM = fvl_m

       end if


       ! store total split flux as the final Rankine–Hugoniot value at Fstar
       tneigh%Fstar(:, k) = FP + FM

       ! call calc_pure_euler2d_flux(Q = UL, gamma = elem%gamma &
       !      , F = FP, G = FM)

! if (all(elem%edgs(:)%tag .eq. 0)) then
! print *, '226644', (FP*nx + FM*ny - tneigh%Fstar(:, k))
! end if

! print *, 'FP = ', FP
! print *, 'FM = ', FM

       if ( elem%local_props%is_viscous ) then

          ! calc left viscous flux FvL
          !
          telem => elem
          ! evaluate the gradient of conservative vars
          r = tneigh%xloc_in(1, k); s = tneigh%xloc_in(2, k)
          call telem%comp_u(r, s, UL) 
          call telem%comp_Wij_point_dg(r, s, Wk)
          ! now use Wk and Hs= ubar to compute viscous fluxes at edges 
          call calc_Fv(telem%local_props%mu, telem%local_props%lambda &
               , telem%gamma, telem%local_props%Pr &
               , Wk, UL, FvL, edge_adiabatic)

          ! calc right viscous flux FvR
          !
          if ( wspace%bcs(tedg%tag)%name .eq. 'interior') then
             telem => wspace%elems(tneigh%elnum)
             ! evaluate the gradient of conservative vars
             r = tneigh%xloc_out(1, k); s = tneigh%xloc_out(2, k)
             call telem%comp_u(r, s, UR)  
             call telem%comp_Wij_point_dg(r, s, Wk)
             ! now use Wk and Hs= ubar to compute viscous fluxes at edges 
             call calc_Fv(telem%local_props%mu, telem%local_props%lambda &
                  , telem%gamma, telem%local_props%Pr &
                  , Wk, UR, FvR, edge_adiabatic)

             ! take average
             FvL = 0.5d0 * (FvL + FvR)
          end if

          nx = tneigh%n(1, k); ny = tneigh%n(2, k)
          tneigh%Fstar(:, k) = tneigh%Fstar(:, k) - FvL(:, 1) * nx - FvL(:, 2) * ny 

       end if

    end do ! next gauss point per neighboring element (on shared segment)

    ! done here
  end subroutine comp_Fstar

  ! computes the right hand side for an element
  ! and stores in rhs(1:neqs, 1:npe).
  ! HINT : needs the element to be initialized with
  !        elem%U before this.
  subroutine comp_elem_rhs(wspace, elem, rhs)
    implicit none
    class(dg_wspace), intent(inout) :: wspace
    class(element_dg2d), intent(inout), target :: elem
    real*8, dimension(:, :), intent(out) :: rhs

    ! local vars
    integer :: iedg, ineigh
    real*8, dimension(size(rhs, 1), size(rhs, 2)) :: tmp
    type(edg_dg), pointer :: tedg => null()
    type(neigh_dg), pointer :: tneigh => null()

    ! HARD reset
    rhs = 0.0d0

    ! compute interior fluxes
    call elem%comp_flux_interior()

print *, 'gulgulFk', elem%Fk
    ! accumulate interior integral to rhs
    ! comp_inter_integral(elem, integ)
    call elem%comp_int_integ(tmp)
    rhs = rhs + tmp
print *, 'gulgulFktmp', rhs

    ! compute boundary integral
    tmp = 0.0d0
    do iedg = 1, elem%nedgs ! loop over edges

       tedg => elem%edgs(iedg) !pick this edge
! print *, 'iedg = ', iedg

       do ineigh = 1, size(tedg%neighs) ! loop over neigh segments on that edge

          tneigh => tedg%neighs(ineigh) ! pick this neighboring segment

          ! first compute the edge flux per that neigh segment
          call wspace%comp_Fstar(elem, tedg, tneigh)

          ! then compute boundary integral over that little segment
          ! and accumulate to tmp
          call elem%comp_bnd_integ(tneigh, tmp) 
! print *, 'tmp = ', tmp
! print *, 'tneigh%Fstar = ', tneigh%Fstar
       end do !segments per that edge

    end do ! edges per that element

    ! finally merge them together to form final rhs
    rhs = rhs - tmp

    ! if (all(elem%edgs(:)%tag .eq. 0) ) then
    !    print *, '888elem # = ' , elem%number, 'rhs_tot = ', rhs
    ! end if

    ! done here
  end subroutine comp_elem_rhs

  ! initializes all elements in the domain
  ! with the given primitive variables
  subroutine init_field(wspace, rho, u, v, P)
    implicit none
    class(dg_wspace), intent(inout), target :: wspace
    real*8, intent(in) :: rho, u, v, P

    ! local vars
    integer :: i
    type(element_dg2d), pointer :: elem => null() 

    ! first init wall nodes which they are important in
    ! viscous boundary conditions 
    call init_wall_nodes(wspace = wspace, tol = 1.d-10)

    do i = 1, size(wspace%elems)

       elem => wspace%elems(i)
       call elem%init_elem_U(rho, u, v, P)

       ! ! apply non-slip viscous wall boundary condition if needed
       ! if ( elem%local_props%is_viscous ) then
       !    call elem%apply_non_slip_wall_bc(elem%U)

       !    ! compute imposed mass matrix too
       !    call elem%create_imposed_mass_matrix()

       ! end if

       call elem%init_mms()

! if ( (elem%number .eq. 3) .or. (elem%number .eq. 4) ) elem%U(3, :) = elem%U(3, :) + 1.4d0 
    end do


    ! done here
  end subroutine init_field

  ! marchs the field by "itrs" number of steps
  ! with fixed time step "dt"
  subroutine march_field(wspace, dt, itrs)
    implicit none
    class(dg_wspace), intent(inout) :: wspace
    real*8, intent(in) :: dt
    integer, intent(in) :: itrs

    ! local vars
    integer :: itr, i

    do itr = 1, itrs ! time step loop

       ! compute the gradient of conservatives vars required
       ! for viscous terms if elements are selected to be viscous
       if ( all(wspace%elems(:)%local_props%is_viscous) ) then
          print *, 'compute grad of U at iteration ', itr
          call comp_Wij(wspace)

          ! ! echo
          ! do i = 1, size(wspace%elems) 
          !    print *, 'gulgul element # ', i, 'Wij = ', minval(wspace%elems(i)%Wij)
          ! end do
          ! stop
       end if

       ! loop over all elements and find rhs
       do i = 1, size(wspace%elems) 
          ! compute rhs
          call wspace%comp_elem_rhs(wspace%elems(i), wspace%elems(i)%rhs)
          print*, 'gulrhs element#', i,  'RHS = ', wspace%elems(i)%rhs(4, :)

       end do
       ! stop

       do i = 1, size(wspace%elems) 
          ! write rhs to U for visualization
          wspace%elems(i)%U = wspace%elems(i)%rhs 
       end do

       return

       ! loop over all elements and update elem%U using elem%rhs 
       do i = 1, size(wspace%elems) 
          ! update elem%U
! if (all(wspace%elems(i)%edgs(:)%tag .eq. 0) ) then
          call wspace%elems(i)%update_explicit_euler(wspace%elems(i)%rhs, dt)
! end if 
       end do

print *, 'itr = ', itr
    end do

    ! done here
  end subroutine march_field

  ! roughly converts DG solution to a PG 
  ! representation for plotting purpuse
  ! using already developed high-order
  ! visulization subroutines for PG
  !
  ! NOTE : This must be replaced with
  !        visualization subroutines
  !        designed for DG solutions
  !
  subroutine udg2upg(wspace, upg)
    implicit none
    class(dg_wspace), intent(in), target :: wspace
    real*8, dimension(:, :), intent(out) :: upg

    ! local vars
    integer :: i, j, tpt
    class(element_dg2d), pointer :: telem => null()

    do i = 1, size(wspace%elems)

       telem => wspace%elems(i)
       do j = 1, telem%npe
          tpt = wspace%grd%icon(i, j)
          upg(:, tpt) = telem%U(:, j)
       end do

    end do

    ! done here
  end subroutine udg2upg

  ! scatter plot very handy for visualization of
  ! curves boundaries using additional scattered 
  ! points interpolated on boundary parametrization.
  !
  ! this is only for DG grids  
  !
  subroutine scatter_tecplot_dg(wspace, outfile, append_flag, title)
    implicit none
    class(dg_wspace), intent(in), target :: wspace
    character(len = *), intent(in) :: outfile, title
    logical, intent(in) :: append_flag

    ! local vars
    integer i, j, k, neqs, nnodesg
    real*8, dimension(:, :), pointer :: u => null()
    real*8, dimension(:), pointer :: x => null(), y => null()

    ! init
    neqs = wspace%elems(1)%neqs
    ! count total nodes in all elements
    nnodesg = 0
    do i = 1, size(wspace%elems)
       nnodesg = nnodesg + wspace%elems(i)%npe
    end do

    ! opening the output file
    if (append_flag) then
       ! append
       open(unit = 10, file = outfile, position = 'append')
    else !or just open a new file
       open(unit = 10, file = outfile)
    end if

    ! write header
    write(10, *) 'title = "',title,'"'

    write(10, '(A)', advance = 'no') 'variables = "x", "y"'
    do i = 1, neqs
       write(10, '(A, I1, A)', advance = 'no') ', "u', i,'"'
    end do

    write(10,*)
    write(10, '(A, I7, A)', advance = 'no') 'zone I = ' &
         , nnodesg, ', datapacking = point'

    write(10,*)
    ! write coordinates and values of the vector field [u]
    do k = 1, size(wspace%elems)
       x => wspace%elems(k)%x(1, :)
       y => wspace%elems(k)%x(2, :)
       u => wspace%elems(k)%U
       do i = 1, size(x)
          write(10, '(F30.17, A, F30.17)', advance = 'no') x(i), ' ',  y(i)
          do j = 1, neqs
             write(10, '(A, F30.17)', advance = 'no') ' ',  u(j,i)
          end do
          write(10,*) 
       end do
    end do

    ! shut down the output
    close(10)

    ! done here
  end subroutine scatter_tecplot_dg

  ! computes the Jacobian of the split fluxes and store
  ! them in tneigh%dFpm so later we can readily evaluate
  ! the boundary integrals.
  !
  subroutine comp_dFpm(wspace, elem, tedg, tneigh)
    implicit none
    class(dg_wspace) :: wspace ! work space
    class(element_dg2d) :: elem
    type(edg_dg) :: tedg
    type(neigh_dg) :: tneigh

    ! local vars
    integer :: k
    real*8 :: r, s, nx, ny
    real*8, dimension(elem%neqs) :: UL, UR
    logical, dimension(4) :: f_select
    real*8, dimension(elem%neqs) :: fvl_p, fvl_m
    real*8, dimension(elem%neqs, elem%neqs) :: d_fvl_p, d_fvl_m


    ! loop over Gauss points per this neighboring segment
    do k = 1, tneigh%ngpseg

       ! HARD reset
       UL = 0.0d0; UR = 0.0d0

       ! find outgoing unit normals at that gauss points
       nx = tneigh%n(1, k); ny = tneigh%n(2, k)

       !
       !          <<< compute UL and UR procedure >>>
       ! 
       ! evaluate UL at the current Gauss point using interior (r, s) 
       ! coordinates and the basis function of this element
       r = tneigh%xloc_in(1, k); s = tneigh%xloc_in(2, k) 
       call elem%comp_u(r, s, UL)

       ! now decide on UR ... 
       select case (wspace%bcs(tedg%tag)%name)

       case ('interior') ! then UR is in the other neighboring element

          ! evaluate UR at the current Gauss point using 
          ! the neighboring element's local coordinates (r, s)
          ! and the basis functions of the neighboring element
          r = tneigh%xloc_out(1, k); s = tneigh%xloc_out(2, k)  
          call wspace%elems(tneigh%elnum)%comp_u(r, s, UR)

       case ('inflow', 'outflow') ! boundary edge; then UR is preset in bcs value

          UR = wspace%bcs(tedg%tag)%val
          
       case ('wall')

          ! do nothing for now!

       case default

          print *, 'could not decide on UR in comp_dFpm! stop'
          stop

       end select

       !        <<< compute Flux Jacobian procedure >>>
       ! 
       !      special situation : Wall treatment
       if (wspace%bcs(tedg%tag)%name .eq. 'wall') then

          call calc_wall_flux(UL, elem%neqs, elem%gamma, nx, ny, fvl_p, d_fvl_p)

          ! store wall flux Jacobian
          tneigh%dFpm(:, :, 1, k) = d_fvl_p
          tneigh%dFpm(:, :, 2, k) = 0.0d0

       else

          ! compute dF+/du using UL
          f_select = (/ .false. , .false., .true., .false. /)
          call calc_van_leer(UL, elem%neqs, elem%gamma, nx, ny, f_select &
               , fvl_p, fvl_m, d_fvl_p, d_fvl_m)

          ! store it
          tneigh%dFpm(:, :, 1, k) = d_fvl_p

          ! compute dF-/du using UR
          f_select = (/ .false. , .false., .false., .true. /)
          call calc_van_leer(UR, elem%neqs, elem%gamma, nx, ny, f_select &
               , fvl_p, fvl_m, d_fvl_p, d_fvl_m)

          ! store it
          tneigh%dFpm(:, :, 2, k) = d_fvl_m

       end if

    end do ! next gauss point per neighboring element (on shared segment)

    ! done here
  end subroutine comp_dFpm

  ! computes the boundary integral of split Jacobians
  subroutine comp_jac_bnd_integral(wspace, elem, tedg, tneigh, integ)
    implicit none
    class(dg_wspace) :: wspace ! work space
    class(element_dg2d) :: elem
    type(edg_dg) :: tedg
    type(neigh_dg) :: tneigh
    real*8, dimension(:,:), intent(inout) :: integ 

    ! local vars
    integer :: i, k
    real*8 :: r, s
    real*8, dimension(elem%neqs) :: UL, UR, dFXu

    ! loop over Gauss points per this neighboring segment
    do k = 1, tneigh%ngpseg

       ! HARD reset
       UL = 0.0d0; UR = 0.0d0

       !
       !          <<< compute UL and UR procedure >>>
       ! 
       ! evaluate UL at the current Gauss point using interior (r, s) 
       ! coordinates and the basis function of this element
       r = tneigh%xloc_in(1, k); s = tneigh%xloc_in(2, k) 
       call elem%comp_u(r, s, UL)

       ! now decide on UR ... 
       select case (wspace%bcs(tedg%tag)%name)

       case ('interior') ! then UR is in the other neighboring element

          ! evaluate UR at the current Gauss point using 
          ! the neighboring element's local coordinates (r, s)
          ! and the basis functions of the neighboring element
          r = tneigh%xloc_out(1, k); s = tneigh%xloc_out(2, k)  
          call wspace%elems(tneigh%elnum)%comp_u(r, s, UR)

       case ('inflow', 'outflow') ! boundary edge; then UR is preset in bcs value

          UR = 0.0d0 !wspace%bcs(tedg%tag)%val
       UL = 0.0d0; UR = 0.0d0          
       case ('wall')

          ! do nothing for now!
       UL = 0.0d0; UR = 0.0d0
       case default

          print *, 'could not decide on UR in comp_jac_bnd_integral! stop'
          stop

       end select

       ! compute the term under the integral
       dFXu = matmul(tneigh%dFpm(:, :, 1, k), UL) &
            + matmul(tneigh%dFpm(:, :, 2, k), UR)

       ! add to boundary integral
       do i = 1, elem%npe
          integ(:, i) = integ(:, i) &
               + tneigh%psi_in(i, k) * dFXu  &
               * tneigh%s(k) * tneigh%W(k)
       end do


    end do ! next gauss point per neighboring element (on shared segment)

    ! done here
  end subroutine comp_jac_bnd_integral

  ! computes matrix-free matrix-vector product
  ! of the diagonal and off diagonal Jacobian times
  ! the current solution in the Krylov subspace. 
  subroutine comp_Ax(wspace, elem, integ)
    implicit none
    class(dg_wspace) :: wspace ! work space
    class(element_dg2d), target :: elem
    real*8, dimension(:, :), intent(out) :: integ

    ! locals vars
    integer :: iedg, ineigh
    type(edg_dg), pointer :: tedg => null()
    type(neigh_dg), pointer :: tneigh => null()
    real*8, dimension((elem%neqs * elem%npe), 1) :: MU

    integ = 0.0d0

    ! add contribution of the main diagonal (element itself)
    !
    !        ------------  DIAGONAL -------------------
    !
    ! interior jacobian integral
    call elem%comp_int_jac_integ(integ)
    ! contribution due to 1/dt * Mass matrix
    MU = matmul(elem%Mass, reshape(elem%U, (/ (elem%neqs * elem%npe), 1 /)))
    integ =  (1.0d0 / wspace%dt) * reshape(MU, (/elem%neqs, elem%npe/) ) - integ
    !
    !
    ! add contribution of off diagonal terms (edge integrals)
    !
    !        ----------  OFF DIAGONAL ------------------
    !
    !
    do iedg = 1, elem%nedgs ! loop over edges

       tedg => elem%edgs(iedg) !pick this edge

       do ineigh = 1, size(tedg%neighs) ! loop over neigh segments on that edge

          tneigh => tedg%neighs(ineigh) ! pick this neighboring segment

          ! accumulate the split Jacobian integral 
          call wspace%comp_jac_bnd_integral(elem, tedg, tneigh, integ)

       end do !segments per that edge

    end do ! edges per that element

    ! done here
  end subroutine comp_Ax

  ! computes full matrix-vector product
  ! of the diagonal and off diagonal Jacobian times
  ! the current solution in the Krylov subspace. 
  ! where A*x is defined as:
  !
  ! A*x = (-dF_newton/du^(n+1,s)) * elem%U
  !
  !
  subroutine comp_full_Ax(wspace)
    implicit none
    class(dg_wspace), target :: wspace ! work space

    ! locals vars
    integer :: ielem
    class(element_dg2d), pointer :: elem => null()

    do ielem = 1, size(wspace%elems)
       elem => wspace%elems(ielem)

       ! compute a complete row (diagonal+off diagonal) of
       ! matrix-vector product corresponding to an element
       ! on that row and store the result in elem%Ax.
       call comp_Ax(wspace, elem, elem%Ax)

       ! change the sign
       elem%Ax = -1.0d0 * elem%Ax
 
    end do

    ! done here
  end subroutine comp_full_Ax

  ! sums the elemental number of unknowns
  ! for all elements in the wspace
  !
  ! NOTE : includes fixed and free nodes together
  !
  subroutine comp_total_num_unknowns(wspace, ntot)
    implicit none
    class(dg_wspace), intent(in) :: wspace
    integer, intent(out) :: ntot

    ! local vars
    integer :: i

    ! HARD reset
    ntot = 0

    do i = 1, size(wspace%elems)

       ! add to total number of nodes
       ntot = ntot + size(wspace%elems(i)%U)

    end do

    ! done here
  end subroutine comp_total_num_unknowns

  ! put all elements%Ax in one contigous
  ! 1d array "contig" and pops it outttt!
  subroutine pop_Ax(wspace, contig)
    implicit none
    class(dg_wspace), intent(in), target :: wspace
    real*8, dimension(:), intent(out) :: contig

    ! local vars
    integer :: ielem, ieq, inpe, jj
    class(element_dg2d), pointer :: elem => null()

    jj = 1
    do ielem = 1, size(wspace%elems)
       elem => wspace%elems(ielem)
       do inpe = 1, elem%npe
          do ieq = 1, elem%neqs
             contig(jj) = elem%Ax(ieq, inpe)
             jj = jj + 1
          end do
       end do
    end do

    ! done here
  end subroutine pop_Ax

  function norm(x)
    implicit none
    real*8, dimension(:), intent(in) :: x
    real*8 :: norm

    ! locals
    integer :: i

    norm = 0.0d0 !reset
    do i = 1, size(x)
       norm = norm + x(i)* x(i)
    end do

    norm = sqrt(norm)

    ! done
  end function norm

  ! puts values of the contigous 1d array "contig"
  ! into elements%U while preserving the order 
  ! of the data
  subroutine push_U(wspace, contig)
    implicit none
    class(dg_wspace), target :: wspace
    real*8, dimension(:), intent(in) :: contig

    ! local vars
    integer :: ielem, ieq, inpe, jj
    class(element_dg2d), pointer :: elem => null()

    jj = 1
    do ielem = 1, size(wspace%elems)
       elem => wspace%elems(ielem)
       do inpe = 1, elem%npe
          do ieq = 1, elem%neqs
             elem%U(ieq, inpe) = contig(jj)
             jj = jj + 1
          end do
       end do
    end do

    ! done here
  end subroutine push_U

  subroutine gmres_raw(wspace, epsil, num, nrst, res)
    implicit none
    ! inputs
    class(dg_wspace), target :: wspace ! work space
    real*8, intent(in) :: epsil
    integer, intent(in) :: num, nrst

    !outputs
    real*8, dimension(:), allocatable :: res

    ! locals
    integer :: k, j, i
    integer :: finish_flag, jrst
    real*8, dimension(:,:), allocatable :: v, h
    real*8, dimension(:), allocatable :: r0, g, yk, uk, s, c, zk, alpha, tmp
    real*8 :: norm_tmp, delta, gamma, resi
    integer :: nrows, nvars, len_res, ielem
    class(element_dg2d), pointer :: elem => null()

    ! resetting the flag
    finish_flag = 0
    resi = 0.0d0
    call comp_total_num_unknowns(wspace, nrows)
    !nrows = size(b)
    nvars = nrows
    !nvars = size(mf%free)
    len_res = 0
    if ( allocated(res) ) deallocate(res)

    ! we don't use dynamic mem allocation in the gmres inner-outer loops 
    ! to increase spped. this might not be memory efficient though.
    allocate(v(nvars,num+1), r0(nrows), g(num+1), s(num), c(num))
    allocate(yk(nrows), uk(nvars), zk(nvars) )
    allocate(h(num+1,num), alpha(num))

    ! do j number of restarts
    do jrst = 1, nrst    
       ! reset ------------
       v = 0.0d0; h = 0.0d0
       c = 0.0d0; s = 0.0d0
       zk = 0.0d0; g = 0.0d0
       alpha = 0.0d0
       yk = 0.0d0
       ! ------------------
       ! ***********************
       ! sol = elem%U
       ! r0 = elem%Ax
       call wspace%comp_full_Ax()
       ! call mf%Ax(sol, r0)
       ! b = elem%rhs
       do ielem = 1, size(wspace%elems)
          wspace%elems(ielem)%Ax = wspace%elems(ielem)%rhs - wspace%elems(ielem)%Ax 
       end do
       !r0 = b - r0
       ! r0 = b - matmul(A,sol)
       ! ***********************
       call pop_Ax(wspace, v(:,1))
       !v(:,1) = r0(mf%free)
       norm_tmp = norm(v(:,1))
       v(:,1) = v(:,1)/norm_tmp
       g(1) = norm_tmp
       ! Arnoldi iterations    
       do k = 1, num
          g(k+1) = 0.0d0
          ! ***********************
          call push_U(wspace, v(:,k))
          ! yk = 0.0d0
          call solve_prec_elemental(wspace)
          ! yk(mf%free) = v(:,k) / N
          call wspace%comp_full_Ax()
          call pop_Ax(wspace, uk)

          ! call mf%Ax(yk, r0)
          ! uk = r0(mf%free)
          ! yk = v(:,k) / N
          ! uk = matmul(A, yk)
          ! ***********************
          do j = 1, k
             h(j,k) = dot_product(v(:,j),uk)
             !grahm-schmidt orthogonalization
             uk = uk - h(j,k) * v(:,j) 
          end do
          h(k+1,k) = norm(uk)        
          v(:,k+1) = uk/h(k+1,k)
          ! applying givens
          do j = 1, (k-1)
             delta = h(j,k)
             h(j,k) = c(j)*delta + s(j)*h(j+1,k)
             h(j+1,k) = -s(j)*delta + c(j) * h(j+1,k)
          end do
          gamma = sqrt(h(k,k)**(2.0d0) + h(k+1,k)**(2.0d0))
          c(k) = h(k,k) / gamma
          s(k) = h(k+1,k) / gamma
          h(k,k) = gamma
          h(k+1,k) = 0.0d0
          delta = g(k)
          g(k) = c(k) * delta + s(k) * g(k+1)
          g(k+1) = -s(k) * delta + c(k) * g(k+1)
          resi = abs(g(k+1))
          ! add to residual history
          len_res = len_res + 1
          allocate(tmp(len_res))
          if(len_res > 1) tmp(1:(len_res-1)) = res(1:(len_res-1)) 
          tmp(len_res) = resi
          call move_alloc(tmp, res)
          !
          if( resi <= epsil) then
             finish_flag = 1
             goto 100
          end if

       end do
       k = k - 1
       ! solving backward for alpha
100    alpha(k) = g(k) / h(k,k)
       do i = k-1, 1, -1
          alpha(i) = g(i) / h(i,i)
          do j = i+1, k      
             alpha(i) = alpha(i) - h(i,j)/h(i,i) * alpha(j)     
          end do
       end do
       ! compute the final directional vector using
       ! the combination of search vectors 
       do j = 1, k
          zk = zk + alpha(j)*v(:,j)
       end do
       ! updating solution
       ! ***********************
       !zk = zk/N
       !r0 = 0.0d0
       ! push zk=deltaU = delta_sol into elems%U 
       call push_U(wspace, zk)
       call solve_prec_elemental(wspace)
       ! r0(mf%free) = zk
       ! add U0 to U=DU to get new U
       ! do ielem = 1, size(wspace%elems)
       !    wspace%elems(ielem)%U = wspace%elems(ielem)%Us + wspace%elems(ielem)%U  
       ! end do

       do ielem = 1, size(wspace%elems)
          elem => wspace%elems(ielem)
          elem%Us = elem%Us + elem%U
       end do

       !sol = sol + r0     
       ! sol = sol + zk/N    
       ! ***********************
       if(finish_flag == 1) then 
          goto 200        
       end if


    end do
    ! subroutine final clean-ups
200 deallocate(v, r0, g, s, c)
    deallocate(yk, uk, zk, h, alpha)

    return !done <here>

  end subroutine gmres_raw

  ! computes Fs in newton iterations
  ! and store it in elems%rhs
  subroutine comp_newton_Fs(wspace)
    implicit none
    class(dg_wspace), target :: wspace

    ! local vars
    integer :: ielem
    class(element_dg2d), pointer :: elem => null()
    real*8, dimension(:, :), allocatable :: MdU

    do ielem = 1, size(wspace%elems)

       elem => wspace%elems(ielem)
       call wspace%comp_elem_rhs(elem, elem%rhs)

       if(allocated(MdU)) deallocate(MdU)
       allocate(MdU((elem%neqs * elem%npe), 1))

       MdU = matmul(     elem%Mass &
            , reshape( (elem%U - elem%U0) , (/ (elem%neqs * elem%npe), 1 /) )  )

       elem%rhs =  (1.0d0 / wspace%dt) * reshape(MdU, (/elem%neqs, elem%npe/) ) - elem%rhs

    end do

    ! little clean-up
    if(allocated(MdU)) deallocate(MdU)

    ! done here
  end subroutine comp_newton_Fs

  ! refreshes the matrix free storage
  ! that holds the entire Jacobian matrix
  ! and reinitializes it with the current elems%U 
  subroutine refresh_jacobian(wspace)
    implicit none
    class(dg_wspace), target :: wspace

    ! local vars
    integer :: ielem, iedg, ineigh
    class(element_dg2d), pointer :: elem => null()
    type(edg_dg), pointer :: tedg => null()
    type(neigh_dg), pointer :: tneigh => null()


    do ielem = 1, size(wspace%elems)

       elem => wspace%elems(ielem)

       call elem%comp_pure_flux_jac()

       do iedg = 1, elem%nedgs ! loop over edges

          tedg => elem%edgs(iedg) !pick this edge

          do ineigh = 1, size(tedg%neighs) ! loop over neigh segments on that edge

             tneigh => tedg%neighs(ineigh) ! pick this neighboring segment

             call wspace%comp_dFpm(elem, tedg, tneigh)

          end do !segments per that edge

       end do ! edges per that element

    end do

    ! done here
  end subroutine refresh_jacobian

  ! marches the workspace using 1st-order euler implicit
  ! method and update elems%U in place
  subroutine march_euler_implicit(wspace, dt, itrs, inewtons, num, nrst, epsil)
    implicit none
    class(dg_wspace), target :: wspace
    real*8, intent(in) :: dt, epsil
    integer, intent(in) :: itrs, inewtons, num, nrst

    ! local vars
    integer :: itr, ielem, inewton
    class(element_dg2d), pointer :: elem => null()
    real*8, dimension(:), allocatable :: res

    ! resets
    wspace%dt = dt

    do itr = 1, itrs

       ! take a copy of current solution as
       ! the initial solution
       do ielem = 1, size(wspace%elems)
          elem => wspace%elems(ielem)
          elem%U0 = elem%U
       end do

       ! begin newton iterations
       do inewton = 1, inewtons

          ! refresh Jacobian matrix
          call refresh_jacobian(wspace)

          ! compute RHS
          call comp_newton_Fs(wspace)

          ! take a copy of current solution as
          ! the current Us_n+1 solution before
          ! Krylov subspace-based update 
          do ielem = 1, size(wspace%elems)
             elem => wspace%elems(ielem)
             elem%Us = elem%U
             elem%U = 0.0d0
          end do

          ! solve for U in place
          call comp_prec_elemental(wspace)
          call gmres_raw(wspace, epsil, num, nrst, res)

          do ielem = 1, size(wspace%elems)
             elem => wspace%elems(ielem)
             elem%U = elem%Us
          end do

          ! show progress
          print *, 'itr = ', itr, 'inewton = ', inewton &
               , 'res_max = ', maxval(abs(res)), 'res_min = ', minval(abs(res))

       end do ! end of newton itrs


    end do ! end of physical time step

    ! done here
  end subroutine march_euler_implicit

  ! computes the elemental preconditioner
  ! based on initially provided elem%Mass and addes A*Xu due to
  ! interior integral of pure flux jacobian
  ! and boundary integral of dF+ and write
  ! it to elem%LUmass_imp. Then overwrite elem%LUmass_imp
  ! with its LU decomposition in place using
  ! LAPACK routines. 
  subroutine comp_prec_elemental(wspace)
    implicit none
    class(dg_wspace), target :: wspace

    ! local vars
    integer :: ielem, k, m, l, q, i, j, ii, qq
    integer :: iedg, ineigh
    real*8 :: tmp
    class(element_dg2d), pointer :: elem => null()
    class(edg_dg), pointer :: tedg => null()
    class(neigh_dg), pointer :: tneigh => null()

    ! LAPACK LU temp. vars
    integer :: LDA, INFO

    do ielem = 1, size(wspace%elems)
       elem => wspace%elems(ielem)

       ! 
       !        NOTE
       !
       ! assuming elem%Mass is already allocated and initialized
       !
       elem%LUmass_imp = 1.0d0 / wspace%dt * elem%Mass

       ! then add contributions from interior and boundary here ...
       ! 1- contribution from interior integral of pure Jacobian matrix
       do k = 1, elem%ngauss
          do m = 1, elem%npe
             do l = 1, elem%npe
                do q = 1, elem%neqs
                   do i = 1, elem%neqs

                      ! comp temp val
                      tmp = 0.0d0

                      do j = 1, 2 !ndim
                         tmp = tmp - elem%d_psi_d_x(l, k, j) * elem%dFk(i,q,k, j) &
                              * elem%psi(m, k) &
                              * elem%coeff * elem%JJ(k) * elem%W(k)
                      end do

                      ! add it to correct location to diagonal matrix
                      ii = (l-1) * elem%neqs + i
                      qq = (m-1) * elem%neqs + q
                      elem%LUmass_imp(ii, qq) = elem%LUmass_imp(ii, qq) + tmp

                   end do
                end do
             end do
          end do
       end do

       ! 2- add contributions due boundary integral of split
       !    jacobians
       do iedg = 1, elem%nedgs ! loop over edges

          tedg => elem%edgs(iedg) !pick this edge

          do ineigh = 1, size(tedg%neighs) ! loop over neigh segments on that edge

             tneigh => tedg%neighs(ineigh) ! pick this neighboring segment

             ! loop over Gauss points per this neighboring segment
             do k = 1, tneigh%ngpseg

                ! add to boundary integral matrix
                do m = 1, elem%npe
                   do l = 1, elem%npe
                      do q = 1, elem%neqs
                         do i = 1, elem%neqs
                            ! find loc in the maxx matrix
                            ii = (l-1) * elem%neqs + i
                            qq = (m-1) * elem%neqs + q
                            ! compute the entry
                            elem%LUmass_imp(ii, qq) = elem%LUmass_imp(ii, qq) &
                                 + tneigh%psi_in(l, k) * tneigh%dFpm(i, q, 1, k) &
                                 * tneigh%psi_in(m, k) * tneigh%s(k) * tneigh%W(k)
                         end do
                      end do
                   end do
                end do

             end do ! next gauss point per neighboring element (on shared segment)

          end do !segments per that edge

       end do ! edges per that element





       ! finally
       ! compute and store LU of the mass matrix of implicit formulation
       LDA = size(elem%Mass, 1)
       call dgetrf(LDA, LDA, elem%LUmass_imp, LDA, elem%IPIVmass_imp, INFO)
       if ( INFO .ne. 0) then
          print *, 'something is wrong in LU factorization of ' &
               , ' mass matrix of implicit formulation! stop'
          stop
       end if

    end do ! next element in the workspace

    ! done here
  end subroutine comp_prec_elemental

  ! applys diagonal preconditioer by
  ! solving the following
  ! 
  ! elem%LUmass_imp * x = elem%U
  !
  ! and then puts "x" in elem%U
  ! this is done for all elements
  ! in the given workspace
  !
  subroutine solve_prec_elemental(wspace)
    implicit none
    class(dg_wspace), target :: wspace

    ! local vars
    integer :: ielem
    class(element_dg2d), pointer :: elem => null()

    ! LAPACK LU temp. vars
    integer :: N, INFO
    real*8, dimension(:, :), allocatable :: rhs_lapack

    do ielem = 1, size(wspace%elems)
       elem => wspace%elems(ielem)

       ! init
       N = elem%neqs * elem%npe
       allocate(rhs_lapack(N, 1))

       ! solve using already stored LU
       rhs_lapack = reshape(elem%U, (/ N, 1 /) )

       CALL DGETRS( 'No transpose', N, 1, elem%LUmass_imp &
            , N, elem%IPIVmass_imp, rhs_lapack, N, INFO )

       if ( INFO .ne. 0) then
          print *, 'something is wrong in LU solve in' &
               , ' elemental preconditioner! stop'
          stop
       end if

       ! update
       elem%U = reshape( rhs_lapack, (/ elem%neqs, elem%npe /))

       ! clean ups
       deallocate(rhs_lapack)

    end do ! next element in the workspace


    ! done here
  end subroutine solve_prec_elemental

  ! computes the on-edge U-average 1/2(u- + u+)
  ! and stote it in tneigh%Hs
  subroutine comp_Hs(wspace, elem, tedg, tneigh)
    implicit none
    class(dg_wspace) :: wspace ! work space
    class(element_dg2d) :: elem
    type(edg_dg) :: tedg
    type(neigh_dg) :: tneigh

    ! local vars
    integer :: k
    real*8 :: r, s
    real*8, dimension(elem%neqs) :: UL, UR

    ! loop over Gauss points per this neighboring segment
    do k = 1, tneigh%ngpseg

       ! HARD reset
       UL = 0.0d0; UR = 0.0d0

       !
       !          <<< compute UL and UR procedure >>>
       ! 
       ! evaluate UL at the current Gauss point using interior (r, s) 
       ! coordinates and the basis function of this element
       r = tneigh%xloc_in(1, k); s = tneigh%xloc_in(2, k) 
       call elem%comp_u(r, s, UL)

       ! now decide on UR ... 
       select case (wspace%bcs(tedg%tag)%name)

       case ('interior') ! then UR is in the other neighboring element

          ! evaluate UR at the current Gauss point using 
          ! the neighboring element's local coordinates (r, s)
          ! and the basis functions of the neighboring element
          r = tneigh%xloc_out(1, k); s = tneigh%xloc_out(2, k)  
          call wspace%elems(tneigh%elnum)%comp_u(r, s, UR)

       case ('inflow', 'outflow') ! boundary edge; then UR is preset in bcs value

          !UR = wspace%bcs(tedg%tag)%val
UR = UL
       case ('wall') ! NOTE : ONLY Viscous

       ! HARD reset
       ! UL(2:3) = 0.0d0; UR = 0.0d0
       !UL(2:3) = 0.0d0
          UR = UL

       case default

          print *, 'could not decide on UR in comp_Hs! stop'
          stop

       end select

       ! store the avarage u on this neigh at tneigh%Hs
       tneigh%Hs(:, k) = 0.5d0 * (UL + UR) 

    end do ! next gauss point per neighboring element (on shared segment)

    ! done here
  end subroutine comp_Hs

  ! computes the Wij which is the derivative of
  ! conservative vars and store it per element.
  ! later we can use W2w transform to obtain the
  ! derivative of primitive vars.
  !
  subroutine comp_Wij(wspace)
    implicit none
    class(dg_wspace), target :: wspace

    ! local vars
    integer :: ielem, iedg, ineigh
    class(element_dg2d), pointer :: elem => null()
    type(edg_dg), pointer :: tedg => null()
    type(neigh_dg), pointer :: tneigh => null()

    ! LAPACK vars
    integer :: NN, INFO

    do ielem = 1, size(wspace%elems)
       elem => wspace%elems(ielem)

       ! 1- compute -int_Omega(w,j * ui dOmega)
       !
       ! NOTE : also HARD reset elem%Wij back to zero
       !
       call elem%comp_inter_ui_integral(elem%Wij)  

       ! 2- accumulate the boundary integral 
       do iedg = 1, elem%nedgs ! loop over edges

          tedg => elem%edgs(iedg) !pick this edge


          do ineigh = 1, size(tedg%neighs) ! loop over neigh segments on that edge

             tneigh => tedg%neighs(ineigh) ! pick this neighboring segment

             ! first compute the u-average per that neigh segment
             call comp_Hs(wspace, elem, tedg, tneigh)


             ! then compute boundary integral over that little segment
             call elem%comp_bnd_ubar_integral(tneigh, elem%Wij)

          end do !segments per that edge

       end do ! edges per that element

       ! 3- now this element is ready for and LU solve
       ! ----------------  LAPACK SOLVE  ------------------------
       ! init
       NN = elem%neqs * elem%npe

       ! solve using already stored LU of elemental mass matrix

       CALL DGETRS( 'No transpose', NN, 2, elem%LUmass &
            , NN, elem%IPIVmass, elem%Wij, NN, INFO )

       if ( INFO .ne. 0) then
          print *, 'something is wrong in LU solve in computing Wij! stop'
          stop
       end if
       ! --------------------------------------------------------

    end do ! end of elements

    ! done here
  end subroutine comp_Wij


  ! just call it once per grid configuration
  ! 
  subroutine init_wall_nodes(wspace, tol)
    implicit none
    class(dg_wspace), target :: wspace
    real*8, intent(in) :: tol

    ! local vars
    integer :: ielem, iedg, ipt, icurve, pt, ii
    class(element_dg2d), pointer :: elem => null()
    class(edg_dg), pointer :: tedg => null()
    real*8 :: xx, yy, tt

    do ielem = 1, size(wspace%elems)

       elem => wspace%elems(ielem)

       do iedg = 1, size(elem%edgs)

          tedg => elem%edgs(iedg)

          ! first add the on-edge wall nodes if this edge is tagged "Wall"
          if (wspace%bcs(tedg%tag)%name .eq. 'wall') then
             do ii = 1, size(tedg%pts)
                call add_wall_node(elem%wall_nodes, tedg%pts(ii))
             end do
          end if

          ! for each nodes on edge see if it is on any wall???
          do ipt = 1, size(tedg%pts)

             pt = tedg%pts(ipt)
             xx = elem%x(1, pt)
             yy = elem%x(2, pt)

             ! loop over wall boundaries
             ! NOTE : curve number is the tag
             do icurve = 1, size(wspace%grd%bn_curves)

                ! if that boundary curve is not wall then go to
                ! the next boundary curve 
                if ( wspace%bcs(icurve)%name .ne. 'wall' ) cycle

                ! then find that node on this wall
                tt = 0.0d0
                call find_t(wspace%grd, icurve, xx, yy, tol, tt)

                if ( tt .eq. -1.0d0 ) then

                   cycle ! not found on that wall curve

                else ! add it to the wall_nodes of this element

                   call add_wall_node(elem%wall_nodes, pt)

                end if

             end do ! next wall curve

          end do ! next node on this edge

       end do ! next edge of this element

    end do ! next element in the workspace

    ! done here
  end subroutine init_wall_nodes

  ! add an integer "pt" to the end of integer array "idat"
  subroutine add_wall_node(idat, pt)
    implicit none
    integer, dimension(:), allocatable :: idat
    integer, intent(in) :: pt

    ! local vars
    integer :: n, ii
    integer, dimension(:), allocatable :: itmp

    if ( .not. allocated(idat) ) then
       n = 1
    else

       ! check if it already exist then quickly return
       do ii = 1, size(idat)
          if ( idat(ii) .eq. pt ) return
       end do

       ! if not then continue adding it ...
       n = size(idat) + 1

    end if

    allocate(itmp(n))

    if ( n > 1 ) itmp(1:(n-1)) = idat(1:(n-1))

    itmp(n) = pt

    call move_alloc(itmp, idat)

    if (allocated(itmp)) deallocate(itmp)

    ! done here
  end subroutine add_wall_node

end module dg_workspace
