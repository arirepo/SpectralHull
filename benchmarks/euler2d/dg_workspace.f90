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
     type(grid), public :: grd
     type(element_dg2d), dimension(:), allocatable :: elems
     type(bc), dimension(:), allocatable :: bcs ! zero-based

   contains
     procedure, public :: init => init_wspace 
     procedure, public :: init_edg_quadrat => init_elem_edg_quadratures
     procedure, public :: comp_Fstar

  end type dg_wspace


  public :: dg_wspace

contains

  ! initilizes the workspace
  ! HINT : the grid in wspace%grd should already
  ! been allocated and initialized properly
  subroutine init_wspace(wspace, neqs, gamma, bc_names, bc_vals, tol)
    implicit none
    class(dg_wspace), intent(inout) :: wspace
    integer, intent(in) :: neqs
    real*8, intent(in) :: gamma
    character(len = 128), dimension(:), intent(in) :: bc_names
    real*8, dimension(:, :), intent(in) :: bc_vals
    real*8, intent(in) :: tol

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

    ! proceed to initialize the elements one by one
    do i = 1, size(wspace%elems)
       call wspace%elems(i)%init(wspace%elems(i)%element &
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
       r = nint(1.5d0 * dble(max_npedg) / dble(2)) ! 1.5 is safety factor :)
       tneigh%ngpseg = r

       ! allocate neigh specific arrays
       allocate(tneigh%xi(r), tneigh%W(r))
       allocate(tneigh%xloc_in(2, r), tneigh%xloc_out(2, r))
       allocate(tneigh%x(2, r), tneigh%dx(2, r))
       allocate(tneigh%n(2, r), tneigh%s(r))
       allocate(tneigh%psi_in(elem%npe, r))
       ! alloc/init Fstar(1:neqs, 1:ngpseg) 
       allocate(tneigh%Fstar(elem%neqs, r))

       tneigh%xi = 0.0d0; tneigh%W = 0.0d0
       tneigh%xloc_in = 0.0d0; tneigh%xloc_out = 0.0d0
       tneigh%x = 0.0d0; tneigh%dx = 0.0d0
       tneigh%n = 0.0d0; tneigh%s = 0.0d0
       tneigh%psi_in = 0.0d0
       tneigh%Fstar = 0.0d0

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
    class(dg_wspace) :: wspace ! work space
    class(element_dg2d) :: elem
    type(edg_dg) :: tedg
    type(neigh_dg) :: tneigh

    ! local vars
    integer :: k
    real*8 :: r, s, nx, ny
    real*8, dimension(elem%neqs) :: UL, UR, FP, FM
    logical, dimension(4) :: f_select
    real*8, dimension(elem%neqs) :: fvl_p, fvl_m
    real*8, dimension(elem%neqs, elem%neqs) :: d_fvl_p, d_fvl_m

    ! loop over Gauss points per this neighboring segment
    do k = 1, tneigh%ngpseg

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

       case default

          ! do nothing for now!

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

    end do ! next gauss point per neighboring element (on shared segment)

    ! done here
  end subroutine comp_Fstar

end module dg_workspace
