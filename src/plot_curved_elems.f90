module plot_curved_elems
  use grid_opt
  use element_opt
  use lag_basis
  use spem_2d
  use trimesher
  use spline
  implicit none


  private

  public :: vis_curved_grid, adaptive_vis_curved_grid

contains

  ! computes Jacobian of the transformation "jac" and 
  ! "Jstar" or the inverse of the Jacobian
  ! of the transformation at the point "(r,s)" 
  ! within the element with number "ielem" in grid "grd". 
  ! The determinant of Jacobian of the transformation 
  ! is returned in "JJ".

  subroutine comp_Jstar_p1(grd, ielem, r, s, jac, Jstar, JJ)
    implicit none
    type(grid), intent(in) :: grd
    integer, intent(in) :: ielem
    real*8, intent(in) :: r, s
    real*8, dimension(:,:), intent(out) :: jac, Jstar
    real*8, intent(out) :: JJ

    ! local vars
    integer :: i
    real*8 :: val, xi, yi
    real*8, dimension(2), save :: der
    real*8 :: dx_dr, dx_ds, dy_dr, dy_ds

    ! hard reset
    jac = 0.0d0; Jstar = 0.0d0; JJ = 0.0d0
    val = 0.0d0; xi = 0.0d0; yi = 0.0d0; der = 0.0d0
    dx_dr = 0.0d0; dx_ds = 0.0d0; dy_dr= 0.0d0; dy_ds = 0.0d0

    ! compute the components of jac
    do i = 1, 3 ! assuming iso-geometric expansion for x, y
       call psi(1, i, r, s, val, der)
       xi = grd%x(grd%icon(ielem,i))
       yi = grd%y(grd%icon(ielem,i))
       dx_dr = dx_dr + xi * der(1)
       dx_ds = dx_ds + xi * der(2)
       dy_dr = dy_dr + yi * der(1)
       dy_ds = dy_ds + yi * der(2) 
    end do

    ! compute jac
    jac(1,1) = dx_dr; jac(1,2) = dy_dr
    jac(2,1) = dx_ds; jac(2,2) = dy_ds

    ! comp JJ, i.e. det of jac
    JJ = dx_dr * dy_ds - dx_ds * dy_dr

    ! check it, should be valid grid!!!
    if (JJ <= 0.0d0) then 
       print *, 'error : negative Jacobian at element (',ielem,')'
       stop
    end if

    ! comp inverse of jac and store it in Jstar
    Jstar(1,1) = 1.0d0 / JJ * jac(2,2)
    Jstar(2,2) = 1.0d0 / JJ * jac(1,1)
    Jstar(1,2) = -1.0d0 / JJ * jac(1,2)
    Jstar(2,1) = -1.0d0 / JJ * jac(2,1)

    ! done here
  end subroutine comp_Jstar_p1



  ! computes the gradient at the interior  
  ! nodes of an element using the values of
  ! the function given at nodes and the 
  ! gradient of the basis function
  !
  !
  ! i.e. for element i Gauss nodej, we have:
  !
  !    grad u(Gauss node j) = sum_k (u_k *grad psi_k(rj, sj) )
  !
  !

  subroutine comp_grad_u(u, grd, dudx, dudy)
    implicit none
    real*8, dimension(:,:), intent(in) :: u
    type(grid), intent(in) :: grd
    ! dimension(elem%neqs, grd%ntri)
    real*8, dimension(:,:), intent(out) :: dudx, dudy

    ! local vars
    integer :: i, k, ptk
    real*8 :: r, s, val, JJ
    real*8, dimension(2), save :: der
    real*8, dimension(2,2), save :: jac, Jstar

    ! hard reset
    dudx = 0.0d0
    dudy = 0.0d0

    r = 0.25d0; s = 0.25d0
    do i = 1, grd%ntri
       ! subroutine comp_Jstar_p1(grd, ielem, r, s, jac, Jstar, JJ)
       call comp_Jstar_p1(grd, i, r, s, jac, Jstar, JJ)
       do k = 1, 3
          ptk = grd%icon(i,k)
          val = 0.0d0; der = 0.0d0
          call psi(1, k, r, s, val, der)
          der = matmul(Jstar, der)
          dudx(:,i) = dudx(:,i) + u(:,ptk) * der(1)
          dudy(:,i) = dudy(:,i) + u(:,ptk) * der(2)
       end do
    end do

    ! done here
  end subroutine comp_grad_u

  !              only P1
  ! computes the interpolated value of the  
  ! primary variable  <u> using the values of
  ! <u> given at interpolation nodes
  ! and the value of the basis functions at those points.
  !
  !
  ! i.e. for element i at nodej, we have:
  !
  !    u(node j) = sum_k (u_k * psi_k(rj, sj) )
  !
  ! where k = 1..3
  !

  subroutine comp_u_p1(u, r, s, u_out)
    implicit none
    real*8, dimension(:,:), intent(in) :: u
    real*8, intent(in) :: r, s
    ! dimension(elem%neqs)
    real*8, dimension(:), intent(out) :: u_out

    ! local vars
    integer :: k
    real*8 :: val
    real*8, dimension(2), save :: der

    ! bullet proofing
    if(size(u, 2) > 3) then
       print *, 'u should be u(1..neqs, 1..3) for P1' &
            , ' subtriangle interpolation. stop.'
       stop
    end if

    ! hard reset
    u_out = 0.0d0
    do k = 1, 3
       val = 0.0d0
       call psi(1, k, r, s, val, der)
       u_out = u_out + u(:,k) * val
    end do

    ! done here
  end subroutine comp_u_p1

  !
  ! visualizes the curved element grid
  ! both values of the function and gradient-based values
  !
  subroutine vis_curved_grid(fem, outfile)
    implicit none
    type(fem_struct), intent(in) :: fem
    character(len = *), intent(in) :: outfile

    ! local vars
    integer :: ii, neqs, jj
    type(grid) :: grd2
    real*8, dimension(:,:), allocatable :: utmp, d_utmp_dx, d_utmp_dy, u_cell

    ! init
    neqs = fem%neqs

    do ii = 1, fem%grd%ncellsg

       ! produce sub-triangle for a single curved SEM element
       ! just to use in the visuallization procedure
       select case (fem%grd%elname(ii))

       case (GEN_TRIANGLE)
          call mesh_a_triangle(grd_in = fem%grd, ielem = ii, grd_out = grd2)
       case (GEN_QUADRI)
          call mesh_a_quadri(grd_in = fem%grd, ielem = ii, grd_out = grd2)
       case default
          print *, 'unknown element name in vis_curved_grid(...)! stop'
          stop

       end select

       allocate(utmp(neqs, grd2%nnodesg))
       allocate(d_utmp_dx(neqs, grd2%ntri), d_utmp_dy(neqs, grd2%ntri))
       allocate(u_cell(neqs, grd2%ntri))

       ! project fem%u on sub-element triangulation
       do jj = 1, grd2%nnodesg
          utmp(:, jj) = fem%u(:, fem%grd%icon(ii, jj))
       end do

       ! append the sub-element triangulation as a 
       ! new zone to already generated zonal grid 
       if ( ii .eq. 1 ) then
          call write_u_tecplot(outfile, grd2, utmp)
       else
          call write_u_tecplot(outfile, grd2, utmp, .true.)
       end if

       ! now proceed to write gradient-based values
       !
       ! first compute the gradient for all sub-element triangles
       call comp_grad_u(u = utmp, grd = grd2, dudx = d_utmp_dx, dudy = d_utmp_dy)
       ! then call the function like CP here
       u_cell = 1.0d0 - (d_utmp_dx**2.0d0 + d_utmp_dy**2.0d0) 
       ! then append
       if ( ii .eq. 1 ) then
          call write_u_cell_centered_tecplot('grad_'//outfile, grd2, u_cell)
       else
          call write_u_cell_centered_tecplot('grad_'//outfile, grd2, u_cell, .true.)
       end if

       ! little clean-up
       deallocate(utmp)
       deallocate(d_utmp_dx, d_utmp_dy)
       deallocate(u_cell)

    end do

    ! done here
  end subroutine vis_curved_grid

  ! recursively generates adapted subtriangulation
  ! such that all gradients are resolved in one
  ! very high-order element using a lot of small P1
  ! elements
  recursive subroutine gen_tris(u, grd, elem, ielem, tris, x, y, l, nnodes &
                                   , ntris, tol, tol_sub, pt_tags, tt, func, n_rec &
                                   , min_area, critic_area, npt_err_integ)
    implicit none
    real*8, dimension(:, :), intent(in) :: u
    type(grid), target, intent(in) :: grd
    type(element), intent(inout) :: elem
    integer, intent(in) :: ielem
    integer, dimension(:,:), allocatable, intent(inout) :: tris
    real*8, dimension(:), allocatable, intent(inout) :: x, y
    integer, intent(inout) :: l, nnodes, ntris 
    real*8, intent(in) :: tol ! tol of snapping 
    real*8, dimension(:), intent(in) :: tol_sub !tol of sub triangulation
    integer, dimension(:), allocatable, intent(inout) :: pt_tags
    real*8, dimension(:), allocatable, intent(inout) :: tt
    interface
       subroutine func(xx, yy, val)
         implicit none
         real*8, intent(in) :: xx, yy
         real*8, intent(out) :: val
       end subroutine func
    end interface
    integer, dimension(3), intent(inout) :: n_rec
    real*8, intent(in) :: min_area, critic_area
    character(len = *), intent(in) :: npt_err_integ

    ! local vars
    integer :: i, ptA, ptB, ptC
    real*8, dimension(3) :: xi, yi, this_tt, rr, ss
    integer, dimension(3) :: pt, tg
    real*8 :: xA, xB, xC, yA, yB, yC, tt0
    real*8, dimension(:), allocatable :: tmp
    integer, dimension(:, :), allocatable :: itmp
    integer, dimension(:), allocatable :: tag_tmp
    type(curve), pointer :: tcurve=>null()
    real*8, dimension(:), pointer :: t=>null(), xs=>null(), ys=>null()

    ! logic
    integer :: npt_mid
    real*8 :: r1x, r1y , r2x, r2y, area, target_area, tmp_r, tmp_s
    real*8 :: xmid, ymid, mid_JJ
    real*8, dimension(3) :: rmid, smid, Wmid
    real*8, dimension(2, 2) :: mid_jac, mid_Jstar
    real*8, dimension(size(u,1)) :: uP1_mid, u_spectral_mid, int_error, int_u
    real*8, dimension(size(u,1)) :: rel_error
    real*8, dimension(size(u,1), 3) :: uP1
    real*8, dimension(1,3) :: xi_mid, yi_mid
    real*8, dimension(1) :: mid_tmp

    ! bullet proofing
    if ( l .eq. 1 ) then ! only initially!
       if (      (.not. allocated(tris)) & 
            .or. (.not. allocated(x)   ) &
            .or. (.not. allocated(y)   )  )   then
          print *, 'fatal : first allocate mother triangle '&
               , ' in order to do sub triangulation. stop.'
          stop
       elseif ( (size(x) .ne. size(y))  &
           .or. (size(x) .ne. nnodes)    ) then
          print *, 'number of input nnodes does not matches '&
               , '  with size of <x> or <y> vectors. stop'
          stop
       elseif ( size(tris, 1) .ne. ntris ) then
          print *, 'number of triangles is not equal to the '&
               , ' number of triangles in connectivity matrix. stop'
          stop
       end if
    end if
    if (grd%ntri <= 0 ) then
       print *, 'the input grid needs to be triangulated' &
            , '  first befor subtriangulation. stop'
       stop
    end if
    if (.not. allocated(grd%bn_curves)) then
       print *, 'fatal: the input grid does not have parametric' &
            , ' boundary curves! stop'
       stop
    end if


    ! finding three vetrices of the mother triangle
    do i = 1, 3
       pt(i) = tris(l, i)
       xi(i) = x(pt(i)); yi(i) = y(pt(i))
       tg(i) = pt_tags(pt(i)) 
    end do


    ! computing mid points
    this_tt = -1.0d0 ! only for internal edges
    ! straight edge by default
    xA = 0.5d0 * (xi(1) + xi(2)); yA = 0.5d0 * (yi(1) + yi(2))
    ! snap it to curved boundary if required
    if ( (tg(1) .eq. tg(2)) .and. (tg(1) .ne. 0) ) then ! is bn edge
       tt0 = 0.5d0 * (tt(pt(1)) + tt(pt(2)))
       ! find new snapped coords on the curve
       tcurve => grd%bn_curves(tg(1))
       t => tcurve%t
       xs => tcurve%x
       ys => tcurve%y
       call spline_nurbs_eval2(tt0, xs, ys, tcurve%a, tcurve%b, tcurve%c, &
                    & tcurve%Mx, tcurve%My, t, xA, yA, 'interp', tcurve%btype)
!      call spline_eval2(t, xs, tt0, xA, 'interp', tcurve%Mx)
!      call spline_eval2(t, ys, tt0, yA, 'interp', tcurve%My)
       this_tt(1) = tt0
    end if
    ! straight edge by default
    xB = 0.5d0 * (xi(2) + xi(3)); yB = 0.5d0 * (yi(2) + yi(3))
    ! snap it to curved boundary if required
    if ( (tg(2) .eq. tg(3)) .and. (tg(2) .ne. 0) ) then ! is bn edge
       tt0 = 0.5d0 * (tt(pt(2)) + tt(pt(3)))
       ! find new snapped coords on the curve
       tcurve => grd%bn_curves(tg(2))
       t => tcurve%t
       xs => tcurve%x
       ys => tcurve%y
       call spline_nurbs_eval2(tt0, xs, ys, tcurve%a, tcurve%b, tcurve%c, &
                    & tcurve%Mx, tcurve%My, t, xB, yB, 'interp', tcurve%btype)
!      call spline_eval2(t, xs, tt0, xB, 'interp', tcurve%Mx)
!      call spline_eval2(t, ys, tt0, yB, 'interp', tcurve%My)
       this_tt(2) = tt0
    end if
    ! straight edge by default
    xC = 0.5d0 * (xi(3) + xi(1)); yC = 0.5d0 * (yi(3) + yi(1))
    ! snap it to curved boundary if required
    if ( (tg(3) .eq. tg(1)) .and. (tg(3) .ne. 0) ) then ! is bn edge
       tt0 = 0.5d0 * (tt(pt(3)) + tt(pt(1)))
       ! find new snapped coords on the curve
       tcurve => grd%bn_curves(tg(3))
       t => tcurve%t
       xs => tcurve%x
       ys => tcurve%y
       call spline_nurbs_eval2(tt0, xs, ys, tcurve%a, tcurve%b, tcurve%c, &
                    & tcurve%Mx, tcurve%My, t, xC, yC, 'interp', tcurve%btype)
!      call spline_eval2(t, xs, tt0, xC, 'interp', tcurve%Mx)
!      call spline_eval2(t, ys, tt0, yC, 'interp', tcurve%My)
       this_tt(3) = tt0
    end if

    ! put the logic for refinement here between two lines
    ! ---------------------------------
    ! compute the area of this sub-division
    r1x = xi(2) - xi(1)
    r1y = yi(2) - yi(1)
    r2x = xi(3)- xi(1)
    r2y = yi(3)- yi(1)
    area = 0.5d0 * abs(r1x * r2y - r2x * r1y)
    ! call func((1.0d0 / 3.0d0 * sum(xi)), (1.0d0 / 3.0d0 * sum(yi)), target_area)
    ! print *, 'hey abs(area - target_area) = ', abs(area - target_area) &
    !      , 'area = ', area, 'target area = ', target_area

    ! interpolate the value of <u> from mother high-order element
    ! on the vertices of current P1 : (xi, yi) element 
    do i = 1, 3
       call xy2rs(xi(i), yi(i), elem, ielem, grd, 40, tol, rr(i), ss(i))
       call comp_u_point(u, grd, elem, ielem, rr(i), ss(i), uP1(:,i))
    end do

    ! deciding on the order of error integration formulae.
    ! simply decides on the input string for the number 
    ! of desired Gaussian quadrature points for evaluating
    ! the error integral
    if ( npt_err_integ .eq. '1pt') then
       npt_mid = 1
       rmid = (/ 1.0d0 / 3.0d0, 0.0d0, 0.0d0 /) 
       smid = (/ 1.0d0 / 3.0d0, 0.0d0, 0.0d0 /)
       Wmid = (/ 1.0d0 / 2.0d0, 0.0d0, 0.0d0 /) 
    elseif ( npt_err_integ .eq. '3pt') then
       npt_mid = 3
       rmid = (/ 1.0d0 / 6.0d0, 2.0d0 / 3.0d0, 1.0d0 / 6.0d0 /) 
       smid = (/ 1.0d0 / 6.0d0, 1.0d0 / 6.0d0, 2.0d0 / 3.0d0 /)
       Wmid = (/ 1.0d0 / 6.0d0, 1.0d0 / 6.0d0, 1.0d0 / 6.0d0 /) 
    else
       print *, 'unknown quadrature formulae for evaluating' &
            , '  error integral in tri-subdivision! stop'
       stop
    end if

    ! starting computing error integral ...
    xi_mid(1, :) = xi; yi_mid(1, :) = yi
    int_error = 0.0d0; int_u = 0.0d0
    
    do i = 1, npt_mid
       ! compute the physical (x, y) of the current Gauss point
       call comp_u_p1(xi_mid, rmid(i), smid(i), mid_tmp)
       xmid = mid_tmp(1)  
       call comp_u_p1(yi_mid, rmid(i), smid(i), mid_tmp)
       ymid = mid_tmp(1)  

       ! compute uP1 at the current Gauss point
       call comp_u_p1(uP1, rmid(i), smid(i), uP1_mid)
       ! compute the spectral <u> of mother triangle at the current Gauss point
       call xy2rs(xmid, ymid, elem, ielem, grd, 40, tol, tmp_r, tmp_s)
       call comp_u_point(u, grd, elem, ielem, tmp_r, tmp_s, u_spectral_mid)
       ! compute the spectral Jacobian on the mother triangle at the current GP
       call comp_Jstar_point(grd, elem, ielem, tmp_r, tmp_s &
            , mid_jac, mid_Jstar, mid_JJ)

       ! finally, computing the error integrals
       int_error = int_error + abs(u_spectral_mid - uP1_mid) * Wmid(i)   
       int_u = int_u + abs(u_spectral_mid) * Wmid(i) * mid_JJ
    end do

    ! find the relative error
    ! rel_error = int_error / int_u
    rel_error = int_error
    ! print *, 'rel_error = ', rel_error
 
    if ( all(rel_error <= tol_sub) .and. &
         (n_rec(1) >= n_rec(2)) .and. &
         (area <= min_area) ) then
       l = l + 1
       return
    end if

    if ( (area <= critic_area) .or. (n_rec(1) >= n_rec(3)) ) then
       print *, 'Warning : in adaptive visuallization, the minimum area ' &
              , 'or number of recursive levels passed the ciritical ' &
              , 'limits! return to coarser level ...'
       l = l + 1
       return
    end if

    ! ! PRESSURE VALVE : prevent infinite recurstion!!!!
    ! if ( n_rec(1) >= n_rec(3) ) then
    !    print *, 'warning : maximum recursion level in sub-triangulation happened!'
    !    l = l + 1
    !    return
    ! end if

    ! if ( area  <= (tol + target_area) ) then
    !    l = l + 1
    !    return
    ! end if
    ! ---------------------------------

    ! prepare for recursive adaptation
    ! 
    ! first add new nodes ...
    ptA = nnodes+1; ptB = nnodes+2; ptC = nnodes+3

    ! add "x" of new nodes
    allocate( tmp(nnodes+3) )
    tmp(1:nnodes) = x(1:nnodes)
    tmp(ptA) = xA
    tmp(ptB) = xB 
    tmp(ptC) = xC
    call move_alloc( tmp, x)
    ! add "y" of new nodes
    allocate( tmp(nnodes+3) )
    tmp(1:nnodes) = y(1:nnodes)
    tmp(ptA) = yA
    tmp(ptB) = yB 
    tmp(ptC) = yC
    call move_alloc( tmp, y)

    ! add <tt> of new nodes
    allocate( tmp(nnodes+3) )
    tmp(1:nnodes) = tt(1:nnodes)
    tmp(ptA) = this_tt(1)
    tmp(ptB) = this_tt(2) 
    tmp(ptC) = this_tt(3)
    call move_alloc( tmp, tt)

    ! add point tags
    allocate( tag_tmp(nnodes+3) )
    tag_tmp(1:nnodes) = pt_tags(1:nnodes)
    ! add tag for point A
    if ( (tg(1) .eq. tg(2)) .and. (tg(1) .ne. 0) ) then ! is boundary edge
       tag_tmp(ptA) = tg(1)
    elseif ( (tg(1) .ne. tg(2)) .and. (tg(1) .ne. 0) .and. (tg(2) .ne. 0) ) then
       print *, 'misleading point tag for point A in sub triangulation! stop'
       stop
    else  
       tag_tmp(ptA) = 0
    end if
    ! add tag for point B
    if ( (tg(2) .eq. tg(3)) .and. (tg(2) .ne. 0) ) then ! is boundary edge
       tag_tmp(ptB) = tg(2)
    elseif ( (tg(2) .ne. tg(3)) .and. (tg(2) .ne. 0) .and. (tg(3) .ne. 0) ) then
       print *, 'misleading point tag for point B in sub triangulation! stop'
       stop
    else  
       tag_tmp(ptB) = 0
    end if 
    ! add tag for point C
    if ( (tg(3) .eq. tg(1)) .and. (tg(3) .ne. 0) ) then ! is boundary edge
       tag_tmp(ptC) = tg(3)
    elseif ( (tg(3) .ne. tg(1)) .and. (tg(3) .ne. 0) .and. (tg(1) .ne. 0) ) then
       print *, 'misleading point tag for point C in sub triangulation! stop'
       stop
    else  
       tag_tmp(ptC) = 0
    end if 

    call move_alloc(tag_tmp, pt_tags)

   ! update nodal index
    nnodes = nnodes + 3


    ! now, add sub triangles
    ! 
    allocate(itmp(ntris + 3, 3))
    itmp(1:l, :) = tris(1:l, :) 
    if (ntris > 1) then
       itmp((l+4):(ntris + 3), :) = tris((l+1):ntris, :)
    end if

    itmp(  l  , :) = (/ pt(1), ptA, ptC /)
    itmp((l+1), :) = (/ ptA, pt(2), ptB /)
    itmp((l+2), :) = (/ ptC, ptB, pt(3) /)
    itmp((l+3), :) = (/ ptA, ptB, ptC /)
    call move_alloc(itmp, tris)
    ntris = ntris + 3

    ! one level of recursion is about to start so
    ! increase its indicator
    n_rec(1) = n_rec(1) + 1

    do i = 1, 4
       call gen_tris(u, grd, elem, ielem, tris, x, y, l, nnodes &
            , ntris, tol, tol_sub, pt_tags, tt, func, n_rec &
            , min_area, critic_area, npt_err_integ)
    end do

    ! done here
  end subroutine gen_tris
  
  subroutine sample_func(xx, yy, val)
    implicit none
    real*8, intent(in) :: xx, yy
    real*8, intent(out) :: val

    ! local cars
    real*8 :: r2

    ! r2 = (xx - 1.5d0)**2 +  (yy - 1.0d0)**2

    ! val = 1.0d-3 * abs(cos( 4.0d0 * exp(sqrt(r2)) ) )
    val = 0.05d0

    ! done here
  end subroutine sample_func

  ! a little tester for subtriangulation subroutine
  ! zone 98 = straight triangle
  ! zone 89 = curved
  subroutine test_gen_tris(grd, u, elem, ielem, tol, tol_sub &
       , outfile, appendit, n_rec, min_area, critic_area, npt_err_integ)
    implicit none
    type(grid), intent(in) :: grd
    real*8, dimension(:,:), intent(in) :: u
    type(element), intent(inout) :: elem
    integer, intent(in) :: ielem
    real*8, intent(in) :: tol
    real*8, dimension(:), intent(in) :: tol_sub
    character(len = *), intent(in) :: outfile
    logical, intent(in) :: appendit
    integer, dimension(3), intent(inout) :: n_rec
    real*8, intent(in) :: min_area, critic_area
    character(len = *), intent(in) :: npt_err_integ

    ! local vars
    integer , dimension(:,:), allocatable :: tris
    real*8, dimension(:), allocatable :: x, y, tt
    real*8, dimension(:,:), allocatable:: u_out
    integer, dimension(:), allocatable :: pt_tags
    type(grid) :: grd_out

    integer :: i, l, nnodes, ntris, edgtag, edgnum, edgloc
    real*8 :: tt1, tt2, rr, ss


    ! a sample test case
    allocate(x(3), y(3) )  
    do i = 1, 3
       x(i) = grd%x(grd%icon(ielem, i))
       y(i) = grd%y(grd%icon(ielem, i))
    end do

    allocate(pt_tags(3), tt(3))

    edgtag = grd%el2bn(ielem, 1) 

    if (  edgtag .eq. 0 ) then ! interior element
       pt_tags = 0
       tt = -1.0d0
    else
       edgnum = grd%el2bn(ielem, 2)
       edgloc = grd%ibedgeELEM_local_edg(edgnum)

       select case (edgloc)
       case (1)
          pt_tags(1:2) = edgtag; pt_tags(3) = 0
          call find_t(grd, edgtag, x(1), y(1), tol, tt1)
          if ( tt1 .eq. -1.0d0) then
             print *, 'error : can not evaluate tt1. stop.'
             stop
          end if
          call find_t(grd, edgtag, x(2), y(2), tol, tt2)
          if ( tt2 .eq. -1.0d0) then
             print *, 'error : can not evaluate tt2. stop.'
             stop
          end if
          tt = (/ tt1, tt2, -1.0d0 /)

       case (2)
          pt_tags(2:3) = edgtag; pt_tags(1) = 0
          call find_t(grd, edgtag, x(2), y(2), tol, tt1)
          if ( tt1 .eq. -1.0d0) then
             print *, 'error : can not evaluate tt1. stop.'
             stop
          end if
          call find_t(grd, edgtag, x(3), y(3), tol, tt2)
          if ( tt2 .eq. -1.0d0) then
             print *, 'error : can not evaluate tt2. stop.'
             stop
          end if
          tt = (/ -1.0d0, tt1, tt2 /)

       case (3)
          pt_tags(3) = edgtag; pt_tags(1) = edgtag; pt_tags(2) = 0
          call find_t(grd, edgtag, x(3), y(3), tol, tt1)
          if ( tt1 .eq. -1.0d0) then
             print *, 'error : can not evaluate tt1. stop.'
             stop
          end if
          call find_t(grd, edgtag, x(1), y(1), tol, tt2)
          if ( tt2 .eq. -1.0d0) then
             print *, 'error : can not evaluate tt2. stop.'
             stop
          end if
          tt = (/ tt2, -1.0d0, tt1 /)
 
       case default
          print *, 'error: unknown local edge in test_gen_tris(...). stop'
          stop
       end select

    end if
    ! print *, 'pt_tags = ', pt_tags
    ! print *, 'tt = ', tt
    ! stop
    allocate( tris(1, 3) )
    tris(1, :) = (/ 1, 2, 3 /)
    ntris = 1
    nnodes = 3
    l = 1

    call gen_tris(u, grd, elem, ielem, tris, x, y, l, nnodes &
         , ntris, tol, tol_sub, pt_tags, tt, sample_func, n_rec & 
         , min_area, critic_area, npt_err_integ)

    ! write output as Tecplot file
    allocate(grd_out%x(nnodes), grd_out%y(nnodes))
    allocate(grd_out%icon(ntris, 3))
    allocate(u_out(size(u,1), nnodes))  
    grd_out%x = x
    grd_out%y = y
    grd_out%icon = tris
    grd_out%nnodesg = nnodes
    grd_out%ntri = ntris
    grd_out%nquad4 = 0
    grd_out%ncellsg = ntris

    ! computing <u> at the generated sub triangle grid
    do i = 1, nnodes
       call xy2rs(grd_out%x(i), grd_out%y(i), elem, ielem, grd, 40, tol, rr, ss)
       call comp_u_point(u, grd, elem, ielem, rr, ss, u_out(:,i))
    end do
    ! do i = 1, nnodes
    !    call sample_func(grd_out%x(i), grd_out%y(i), u(1, i) )
    ! end do

    call write_u_tecplot(outfile, grd_out, u_out, appendit)

    ! clean up    
    deallocate(tris, x, y, tt, pt_tags, u_out)

    ! done here
  end subroutine test_gen_tris

  !
  ! visualizes the curved element grid
  ! using adaptive sub triangulation
  !
  subroutine adaptive_vis_curved_grid(fem, outfile, tol, tol_sub &
       , n_rec, min_area, critic_area, npt_err_integ)
    implicit none
    type(fem_struct), intent(inout) :: fem
    character(len = *), intent(in) :: outfile
    real*8, intent(in) :: tol, tol_sub
    integer, dimension(2), intent(in) :: n_rec !(/ min_num_recursion, max_num_rec/)
    real*8, intent(in) :: min_area, critic_area
    character(len = *), intent(in) :: npt_err_integ

    ! local vars
    integer :: ii
    logical :: appendit_flag
    real*8, dimension(fem%neqs) :: tol_sub_fin
    integer, dimension(3) :: n_rec_fin

    appendit_flag = .false.

    ! finding the average of the solution field
    ! such that it is multiplied by tol_sub to get
    ! the final tolerance of sub triangulaition
    do ii = 1, fem%neqs
       tol_sub_fin(ii) = tol_sub * sqrt(sum(fem%u(ii, :)**2.0d0))
    end do


    do ii = 1, fem%grd%ntri

       ! filling the number of recursion information array
       n_rec_fin(2:3) = n_rec ! (min possible, max possible)
       n_rec_fin(1) = 0 ! current level

       call test_gen_tris(fem%grd, fem%u, fem%elems(ii), ii, tol &
            , tol_sub_fin, outfile, appendit_flag, n_rec_fin &
            , min_area, critic_area, npt_err_integ)
       appendit_flag = .true.
       print *, 'element ', ii, ' sub-divided successfully!'

    end do

    ! done here
  end subroutine adaptive_vis_curved_grid



end module plot_curved_elems
