module bn_integral
  use grid_opt
  use element_opt
  use spline
  implicit none

  private


  public :: comp_elem_bn_integral

contains

  subroutine comp_elem_bn_integral(grd, elem, ielem, loc_edg, tol, func, sum)
    implicit none
    type(grid), target, intent(in) :: grd
    type(element), intent(inout) :: elem
    integer, intent(in) :: ielem, loc_edg
    real*8, intent(in) :: tol
    interface
       ! custom subroutine to compute the integrand over the edge
       subroutine func(elem, x0, y0, x, y, val)
         use element_opt
         implicit none
         type(element), intent(inout) :: elem
         real*8, intent(in) :: x0, y0, x, y
         real*8, dimension(:), intent(out) :: val
       end subroutine func
    end interface
    real*8, dimension(:), intent(out) :: sum

    ! local vars
    integer :: i, tag, p, edgnum, pt1, pt2, pt3, r
    type(curve), pointer :: tcurve=>null()
    real*8, dimension(:), allocatable :: xi, W, derjac, xq, yq
    real*8, dimension(:), pointer :: t=>null(), xs=>null(), ys=>null()
    real*8 :: tt, x1, y1, x2, y2, alpha, beta, xdot, ydot, x0, y0
    logical :: do_snap
    real*8, dimension(size(sum)) :: val
    real*8 :: t1, t2

print *, 'size(val) = ', size(val)
! stop

    ! hard RESET
    sum = 0.0d0
    do_snap = .false.
    p = grd%p(ielem)
    r = int(2.0d0 * p)
    allocate(xi(r), W(r), derjac(r), xq(r), yq(r))

    alpha = 0.0d0
    beta  = 0.0d0

    ! computing Legendre-Gauss-Jacobi points for integration
    ! and corresponding weight functions
    call ZEJAGA(r, alpha, beta, xi, derjac)
    call WEJAGA(r, alpha, beta, xi, derjac, W)

    ! three points on the corners of this triangle
    pt1 = grd%icon(ielem, 1)
    pt2 = grd%icon(ielem, 2)
    pt3 = grd%icon(ielem, 3)
    if     (loc_edg .eq. 1) then
       x1 = grd%x(pt1); x2 = grd%x(pt2)
       y1 = grd%y(pt1); y2 = grd%y(pt2)
    elseif (loc_edg .eq. 2) then
       x1 = grd%x(pt2); x2 = grd%x(pt3)
       y1 = grd%y(pt2); y2 = grd%y(pt3)
    elseif (loc_edg .eq. 3) then
       x1 = grd%x(pt3); x2 = grd%x(pt1)
       y1 = grd%y(pt3); y2 = grd%y(pt1)
    else
       print *, 'wrong local edg number! stop'
       stop
    end if

    ! find element tag
    tag =    grd%el2bn(ielem, 1)
    edgnum = grd%el2bn(ielem, 2)
    print *, 'tag =', tag, 'edgnum =', edgnum 
    if ( tag .ne. 0 ) then
       if (grd%ibedgeELEM_local_edg(edgnum) .eq. loc_edg) then
          ! is a boundary element 
          do_snap = .true.
          tcurve => grd%bn_curves(tag)
          t => tcurve%t
          xs => tcurve%x
          ys => tcurve%y
       end if
    end if

    ! find the constant derivatives
    xdot = 1.0d0 / 2.0d0 * (x2 - x1)
    ydot = 1.0d0 / 2.0d0 * (y2 - y1)

    do i = 1, size(xi)

       ! Mapping of coordinates of Gauss-Legendre quadrature
       ! for straight edges to physical space
       xq(i) = x1 + (xi(i) + 1.0d0) / 2.0d0 * (x2 - x1) 
       yq(i) = y1 + (xi(i) + 1.0d0) / 2.0d0 * (y2 - y1) 

       ! if boundary edge also snapp it
       if ( do_snap ) then
          call find_t(grd, tag, x1, y1, tol, t1)
          if ( t1 .eq. -1.0d0) then
             print *, 'error : the parameter t can not be computed' &
                  , ' for point (', x1, ',', y1,') on ' &
                  , 'edg # ', edgnum,'. stop.'
             stop
          end if
          call find_t(grd, tag, x2, y2, tol, t2)
          if ( t2 .eq. -1.0d0) then
             print *, 'error : the parameter t can not be computed' &
                  , ' for point (', x2, ',', y2,') on ' &
                  , 'edg # ', edgnum,'. stop.'
             stop
          end if

          call find_t(grd, tag, xq(i), yq(i), tol, tt)
          if ( tt .eq. -1.0d0) then
             print *, 'error : the parameter t can not be computed' &
                  , ' for point (', xq(i), ',', yq(i),') on ' &
                  , 'edg # ', edgnum,'. stop.'
             stop
          end if
          ! find new snapped xq(i) and yq(i) on the curve
          call spline_nurbs_eval2(tt, xs, ys, tcurve%a, tcurve%b, tcurve%c, &
                    & tcurve%Mx, tcurve%My, t, xdot, ydot, 'interp', tcurve%btype)
!         call spline_eval2(t, xs, tt, xdot, 'interp', tcurve%Mx)
!         call spline_eval2(t, ys, tt, ydot, 'interp', tcurve%My)
          ! print *, 'xq(i) = ', xq(i), 'xdot = ', xdot
          ! print *, 'yq(i) = ', yq(i), 'ydot = ', ydot
          ! print *, 'diff_x = ', (xq(i) - xdot), 'diff_y = ', (yq(i) - ydot)

          xq(i) = xdot
          yq(i) = ydot
          ! computing derivatives xdot and ydot ...
          call spline_nurbs_eval2(tt, xs, ys, tcurve%a, tcurve%b, tcurve%c, &
                     & tcurve%Mx, tcurve%My, t, xdot, ydot, 'diff1', tcurve%btype)
!         call spline_eval2(t, xs, tt, xdot, 'diff1', tcurve%Mx)
!         call spline_eval2(t, ys, tt, ydot, 'diff1', tcurve%My)
          print *, 'tag = ', tag, 'xdot = ', xdot, 'ydot = ', ydot
          print *, 'size(t) = ', size(t)
          xdot = xdot * 0.5d0 * (t2 - t1)
          ydot = ydot * 0.5d0 * (t2 - t1)

       end if

       ! convert physical (xq(i), yq(i)) to 0<=(x0, y0)<=1 in the master element
       call xy2rs(xq(i) , yq(i), elem, ielem, grd, 40, tol, x0, y0)

       ! evaluate the integrand on point 0<=(x0, y0)<=1 in the master element 
       print *, 'ielem = ', ielem
       call func(elem, x0, y0, xq(i) , yq(i), val)

       ! accumulating to the boundary integral
       sum = sum + val * sqrt(xdot**2.0d0 + ydot**2.0d0) * W(i)

    end do

    ! little clean-up
    deallocate(xi, W, derjac, xq, yq)

    ! done here
  end subroutine comp_elem_bn_integral


end module bn_integral
