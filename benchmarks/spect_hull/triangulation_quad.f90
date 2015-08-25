module triangulation_quad
  use grid_opt
  use ps
  use trimesher
  use lag_basis
  use dunavant
  implicit none

  private



  public :: gen_triang_quad

contains

  ! triangulates a polygon and returns 
  ! connectivity array "icon" and coordinates
  ! of the points "x" and "y"
  !
  subroutine triangulate_polygon(ar, icon, x, y)
    implicit none
    real*8, dimension(:, :), intent(in) :: ar
    integer, dimension(:, :), allocatable :: icon
    real*8, dimension(:), allocatable :: x, y

    ! local vars
    type(grid) :: grd_out
    integer :: npt, ii, pt1, pt2, indx
    type(point), dimension(:), allocatable :: pts
    type(connector), dimension(:), allocatable :: cons 
    type(point) :: holes(1)
    integer :: report_before_in, report_after_in
    integer, dimension(2) :: pts_list

    ! init.
    report_before_in = 0 
    report_after_in = 14

    ! allocate space for all points
    npt = size(ar, 1) - 1 ! one duplicate at the end!
    allocate(pts(npt)) 

    ! now fill the x/y of thoses skeleton points
    do ii = 1, npt ! loop over all edges
       ! add the point of the start of this edge
       pts(ii)%x = ar(ii, 1)
       pts(ii)%y = ar(ii, 2)
       pts(ii)%tag = ii
    end do

    ! then dimension all edges
    allocate(cons(npt))
    indx = 1

    ! proceeding to fill connectors array
    do ii = 1, npt ! loop over all edges
       pt1 = ii
       pt2 = ii + 1
       if ( ii .eq. npt ) pt2 = 1 !loop closed

       pts_list = (/ pt1, pt2 /)

       call add_to_cons(pts, pts_list, indx, cons)

    end do

    ! setting up the holes which there is no hole
    ! here in this case
    holes(1)%tag = -1
    holes(1)%x = 1.0d20 ! some point outside the domain
    holes(1)%y = 1.0d20 ! some point outside the domain

    ! perform triangulation
    call meshit_visual('j', pts, cons, holes, report_before_in &
         , report_after_in, grd_out)

    ! setting up the outputs
    !
    ! set "icon"
    if ( allocated(icon) ) deallocate(icon)
    allocate( icon(size(grd_out%icon, 1), size(grd_out%icon, 2)) )
    icon = grd_out%icon

    ! set "x"
    if ( allocated(x) ) deallocate(x)
    allocate(x(size(grd_out%x)))
    x = grd_out%x

    ! set "y"
    if ( allocated(y) ) deallocate(y)
    allocate(y(size(grd_out%x)))
    y = grd_out%y

    ! clean ups
    deallocate(pts, cons)

    ! done here
  end subroutine triangulate_polygon


  ! main subroutine for generating Gauss-legendre
  ! quadrature for polygons based on triangulation
  ! of the polygon and using Gauss-Legendre on each
  ! sub-triangle of the polygon
  !
  subroutine gen_triang_quad(order, ar, xy, W)
    implicit none
    integer, intent(in) :: order
    real*8, dimension(:, :), intent(in) :: ar
    real*8, dimension(:, :), allocatable :: xy
    real*8, dimension(:), allocatable :: W

    ! local vars
    integer :: nt, degree, ngauss, indx, ii, jj, kk, tpt
    real*8 :: psi_kk, der_kk(2), area
    integer, dimension(:, :), allocatable :: icon
    real*8, dimension(:), allocatable :: x, y
    real*8, dimension(:, :), allocatable :: xy0
    real*8, dimension(:), allocatable :: W0
    real*8 :: node_xy(2,3)

    ! first perform triangulation and obtain
    ! a list of triangles and their coordinates
    ! -----------------------------------------------
    call triangulate_polygon(ar, icon, x, y)
    nt = size(icon, 1) ! number of triangles

    ! then generate ONE refrence Gauss-Legendre
    ! quad for the given "order" for a refrence triangle 
    ! defined on 0 <= r <= 1 and 0 <= s <= 1
    ! -----------------------------------------------
    ! check to see the order of exactness is available
    ! in the tables for the given rule
    call dunavant_degree ( order, degree )

    ! compute the number of required Gauss points
    call dunavant_order_num( order, ngauss )

    ! allocate space for that
    allocate( xy0(2,ngauss), W0(ngauss) )

    ! compute the absicca and weights for that rule
    call dunavant_rule( order, ngauss, xy0, W0 )

    ! then allocate the output arrays based on
    ! total number of GL quadrature points
    ! -----------------------------------------------
if ( allocated( xy) ) deallocate(xy)
if ( allocated( W) ) deallocate(W)

    allocate(xy(2, (nt * ngauss) ), W(nt * ngauss) )
    xy = 0.0d0
    W = 0.0d0
    indx = 1

    ! then fill the output arrays using standard FEM P1
    ! transformation method
    do ii = 1, nt

       ! subtriangle operations
       do jj = 1, ngauss

          do kk = 1, 3

             call psi(etype= 1, i= kk, r = xy0(1, jj), s= xy0(2, jj) &
                  , val= psi_kk, der= der_kk)

             tpt = icon(ii, kk)

             ! capture nodes of this triangle
             node_xy(1, kk) = x(tpt)
             node_xy(2, kk) = y(tpt)

             ! transform x, y
             xy(1, indx) = xy(1, indx) + psi_kk * x(tpt)
             xy(2, indx) = xy(2, indx) + psi_kk * y(tpt)

          end do !loop over vertices

          !
          call triangle_area ( node_xy, area )

          ! Bullet proofing
          if ( area <= 0.0d0 ) then
             print *, 'winding of subtriangle in gen_triang_quad(...)' &
                  , ' is incorrect! stop'
             stop
          end if

          ! transform W
          W(indx) = area * W0(jj)

          ! next total Gauss point 
          indx = indx + 1

       end do ! loop over each subtriangle Gauss points

    end do ! loop over all sub triangles

    ! cleanups
    if ( allocated(icon) ) deallocate(icon)
    if ( allocated(x) ) deallocate(x)
    if ( allocated(y) ) deallocate(y)
    if ( allocated(xy0) ) deallocate(xy0)
    if ( allocated(W0) ) deallocate(W0)

    ! done here
  end subroutine gen_triang_quad

end module triangulation_quad

! program tester
!   use triangulation_quad
!   implicit none

!   ! local vars
!   integer :: order
!   real*8, dimension(5, 2) :: ar
!   real*8, dimension(:,:), allocatable :: xy
!   real*8, dimension(:), allocatable :: W

!   ! set the polygon and order
!   ar(1, :) = (/ 0.0d0, 0.0d0/)
!   ar(2, :) = (/ 1.0d0, 0.0d0/)
!   ar(3, :) = (/ 1.0d0, 1.0d0/)
!   ar(4, :) = (/ 0.0d0, 1.0d0/)
!   ar(5, :) = (/ 0.0d0, 0.0d0/)
!   order = 4

!   ! compute quads
!   call gen_triang_quad(order, ar, xy, W)

!   !display quads
!   print *, 'xy(1, :) = ', xy(1, :)
!   print *, 'xy(2, :) = ', xy(2, :)
!   print *, 'W = ', W
!   print *, 'sum(W) = ', sum(W)

!   ! cleanups
!   if ( allocated( xy ) ) deallocate(xy)
!   if ( allocated( W ) ) deallocate(W)

!   ! done here
! end program tester
