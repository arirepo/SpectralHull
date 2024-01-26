module quadri_elem
  use globals
  implicit none


  private


  public :: gen_master_quadri, pt_dist_1d, fill_elem_coords

contains

  !> @details Generates 1D point distribution for sides of
  !! a spectral element.
  subroutine pt_dist_1d(n, xmin, xmax, method, x)
    implicit none
    integer, intent(in) :: n	!< Number of desired points (input).
    real*8, intent(in) :: xmin, xmax	!< Max and min position in 1D (input).
    character(len = *), intent(in) :: method	!< Method of interest to generate points (input).
    real*8, dimension(:), intent(out) :: x	!< Output: array of the generated points.

    ! local vars
    integer :: i

    ! bullet proofing
    if ( size(x) .ne. n ) then
       print *, 'incorrect size of 1d point distribution! stop'
       stop
    end if

    select case (method)

    case ('equal_space')

       x = (/ ( xmin + dble(i-1) / dble(n - 1) * (xmax - xmin) , i = 1, n ) /)

    case ('cheby')

       x = cos( PI * (/ ( dble(i), i = (n-1), 0, -1) /) / dble(n - 1) )

       case default 
 
         print *, 'Can not recognize the method for 1d point distribution. stop!!'
         stop
 
     end select
 
 
     ! done here
   end subroutine pt_dist_1d
 
  !> @details Packs coordinates of the 1D distributed points x1, x2
  !> into two-dimensional quadrilateral 
  !>
  !>        edge3
  !>  4 --------------- 3
  !>  |                 |
  !>  |                 |
  !>  | edge4           | edge2
  !>  |                 |
  !>  1 ----------------2
  !>         edge1
  !>  
  subroutine fill_elem_coords(x1, x2, x, y)
    implicit none
    real*8, dimension(:), intent(in), target :: x1  !< Input array in the 1st dimension.
	real*8, dimension(:), intent(in), target :: x2  !< Input array in the 2nd dimension.
    real*8, dimension(:), allocatable :: x  !< Output array in the 1st dimension.
	real*8, dimension(:), allocatable :: y  !< Output array in the 2ns dimension.

    ! local vars
    integer :: i, j, n1, n2
    real*8, dimension(:), pointer :: edg => null()
    real*8, dimension(:), allocatable :: tmp

    ! bullet proofing
    if (allocated(x)) then
       if (size(x) > 0) then 
          print *, 'Cannot fill quad4 element because <x> might has' &
               , ' been already filled! try deallocating x! stop!'
          stop
       end if
    end if
    if (allocated(y)) then
       if (size(y) > 0) then 
          print *, 'Cannot fill quad4 element because <y> might has' &
               , ' been already filled! Try deallocating y! stop!'
          stop
       end if
    end if

    ! init
    n1 = size(x1); n2 = size(x2)

    ! add pt1
    call concat(x, x1(1:1));   call concat(y, x2(1:1))
    ! add pt2
    call concat(x, x1(n1:n1)); call concat(y, x2(1:1))
    ! add pt3
    call concat(x, x1(n1:n1)); call concat(y, x2(n2:n2))
    ! add pt4
    call concat(x, x1(1:1));   call concat(y, x2(n2:n2))

    ! add edg1
    edg => x1(2:(n1-1))
    allocate(tmp(size(edg))); tmp = x2(1)
    call concat(x, edg); call concat(y, tmp)
    deallocate(tmp)
    ! add edg2
    edg => x2(2:(n2-1))
    allocate(tmp(size(edg))); tmp = x1(n1)
    call concat(x, tmp); call concat(y, edg)
    deallocate(tmp)
    ! add edg3
    edg => x1((n1-1):2:(-1))
    allocate(tmp(size(edg))); tmp = x2(n2)
    call concat(x, edg); call concat(y, tmp)
    deallocate(tmp)
    ! add edg4
    edg => x2((n2-1):2:(-1))
    allocate(tmp(size(edg))); tmp = x1(1)
    call concat(x, tmp); call concat(y, edg)
    deallocate(tmp)

    ! add interior points
    do i = 2, (n1-1)
       do j = 2, (n2-1)
          call concat(x, x1(i:i))
          call concat(y, x2(j:j))
       end do
    end do

    if ((size(x) .ne. (n1*n2)) .or. (size(y) .ne. (n1*n2))) then
       print *, 'the length of coord vectors x or y is wrong! stop'
       stop
    end if

    ! done here
  end subroutine fill_elem_coords

  !> @details Concatenates the input array <b> to the end of
  !> array <a>
  subroutine concat(a, b)
    implicit none
    real*8, dimension(:), allocatable :: a  !< The host array
    real*8, dimension(:), intent(in)  :: b  !< The array to be concatenated.

    ! local vars
    integer :: na, nb, total
    real*8, dimension(:), allocatable :: tmp

    ! bulletproofing
    if ( allocated(a) ) then
       na = size(a)
    else
       na = 0
    end if

    nb = size(b)
    total = na + nb
    allocate(tmp(total))

    if ( na > 0 ) tmp(1:na) = a
    tmp((na + 1):total) = b

    call move_alloc(tmp, a)

    ! done here
  end subroutine concat

  !> @details Generates arbitrary high-order points for a master quadrilateral element
  subroutine gen_master_quadri(n1, n2, xmin, xmax, ymin, ymax, method, x, y)
    implicit none
    integer, intent(in) :: n1  !< Number of high-order points along x direction.
	integer, intent(in) :: n2  !< Number of high-order points along y direction.
    real*8, intent(in) :: xmin  !< Lower-bound of the domain in the x-direction.
	real*8, intent(in) :: xmax  !< Upper-bound of the domain in the x-direction.
	real*8, intent(in) :: ymin  !< Lower-bound of the domain in the y-direction.
	real*8, intent(in) :: ymax  !< Upper-bound of the domain in the y-direction.
    character(len = *), intent(in) :: method !< Method of interest to generate points.
    real*8, dimension(:), allocatable :: x  !< Output: coordinates of the points in the x-direction.
	real*8, dimension(:), allocatable :: y  !< Output: coordinates of the points in the y-direction.

    ! local vars
    real*8 :: x1d(n1), y1d(n2)

    ! generate 1d point distribution on sides
    call pt_dist_1d(n1, xmin, xmax, method, x1d)
    call pt_dist_1d(n2, ymin, ymax, method, y1d)

    ! fill out the final master element using
    ! cartesian tensor product 
    call fill_elem_coords(x1d, y1d, x, y)

    ! done here
  end subroutine gen_master_quadri

end module quadri_elem

! program tester 
!   use gnufor2
!   use quadri_elem
!   implicit none

!   ! local vars
!   integer :: i, nx, ny
!   real*8, dimension(:), allocatable :: x, y

!   nx = 10; ny = 20
!   call  gen_master_quadri(nx, ny, -1.0d0 &
!                         , 1.0d0, -1.0d0 &
!                         , 1.0d0, 'cheby', x, y)

 
!   ! write to file
!   open(unit = 1, file = 'pts.dat')
!   do i = 1, size(x)
!      write(1, *) x(i), y(i)
!   end do
!   close(1)

!   ! clean ups
!   if (allocated(x)) deallocate(x)
!   if (allocated(y)) deallocate(y)

!   ! done here
! end program tester

