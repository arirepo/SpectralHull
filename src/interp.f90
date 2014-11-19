module interp
  ! for interpolation
  implicit none

  private

  real*8, dimension(:), allocatable :: wx, wy
  real*8, dimension(:), allocatable :: ztmp

  public :: lag1d, lag2d

contains

  ! statically computes the weights and store them
  ! in the local wx and wy arrays and save them for
  ! any future use
  subroutine comp_weights(x, w)
    implicit none
    real*8, dimension(:), intent(in) :: x
    real*8, dimension(:), allocatable, intent(out) :: w

    ! local vars
    integer :: n, m, j
    real*8 :: tmp

    n = size(x)
    allocate(w(n))

    do j = 1, n

       tmp = 1.0d0

       do m = 1, n

          if (m .eq. j) cycle
          tmp = tmp * (x(j) - x(m))

       end do

       w(j) = 1.0d0 / tmp

    end do

    ! done here
  end subroutine comp_weights

  subroutine lag1d(x, z, xx, zz, w)
    implicit none
    real*8, dimension(:), intent(in) :: x,z
    real*8, intent(in) :: xx
    real*8, intent(out) :: zz
    real*8, dimension(:), allocatable :: w ! given weights

    ! local vars
    integer :: j, n 
    real*8 :: bot, top, coeff

    n = size(x)

    ! first explicitly check to see
    ! if xx is in x then no need
    ! for interpolation and return value 
    do j = 1, n

       if ( x(j) .eq. xx ) then
          zz = z(j)
          return
       end if

    end do
 
    if ( .not. allocated(w) ) call comp_weights(x, w)

    bot = 0.0d0
    top = 0.0d0

    do j = 1, n
       coeff = w(j) / (xx - x(j))
       bot = bot + coeff 
       top = top + coeff * z(j)
    end do

    if ( bot .eq. 0.0d0 ) then
       print *, 'the bottom of 1d Lagrange interpolation is zero! stop'
       stop
    end if

    zz = top / bot

    ! done here
  end subroutine lag1d

  subroutine lag2d(x, y, z, xx, yy, zz)
    implicit none
    real*8, dimension(:), intent(in) :: x, y
    real*8, dimension(:,:), intent(in) :: z

    real*8, intent(in) :: xx, yy
    real*8, intent(out) :: zz

    ! local vars
    integer :: n1, n2, j

    n1 = size(x); n2 = size(y) 

    if ( .not. allocated(ztmp) ) allocate(ztmp(n2))

    ! interpolation in the first direction
    do j = 1, n2
       call lag1d(x, z(:,j), xx, ztmp(j), wx)
    end do

    ! interpolation in the second direction
    call lag1d(y, ztmp, yy, zz, wy)

    ! done here
  end subroutine lag2d

  ! HARD reset
  subroutine dealloc_interp_module
    implicit none

    if ( allocated(wx) ) then 
       deallocate(wx)
    else
       print *, 'warning : wx already deallocated!'
    end if

    if ( allocated(wy) ) then 
       deallocate(wy)
    else
       print *, 'warning : wy already deallocated!'
    end if

    if ( allocated(ztmp) ) then 
       deallocate(ztmp)
    else
       print *, 'warning : ztmp already deallocated!'
    end if

    ! done here
  end subroutine dealloc_interp_module

end module interp

! program test_interp
!   use gnufor2
!   use interp
!   implicit none

!   ! local vars
!   integer :: i, j
!   integer :: n1, n2, nn1, nn2
!   real*8, dimension(:), allocatable :: x, y, xx, yy, s1, s2, ss1, ss2
!   real*8, dimension(:,:), allocatable :: z, zz
!   real*8 :: xmin, xmax, ymin, ymax
!   real*8, parameter :: PI = 4.D0*DATAN(1.D0)

!   xmin = -4.0d0; xmax = 4.0d0; ymin = -4.0d0; ymax = 4.0d0
 
!   n1 = 50; n2 = 50
!   allocate(x(n1), y(n2), s1(n1), s2(n2), z(n1, n2))

!   s1 = cos(PI * (/ ( dble(i) , i = 0, (n1 - 1)) /) / dble(n1 - 1) )
!   s1 = s1(n1:1:-1)
!   s2 = cos(PI * (/ ( dble(i) , i = 0, (n2 - 1)) /) / dble(n2 - 1) )
!   s2 = s2(n2:1:-1)
  
!   x = (/ (xmin + (xmax - xmin) * (s1(j) + 1.0d0) / 2.0d0, j = 1, n1) /)
!   y = (/ (ymin + (ymax - ymin) * (s2(j) + 1.0d0) / 2.0d0, j = 1, n2) /)

!   ! x = (/ (xmin + (xmax - xmin) / dble(n1 - 1) * dble(j - 1), j = 1, n1) /)
!   ! y = (/ (ymin + (ymax - ymin) / dble(n2 - 1) * dble(j - 1), j = 1, n2) /)

!   do i = 1, n1
!      do j = 1, n2
!         z(i,j) = f(x(i), y(j))   
!      end do
!   end do

!   nn1 = 100; nn2 = 100
!   allocate(xx(nn1), yy(nn2), zz(nn1, nn2))
!   xx = (/ (xmin + (xmax - xmin) / dble(nn1 - 1) * dble(j - 1), j = 1, nn1) /)
!   yy = (/ (ymin + (ymax - ymin) / dble(nn2 - 1) * dble(j - 1), j = 1, nn2) /)

!   do i = 1, nn1
!      do j = 1, nn2
!         call lag2d(x, y, z, xx(i), yy(j), zz(i,j))
!      end do
!   end do

!   call surf(x, y, z, contour = 'xy')
!   call surf(xx, yy, zz, contour = 'xy')




!   ! done here
! contains

!   function f(x,y)
!     implicit none
!     real*8, intent(in) :: x,y
!     real*8 :: f

!     f = cos(2.0d0 * x**2.0d0 - 2.0d0 * x + 1.0d0) + sin(3.0d0 * y )

!   end function f


! end program test_interp
