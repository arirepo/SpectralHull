module spline
  ! implements variants of third-order spline 
  ! for interpolation, 1st and 2nd differentiation 

  use globals
  ! use utils
  ! use diff ! optional just for tester subroutine
  ! use gnufor2 ! optional just for tester subroutine

  implicit none

  private


  public :: spline_eval
  ! public :: spline_tester
  public :: comp_spline_weights
  public :: spline_eval2

contains
  
  !---------------SCALAR VERSION------------------------
  ! from 
  ! <http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm>
  ! with modifications (see below)
  subroutine scalar_solve_tridiag(a,b,c,d,x,n, cp, dp)
    implicit none
    !        a - sub-diagonal (means it is the diagonal below the main diagonal)
    !        b - the main diagonal
    !        c - sup-diagonal (means it is the diagonal above the main diagonal)
    !        d - right part d(z,n) <for n nodes with z number of eqs>
    !        x - the answer x(z,n)
    !        n - number of equations in the tridiagonal system
    !        cp(n), dp(z,n) are temporary real arrays that are fed into
    !        the subroutine to prevent dynamic allocation and increase
    !        the performance.
    integer,intent(in) :: n
    real(rk),dimension(n),intent(in) :: a,b,c,d
    real(rk),dimension(n),intent(out) :: x
    real(rk),dimension(n) :: cp, dp
    real(rk) :: m
    integer i

    ! initialize c-prime and d-prime
    cp(1) = c(1)/b(1) !c-prime
    dp(1) = d(1)/b(1) !d-prime
    ! solve for vectors c-prime and d-prime
    do i = 2,n
       m = b(i)-cp(i-1)*a(i)
       cp(i) = c(i)/m
       dp(i) = (d(i)-dp(i-1)*a(i))/m
    enddo
    ! initialize x
    x(n) = dp(n)
    ! solve for x from the vectors c-prime and d-prime
    do i = n-1, 1, -1
       x(i) = dp(i)-cp(i)*x(i+1)
    end do

  end subroutine scalar_solve_tridiag
  !-----------------------------------------------------

  ! computes the spline weight M for nonequally spaced distribution 
  ! y(x_i)
  subroutine comp_spline_weights(y, x , M, a, b, c, d, cp, dp, mode)
    implicit none

    ! The original data y on 
    ! uneven spaced grid x
    real(rk), dimension(:), intent(in) :: y, x
    ! the spline weights
    real(rk), dimension(:), intent(out) :: M
    ! temp buffers
    real(rk), dimension(:), intent(inout) :: a,b,c,d, cp, dp
    character(len = *), intent(in) :: mode ! the type of spline

    ! local vars
    integer :: i, n
    real(rk) :: hi, him1
    ! determing the size of input data
    n = size(y)

    ! select the type of spline 
    select case (mode)
    case ('natural')
       ! BC 1
       a(1) = 0._rk
       b(1) = 1._rk
       c(1) = 0._rk
       d(1) = 0._rk
       ! BC 2
       a(n) = 0._rk
       b(n) = 1._rk
       c(n) = 0._rk
       d(n) = 0._rk
    case ('parabolic')
       ! BC 1
       a(1) = 0._rk
       b(1) = 1._rk
       c(1) = -1._rk
       d(1) = 0._rk
       ! BC 2
       a(n) = -1._rk
       b(n) = 1._rk
       c(n) = 0._rk
       d(n) = 0._rk

       ! =============================
       ! add more modes here in future 
       ! =============================
    case default
       print *, 'fatal: no recognized mode for computing spline weights!'
       stop
    end select

    ! fill out a, b, c and d for tridiagonal system for
    ! spline weights M(i)
    do i = 2, (n-1)
       ! compute non-even spacing 
       hi = x(i+1) - x(i)
       him1 = x(i) - x(i-1)
       ! fill a, b, c and d
       a(i) = 1._rk / 6._rk * him1
       b(i) = 1._rk / 3._rk * (him1 + hi)
       c(i) = 1._rk / 6._rk * hi
       d(i) = (y(i+1) - y(i)) / hi + (y(i-1) - y(i)) / him1
    end do

    call scalar_solve_tridiag(a,b,c,d,M, n, cp, dp)
    ! Now the array <M> contains the weights

    ! done here
  end subroutine comp_spline_weights
  
  ! evaluates the interpolation, 1st or 2nd derivatives
  ! of the given non-equally spaced data y(x) at locations
  ! provided in array xx.
  ! NOTE FOR DOCUMENTATION :
  ! please see the notes for spline in ../doc folder
  ! which uses same notation and equations.  
  subroutine spline_eval(x, y, xx, yy, opt, mode)
    implicit none

    real(rk), dimension(:) , intent(in) :: x,y,xx
    real(rk), dimension(:) , intent(out) :: yy
    ! opt is the operation type : 'interp' 'diff1' 'diff2'
    ! mode is spline mode 'natural' 'parabolic' etc. 
    character(len = *), intent(in) :: opt, mode 

    ! local vars
    integer :: n, i, j, nxx
    real(rk), dimension(:), allocatable :: a,b,c,d,M,cp,dp
    real(rk) :: hi, ai, bi, ci, di

    ! acquiring the size of two vectors
    n = size(y)
    nxx = size(xx)

    ! allocating all allocatable arrays
    allocate(a(n), b(n), c(n), d(n), M(n), cp(n), dp(n)) 

    ! computing the weights regardless of the operation type
    call comp_spline_weights(y, x , M, a, b, c, d, cp, dp, mode)

    ! evalute the spline operation of y(x) on yy(xx)
    do j = 1, nxx
       ! first find the location of xx(j) 
       !in x and that's the "i" location.
       i = find_i(x, xx(j))

       ! compute the spline coefficients
       hi = x(i+1) - x(i)
       bi = M(i) / 2._rk
       ai = (M(i+1) - M(i)) / (6._rk * hi)
       ci = (y(i+1) - y(i)) / hi - M(i)*hi/2._rk + (M(i) - M(i+1))/ 6._rk * hi
       di = y(i)
       ! compute the value of  interpolation, 
       ! or differentiation based on "opt"
       select case (opt) ! type of operation

       case ('interp') ! interpolation
          yy(j) = ai * (xx(j) - x(i))**3._rk &
                + bi * (xx(j) - x(i))**2._rk &
                + ci * (xx(j) - x(i)) + di

       case ('diff1') ! first derivative
          yy(j) = 3._rk * ai * (xx(j) - x(i))**2._rk &
                + 2._rk * bi * (xx(j) - x(i)) &
                + ci

       case ('diff2') ! second derivative
          yy(j) = 6._rk * ai * (xx(j) - x(i)) &
                + 2._rk * bi

       case default ! not specified, would be dangerous
          print *, 'fatal : no recognized option is provided to spline_eval routine'
          stop 
       end select

    end do  !end of spline computations for subgrid points

    ! lil clean-up
    deallocate(a, b, c, d, M, cp, dp) 

    ! done evaluation of spline
  end subroutine spline_eval
  
  !finds the given location xxk in <x> array
  !such that x(find_i) <= xxk <= x(find_i+1) 
  ! the returned result contains find_i   
  function find_i(x, xxk)
    implicit none
    real(rk), dimension(:), intent(in) :: x
    real(rk), intent(in) :: xxk
    integer :: find_i

    ! local vars
    integer :: n, loc, j 
    logical :: found_flag

    found_flag = .false. ! not initially found!
    n = size(x)

    do j = 1, (n-1) !brute force search <NOT GOOD!>
       if( (xxk >= (x(j) - 1.0d-15) ) .and. (xxk <= (x(j+1) + 1.0d-15)) ) then
          loc = j ! this will be returned!
          found_flag = .true.
          exit
       else
          continue
       end if
    end do
    ! decide ...
    if(found_flag) then ! we are fine
       find_i = loc 
    else
       print *, 'fatal : the <xx> data not in range of <x>. error is spline! exit'
       print *, 'x = ', x
       print *, 'xxk = ', xxk
       stop
    end if
    ! done this fuction
  end function find_i

  ! assumes weights are already computed/stored
  subroutine spline_eval2(x, y, xx, yy, opt, M)
    implicit none

    real(rk), dimension(:) , intent(in) :: x,y
    real(rk), intent(in) :: xx
    real(rk), intent(out) :: yy
    ! opt is the operation type : 'interp' 'diff1' 'diff2'
    character(len = *), intent(in) :: opt 
    real(rk), dimension(:), intent(in) :: M

    ! local vars
    integer :: n, i
    real(rk) :: hi, ai, bi, ci, di

    ! acquiring the size of two vectors
    n = size(y)


    ! first find the location of xx 
    !in x and that's the "i" location.
    i = find_i(x, xx)

    ! compute the spline coefficients
    hi = x(i+1) - x(i)
    bi = M(i) / 2._rk
    ai = (M(i+1) - M(i)) / (6._rk * hi)
    ci = (y(i+1) - y(i)) / hi - M(i)*hi/2._rk + (M(i) - M(i+1))/ 6._rk * hi
    di = y(i)
    ! compute the value of  interpolation, 
    ! or differentiation based on "opt"
    select case (opt) ! type of operation

    case ('interp') ! interpolation
       yy = ai * (xx - x(i))**3._rk &
            + bi * (xx - x(i))**2._rk &
            + ci * (xx - x(i)) + di

    case ('diff1') ! first derivative
       yy = 3._rk * ai * (xx - x(i))**2._rk &
            + 2._rk * bi * (xx - x(i)) &
            + ci

    case ('diff2') ! second derivative
       yy = 6._rk * ai * (xx - x(i)) &
            + 2._rk * bi

    case default ! not specified, would be dangerous
       print *, 'fatal : no recognized option is provided to spline_eval routine'
       stop 
    end select

    ! done evaluation of spline
  end subroutine spline_eval2

  ! ! This is a tester subroutine for spline
  ! ! just use it to make sure that
  ! ! spline implementation is correct  
  ! subroutine spline_tester(nn)
  !   implicit none

  !   ! number of tester collocation points
  !   integer, intent(in) :: nn 

  !   ! local vars for the lil tester script
  !   real(rk), dimension(:), allocatable :: x0, xi0, f, df, df_spline, df_FD, tmp, df_cheb
  !   real(rk), dimension(:,:), allocatable :: DD
  !   integer :: j

  !   ! local allocation
  !   allocate(x0(nn), xi0(nn), f(nn), df(nn), df_spline(nn), df_FD(nn), tmp(nn), DD(nn,nn), df_cheb(nn))

  !   ! fill DD and xi0
  !   call cheby_diff(nn-1, 'matlab', DD, xi0)

  !   ! initializing the Chebyshev and equally spaced meshes ...
  !   ! NOTE: both intervals are set on [-1,1].
  !   x0 = (/ (real(j, rk), j = 0,(nn-1)) /)
  !   if (any(xi0 /= cos( x0 * PI / real(nn-1,rk)  ))) then
  !      print *, 'fatal : the generated Chebyshev grid is ' & 
  !           , 'not consistent with cheby_diff() subroutine!'
  !      stop
  !   end if
  !   x0 = 2._rk * x0 / real(nn-1, rk) - 1._rk

  !   ! sample function for diff. and interp.
  !   f = x0**(3._rk) ! the function itself
  !   df = 3._rk * (x0**2._rk) !its analytical der.

  !   ! finding the transformed derivative using 2nd-order central FD
  !   ! note assuming arbitrary Dxi like Dxi = 1._rk
  !   call test_diff(f, df_FD, 1._rk )
  !   call test_diff(x0, tmp, 1._rk )
  !   df_FD = df_FD / tmp

  !   ! finding the transformed derivative using Chebyshev differentiation matrix
  !   df_cheb = matmul(DD,f)/  matmul(DD,x0) 
  !   ! the following is semi-analytical
  !   ! df_cheb = matmul(DD, f) / (-2._rk/(PI*sqrt(1._rk - xi0**2._rk)))

  !   ! finding the transformed derivative using spline
  !   xi0 = xi0(nn:1:-1) ! first rotate the Cheby grid
  !   call spline_eval(xi0, f, xi0, df_spline, 'diff1', 'parabolic')
  !   call spline_eval(xi0, x0, xi0, tmp, 'diff1', 'parabolic')
  !   df_spline = df_spline / tmp
  !   xi0 = xi0(nn:1:-1) ! rotate it back to original state

  !   ! ! plot the results
  !   ! call plot(x0, df, x0, df_FD) ! finite difference
  !   ! call plot(x0, df, x0, df_cheb) ! Chebyshev
  !   ! call plot(x0, df, x0, df_spline) ! spline

  !   ! plot the interpolation of function on Chebyshev grid
  !   call spline_eval(x0, f, xi0, tmp, 'interp', 'parabolic')
  !   ! call plot(x0, f, xi0, tmp) ! spline interpolation
  !   ! call plot(x0, df, xi0, matmul(DD, tmp)) ! Chebyshev diff on interpolated data

  !   ! local deallocation
  !   deallocate(x0, xi0, f, df, df_spline, df_FD, tmp, DD, df_cheb)
    
  !   ! subroutine ends up here!

  ! contains ! local subroutines 

  !   ! sample FD central diff with upwind schemes near boundaries
  !   ! for scalar case
  !   subroutine test_diff(f, df, dx)
  !     implicit none
  !     real(rk), dimension(:), intent(in) :: f
  !     real(rk), dimension(:), intent(out) :: df
  !     real(rk), intent(in) :: dx

  !     ! local vars
  !     integer :: i, n

  !     n = size(f)
  !     df(1) = (-1.5_rk * f(1) + 2.0_rk * f(2) - 0.5_rk * f(3)) / dx
  !     do i = 2, (n-1)
  !        df(i) = (f(i+1) - f(i-1)) / (2._rk * dx)
  !     end do
  !     df(n) = (1.5_rk * f(n) - 2.0_rk * f(n-1) + 0.5_rk * f(n-2)) / dx

  !     ! done
  !   end subroutine test_diff

  ! end subroutine spline_tester

  ! done this module
end module spline
