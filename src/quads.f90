module quads
  ! implements tables of quadratures
  ! for FEM integration
  implicit none

  private

  public :: p2n, get_quad

contains

  ! converts the polynomial order of the 
  ! integrand "p" to the required number
  ! of Gauss points "n" needed for quadrature
  ! to be exact. This depends on element
  ! type also.

  subroutine p2n(etype, p, n)
    implicit none
    integer, intent(in) :: etype ! element type
    integer, intent(in) :: p ! order of the integrand
    integer, intent(out) :: n ! number of Gauss points 

    ! local vars
    integer :: nx 

    select case (etype) ! what element ?

    case (1) ! triangle

       if(p .eq. 1) then 
          n = 1
       elseif( p .eq. 2 ) then 
          n = 3
       elseif ( p .eq. 3) then
          n = 4
       else
          print *, 'p > 3 is not implemented for Gauss points of a triangle' &
               , ' at this version. stop'
          stop
       end if

    case (2) ! quad 

       nx = ceiling( dble(p+3)/(2.0D0) )
       n = nx * nx

    case default

       print *,' unable to recognize element type in p2n(...). stop.'
       stop

    end select

    ! done here
  end subroutine p2n

  
  ! generates the location (r_i,s_i) and weights W_i
  ! of Jacobi-Gauss-Legendre quadrature 
  ! corresponding to the master element for numerically
  ! integrating surface integral on triangles and quadrilaterals.
  
  subroutine get_quad(etype, n, alpha, beta, r, s, W)
    implicit none
    integer, intent(in) :: etype ! element type
    integer, intent(in) :: n ! number of Gauss points
    real*8, intent(in) :: alpha, beta ! coeff. of Jacobi polynomials  
    ! the coords of Gauss point (r,s) in master elem and 
    ! the corresponding Weight at that point
    real*8, dimension(n), intent(out) :: r,s, W 


    ! local vars
    integer :: nx
    real*8, dimension(:), allocatable :: xi, wi, derjac
    integer :: i, i1, i2


    select case (etype) ! what element?

    case(1) ! a triangle master element

       select case (n) 

       case (1) ! one Gauss point.

          r(1) = 1.0d0 / 3.0d0; s(1) = 1.0d0 / 3.0d0; W(1) = 1.0d0 / 2.0d0

       case (3) ! three Gauss points.

          r(1) = 1.0d0 / 6.0d0; s(1) = 1.0d0 / 6.0d0; W(1) = 1.0d0 / 6.0d0
          r(2) = 2.0d0 / 3.0d0; s(2) = 1.0d0 / 6.0d0; W(2) = 1.0d0 / 6.0d0
          r(3) = 1.0d0 / 6.0d0; s(3) = 2.0d0 / 3.0d0; W(3) = 1.0d0 / 6.0d0

       case (4) ! four Gauss points.

          r(1) = 1.0d0 / 3.0d0; s(1) = 1.0d0 / 3.0d0; W(1) = -9.0d0 / 32.0d0
          r(2) = 3.0d0 / 5.0d0; s(2) = 1.0d0 / 5.0d0; W(2) = 25.0d0 / 96.0d0
          r(3) = 1.0d0 / 5.0d0; s(3) = 3.0d0 / 5.0d0; W(3) = 25.0d0 / 96.0d0
          r(4) = 1.0d0 / 5.0d0; s(4) = 1.0d0 / 5.0d0; W(4) = 25.0d0 / 96.0d0

       case default ! not found

          print *, 'unknown number of Gauss points "n" for a triangle.'
          print *, 'error happened in get_quad(...). stop.'
          stop 

       end select


    case(2) ! a quad

       ! computing quadrature points in the "x" direction
       nx = int(sqrt(dble(n)))

       if ( (nx * nx) .ne. n ) then 
          print *, 'the number of Gauss points inside quad should be n = nx*nx' &
               , ' where nx is some integer. Refer to get_quad(...) subroutine' &
               , ' for more info. stopped!'
          stop
       end if


       ! computing the Jacobi-Gauss-Legendre coordinates
       ! and corresponding weight functions
       allocate( xi(nx), wi(nx), derjac(nx) )
       call ZEJAGA(nx, alpha, beta, xi, derjac)
       call WEJAGA(nx, alpha, beta, xi, derjac, wi)

       ! filling the output arrays
       do i = 1, nx

          i1 = (i-1) * nx + 1
          i2 = i * nx

          r(i1:i2) = xi
          s(i1:i2) = xi(i)
          W(i1:i2) = wi(i) * wi

       end do

       deallocate( xi, wi, derjac )


    case default ! element type not found

       print *, 'undefined element type in get_quad(...). stop.'
       stop

    end select


    ! done here
  end subroutine get_quad



end module quads
