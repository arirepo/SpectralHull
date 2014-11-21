module nurbs
  ! creates nurbs curves of order p to
  ! interpolate a given set of points x(np), y(np) 

  use globals

  implicit none

  private

  public :: parametrize_points
! arguments: real e, real x(:), real y(:), real t(:)

  public :: nurbs_eval
! arguments: real Cx(:), real Cy(:), integer np, real u,
!            real x, real y, real w(:), real knot(:), integer prime

  public :: comp_nurbs
! arguments: real x(:), real y(:), real u(:), 
!            real Cx(:), real Cy(:), real w(:), real knot(:)

  public :: print_vals
! arguments: real Cx(:), real Cy(:), real w(:), real knot(:)

  public :: tester

  contains
    subroutine nurbs_eval(Cx, Cy, np, u, x, y, w, knot, prime)
      real(rk), dimension(:), intent(in) :: Cx
      real(rk), dimension(:), intent(in) :: Cy
      integer,                intent(in) :: np
      real(rk),               intent(in) :: u
      real(rk), dimension(:), intent(in) :: w
      real(rk), dimension(:), intent(in) :: knot
      real(rk),               intent(out) :: x, y
      integer,                intent(in) :: prime

      real(rk), dimension(:), allocatable :: N
      real(rk) :: denom
      integer :: i, p, nc

      nc = size(Cx)
      if((size(Cy) < nc) .or. (size(Cy) > nc))then
        write(*,*)'Error in nurbs routine nurbs_eval: ',&
                & 'size of control points vectors are different'
      end if
      p = size(knot) - nc - 1

      allocate( N(nc) )

      call evaluate_bases(p, N, knot, u, prime)
      denom = dot_product(w, N)

      x = 0.0_rk
      y = 0.0_rk
      do i = 1, nc
        x = x + N(i) * w(i) * Cx(i)
        y = y + N(i) * w(i) * Cy(i)
      end do
      if( abs(x) > 0 )then
        x = x / denom
      end if
      if( abs(y) > 0 )then
        y = y / denom
      end if

      if(prime > 0)then
        x = 1.0_rk
      end if
      deallocate(N)
    end subroutine nurbs_eval

    ! evaluates nc b-spline basis functions of order p at point u.
    ! the results are stored in N(1:nc)
    subroutine evaluate_bases(p, N, knot, u, prime)
      integer, intent(in) :: p
      real(rk), dimension(:), intent(in out) :: N
      real(rk), dimension(:), intent(in) :: knot
      real(rk) :: u
      integer, intent(in) :: prime

      integer :: i, j, k, nc
      real(rk) :: v1, v2, d1, d2
      real(rk), dimension(:), allocatable :: N1

      nc = size(N)
      if(u > 1.0_rk .or. u < 0.0_rk)then
        write(*,*)'Error in subroutine evaluate_bases: ',&
                & 'input value u not within range [0, 1]'
      else if(u < knot(nc+p+1))then
        allocate( N1(nc + 1) )
        if(prime < 1)then
          N1(:) = 0.0_rk
          do i = 1, nc
            if( (u >= knot(i)) .and. (u < knot(i + 1)) )then
              N1(i) = 1.0_rk
            end if
          end do
          do k = 1, p
            do i = 1, nc
              v1 = 0.0_rk
              v2 = 0.0_rk
              d1 = knot(i + k) - knot(i)
              d2 = knot(i + k + 1) - knot(i + 1)
              if( (abs(N1(i + 1)) > 0) .and. (abs(knot(i + k + 1) - u) > 0) )then
                v2 = N1(i + 1) * ( knot(i + k + 1) - u ) / d2
              end if
              if( (abs(N1(i)) > 0) .and. (abs(u - knot(i)) > 0) )then
                v1 = N1(i) * ( u - knot(i) ) / d1
              end if
              N(i) = v1 + v2
            end do
            N1(:nc) = N(:nc)
          end do
        else 

        end if
        deallocate( N1 )
      else if(.not. (u < knot(nc + p + 1)))then
        N(:) = 0.0_rk
        N(nc) = 1.0_rk
      end if
    end subroutine evaluate_bases  

    ! this subroutine creates a nurbs curve to approximate a given set of points x, y. 
    ! The resulting curve is described by Cx, Cy, w, knot.
    subroutine comp_nurbs(x, y, u, Cx, Cy, w, knot)
      real(rk), dimension(:), intent(in out) :: x
      real(rk), dimension(:), intent(in out) :: y
      real(rk), dimension(:), intent(in) :: u
      real(rk), dimension(:), intent(in out) :: Cx
      real(rk), dimension(:), intent(in out) :: Cy
      real(rk), dimension(:), intent(in out) :: w
      real(rk), dimension(:), intent(in out) :: knot

      real(rk) :: D, nx, ny
      integer :: nc, np, p
      integer :: i

      nc = size(Cx)
      np = size(x)
      p = size(knot) - nc - 1
      w(:) = 1.0_rk

      ! create knot vector
      call generate_knots(p, nc, knot, u)
      ! solve a system of equations for simplified nurbs curve with all weights w = 1
      call solve_system(x, y, Cx, Cy, knot, u)

      D = 0.0;
      do i = 1, np
        call nurbs_eval(Cx, Cy, np, u(i), nx, ny, w, knot, 0)
        D = MAX(D, sqrt((nx - x(i))**2 + (ny - y(i))**2))
        x(i) = nx; y(i) = ny
      end do
      write(*,*)
      write(*,*)'maximum error between boundary points and computed nurbs curve:', D
      write(*,*)
    end subroutine comp_nurbs

    ! given np long vector u, 0 <= u_i <= 1 generate a knot vector 
    ! nc + p long for creating b-spline basis functions.
    subroutine generate_knots(p, nc, knot, u)
      integer, intent(in) :: p
      integer, intent(in) :: nc
      real(rk), dimension(:), intent(in out) :: knot
      real(rk), dimension(:), intent(in) :: u

      integer :: np, i, j

      np = size(u)

      knot(1:p + nc + 1) = 0.0_rk
      do i = 2, nc - p
        do j = i, i + p - 1
          knot(i + p) = knot(i + p) + u(j)
        end do
        knot(i + p) = knot(i + p) / real(p)
      end do
      knot(nc + 1:p + nc + 1) = 1.0_rk
    end subroutine generate_knots

    ! solve a system of equations X = N C for C where X are geometry points, N are
    ! basis functions evaluated at points X and C are control points.
    ! The result is a simplified NURBS curve that assumes that all weights w are 1 and 
    ! interpolates the points X
    subroutine solve_system(x, y, Cx, Cy, knot, t)
      real(rk), dimension(:), intent(in) :: x
      real(rk), dimension(:), intent(in) :: y
      real(rk), dimension(:), intent(in out) :: Cx
      real(rk), dimension(:), intent(in out) :: Cy
      real(rk), dimension(:), intent(in) :: knot
      real(rk), dimension(:), intent(in) :: t

      integer :: np, lda, ldb, info_x, info_y
      real(rk), dimension(:,:), allocatable :: N
      real(rk), dimension(:),   allocatable :: N1
      integer,  dimension(:,:), allocatable :: IPIV
      integer :: i, p

      np = size(x)
      if( abs(size(y) - np) > 0 )then
        write(*,*)'Error in nurbs solve_system routine: ',&
                & 'x, y arrays of different dimensions'
        stop
      else if( abs(size(t) - np) > 0 )then
        write(*,*)'Error in nurbs solve_system routine: ',&
                & 'x, t arrays of different dimensions'
        stop
      else if( abs(size(Cy) - np) > 0 )then
        write(*,*)'Error in nurbs solve_system routine: ',&
                & 'x, Cy arrays of different dimensions'
        stop
      else if( abs(size(Cx) - np) > 0 )then
        write(*,*)'Error in nurbs solve_system routine: ',&
                & 'x, Cx arrays of different dimensions'
        stop
      end if

      p = size(knot) - np - 1
      lda = np
      ldb = np

      allocate(N(np,np))
      allocate(N1(np))
      allocate(IPIV(np,np))

      do i = 1, np
        call evaluate_bases(p, N1, knot, t(i), 0)
        N(i,:) = N1(:)
      end do
      Cx(:) = x(:)
      call dgesv(np, 1, N, lda, IPIV, Cx, ldb, info_x)

      do i = 1, np
        call evaluate_bases(p, N1, knot, t(i), 0)
        N(i,:) = N1(:)
      end do
      Cy(:) = y(:)
      call dgesv(np, 1, N, lda, IPIV, Cy, ldb, info_y)

      if(abs(info_x) > 0 .or. abs(info_y) > 0)then
        write(*,*)'Error in subroutine solve_system: linear solve failed'
!       stop
      end if

      deallocate(N)
      deallocate(N1)
      deallocate(IPIV)

!     real(rk), dimension(:,:), allocatable :: N  ! matrix of nc basis functions evaluated at np points
!     real(rk), dimension(:), allocatable :: N1   ! nc basis functions evaluated at a given point 
!     real(rk), dimension(:,:), allocatable :: U  ! SVD matrix
!     real(rk), dimension(:,:), allocatable :: DI ! diagonal SVD matrix
!     real(rk), dimension(:), allocatable :: D    ! vector of diagonals DI
!     real(rk), dimension(:,:), allocatable :: V  ! SVD matrix
!     real(rk), dimension(:), allocatable :: work ! lapack DGESVD vector

!     integer :: np, nc, p, i
!     integer :: lda, ldu, ldvt, lwork, info

!     character :: jobu = 'A'
!     character :: jobvt = 'A'

!     np = size(x)
!     if( (size(y) > np) .or. (size(y) < np) )then
!       write(*,*)'Error in nurbs solve_system routine: x, y arrays of different dimensions'
!     else if( (size(t) > np) .or. (size(t) < np) )then
!       write(*,*)'Error in nurbs solve_system routine: x, t arrays of different dimensions'
!     end if
!     nc = size(Cx)
!     if( (size(Cy) > nc) .or. (size(Cy) < nc) )then
!       write(*,*)'Error in nurbs solve_system routine: arrays of different dimensions'
!     end if
!     p = size(knot) - nc - 1
!     lda = np
!     ldu = np
!     ldvt = nc
!     lwork = MAX(1, 3 * MIN(nc, np) + MAX(nc, np), 5 * MIN(nc, np))

!     allocate( N(np, nc) )
!     allocate( N1(nc) )
!     allocate( U(np, np) )
!     allocate( DI(nc, np) )
!     allocate( D( MIN(np, nc) ) )
!     allocate( V(nc, nc) )
!     allocate( work( lwork ) )

!     do i = 1, np
!       call evaluate_bases(p, N1, knot, t(i), 0)
!       N(i,:) = N1(:)
!     end do

!     !http://www.netlib.org/templates/single/GMRES.f
!     call DGESVD(jobu, jobvt, np, nc, N, lda, D, U, ldu, V, ldvt, work, lwork, info)

!     if(abs(info) > 0) then
!       write(*,*)
!       write(*,*)'SVD Error: Lapack failed to decompose matrix for nurbs point fitting'
!       write(*,*)
!     end if

!     DI(:,:) = 0.0_rk
!     do i = 1, MIN(nc, np)
!       if( D(i) > 0 )then
!         DI(i,i) = 1.0_rk / D(i)
!       else
!         DI(i,i) = 0.0_rk
!       end if
!     end do

!     Cx = matmul(matmul(matmul(trans(V, nc, nc), DI), trans(U, np, np)), x)
!     Cy = matmul(matmul(matmul(trans(V, nc, nc), DI), trans(U, np, np)), y)

!     Cx(1) = x(1); Cx(nc) = x(np)
!     Cy(1) = y(1); Cy(nc) = y(np)
!     deallocate( N )
!     deallocate( N1 )
!     deallocate( U )
!     deallocate( DI )
!     deallocate( D )
!     deallocate( V )
!     deallocate( work )
    end subroutine solve_system

    ! returns the transpose of a mxn matrix Mat
    function trans(Mat, m, n)
      real(rk), dimension(:,:), intent(in) :: Mat
      integer, intent(in) :: m
      integer, intent(in) :: n
      real(rk), dimension(n, m) :: trans
      integer :: i
      do i = 1, m
        trans(:n, i) = Mat(i, :n)
      end do
    end function trans

    ! prints the points x, y for a NURBS curve with
    ! control points Cx, Cy, weights w, and knot vector knot
    ! evaluated at 0 <= t <= 1.0 
    subroutine print_vals(num, Cx, Cy, w, knot)
      integer, intent(in) :: num
      real(rk), dimension(:), intent(in) :: Cx
      real(rk), dimension(:), intent(in) :: Cy
      real(rk), dimension(:), intent(in) :: w
      real(rk), dimension(:), intent(in) :: knot

      real(rk), dimension(:), allocatable :: N
      real(rk) :: xv, denom, xnum, ynum, v1, v2
      integer :: i, j, nc, p, fnum

      character(len=1024) :: fname
      character(len=1024) :: format_string

      format_string = "(A12, I0, A4)"
      write(fname, format_string) "nurbs_points", num, ".dat"

      nc = size(Cx)
      p = size(knot) - nc - 1

      allocate(N(nc))
      fnum = 30 + num
      open(fnum,file=trim(fname))
      xv = 0.0_rk
      do i = 1, 1000
        call evaluate_bases(p, N, knot, xv, 0)
        denom = dot_product(w, N)
        xnum = 0.0_rk
        ynum = 0.0_rk
        do j = 1, nc
          xnum = xnum + N(j) * w(j) * Cx(j)
          ynum = ynum + N(j) * w(j) * Cy(j)
        end do
        if(abs(xnum) > 0)then
          xnum = xnum / denom
        end if
        if(abs(ynum) > 0)then
          ynum = ynum / denom
        end if
        write(fnum,*) xnum, ynum
        xv = xv + 0.001_rk
      end do
      close(fnum)
      deallocate(N)
      open(18,file='control_points.dat')
      do i = 1, nc
        write(18,*)Cx(i),Cy(i)
      end do
      close(18)
    end subroutine print_vals


    ! given np points x, y, create a vector t np long where, 0 <= t_i <= 1
    subroutine parametrize_points(e, x, y, t)
      real(rk) :: e
      real(rk), dimension(:), intent(in) :: x
      real(rk), dimension(:), intent(in) :: y
      real(rk), dimension(:), intent(in out) :: t

      integer :: i, j, np
      np = size(x)

      t(:) = 0.0_rk
      do i = 1, np
        do j = 2, i
          t(i) = t(i) + ( (x(j) - x(j - 1))**2 + &
                        & (y(j) - y(j - 1))**2 )**(0.5 * e)
        end do
      end do
      t = t / t(np)
    end subroutine parametrize_points


    subroutine tester()
      integer :: np 
      integer :: nc
      integer :: p = 3
      real(rk), dimension(:), allocatable :: x
      real(rk), dimension(:), allocatable :: y
      real(rk), dimension(:), allocatable :: Cx
      real(rk), dimension(:), allocatable :: Cy
      real(rk), dimension(:), allocatable :: t
      real(rk), dimension(:), allocatable :: w
      real(rk), dimension(:), allocatable :: knot
      real(rk) :: e = 0.5_rk, nnx, nny
      integer :: nseg, ssize, i, rem, j, k
      logical :: found
      real(rk), dimension(:), allocatable :: nx
      real(rk), dimension(:), allocatable :: ny
      real(rk), dimension(:), allocatable :: nCx
      real(rk), dimension(:), allocatable :: nCy
      real(rk), dimension(:), allocatable :: nt
      real(rk), dimension(:), allocatable :: nw
      real(rk), dimension(:), allocatable :: nknot
      real(rk), dimension(:), allocatable :: tmp

      write(*,*)'enter number of points to read'
      read(*,*)np
      nc = np
      if((np - p) < 1)then
        write(*,*)'Error: must have at least p + 1 points'
        write(*,*)'changing p from',p,'to',p - 1
        p = p - 1
      end if

      allocate(x(np))
      allocate(y(np))
      allocate(Cx(nc))
      allocate(Cy(nc))
      allocate(t(np))
      allocate(w(nc))
      allocate(knot(nc + p + 1))
      call read_file(x, y)

      call parametrize_points(e, x, y, t)
      call comp_nurbs(x, y, t, Cx, Cy, w, knot)
      call print_vals(1, Cx, Cy, w, knot)

      deallocate(x)
      deallocate(y)
      deallocate(Cx)
      deallocate(Cy)
      deallocate(t)
      deallocate(w)
      deallocate(knot)
    end subroutine tester

    subroutine read_file(x, y)
      real(rk), dimension(:), intent(in out) :: x
      real(rk), dimension(:), intent(in out) :: y
      integer :: np, i

      np = size(x)
      open(100, file='parabolapoints')
      do i = 1, np
        read(100,*) x(i), y(i)
      end do
      close(100)
    end subroutine read_file
end module nurbs
