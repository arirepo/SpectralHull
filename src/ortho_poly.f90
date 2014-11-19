module ortho_poly
  ! implements all sorts of crazy things
  ! with orthogonal polynomials
  !
  ! ATTENTION : 
  ! NEEDS LAPACK AND BLAS TO LINK
  !

  use globals

  implicit none

  private

  public :: comp_jacobi_diff_matrix, comp_ortho_diff_matrix
  public :: comp_jacobi_integ_matrix
  public :: comp_gauss_integ_matrix, test_gauss_integ_matrix

contains

  ! computes the differentiation matrix comming
  ! from general Jacobi polynomials P_n^{alpha, beta}(x)
  ! and also the corresponding point distribution
  ! on the given open interval x = (x0 x1)

  subroutine comp_jacobi_diff_matrix(alpha, beta, n, D, DD, x)
    implicit none
    real(rk), intent(in) :: alpha, beta
    integer , intent(in) :: n
    real(rk), dimension(n+1, n+1), intent(out) :: D, DD
    real(rk), dimension(n+1), intent(out) :: x

    ! local vars
    real(rk), dimension(:,:), allocatable :: uixk, duixk, Psi, dduixk 
    integer i, k

    ! LAPACK vars
    real(rk), dimension(n+1) :: work
    integer , dimension(n+1) :: ipiv
    integer nnn, info

    ! external subroutines
    external DGETRF
    external DGETRI

    ! local allocations
    allocate( uixk(n+1, n+1), duixk(n+1, n+1) &
            , Psi(n+1, n+1), dduixk(n+1, n+1) )

    ! compute point distribution (relatives)
    call ZEJAGL(n, alpha, beta, x, uixk(:,1))       
    
    ! check to see the jacobi polynomial is not zero 
    ! at any points among the relative points
    if( any(uixk(:,1) == 0.0_rk) ) then
       print *, ' in comp_jacobi_diff_matrix(...)'
       print *, 'uixk(:,1) = 0 somewhere!!! exit.'
       stop
    end if
    ! compute "uixk" and "duixk" matrices by evaluating
    ! the value and the derivative of Jacobi polynomials
    ! for orders i = 0:n at xk collocation points
    ! given by x-distribution for k = 1:(n+1)
    ! the results should be (n+1)x(n+1) matrices.
    do i = 0 , n
       do k = 1 , (n+1)
          call VAJAPO(i, alpha, beta, x(k), uixk(k,i+1) &
                       , duixk(k,i+1), dduixk(k,i+1))
       end do
    end do
    ! initialize Psi with uixk where it will be overwritten
    ! by LU factors computed by LAPACK 
    Psi = uixk
    nnn = n+1 ! size n in LAPACK
    ! do LU factorization using LAPACK
    call DGETRF(nnn, nnn, Psi, nnn, ipiv, info)
    if( info /= 0 ) then 
       print *,  'uixk is numerically singular or ' &
              ,  'LU can not be computed!'
       stop
    end if
    ! compute inverse using already computed LU
    call DGETRI(nnn, Psi, nnn, ipiv, work, nnn, info)
    if( info /= 0 ) then 
       stop 'inversion of <uixk> failed!'
    end if

    ! computing the final Jacobi 1st differentiation matrix
    D  = matmul( duixk, Psi)
    ! computing the final Jacobi 2nd differentiation matrix
    DD = matmul(dduixk, Psi)

    ! clean ups
    deallocate( uixk, duixk, Psi, dduixk)

    ! done here
  end subroutine comp_jacobi_diff_matrix

  ! computes the differentiation matrix arising
  ! from general orthogonal polynomials
  ! including:
  !
  ! polclass = "jacobi"
  ! alpha, beta = something    =>   Jacobi
  ! polclass = "legendre"
  ! alpha =  0.,  beta = 0.    =>   Legendre
  ! polclass = "cheby"
  ! alpha =-0.5, beta =-0.5    =>   Chebyshev
  ! polclass = "slaguerre"
  ! alpha = something > -1 for  scaled Laguerre
  !
  ! and also the corresponding point distribution
  ! on the given open interval x = (x0 x1)
  ! which depends on the class of the polynomials

  subroutine comp_ortho_diff_matrix(polclass, alpha, beta, n, D, DD, x)
    implicit none
    character(len=*), intent(in) :: polclass
    real(rk), intent(in) :: alpha, beta
    integer , intent(in) :: n
    real(rk), dimension(n+1, n+1), intent(out) :: D, DD
    real(rk), dimension(n+1), intent(out) :: x

    ! local vars
    real(rk), dimension(:,:), allocatable :: uixk, duixk, Psi, dduixk 
    real(rk), dimension(:), allocatable :: Lhat  
    integer i, k, nlagu

    ! LAPACK vars
    real(rk), dimension(n+1) :: work
    integer , dimension(n+1) :: ipiv
    integer nnn, info

    ! external subroutines
    external DGETRF
    external DGETRI

    ! local allocations
    allocate( uixk(n+1, n+1), duixk(n+1, n+1) &
            , Psi(n+1, n+1), dduixk(n+1, n+1), Lhat(n+1) )

    ! HARD reset !
    uixk = 0._rk; duixk = 0._rk
    Psi = 0._rk; dduixk = 0._rk
    nlagu = n + 1

    ! compute point distribution (relatives)
    select case (polclass)
    case ('jacobi')
       call ZEJAGL(n, alpha, beta, x, uixk(:,1))       
    case ('legendre')
       call ZELEGL(n, x, uixk(:,1))
    case ('cheby')
       call ZECHGL(n, x)                                           
    case ('slaguerre')
       !call ZELAGR(nlagu, alpha, x, uixk(:,1))
       call ZELAGR(nlagu, alpha, x, Lhat)
       !print *, ' x    = ', x
       ! print *, ' Lhat = ', Lhat
       ! stop
       uixk(:,1) = Lhat
    case default
       print *, 'unknown type of orthogonal polynomials in'
       print *, 'comp_ortho_diff_matrix(...) subroutine.'
       print *, 'preparing to exit ...'
       stop
    end select
    
    ! check to see if the polynomial is not zero 
    ! at any points among the relative points
    if( any(uixk(:,1) == 0.0_rk) ) then
       print *, ' in comp_ortho_diff_matrix(...)'
       print *, 'uixk(:,1) = 0 somewhere!!! exit.'
       print *, 'uixk(:,1) = ', uixk(:,1)
       stop
    end if
    ! compute "uixk" and "duixk" matrices by evaluating
    ! the value and the derivative of the polynomials
    ! for orders i = 0:n at xk collocation points
    ! given by x-distribution for k = 1:(n+1)
    ! the results should be (n+1)x(n+1) matrices.
    select case (polclass)

    case ('jacobi')

       do i = 0 , n
          do k = 1 , (n+1)
             call VAJAPO(i, alpha, beta, x(k), uixk(k,i+1) &
                  , duixk(k,i+1), dduixk(k,i+1))
          end do
       end do

    case ('legendre')

       do i = 0 , n
          do k = 1 , (n+1)
             call VALEPO(i, x(k), uixk(k,i+1) &
                  , duixk(k,i+1), dduixk(k,i+1))
          end do
       end do

    case ('cheby')

       do i = 0 , n
          do k = 1 , (n+1)
             call VACHPO(i, x(k), uixk(k,i+1) &
                  , duixk(k,i+1), dduixk(k,i+1))
          end do
       end do

    case ('slaguerre')

       do i = 0 , n
          do k = 1 , (n+1)
             call VALASF(i, alpha, x(k), uixk(k,i+1) &
                  , duixk(k,i+1))
             ! call VALAPO(i, alpha, x(k), uixk(k,i+1) &
             !      , duixk(k,i+1), dduixk(k,i+1))

          end do
       end do

    end select
    ! initialize Psi with uixk where it will be overwritten
    ! by LU factors computed by LAPACK 
    Psi = uixk
    nnn = n+1 ! size n in LAPACK
    ! do LU factorization using LAPACK
    call DGETRF(nnn, nnn, Psi, nnn, ipiv, info)
    if( info /= 0 ) then 
       print *,  'uixk is numerically singular or ' &
              ,  'LU can not be computed!'
       stop
    end if
    ! compute inverse using already computed LU
    call DGETRI(nnn, Psi, nnn, ipiv, work, nnn, info)
    if( info /= 0 ) then 
       stop 'inversion of <uixk> failed!'
    end if

    ! finalize ...    
    if ( polclass == 'slaguerre') then
       call comp_exact_D_scaled_laguerre(nlagu, alpha, Lhat, x, D)
       DD = matmul(D, D)
    else
       ! computing the final 1st differentiation matrix
       D  = matmul( duixk, Psi)
       ! computing the final 2nd differentiation matrix
       DD = matmul(dduixk, Psi)
    end if

    ! clean ups
    deallocate( uixk, duixk, Psi, dduixk, Lhat)

    ! done here
  end subroutine comp_ortho_diff_matrix

  !
  ! the following fills the scaled Laguerre
  ! first differentiation matrix by using
  ! analytical relations given in 
  ! eq. (7.5.3) page 147 of Daniele Funaro's book.
  !
  ! the second derivative DD is simply obtained
  ! by squaring D as suggested in page 147 of
  ! that book.
  subroutine comp_exact_D_scaled_laguerre(n, alpha, L, eta, D)
    implicit none
    integer, intent(in) :: n
    real(rk), intent(in) :: alpha
    real(rk), dimension(0:(n-1)), intent(in) :: L, eta
    real(rk), dimension(0:(n-1), 0:(n-1)), intent(out) :: D

    ! local vars
    integer :: i,j

    ! begin filling the matrix
    D(0,0) = - dble(n-1)/ (alpha + 2._rk)
    
    forall (i = 1:(n-1), j = 0:0)
       D(i,j) = (alpha + 1._rk) * L(i) / eta(i)
    end forall

    forall (i = 0:0, j = 1: (n-1) )
       D(i,j) = -1._rk / ((alpha + 1._rk) * eta(j) * L(j))
    end forall

    forall ( i = 1:(n-1), j = 1:(n-1), i == j )
       D(i,j) = (eta(i) - alpha) / (2._rk * eta(i))
    end forall

    forall ( i = 1:(n-1), j = 1:(n-1), i /= j )
       D(i,j) = L(i) / (L(j) * (eta(i) - eta(j))) 
    end forall
    ! end filling the matrix

    ! done here
  end subroutine comp_exact_D_scaled_laguerre

  ! computes the integration matrix (operator) arising
  ! from general Jacobi polynomials P_n^{alpha, beta}(x)
  ! and also the corresponding point distribution
  ! on the given open interval x = (x0 x1) (n points)

  subroutine comp_jacobi_integ_matrix(alpha, beta, n, S, SS, x)
    implicit none
    real(rk), intent(in) :: alpha, beta
    integer , intent(in) :: n
    real(rk), dimension(n, n), intent(out) :: S, SS
    real(rk), dimension(n), intent(out) :: x

    ! local vars
    real(rk), dimension(:,:), allocatable :: D, DD
    real(rk), dimension(:), allocatable :: xx

    ! LAPACK vars
    real(rk), dimension(n+1) :: work
    integer , dimension(n+1) :: ipiv
    integer nnn, info

    ! external subroutines
    external DGETRF
    external DGETRI

    ! local allocations
    allocate( D(n+1, n+1), DD(n+1, n+1) &
            , xx(n+1))

    ! compute the Jacobi differentiation matrix
    call comp_jacobi_diff_matrix(alpha, beta, n, D, DD, xx)
    ! return the integration grid
    x = xx(2:(n+1))

    !============================================================
    ! FIRST INTEGRATION OPERATOR :  S = D(2:(n+1), 2:(n+1))^(-1)
    !============================================================
    ! initialize S with D(2:n+1, 2:n+1) where it will be overwritten
    ! by LU factors computed by LAPACK 
    S = D(2:(n+1), 2:(n+1))
    nnn = n ! size n in LAPACK
    ! do LU factorization using LAPACK
    call DGETRF(nnn, nnn, S, nnn, ipiv, info)
    if( info /= 0 ) then 
       print *,  ' D(2:(n+1), 2:(n+1)) is numerically singular or ' &
              ,  'LU can not be computed!'
       stop
    end if
    ! compute inverse using already computed LU
    call DGETRI(nnn, S, nnn, ipiv, work, nnn, info)
    if( info /= 0 ) then 
       stop 'inversion of <D> failed!'
    end if

    !============================================================
    ! SECOND INTEGRATION OPERATOR: SS = DD(2:(n+1), 2:(n+1))^(-1)
    !============================================================
    ! ! initialize SS with DD(2:n+1, 2:n+1) where it will be overwritten
    ! ! by LU factors computed by LAPACK 
    ! SS = DD(2:(n+1), 2:(n+1))
    ! nnn = n ! size n in LAPACK
    ! ! do LU factorization using LAPACK
    ! call DGETRF(nnn, nnn, SS, nnn, ipiv, info)
    ! if( info /= 0 ) then 
    !    print *,  ' DD(2:(n+1), 2:(n+1)) is numerically singular or ' &
    !           ,  'LU can not be computed!'
    !    stop
    ! end if
    ! ! compute inverse using already computed LU
    ! call DGETRI(nnn, SS, nnn, ipiv, work, nnn, info)
    ! if( info /= 0 ) then 
    !    stop 'inversion of <DD> failed!'
    ! end if

    SS = matmul(S, S) ! NOT VERY CONSISTENT and ACCURATE!

    ! clean ups
    deallocate( D, DD, xx)

    ! done here
  end subroutine comp_jacobi_integ_matrix
  
  ! computes the Gauss integration matrix [S]
  ! when multiplied to the vector [f] which 
  ! is the smapled values of function f(x)
  ! at "n" collocated points "x ={xk}", it gives
  ! the initial function evaluated on a
  ! set of collocation points "x" such that
  ! 
  !
  !   [S] * [f] =  int _{-1}^{xk} f(xi) dxi + O(dx^n)
  ! 
  ! this formula is exact for case where f is 
  ! a polynomial of degree 2n-1.
  !    
  subroutine comp_gauss_integ_matrix(alpha, beta, n, S, x)
    implicit none
    real(rk), intent(in) :: alpha, beta
    integer , intent(in) :: n
    real(rk), dimension(n, n), intent(out) :: S
    real(rk), dimension(n), intent(out) :: x

    !local vars
    real(rk), dimension(n+1, n+1) :: D, uixk, duixk, Psi
    real(rk), dimension(n+1) :: xx, dxx
    real(rk) ::  tmp 
    integer :: i, k

    ! LAPACK vars
    real(rk), dimension(n+1) :: work
    integer , dimension(n+1) :: ipiv
    integer nnn, info

    ! external subroutines
    external DGETRF
    external DGETRI

    ! first compute the zeros of Jacobi polynomial
    ! of degree n
    call ZEJAGA(n, alpha, beta, xx, dxx)
    ! add -1 to the first of it
    xx(2:(n+1)) = xx(1:n)
    xx(1) = -1._rk  
    x = xx(2:(n+1)) ! return this

    ! Then find the Fourier coefficients
    ! compute "uixk" and "duixk" matrices by evaluating
    ! the value and the derivative of Jacobi polynomials
    ! for orders i = 0:n at xk collocation points
    ! given by x-distribution for k = 1:(n+1)
    ! the results should be (n+1)x(n+1) matrices.
    do i = 0 , n
       do k = 1 , (n+1)
          call VAJAPO(i, alpha, beta, xx(k), uixk(k,i+1) &
               , duixk(k,i+1), tmp)
       end do
    end do
    ! initialize Psi with uixk where it will be overwritten
    ! by LU factors computed by LAPACK 
    Psi = uixk
    nnn = n+1 ! size n in LAPACK
    ! do LU factorization using LAPACK
    call DGETRF(nnn, nnn, Psi, nnn, ipiv, info)
    if( info /= 0 ) then 
       print *,  'uixk is numerically singular or ' &
            ,  'LU can not be computed!' &
            , ' error in comp_gauss_integ_matrix(...)'
       stop
    end if
    ! compute inverse using already computed LU
    call DGETRI(nnn, Psi, nnn, ipiv, work, nnn, info)
    if( info /= 0 ) then 
       stop 'inversion of <uixk> failed!'
    end if

    ! computing the final Jacobi 1st differentiation matrix
    ! on Guass-Jacobi points
    D  = matmul( duixk, Psi)

    ! Proceeding to compute the integration operator ...
    S = D(2:(n+1), 2:(n+1))
    nnn = n ! size n in LAPACK
    ! do LU factorization using LAPACK
    call DGETRF(nnn, nnn, S, nnn, ipiv, info)
    if( info /= 0 ) then 
       print *,  'S is numerically singular or ' &
            ,  'LU can not be computed!'
       stop
    end if
    ! compute inverse using already computed LU
    call DGETRI(nnn, S, nnn, ipiv, work, nnn, info)
    if( info /= 0 ) then 
       stop 'inversion of <S> failed!'
    end if

    ! done here
  end subroutine comp_gauss_integ_matrix

  ! finds the norm-inf error between the analytically
  ! computed weights for Jacobi-Gauss quadrature
  ! and the matrix-based computation.
  ! the final difference between these values
  ! is returned in "err" variable.  
  subroutine test_gauss_integ_matrix(alpha, beta, n, err)
    implicit none
    real(rk), intent(in) :: alpha, beta
    integer , intent(in) :: n
    real(rk), intent(out) :: err

    !local vars
    real(rk), dimension(n, n) :: S
    real(rk), dimension(n) :: x, dx, ww

    ! first find zeros and derivative at zeros
    ! of the Jacobi polynomial
    call ZEJAGA(n, alpha, beta, x, dx)
    ! then compute the Jacobi-Gauss integration operator
    ! NOTE that "x" remains the same in this case
    call comp_gauss_integ_matrix(alpha, beta, n, S, x)
    ! then compute analytical weights
    call WEJAGA(n, alpha, beta, x, dx, ww)
    ! compute inf-norm of error
    err = maxval(abs(S(n,:) - ww))

    ! done here
  end subroutine test_gauss_integ_matrix

end module ortho_poly
