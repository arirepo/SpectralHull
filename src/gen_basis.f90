module gen_basis
  use grid_opt
  ! generates arbitrary 2D basis functions
  ! using Vandermonde matrix operations.
  implicit none

  private

  type basis
     private

     real*8, dimension(:,:), allocatable :: MT, M0
     integer, dimension(:), allocatable :: piv
     integer :: d
     integer :: elname
     procedure(gen_eval_M_pt), private, pointer, nopass  :: eval_M_pt

   contains

     procedure, public :: init => initialize
     procedure, public :: eval => evaluate 

  end type basis
  
  abstract interface
     subroutine gen_eval_M_pt(xx, yy, d, op, MT)
       implicit none
       real*8, intent(in) :: xx, yy
       integer, intent(in) :: d, op  
       real*8, dimension(:), intent(out) :: MT
     end subroutine gen_eval_M_pt
  end interface

  public :: basis

contains

  ! computes terms in pascal triangle and also their derivatives
  ! with respect to x, y.
  !
  ! inputs:
  !
  ! the terms are evaluated at PHYSICAL point (xx, yy) which can
  ! be inside master element as a special case.
  ! d is the degree of the 2d basis functions  
  ! the operation: (0=pascal), (1=dpascal/dx), (2=dpascal/dy)
  !
  ! outputs:
  !
  ! MT is a real vector containing all pascal terms or
  ! their derivatives d/dx or d/dy at point (xx, yy)
  !
  ! MT = [sum_{i=0}^d sum_{j=0}^i xx**(i-j) * yy**j]
  !
  subroutine comp_pascal_tri(xx, yy, d, op, MT)
    implicit none
    real*8, intent(in) :: xx, yy
    integer, intent(in) :: d, op  
    real*8, dimension(:), intent(out) :: MT

    ! local vars
    integer :: id, i, j, jj

    id = size(MT) !

    ! bug checking
    if ( id .ne. ( (d+1)*(d+2)/2 ) ) then
       print *, 'the length of MT does not match with' &
            , ' order of requested polynomial. stop.'
       stop
    end if
    if ( d < 1 ) then
       print *, 'the degree of requested polynomials should be >= 1'
       stop
    end if

    ! fill it!    
    jj = 1
    do i = 0, d
       do j = 0, i

          if     ( op .eq. 0 ) then !pascal terms
             MT(jj) = (xx**(i-j)) * (yy**j)  
          elseif ( op .eq. 1 ) then !d/dx of pascal terms
             if ((i-j-1) < 0) then
                MT(jj) = 0.0d0
             else
                MT(jj) = dble(i-j) * (xx**(i-j-1)) * (yy**j) 
             end if
          elseif ( op .eq. 2 ) then !d/dy of pascal terms
             if ( (j-1) < 0) then
                MT(jj) = 0.0d0
             else
                MT(jj) = (xx**(i-j)) * dble(j) * (yy**(j-1)) 
             end if
          else
             print *, 'unknown operation in comp_pascal_tri! stop.'
             stop 
          end if

          jj = jj + 1
       end do
    end do

    ! done here
  end subroutine comp_pascal_tri

  ! computes basic polynomials for quad elem and also their derivatives
  ! with respect to x, y.
  !
  ! inputs:
  !
  ! the terms are evaluated at PHYSICAL point (xx, yy) which can
  ! be inside master element as a special case.
  ! d is the degree of the 2d basis functions  
  ! the operation: (0=poly), (1=dpoly/dx), (2=dpoly/dy)
  !
  ! outputs:
  !
  ! MT is a real vector containing all polynomial terms or
  ! their derivatives d/dx or d/dy at point (xx, yy)
  !
  ! MT = [sum_{i=0}^d sum_{j=0}^d xx**i * yy**j]
  !
  subroutine comp_poly_quad(xx, yy, d, op, MT)
    implicit none
    real*8, intent(in) :: xx, yy
    integer, intent(in) :: d, op  
    real*8, dimension(:), intent(out) :: MT

    ! local vars
    integer :: id, i, j, jj

    id = size(MT) !

    ! bug checking
    if ( id .ne. ( (d+1)**2 ) ) then
       print *, 'the length of MT does not match with' &
            , ' order of requested polynomials for quad! stop.'
       stop
    end if
    if ( d < 1 ) then
       print *, 'the degree of requested polynomials should be >= 1 for quad!'
       stop
    end if

    ! fill it!    
    jj = 1
    do i = 0, d
       do j = 0, d

          if     ( op .eq. 0 ) then !poly terms
             MT(jj) = (xx**i) * (yy**j)  
          elseif ( op .eq. 1 ) then !d/dx of poly
             if (i .eq. 0) then
                MT(jj) = 0.0d0
             else
                MT(jj) = dble(i) * (xx**(i-1)) * (yy**j) 
             end if
          elseif ( op .eq. 2 ) then !d/dy of poly
             if ( j .eq. 0) then
                MT(jj) = 0.0d0
             else
                MT(jj) = (xx**i) * dble(j) * (yy**(j-1)) 
             end if
          else
             print *, 'unknown operation in comp_poly_quad! stop.'
             stop 
          end if

          jj = jj + 1
       end do
    end do

    ! done here
  end subroutine comp_poly_quad

  ! --------------------------------------------
  ! NOTE : The matrix "M" and its transpose "MT" 
  ! used below, is in fact the matrix "V" used 
  ! in module "approx_fekete". It is just different
  ! symbols for the same thing! 
  ! --------------------------------------------
  ! creates MT matrix and then
  ! performs LU factorization of MT matrix.
  ! later this will be used to solve the system and
  ! evaluate basis functions and their 
  ! derivatives at a point 
  subroutine comp_lu_MT(this, x, y)
    implicit none
    class(basis), target, intent(inout) :: this
    real*8, dimension(:), intent(in) :: x, y

    ! local vars
    integer :: i, id, INFO, d, op
    real*8, dimension(:,:), pointer :: MT => null()
    integer, dimension(:), pointer :: piv => null()

    MT  => this%MT
    piv => this%piv
    d = this%d
    op = 0 ! only interpolation

    id = size(x) ! number of points

    ! bug checking
    if ( (id .ne. size(MT,1)) .or. & 
         (id .ne. size(MT,2)) ) then
       print *, 'MT matrix is not (numberofpoints*numberofpoints). stop'
       stop
    end if
    if ( size(piv) .ne. id ) then
       print *, 'the length of pivot vector is not numberofpoints. stop'
       stop
    end if

    ! start filling MT column wise
    do i = 1, id
       call this%eval_M_pt(x(i), y(i), d, op, MT(:,i))
    end do

    ! now, perform lu of MT and store in place
    call DGETRF( id, id, MT, id, piv, INFO )
    if (INFO .ne. 0) then
       print *, 'somethig is wrong in LU decomposition in MT! stop.'
       stop
    end if

    ! done here
  end subroutine comp_lu_MT

  
  ! initializes the basis data type
  subroutine initialize(this, x, y, elname)
    implicit none
    class(basis), intent(inout) :: this
    real*8, dimension(:), intent(in) :: x, y
    integer, intent(in) :: elname

    ! local vars
    integer :: id
    real*8  :: delta

    id = size(x)

    ! bulletproofing
    if ( allocated(this%MT) ) then
       print *, 'this object for basis function is already initialized! stop'
       stop
    end if

    !
    allocate(this%MT(id, id), this%M0(id, 1), this%piv(id))
    this%M0 = 0.0d0
    this%elname = elname

    ! obtaining the degree of required polynomial
    ! from the given number of points and performing
    ! element specific initialization
    select case (this%elname)
    case (GEN_TRIANGLE ) ! do it the old way

       delta = sqrt(1.0d0 + 8.0d0 * dble(id))
       this%d = maxval( (/ nint(-1.5d0 + 0.5d0 * delta) &
            , nint(-1.5d0 - 0.5d0 * delta) /) )

       if ( this%d <= 0 ) then
          print *, 'degree of basis is this%d <= 0! stop.'
          stop
       end if

       ! initialize generic evaluation procedure pointer
       this%eval_M_pt => comp_pascal_tri
 
    case ( GEN_QUADRI )

       this%d = nint(sqrt(dble(id))) - 1

       ! initialize generic evaluation procedure pointer
       this%eval_M_pt => comp_poly_quad

    case default

       print *, 'unregognized element name in init. basis func.! stop'
       stop

    end select

    ! fill MT matrice
    call comp_lu_MT(this, x, y)

    ! done here
  end subroutine initialize

  ! evaluates the basis function or
  ! it derivatives d/dx, d/dy at 
  ! some point (x0, y0)
  ! op = 0 => evaluate the basis function
  ! op = 1 => evaluate the d/dx of basis
  ! op = 2 => evaluate the d/dy of basis
  subroutine evaluate(this, x0, y0, op, val)
    implicit none
    class(basis), intent(inout) :: this
    real*8, intent(in) :: x0, y0
    integer, intent(in) :: op
    real*8, dimension(:), intent(out) :: val

    ! local vars
    integer :: id, INFO

    ! bullet proofing
    if (.not. allocated(this%MT)) then
       print *, 'fatal: please first initialize this basis object' &
            , ' before evaluating. stop.'
       stop
    end if

    id = size(this%MT, 1)

    if ( size(val) .ne. id) then
       print *, 'in evaluating basis functions, size of' &
            , ' output array size<val> = ', size(val), ' is not equal to' &
            , ' the number of basis functions (numberofpoints = ', id, '). stop.'
       stop
    end if

    ! first fill this%M0 column vector with
    ! pascal terms at point (x0, y0) according to the requested operation
    call this%eval_M_pt(x0, y0, this%d, op, val)

    ! solve psi = basis = MT\M0 using already computed LU
    CALL DGETRS( 'No transpose', id, 1, this%MT, id, this%piv, val, id, INFO )
    if ( INFO .ne. 0) then
       print *, 'could not solve to evaluate basis at point (x0, y0)'
       stop
    end if

    ! val = this%M0(:,1)

    ! done here
  end subroutine evaluate

end module gen_basis

! program tester
!   use gen_basis
!   implicit none

!   integer :: i, npelem
!   real*8 :: x0, y0
!   real*8, dimension(:), allocatable :: xi, eta

!   ! first order element
!   print *, 'first order element :'
!   npelem = 3
!   allocate( xi(npelem), eta(npelem))

!   xi  = (/ 0.0d0, 1.0d0, 0.0d0 /)
!   eta = (/ 0.0d0, 0.0d0, 1.0d0 /)  
  
!   do i = 1, size(xi)

!      x0 = xi(i); y0 = eta(i)
!      call print_basis(xi, eta, x0, y0)

!   end do

!   deallocate(xi, eta)

!   ! second order element
!   print *, 'second order element:'
!   npelem = 6
!   allocate( xi(npelem), eta(npelem))

!   xi  = (/ 0.0d0, 0.5d0, 1.0d0, 0.5d0, 0.0d0, 0.0d0 /)
!   eta = (/ 0.0d0, 0.0d0, 0.0d0, 0.5d0, 1.0d0, 0.5d0 /)  
  
!   do i = 1, size(xi)

!      x0 = xi(i); y0 = eta(i)
!      call print_basis(xi, eta, x0, y0)

!   end do

!   deallocate(xi, eta)

!   ! done here

! contains

!   subroutine print_basis(xi, eta, x0, y0)
!     implicit none
!     real*8, dimension(:), intent(in) :: xi, eta
!     real*8, intent(in) :: x0, y0

!     ! local vars
!     type(basis) :: tbasis
!     real*8, dimension(size(xi)) :: val

!     call tbasis%init(xi, eta)

!     print *, 'at point (x0, y0) = ', x0, ',', y0
!     call tbasis%eval(x0, y0, 0, val)
!     print *, 'psi = ', val
!     call tbasis%eval(x0, y0, 1, val)
!     print *, 'd psi/dx = ', val
!     call tbasis%eval(x0, y0, 2, val)
!     print *, 'd psi/dy = ', val

!     ! done here
!   end subroutine print_basis

! end program tester
