module precond_opt
  use sparse_opt
  implicit none

  private


  ! preconditioner
  type precond

     ! band LU Lapack data struct. section
     integer :: n, kl, ku, ldab, INFO
     integer, dimension(:), allocatable :: ipiv
     real*8, dimension(:, :), allocatable :: AB, BB
     logical :: enabled = .false.

   contains
     procedure :: init => init_precond
     procedure :: solve => solve_precond
  end type precond


  public :: precond

contains

  function get_global_scalar(Msp, zeta, Ksp, ifree, neqs, ntot, i, j)
    implicit none
    class(smat), intent(in) :: Msp
    real*8, dimension(:, :), intent(in) :: zeta
    class(smat), intent(in) :: Ksp
    integer, dimension(:), intent(in) :: ifree
    integer, intent(in) :: neqs, ntot, i, j
    real*8 :: get_global_scalar

    ! locar vars
    integer :: itim, jtim, ii, jj
    integer :: inode, jnode, iscalar, jscalar
    real*8, dimension(size(Ksp%a, 1), size(Ksp%a, 1)) :: tmp
    real*8 :: Mval, Kval

    itim = ceiling(dble(i) / dble(ntot)) 
    jtim = ceiling(dble(j) / dble(ntot))

    ii = i - (itim - 1) * ntot
    jj = j - (jtim - 1) * ntot

    inode = ceiling(dble(ii) / dble(neqs)) 
    jnode = ceiling(dble(jj) / dble(neqs))

    iscalar = ii - (inode - 1) * neqs
    jscalar = jj - (jnode - 1) * neqs

    tmp = Msp%get(ifree(inode),ifree(jnode))  
    Mval = tmp(iscalar, jscalar)
    tmp = Ksp%get(ifree(inode),ifree(jnode))  
    Kval = tmp(iscalar, jscalar)

    get_global_scalar = Mval + zeta(itim, jtim) * Kval

    ! done here
  end function get_global_scalar

  subroutine init_precond(this, Msp, zeta, Ksp, ifree, kl, ku)
    implicit none
    class(precond), intent(inout) :: this
    class(smat), intent(in) :: Msp
    real*8, dimension(:, :), intent(in) :: zeta
    class(smat), intent(in) :: Ksp
    integer, dimension(:), intent(in) :: ifree
    integer, intent(in) :: kl, ku

    ! local vars
    integer :: neqs, w, ntot , nS
    integer :: i, j, i1, i2

    ! inits
    neqs = size(Ksp%a, 1)
    w = size(ifree)
    ntot = w * neqs
    nS = size(zeta, 1)
    this%n = ntot * nS 
    this%kl = kl
    this%ku = ku
    this%ldab = 2*kl+ku+1
    this%INFO = 0

    if( allocated(this%ipiv) ) deallocate(this%ipiv)
    allocate(this%ipiv(this%n))
    if( allocated(this%AB) ) deallocate(this%AB)
    allocate(this%AB(this%ldab, this%n))
    this%AB = 0.0d0
    if( allocated(this%BB) ) deallocate(this%BB)
    allocate(this%BB(this%n, 1))
    this%BB = 0.0d0

    ! pack space-time sparse matrix into LAPACK band storage matrix
    do j = 1, this%n
       i1 = maxval((/ 1, j - ku /))
       i2 = minval((/ this%n, j + kl /))
       do i = i1, i2
          this%AB(kl+ku+1+i-j, j) = get_global_scalar(Msp, zeta, Ksp, ifree, neqs, ntot, i, j)
       end do
    end do

    ! perform band LU and store it
    call DGBTRF( this%n, this%n, kl, ku, this%AB, this%ldab, this%ipiv, this%INFO )
    if ( this%INFO .ne. 0 ) then
       print *, 'something is wrong in band LU decomposition! INFO = ', this%INFO
       stop
    end if

    ! lock it then 
    this%enabled = .true.

    ! done here
  end subroutine init_precond


  subroutine solve_precond(this, b, x)
    implicit none
    class(precond), intent(inout) :: this
    real*8, dimension(:), intent(in) :: b
    real*8, dimension(:), intent(out) :: x

    if ( .not. this%enabled ) then
       x = b
       return
    end if
 
    this%BB(:, 1) = b

    call DGBTRS( 'No transpose', this%n, this%kl, this%ku, 1, this%AB &
         , this%ldab, this%ipiv, this%BB, this%n, this%INFO )

    if ( this%INFO .ne. 0 ) then
       print *, 'something is wrong in band LU solver DGBTRS! INFO = ', this%INFO
       stop
    end if

    x = this%BB(:, 1)

    ! done here
  end subroutine solve_precond

end module precond_opt
