module sparse_opt
  implicit none

  private

  type ints
     integer, dimension(:), allocatable :: ja
  end type ints

  ! implements block sparse storage for matrix
  ! with double precision entries and overloads
  ! some important operators
  type smat
     integer :: n, m, nnz
     integer, dimension(:), allocatable :: ia, ja
     ! (1..neqs, 1..neqs, 1..nnz)
     real*8, dimension(:, :, :), allocatable :: a

   contains
     procedure :: init => init_smat
     procedure :: smatmul => sparse_matmul
     procedure :: set => set_smat
     procedure :: get => get_smat
     procedure :: get_scalar => get_scalar_smat
     procedure :: sp2full => convert_to_full_matrix

  end type smat

  public :: ints, smat

contains

  subroutine init_smat(this, rows, neqs)
    implicit none
    class(smat), intent(inout) :: this
    type(ints), dimension(:), intent(in) :: rows
    integer, intent(in) :: neqs

    ! local vars
    integer :: i, j, ii
    
    this%n = size(rows)
    this%m = 1 ! at least has one column 
    do i = 1, size(rows)
       this%m = maxval( (/ this%m, rows(i)%ja(:)/) )
    end do

    this%nnz = 0
    do i = 1, size(rows)
       this%nnz = this%nnz + size(rows(i)%ja)
    end do

    ! compute ia
    if( allocated(this%ia) ) deallocate(this%ia)
    allocate(this%ia(this%n + 1))
    this%ia(1) = 1
    do i = 2, (this%n + 1)
       this%ia(i) = this%ia(i-1) + size(rows(i-1)%ja)
    end do

    ! compute ja
    if(allocated(this%ja)) deallocate(this%ja)
    allocate(this%ja(this%nnz))
    ii = 1
    do i = 1, this%n
       do j = 1, size(rows(i)%ja)
          this%ja(ii) = rows(i)%ja(j)
          ii = ii + 1
       end do
    end do

    if (allocated(this%a)) deallocate(this%a)
    allocate(this%a(neqs, neqs, this%nnz))
    this%a = 0.0d0

    ! done here
  end subroutine init_smat

  subroutine sparse_matmul(this, x, xx)
    implicit none
    class(smat), intent(in) :: this
    real*8, dimension(size(this%a, 1), this%m) , intent(in) :: x
    real*8, dimension(size(this%a, 1), this%m) , intent(out) :: xx

    ! local vars
    integer :: i, j, j1, j2

    ! Hard RESET
    xx = 0.0d0

    do i = 1, this%n
       j1 = this%ia(i)
       j2 = this%ia(i+1) - 1
       do j = j1, j2
          xx(:, i) = xx(:, i) + matmul( this%a(:, :, j) , x(:, this%ja(j)) )
       end do
    end do

    ! done here
  end subroutine sparse_matmul

  subroutine set_smat(this, i, j, val)
    implicit none
    class(smat), intent(inout), target :: this
    integer, intent(in) :: i, j
    real*8, dimension(size(this%a,1), size(this%a,1)), intent(in) :: val

    ! local vars
    integer :: j1, j2, jj, loc
    integer, dimension(:), pointer :: sec => null()

    j1 = this%ia(i)
    j2 = this%ia(i+1) - 1
    sec => this%ja(j1:j2)
    jj = minloc(array = sec, dim = 1, mask = (sec .eq. j) )
    if ( jj .eq. 0 ) then
       ! do nothing
    else
       loc = j1 + jj - 1 
       this%a(:, :, loc) = val
    end if

    ! done here
  end subroutine set_smat

  function get_smat(this, i, j)
    implicit none
    class(smat), intent(in), target :: this
    integer, intent(in) :: i, j
    real*8, dimension(size(this%a,1), size(this%a,1)) :: get_smat

    ! local vars
    integer :: j1, j2, jj, loc
    integer, dimension(:), pointer :: sec => null()

    j1 = this%ia(i)
    j2 = this%ia(i+1) - 1
    sec => this%ja(j1:j2)
    jj = minloc(array = sec, dim = 1, mask = (sec .eq. j) )
    if ( jj .eq. 0 ) then
       get_smat = 0.0d0
    else
       loc = j1 + jj - 1 
       get_smat = this%a(:, :, loc)
    end if

    ! done here
  end function get_smat

  function get_scalar_smat(this, i, j)
    implicit none
    class(smat), intent(in) :: this
    integer, intent(in) :: i, j
    real*8 :: get_scalar_smat

    ! local vars
    integer :: neqs, blk_i, blk_j, di, dj
    real*8, dimension(size(this%a, 1), size(this%a, 1)) :: blk

    neqs = size(this%a, 1)

    blk_i = ceiling(dble(i) / dble(neqs))
    blk_j = ceiling(dble(j) / dble(neqs))
    blk = this%get(blk_i, blk_j)

    di = i - (blk_i - 1) * neqs - 1
    dj = j - (blk_j - 1) * neqs - 1
    get_scalar_smat = blk( 1 + di, 1 + dj)

    ! done here
  end function get_scalar_smat

  subroutine convert_to_full_matrix(this, A)
    implicit none
    class(smat), intent(in) :: this
    real*8, dimension(:, :), allocatable :: A

    ! local vars
    integer :: i, j, neqs
    integer :: i1, i2, j1, j2
    real*8, dimension(size(this%a, 1), size(this%a, 1)) :: tmp

    ! init
    neqs = size(this%a, 1)

    ! bullet proofing ...
    if ( allocated(A) ) deallocate(A)
    allocate(A(neqs * this%n, neqs * this%m))
 
    do i = 1, this%n
       do j = 1, this%m
          i1 = (i-1) * neqs + 1
          i2 = i * neqs
          j1 = (j-1) * neqs + 1
          j2 = j * neqs

          tmp = this%get(i, j)

          A(i1:i2, j1:j2) = tmp
       end do
    end do

    ! done here
  end subroutine convert_to_full_matrix

end module sparse_opt

! program tester
!   use sparse_opt
!   implicit none

!   ! local vars
!   integer :: neqs, i, j
!   type(ints), dimension(:), allocatable :: rows
!   type(smat), target :: tsmat
!   real*8, dimension(:, :), allocatable :: x, xx, val
!   real*8, dimension(:, :, :, :), allocatable :: full

!   ! sample matrix
!   ! 3 - 4- 5
!   ! |   | /
!   ! 2 - 1
!   !
!   ! * * - * *
!   ! * * * - -
!   ! - * * * -
!   ! * - * * *
!   ! * - - * *
!   allocate(rows(5))
!   allocate(rows(1)%ja(4))
!   rows(1)%ja = (/ 1, 2, 4, 5 /)
!   allocate(rows(2)%ja(3))
!   rows(2)%ja = (/ 1, 2, 3/)
!   allocate(rows(3)%ja(3))
!   rows(3)%ja = (/ 2, 3, 4 /)
!   allocate(rows(4)%ja(4))
!   rows(4)%ja = (/ 1, 3, 4, 5 /)
!   allocate(rows(5)%ja(3))
!   rows(5)%ja = (/ 1, 4, 5 /)
!   neqs = 1

!   call tsmat%init(rows, neqs)
!   tsmat%a = 1.0d0

!   allocate(x(neqs, 5), xx(neqs, 5))
!   x = 1.0d0

!   call tsmat%smatmul(x, xx)
!   print *, 'xx = ', xx

!   allocate(val(1,1))
!   val = 19.0d0
!   call tsmat%set(1,5, val)

!   print *, 'tsmat%n = ', tsmat%n
!   print *, 'tsmat%m = ', tsmat%m
!   print *, 'tsmat%nnz = ', tsmat%nnz
!   print *, 'tsmat%ia = ', tsmat%ia
!   print *, 'tsmat%ja = ', tsmat%ja
!   print *, 'tsmat%a = ', tsmat%a 

!   print *, 'sh = ', shape(tsmat%get(1,1)), tsmat%get(1,3)

!   ! convert to full matrix
!   allocate(full(1, 1, tsmat%n, tsmat%n))
!   do i = 1, tsmat%n
!      do j = 1, tsmat%n
!         full(:, :, i, j) = tsmat%get_scalar(i,j) !tsmat%get(i, j)
!      end do
!   end do

!   print *, 'here is full matrix :' 

!   do i = 1, tsmat%n
!      print *, (/ ( full(1,1, i, j) , j = 1, tsmat%n ) /)
!   end do

!   ! test get_scalar
!   tsmat%a(1, 1, :) = (/ ( dble(i), i = 1, tsmat%nnz) /)

!   print *, tsmat%a
!   print *, 'get_scalar(1,3) = ', tsmat%get_scalar(1,5)

!   ! done here
! end program tester

