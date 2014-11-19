module fem_reordering
  implicit none

  private

  public :: fem_reorder, arash_sort 
contains

  ! sorts the real*8 vector <x> without touching it!
  ! instead, returns the permutation vector <p> that
  ! corrresponds to the sort operations based on:
  !
  ! typ = '+'  ====> sorts smaller to larger values
  ! typ = '-'  ====> sorts larger to smaller values
  ! 
  subroutine arash_sort(x, p, typ)
    implicit none
    real*8, dimension(:), intent(in) :: x
    integer, dimension(size(x)), intent(out) :: p
    character, intent(in) :: typ

    ! local vars
    integer :: n, i, j, tmp

    ! bulletproofing
    select case (typ)
    case ('+', '-')
       ! do nothing
    case default
       print *, 'fatal: unknown sorting operation! stop'
       stop
    end select

    n = size(x)

    ! initial permutation vector
    p = (/ (i , i = 1, n) /)

    do i = 1, (n-1)
       do j = (i+1), n

          if (  ( ( x(p(i)) > x(p(j)) ) & 
               .and. (    typ .eq. '+'   ) ) &
               .or.  ( ( x(p(i)) < x(p(j)) ) & 
               .and. (    typ .eq. '-'   ) ) ) then
             tmp = p(j)
             p(j) = p(i)
             p(i) = tmp
          end if

       end do
    end do


    ! done here
  end subroutine arash_sort

  !
  ! reorders the tuple:
  ! (xi, eta) :: xy(1,:) = xi, xy(2,:) = eta
  ! of the master element in the desired
  ! counterclockwise FEM orientation described below:
  !
  ! 
  !     [3]
  !      |  \
  !     (10) (9)
  !      |      \
  !edg3 (11)     (8) <edg2>
  !      |   no     \
  !     (12)  order  (7)
  !      |             \
  !     [1]-(4)-(5)-(6)-[2]
  !         <edg1>
  ! interior points are left untouched on 
  ! and have the original ordering.

  subroutine fem_reorder(xy, tol)
    implicit none
    real*8, dimension(:,:), intent(inout) :: xy
    real*8, intent(in) :: tol

    ! local vars
    integer :: i, n, indx, jj, n_tmp
    logical,      dimension(size(xy, 2))          :: reordered 
    real*8, dimension(size(xy,1), size(xy, 2))    :: xy_tmp
    real*8, dimension(:, :), allocatable    :: tmp, tmp2
    integer, dimension(:), allocatable :: p !permutation

    ! bulletproofing
    if ( size(xy, 1) .ne. 2 ) then
       print *, ' currently <xy> must be xy(2, :). The first' &
            , ' dimension did not match! stop'
       stop
    end if

    ! init section
    ! ------------------------
    ! total number of points in a master element 
    n = size(xy, 2)
    reordered = .false.
    xy_tmp = 0.0d0
    indx = 1

    ! find and reorder point 1, 2, 3
    allocate(tmp(2,3))
    tmp(1,:) = (/ 0.0d0, 1.0d0, 0.0d0 /)
    tmp(2,:) = (/ 0.0d0, 0.0d0, 1.0d0 /)
    call add2list(xy, tmp, tol, xy_tmp, indx, reordered)
    ! print *, 'xy_tmp = [', xy_tmp, ']'
    ! stop 
    deallocate(tmp)

    ! 1- find edg1 points
    jj = 1
    do i = 1, n
       if(    (abs(xy(2, i) - 0.0d0) <= tol ) &
            .and. (        .not. reordered(i)     ) ) then
          allocate( tmp2(2, jj) )
          if ( jj > 1) tmp2(:,1:(jj-1) ) = tmp(:,1:(jj-1) )
          tmp2(:,jj) = xy(:, i)
          call move_alloc(tmp2, tmp)
          jj = jj + 1
       end if
    end do
    ! 2- sort edg1 points
    n_tmp = size(tmp, 2)
    if(allocated(tmp) .and. (n_tmp >= 1) ) then 

       allocate(p(n_tmp), tmp2(2, n_tmp))
       call arash_sort(tmp(1,:), p, '+')
       tmp2(:,1:n_tmp) = tmp(:,p)
       tmp = tmp2
       deallocate(p, tmp2) 
       ! 3- add reordered edg1 points to final list    
       call add2list(xy, tmp, tol, xy_tmp, indx, reordered)
       deallocate(tmp) 
    end if

    ! 1- find edg2 points
    jj = 1
    do i = 1, n
       if(    (abs(xy(1, i) + xy(2, i) - 1.0d0) <= tol ) &
            .and. (        .not. reordered(i)     ) ) then
          allocate( tmp2(2, jj) )
          if ( jj > 1) tmp2(:,1:(jj-1) ) = tmp(:,1:(jj-1) )
          tmp2(:,jj) = xy(:, i)
          call move_alloc(tmp2, tmp)
          jj = jj + 1
       end if
    end do
    ! 2- sort edg2 points
    n_tmp = size(tmp, 2)
    if(allocated(tmp) .and. (n_tmp >= 1) ) then 

       allocate(p(n_tmp), tmp2(2, n_tmp))
       call arash_sort(tmp(1,:), p, '-')
       tmp2(:,1:n_tmp) = tmp(:,p)
       tmp = tmp2
       deallocate(p, tmp2) 
       ! 3- add reordered edg2 points to final list    
       call add2list(xy, tmp, tol, xy_tmp, indx, reordered)
       deallocate(tmp) 
    end if

    ! 1- find edg3 points
    jj = 1
    do i = 1, n
       if(    (abs(xy(1, i) - 0.0d0) <= tol ) &
            .and. (        .not. reordered(i)     ) ) then
          allocate( tmp2(2, jj) )
          if ( jj > 1) tmp2(:,1:(jj-1) ) = tmp(:,1:(jj-1) )
          tmp2(:,jj) = xy(:, i)
          call move_alloc(tmp2, tmp)
          jj = jj + 1
       end if
    end do
    ! 2- sort edg3 points
    n_tmp = size(tmp, 2)
    if(allocated(tmp) .and. (n_tmp >= 1) ) then 

       allocate(p(n_tmp), tmp2(2, n_tmp))
       call arash_sort(tmp(2,:), p, '-')
       tmp2(:,1:n_tmp) = tmp(:,p)
       tmp = tmp2
       deallocate(p, tmp2) 
       ! 3- add reordered edg3 points to final list    
       call add2list(xy, tmp, tol, xy_tmp, indx, reordered)
       deallocate(tmp) 
    end if


    ! add interior points by filling tmp first (like above)
    ! and just adding without sorting
    ! 1- find interior points
    jj = 1
    do i = 1, n
       if(    .not. reordered(i)    ) then
          allocate( tmp2(2, jj) )
          if ( jj > 1) tmp2(:,1:(jj-1) ) = tmp(:,1:(jj-1) )
          tmp2(:,jj) = xy(:, i)
          call move_alloc(tmp2, tmp)
          jj = jj + 1
       end if
    end do
    ! 2- add interior points without sorting
    n_tmp = size(tmp, 2)
    if(allocated(tmp) .and. (n_tmp >= 1) ) then 
       call add2list(xy, tmp, tol, xy_tmp, indx, reordered)
       deallocate(tmp) 
    end if

    ! final check ups to see the "indx" matches 
    ! with the size of original array
    if ( (indx - 1) .ne. n ) then
       print *, 'fatal error in fem_reorder: not all points' &
            , '  are reordered or additional points are reordered! stop'
       stop
    end if

    ! return the sorted coordinates
    xy = xy_tmp

    ! done here
  end subroutine fem_reorder

  subroutine add2list(xy, xy_favor, tol, xy_out, indx, reordered)
    implicit none
    real*8, dimension(:,:), intent(in) :: xy, xy_favor
    real*8, intent(in) :: tol
    real*8, dimension(:,:), intent(inout) :: xy_out
    integer, intent(inout) :: indx
    logical, dimension(:), intent(inout) :: reordered

    ! local vars
    integer :: i, j, n, n_favor
    logical :: found 

    ! init
    n_favor = size(xy_favor, 2)
    n = size(xy, 2)

    do j = 1, n_favor 
       found = .false.
       do i = 1, n
          if (  (maxval( abs(xy(:, i) - xy_favor(:, j)) ) <= tol) & 
               .and. (.not. reordered(i))  ) then ! found
             found = .true.
             xy_out(:, indx) = xy(:, i)
             indx = indx + 1
             reordered(i) = .true.
             exit
          end if
       end do

       if (.not. found) then
          print *, 'fatal : point (', xy_favor(1, j), xy_favor(2, j), ') in' &
               , '  the master element list not found! stop'
          stop
       end if
    end do


    ! done here
  end subroutine add2list

end module fem_reordering

! ! little tester for sorting and reordering
! program tester
!   use fem_reordering
!   implicit none

!   real*8, dimension(:), allocatable :: x
!   integer, dimension(:), allocatable :: p
!   real*8, dimension(:,:), allocatable :: xy

!   ! test case (A)
!   x = [2.0d0, -1.0d0, 3.0d0, 3.1d0]
!   allocate(p(size(x)))
!   print *, 'before sorting'
!   print *, ' x = [', x, ']'
!   print *, 'sort smaller to larger:'
!   call arash_sort(x, p, '+')
!   print *, 'x(p) = [', x(p), ']'

!   print *, 'sort larger to smaller:'
!   call arash_sort(x, p, '-')
!   print *, 'x(p) = [', x(p), ']'
!   deallocate(p)

!   ! test case (B)
!   x = [1.0d0, 3.0d0, 4.0d0, 8.1d0]
!   allocate(p(size(x)))
!   print *, 'before sorting'
!   print *, ' x = [', x, ']'
!   print *, 'sort smaller to larger:'
!   call arash_sort(x, p, '+')
!   print *, 'x(p) = [', x(p), ']'
!   print *, ' p = [', p, ']'

!   print *, 'sort larger to smaller:'
!   call arash_sort(x, p, '-')
!   print *, 'x(p) = [', x(p), ']'
!   print *, ' p = [', p, ']'

!   deallocate(p)

!   ! case study (A) for FEM reordering
!   allocate(xy(2,6))
!   xy(1,:) = (/ 1.0d0, 0.5d0, 0.0d0, 0.0d0, 0.0d0, 0.5d0 /)
!   xy(2,:) = (/ 0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.5d0, 0.5d0 /)
!   call fem_reorder(xy, 1.0d-15 )
!   print *, 'xy = [', xy, ']'
!   deallocate( xy )  

!   ! case study (B) for FEM reordering
!   allocate(xy(2,3))
!   xy(1,:) = (/ 1.0d0, 0.0d0, 0.0d0 /)
!   xy(2,:) = (/ 0.0d0, 0.0d0, 1.0d0 /)
!   call fem_reorder(xy, 1.0d-15 )
!   print *, 'xy = [', xy, ']'
!   deallocate( xy )  

!   ! case study (C) for FEM reordering
!   allocate(xy(2,15))
!   xy(1,:) = (/ 0.0d0, 0.25d0, 0.5d0, 0.75d0, 1.0d0, 0.00d0, 0.25d0, 0.50d0, 0.75d0, 0.0d0, 0.25d0, 0.5d0, 0.00d0, 0.25d0, 0.0d0 /)
!   xy(2,:) = (/ 0.0d0, 0.00d0, 0.0d0, 0.00d0, 0.0d0, 0.25d0, 0.25d0, 0.25d0, 0.25d0, 0.5d0, 0.50d0, 0.5d0, 0.75d0, 0.75d0, 1.0d0 /)
!   call fem_reorder(xy, 1.0d-15 )
!   print *, 'xy = [', xy, ']'
!   deallocate( xy )  

!   ! done here
! end program tester
