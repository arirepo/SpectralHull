module reorder_opt
  use grid_opt
  use fem_reordering
  use rev_cuthill_mckee
  use sparse_opt

  implicit none

  private

  type reorder
     type(ints), dimension(:), allocatable :: rows
     integer, private :: node_num, adj_num
     integer, dimension(:), allocatable, private :: adj_row, adj, perm, perm_inv
     logical, private :: pattern_found = .false., already_reordered = .false.

   contains
     procedure, public :: init => init_reorder
     procedure, public :: spy => write_sparsity_pattern
     procedure, public :: reorder => reorder_grid 
  end type reorder

  public :: reorder

contains

  subroutine find_sparsity_pattern(grd, rows)
    implicit none
    type(grid), intent(in) :: grd
    type(ints), dimension(:), allocatable :: rows

    ! local vars
    integer :: i, j, k, ptj, ptk, tsize
    integer, dimension(:), allocatable :: itmp, itmp2

    ! refresh the pattern if already exists
    if ( allocated(rows) ) deallocate(rows)
    allocate(rows(grd%nnodesg))

    ! create the pattern
    do i = 1, grd%ncellsg
       do j = 1, grd%npe(i)
          ptj = grd%icon(i, j)
          do k = 1, grd%npe(i)
             ptk = grd%icon(i, k)
             if ( .not. allocated(rows(ptj)%ja) ) then
                allocate(rows(ptj)%ja(1))
                rows(ptj)%ja(1) = ptk
             else
                if ( .not. any( rows(ptj)%ja .eq. ptk ) ) then
                   tsize = size(rows(ptj)%ja)
                   allocate( itmp( tsize + 1) )
                   itmp(1:tsize) = rows(ptj)%ja(1:tsize)
                   itmp(tsize + 1) = ptk
                   call move_alloc(itmp, rows(ptj)%ja)    
                end if
             end if
          end do
       end do
    end do

    ! sort rows according to increasing column index
    do i = 1, grd%nnodesg
       allocate(itmp(size(rows(i)%ja)))
       allocate(itmp2(size(itmp)))
       call arash_sort(dble(rows(i)%ja), itmp, '+')
       itmp2 = rows(i)%ja(itmp) 
       rows(i)%ja = itmp2
       deallocate(itmp, itmp2)
    end do

    ! done here
  end subroutine find_sparsity_pattern

  subroutine write_sparsity_pattern(this, outfile)
    implicit none
    class(reorder), intent(in) :: this
    character(len = *), intent(in) :: outfile

    ! local vars
    integer :: i, j 

    ! bullet proofing ...
    if ( .not. this%pattern_found ) then
       print *, 'no sparsity pattern has been found yet! can not write it! stop.'
       stop
    end if
 
    open(unit = 7, file = outfile, status = 'replace')

    do i = 1, size(this%rows)
       do j = 1, size(this%rows(i)%ja)
          write(7, *) i, this%rows(i)%ja(j), '1.0'
       end do
    end do

    close(7)

    ! done here
  end subroutine write_sparsity_pattern

  subroutine find_adj(rows, node_num, adj_num, adj_row, adj)
    implicit none
    type(ints), dimension(:), intent(in) :: rows
    integer, intent(out) :: node_num, adj_num
    integer, dimension(:), allocatable :: adj_row, adj

    ! local vars
    integer :: i, j, ii

    ! refresh the storage, if it is already filled
    if ( allocated(adj_row) ) deallocate(adj_row)
    if ( allocated(adj) ) deallocate(adj)

    ! comp total number of nodes
    node_num = size(rows)
    ! compute total number of adj
    adj_num = 0
    do i = 1, node_num
       adj_num = adj_num + (size(rows(i)%ja) - 1) ! minus diagonal
    end do

    ! allocate outputs
    allocate(adj_row(node_num + 1), adj(adj_num))

    ! compute adj_row
    adj_row(1) = 1
    do i = 2, (node_num+1)
       adj_row(i) = adj_row(i-1) + size(rows(i-1)%ja) - 1
    end do

    ! compute adj
    ii = 1
    do i = 1, node_num
       do j = 1, size(rows(i)%ja)
          if ( i .eq. rows(i)%ja(j)) cycle

          adj(ii) = rows(i)%ja(j)
          ii = ii + 1
       end do
    end do

    ! done here
  end subroutine find_adj
  

  subroutine init_reorder(this, grd)
    implicit none
    class(reorder), intent(inout) :: this
    type(grid), intent(in) :: grd

    ! local vars
    integer :: i

    ! init.=- HARD reset double ckeck to make sure
    this%pattern_found = .false.
    this%already_reordered = .false.

    ! first capture sparsity pattern from input grid
    call find_sparsity_pattern(grd, this%rows)
    ! lock it
    this%pattern_found = .true.

    ! then find adjacency info and store in this object
    call find_adj(this%rows, this%node_num, this%adj_num, this%adj_row, this%adj)

    ! treat the rest of allocatable vars
    if (allocated(this%perm)) deallocate(this%perm)
    if (allocated(this%perm_inv)) deallocate(this%perm_inv)
    allocate(this%perm(this%node_num), this%perm_inv(this%node_num))

    ! init permutation and inverse of permutation vectors
    this%perm = (/ (i, i = 1, this%node_num) /)
    this%perm_inv = this%perm
    
    ! done here
  end subroutine init_reorder

  ! reorders a simple integer array <a>
  ! using inverse permutation vector <pinv>
  !
  ! a = an integer vector of arbitrary size
  !     containing a set of original nodes 
  !     with original numbers between 1 and grd%nnodesg.
  !
  ! pinv = an integer array of length (1:grd%nnodesg)
  !        containing the inverse of permutation
  !        obtained by applying reverse cuthill-mckee
  !        algorithm.
  ! 
  subroutine reorder_array(a, pinv)
    implicit none
    integer, dimension(:), intent(inout) :: a
    integer, dimension(:), intent(in) :: pinv

    ! local vars
    integer :: i

    do i = 1, size(a)
       a(i) = pinv(a(i))
    end do

    ! done here
  end subroutine reorder_array

  ! 
  subroutine reorder_coord(x, pinv)
    implicit none
    real*8, dimension(:), intent(inout) :: x
    integer, dimension(:), intent(in) :: pinv

    ! local vars
    integer :: i
    real*8, dimension(size(x)) :: tmp

    do i = 1, size(x)
       tmp(pinv(i)) = x(i)
    end do

    x = tmp

    ! done here
  end subroutine reorder_coord

  subroutine reorder_grid(this, grd)
    implicit none
    class(reorder), intent(inout) :: this
    type(grid), intent(inout) :: grd

    ! local vars
    integer :: i

    ! bullet proofing ...
    if ( this%already_reordered ) then
       print *, 'the given grid is already reordered! no need for that! stop'
       stop
    end if

    ! then run reverse cuthill-mckee algorithm and store permutation
    print *, 'reordering using reversed Cuthill-Mckee algorithm ...'
    call genrcm( this%node_num, this%adj_num, this%adj_row, this%adj, this%perm )
    print *, 'reordering done!'

    ! compute the inverse of permutation matrix
    call perm_inverse3 ( this%node_num, this%perm, this%perm_inv )

    !
    ! Now, apply reordering to the given grid ...
    !
    ! 1- reorder coordinates
    call reorder_coord(grd%x, this%perm_inv)
    call reorder_coord(grd%y, this%perm_inv)

    ! 2- reorder icon (connectivity matrix)
    do i = 1, size(grd%icon, 1)
       call reorder_array(grd%icon(i, :), this%perm_inv)
    end do

    ! 3- reorder ibedge  
    do i = 1, size(grd%ibedge, 1)
       call reorder_array(grd%ibedge(i, :), this%perm_inv)
    end do

    ! 4- reorder duplicate nodes (if available)
    if( allocated(grd%dup_nodes) ) then
       if ( size(grd%dup_nodes) > 0 ) then
          call reorder_array(grd%dup_nodes, this%perm_inv)
       end if
    end if

    ! 5- reorder el2edg (nodes on edges of element)
    do i = 1, size(grd%icon, 1)
       if (allocated(grd%el2edg(i)%edg1)) then
          if ( size(grd%el2edg(i)%edg1) > 0 ) then
             call reorder_array(grd%el2edg(i)%edg1, this%perm_inv)
             call reorder_array(grd%el2edg(i)%edg2, this%perm_inv)
             call reorder_array(grd%el2edg(i)%edg3, this%perm_inv)
          end if
       end if
    end do

    ! reinitialize the newly reordered grid
    call this%init(grd)

    ! finally, lock it
    this%already_reordered = .true.

    ! done here
  end subroutine reorder_grid

end module reorder_opt
