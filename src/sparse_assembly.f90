module sparse_assembly
  use grid_opt
  use element_opt
  use sparse_opt
  implicit none
  private



  public :: sparse_assem, sparse_assem_M

contains

  subroutine sparse_assem(grd, elems, A)
    implicit none
    type(grid), intent(in) :: grd
    type(element), dimension(:), intent(inout), target :: elems
    class(smat), intent(inout) :: A

    ! local vars
    integer :: i

    ! HARD reset
    A%a = 0.0d0

    ! assemble stiffness matrix
    do i = 1, size(elems)
       elems(i)%A => elems(i)%K 
       call partial_assem(grd, i, elems(i), A)
    end do


    ! add more partial accumulations
    ! including mass matrix, damping matrix and etc ...
    ! HERE :
    ! :
    ! : 
    ! :
    !

    ! done here
  end subroutine sparse_assem

  subroutine sparse_assem_M(grd, elems, A)
    implicit none
    type(grid), intent(in) :: grd
    type(element), dimension(:), intent(inout), target :: elems
    class(smat), intent(inout) :: A

    ! local vars
    integer :: i

    ! HARD reset
    A%a = 0.0d0

    ! assemble mass matrix
    do i = 1, size(elems)
       elems(i)%A => elems(i)%M 
       call partial_assem(grd, i, elems(i), A)
    end do


    ! done here
  end subroutine sparse_assem_M

  ! assembles and accumulates one element 
  ! into global sparse matrix A
  subroutine partial_assem(grd, ielem, elem, A)
    implicit none
    type(grid), intent(in) :: grd
    integer, intent(in) :: ielem
    type(element), intent(inout) :: elem
    class(smat), intent(inout) :: A

    ! local vars
    integer :: i, j, pti, ptj 
    integer :: npe, neqs
    real*8, dimension(size(elem%A, 1), size(elem%A, 1)) :: tmp

    ! init
    npe = size(elem%A, 3)
    neqs = size(elem%A, 1)

    do i = 1, npe
       pti = grd%icon(ielem, i)
       do j = 1, npe
          ptj = grd%icon(ielem, j)

          tmp = A%get(pti, ptj)
          tmp = tmp + elem%A(:, :, i, j)
          call A%set(pti, ptj, tmp)
       end do
    end do


    ! done here
  end subroutine partial_assem

end module sparse_assembly
