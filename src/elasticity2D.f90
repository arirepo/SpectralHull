module elasticity2D
  ! includes data types for abstract representation of each
  ! element and various subroutines to compute
  ! and manipulate elements
  use grid_opt
  use element_opt
  use solidMat

  implicit none

  private

  public :: comp_2D_elem_K

contains

  !> @brief Computes the stiffness matrix [K]
  !! for a given 2D element and plane stress(0) or plane strain (1)
  !! condition.
  !> @note the element must be initialized before this.

  subroutine comp_2D_elem_K(elem, ielem, grd, mater)
    implicit none
    type(element), intent(inout) :: elem
    integer, intent(in) :: ielem ! element number
    type(grid), intent(in) :: grd
    type(isoMat2D), intent(in) :: mater

    ! local vars
    integer :: i, j, k, l, m, npe
    real*8, dimension(2) :: der1, der2
    real*8 :: coeff, C11,C12, C21, C22, C66

    ! check this element must be initialized before
    if (elem%init_lock .ne. 1221360) then
      print *, 'Element ', ielem,' is not initialized! stop.'
      stop
    end if

    ! check this material must be initialized before
    if (mater%init_lock .ne. 1221360) then
      print *, 'Material is not initialized! stop.'
      stop
    end if

    npe = grd%npe(ielem)

    C11 = mater%C(1)
    C12 = mater%C(2)
    C21 = C12
    C22 = mater%C(3)
    C66 = mater%C(4)

    ! HARD reset
    elem%K = 0.0d0

    ! select the coefficient of the Jacobian of the transformation
    select case (grd%elname(ielem))
      case (GEN_QUADRI)
        coeff = 1.0d0
      case (GEN_TRIANGLE)
        coeff = 0.5d0
    end select

    ! fill out stiffness matrix
    ! do l = 1, elem%neqs
    ! do m = 1, elem%neqs
    do i = 1, npe
      do j = 1, npe
        do k = 1, elem%ngauss
          ! get grad of basis functions
          der1(1) = elem%d_psi_d_xi(i,k); der1(2) = elem%d_psi_d_eta(i,k)
          der2(1) = elem%d_psi_d_xi(j,k); der2(2) = elem%d_psi_d_eta(j,k)
          ! transform computational grads to physical grads
          der1 = matmul(elem%Jstar(:,:,k), der1)
          der2 = matmul(elem%Jstar(:,:,k), der2)
          ! accumulate to the stiffness matrix of this element
          !elem%K(l,m,i,j) = elem%K(l,m,i,j) &
            !     + sum(der1 * der2) * coeff * elem%JJ(k) * elem%W(k)

          ! accumulate to the stiffness matrix of this element
          elem%K(1,1,i,j) = elem%K(1,1,i,j) &
            + ( C11 * der1(1) * der2(1) + C66 * der1(2) * der2(2) ) * coeff * elem%JJ(k) * elem%W(k)

          elem%K(1,2,i,j) = elem%K(1,2,i,j) &
            + ( C12 * der1(1) * der2(2) + C66 * der1(2) * der2(1) ) * coeff * elem%JJ(k) * elem%W(k)

          elem%K(2,1,i,j) = elem%K(2,1,i,j) &
            + ( C12 * der1(2) * der2(1) + C66 * der1(1) * der2(2) ) * coeff * elem%JJ(k) * elem%W(k)

          elem%K(2,2,i,j) = elem%K(2,2,i,j) &
            + ( C22 * der1(2) * der2(2) + C66 * der1(1) * der2(1) ) * coeff * elem%JJ(k) * elem%W(k)

        end do
      end do
    end do
    ! end do
    ! end do

    ! done here
  end subroutine comp_2D_elem_K

end module elasticity2D
