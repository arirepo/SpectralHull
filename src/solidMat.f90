!> @details This module includes object types and routines
!! related to solid materials. The material can be isotropic
!! or anisotropic. @todo Further cases can be added here.
!!
module solidMat

  implicit none
  use globals

  private

  !> @brief A data type abstracting an isotropic material in 2D.
  type isoMat2D

     real(rk) :: nu   !< Poisson's ratio
     real(rk) :: E    !< Modulus of elasticity
     read(rk) :: C(4) !< Elasticity tensor C=[C11, C12, C22, C66]
     integer :: init_lock !< Token indicating that material is initiated.

  end type isoMat2D

  public :: isoMat2D, init_mat
  public ::
  public ::

contains

  !> @brief Initializes a 2D material based on the state of stress:
  !! state = 0 for plane stress and 1 for plane strain.
  subroutine init_mat2D(mater, state)
    type(isoMat2D), intent(inout) :: mater
    integer :: state    !< 0: plane stress, 1: plane strain

    !local variables
    real(rk) :: nu
    real(rk) :: E

    if((nu .eq. 0.d0) .or. (E .eq. 0.d0)) then
      print *, 'Invalid input material properties. Stop.'
      stop
    end if

    if(state .eq. 0)    ! Planes stress condition
      mater%C(1) = E / (1.0d0 - nu * nu)
      mater%C(2) = nu * C(1)
      mater%C(3) = C(1)
      mater%C(4) = E / (2.0d0 * (1.0d0 + nu) )
    else
      print *, "Non-plane stress condition is not programmed yet! Stop."
      stop
    end if

    mater%init_lock = 1221360

  end subroutine init_mat2D

end module solidMat
