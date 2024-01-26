module element_opt
  ! includes data types for abstract representation of each
  ! element and various subroutines to compute
  ! and manipulate elements
  use grid_opt
  use gen_basis
  use dunavant
  use approx_fekete, only : fekete

  implicit none

  private

  ! a datatype abstracting an element
  type element

     ! header ----
     integer :: neqs
     ! the following are local coords of Gauss points.
     ! r(k), s(k), W(k); <k = 1...ngauss> are coords and
     ! weights of Gauss-Legendre points on a triangle
     integer :: ngauss
     real*8, dimension(:), allocatable :: r, s, W
     ! generalized basis functions
     type(basis) :: tbasis
     ! psi and d_psi have dimensions = (1..npe, 1..ngauss)
     real*8, dimension(:,:), allocatable :: psi, d_psi_d_xi, d_psi_d_eta
     ! Jacobian of transformation <jac> and its inverse called Jstar
     ! have the dimension(1..2, 1..2, 1..ngauss)
     real*8, dimension(:,:,:), allocatable :: jac, Jstar
     ! the value of Jacobian of the element at each interpolation point
     real*8, dimension(:), allocatable :: JJ

     integer :: init_lock ! would be 1221360 if initialized!
     ! -----

     ! K(neqs, neqs, npe, npe)
     real*8, dimension(:,:,:,:), allocatable :: K, M
     ! Q,f(neqs, npe)
     real*8, dimension(:,:), allocatable :: Q, f

     ! generic pointer array
     ! used for matrix-free Ax
     real*8, dimension(:,:,:,:), pointer :: A => null()

  end type element

  public :: element, init_elem, comp_elem_KQf, assemble_all
  public :: comp_u_point, comp_grad_u_point, xy2rs
  public :: comp_Jstar_point

contains

  ! initializes each element according to the
  ! grid specification
  subroutine init_elem(elem, ielem, grd, neqs)
    implicit none
    type(element), intent(inout) :: elem
    integer, intent(in) :: ielem ! element number
    type(grid), intent(inout) :: grd
    integer, intent(in) :: neqs

    ! local vars
    integer :: i, p, rule, order_num, npe, degree
    real*8 :: x0, y0
    real*8, dimension(:,:), allocatable :: xy
    type(fekete)  :: tfekete

    if( elem%init_lock .eq. 1221360 ) then
       print *, 'error : the element #', ielem, ' is already initialized! stop'
       stop
    end if

    p = grd%p(ielem)
    npe = grd%npe(ielem)

    ! change logic here for general mixed element
    ! higher order elements
    ! -------------------------------------------
    elem%neqs = neqs

    ! compute the "order_num" or number of Gauss points
    ! rule = maxval((/ ceiling(dble(p+5)/dble(2) ), p /))
    rule = int(2.0d0 * p) ! infact it should be "p" for linear
                 ! Laplace equation for constant shape elements
                 ! but since we've got rational function
                 ! for curvilinear elements, then we put
                 ! it 2*p for safety.
                 ! NOTE : in general this is not correct
                 !
    ! rule = p ! infact it should be "p" for linear

    select case (grd%elname(ielem))

    case (GEN_TRIANGLE)
       ! check to see the order of exactness is available
       ! in the tables for the given rule
       call dunavant_degree ( rule, degree )

       ! compute the number of required Gauss points
       call dunavant_order_num( rule, order_num )
       elem%ngauss = order_num

       ! allocate space for that
       allocate( xy(2,order_num))
       allocate(elem%r(order_num), elem%s(order_num), elem%W(order_num))

       ! compute the absicca and weights for that rule
       call dunavant_rule( rule, order_num, xy, elem%W )
       elem%r = xy(1,:)
       elem%s = xy(2,:)
       deallocate(xy)

    case (GEN_QUADRI)

       call grd%tfekete_table%lookup(d = rule, name = 'quadri', fekete_out = tfekete &
            , spacing = 'equal_space', s = 3, echo = .true.)

       ! call tfekete%init(d = rule, name = 'quadri' &
       !      , spacing = 'equal_space', s = 3, echo = .true.)
       elem%ngauss = size(tfekete%w)

       ! allocate space for that
       allocate(elem%r(elem%ngauss), elem%s(elem%ngauss), elem%W(elem%ngauss))

       elem%r = tfekete%fin_coords(1, :)
       elem%s = tfekete%fin_coords(2, :)
       elem%W = tfekete%w

       ! ! deallocate stuff in tfekete object
       ! call tfekete%clean()

    case default
       print *, 'unknown elname in obtaining quadrature rules! stop'
       stop

    end select

    ! initializing the basis functions and their derivatives
    call elem%tbasis%init(grd%maselem(ielem)%xi, grd%maselem(ielem)%eta &
         , grd%elname(ielem))

    allocate(elem%psi(npe,elem%ngauss) &
           , elem%d_psi_d_xi(npe,elem%ngauss) &
           , elem%d_psi_d_eta(npe,elem%ngauss) )

    ! evaluating the basis function and their derivatives at
    ! Gauss - Legendre (quadrature) points and then storing
    do i = 1, elem%ngauss
       x0 = elem%r(i)
       y0 = elem%s(i)
       call elem%tbasis%eval(x0, y0, 0,  elem%psi(:,i)        )
       call elem%tbasis%eval(x0, y0, 1,  elem%d_psi_d_xi(:,i) )
       call elem%tbasis%eval(x0, y0, 2,  elem%d_psi_d_eta(:,i))
    end do

    ! allocate Jacobian of transformation
    allocate(elem%jac(2,2, elem%ngauss), elem%Jstar(2,2, elem%ngauss) )
    allocate(elem%JJ(elem%ngauss) )

    ! compute the Jacobian of transformation matrix , its inverse and
    ! the value of Jacobian at each Gauss-Legendre point
    do i = 1, elem%ngauss
       call comp_Jstar(grd, elem, ielem, i, elem%jac(:,:,i) &
            , elem%Jstar(:,:,i), elem%JJ(i) )
    end do

    ! initialize the stiffness matrix
    allocate(elem%K(elem%neqs,elem%neqs, npe, npe))
    allocate(elem%M(elem%neqs,elem%neqs, npe, npe))
    elem%K = 0.0d0
    elem%M = 0.0d0

    ! allocating and initializing rhs
    allocate(elem%Q(elem%neqs,npe))
    allocate(elem%f(elem%neqs,npe))
    elem%Q = 0.0d0
    elem%f = 0.0d0

    ! initialization terminated successfully!
    elem%init_lock = 1221360 !locked
    ! -------------------------------------------

    ! done here
  end subroutine init_elem


  ! computes Jacobian of the transformation "jac" and
  ! "Jstar" or the inverse of the Jacobian
  ! of the transformation at the Gauss point "(r(k),s(k))"
  ! within the element "elem" with number "ielem" in grid "grd".
  ! The determinant of Jacobian of the transformation
  ! is returned in "JJ".

  subroutine comp_Jstar(grd, elem, ielem, k, jac, Jstar, JJ)
    implicit none
    type(grid), intent(in) :: grd
    type(element), intent(in) :: elem
    integer, intent(in) :: ielem, k
    real*8, dimension(:,:), intent(out) :: jac, Jstar
    real*8, intent(out) :: JJ

    ! local vars
    integer :: i
    real*8 :: xi, yi
    real*8, dimension(2) :: der
    real*8 :: dx_dr, dx_ds, dy_dr, dy_ds
    integer :: npe

    ! hard reset
    jac = 0.0d0; Jstar = 0.0d0; JJ = 0.0d0
    xi = 0.0d0; yi = 0.0d0; der = 0.0d0
    dx_dr = 0.0d0; dx_ds = 0.0d0; dy_dr= 0.0d0; dy_ds = 0.0d0
    npe = grd%npe(ielem)
    ! print *, '============================================='
    ! print *, '============================================='
    ! print *, 'computing jacobian for elem #', ielem
    ! print *, '============================================='
    ! print *, '============================================='

    ! compute the components of jac
    do i = 1, npe ! assuming iso-geometric expansion for x, y
       ! print *, 'at point #', grd%icon(ielem,i)
       der(1) =  elem%d_psi_d_xi(i,k)
       der(2) = elem%d_psi_d_eta(i,k)
       ! print *, 'der = [', der, ']'
       xi = grd%x(grd%icon(ielem,i))
       yi = grd%y(grd%icon(ielem,i))
       ! print *, 'xi = ', xi, 'yi = ', yi
       dx_dr = dx_dr + xi * der(1)
       dx_ds = dx_ds + xi * der(2)
       dy_dr = dy_dr + yi * der(1)
       dy_ds = dy_ds + yi * der(2)
       ! print *, ' total dx_dr = ', dx_dr
       ! print *, ' total dx_ds = ', dx_ds
       ! print *, ' total dy_dr = ', dy_dr
       ! print *, ' total dy_ds = ', dy_ds

    end do


    ! compute jac
    jac(1,1) = dx_dr; jac(1,2) = dy_dr
    jac(2,1) = dx_ds; jac(2,2) = dy_ds

    ! comp JJ, i.e. det of jac
    JJ = dx_dr * dy_ds - dx_ds * dy_dr

    ! check it, should be valid grid!!!
    if (JJ <= 0.0d0) then
       print *, 'error : negative or zero Jacobian at element (',ielem,')'
       ! print *, ' JJ = ', JJ
       stop
    ! else
    !    print *, 'OKKKK : positive Jacobian at element (',ielem,')'
    !    print *, ' JJ = ', JJ
    end if

    ! comp inverse of jac and store it in Jstar
    Jstar(1,1) =  1.0d0 / JJ * jac(2,2)
    Jstar(2,2) =  1.0d0 / JJ * jac(1,1)
    Jstar(1,2) = -1.0d0 / JJ * jac(1,2)
    Jstar(2,1) = -1.0d0 / JJ * jac(2,1)

    ! done here
  end subroutine comp_Jstar


  !> @brief computes mass matrix [K] and rhs vectors
  !! for the given element
  !> @note the element must be initialized before this.

  subroutine comp_elem_KQf(elem, ielem, grd)
    implicit none
    type(element), intent(inout) :: elem
    integer, intent(in) :: ielem ! element number
    type(grid), intent(in) :: grd

    ! local vars
    integer :: i, j, k, l, m, npe
    real*8, dimension(2) :: der1, der2
    real*8 :: coeff

    ! check this element must be initialized before
    if (elem%init_lock .ne. 1221360) then
       print *, 'element ', ielem,' is not initialized! stop.'
       stop
    end if

    npe = grd%npe(ielem)
    ! print *, 'the area of elem #', ielem , sum(elem%JJ * elem%W)

    ! HARD reset
    elem%K = 0.0d0
    elem%M = 0.0d0

    ! select the coefficient of the Jacobian of the transformation
    select case (grd%elname(ielem))
    case (GEN_QUADRI)
       coeff = 1.0d0
    case (GEN_TRIANGLE)
       coeff = 0.5d0
    end select

    ! fill out stiffness matrix
    do l = 1, elem%neqs
       do m = 1, elem%neqs
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
                   elem%K(l,m,i,j) = elem%K(l,m,i,j) &
                        + sum(der1 * der2) * coeff * elem%JJ(k) * elem%W(k)
                   ! accumulate to the mass matrix of this element
                   elem%M(l,m,i,j) = elem%M(l,m,i,j) &
                        + elem%psi(i,k) * elem%psi(j,k) * coeff * elem%JJ(k) * elem%W(k)

                end do
             end do
          end do
       end do
    end do

    ! HARD reset
    elem%Q = 0.0d0
    elem%f = 0.0d0

    ! done here
  end subroutine comp_elem_KQf

  !< @brief
  !< assembles the entire stiffness matrix
  subroutine assemble_all( elems, grd, KK, rhs)
    implicit none
    type(element), dimension(:), target, intent(in) :: elems
    type(grid), target, intent(in) :: grd
    real*8, dimension(:,:,:,:), intent(out) :: KK
    real*8, dimension(:,:), intent(out) :: rhs

    ! local vars
    integer :: k, i, j, npe
    integer, dimension(:,:), pointer :: icon
    type(element), pointer :: elem
    integer :: pti, ptj

    ! load connectivity matrix
    icon => grd%icon

    ! hard reset
    KK = 0.0d0
    rhs = 0.0d0

    ! start assembling
    do k = 1, grd%ncellsg ! loop over all cells
       elem => elems(k) ! load this element
       npe = grd%npe(k)
       do i = 1, npe
          pti = icon(k,i)
          rhs(:,pti) = rhs(:,pti) + elem%Q(:, i) + elem%f(:, i)

          do j = 1, npe
             ptj = icon(k,j)
             KK(:,:,pti, ptj) = KK(:,:,pti, ptj) + elem%K(:,:,i,j)
          end do

       end do

    end do

    ! done here
  end subroutine assemble_all

  !< @details
  !< computes Jacobian of the transformation "jac" and
  !! "Jstar" or the inverse of the Jacobian
  !! of the transformation at any point "(r,s)"
  !! within the element "elem" with number "ielem" in grid "grd".
  !< The determinant of Jacobian of the transformation
  !! is returned in "JJ".

  subroutine comp_Jstar_point(grd, elem, ielem, r, s, jac, Jstar, JJ)
    implicit none
    type(grid), intent(in) :: grd
    type(element), intent(inout) :: elem
    integer, intent(in) :: ielem
    real*8, intent(in) :: r, s
    real*8, dimension(:,:), intent(out) :: jac, Jstar
    real*8, intent(out) :: JJ

    ! local vars
    integer :: i
    real*8 :: xi, yi
    real*8, dimension(2), save :: der
    real*8 :: dx_dr, dx_ds, dy_dr, dy_ds
    integer :: npe
    real*8, dimension(grd%npe(ielem)) :: d_psi_d_xi, d_psi_d_eta

    ! hard reset
    jac = 0.0d0; Jstar = 0.0d0; JJ = 0.0d0
    xi = 0.0d0; yi = 0.0d0; der = 0.0d0
    dx_dr = 0.0d0; dx_ds = 0.0d0; dy_dr= 0.0d0; dy_ds = 0.0d0
    npe = grd%npe(ielem)

    ! computing the derivative of basis function at that point
    call elem%tbasis%eval(r, s, 1, d_psi_d_xi)
    call elem%tbasis%eval(r, s, 2, d_psi_d_eta)

    ! compute the components of jac
    do i = 1, npe ! assuming iso-geometric expansion for x, y

       der(1) =  d_psi_d_xi(i)
       der(2) = d_psi_d_eta(i)
       xi = grd%x(grd%icon(ielem,i))
       yi = grd%y(grd%icon(ielem,i))

       dx_dr = dx_dr + xi * der(1)
       dx_ds = dx_ds + xi * der(2)
       dy_dr = dy_dr + yi * der(1)
       dy_ds = dy_ds + yi * der(2)

    end do


    ! compute jac
    jac(1,1) = dx_dr; jac(1,2) = dy_dr
    jac(2,1) = dx_ds; jac(2,2) = dy_ds

    ! comp JJ, i.e. det of jac
    JJ = dx_dr * dy_ds - dx_ds * dy_dr

    ! check it, should be valid grid!!!
    if (JJ <= 0.0d0) then
       print *, 'error in comp_Jstar_point: negative or zero' &
            , ' Jacobian at element (',ielem,')'
       stop
    end if

    ! comp inverse of jac and store it in Jstar
    Jstar(1,1) =  1.0d0 / JJ * jac(2,2)
    Jstar(2,2) =  1.0d0 / JJ * jac(1,1)
    Jstar(1,2) = -1.0d0 / JJ * jac(1,2)
    Jstar(2,1) = -1.0d0 / JJ * jac(2,1)

    ! done here
  end subroutine comp_Jstar_point

  !
  !< @details
  !< computes the primary variable <u> at a
  !! node of an element using the values of the basis function
  !! given at that point
  !
  ! computes the primary variable <u> at a  
  ! node of an element using the values of the basis function
  ! given at that point
  !
  !< i.e. for element "i" at point (r,s) we have:
  !
  !<    u(r,s) = sum_k (u_k psi_k(r, s) )
  !
  !

  subroutine comp_u_point(u, grd, elem, ielem, r, s, u_out)
    implicit none
    real*8, dimension(:,:), intent(in) :: u
    type(grid), intent(in) :: grd
    type(element), intent(inout) :: elem
    integer, intent(in) :: ielem
    real*8, intent(in) :: r, s
    ! dimension(elem%neqs)
    real*8, dimension(:), intent(out) :: u_out

    ! local vars
    integer :: k, ptk, npe
    real*8, dimension(grd%npe(ielem)) :: psi

    ! hard reset
    u_out = 0.0d0
    npe = grd%npe(ielem)

    call elem%tbasis%eval(r, s, 0, psi)

    do k = 1, npe
       ptk = grd%icon(ielem,k)
       u_out(:) = u_out(:) + u(:,ptk) * psi(k)
    end do

    ! done here
  end subroutine comp_u_point

  !
  !< @details
  !! Computes the gradient at an interior
  !! node of an element using the values of
  !! the gradient of the basis functions
  !! given at that point
  !
  !
  !< i.e. for element "i" at point (r,s) we have:
  !!
  !!   /f$ \grad u(r,s) = sum_{k} (u_{k} * \grad \psi_{k(r, s)}) /f$
  !
  !

  subroutine comp_grad_u_point(u, grd, elem, ielem, r, s, dudx, dudy)
    implicit none
    real*8, dimension(:,:), intent(in) :: u !< The value of solution
    type(grid), intent(in) :: grd   !< Grid under calculation
    type(element), intent(inout) :: elem  !< Current element
    integer, intent(in) :: ielem  !< Current element number in the grid.
    real*8, intent(in) :: r, s    !< Point coordinates in the master element
    ! dimension(elem%neqs)
    !< @todo This may need to be adjusted for 2D -3D cases or to be generalized.
    real*8, dimension(:), intent(out) :: dudx, dudy !< Gradients of U in x and y directions.

    ! local vars
    integer :: k, ptk, npe
    real*8 :: JJ
    real*8, dimension(2), save :: der
    real*8, dimension(2,2), save :: jac, Jstar
    real*8, dimension(grd%npe(ielem)) :: d_psi_d_xi, d_psi_d_eta

    ! hard reset
    dudx = 0.0d0
    dudy = 0.0d0
    npe = grd%npe(ielem)


    call comp_Jstar_point(grd, elem, ielem, r, s, jac, Jstar, JJ)
    call elem%tbasis%eval(r, s, 1, d_psi_d_xi)
    call elem%tbasis%eval(r, s, 2, d_psi_d_eta)

    do k = 1, npe
       ptk = grd%icon(ielem,k)

       der(1) = d_psi_d_xi(k)
       der(2) = d_psi_d_eta(k)
       der = matmul(Jstar, der)
       dudx(:) = dudx(:) + u(:,ptk) * der(1)
       dudy(:) = dudy(:) + u(:,ptk) * der(2)
    end do

    ! done here
  end subroutine comp_grad_u_point

  !< @details
  !! computes local coords (r,s) given the
  !! global coordinates (x,y) using a robust
  !! Newton method.
  !<
  !! This should be applied consistently to
  !! higher-order elements also.

  subroutine xy2rs(x, y, elem, ielem, grd, maxitr, tolrs, r, s)
    implicit none
    real*8, intent(in) :: x, y
    type(element), intent(inout) :: elem
    integer, intent(in) :: ielem
    type(grid), intent(in) :: grd
    integer, intent(in) :: maxitr
    real*8, intent(in) :: tolrs
    real*8, intent(out) :: r, s

    ! local vars
    integer :: i, j, npt, pt ! points per triangle
    real*8 :: val
    real*8, dimension(grd%npe(ielem)) :: psi, d_psi_d_xi, d_psi_d_eta
    real*8, dimension(2), save :: der, F, delta
    real*8, dimension(2,2), save:: JJ, Jinv


    ! initial guess
    r = 0.25d0; s =0.25d0
    npt = grd%npe(ielem)

    do j = 1, maxitr ! main Newton iteration loop

       ! setting up Newton functional
       F = (/ x, y /)

       ! evaluating basis functions
       call elem%tbasis%eval(r, s, 0,  psi        )
       call elem%tbasis%eval(r, s, 1,  d_psi_d_xi )
       call elem%tbasis%eval(r, s, 2,  d_psi_d_eta)

       ! HARD reset
       JJ = 0.0d0; Jinv = 0.0d0

       do i = 1, npt

          pt = grd%icon(ielem, i)
          val = psi(i)
          der(1) = d_psi_d_xi(i); der(2) = d_psi_d_eta(i)

          F = F - (/ (grd%x(pt) * val) , (grd%y(pt) * val) /)
          JJ(1,1) =  JJ(1,1) + grd%x(pt) * der(1)
          JJ(1,2) =  JJ(1,2) + grd%x(pt) * der(2)
          JJ(2,1) =  JJ(2,1) + grd%y(pt) * der(1)
          JJ(2,2) =  JJ(2,2) + grd%y(pt) * der(2)

       end do

       Jinv(1,1) = JJ(2,2)
       Jinv(2,2) = JJ(1,1)
       Jinv(2,1) = -JJ(2,1)
       Jinv(1,2) = -JJ(1,2)
       Jinv = 1.0d0 /(JJ(1,1) * JJ(2,2) - JJ(2,1) * JJ(1,2)) * Jinv
       delta = matmul(Jinv, F)

       r = r + delta(1)
       s = s + delta(2)

       if ( sqrt(sum(delta*delta)) .le. tolrs ) then
          ! robust bug checking before returning
          if ( (npt .eq. 3) .and. (j > 2) ) then
             print *, 'number of itr is ', j, ' and is greater than 2 ' &
                    , 'for linear element with number', ielem, '! stop.'
             stop
          else if (  (r > (1.0d0 + tolrs)) .or. (s > (1.0d0 + tolrs)) &
               .or.  (r < (0.0d0 - tolrs)) .or. (s < (0.0d0 - tolrs)) ) then
             ! print *, '(r,s) = ', r, s ,' out of range! stop.'
             ! stop
          end if
          ! print *, 'j = ', j , 'error = ', sqrt(sum(delta*delta))
          return
       end if

    end do

    print *, 'xy2rs(...) did not converge' &
         , ' to desired tolerance in ', maxitr, 'iterations' &
         , ' at point (', x, y, ')! stop.'
    stop

    ! done here
  end subroutine xy2rs

end module element_opt
