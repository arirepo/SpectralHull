module element_opt_dg2d
  use element_opt
  use grid_opt
  use euler2d_eqs
  implicit none

  private

  type neigh_dg
     ! the neighbor element number and number of 
     !gauss points per this common segment
     integer :: elnum, ngpseg
     ! the element-wise location of 1d Gauss legendre quad rule on the segment
     ! dim = 1:ngpseg
     real*8, dimension(:), allocatable :: r_in, s_in, r_out, s_out, W
     ! Riemann flux on the edge stored at Gauss points (1:neqs, 1:ngpseg) 
     real*8, dimension(:, :), allocatable :: Fstar
  end type neigh_dg

  type edg_dg
     integer :: tag !can be used for BCs
     integer, dimension(:) , allocatable :: pts
     type(neigh_dg), dimension(:) , allocatable :: neighs
  end type edg_dg

  type, extends(element) :: element_dg2d
     private
     integer :: number, npe, elname, p, eltype, npedg, nedgs
     ! (1..2, 1..npe) (1, npe) = x <> (2, npe) = y
     real*8, dimension(:, :), allocatable :: x
     real*8 :: gamma
     ! (neqs * npe, neqs * npe)
     real*8, dimension(:, :), allocatable :: Mass
     ! elemental solution : Uij (i=1,neqs <> j=1,npe) 
     real*8, dimension(:, :), allocatable :: U
     ! physical derivatives of basis functions at Gauss points
     ! (1..npe, 1..ngauss, 1..2)
     real*8, dimension(:, :, :), allocatable :: d_psi_d_x
     ! fluxes evaluated at Gauss points (1..neqs, 1..ngauss, 1..2) 
     real*8, dimension(:, :, :), allocatable :: Fk
     real*8 :: coeff

     !
     type(edg_dg), dimension(:), allocatable :: edgs

  contains

     procedure, nopass, public :: init => init_elem_dg2d
     procedure, public :: comp_mass => comp_mass_mat_dg2d
     procedure, public :: comp_u => comp_u_point_dg
     procedure, public :: comp_flux => comp_flux_interior
     procedure, public :: comp_metric_jacobian => comp_Jstar_point_dg
     procedure, private :: comp_dpsi => comp_d_psi_d_x
     procedure, public :: comp_int_integ => comp_inter_integral

  end type element_dg2d


  public :: element_dg2d


contains
  
  subroutine init_elem_dg2d(elem, ielem, grd, neqs, gamma)
    implicit none
    class(element), intent(inout) :: elem
    integer, intent(in) :: ielem ! element number
    type(grid), intent(inout) :: grd
    integer, intent(in) :: neqs
    real*8, intent(in) :: gamma

    ! local vars
    integer :: npe, ii, jj, edgnum, local_edg, pt1, pt2, last_pt
    integer :: int1, int2
    integer, dimension(:), allocatable :: pts

    ! init the base (parent) class
    call init_elem(elem, ielem, grd, neqs)  
    npe = size(elem%psi, 1)

    ! do additional init for dg2d element
    select type (elem)

       class is (element)

          ! do nothing!

       class is (element_dg2d) !further init/alloc for DG element

          elem%number = ielem
          elem%npe = npe
          elem%elname = grd%elname(ielem)
          elem%p = grd%p(ielem)
          elem%eltype = grd%eltype(ielem)
          elem%npedg = 0 ! init p=0,1

          if( allocated(elem%x) ) deallocate(elem%x)
          allocate( elem%x(2, npe) )
          elem%x(1, :) = grd%x( grd%icon(ielem, 1:grd%npe(ielem)) ) 
          elem%x(2, :) = grd%y( grd%icon(ielem, 1:grd%npe(ielem)) ) 

          
          elem%gamma = gamma

          if ( allocated(elem%Mass) ) deallocate(elem%Mass)       
          allocate(elem%Mass(neqs * npe, neqs * npe))
          elem%Mass = 0.0d0

          if ( allocated(elem%U) ) deallocate(elem%U)
          allocate(elem%U(neqs, npe))
          elem%U = 0.0d0

          if ( allocated(elem%d_psi_d_x) ) deallocate(elem%d_psi_d_x)
          allocate(elem%d_psi_d_x(npe, elem%ngauss, 2))
          elem%d_psi_d_x = 0.0d0
          call elem%comp_dpsi()

          if ( allocated(elem%Fk) ) deallocate(elem%Fk)
          allocate(elem%Fk(neqs, elem%ngauss, 2))
          elem%Fk = 0.0d0

          ! select the coefficient of the Jacobian of the transformation
          select case (elem%elname)
          case ( GEN_QUADRI)
             elem%coeff = 1.0d0
             elem%nedgs = 4
          case ( GEN_TRIANGLE)
             elem%coeff = 0.5d0
             elem%nedgs = 3
          end select

          !
          ! allocate and init edges
          !
          if ( allocated(elem%edgs) ) deallocate(elem%edgs)
          allocate(elem%edgs(elem%nedgs))

          ! adding neighbors
          select case (elem%elname)
          case ( GEN_QUADRI)

             ! allocate/init neighbors (initially one neighbor!)
             do ii = 1, 4
                allocate(elem%edgs(ii)%neighs(1))
             end do
             elem%edgs(2)%neighs(1)%elnum = grd%e2e(ielem, 1)
             elem%edgs(3)%neighs(1)%elnum = grd%e2e(ielem, 2)
             elem%edgs(4)%neighs(1)%elnum = grd%e2e(ielem, 3)
             elem%edgs(1)%neighs(1)%elnum = grd%e2e(ielem, 4)

          case ( GEN_TRIANGLE)

             ! allocate/init neighbors (initially one neighbor!)
             do ii = 1, 3
                allocate(elem%edgs(ii)%neighs(1))
             end do
             elem%edgs(2)%neighs(1)%elnum = grd%e2e(ielem, 1)
             elem%edgs(3)%neighs(1)%elnum = grd%e2e(ielem, 2)
             elem%edgs(1)%neighs(1)%elnum = grd%e2e(ielem, 3)

          end select

          ! adding points per edge
          if ( elem%p > 0 ) then ! we have points on the edges 
             elem%npedg = size(grd%el2edg(ielem)%edg1) + 2 
             ! end points are included!

             last_pt = elem%nedgs + 1
             allocate(pts(elem%npedg))

             do ii = 1, size(elem%edgs)
                pt1 = ii
                pt2 = ii + 1
                if ( ii .eq. size(elem%edgs) ) pt2 = 1

                if ( elem%p >= 2 ) then !higher order elements
                   int1 = last_pt
                   int2 = last_pt + elem%npedg - 3
                   pts = (/ pt1, (/ (jj, jj = int1, int2) /) , pt2 /)
                   last_pt = int2 + 1
                else
                   pts = (/ pt1, pt2 /)
                end if

                allocate(elem%edgs(ii)%pts(elem%npedg))
                elem%edgs(ii)%pts = pts                
             end do

             deallocate(pts)
          end if

          ! specify the tag of the edges
          elem%edgs(:)%tag = 0 ! initially everything is interior edge
          if ( grd%el2bn(ielem, 1) .ne. 0 ) then ! there is a boundary edge
             edgnum = grd%el2bn(ielem, 2)
             local_edg = grd%ibedgeELEM_local_edg(edgnum)
             elem%edgs(local_edg)%tag =  grd%ibedgeBC(edgnum)
          end if

       class default

          print *, 'unknown type of element object! stop'
          stop

    end select


    ! done here
  end subroutine init_elem_dg2d

  ! computes mass matrix for time marching of dg2d element
  subroutine comp_mass_mat_dg2d(elem)
    implicit none
    class(element_dg2d), intent(inout) :: elem

    ! local vars
    integer :: i, j, k, l, ii, jj 
    real*8, dimension(:, :), allocatable :: MM


    allocate(MM(elem%neqs * elem%npe, elem%neqs * elem%npe)) 

    ! HARD reset
    MM = 0.0d0
    elem%Mass = 0.0d0

    ! fill out mass matrix for scalar eqs.
    do i = 1, elem%npe
       do j = 1, elem%npe
          do k = 1, elem%ngauss

             ! accumulate to the mass matrix of this element
             MM(i,j) = MM(i,j) &
                  + elem%psi(i,k) * elem%psi(j,k) &
                  * elem%coeff * elem%JJ(k) * elem%W(k)

          end do
       end do
    end do

    ! distribute diagonally to the mass matrix of element
    ! containing system of equations
    do i = 1, elem%npe
       do j = 1, elem%npe
          do l = 1, elem%neqs
             ii = (i-1) * elem%neqs + l
             jj = (j-1) * elem%neqs + l
             elem%Mass(ii, jj) = MM(i, j)
          end do
       end do
    end do

    ! clean ups
    deallocate(MM)

    ! done here
  end subroutine comp_mass_mat_dg2d

  ! evaluates U at a local point (r, s)
  ! in this discontinious element
  subroutine comp_u_point_dg(elem, r, s, u)
    implicit none
    class(element_dg2d), intent(inout) :: elem
    real*8, intent(in) :: r, s
    real*8, dimension(:), intent(out) :: u

    ! local vars
    integer :: k
    real*8, dimension(size(elem%psi, 1)) :: psi

    ! hard reset
    u = 0.0d0

    ! eval basis funcs at point (r, s)
    call elem%tbasis%eval(r, s, 0, psi)

    ! find u ...
    do k = 1, elem%npe
       u = u + elem%U(:,k) * psi(k)
    end do

    ! done here
  end subroutine comp_u_point_dg
 
  ! evaluates the flux at gauss points
  ! and store them in elem%Fk
  subroutine comp_flux_interior(elem)
    implicit none
    class(element_dg2d), intent(inout) :: elem

    ! local vars
    integer :: k
    real*8, dimension(elem%neqs) :: Uk


    do k = 1, elem%ngauss

       ! evaluate U at the current Gauss point 
       call elem%comp_u(elem%r(k), elem%s(k), Uk)

       ! evaluate Fluxes and store them
       call calc_pure_euler2d_flux(Q = Uk, gamma = elem%gamma &
            , F = elem%Fk(:,k, 1), G = elem%Fk(:,k, 2))

    end do

    ! done here
  end subroutine comp_flux_interior

  ! computes Jacobian of the transformation "jac" and 
  ! "Jstar" or the inverse of the Jacobian
  ! of the transformation at any local point "(r,s)" 
  ! within the element "elem". 
  ! The determinant of Jacobian of the transformation 
  ! is returned in "JJ".  
  subroutine comp_Jstar_point_dg(elem, r, s, jac, Jstar, JJ)
    implicit none
    class(element_dg2d), intent(inout) :: elem
    real*8, intent(in) :: r, s
    real*8, dimension(:,:), intent(out) :: jac, Jstar
    real*8, intent(out) :: JJ

    ! local vars
    integer :: i
    real*8 :: xi, yi
    real*8, dimension(2), save :: der
    real*8 :: dx_dr, dx_ds, dy_dr, dy_ds
    real*8, dimension(size(elem%psi,1)) :: d_psi_d_xi, d_psi_d_eta

    ! hard reset
    jac = 0.0d0; Jstar = 0.0d0; JJ = 0.0d0
    xi = 0.0d0; yi = 0.0d0; der = 0.0d0
    dx_dr = 0.0d0; dx_ds = 0.0d0; dy_dr= 0.0d0; dy_ds = 0.0d0

    ! computing the derivative of basis function at that point
    call elem%tbasis%eval(r, s, 1, d_psi_d_xi)
    call elem%tbasis%eval(r, s, 2, d_psi_d_eta)

    ! compute the components of jac
    do i = 1, elem%npe ! assuming iso-geometric expansion for x, y

       der(1) =  d_psi_d_xi(i)
       der(2) = d_psi_d_eta(i)
       xi = elem%x(1,i)
       yi = elem%x(2,i)

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
       print *, 'error in comp_Jstar_point_dg: negative or zero' &
            , ' Jacobian at element (',elem%number,')'
       stop
    end if

    ! comp inverse of jac and store it in Jstar
    Jstar(1,1) =  1.0d0 / JJ * jac(2,2)
    Jstar(2,2) =  1.0d0 / JJ * jac(1,1)
    Jstar(1,2) = -1.0d0 / JJ * jac(1,2)
    Jstar(2,1) = -1.0d0 / JJ * jac(2,1)

    ! done here
  end subroutine comp_Jstar_point_dg

  subroutine comp_d_psi_d_x(elem)
    implicit none
    class(element_dg2d), intent(inout) :: elem

    ! local vars
    integer :: i, k
    real*8, dimension(2) :: der

    ! compute and store ...
    do i = 1, elem%npe
       do k = 1, elem%ngauss
          ! get grad of basis functions
          der(1) = elem%d_psi_d_xi(i,k); der(2) = elem%d_psi_d_eta(i,k)
          ! transform computational grads to physical grads
          der = matmul(elem%Jstar(:,:,k), der)
          ! store
          ! (1..elem%npe, 1..ngauss, 1..2)
          elem%d_psi_d_x(i, k, :) = der
       end do
    end do

    ! done here
  end subroutine comp_d_psi_d_x

  ! computes int_Omega(w,j Fij dOmega) in the interior 
  ! of the DG element.
  subroutine comp_inter_integral(elem, integ)
    implicit none
    class(element_dg2d), intent(in) :: elem
    real*8, dimension(:,:), intent(out) :: integ 

    ! local vars
    integer :: i, k, idim

    ! HARD reset
    integ = 0.0d0

    do i = 1, elem%npe
       do k = 1, elem%ngauss
          do idim = 1, 2 !2d case
             integ(:, i) = integ(:, i) &
                  + elem%d_psi_d_x(i, k, idim) * elem%Fk(:,k, idim) &
                  * elem%coeff * elem%JJ(k) * elem%W(k)
          end do
       end do
    end do

    ! done here
  end subroutine comp_inter_integral

end module element_opt_dg2d