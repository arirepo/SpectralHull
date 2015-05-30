module element_opt_dg2d
  use element_opt
  use grid_opt
  use euler2d_eqs
  implicit none

  private

  type phys_props

     ! viscous props
     logical :: is_viscous = .false.
     real*8 :: mu, lambda, Pr
     logical :: adia

  end type phys_props

  type neigh_dg
     ! the neighbor element number and number of 
     !gauss points per this common segment
     integer :: elnum, ngpseg
     ! physical coords (x,y) of start and end of this segment shared with neighbor
     real*8 :: xs(2), xe(2) 
     ! the element-wise location of 1d Gauss legendre quad rule on the segment
     ! dim = 1:ngpseg
     real*8, dimension(:), allocatable :: xi, W
     ! local coords in the element (xloc_in) 
     ! and in the neighboring element (xloc_out)
     ! (1, 1:ngpseg) = xi, (2, 1:ngpseg) = eta 
     real*8, dimension(:,:), allocatable :: xloc_in, xloc_out
     ! physical coords (x(t), y(t)) in parametric space <t> 
     ! and derivatives (dx/dt, dy/dt) at the gauss points 
     ! (1, 1:ngpseg) = x, (2, 1:ngpseg) = y 
     real*8, dimension(:,:), allocatable :: x, dx
     ! unit normal vector at the gauss points of this segment
     ! (1, 1:ngpseg) = nx , (2, 1:ngpseg) = ny
     real*8, dimension(:, :), allocatable :: n
     ! arc length at the gauss points of this segment
     ! s(k) = sqrt(xdot(k)^2 + ydot(k)^2) , k = 1..ngpseg 
     real*8, dimension(:), allocatable :: s 
     ! all basis functions evaluated at all gauss points per
     ! this shared segment with neighbor element.
     ! (1:npe, 1:ngpseg)
     real*8, dimension(:, :), allocatable :: psi_in

     ! Riemann flux on the edge stored at Gauss points (1:neqs, 1:ngpseg) 
     real*8, dimension(:, :), allocatable :: Fstar

     ! derivative of split flux dF(+,-)/du on the edge 
     ! stored at Gauss points (1:neqs, 1:neqs, 1=plus;2=minus, 1:ngpseg) 
     real*8, dimension(:, :, :, :), allocatable :: dFpm

     ! vars related to viscous terms
     ! u-average stored at Gauss points (1:neqs, 1:ngpseg) 
     real*8, dimension(:, :), allocatable :: Hs

  end type neigh_dg

  type edg_dg
     integer :: tag !can be used for BCs
     integer, dimension(:) , allocatable :: pts
     type(neigh_dg), dimension(:) , allocatable :: neighs
  end type edg_dg

  type, extends(element) :: element_dg2d
     private
     integer, public :: number, npe, elname, p, eltype, npedg, nedgs
     ! (1..2, 1..npe) (1, npe) = x <> (2, npe) = y
     real*8, dimension(:, :), allocatable, public :: x !physic. coords
     ! (1..2, 1..npe) (1, npe) = xi <> (2, npe) = eta
     real*8, dimension(:, :), allocatable :: x_loc !local coords
     real*8, public :: gamma
     ! (neqs * npe, neqs * npe)
     real*8, dimension(:, :), allocatable, public :: Mass
     real*8, dimension(:, :), allocatable, public :: LUmass, LUmass_imp, LU_imposed_mass 
     ! the pivot matrix (neqs * npe)
     integer, dimension(:), allocatable, public :: IPIVmass, IPIVmass_imp, IPIV_LU_imposed_mass
     ! elemental solution : Uij (i=1,neqs <> j=1,npe) 
     real*8, dimension(:, :), allocatable, public :: U, Us, rhs, So
     ! physical derivatives of basis functions at Gauss points
     ! (1..npe, 1..ngauss, 1..2)
     real*8, dimension(:, :, :), allocatable, public :: d_psi_d_x
     ! fluxes evaluated at Gauss points (1..neqs, 1..ngauss, 1..2) 
     real*8, dimension(:, :, :), allocatable, public :: Fk
     real*8, public :: coeff
     ! pure flux jac. evaluated at Gauss points (1..neqs, 1..neqs, 1..ngauss, 1..2) 
     real*8, dimension(:, :, :, :), allocatable, public :: dFk

     ! data struct related to matrix-free iterative procedures 
     real*8, dimension(:, :), allocatable, public :: U0, Ax

     ! data struct related to viscous terms and Navier-Stokes eqs
     ! (neqs*npe, ndim=1..2)
     real*8, dimension(:, :), allocatable, public :: Wij
     type(phys_props), public :: local_props
     integer, dimension(:), allocatable, public :: wall_nodes

     !
     type(edg_dg), dimension(:), allocatable, public :: edgs

  contains

     procedure, nopass, public :: init => init_elem_dg2d
     procedure, public :: comp_mass => comp_mass_mat_dg2d
     procedure, public :: comp_u => comp_u_point_dg
     procedure, public :: comp_flux_interior
     procedure, public :: comp_metric_jacobian => comp_Jstar_point_dg
     procedure, private :: comp_dpsi => comp_d_psi_d_x
     procedure, public :: comp_int_integ => comp_inter_integral
     procedure, public :: xy2rs => xy2rs_dg
     procedure, public :: comp_bnd_integ => comp_bnd_integral
     procedure, public :: init_elem_U
     procedure, public :: update_explicit_euler
     procedure, public :: comp_x => comp_x_point_dg
     procedure, public :: comp_source => comp_elem_source_term
     procedure, public :: init_mms => init_elem_mms
     procedure, public :: comp_pure_flux_jac => comp_pure_flux_jac_interior
     procedure, public :: comp_int_jac_integ => comp_inter_jac_integral
     procedure, public :: comp_inter_ui_integral
     procedure, public :: comp_bnd_ubar_integral
     procedure, public :: comp_Wij_point_dg
     procedure, public :: apply_non_slip_wall_bc
     procedure, public :: create_imposed_mass_matrix
     procedure, public :: comp_elem_source_term_viscous

  end type element_dg2d



  public :: element_dg2d, edg_dg, neigh_dg, phys_props


contains
  
  subroutine init_elem_dg2d(elem, ielem, grd, neqs, gamma)
    implicit none
    class(element), intent(inout) :: elem
    integer, intent(in) :: ielem ! element number
    type(grid), intent(inout) :: grd
    integer, intent(in) :: neqs
    real*8, intent(in) :: gamma

    ! local vars
    integer :: npe, ii, jj, pt1, pt2, last_pt
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
          elem%npedg = elem%p + 1 ! init p=0,1

          if( allocated(elem%x) ) deallocate(elem%x)
          allocate( elem%x(2, npe) )
          elem%x(1, :) = grd%x( grd%icon(ielem, 1:grd%npe(ielem)) ) 
          elem%x(2, :) = grd%y( grd%icon(ielem, 1:grd%npe(ielem)) ) 

          if( allocated(elem%x_loc) ) deallocate(elem%x_loc)
          allocate( elem%x_loc(2, npe) )
          elem%x_loc(1, :) = grd%maselem(ielem)%xi
          elem%x_loc(2, :) = grd%maselem(ielem)%eta

          
          elem%gamma = gamma

          if ( allocated(elem%Mass) ) deallocate(elem%Mass)       
          allocate(elem%Mass(neqs * npe, neqs * npe))
          elem%Mass = 0.0d0

          if ( allocated(elem%LUmass) ) deallocate(elem%LUmass)       
          allocate(elem%LUmass(neqs * npe, neqs * npe))
          elem%LUmass = 0.0d0

          if ( allocated(elem%LUmass_imp) ) deallocate(elem%LUmass_imp)       
          allocate(elem%LUmass_imp(neqs * npe, neqs * npe))
          elem%LUmass_imp = 0.0d0

          if ( allocated(elem%LU_imposed_mass) ) deallocate(elem%LU_imposed_mass)       
          allocate(elem%LU_imposed_mass(neqs * npe, neqs * npe))
          elem%LU_imposed_mass = 0.0d0

          if ( allocated(elem%IPIVmass) ) deallocate(elem%IPIVmass)       
          allocate(elem%IPIVmass(neqs * npe))
          elem%IPIVmass = (/ (ii, ii = 1, (neqs * npe) ) /)

          if ( allocated(elem%IPIVmass_imp) ) deallocate(elem%IPIVmass_imp)       
          allocate(elem%IPIVmass_imp(neqs * npe))
          elem%IPIVmass_imp = (/ (ii, ii = 1, (neqs * npe) ) /)

          if ( allocated(elem%IPIV_LU_imposed_mass) ) deallocate(elem%IPIV_LU_imposed_mass)
          allocate(elem%IPIV_LU_imposed_mass(neqs * npe))
          elem%IPIV_LU_imposed_mass = (/ (ii, ii = 1, (neqs * npe) ) /)

          if ( allocated(elem%U) ) deallocate(elem%U)
          allocate(elem%U(neqs, npe))
          elem%U = 0.0d0

          if ( allocated(elem%Us) ) deallocate(elem%Us)
          allocate(elem%Us(neqs, npe))
          elem%Us = 0.0d0

          if ( allocated(elem%rhs) ) deallocate(elem%rhs)
          allocate(elem%rhs(neqs, npe))
          elem%rhs = 0.0d0

          if ( allocated(elem%So) ) deallocate(elem%So)
          allocate(elem%So(neqs, elem%ngauss))
          elem%So = 0.0d0

          if ( allocated(elem%U0) ) deallocate(elem%U0)
          allocate(elem%U0(neqs, npe))
          elem%U0 = 0.0d0

          if ( allocated(elem%Ax) ) deallocate(elem%Ax)
          allocate(elem%Ax(neqs, npe))
          elem%Ax = 0.0d0

          if ( allocated(elem%Wij) ) deallocate(elem%Wij)
          allocate(elem%Wij((neqs*npe), 2))
          elem%Wij = 0.0d0

          if ( allocated(elem%d_psi_d_x) ) deallocate(elem%d_psi_d_x)
          allocate(elem%d_psi_d_x(npe, elem%ngauss, 2))
          elem%d_psi_d_x = 0.0d0
          call elem%comp_dpsi()

          if ( allocated(elem%Fk) ) deallocate(elem%Fk)
          allocate(elem%Fk(neqs, elem%ngauss, 2))
          elem%Fk = 0.0d0

          if ( allocated(elem%dFk) ) deallocate(elem%dFk)
          allocate(elem%dFk(neqs, neqs, elem%ngauss, 2))
          elem%dFk = 0.0d0

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

          !
          ! allocate/init neighbors (initially one neighbor!)
          !
          ! determine the start and the end of the segments shared with neighbors
          do ii = 1, elem%nedgs
             allocate(elem%edgs(ii)%neighs(1))
             pt1 = ii
             pt2 = ii + 1
             if ( ii .eq. elem%nedgs ) pt2 = 1
             elem%edgs(ii)%neighs(1)%xs = elem%x(:, pt1)
             elem%edgs(ii)%neighs(1)%xe = elem%x(:, pt2) 
          end do

          ! adding neighbors
          select case (elem%elname)
          case ( GEN_QUADRI)

             elem%edgs(2)%neighs(1)%elnum = grd%e2e(ielem, 1)
             elem%edgs(3)%neighs(1)%elnum = grd%e2e(ielem, 2)
             elem%edgs(4)%neighs(1)%elnum = grd%e2e(ielem, 3)
             elem%edgs(1)%neighs(1)%elnum = grd%e2e(ielem, 4)

          case ( GEN_TRIANGLE)

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
          do ii = 1, size(elem%edgs)
             elem%edgs(ii)%tag = grd%local_edg_bc(elem%number, ii)
          end do

          ! ! specify the tag of the edges
          ! elem%edgs(:)%tag = 0 ! initially everything is interior edge
          ! if ( grd%el2bn(ielem, 1) .ne. 0 ) then ! there is a boundary edge
          !    edgnum = grd%el2bn(ielem, 2)
          !    local_edg = grd%ibedgeELEM_local_edg(edgnum)
          !    elem%edgs(local_edg)%tag =  grd%ibedgeBC(edgnum)
          ! end if

          ! compute and store mass matrix and its LU info
          call elem%comp_mass()

          ! compute source term for MMS
          ! call elem%comp_source()
          call elem%comp_elem_source_term_viscous()

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

    ! LAPACK LU temp. vars
    integer :: LDA, INFO

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

    ! compute and store LU of the mass matrix
    elem%LUmass = elem%Mass
    LDA = size(elem%Mass, 1)
    call dgetrf(LDA, LDA, elem%LUmass, LDA, elem%IPIVmass, INFO)
    if ( INFO .ne. 0) then
       print *, 'something is wrong in LU factorization of mass matrix! stop'
       stop
    end if

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
    integer :: k, idim
    real*8, dimension(elem%neqs) :: Uk
    real*8, dimension(elem%neqs, 2) :: Wk, Fv

    do k = 1, elem%ngauss

       ! evaluate U at the current Gauss point 
       call elem%comp_u(elem%r(k), elem%s(k), Uk)

       ! evaluate Fluxes and store them
       call calc_pure_euler2d_flux(Q = Uk, gamma = elem%gamma &
            , F = elem%Fk(:,k, 1), G = elem%Fk(:,k, 2))

       if ( elem%local_props%is_viscous ) then ! add viscous fluxes

          ! evaluate the gradient of conservative vars
          call elem%comp_Wij_point_dg(elem%r(k), elem%s(k), Wk)

          ! compute viscous fluxes
          call calc_Fv(elem%local_props%mu, elem%local_props%lambda &
               , elem%gamma, elem%local_props%Pr &
               , Wk, Uk, Fv, .false.)

          ! finally add viscous fluxes
          do idim = 1, 2
             elem%Fk(:,k, idim) = elem%Fk(:,k, idim) - Fv(:, idim)
          end do

       end if

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
! if (elem%number .eq. 3) print *, 'i = ', i, 'derder = ', der
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
          ! add MMS
          integ(:, i) = integ(:, i) &
               + elem%psi(i, k) * elem%So(:,k) &
               * elem%coeff * elem%JJ(k) * elem%W(k)
          ! if ( all (elem%edgs(:)%tag .eq. 0) ) then
          !    print *, 'mms ', elem%number, elem%coeff, elem%JJ(k), elem%W(k), elem%ngauss
          ! end if

       end do
    end do

    ! done here
  end subroutine comp_inter_integral

  ! computes local coords (r,s) given the
  ! global coordinates (x,y) using a robust
  ! Newton method.
  ! 
  ! This should be applied consistently to
  ! higher-order elements also.
  
  subroutine xy2rs_dg(elem, x, y, maxitr, tolrs, r, s)
    implicit none
    class(element_dg2d), intent(inout) :: elem
    real*8, intent(in) :: x, y
    integer, intent(in) :: maxitr
    real*8, intent(in) :: tolrs
    real*8, intent(out) :: r, s

    ! local vars
    integer :: i, j, npt, pt 
    real*8 :: val
    real*8, dimension(elem%npe) :: psi, d_psi_d_xi, d_psi_d_eta
    real*8, dimension(2), save :: der, F, delta
    real*8, dimension(2,2), save:: JJ, Jinv

    ! initial guess
    r = 0.25d0; s =0.25d0
    npt = elem%npe

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

          pt = i
          val = psi(i)
          der(1) = d_psi_d_xi(i); der(2) = d_psi_d_eta(i)

          F = F - (/ (elem%x(1, pt) * val) , (elem%x(2, pt) * val) /)
          JJ(1,1) =  JJ(1,1) + elem%x(1, pt) * der(1) 
          JJ(1,2) =  JJ(1,2) + elem%x(1, pt) * der(2)
          JJ(2,1) =  JJ(2,1) + elem%x(2, pt) * der(1)
          JJ(2,2) =  JJ(2,2) + elem%x(2, pt) * der(2)

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
          if ( (elem%p < 2) .and. (j > 2) ) then
             print *, 'number of itr is ', j, ' and is greater than 2 ' &
                    , 'for linear element with number', elem%number, '! stop.'
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

    print *, 'xy2rs_dg(...) did not converge' &
         , ' to desired tolerance in ', maxitr, 'iterations' &
         , ' at point (', x, y, ')! stop.'
    stop

    ! done here
  end subroutine xy2rs_dg

  ! computes int_dOmega(w Fijnj dn) on the current <tneigh>
  ! shared segment and accumulates the result to <integ>!
  !
  ! HINT : be careful! freez the initial value of <integ>
  ! before looping over edges of the element and then looping
  ! on the neighbors per edge!
  subroutine comp_bnd_integral(elem, tneigh, integ)
    implicit none
    class(element_dg2d), intent(in) :: elem
    type(neigh_dg), intent(in) :: tneigh
    real*8, dimension(:,:), intent(inout) :: integ 

    ! local vars
    integer :: i, k

    ! HINT : index <i> can be restricted to nodes 
    ! on the current edg to improve performance!
    do k = 1, tneigh%ngpseg

       ! add to boundary integral
       do i = 1, elem%npe
          integ(:, i) = integ(:, i) &
               + tneigh%psi_in(i, k) * tneigh%Fstar(:, k) &
               * tneigh%s(k) * tneigh%W(k)
       end do

! if (elem%number .eq. 3) then
!    print *, 'tneigh%s(k) = ', tneigh%s(k), 'tneigh%W(k)=', tneigh%W(k) &
!         , 'elem', elem%number, 'xs=', tneigh%xs, 'xe=', tneigh%xe &
!         , 'tneigh%ngpseg = ', tneigh%ngpseg, 'tneigh%n(:,k)=' &
!         , tneigh%n(:,k), 'tneigh%elnum = ', tneigh%elnum &
!         , 'xloc_in = ', tneigh%xloc_in, 'xloc_out = ', tneigh%xloc_out 
! end if

    end do ! next gauss point per neighbor segment

    ! done here
  end subroutine comp_bnd_integral

  ! initializes the element's conservative vars elem%U
  ! with the given constant primitive vars
  subroutine init_elem_U(elem, rho, u, v, P)
    implicit none
    class(element_dg2d), intent(inout) :: elem
    real*8, intent(in) :: rho, u, v, P

    ! local vars
    integer :: i
! real*8 :: xx, yy, rr, uu, vv, CC

    do i = 1, elem%npe

! xx = elem%x(1, i); yy = elem%x(2, i)
! rr = sqrt(xx**2.0d0 + yy**2.0d0)
! CC = -2.0d0

! uu = (1.0d0 - exp( CC * (rr - 0.5d0)**2.0d0 )) * u
! vv = (1.0d0 - exp( CC * (rr - 0.5d0)**2.0d0 )) * v

       ! u2U(rho, u, v, P, gamma, UU)
       call u2U(rho, u, v, P, elem%gamma, elem%U(:, i))

    end do

    ! done here
  end subroutine init_elem_U

  ! updates the conservative solution vector
  ! i.e. elem%U by LU solve of mass matrix on the
  ! given residual  
  subroutine update_explicit_euler(elem, rhs, dt)
    implicit none
    class(element_dg2d), intent(inout) :: elem
    real*8, dimension(:, :), intent(inout) :: rhs
    real*8, intent(in) :: dt

    ! local vars
    integer :: N, INFO
    real*8, dimension(:, :), allocatable :: rhs_lapack, rhs_new

    ! init
    N = elem%neqs * elem%npe
    allocate(rhs_lapack(N, 1))
    allocate(rhs_new(elem%neqs, elem%npe))

    ! apply zero-rhs if imposed wall exist in viscous flow
    if ( elem%local_props%is_viscous ) then
       call elem%apply_non_slip_wall_bc(rhs)
    end if

    ! solve using already stored LU
    rhs_lapack = reshape(rhs, (/ N, 1 /) )

    if ( elem%local_props%is_viscous ) then
       if ( allocated(elem%wall_nodes) ) then

          CALL DGETRS( 'No transpose', N, 1, elem%LU_imposed_mass &
               , N, elem%IPIV_LU_imposed_mass, rhs_lapack, N, INFO )

       else

          CALL DGETRS( 'No transpose', N, 1, elem%LUmass &
               , N, elem%IPIVmass, rhs_lapack, N, INFO )

       end if
    else

       CALL DGETRS( 'No transpose', N, 1, elem%LUmass &
            , N, elem%IPIVmass, rhs_lapack, N, INFO )

    end if

    if ( INFO .ne. 0) then
       print *, 'something is wrong in LU solve in explicit euler update! stop'
       stop
    end if

    rhs_new = reshape( rhs_lapack, (/ elem%neqs, elem%npe /))

    ! update
    elem%U = elem%U + dt * rhs_new   


    ! clean ups
    deallocate(rhs_lapack, rhs_new)

    ! done here
  end subroutine update_explicit_euler

  ! evaluates coords x at a local point (r, s)
  ! in this discontinious element
  subroutine comp_x_point_dg(elem, r, s, x)
    implicit none
    class(element_dg2d), intent(inout) :: elem
    real*8, intent(in) :: r, s
    real*8, dimension(:), intent(out) :: x

    ! local vars
    integer :: k
    real*8, dimension(elem%npe) :: psi

    ! hard reset
    x = 0.0d0

    ! eval basis funcs at point (r, s)
    call elem%tbasis%eval(r, s, 0, psi)

    ! find x ...
    do k = 1, elem%npe
       x = x + elem%x(:,k) * psi(k)
    end do

    ! done here
  end subroutine comp_x_point_dg

  ! computes element's sorce term for MMS
  !
  subroutine comp_elem_source_term(elem)
    implicit none
    class(element_dg2d), intent(inout) :: elem

    ! local vars
    integer :: k
    real*8, dimension(2) :: xx
    real*8 :: x, y

    do k = 1, elem%ngauss

       call elem%comp_x(elem%r(k), elem%s(k), xx)
       x = xx(1); y = xx(2)

if (elem%number .eq. 4) print *, 'xso = ', x, 'yso = ', y 
       elem%So(1, k) =  x ** 2 + 2.0d0 * x * y + 1.0d0

       elem%So(2, k) = x ** 3 + 2.0d0 * x ** 2 * y + 2.0d0 * x * y ** 2 + x + 2.0d0 * y

       elem%So(3, k) = 2.0d0 * x ** 3 + 5.0d0 * x ** 2 * y + 2.0d0 * x * y ** 2 + 2.0d0 * x + 5.0d0 * y

       elem%So(4, k) = dble((3.0d0 * elem%gamma * x ** 4 + 12.0d0 * elem%gamma * x ** 3 * y + 12.0d0 *& 
     elem%gamma * x ** 2 * y ** 2 + 4.0d0 * elem%gamma * x * y ** 3 - 3.0d0 * x ** 4 - &
     12.0d0 * x ** 3 * y - 12.0d0 * x ** 2 * y ** 2 - 4.0d0 * x * y ** 3 + 3.0d0 * elem%gamma * &
     x ** 2 + 14.0d0 * elem%gamma * x * y + 14.0d0 * y ** 2 * elem%gamma - 3.0d0 * x ** 2 &
     - 10.0d0 * x * y - 8.0d0 * y ** 2 + 2.0d0 * elem%gamma) / (elem%gamma - 1.0d0)) / 0.2D1

    end do

    ! done here
  end subroutine comp_elem_source_term

  !
  ! method of manufacturing
  !
  subroutine init_elem_mms(elem)
    implicit none
    class(element_dg2d), intent(inout) :: elem

    ! local vars
    integer :: i
    real*8 :: rho, u, v, P, x, y

    do i = 1, elem%npe

       x = elem%x(1, i)
       y = elem%x(2, i)

       rho = 1.0d0 + x**2  
       u = y
       v = x+y
       P = 1.0d0 + y**2

       ! u2U(rho, u, v, P, gamma, UU)
       call u2U(rho, u, v, P, elem%gamma, elem%U(:, i))

    end do

    ! done here
  end subroutine init_elem_mms

  ! evaluates the pure flux jac. at gauss points
  ! and store them in elem%dFk
  subroutine comp_pure_flux_jac_interior(elem)
    implicit none
    class(element_dg2d), intent(inout) :: elem

    ! local vars
    integer :: k
    real*8, dimension(elem%neqs) :: Uk

    do k = 1, elem%ngauss

       ! evaluate U at the current Gauss point 
       call elem%comp_u(elem%r(k), elem%s(k), Uk)

       ! evaluate pure flux jac. and store them
       call calc_pure_euler2d_flux(Q = Uk, gamma = elem%gamma &
            , F = elem%Fk(:,k, 1), G = elem%Fk(:,k, 2))

       call calc_pure_euler2d_flux_jac(Q = Uk, gamma = elem%gamma &
            , dF = elem%dFk(:, :, k, 1), dG = elem%dFk(:, :, k, 2))

    end do

    ! done here
  end subroutine comp_pure_flux_jac_interior

  ! computes int_Omega(w,j dFij/du * Xu dOmega) in the interior 
  ! of the DG element.
  subroutine comp_inter_jac_integral(elem, integ)
    implicit none
    class(element_dg2d) :: elem
    real*8, dimension(:,:), intent(out) :: integ 

    ! local vars
    integer :: i, k, idim
    real*8, dimension(elem%neqs) :: Xu

    ! HARD reset
    integ = 0.0d0

    
    do k = 1, elem%ngauss

       ! evaluate U at the current Gauss point
       ! at the current Arnoldi iteration 
       call elem%comp_u(elem%r(k), elem%s(k), Xu)

       do i = 1, elem%npe
          do idim = 1, 2 !2d case
             integ(:, i) = integ(:, i) &
                  + elem%d_psi_d_x(i, k, idim) * matmul(elem%dFk(:,:,k, idim), Xu) &
                  * elem%coeff * elem%JJ(k) * elem%W(k)
          end do
       end do
    end do

    ! done here
  end subroutine comp_inter_jac_integral

  ! computes -int_Omega(w,j * ui dOmega) in the interior 
  ! of the DG element.
  !
  ! NOTE : integ(neqs*npe, ndim)
  !
  subroutine comp_inter_ui_integral(elem, integ)
    implicit none
    class(element_dg2d) :: elem
    real*8, dimension(:, :), intent(out) :: integ 

    ! local vars
    integer :: i, k, idim, i1, i2
    real*8, dimension(elem%neqs) :: u

    ! HARD reset
    integ = 0.0d0

    
    do k = 1, elem%ngauss

       ! evaluate U at the current Gauss point
       call elem%comp_u(elem%r(k), elem%s(k), u)

       do i = 1, elem%npe

          ! select the range
          i1 = (i-1)*elem%neqs + 1
          i2 = i * elem%neqs

          do idim = 1, 2 !2d case

             ! add contribution
             integ(i1:i2, idim) = integ(i1:i2, idim) &
                  - elem%d_psi_d_x(i, k, idim) * u &
                  * elem%coeff * elem%JJ(k) * elem%W(k)
          end do

       end do

    end do

    ! done here
  end subroutine comp_inter_ui_integral

  ! computes int_dOmega(w ubar_i nj) on the current <tneigh>
  ! shared segment and accumulates the result to <integ>!
  !
  ! integ(neqs*npe, ndim)
  !
  subroutine comp_bnd_ubar_integral(elem, tneigh, integ)
    implicit none
    class(element_dg2d), intent(in) :: elem
    type(neigh_dg), intent(in) :: tneigh
    real*8, dimension(:,:), intent(inout) :: integ 

    ! local vars
    integer :: i, k, i1, i2, idim

    ! HINT : index <i> can be restricted to nodes 
    ! on the current edg to improve performance!
    do k = 1, tneigh%ngpseg

       ! add to boundary integral
       do i = 1, elem%npe

          ! select range
          i1 = (i-1)*elem%neqs + 1
          i2 = i * elem%neqs

          do idim = 1, 2
             ! add contribution
             integ(i1:i2, idim) = integ(i1:i2, idim) &
                  + tneigh%psi_in(i, k) * tneigh%Hs(:, k) * tneigh%n(idim, k) &
                  * tneigh%s(k) * tneigh%W(k)
          end do ! next dimension

       end do ! next point per element


    end do ! next gauss point per neighbor segment

    ! done here
  end subroutine comp_bnd_ubar_integral

  ! evaluates Wij at a local point (r, s)
  ! in this discontinious element
  !
  ! W(neqs, ndim)
  !
  subroutine comp_Wij_point_dg(elem, r, s, W)
    implicit none
    class(element_dg2d), intent(inout) :: elem
    real*8, intent(in) :: r, s
    real*8, dimension(:, :), intent(out) :: W

    ! local vars
    integer :: k, i1, i2, idim
    real*8, dimension(size(elem%psi, 1)) :: psi

    ! hard reset
    W = 0.0d0

    ! eval basis funcs at point (r, s)
    call elem%tbasis%eval(r, s, 0, psi)

    ! find W ...
    do k = 1, elem%npe

       ! find the range
       i1 = (k-1) * elem%neqs + 1
       i2 = k * elem%neqs

       ! add to W at point
       do idim = 1, 2
          W(:, idim) = W(:, idim) + elem%Wij(i1:i2, idim) * psi(k)
       end do

    end do

    ! done here
  end subroutine comp_Wij_point_dg

  !
  ! applys non-slip wall boundary conditions
  ! if viscous terms are included
  !
  subroutine apply_non_slip_wall_bc(elem, rhs)
    implicit none
    class(element_dg2d), intent(in) :: elem
    real*8, dimension(:, :), intent(inout) :: rhs

    ! local vars
    integer :: inode

    if ( allocated(elem%wall_nodes) ) then

       do inode = 1, size(elem%wall_nodes)
          rhs(2:3, elem%wall_nodes(inode)) = 0.0d0
       end do

    end if

    ! done here
  end subroutine apply_non_slip_wall_bc

  ! compute the imposed mass matrix and store
  ! its LU factorization for future solves
  ! This is important if viscous wall BCs exist
  ! 
  subroutine create_imposed_mass_matrix(elem)
    implicit none
    class(element_dg2d) :: elem

    ! local vars
    integer :: inode, ieq, ii, ind

    ! LAPACK VARS
    integer :: LDA, INFO 

    ! apply wall bcs (if any)
    if( allocated(elem%wall_nodes) ) then

       ! first take a fresh copy of already computed mass matrix
       elem%LU_imposed_mass = elem%Mass

       ! then impose fixed BCs
       do inode = 1, size(elem%wall_nodes)
          do ieq = 2, 3
             ! compute the range
             ind = elem%wall_nodes(inode)
             ii = (ind - 1) * elem%neqs + ieq
             ! freez the row
             elem%LU_imposed_mass(ii, :) = 0.0d0
             ! put one on the main diagonal
             elem%LU_imposed_mass(ii, ii) = 1.0d0
          end do
       end do

       ! then store LU
       LDA = size(elem%LU_imposed_mass, 1)
       call dgetrf(LDA, LDA, elem%LU_imposed_mass, LDA, elem%IPIV_LU_imposed_mass, INFO)
       if ( INFO .ne. 0) then
          print *, 'something is wrong in LU factorization of imposed mass matrix! stop'
          stop
       end if

    end if

    ! done here
  end subroutine create_imposed_mass_matrix

  ! computes element's sorce term for MMS
  ! of viscous flow
  subroutine comp_elem_source_term_viscous(elem)
    implicit none
    class(element_dg2d), intent(inout) :: elem

    ! local vars
    integer :: k
    real*8, dimension(2) :: xx
    real*8 :: x, y, cg, Pr0, gamma0, mu0

    cg = elem%local_props%lambda
    gamma0 = elem%gamma
    Pr0 = elem%local_props%Pr
    mu0 = elem%local_props%mu

    do k = 1, elem%ngauss

       call elem%comp_x(elem%r(k), elem%s(k), xx)
       x = xx(1); y = xx(2)

       ! fill rhs vector
       elem%So(1, k) =  0.2D1 * x * y + 0.1D1 + x ** 2

       elem%So(2, k) = 0.2D1 * x * y ** 2 + x + 0.2D1 * y + x ** 3 + 0.2D1 * x ** 2 * y

       elem%So(3, k) = 0.5D1 * x ** 2 * y + 0.2D1 * x * y ** 2 + 0.5D1 * y + 0.2D1 * x &
            + 0.2D1 * x ** 3

       elem%So(4, k) = -0.5000000000D0 * (-0.12D2 * y ** 3 * Pr0 * x ** 5 * gamma0 &
            - 0.4D1 * y ** 3 * Pr0 * x ** 7 * gamma0 - 0.50D2 * y ** 2 * Pr0 * &
            x ** 6 * gamma0 - 0.12D2 * y ** 2 * Pr0 * x ** 8 * gamma0 - 0.4D1 &
            * y ** 3 * Pr0 * x * gamma0 - 0.54D2 * y ** 2 * Pr0 * x ** 2 * gamma0 &
            - 0.12D2 * y ** 3 * Pr0 * x ** 3 * gamma0 - 0.78D2 * y ** 2 * &
            Pr0 * x ** 4 * gamma0 - 0.54D2 * y * Pr0 * x ** 3 * gamma0 - 0.14D2 &
            * y * Pr0 * x * gamma0 - 0.78D2 * y * Pr0 * x ** 5 * gamma0 - 0.50D2 &
            * y * Pr0 * x ** 7 * gamma0 - 0.12D2 * y * Pr0 * x ** 9 * gamma0 &
            + 0.36D2 * mu0 * Pr0 * x ** 2 * gamma0 + 0.36D2 * mu0 * Pr0 * &
            x ** 4 * gamma0 + 0.12D2 * mu0 * Pr0 * x ** 6 * gamma0 + 0.12D2 * &
            mu0 * x ** 2 * y ** 2 * gamma0 + 0.2D1 * mu0 * cg * Pr0 * gamma0 & 
            - 0.6D1 * mu0 * cg * Pr0 * x ** 2 - 0.6D1 * mu0 * cg * Pr0 * x ** 4 &
            - 0.2D1 * mu0 * cg * Pr0 * x ** 6 + 0.6D1 * mu0 * cg * Pr0 * x ** 2 &
            * gamma0 + 0.6D1 * mu0 * cg * Pr0 * gamma0 * x ** 4 + 0.2D1 * & 
            mu0 * cg * Pr0 * x ** 6 * gamma0 - 0.9D1 * Pr0 * x ** 2 * gamma0 - &
            0.18D2 * Pr0 * x ** 4 * gamma0 - 0.20D2 * Pr0 * x ** 6 * gamma0 + &
            0.36D2 * y ** 2 * Pr0 * x ** 2 + 0.42D2 * y * Pr0 * x ** 3 + 0.4D1 &
            * y ** 3 * Pr0 * x + 0.10D2 * y * Pr0 * x + 0.66D2 * y * Pr0 * x** 5 &
            + 0.46D2 * y * Pr0 * x ** 7 + 0.12D2 * y * Pr0 * x ** 9 - 0.14D2 &
            * y ** 2 * Pr0 * gamma0 + 0.12D2 * y ** 3 * Pr0 * x ** 3 + 0.60D2 * &
            y ** 2 * Pr0 * x ** 4 + 0.12D2 * y ** 3 * Pr0 * x ** 5 + 0.4D1 &
            * y ** 3 * Pr0 * x ** 7 + 0.44D2 * y ** 2 * Pr0 * x ** 6 + 0.12D2 &
            * y ** 2 * Pr0 * x ** 8 + 0.12D2 * mu0 * Pr0 * gamma0 - 0.36D2 &
            * mu0 * Pr0 * x ** 2 - 0.36D2 * mu0 * Pr0 * x ** 4 - 0.12D2 * mu0 &
            * Pr0 * x ** 6 + 0.2D1 * mu0 * x ** 2 * gamma0 - 0.4D1 * mu0 * y ** 2 &
            * gamma0 + 0.18D2 * mu0 * x ** 2 * gamma0 ** 2 + 0.18D2 * mu0 &
            * x ** 4 * gamma0 ** 2 - 0.14D2 * mu0 * x ** 4 * gamma0 - 0.6D1 * &
            mu0 * x ** 6 * gamma0 + 0.6D1 * mu0 * x ** 6 * gamma0 ** 2 - 0.12D2 &
            * Pr0 * x ** 8 * gamma0 - 0.3D1 * Pr0 * x ** 10 * gamma0 - 0.2D1 &
            * Pr0 * gamma0 + 0.3D1 * Pr0 * x ** 2 + 0.12D2 * Pr0 * x ** 4 + &
            0.18D2 * Pr0 * x ** 6 + 0.8D1 * y ** 2 * Pr0 - 0.12D2 * mu0 * Pr0 & 
            - 0.6D1 * mu0 * gamma0 + 0.6D1 * mu0 * gamma0 ** 2 + 0.12D2 * Pr0 & 
            * x ** 8 + 0.3D1 * Pr0 * x ** 10 - 0.2D1 * mu0 * cg * Pr0) / Pr0 / &
            (gamma0 - 0.1D1) / (0.1D1 + x ** 2) ** 3

    end do

print *, 'gulgulsource', elem%local_props%lambda, elem%local_props%mu, elem%local_props%Pr
    ! done here
  end subroutine comp_elem_source_term_viscous

end module element_opt_dg2d
