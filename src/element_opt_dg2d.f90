module element_opt_dg2d
  use element_opt
  use grid_opt
  use euler2d_eqs
  use spline
  implicit none

  private

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

   contains

     procedure, private :: dealloc_neigh

  end type neigh_dg

  type edg_dg
     integer :: tag !can be used for BCs
     integer, dimension(:) , allocatable :: pts
     type(neigh_dg), dimension(:) , allocatable :: neighs

   contains

     procedure, private :: dealloc_edg

  end type edg_dg

  type, extends(element) :: element_dg2d
     private
     integer, public :: number, npe, elname, p, eltype, npedg, nedgs
     ! (1..2, 1..npe) (1, npe) = x <> (2, npe) = y
     real*8, dimension(:, :), allocatable, public :: x !physic. coords
     ! (1..2, 1..npe) (1, npe) = xi <> (2, npe) = eta
     real*8, dimension(:, :), allocatable, public :: x_loc !local coords
     real*8, public :: gamma
     ! (neqs * npe, neqs * npe)
     real*8, dimension(:, :), allocatable, public :: Mass
     real*8, dimension(:, :), allocatable, public :: LUmass, LUmass_imp 
     ! the pivot matrix (neqs * npe)
     integer, dimension(:), allocatable, public :: IPIVmass, IPIVmass_imp
     ! elemental solution : Uij (i=1,neqs <> j=1,npe) 
     real*8, dimension(:, :), allocatable, public :: U, Us, rhs, So, Un
     real*8, dimension(:, :, :), allocatable, public :: Urk ! (neqs, npe, RKstages)
     ! physical derivatives of basis functions at Gauss points
     ! (1..npe, 1..ngauss, 1..2)
     real*8, dimension(:, :, :), allocatable, public :: d_psi_d_x
     ! fluxes evaluated at Gauss points (1..neqs, 1..ngauss, 1..2) 
     real*8, dimension(:, :, :), allocatable :: Fk
     real*8, public :: coeff
     ! pure flux jac. evaluated at Gauss points (1..neqs, 1..neqs, 1..ngauss, 1..2) 
     real*8, dimension(:, :, :, :), allocatable, public :: dFk

     ! data struct related to matrix-free iterative procedures 
     real*8, dimension(:, :), allocatable, public :: U0, Ax

     !
     type(edg_dg), dimension(:), allocatable, public :: edgs

  contains

     procedure, nopass, public :: init => init_elem_dg2d
     procedure, public :: comp_mass => comp_mass_mat_dg2d
     procedure, public :: comp_u => comp_u_point_dg
     procedure, public :: comp_flux_interior
     procedure, public :: comp_metric_jacobian => comp_Jstar_point_dg
     procedure, public :: comp_dpsi => comp_d_psi_d_x
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
     procedure, public :: comp_dt_Minv_rhs
     procedure, public :: transform_tri
     procedure, public :: transform_quadri
     procedure, public :: alloc_init_loc_matrices
     procedure, public :: dealloc_elem

  end type element_dg2d



  public :: element_dg2d, edg_dg, neigh_dg


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

          ! allocate and initialize the local matrices
          ! subroutine alloc_init_loc_matrices(elem, npe, neqs)
          call elem%alloc_init_loc_matrices(npe, neqs)

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
          call elem%comp_source()

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
!           ! add MMS
!           integ(:, i) = integ(:, i) &
!                + elem%psi(i, k) * elem%So(:,k) &
!                * elem%coeff * elem%JJ(k) * elem%W(k)
! if ( all (elem%edgs(:)%tag .eq. 0) ) then
! print *, 'mms ', elem%number, elem%coeff, elem%JJ(k), elem%W(k), elem%ngauss
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

! print *, ' size(elem%edgs) = ', size(elem%edgs)

! if ( elem%elname .eq. GEN_SPHULL ) then
! r = x
! s = y
! return
! end if

    ! initial guess
!     r = 0.25d0; s =0.25d0
    r = sum(elem%r)/size(elem%r); s = sum(elem%s)/size(elem%s)
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
             !stop
          else if (  (r > (1.0d0 + tolrs)) .or. (s > (1.0d0 + tolrs)) &
               .or.  (r < (0.0d0 - tolrs)) .or. (s < (0.0d0 - tolrs)) ) then
             ! print *, '(r,s) = ', r, s ,' out of range! stop.'
             ! stop 
          end if
!           print *, 'jqjqq = ', j , 'error = ', sqrt(sum(delta*delta)), 'elem%number = ', elem%number, 'tol = ', tolrs
          return
       end if

    end do

    print *, 'xy2rs_dg(...) did not converge' &
         , ' to desired tolerance in ', maxitr, 'iterations' &
         , ' at point (', x, y, ')! the residual is : ', sqrt(sum(delta*delta)) &
         , ' for element name ', elem%elname, ' size(elem%edgs) = ' &
         , size(elem%edgs), 'elem%number = ', elem%number, ' . stop.'
print *, 'elem%x(1, :) ' , elem%x(1, :)
print *, 'elem%x(2, :) ' , elem%x(2, :)
print *, 'elem%x_loc(1, :) ' , elem%x_loc(1, :)
print *, 'elem%x_loc(2, :) ' , elem%x_loc(2, :)
do i = 1, size(elem%edgs)
print *, 'elem%edgs(i)%neighs(1)%elnum = ', elem%edgs(i)%neighs(1)%elnum
end do
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
    real*8 :: r, rho0

    do i = 1, elem%npe

r = (elem%x(1, i) - 1.0d0)**2 + (elem%x(2, i) - 0.5d0)**2
rho0 = (1 + 0.1d0 * exp(-38.0d0 * r )) * rho
       ! u2U(rho, u, v, P, gamma, UU)
       call u2U(rho0, u, v, P, elem%gamma, elem%U(:, i))

    end do

    ! done here
  end subroutine init_elem_U

  ! updates the conservative solution vector
  ! i.e. elem%U by LU solve of mass matrix on the
  ! given residual  
  subroutine update_explicit_euler(elem, rhs, dt)
    implicit none
    class(element_dg2d), intent(inout) :: elem
    real*8, dimension(:, :), intent(in) :: rhs
    real*8, intent(in) :: dt

    ! local vars
    integer :: N, INFO
    real*8, dimension(:, :), allocatable :: rhs_lapack, rhs_new

    ! init
    N = elem%neqs * elem%npe
    allocate(rhs_lapack(N, 1))
    allocate(rhs_new(elem%neqs, elem%npe))

    ! solve using already stored LU
    rhs_lapack = reshape(rhs, (/ N, 1 /) )

    CALL DGETRS( 'No transpose', N, 1, elem%LUmass &
         , N, elem%IPIVmass, rhs_lapack, N, INFO )

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

  ! computes dt * M^-1 * rhs and stores the results
  ! in rhs itself
  !
  ! used in TVD-RK scheme
  !  
  subroutine comp_dt_Minv_rhs(elem, rhs, dt)
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

    ! solve using already stored LU
    rhs_lapack = reshape(rhs, (/ N, 1 /) )

    CALL DGETRS( 'No transpose', N, 1, elem%LUmass &
         , N, elem%IPIVmass, rhs_lapack, N, INFO )

    if ( INFO .ne. 0) then
       print *, 'something is wrong in LU solve in TVDRK update! stop'
       stop
    end if

    rhs_new = reshape( rhs_lapack, (/ elem%neqs, elem%npe /))

    ! update
    rhs = dt * rhs_new   

    ! clean ups
    deallocate(rhs_lapack, rhs_new)

    ! done here
  end subroutine comp_dt_Minv_rhs

  ! general master element r : 0->1, s : 0->1
  ! to physical coordinate transformation
  ! for linear and curved (nonlinear) DG triangle
  !
  ! NOTE : the grd structure is only used for mapping
  ! curved boundaries (if user decides to do so).
  ! if grd%linear_boundaries = .false.; which is default,
  ! then the boundaries are treated as either NURB or spline
  ! curves otherwise if grd%linear_boundaries = .true. then
  ! the boundaries are simply linear.
  !
  subroutine transform_tri(elem, grd, xi, eta, x, y, tol)
    implicit none
    class(element_dg2d), intent(inout) :: elem
    type(grid), intent(in), target :: grd
    real*8, intent(in) :: xi, eta
    real*8, intent(out) :: x, y
    real*8, intent(in) :: tol

    ! local vars
    integer :: nc, tag, pt1, pt2, pt3
    real*8 :: x1_x, x1_y, x2_x, x2_y, x3_x, x3_y
    real*8 :: Cx, Cy
    real*8 :: tt, t1, t2
    type(curve), pointer :: tcurve
    real*8, dimension(:), pointer :: t, xs, ys
    integer :: i, edgnum
    real*8 :: tangx, tangy, norm_tang, tmpx, tmpy
    real*8 :: xtmp(2,3)

    ! determine tag
    if ( all(elem%edgs(:)%tag .eq. 0) .or. grd%linear_boundaries ) then
       tag = 0 ! interior
    elseif ( count(elem%edgs(:)%tag .ne. 0) > 1 ) then
       print *, 'two edges in one element has tag(bc) of non-interior! stop'
       stop
    else ! we are good to go
       do i = 1, size(elem%edgs)
          tag = elem%edgs(i)%tag
          if ( tag .ne. 0 ) then ! we found it!
             edgnum = i
             exit
          else
             cycle
          end if
       end do
    end if

    ! first determine points 1, 2, 3
    if ( (tag .eq. 0) .or. grd%linear_boundaries ) then ! interior
       pt1 = 1
       pt2 = 2
       pt3 = 3
    else ! boundary

       select case (edgnum)
       case (1)
          pt1 = 1; pt2 = 2; pt3 = 3
       case (2)
          pt1 = 2; pt2 = 3; pt3 = 1
       case (3)
          pt1 = 3; pt2 = 1; pt3 = 2
       case default
          print *, 'could not find the point orientation of' &
               , ' boundary curved triangle transformation! stop'
          stop
       end select

       xtmp(:, 1) = elem%x(:, pt1)
       xtmp(:, 2) = elem%x(:, pt2)
       xtmp(:, 3) = elem%x(:, pt3)
       ! swap
       elem%x(1:2, 1:3) = xtmp

    end if

    ! find coordinates
    x1_x =  elem%x(1, 1)
    x1_y =  elem%x(2, 1)
    x2_x =  elem%x(1, 2)
    x2_y =  elem%x(2, 2)
    x3_x =  elem%x(1, 3)
    x3_y =  elem%x(2, 3)

    ! find tangential vector to side y1
    norm_tang = sqrt( (x2_x - x1_x)**2.0d0 + (x2_y - x1_y)**2.0d0 )
    tangx = (x2_x - x1_x) / norm_tang
    tangy = (x2_y - x1_y) / norm_tang

    ! compute Cx, Cy
    if( (tag .eq. 0) .or. grd%linear_boundaries ) then ! no curve bn
       Cx = x1_x + xi * (x2_x - x1_x)
       Cy = x1_y + xi * (x2_y - x1_y)
    else ! do boundary curve interpolation

       tcurve => grd%bn_curves(tag)
       xs => tcurve%x
       ys => tcurve%y
       t  => tcurve%t
       nc = size(xs)

       ! find_t(grd, tag, x, y, tol, t)
       tmpx = x1_x + piecewise_tol * tangx
       tmpy = x1_y + piecewise_tol * tangy 
       call find_t(grd, tag, tmpx, tmpy, tol, t1)
       tmpx = x2_x - piecewise_tol * tangx
       tmpy = x2_y - piecewise_tol * tangy 
       call find_t(grd, tag, tmpx, tmpy, tol, t2)
       t1 = dble(nint(t1))
       t2 = dble(nint(t2))

       tt = t1 + xi * (t2 - t1)

       call spline_nurbs_eval2(tt, xs, ys, tcurve%a, tcurve%b, tcurve%c &
            , tcurve%Mx, tcurve%My, t, Cx, Cy, 'interp', tcurve%btype)

    end if

    ! finalize the transformation
    if ( abs(xi - 1.0d0) <= tol ) then
       x = x2_x 
       y = x2_y
    else
       x = (1.0d0 - xi - eta) / (1.0d0 - xi) * Cx &
            + (xi * eta) / (1.0d0 - xi) * x2_x + eta * x3_x
       y = (1.0d0 - xi - eta) / (1.0d0 - xi) * Cy &
            + (xi * eta) / (1.0d0 - xi) * x2_y + eta * x3_y
    end if


    ! done here
  end subroutine transform_tri

  ! =================================
  ! ============  MAP ===============
  ! =================================
  !                <edge : y3>
  !      (pt4) *-----------------* (pt3)
  !            |                 |
  !            |                 |
  !            |                 |
  !<edge : y4> |                 | <edge : y2>
  !            | --     ------   |
  !            |/  \   /      \  |
  !            *   ----        --*
  !         (pt1)  <edge : y1>  (pt2)
  ! =================================
  ! =================================

  ! the edge y1 corresponding to the only curved side 
  ! of the quad. Matches this curved side of the boundary
  ! quads by replacing this with bn_curves structure.
  ! Otherwise, if the quad element is an interior quad, use straight 
  ! line definition. This is DG implementation ONLY. 
  !
  ! NOTE : the grd structure is only used for mapping
  ! curved boundaries (if user decides to do so).
  ! if grd%linear_boundaries = .false.; which is default,
  ! then the boundaries are treated as either NURB or spline
  ! curves otherwise if grd%linear_boundaries = .true. then
  ! the boundaries are simply linear.
  !

  subroutine transform_quadri(elem, grd, xi, eta, x, y, tol)
    implicit none
    class(element_dg2d), intent(inout) :: elem
    type(grid), target, intent(in) :: grd
    real*8, intent(in) :: xi, eta
    real*8, intent(out) :: x, y
    real*8, intent(in) :: tol

    ! local vars
    integer :: nc, tag
    integer, dimension(4) :: pt
    real*8, dimension(4) :: xx, yy
    real*8 :: tt, t1, t2
    type(curve), pointer :: tcurve => null()
    real*8, dimension(:), pointer :: t => null(), xs => null(), ys => null()
    integer :: i, edgnum, j, k
    logical :: duplicate
    real*8 :: y1x, y1y, y3x, y3y
    real*8 :: area
    integer :: twist, pt_tmp
    real*8 :: tangx, tangy, norm_tang, tmpx, tmpy
    real*8 :: xtmp(2,4)

    ! determine tag
    if ( all(elem%edgs(:)%tag .eq. 0) .or. grd%linear_boundaries ) then
       tag = 0 ! interior or linear boundaries
    elseif ( count(elem%edgs(:)%tag .ne. 0) > 1 ) then
       print *, 'two edges in one quad element has tag(bc) of non-interior! stop'
       stop
    else ! we are good to go
       do i = 1, size(elem%edgs)
          tag = elem%edgs(i)%tag
          if ( tag .ne. 0 ) then ! we found it!
             edgnum = i
             exit
          else
             cycle
          end if
       end do
    end if

    ! first determine points 1, 2, 3 and 4
    if ( (tag .eq. 0) .or. grd%linear_boundaries ) then ! interior
       pt = (/ 1, 2, 3, 4 /)
    else ! boundary

       select case (edgnum)
       case (1)
          pt = (/ 1, 2, 3, 4 /)
       case (2)
          pt = (/ 2, 3, 4, 1 /)
       case (3)
          pt = (/ 3, 4, 1, 2 /)
       case (4)
          pt = (/ 4, 1, 2, 3 /)
       case default
          print *, 'could not find the point orientation of' &
               , ' boundary curved quadrilateral transformation! stop'
          stop
       end select

       xtmp(:, 1) = elem%x(:, pt(1))
       xtmp(:, 2) = elem%x(:, pt(2))
       xtmp(:, 3) = elem%x(:, pt(3))
       xtmp(:, 4) = elem%x(:, pt(4))
       ! swap
       elem%x(1:2, 1:4) = xtmp

    end if

    ! see if the chosen points lead to
    ! a good or twisted element, if
    ! twisted then swap points 3, 4
    ! to make it good element
    !
    ! subroutine comp_quad4_area(x, y, area, twist)
    !
    call comp_quad4_area(elem%x(1, :), elem%x(2, :), area, twist)
    if ( twist .eq. 1 ) then ! twisted!

       !swapping ...
       pt = (/ 1, 2, 3, 4 /)
       pt_tmp = pt(4)
       pt(4) = pt(3)
       pt(3) = pt_tmp
       xtmp(:, 1) = elem%x(:, pt(1))
       xtmp(:, 2) = elem%x(:, pt(2))
       xtmp(:, 3) = elem%x(:, pt(3))
       xtmp(:, 4) = elem%x(:, pt(4))
       ! swap
       elem%x(1:2, 1:4) = xtmp

       ! double check
       call comp_quad4_area(elem%x(1, :), elem%x(2, :), area, twist)
       if (twist .eq. 1) then ! still twisted wow !!!
          print *, 'the quad element #', elem%number, 'is twisted' &
               , ' and cant be fixed! stop'
          stop
       end if

    end if

    ! find coordinates
print *, 'elem%x(1, :) = ',  elem%x(1, :)
print *, 'elem%x(2, :) = ',  elem%x(2, :)
    xx = elem%x(1, 1:4)
    yy = elem%x(2, 1:4)

    ! find tangential vector to side y1
norm_tang = sqrt( (xx(2) - xx(1))**2.0d0 + (yy(2) - yy(1))**2.0d0 )
tangx = (xx(2) - xx(1)) / norm_tang
tangy = (yy(2) - yy(1)) / norm_tang

    ! compute y1 or the possibly curved side
    if( tag .eq. 0 ) then ! no curve bn
       y1x = xx(1) + 0.5d0 * (xi + 1.0d0) * (xx(2) - xx(1))
       y1y = yy(1) + 0.5d0 * (xi + 1.0d0) * (yy(2) - yy(1))

    else ! do boundary curve interpolation

       tcurve => grd%bn_curves(tag)
       xs => tcurve%x
       ys => tcurve%y
       t  => tcurve%t
       nc = size(xs)

       ! find_t(grd, tag, x, y, tol, t)
tmpx = xx(1) + piecewise_tol * tangx
tmpy = yy(1) + piecewise_tol * tangy 
       call find_t(grd, tag, tmpx, tmpy, tol, t1)
tmpx = xx(2) - piecewise_tol * tangx
tmpy = yy(2) - piecewise_tol * tangy 
       call find_t(grd, tag, tmpx, tmpy, tol, t2)
t1 = dble(nint(t1))
t2 = dble(nint(t2))
       tt = t1 + 0.5d0 * (xi + 1.0d0) * (t2 - t1)

       call spline_nurbs_eval2(tt, xs, ys, tcurve%a, tcurve%b, tcurve%c &
                   , tcurve%Mx, tcurve%My, t, y1x, y1y, 'interp', tcurve%btype)
    end if

    ! finalize the transformation
    y3x = xx(4) + 0.5d0 * (xi + 1.0d0) * (xx(3) - xx(4))
    y3y = yy(4) + 0.5d0 * (xi + 1.0d0) * (yy(3) - yy(4))

    x = (eta + 1.0d0) / 2.0d0 * y3x - (eta - 1.0d0) / 2.0d0 * y1x
    y = (eta + 1.0d0) / 2.0d0 * y3y - (eta - 1.0d0) / 2.0d0 * y1y
  
    ! done here
  end subroutine transform_quadri

  ! 
  subroutine alloc_init_loc_matrices(elem, npe, neqs)
    implicit none
    class(element_dg2d), intent(inout) :: elem
    integer, intent(in) :: npe, neqs 

    ! local vars
    integer :: ii

    if ( allocated(elem%Mass) ) deallocate(elem%Mass)       
    allocate(elem%Mass(neqs * npe, neqs * npe))
    elem%Mass = 0.0d0

    if ( allocated(elem%LUmass) ) deallocate(elem%LUmass)       
    allocate(elem%LUmass(neqs * npe, neqs * npe))
    elem%LUmass = 0.0d0

    if ( allocated(elem%LUmass_imp) ) deallocate(elem%LUmass_imp)       
    allocate(elem%LUmass_imp(neqs * npe, neqs * npe))
    elem%LUmass_imp = 0.0d0

    if ( allocated(elem%IPIVmass) ) deallocate(elem%IPIVmass)       
    allocate(elem%IPIVmass(neqs * npe))
    elem%IPIVmass = (/ (ii, ii = 1, (neqs * npe) ) /)

    if ( allocated(elem%IPIVmass_imp) ) deallocate(elem%IPIVmass_imp)       
    allocate(elem%IPIVmass_imp(neqs * npe))
    elem%IPIVmass_imp = (/ (ii, ii = 1, (neqs * npe) ) /)

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

    if ( allocated(elem%Un) ) deallocate(elem%Un)
    allocate(elem%Un(neqs, npe))
    elem%Un = 0.0d0

    if ( allocated(elem%Urk) ) deallocate(elem%Urk)
    allocate(elem%Urk(neqs, npe, 2))
    elem%Urk = 0.0d0

    if ( allocated(elem%U0) ) deallocate(elem%U0)
    allocate(elem%U0(neqs, npe))
    elem%U0 = 0.0d0

    if ( allocated(elem%Ax) ) deallocate(elem%Ax)
    allocate(elem%Ax(neqs, npe))
    elem%Ax = 0.0d0

    if ( allocated(elem%d_psi_d_x) ) deallocate(elem%d_psi_d_x)
    allocate(elem%d_psi_d_x(npe, elem%ngauss, 2))
    elem%d_psi_d_x = 0.0d0


    if ( allocated(elem%Fk) ) deallocate(elem%Fk)
    allocate(elem%Fk(neqs, elem%ngauss, 2))
    elem%Fk = 0.0d0

    if ( allocated(elem%dFk) ) deallocate(elem%dFk)
    allocate(elem%dFk(neqs, neqs, elem%ngauss, 2))
    elem%dFk = 0.0d0

    ! done here
  end subroutine alloc_init_loc_matrices

  ! completely cleans/deallocates
  ! a neighbor struct
  !
  subroutine dealloc_neigh(tneigh)
    implicit none
    class(neigh_dg), intent(inout) :: tneigh

    if(allocated(tneigh%xi)) deallocate(tneigh%xi)
    if(allocated(tneigh%W)) deallocate(tneigh%W)
    if(allocated(tneigh%xloc_in)) deallocate(tneigh%xloc_in)
    if(allocated(tneigh%xloc_out)) deallocate(tneigh%xloc_out)
    if(allocated(tneigh%x)) deallocate(tneigh%x)
    if(allocated(tneigh%dx)) deallocate(tneigh%dx)
    if(allocated(tneigh%n)) deallocate(tneigh%n)
    if(allocated(tneigh%s)) deallocate(tneigh%s)
    if(allocated(tneigh%psi_in)) deallocate(tneigh%psi_in)
    if(allocated(tneigh%Fstar)) deallocate(tneigh%Fstar)
    if(allocated(tneigh%dFpm)) deallocate(tneigh%dFpm)

    ! done here
  end subroutine dealloc_neigh

  ! completely cleans/deallocates
  ! an edg_dg struct
  !
  subroutine dealloc_edg(tedg)
    implicit none
    class(edg_dg), intent(inout) :: tedg

    ! local vars
    integer :: i

    if(allocated(tedg%pts)) deallocate(tedg%pts)

    if (allocated(tedg%neighs)) then
       do i = 1, size(tedg%neighs)
          call tedg%neighs(i)%dealloc_neigh()
       end do
    end if

    if (allocated(tedg%neighs)) deallocate(tedg%neighs)

    ! done here
  end subroutine dealloc_edg

  ! completely deallocates/cleans
  ! an element/hull
  !
  subroutine dealloc_elem(elem)
    implicit none
    class(element_dg2d), intent(inout) :: elem

    ! local vars
    integer :: i

    if(allocated(elem%r)) deallocate(elem%r)
    if(allocated(elem%s)) deallocate(elem%s)
    if(allocated(elem%W)) deallocate(elem%W)

    call elem%tbasis%dealloc_basis()

    if(allocated(elem%psi)) deallocate(elem%psi)
    if(allocated(elem%d_psi_d_xi)) deallocate(elem%d_psi_d_xi)
    if(allocated(elem%d_psi_d_eta)) deallocate(elem%d_psi_d_eta)

    if(allocated(elem%jac)) deallocate(elem%jac)
    if(allocated(elem%Jstar)) deallocate(elem%Jstar)
    if(allocated(elem%JJ)) deallocate(elem%JJ)

    if(allocated(elem%K)) deallocate(elem%K)
    if(allocated(elem%M)) deallocate(elem%M)
    if(allocated(elem%Q)) deallocate(elem%Q)
    if(allocated(elem%f)) deallocate(elem%f)

    if ( associated(elem%A)) nullify(elem%A)

    if(allocated(elem%x)) deallocate(elem%x)
    if(allocated(elem%x_loc)) deallocate(elem%x_loc)

    if(allocated(elem%Mass)) deallocate(elem%Mass)
    if(allocated(elem%LUmass)) deallocate(elem%LUmass)
    if(allocated(elem%LUmass_imp)) deallocate(elem%LUmass_imp)

    if(allocated(elem%IPIVmass)) deallocate(elem%IPIVmass)
    if(allocated(elem%IPIVmass_imp)) deallocate(elem%IPIVmass_imp)

    if(allocated(elem%U)) deallocate(elem%U)
    if(allocated(elem%Us)) deallocate(elem%Us)
    if(allocated(elem%rhs)) deallocate(elem%rhs)
    if(allocated(elem%So)) deallocate(elem%So)
    if(allocated(elem%Un)) deallocate(elem%Un)
    if(allocated(elem%Urk)) deallocate(elem%Urk)
    if(allocated(elem%d_psi_d_x)) deallocate(elem%d_psi_d_x)
    if(allocated(elem%Fk)) deallocate(elem%Fk)
    if(allocated(elem%dFk)) deallocate(elem%dFk)
    if(allocated(elem%U0)) deallocate(elem%U0)
    if(allocated(elem%Ax)) deallocate(elem%Ax)

    if ( allocated(elem%edgs) ) then
       do i = 1, size(elem%edgs)
          call elem%edgs(i)%dealloc_edg()
       end do
    end if

    if ( allocated(elem%edgs) ) deallocate(elem%edgs)

    ! done here
  end subroutine dealloc_elem

end module element_opt_dg2d
