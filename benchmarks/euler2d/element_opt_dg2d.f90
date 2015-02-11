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
     real*8, dimension(:, :), allocatable :: x !physic. coords
     ! (1..2, 1..npe) (1, npe) = xi <> (2, npe) = eta
     real*8, dimension(:, :), allocatable :: x_loc !local coords
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
     procedure, public :: comp_flux_interior
     procedure, public :: comp_metric_jacobian => comp_Jstar_point_dg
     procedure, private :: comp_dpsi => comp_d_psi_d_x
     procedure, public :: comp_int_integ => comp_inter_integral
     procedure, public :: xy2rs => xy2rs_dg
     procedure, public :: init_edg_quadrat => init_elem_edg_quadratures
     procedure, public :: comp_Fstar
  end type element_dg2d

  type bc

     character(len = 128) :: name
     ! one value per tag at this version
     ! val(1:neqs)
     real*8, dimension(:), allocatable :: val 

  end type bc

  type dg_srtuct

     type(grid) :: grd
     type(element_dg2d), dimension(:), allocatable :: elems
     type(bc), dimension(:), allocatable :: bcs ! zero-based

  end type dg_srtuct


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

  subroutine init_elem_edg_quadratures(elem, wspace, iedg, tol)
    implicit none
    class(element_dg2d), intent(inout), target :: elem
    class(dg_srtuct), intent(inout), target :: wspace ! work space
    integer, intent(in) :: iedg !edge number
    real*8, intent(in) :: tol

    ! local vars
    integer :: i, j, max_npedg, r, tag
    type(edg_dg), pointer :: tedg => null()
    type(neigh_dg), pointer :: tneigh => null()
    real*8, parameter :: alpha = 0.0d0, beta  = 0.0d0
    logical :: do_snap = .false.
    type(curve), pointer :: tcurve=>null()
    real*8, dimension(:), pointer :: t=>null(), xs=>null(), ys=>null()
    real*8 :: tt, x1, y1, x2, y2, xq, yq, xdot, ydot
    type(grid), pointer :: grd => null()
    real*8 :: t1, t2
    real*8, dimension(:), allocatable :: derjac

    ! HARD reset
    do_snap = .false.

    grd => wspace%grd
    tedg => elem%edgs(iedg)
    tag = tedg%tag

    if ( tag .ne. 0 ) then
       ! is a boundary element 
       do_snap = .true.
       tcurve => grd%bn_curves(tag)
       t => tcurve%t
       xs => tcurve%x
       ys => tcurve%y
    end if

    do j = 1, size(tedg%neighs) ! loop over all neighbors to that edge (shared segments)

       tneigh => tedg%neighs(j) ! select this neighbor 

       ! find the maximum number of 1d interpolation (Lagrange) points per this shared 
       ! segment with the j^th neighbor of this edge.
       !
       ! HINT : use npedg info to see how many points per edge we have in this elem
       ! and this neighbor. and take the maximum as 
       ! the highest possible degree of a 1d polynomial that can be defined on that
       ! part of the shared segment 
       ! 
       if (tneigh%elnum .eq. -1)  then !wall
          max_npedg = elem%npedg ! count on element itself
       else
          max_npedg = max(elem%npedg, wspace%elems(tneigh%elnum)%npedg)
       end if

       ! degree of exactness relation for 1d Gauss legendre: 2r - 1 = max_npedg -1
       r = nint(1.5d0 * dble(max_npedg) / dble(2)) ! 1.5 is safety factor :)
       tneigh%ngpseg = r

       ! allocate neigh specific arrays
       allocate(tneigh%xi(r), tneigh%W(r))
       allocate(tneigh%xloc_in(2, r), tneigh%xloc_out(2, r))
       allocate(tneigh%x(2, r), tneigh%dx(2, r))
       allocate(tneigh%n(2, r), tneigh%s(r))
       ! alloc/init Fstar(1:neqs, 1:ngpseg) 
       allocate(tneigh%Fstar(elem%neqs, r))

       tneigh%xi = 0.0d0; tneigh%W = 0.0d0
       tneigh%xloc_in = 0.0d0; tneigh%xloc_out = 0.0d0
       tneigh%x = 0.0d0; tneigh%dx = 0.0d0
       tneigh%Fstar = 0.0d0

       ! computing Legendre-Gauss-Jacobi points for integration
       ! and corresponding weight functions
       allocate(derjac(r))
       call ZEJAGA(r, alpha, beta, tneigh%xi, derjac)
       call WEJAGA(r, alpha, beta, tneigh%xi, derjac, tneigh%W)
       deallocate(derjac)

       ! 1 is the start and 2 is the end of this edge segment; just for convention
       x1 = tneigh%xs(1); x2 = tneigh%xe(1) 
       y1 = tneigh%xs(2); y2 = tneigh%xe(2) 

       ! find the constant derivatives
       xdot = 0.5d0 * (x2 - x1)
       ydot = 0.5d0 * (y2 - y1)

       do i = 1, r

          ! Mapping of coordinates of Gauss-Legendre quadrature
          ! for straight edges to physical space
          xq = x1 + (tneigh%xi(i) + 1.0d0) / 2.0d0 * (x2 - x1) 
          yq = y1 + (tneigh%xi(i) + 1.0d0) / 2.0d0 * (y2 - y1) 

          ! if boundary edge then also snapp <x> and <dx>
          if ( do_snap ) then
             call find_t(grd, tag, x1, y1, tol, t1)
             if ( t1 .eq. -1.0d0) then
                print *, 'error : the parameter t can not be computed' &
                     , ' for point (', x1, ',', y1,') on ' &
                     , 'edg # ', iedg,'. stop.'
                stop
             end if
             call find_t(grd, tag, x2, y2, tol, t2)
             if ( t2 .eq. -1.0d0) then
                print *, 'error : the parameter t can not be computed' &
                     , ' for point (', x2, ',', y2,') on ' &
                     , 'edg # ', iedg,'. stop.'
                stop
             end if

             call find_t(grd, tag, xq, yq, tol, tt)
             if ( tt .eq. -1.0d0) then
                print *, 'error : the parameter t can not be computed' &
                     , ' for point (', xq, ',', yq,') on ' &
                     , 'edg # ', iedg,'. stop.'
                stop
             end if
             ! find new snapped xq and yq on the curve
             xq = 0.0d0; yq = 0.0d0
             call spline_nurbs_eval2(tt, xs, ys, tcurve%a, tcurve%b, tcurve%c, &
                  & tcurve%Mx, tcurve%My, t, xq, yq, 'interp', tcurve%btype)
  
             ! computing derivatives xdot and ydot ...
             xdot = 0.0d0; ydot = 0.0d0
             call spline_nurbs_eval2(tt, xs, ys, tcurve%a, tcurve%b, tcurve%c, &
                  & tcurve%Mx, tcurve%My, t, xdot, ydot, 'diff1', tcurve%btype)

             ! add stuff coming from chain rule for differentiation
             xdot = xdot * 0.5d0 * (t2 - t1)
             ydot = ydot * 0.5d0 * (t2 - t1)

          end if

          ! store physical coords and derivatives
          tneigh%x(1,i) = xq; tneigh%x(2,i) = yq
          tneigh%dx(1,i) = xdot; tneigh%dx(2,i) = ydot
          tneigh%s(i) = sqrt(xdot**2 + ydot**2)
          tneigh%n(1, i) = ydot / tneigh%s(i)
          tneigh%n(2, i) = -xdot / tneigh%s(i) 

          ! convert physical (xq, yq) to local coords
          ! (r,s) in this element store in this element
          !
          ! xy2rs_dg(elem, x, y, maxitr, tolrs, r, s)
          call elem%xy2rs(xq , yq, 40, tol &
               , tneigh%xloc_in(1,i), tneigh%xloc_in(2,i))
          ! (r, s) in the neighbor element stored in this element
          if (tneigh%elnum .ne. -1)  then
             call wspace%elems(tneigh%elnum)%xy2rs(xq , yq &
                  , 40, tol, tneigh%xloc_out(1,i), tneigh%xloc_out(2,i))
          end if

       end do ! next quadrature point per the current segment shared with current neighbor

    end do ! next neighbor per the current edge
 
    ! done here
  end subroutine init_elem_edg_quadratures

  ! computes the on edge Rankine–Hugoniot value approximated 
  ! by either a Rieman solver or flux splitting algorithm
  ! and store it in tneigh%Fstar
  !
  subroutine comp_Fstar(elem, wspace, tedg, tneigh)
    implicit none
    class(element_dg2d) :: elem
    class(dg_srtuct) :: wspace ! work space
    type(edg_dg) :: tedg
    type(neigh_dg) :: tneigh

    ! local vars
    integer :: k
    real*8 :: r, s, nx, ny
    real*8, dimension(elem%neqs) :: UL, UR, FP, FM
    logical, dimension(4) :: f_select
    real*8, dimension(elem%neqs) :: fvl_p, fvl_m
    real*8, dimension(elem%neqs, elem%neqs) :: d_fvl_p, d_fvl_m

    ! loop over Gauss points per this neighboring segment
    do k = 1, tneigh%ngpseg

       ! find outgoing unit normals at that gauss points
       nx = tneigh%n(1, k); ny = tneigh%n(2, k)

       !
       !          <<< compute UL and UR procedure >>>
       ! 
       ! evaluate UL at the current Gauss point using interior (r, s) 
       ! coordinates and the basis function of this element
       r = tneigh%xloc_in(1, k); s = tneigh%xloc_in(2, k) 
       call elem%comp_u(r, s, UL)

       ! now decide on UR ... 
       select case (wspace%bcs(tedg%tag)%name)

       case ('interior') ! then UR is in the other neighboring element

          ! evaluate UR at the current Gauss point using 
          ! the neighboring element's local coordinates (r, s)
          ! and the basis functions of the neighboring element
          r = tneigh%xloc_out(1, k); s = tneigh%xloc_out(2, k)  
          call wspace%elems(tneigh%elnum)%comp_u(r, s, UR)

       case ('inflow', 'outflow') ! boundary edge; then UR is preset in bcs value

          UR = wspace%bcs(tedg%tag)%val

       case default

          ! do nothing for now!

       end select

       !        <<< compute Flux procedure >>>
       ! 
       !      special situation : Wall treatment
       if (wspace%bcs(tedg%tag)%name .eq. 'wall') then

          call calc_wall_flux(UL, elem%neqs, elem%gamma, nx, ny, fvl_p, d_fvl_p)

          ! store F+
          FP = fvl_p
          FM = 0.0d0

       else

          ! compute F+ using UL
          f_select = (/ .true. , .false., .false., .false. /)
          call calc_van_leer(UL, elem%neqs, elem%gamma, nx, ny, f_select &
               , fvl_p, fvl_m, d_fvl_p, d_fvl_m)

          ! store F+
          FP = fvl_p

          ! compute F- using UR
          f_select = (/ .false. , .true., .false., .false. /)
          call calc_van_leer(UR, elem%neqs, elem%gamma, nx, ny, f_select &
               , fvl_p, fvl_m, d_fvl_p, d_fvl_m)

          ! store F-
          FM = fvl_m

       end if


       ! store total split flux as the final Rankine–Hugoniot value at Fstar
       tneigh%Fstar(:, k) = FP + FM

    end do ! next gauss point per neighboring element (on shared segment)

    ! done here
  end subroutine comp_Fstar


end module element_opt_dg2d
