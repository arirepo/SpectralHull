module euler2d_eqs
  implicit none

  private


  public :: calc_van_leer, calc_wall_flux, calc_pure_euler2d_flux
  public :: u2U, calc_pure_euler2d_flux_jac
  public :: comp_char_bcs_euler2d, comp_weakly_proj_wall_U

contains

  ! f_select[1] = .false. -> not compute f+,  = .true. -> compute f+
  ! f_select[2] = .false. -> not compute f-,  = .true. -> compute f-
  ! f_select[3] = .false. -> not compute df+, = .true. -> compute df+
  ! f_select[4] = .false. -> not compute df-, = .true. -> compute df-
  subroutine calc_van_leer(Q, neqs, gamma, nx, ny, f_select &
       , fvl_p, fvl_m, d_fvl_p, d_fvl_m)
    implicit none
    real*8, dimension(:), intent(in) :: Q
    integer, intent(in) :: neqs
    real*8, intent(in) :: gamma, nx, ny
    logical, dimension(:), intent(in) :: f_select
    real*8, dimension(:), intent(out) :: fvl_p, fvl_m
    real*8, dimension(:, :), intent(out) :: d_fvl_p, d_fvl_m


    ! local vars
    integer :: i, j
    real*8 :: sign_pm = 1.0d0, p2 = 0.0d0
    real*8 :: rho, u, v, u_bar, e, P, c
    real*8, dimension(neqs) :: d_rho_dQi, d_u_dQi, d_v_dQi &
         , d_e_dQi, d_P_dQi, d_c_dQi, d_ubar_dQi

    ! HARD reset
    fvl_p = 0.0d0; fvl_m = 0.0d0; d_fvl_p = 0.0d0; d_fvl_m = 0.0d0
    d_rho_dQi = 0.0d0; d_u_dQi = 0.0d0; d_v_dQi = 0.0d0
    d_e_dQi = 0.0d0; d_P_dQi = 0.0d0; d_c_dQi = 0.0d0
    d_ubar_dQi = 0.0d0


    ! step 1: calculating primitive variables
    rho = Q(1)
    u = Q(2) / Q(1)
    v = Q(3) / Q(1)
    e = Q(4)
    u_bar = u * nx + v * ny
    P = (gamma - 1.0d0) * e - 0.5d0 * (gamma - 1.0d0) *rho * ( u*u + v*v)
    c = sqrt(gamma * P/ rho)

    ! step 2: filling the flux van leer (fvl) vector ASAP
    if (u_bar >= c) then

       if(f_select(1) .or. f_select(3)) then

          fvl_p(1) = rho* u_bar
          fvl_p(2) = (rho * u * u_bar + nx * P)
          fvl_p(3) = (rho * v * u_bar + ny * P)
          fvl_p(4) = (e+P)* u_bar
       end if
       if(f_select(2) .or. f_select(4)) then

          fvl_m = 0.0d0
       end if

    elseif (u_bar <= (-c)) then

       if (f_select(2) .or. f_select(4)) then
          fvl_m(1) = rho* u_bar
          fvl_m(2) = (rho * u * u_bar + nx * P)
          fvl_m(3) = (rho * v * u_bar + ny * P)
          fvl_m(4) = (e+P)* u_bar
       end if

       if(f_select(1) .or. f_select(3)) then

          fvl_p= 0.0d0
       end if

    else

       if(f_select(1) .or. f_select(3)) then

          fvl_p(1) = 0.25d0 * rho * c * (u_bar/c + 1.0d0)**2
          fvl_p(2) = fvl_p(1) * (nx/gamma*(-u_bar + 2.0d0 * c) + u)
          fvl_p(3) = fvl_p(1) * (ny/gamma*(-u_bar + 2.0d0 * c) + v)
          fvl_p(4) = fvl_p(1) * ( (-(gamma-1.0d0)*u_bar*u_bar &
               + 2.0d0*(gamma-1.0d0)*u_bar*c+2.0d0*c*c) &
               / (gamma*gamma - 1.0d0) + (u*u + v*v)*0.5d0 )
       end if
       if(f_select(2) .or. f_select(4)) then
          fvl_m(1) = -0.25d0 * rho*c* (u_bar/c - 1.0d0)**2
          fvl_m(2) = fvl_m(1) * (nx/gamma*(-u_bar - 2.0d0 * c) + u)
          fvl_m(3) = fvl_m(1) * (ny/gamma*(-u_bar - 2.0d0 * c) + v)
          fvl_m(4) = fvl_m(1) * ( (-(gamma-1.0d0)*u_bar*u_bar &
               - 2.0d0*(gamma-1.0d0)*u_bar*c+2.0d0*c*c)/ (gamma*gamma - 1.0d0) &
               + (u*u + v*v)*0.5d0 )
       end if

    end if

    !at least one of the Jacobians are needed to be computed
    if(f_select(3) .or. f_select(4)) then
       !step 3: calculating derivatives of primitive variables
       d_rho_dQi(1) = 1.0d0;   d_rho_dQi(2) = 0.0d0;   d_rho_dQi(3) = 0.0d0
       d_rho_dQi(4) = 0.0d0
       d_u_dQi(1) = -Q(2)/(Q(1)**2); d_u_dQi(2) = 1.0d0/Q(1); d_u_dQi(3) = 0.0d0
       d_u_dQi(4) = 0.0d0
       d_v_dQi(1) = -Q(3)/(Q(1)**2); d_v_dQi(2) = 0.0d0; d_v_dQi(3) = 1.0d0/Q(1)
       d_v_dQi(4) = 0.0d0
       d_e_dQi(1) = 0.0d0; d_e_dQi(2) = 0.0d0; d_e_dQi(3) = 0.0d0
       d_e_dQi(4) = 1.0d0

       !step 4: calculating derivatives of physical variables P, c, u_bar

       do i = 1, neqs
          d_P_dQi(i) = (gamma-1.0d0) * (d_e_dQi(i) - 0.50d0 * d_rho_dQi(i) &
               *(u*u + v*v) - rho * (u * d_u_dQi(i) + v * d_v_dQi(i)))
          d_c_dQi(i) = gamma/2.0d0 * (d_P_dQi(i) * rho - d_rho_dQi(i) * P) &
               / (rho * rho * c)
          d_ubar_dQi(i) = d_u_dQi(i) * nx + d_v_dQi(i) * ny
       end do

       ! step 5: using chain-rule to finish the job.
       ! starting filling row by row ..
       if(u_bar >= c) then
          if (f_select(3) ) then
             ! d_fvl_p
             do j = 1, neqs
                d_fvl_p(1, j) = (d_rho_dQi(j)*u_bar + rho*d_ubar_dQi(j))
                d_fvl_p(2, j) = (d_rho_dQi(j) * u * u_bar+ rho * d_u_dQi(j) &
                     * u_bar+ rho * u * d_ubar_dQi(j) + nx * d_P_dQi(j))
                d_fvl_p(3, j) = (d_rho_dQi(j) * v * u_bar+ rho * d_v_dQi(j) &
                     * u_bar+ rho * v * d_ubar_dQi(j) + ny * d_P_dQi(j))
                d_fvl_p(4, j) = ( (d_e_dQi(j)+d_P_dQi(j))* u_bar + (e+P) &
                     * d_ubar_dQi(j) )
             end do
          end if
          if(f_select(4)) then
             ! d_fvl_m
             d_fvl_m = 0.0d0
          end if

       elseif (u_bar <= (-c)) then
          if( f_select(3) ) then
             ! d_fvl_p
             d_fvl_p = 0.0d0
          end if

          if(f_select(4)) then

             ! d_fvl_m
             do j = 1, neqs
                d_fvl_m(1, j) = (d_rho_dQi(j)*u_bar + rho*d_ubar_dQi(j))
                d_fvl_m(2, j) = (d_rho_dQi(j) * u * u_bar+ rho * d_u_dQi(j) &
                     * u_bar+ rho * u * d_ubar_dQi(j) + nx * d_P_dQi(j))
                d_fvl_m(3, j) = (d_rho_dQi(j) * v * u_bar+ rho * d_v_dQi(j) &
                     * u_bar+ rho * v * d_ubar_dQi(j) + ny * d_P_dQi(j))
                d_fvl_m(4, j) = ( (d_e_dQi(j)+d_P_dQi(j))* u_bar + (e+P) &
                     * d_ubar_dQi(j) )
             end do
          end if

       else
          if( f_select(3) ) then
             ! d_fvl_p
             sign_pm = 1.0d0
             p2 = (u_bar / c + sign_pm * 1.0d0)**2
             do j = 1, neqs
                d_fvl_p(1, j) = sign_pm*0.25d0*(d_rho_dQi(j)*c*p2  &
                     + rho * d_c_dQi(j) * p2 + 2.0d0*(rho/c) &
                     * (u_bar/c + sign_pm * 1.0d0) &
                     * (d_ubar_dQi(j)*c-d_c_dQi(j)*u_bar))
                d_fvl_p(2, j) = d_fvl_p(1, j) * (nx/gamma &
                     * (-u_bar + sign_pm*2.0d0*c)+u) + fvl_p(1) &
                     * (nx/gamma*(-d_ubar_dQi(j)+sign_pm* 2.0d0* d_c_dQi(j)) &
                     + d_u_dQi(j))
                d_fvl_p(3, j) = d_fvl_p(1, j) * (ny/gamma &
                     * (-u_bar + sign_pm*2.0d0*c)+v) + fvl_p(1) &
                     * (ny/gamma*(-d_ubar_dQi(j)+sign_pm* 2.0d0 &
                     * d_c_dQi(j))+ d_v_dQi(j))
                d_fvl_p(4, j) = d_fvl_p(1, j) * ((-(gamma-1.0d0) &
                     *u_bar*u_bar+ sign_pm*2.0d0*(gamma-1.0d0)*u_bar*c+2.0d0*c*c) &
                     /(gamma*gamma-1.0d0)+0.50d0*(u*u+v*v)) + fvl_p(1) &
                     * ((-2.0d0*(gamma-1.0d0)*u_bar*d_ubar_dQi(j) &
                     +sign_pm*2.0d0*(gamma-1.0d0)*(d_ubar_dQi(j)*c &
                     +u_bar*d_c_dQi(j))+ 4.0d0*c*d_c_dQi(j)) &
                     /(gamma*gamma-1.0d0)+ u*d_u_dQi(j) + v*d_v_dQi(j))
             end do
          end if

          if( f_select(4) ) then
             ! d_fvl_m
             sign_pm = -1.0d0
             p2 = (u_bar / c + sign_pm * 1.0d0)**2
             do j = 1, neqs
                d_fvl_m(1, j) = sign_pm*0.25d0*(d_rho_dQi(j)*c*p2 &
                     + rho * d_c_dQi(j) * p2 + 2.0d0*(rho/c) &
                     * (u_bar/c + sign_pm * 1.0d0) * (d_ubar_dQi(j) &
                     *c-d_c_dQi(j)*u_bar))
                d_fvl_m(2, j) = d_fvl_m(1, j) * (nx/gamma &
                     * (-u_bar + sign_pm*2.0d0*c)+u) + fvl_m(1) &
                     * (nx/gamma*(-d_ubar_dQi(j)+sign_pm &
                     * 2.0d0* d_c_dQi(j))+ d_u_dQi(j))
                d_fvl_m(3, j) = d_fvl_m(1, j) * (ny/gamma &
                     * (-u_bar + sign_pm*2.0d0*c)+v) + fvl_m(1) &
                     * (ny/gamma*(-d_ubar_dQi(j)+sign_pm* 2.0d0 &
                     * d_c_dQi(j))+ d_v_dQi(j))
                d_fvl_m(4, j) = d_fvl_m(1, j) * ((-(gamma-1.0d0) &
                     *u_bar*u_bar+ sign_pm*2.0d0*(gamma-1.0d0) &
                     *u_bar*c+2.0d0*c*c)/(gamma*gamma-1.0d0) &
                     +0.5d0*(u*u+v*v)) + fvl_m(1) &
                     * ((-2.0d0*(gamma-1.0d0)*u_bar*d_ubar_dQi(j) &
                     +sign_pm*2.0d0*(gamma-1.0d0)*(d_ubar_dQi(j)*c &
                     +u_bar*d_c_dQi(j))+ 4.0d0*c*d_c_dQi(j)) &
                     /(gamma*gamma-1.0d0)+ u*d_u_dQi(j) + v*d_v_dQi(j))
             end do
          end if

       end if

    end if

    ! done here
  end subroutine calc_van_leer


  !computes wall flux and Jacobians
  subroutine calc_wall_flux(Q, neqs, gamma, nx, ny, fw, d_fw)
    implicit none
    real*8, dimension(:), intent(in) :: Q
    integer, intent(in) :: neqs
    real*8, intent(in) :: gamma, nx, ny
    real*8, dimension(:), intent(out) :: fw
    real*8, dimension(:, :), intent(out) :: d_fw

    ! local vars
    integer :: i, j
    real*8 :: rho, u, v, u_bar, e, P, c
    real*8, dimension(neqs) :: d_rho_dQi, d_u_dQi, d_v_dQi, d_e_dQi, d_P_dQi

    ! HARD reset
    fw = 0.0d0; d_fw = 0.0d0
    d_rho_dQi = 0.0d0; d_u_dQi = 0.0d0; d_v_dQi = 0.0d0
    d_e_dQi = 0.0d0; d_P_dQi = 0.0d0

    ! step 1: calculating primitive variables
    rho = Q(1)
    u = Q(2) / Q(1)
    v = Q(3) / Q(1)
    e = Q(4)
    u_bar = u * nx + v * ny
    P = (gamma - 1.0d0) * e - 0.5d0 * (gamma - 1.0d0) *rho * ( u*u + v*v)
    c = sqrt(gamma * P/ rho)

    ! calculate wall flux here
    fw(1) = 0.0d0
    fw(2) = (nx * P)
    fw(3) = (ny * P)
    fw(4) = 0.0d0

    ! calculating derivatives of primitive variables
    d_rho_dQi(1) = 1.0d0;   d_rho_dQi(2) = 0.0d0;   d_rho_dQi(3) = 0.0d0;   d_rho_dQi(4) = 0.0d0
    d_u_dQi(1) = -Q(2)/(Q(1)**2); d_u_dQi(2) = 1.0d0/Q(1); d_u_dQi(3) = 0.0d0; d_u_dQi(4) = 0.0d0
    d_v_dQi(1) = -Q(3)/(Q(1)**2); d_v_dQi(2) = 0.0d0; d_v_dQi(3) = 1.0d0/Q(1); d_v_dQi(4) = 0.0d0
    d_e_dQi(1) = 0.0d0; d_e_dQi(2) = 0.0d0; d_e_dQi(3) = 0.0d0; d_e_dQi(4) = 1.0d0

    ! calculating derivatives of P
    do i = 1, neqs
       d_P_dQi(i) = (gamma-1.0d0) * (d_e_dQi(i) - 0.5d0 * d_rho_dQi(i) *(u*u + v*v) - rho * (u * d_u_dQi(i) + v * d_v_dQi(i)))
    end do

    ! filling d_fw
    do j = 1, neqs
       d_fw(1, j) = 0.0d0
       d_fw(2, j) = nx * d_P_dQi(j)
       d_fw(3, j) = ny * d_P_dQi(j)
       d_fw(4, j) = 0.0d0
    end do

    ! done here
  end subroutine calc_wall_flux

  subroutine calc_pure_euler2d_flux(Q, gamma, F, G)
    implicit none
    real*8, dimension(:), intent(in) :: Q
    real*8, intent(in) :: gamma
    real*8, dimension(:), intent(out) :: F, G

    ! local vars
    real*8 :: rho, u, v, e, P

    ! calculating primitive variables
    rho = Q(1)
    u = Q(2) / Q(1)
    v = Q(3) / Q(1)
    e = Q(4)
    P = (gamma - 1.0d0) * e - 0.5d0 * (gamma - 1.0d0) *rho * ( u*u + v*v)

    ! computing fluxes
    F(1) = rho * u
    F(2) = rho * u * u + P
    F(3) = rho * u * v
    F(4) = (e + P) * u

    G(1) = rho * v
    G(2) = rho * v * u
    G(3) = rho * v * v + P
    G(4) = (e + P) * v

    ! done here
  end subroutine calc_pure_euler2d_flux

  !< @brief
  !! Sets the parameters of Euler equations of compressible fluid dynamics in 2D
  !< @todo Check the description here!
  !! Converts u to U
  !! Here \f$ \rho\ \f$ is the density,
  !! (u,v) is the velocity.
  !! /f$ \gamma \f$ is the ratio of specific heats.
  !< @param \f$ \rho, \gamma, u, v, and P input \f$
  !< @ param UU output array!

  subroutine u2U(rho, u, v, P, gamma, UU)
    implicit none
    ! primitive vars
    real*8, intent(in) :: rho !< \f$ \rho \f$ is the density
    real*8, intent(in) :: u   !< velocity component
    real*8, intent(in) :: v   !< velocity component
    real*8, intent(in) :: P   !<
    real*8, intent(in) :: gamma !< \f$ \gamma \f$ is the ratio of specific heats.

    ! conservative vars
    real*8, dimension(:), intent(out) :: UU !< Output of the subroutine

    UU(1) = rho
    UU(2) = rho * u
    UU(3) = rho * v
    UU(4) = P / (gamma - 1.0d0)  + 0.5d0 * rho * (u**2 + v**2)

    ! done here
  end subroutine u2U

  subroutine calc_pure_euler2d_flux_jac(Q, gamma, dF, dG)
    implicit none
    real*8, dimension(:), intent(in) :: Q
    real*8, intent(in) :: gamma
    real*8, dimension(:, :), intent(out) :: dF, dG

    ! local vars
    real*8 :: rho, u, v, e, zeta, u2

    ! calculating primitive variables
    rho = Q(1)
    u = Q(2) / Q(1)
    v = Q(3) / Q(1)
    e = Q(4)

    ! computing repeated factors
    u2 = u**2 + v**2
    zeta = (gamma - 1.0d0) * u2 - e * gamma / rho

    ! computing pure flux Jacobian dF/dQ
    dF(1, 1) = 0.0d0
    dF(1, 2) = 1.0d0
    dF(1, 3) = 0.0d0
    dF(1, 4) = 0.0d0

    dF(2, 1) = -1.0d0 * (u**2) + 0.5d0 * (gamma - 1.0d0) * u2
    dF(2, 2) = (3.0d0 - gamma) * u
    dF(2, 3) = (1.0d0 - gamma) * v
    dF(2, 4) = gamma - 1.0d0


    dF(3, 1) = -1.0d0 * u * v
    dF(3, 2) = v
    dF(3, 3) = u
    dF(3, 4) = 0.0d0

    dF(4, 1) = zeta * u
    dF(4, 2) = gamma * e / rho + 0.5d0 * (1.0d0 - gamma) * u2 &
         + (1.0d0 - gamma) * (u**2)
    dF(4, 3) = (1.0d0 - gamma) * u * v
    dF(4, 4) = gamma * u

    ! computing pure flux Jacobian dG/dQ
    dG(1, 1) = 0.0d0
    dG(1, 2) = 0.0d0
    dG(1, 3) = 1.0d0
    dG(1, 4) = 0.0d0

    dG(2, 1) = -1.0d0 * u * v
    dG(2, 2) = v
    dG(2, 3) = u
    dG(2, 4) = 0.0d0

    dG(3, 1) = -1.0d0 * (v**2) + 0.5d0 * (gamma - 1.0d0) * u2
    dG(3, 2) = (1.0d0 - gamma) * u
    dG(3, 3) = (3.0d0 - gamma) * v
    dG(3, 4) = gamma - 1.0d0

    dG(4, 1) = zeta * v
    dG(4, 2) = (1.0d0 - gamma) * u * v
    dG(4, 3) = gamma * e / rho + 0.5d0 * (1.0d0 - gamma) * u2 &
         + (1.0d0 - gamma) * (v**2)
    dG(4, 4) = gamma * v

    ! done here
  end subroutine calc_pure_euler2d_flux_jac

  ! compute the value of ghost points (AKA F(-) in DG)
  ! using freezed LODI characteristic decomposition of
  ! two dimensional compressible Euler equations.
  !
  ! The freezed values are evaluated at boundary node
  ! current value at the current time step.
  !
  ! input
  ! Q = [rho, rho * u, rho *v, e] at boundary node
  !     at the current time step
  ! Qinf = [rho_inf, rho_inf * u_inf, rho_inf *v_inf, e_inf]
  !     of the infinity refrence values.
  ! (nx, ny) = unit outward normal vector
  ! output
  ! Qghost = [rho, rho *u, rho *v, e]@ ghost point
  !
  subroutine comp_char_bcs_euler2d(Q, Qinf, gamma, nx, ny, Qghost)
    implicit none
    real*8, dimension(:), intent(in) :: Q, Qinf
    real*8, intent(in) :: gamma, nx, ny
    real*8, dimension(:), intent(out) :: Qghost

    ! local vars
    real*8, dimension(4) :: uu_in, uu_inf, W_in, W_inf, W_star, Lambda
    real*8, dimension(4,4) :: L, R

    ! compute primitive vars at interior and inf
    call Q2q(Q, gamma, uu_in)
    call Q2q(Qinf, gamma, uu_inf)


    ! compute Lambda (eigen values), L (left eigen matrix)
    ! and R (right eigen matrix) freezed about current state
    ! at the interior point adjucent to boundary (+).
    !
    call comp_euler2d_L_R(uu_in, gamma, nx, ny, Lambda, L, R)
    W_in = matmul(L, uu_in)
    W_inf = matmul(L, uu_inf)

    !
    ! decide on the sign of eigen values
    ! how to blend characteristic vars (inlet or outlet)
    ! and compute the final characteristic vars in the
    ! STAR region
    !
    if ( Lambda(2) .eq. 0.0d0 ) then !tangential flow
       W_star = (/ W_inf(1), W_inf(2), W_in(3), W_inf(4) /)
    elseif     ( all(Lambda(2:4) > 0.0d0) .and. (Lambda(1) < 0.0d0 ) ) then
       !subsonic outlet
       W_star = (/ W_inf(1), W_in(2), W_in(3), W_in(4) /)
    elseif ( (Lambda(2) < 0.0d0) .and. (Lambda(3) > 0.0d0 ) ) then !subsonic inlet
       W_star = (/ W_inf(1), W_inf(2), W_in(3), W_inf(4) /)
    elseif ( all(Lambda > 0.0d0) ) then !supersonic outlet
       W_star = W_in
    elseif ( all(Lambda < 0.0d0) ) then !supersonic inlet
       W_star = W_inf
    else
       print *, 'could not find the type of inlet/outlet ' &
            , ' in specifying char-Bcs! stop'
       print *, 'Lambda = ', Lambda
       stop
    end if

    ! multiply by the inverse of left eigen matrix
    ! to obtain the primitive vars from char vars
    uu_in = matmul(R, W_star)

    ! check if density or pressure is negative then exit
    if ( (uu_in(1) <= 0.0d0) .or. (uu_in(4) <= 0.0d0) ) then
       print *, 'something is wrong in specifying char-BCs because ' &
            , ' it yeilds negative density or pressure! stop'
       stop
    end if

    ! convert primitive to conservative vars and finally
    ! save them to output values of Ughost
    ! subroutine u2U(rho, u, v, P, gamma, UU)
    call u2U(uu_in(1), uu_in(2), uu_in(3), uu_in(4), gamma, Qghost)

    ! done here
  end subroutine comp_char_bcs_euler2d

  ! converts conservative vars
  ! Q = [rho, rho *u, rho *v, e]
  ! to primitive vars
  ! qq = [rho, u, v, P]
  !
  subroutine Q2q(Q, gamma, qq)
    implicit none
    real*8, dimension(:), intent(in) :: Q
    real*8, intent(in) :: gamma
    real*8, dimension(:), intent(out) :: qq

    ! local vars
    real*8 :: rho, u, v, P, e

    ! calculating primitive variables
    rho = Q(1)
    u = Q(2) / Q(1)
    v = Q(3) / Q(1)
    e = Q(4)
    P = (gamma - 1.0d0) * e - 0.5d0 * (gamma - 1.0d0) *rho * ( u*u + v*v)

    ! pack and return the values
    qq(1) = rho
    qq(2) = u
    qq(3) = v
    qq(4) = P

    ! done here
  end subroutine Q2q

  ! compute the left and right eigen matrix
  ! of the two dimensional compressible euler
  ! equations in primitive form
  !
  ! input  :
  ! uu = [rho, u, v, P]^T vector of primitive vars
  ! output :
  ! L  = left eigen matrix
  ! R = right eigen matrix; R*L = L*R = I (when nx^2+ny^2=1)
  !
  ! Refer to book "I do like CFD" for more info.
  !
  subroutine comp_euler2d_L_R(uu, gamma, nx, ny, Lambda, L, R)
    implicit none
    real*8, dimension(:), intent(in) :: uu
    real*8, intent(in) :: gamma, nx, ny
    real*8, dimension(:), intent(out) :: Lambda
    real*8, dimension(:, :), intent(out) :: L, R

    ! local vars
    real*8 :: rho, u, v, P, c, lx, ly, qn

    ! bullet proofing ...
    if ( abs(sqrt(nx**2 + ny**2) - 1.0d0) > 1.0d-15 ) then
       print *, 'the vector (nx, ny) is not a unit vector! stop'
       stop
    end if

    ! assign primitive vars
    rho = uu(1)
    u = uu(2)
    v = uu(3)
    P = uu(4)

    ! comp c and init
    c = sqrt(gamma * P / rho)
    lx = -ny
    ly = nx
    qn = u*nx + v*ny

    ! fill Lambda
    Lambda(1) = qn - c
    Lambda(2) = qn
    Lambda(3) = qn + c
    Lambda(4) = qn

    ! fill L
    L(1,1) = 0.0d0
    L(1,2) = -nx * rho / (2.0d0 * c)
    L(1,3) = -ny * rho / (2.0d0 * c)
    L(1,4) = 1.0d0 / (2.0d0 * c**2)

    L(2,1) = 1.0d0
    L(2,2) = 0.0d0
    L(2,3) = 0.0d0
    L(2,4) = -1.0d0 / (c**2)

    L(3,1) = 0.0d0
    L(3,2) = nx * rho / (2.0d0 * c)
    L(3,3) = ny * rho / (2.0d0 * c)
    L(3,4) = 1.0d0 / (2.0d0 * c**2)

    L(4,1) = 0.0d0
    L(4,2) = rho * lx
    L(4,3) = rho * ly
    L(4,4) = 0.0d0

    ! fill R
    R(1,1) = 1.0d0
    R(1,2) = 1.0d0
    R(1,3) = 1.0d0
    R(1,4) = 0.0d0

    R(2,1) = -nx*c / rho
    R(2,2) = 0.0d0
    R(2,3) = nx*c / rho
    R(2,4) = lx / rho

    R(3,1) = -ny*c / rho
    R(3,2) = 0.0d0
    R(3,3) = ny*c / rho
    R(3,4) = ly / rho

    R(4,1) = c*c
    R(4,2) = 0.0d0
    R(4,3) = c*c
    R(4,4) = 0.0d0

    ! done here
  end subroutine comp_euler2d_L_R

  subroutine comp_weakly_proj_wall_U(Uin, gamma, nx, ny, Uout)
    implicit none
    real*8, dimension(:), intent(in) :: Uin
    real*8, intent(in) :: gamma, nx, ny
    real*8, dimension(:), intent(out) :: Uout

    ! local vars
    real*8 :: rho, P
    real*8, dimension(size(Uin)) :: qq
    real*8, dimension(2) :: u0, et, uu

    call Q2q(Uin, gamma, qq)

    rho = qq(1)
    u0(1) = qq(2)
    u0(2) = qq(3)
    P = qq(4)

    et(1) = ny
    et(2) = -nx

    call comp_projected_wall_vel(u0, et, uu)

    Uout(1) = rho
    Uout(2) = rho * uu(1)
    Uout(3) = rho * uu(2)
    Uout(4) = P / (gamma - 1.0d0) + 0.5d0 *rho * ( uu(1)*uu(1) + uu(2)*uu(2))

    ! done here
  end subroutine comp_weakly_proj_wall_U

  subroutine comp_projected_wall_vel(u0, et, u)
    implicit none
    real*8, dimension(:), intent(in) :: u0, et
    real*8, dimension(:), intent(out) :: u

    ! local vars
    real*8 :: u0_dot_et, norm_A1, norm_A2
    real*8, dimension(size(u0)) :: A1, A2


    u0_dot_et = sum(u0*et)

    A1 = -u0 + u0_dot_et * et
    A2 = -u0 - u0_dot_et * et

    norm_A1 = sqrt(sum(A1*A1))
    norm_A2 = sqrt(sum(A2*A2))

    if (norm_A1 < norm_A2 ) then
       u = u0 + A1
    else
       u = u0 + A2
    end if

    ! done here
  end subroutine comp_projected_wall_vel

end module euler2d_eqs

! ! tests euler2d eqs and fluxes etc.
! ! COMMENT after this when use this module
! ! in your project
! program tester
!   use euler2d_eqs
!   use orig_euler2d_flux_wrapper
!   implicit none

!   ! local vars
!   integer :: neqs, i, j, k, m, s, t
!   real*8 :: gamma, nx, ny, length, rho, Minf, alpha, u, v, P, c
!   logical, dimension(4) :: f_select
!   real*8, dimension(4) :: Q, fvl_p, fvl_m, fvl_p_ref, fvl_m_ref, err, fw, fw_ref, err2
!   real*8, dimension(4, 4) :: d_fvl_p, d_fvl_m, d_fvl_p_ref, d_fvl_m_ref, d_fw, d_fw_ref
!   real*8, parameter :: PI = 4.D0*DATAN(1.D0)

!   ! HARD reset
!   err = 0.0d0
!   err2 = 0.0d0

!   ! consts
!   neqs = 4

!   do i = 1, 6
!      do j = 1, 6
!         do k = 1, 6
!            do m = 1, 6
!               do s = 1, 6
!                  do t = 1, 6
!                     ! case
!                     gamma = 1.4d0 + dble(i - 1) / dble(5) * 2.0d0
!                     rho = 1.0d0 + dble(j - 1) / dble(5) * 2.0d0
!                     Minf = 0.01d0 + dble(k - 1) / dble(5) * 2.0d0
!                     alpha = -20.0d0 + dble(m - 1) / dble(5) * 40.0d0

!                     nx = -2.0d0 + dble(s - 1) / dble(5) * 4.0d0
!                     ny = -3.0d0 + dble(t - 1) / dble(5) * 6.0d0
!                     length = sqrt(nx**2 + ny**2)
!                     nx = nx / length
!                     ny = ny / length
!                     f_select = .true.
!                     alpha = alpha * PI / 180.0d0
!                     u = Minf * cos(alpha) * nx
!                     v = Minf * sin(alpha) * ny
!                     c = sqrt(u**2 + v**2) / Minf

!                     ! init. conservative vars
!                     Q(1) = rho
!                     Q(2) = rho * u
!                     Q(3) = rho * v
!                     P = c**2 * rho / gamma
!                     Q(4) = P / (gamma - 1.0d0) + 0.5d0 * rho * (u**2 + v**2)



!                     ! first call FORTRAN subroutine
!                     !
!                     ! f_select[1] = .false. -> not compute f+,  = .true. -> compute f+
!                     ! f_select[2] = .false. -> not compute f-,  = .true. -> compute f-
!                     ! f_select[3] = .false. -> not compute df+, = .true. -> compute df+
!                     ! f_select[4] = .false. -> not compute df-, = .true. -> compute df-
!                     call calc_van_leer(Q, neqs, gamma, nx, ny, f_select &
!                          , fvl_p, fvl_m, d_fvl_p, d_fvl_m)

!                     ! ! print computed vals
!                     ! print *, 'fvl_p = ', fvl_p
!                     ! print *, 'fvl_m = ', fvl_m
!                     ! print *, 'd_fvl_p = ', d_fvl_p
!                     ! print *, 'd_fvl_m = ', d_fvl_m

!                     ! then call oroginal C function and consider the results as reference values
!                     call calc_orig_calc_van_leer(Q, neqs, gamma, nx, ny, f_select &
!                          , fvl_p_ref, fvl_m_ref, d_fvl_p_ref, d_fvl_m_ref)

!                     ! ! print reference vals
!                     ! print *, 'fvl_p_ref = ', fvl_p_ref
!                     ! print *, 'fvl_m_ref = ', fvl_m_ref
!                     ! print *, 'd_fvl_p_ref = ', d_fvl_p_ref
!                     ! print *, 'd_fvl_m_ref = ', d_fvl_m_ref


!                     ! ! print errors
!                     ! print *, 'Linf(fvl_p) = ', maxval(abs(fvl_p - fvl_p_ref))
!                     ! print *, 'Linf(fvl_m) = ', maxval(abs(fvl_m - fvl_m_ref))
!                     ! print *, 'Linf(d_fvl_p) = ', maxval(abs(d_fvl_p - d_fvl_p_ref))
!                     ! print *, 'Linf(d_fvl_m) = ', maxval(abs(d_fvl_m - d_fvl_m_ref))

!                     err(1) = max(err(1), maxval(abs(fvl_p - fvl_p_ref)))
!                     err(2) = max(err(2), maxval(abs(fvl_m - fvl_m_ref)))
!                     err(3) = max(err(3), maxval(abs(d_fvl_p - d_fvl_p_ref)))
!                     err(4) = max(err(4), maxval(abs(d_fvl_m - d_fvl_m_ref)))

!                     ! compute wall flux
!                     call calc_wall_flux(Q, neqs, gamma, nx, ny, fw, d_fw)

!                     ! compute original wall flux from C-file
!                     call calc_orig_wall_flux(Q, neqs, gamma, nx, ny, fw_ref, d_fw_ref)

!                     ! calc errors
!                     err2(1) = max(err2(1), maxval(abs(fw - fw_ref)))
!                     err2(2) = max(err2(2), maxval(abs(d_fw - d_fw_ref)))


!                  end do
!               end do
!            end do
!         end do
!      end do
!   end do

!   ! print total analysis error
!   print *, 'error in the Van Leer flux and analytical Jacobians = ', err
!   print *, 'error in the wall flux and analytical Jacobians = ', err2

!   ! done here
! end program tester
