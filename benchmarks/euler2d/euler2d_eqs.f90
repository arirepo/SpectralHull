module euler2d_eqs
  implicit none

  private


  public :: calc_van_leer, calc_wall_flux, calc_pure_euler2d_flux

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
