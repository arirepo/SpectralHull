module mfree_st
  use grid_opt
  use element_opt
  use mfree
  use ortho_poly
  use precond_opt
  use normal_rand
  implicit none

  private

  type, extends(mf_struct) :: mf_struct_st

     integer :: nS
     real*8 :: alpha, beta, t0, t1, dt
     real*8, dimension(:), allocatable :: x, t, unk, u0, du0, Mu0, Mdu0
     real*8, dimension(:,:), allocatable :: S, SS, zeta, tmp

  contains

    procedure :: init_st => init_mf_struct_st
    procedure :: Ax_st => full_A_x_st
    procedure :: rhs => comp_rhs_st
    procedure :: gmres => gmres_st
    procedure :: window => window_st
    procedure :: stbicg => stab_bi_cg_chen
    procedure :: idrs => idrs_advanced

  end type mf_struct_st

  public :: mf_struct_st

contains

  subroutine init_mf_struct_st(this, grd, elems, tags, values &
       , nS, alpha, beta, t0, t1, u0, du0)
    implicit none
    class(mf_struct_st), intent(inout) :: this
    type(grid), intent(in), target :: grd
    type(element), dimension(:), intent(in), target :: elems
    integer, dimension(:), intent(in) :: tags
    real*8, dimension(:), intent(in) :: values
    integer, intent(in) :: nS
    real*8, intent(in) :: alpha, beta, t0, t1
    real*8, dimension(:), intent(in) :: u0, du0

    ! local vars
    integer :: i, neqs, nnodesg

    ! use the parrent init procedure to init the base
    call this%mf_struct%init(grd, elems, tags, values)

    ! now, init the member of this child class
    neqs = this%elems(1)%neqs
    nnodesg = this%grd%nnodesg
    this%nS = nS
    this%alpha = alpha; this%beta = beta; this%t0 = t0; this%t1 = t1
    this%dt = (t1 - t0)
    allocate(this%x(nS), this%t(nS), this%unk(size(this%free)))
    allocate(this%u0(neqs* nnodesg), this%du0(neqs* nnodesg) &
         , this%Mu0(neqs* nnodesg), this%Mdu0(neqs* nnodesg))
    allocate(this%S(nS, nS), this%SS(nS, nS), this%zeta(nS, nS))
    allocate(this%tmp(neqs* nnodesg, nS))
    ! safe initi.
    this%x = 0.0d0; this%t = 0.0d0; this%unk = 0.0d0
    this%u0 = 0.0d0; this%du0 = 0.0d0; this%Mu0 = 0.0d0; this%Mdu0 = 0.0d0 
    this%S = 0.0d0; this%SS = 0.0d0; this%zeta = 0.0d0
    this%tmp = 0.0d0

    ! fill out the integration operator
    call comp_jacobi_integ_matrix(alpha, beta, nS, this%S, this%SS, this%x)

    ! transform -1 < this%x <= 1 to physical space t0 < this%t <= t1
    this%t = t0 + 0.5d0 * (this%x + 1.0d0) * (t1 - t0)

    ! transforming integration operator this%S
    this%S = 0.5d0 * (t1 - t0) * this%S

    ! computing zeta
    this%SS = 0.0d0
    this%zeta = 0.0d0
    do i = 1, nS
       this%SS(i, i) = this%t(i)
    end do

    ! final formulae
    this%zeta = matmul(this%SS, this%S) - matmul(this%S, this%SS)  

    ! initializing initial conditions
    this%u0 = u0
    this%du0 = du0

    ! computing Mu0 and Mdu0
    ! setting the pointers to the mass matrix
    do i = 1, size(this%elems)
       this%elems(i)%A => this%elems(i)%M 
    end do
    ! perform matrix-vector product
    this%Mu0 = 0.0d0
    call this%mf_struct%Ax_gen(this%u0 , this%Mu0 )
    this%Mdu0 = 0.0d0
    call this%mf_struct%Ax_gen(this%du0, this%Mdu0)

    ! done here
  end subroutine init_mf_struct_st


  subroutine full_A_x_st(this, x, xout)
    implicit none
    class(mf_struct_st), intent(inout) :: this
    real*8, dimension(:,:), intent(in) :: x
    real*8, dimension(:,:), intent(out) :: xout

    ! local vars
    integer :: i, j

    ! HARD reset
    xout = 0.0d0

    ! setting the pointers to the stiffness matrix
    do i = 1, size(this%elems)
       this%elems(i)%A => this%elems(i)%K 
    end do

    ! perform space-time K*x
    do i = 1, this%nS ! for all time splices
       call this%mf_struct%Ax_gen(x(:, i), xout(:, i))
    end do

    ! performing zeta*K*x
    this%tmp = xout
    do i = 1, this%nS
       this%unk = 0.0d0
       do j = 1, this%nS
          this%unk = this%unk + this%zeta(i, j) * xout(this%free, j)
       end do
       this%tmp(this%free, i) =  this%unk
    end do

    ! saving a copy so this%tmp can be used further
    ! so after this xout = zeta *K * x 
    xout = this%tmp

    ! setting the pointers to the mass matrix
    do i = 1, size(this%elems)
       this%elems(i)%A => this%elems(i)%M 
    end do

    ! perform space-time M*x
    this%tmp = 0.0d0
    do i = 1, this%nS ! for all time splices
       call this%mf_struct%Ax_gen(x(:, i), this%tmp(:, i))
    end do

    ! accumulating final result A*x = zeta * K * x + m*x 
    do i = 1, this%nS
       xout(this%free, i) = xout(this%free, i) + this%tmp(this%free, i)
    end do

    ! done here
  end subroutine full_A_x_st

  subroutine comp_rhs_st(this, rhs)
    implicit none
    class(mf_struct_st), intent(inout) :: this
    real*8, dimension(:,:), intent(out) :: rhs

    ! local vars
    integer :: i

    !
    do i = 1, this%nS
       rhs(:, i) = this%Mu0
       rhs(this%free, i) = rhs(this%free, i) + this%Mdu0(this%free) * (this%t(i) - this%t0)
    end do

    ! done here
  end subroutine comp_rhs_st

  subroutine window_st(this, x)
    implicit none
    class(mf_struct_st), intent(inout) :: this
    real*8, dimension(:, :), intent(in) :: x

    ! local vars
    integer :: i, j

    ! window the last time slice solution as initial condition
    ! for the next space-time packet
    this%u0 = x(:, this%nS)

    ! setting the pointers to the mass matrix
    do i = 1, size(this%elems)
       this%elems(i)%A => this%elems(i)%M 
    end do

    ! update this%Mu0
    this%Mu0 = 0.0d0
    call this%mf_struct%Ax_gen(this%u0, this%Mu0)

    ! proceeding to compute Mdu0 ...

    ! setting the pointers to the stiffness matrix
    do i = 1, size(this%elems)
       this%elems(i)%A => this%elems(i)%K 
    end do

    ! perform space-time K*x
    ! HARD reset
    this%tmp = 0.0d0
    do i = 1, this%nS ! for all time splices
       call this%mf_struct%Ax_gen(x(:, i), this%tmp(:, i))
    end do

    ! performing S*K*x at t = t(nS)
    i = this%nS
    this%unk = 0.0d0
    do j = 1, this%nS
       this%unk = this%unk - this%S(i, j) * this%tmp(this%free, j)
    end do

    ! update
    this%Mdu0(this%free) = this%Mdu0(this%free) + this%unk

    this%t0 = this%t1
    this%t1 = this%t0 + this%dt 
    this%t = this%t0 + 0.5d0 * (this%x + 1.0d0) * (this%t1 - this%t0)

    ! computing zeta
    this%SS = 0.0d0
    this%zeta = 0.0d0
    do i = 1, this%nS
       this%SS(i, i) = this%t(i)
    end do

    ! final formulae
    this%zeta = matmul(this%SS, this%S) - matmul(this%S, this%SS)  

    ! done here
  end subroutine window_st

  subroutine gmres_st(mf, b, N, epsil, num, nrst, sol, res)
    implicit none
    ! inputs
    class(mf_struct_st), intent(inout) :: mf
    real*8, dimension(:, :), intent(in) :: b
    type(precond) :: N
    real*8, intent(in) :: epsil
    integer, intent(in) :: num, nrst

    !outputs
    real*8, dimension(:, :), intent(inout) :: sol
    real*8, dimension(:), allocatable :: res

    ! locals
    integer :: k, j, i
    integer :: finish_flag, jrst
    real*8, dimension(:,:), allocatable :: v, h, r0, yk
    real*8, dimension(:), allocatable :: g, uk, s, c, zk, alpha, tmp
    real*8 :: norm_tmp, delta, gamma, resi
    integer :: nrows, nvars, len_res
    real*8, dimension(size(mf%free) * mf%nS) :: buff

    ! resetting the flag
    finish_flag = 0
    resi = 0.0d0
    nrows = size(b, 1) * size(b, 2) 
    nvars = size(mf%free) * mf%nS
    len_res = 0

    ! we don't use dynamic mem allocation in the gmres inner-outer loops 
    ! to increase spped. this might not be memory efficient though.
    allocate(v(nvars,num+1), r0(size(b,1), size(b, 2)), g(num+1), s(num), c(num))
    allocate(yk(size(b,1), size(b, 2)), uk(nvars), zk(nvars) )
    allocate(h(num+1,num), alpha(num))

    ! do j number of restarts
    do jrst = 1, nrst    
       ! reset ------------
       v = 0.0d0; h = 0.0d0
       c = 0.0d0; s = 0.0d0
       zk = 0.0d0; g = 0.0d0
       alpha = 0.0d0
       yk = 0.0d0
       ! ------------------
       ! ***********************
       call mf%Ax_st(sol, r0)
       r0 = b - r0
       ! r0 = b - matmul(A,sol)
       ! ***********************
       v(:,1) = reshape(r0(mf%free, :), (/nvars /))
       norm_tmp = norm(v(:,1))
       v(:,1) = v(:,1)/norm_tmp
       g(1) = norm_tmp
       ! Arnoldi iterations    
       do k = 1, num
          g(k+1) = 0.0d0
          ! ***********************
          yk = 0.0d0
          call N%solve(v(:,k), buff)
          ! buff = v(:,k) / N
          yk(mf%free, :) = reshape(buff, (/ size(mf%free), mf%nS /) )
          ! yk(mf%free, :) = v(:,k) / N
          call mf%Ax_st(yk, r0)
          uk = reshape(r0(mf%free, :), (/ nvars /) )
          ! yk = v(:,k) / N
          ! uk = matmul(A, yk)
          ! ***********************
          do j = 1, k
             h(j,k) = dot_product(v(:,j),uk)
             !grahm-schmidt orthogonalization
             uk = uk - h(j,k) * v(:,j) 
          end do
          h(k+1,k) = norm(uk)        
          v(:,k+1) = uk/h(k+1,k)
          ! applying givens
          do j = 1, (k-1)
             delta = h(j,k)
             h(j,k) = c(j)*delta + s(j)*h(j+1,k)
             h(j+1,k) = -s(j)*delta + c(j) * h(j+1,k)
          end do
          gamma = sqrt(h(k,k)**(2.0d0) + h(k+1,k)**(2.0d0))
          c(k) = h(k,k) / gamma
          s(k) = h(k+1,k) / gamma
          h(k,k) = gamma
          h(k+1,k) = 0.0d0
          delta = g(k)
          g(k) = c(k) * delta + s(k) * g(k+1)
          g(k+1) = -s(k) * delta + c(k) * g(k+1)
          resi = abs(g(k+1))
          ! add to residual history
          len_res = len_res + 1
          allocate(tmp(len_res))
          if(len_res > 1) tmp(1:(len_res-1)) = res(1:(len_res-1)) 
          tmp(len_res) = resi
          call move_alloc(tmp, res)
          !
          if( resi <= epsil) then
             finish_flag = 1
             goto 100
          end if

       end do
       k = k - 1
       ! solving backward for alpha
100    alpha(k) = g(k) / h(k,k)
       do i = k-1, 1, -1
          alpha(i) = g(i) / h(i,i)
          do j = i+1, k      
             alpha(i) = alpha(i) - h(i,j)/h(i,i) * alpha(j)     
          end do
       end do
       ! compute the final directional vector using
       ! the combination of search vectors 
       do j = 1, k
          zk = zk + alpha(j)*v(:,j)
       end do
       ! updating solution
       ! ***********************
       uk = 0.0d0
       call N%solve(zk, uk)
       zk = uk 
       ! zk = zk/N
       r0 = 0.0d0
       r0(mf%free, :) = reshape(zk, (/ size(mf%free), mf%nS /) )
       sol = sol + r0     
       ! sol = sol + zk/N    
       ! ***********************
       if(finish_flag == 1) then 
          goto 200        
       end if
    end do
    ! subroutine final clean-ups
200 deallocate(v, r0, g, s, c)
    deallocate(yk, uk, zk, h, alpha)

    return !done <here>

  end subroutine gmres_st

  ! 
  subroutine stab_bi_cg_chen(mf, b, N, epsil, itrs, sol, res)
    implicit none
    ! inputs
    class(mf_struct_st), intent(inout) :: mf
    real*8, dimension(:, :), intent(in) :: b
    type(precond) :: N
    real*8, intent(in) :: epsil
    integer, intent(in) :: itrs

    !outputs
    real*8, dimension(:, :), intent(inout) :: sol
    real*8, dimension(:), allocatable :: res

    ! locals
    integer :: j, nvars, len_res
    real*8 :: alpha, omega, resi, beta
    real*8, dimension(:), allocatable :: rj, rjp1, rbar0, p
    real*8, dimension(:), allocatable :: ptilde, y, s, stilde
    real*8, dimension(:), allocatable :: y2, dx, tmp
    real*8, dimension(:,:), allocatable :: x0, x, ytmp

    ! init
    nvars = size(mf%free) * mf%nS
    len_res = 0

    ! 
    allocate(rj(nvars), rjp1(nvars), rbar0(nvars), p(nvars))
    allocate(ptilde(nvars), y(nvars), s(nvars), stilde(nvars))
    allocate(y2(nvars), dx(nvars))
    allocate(x0(size(b, 1), size(b, 2)), x(size(b, 1), size(b, 2)))
    allocate(ytmp(size(b, 1), size(b, 2)))

    ! start ...
    call mf%Ax_st(sol, x0)
    x0 = b - x0
    rj = reshape(x0(mf%free, :), (/nvars /))
    !rj = b - A * x0
    rbar0 = rj
    p = rj
    x = sol
    if(allocated(res) ) deallocate(res)
    ytmp = 0.0d0

    do j = 1, itrs

       call N%solve(p, ptilde)
       !
       ytmp(mf%free, :) = reshape(ptilde, (/ size(mf%free), mf%nS /) )
       call mf%Ax_st(ytmp, x0)
       y = reshape(x0(mf%free, :), (/nvars /) )
       ! y = A*ptilde
       alpha = sum(rj * rbar0) / sum(y * rbar0)
       s = rj - alpha * y
       call N%solve(s, stilde)
       !
       ytmp(mf%free, :) = reshape(stilde, (/ size(mf%free), mf%nS /) )
       call mf%Ax_st(ytmp, x0)
       y2 = reshape(x0(mf%free, :), (/ nvars /))
       ! y2 = A*stilde
       omega = sum(y2 * s)/ sum(y2 * y2)
       dx = alpha * ptilde + omega * stilde
       x(mf%free, :) = x(mf%free, :) + reshape(dx, (/ size(mf%free), mf%nS /) )
       resi = sqrt(sum(dx*dx))
       ! add to residual history
       len_res = len_res + 1
       allocate(tmp(len_res))
       if(len_res > 1) tmp(1:(len_res-1)) = res(1:(len_res-1)) 
       tmp(len_res) = resi
       call move_alloc(tmp, res)
       !
       if ( resi <= epsil) then
          goto 500
       end if

       rjp1 = s - omega * y2
       beta = sum(rjp1 * rbar0) / sum(rj * rbar0) * (alpha / omega)
       p = rjp1 + beta * (p - omega * y)
       rj = rjp1
    end do

    print *, 'Warning : Preconditioned Stabilized Bi Conjugate' &
         , ' Gradient Method did not converge ' &
         , 'with the specified tolerance = ', epsil, '!!!!'

500 sol(mf%free, :) = x(mf%free, :)

    ! clean ups
    deallocate(rj, rjp1, rbar0, p)
    deallocate(ptilde, y, s, stilde)
    deallocate(y2, dx)
    deallocate(x0, x)
    deallocate(ytmp)


    ! done here
  end subroutine stab_bi_cg_chen

  !
  !               The Main IDRS Solver
  !
  ! This is the improved version of original SIAM paper
  ! which has smoothing and preconditiong remedies.
  !
  subroutine idrs_advanced(mf, b, s, tol, maxit, NN, x0 &
       , smoothing, x, flag, resvec, iter, relres)
    implicit none
    class(mf_struct_st), intent(inout) :: mf
    real*8, dimension(:, :), intent(in) :: b
    integer, intent(in) :: s
    real*8, intent(in) :: tol
    integer, intent(in) :: maxit
    type(precond) :: NN
    real*8, dimension(:, :), intent(in) :: x0
    integer, intent(in) :: smoothing

    ! outputs
    real*8, dimension(:, :), intent(out) :: x
    integer, intent(out) :: flag
    real*8, dimension(:), allocatable :: resvec
    integer, intent(out) :: iter
    real*8, intent(out) :: relres

    ! local vars
    integer :: k, i, nvars
    integer :: m, n
    real*8 :: angle, mp, normb, tolb, normr, trueres, om, alpha, beta, gamma
    real*8, dimension(size(mf%free) * mf%nS, s) :: P, Q, G, U
    real*8, dimension(s, s) :: MM
    real*8, dimension(size(mf%free) * mf%nS, 1) :: r, x_s, r_s, v, tmp, t
    real*8, dimension(s, 1) :: f, c 
    integer ( kind = 4 ) :: seed
    real*8, dimension(size(b ,1), size(b ,2)) :: x0tmp, uu

    ! HARD init
    x = 0.0d0
    flag = 0
    iter = 0
    relres = 0.0d0
    x0tmp = 0.0d0
    uu = 0.0d0

    nvars = size(mf%free) * mf%nS
    m = nvars
    n = m

    angle = 0.7d0
    seed = 123456789
    call r8mat_normal_01 ( n, s, seed, P )
    call orth(P, Q)
    P = Q

    ! % Number close to machine precision:
    mp = 1.0d3 * epsilon(1.0d0)

    relres = 0.0d0
    relres = relres / relres ! gets NaN by a trick

    ! Compute initial residual:
    x = x0
    normb = norm(reshape( b(mf%free, :), (/ nvars /) ) )
    tolb = tol * normb   ! Relative tolerance

    ! ------------------------------------------------- 
    call mf%Ax_st(x0, x0tmp)
    x0tmp = b - x0tmp
    r(:,1) = reshape(x0tmp(mf%free, :), (/nvars /))
    ! r(:,1) = b - matmul(A, x)
    ! ------------------------------------------------- 

    if (smoothing .ne. 0) then
       x_s(:, 1) = reshape( x0(mf%free, :), (/ nvars /) )
       r_s = r
    end if

    normr = norm(r(:,1))
    if ( allocated(resvec) ) deallocate(resvec) ! just for the 1st time
    call pushit(resvec, normr)
    trueres = 0.0d0

    if (normr <= tolb) then  !Initial guess is a good enough solution
       iter = 0                 
       flag = 0
       relres = normr/normb
       return
    end if
    ! %
    G = 0.0d0; U = 0.0d0; MM = 0.0d0
    do i = 1, s
       MM(i, i) = 1.0d0
    end do
    om = 1.0d0
    !

    ! % Main iteration loop, build G-spaces:
    iter = 0
    do ! while ( normr > tolb && iter < maxit )  

       ! % New righ-hand size for small system:
       f = transpose(matmul(transpose(r), P))
       do k = 1, s 
          !
          ! % Solve small system and make v orthogonal to P:
          call lin_solve(MM(k:s,k:s), f(k:s,1), c(k:s,:))
          v = r - matmul(G(:,k:s), c(k:s,:)) 
          ! %
          ! % Preconditioning:
          call NN%solve(v(:,1), tmp(:,1))
          ! call lin_solve(M1, v(:,1), tmp)
          v = tmp
          ! call lin_solve(M2, v(:,1), tmp)
          ! v = tmp

          ! % Compute new U(:,k) and G(:,k), G(:,k) is in space G_j
          U(:,k) = matmul(U(:,k:s), c(k:s,1)) + om * v(:, 1)
          ! ------------------------------------------------- 
          uu = 0.0d0
          uu(mf%free, :) = reshape( U(:,k), (/ size(mf%free), mf%nS /) )
          call mf%Ax_st(uu, x0tmp)
          G(:,k) = reshape(x0tmp(mf%free, :), (/nvars /))
          ! G(:,k) = matmul(A, U(:,k))
          ! ------------------------------------------------- 

          ! %
          ! % Bi-Orthogonalise the new basis vectors: 
          do i = 1, (k-1)
             alpha =  sum(P(:,i) * G(:,k)) / MM(i,i)
             G(:,k) = G(:,k) - alpha*G(:,i)
             U(:,k) = U(:,k) - alpha*U(:,i)
          end do
          ! %
          ! % New column of MM = P'*G  (first k-1 entries are zero)
          MM(k:s,k:k) = transpose(matmul(transpose(G(:,k:k)), P(:,k:s)))
          if ( MM(k,k) .eq. 0.0d0 ) then
             flag = 3
             return
          end if
          ! %
          ! %  Make r orthogonal to q_i, i = 1..k 
          beta = f(k,1)/MM(k,k)
          r = r - beta*G(:,k:k)
          !
          uu = 0.0d0
          uu(mf%free, :) = reshape(U(:,k), (/ size(mf%free), mf%nS /) )
          x = x + beta * uu
          ! x = x + beta*U(:,k)
          !
          normr = norm(r(:, 1))

          ! %
          ! %  Smoothing:
          if ( smoothing .ne. 0) then
             t = r_s - r
             gamma =  sum(t*r_s)/sum(t*t)
             r_s = r_s - gamma*t
             ! 
             tmp(:, 1) = reshape(x(mf%free, :), (/ nvars /) )
             x_s = x_s - gamma*(x_s - tmp)
             !
             normr = norm(r_s(:,1))
          end if
          call pushit(resvec, normr)
          iter = iter + 1
          if ( (normr < tolb) .or. (iter .eq. maxit ) ) exit
          ! %
          ! % New f = P'*r (first k  components are zero)
          if ( k < s ) then 
             f((k+1):s, 1)   = f((k+1):s, 1) - beta*MM((k+1):s,k)
          end if
       end do
       ! %
       if ( (normr < tolb) .or. (iter .eq. maxit) ) then
          exit
       end if
       ! %
       ! % Now we have sufficient vectors in G_j to compute residual in G_j+1
       ! % Note: r is already perpendicular to P so v = r
       ! %
       ! % Preconditioning:
       v = r
       call NN%solve(v(:,1), tmp(:,1))
       v = tmp
       ! call lin_solve(M1, v(:,1), tmp)
       ! v = tmp
       ! call lin_solve(M2, v(:,1), tmp)
       ! v = tmp

       ! %
       ! % Matrix-vector multiplication:
       ! ------------------------------------------------- 
       uu = 0.0d0
       uu(mf%free, :) = reshape( v(:,1), (/ size(mf%free), mf%nS /) )
       call mf%Ax_st(uu, x0tmp)
       t(:, 1) = reshape(x0tmp(mf%free, :), (/nvars /))
       ! t = matmul(A, v) 
       ! ------------------------------------------------- 


       ! % Computation of a new omega
       om = omega( t(:, 1), r(:, 1), angle )
       if ( om .eq. 0.0d0 ) then
          flag = 3
          return
       end if
       ! %
       r = r - om*t
       ! 
       uu = 0.0d0
       uu(mf%free, :) = reshape(v(:,1), (/ size(mf%free), mf%nS /) )
       x = x + om * uu
       ! x = x + om*v(:,1)
       normr = norm(r(:,1))
       ! %
       ! %     Smoothing:
       if ( smoothing .ne. 0 ) then
          t    = r_s - r
          gamma = sum(t*r_s)/sum(t*t)
          r_s = r_s - gamma*t
          ! 
          tmp(:, 1) = reshape(x(mf%free, :), (/ nvars /) )
          x_s = x_s - gamma*(x_s - tmp)
          ! x_s(:,1) = x_s(:,1) - gamma*(x_s(:,1) - x)
          normr = norm(r_s(:,1))
       end if
       ! %
       call pushit(resvec, normr)
       iter = iter + 1

       if ( (normr > tolb) .and. (iter < maxit) ) then
          cycle ! continue while loop
       else
          exit  ! exit while loop
       end if
    end do ! while

    if ( smoothing .ne. 0) then
       x(mf%free, :) = reshape(x_s(:,1), (/ size(mf%free), mf%nS /) )
    end if

    call mf%Ax_st(x, x0tmp)
    x0tmp = b - x0tmp

    tmp(:, 1) = reshape(x0tmp(mf%free, :), (/ nvars /) )

    relres = norm(tmp(:, 1)) / normb

    if ( relres < tol ) then 
       flag = 0
    else if ( iter .eq. maxit ) then
       flag = 1
    else
       flag = 2
    end if

    ! done here
  end subroutine idrs_advanced

  function omega( t, s, angle ) result(om)
    implicit none
    real*8, dimension(:), intent(in) :: t, s
    real*8, intent(in) :: angle
    real*8 :: om

    ! locals
    real*8 :: ns, nt, ts, rho

    ns = norm(s)
    nt = norm(t)
    ts = sum(t*s)
    rho = abs(ts/(nt*ns))
    om = ts/(nt*nt)
    if ( rho < angle ) then
       om = om*angle/rho
    end if

    ! done here
  end function omega

  ! generates orthonormal basis for range of
  ! a not neccessarily square matrix <A>
  subroutine orth(A, Q)
    implicit none
    real*8, dimension(:, :), intent(in)  :: A
    real*8, dimension(:, :), intent(out) :: Q

    ! local vars
    character*1 :: jobu, jobvt
    real*8, dimension(:, :), allocatable :: u, vt

    integer :: m, n, lda, lwork, INFO
    real*8, dimension(:), allocatable :: s, work

    ! write to temp; do not change original <A>
    Q = A

    ! init
    m = size(A, 1)
    n = size(A, 2)
    jobu = 'O'
    jobvt = 'N'
    lda = m
    allocate(s(min(m,n)))
    allocate(u(1, 1), vt(1, 1)) ! really not referenced!!!
    lwork = max(1,3*min(m,n)+max(m,n),5*min(m,n))
    allocate(work(max(1,lwork)))

    ! call LAPACK double precision generic SVD routine
    ! NOTE : svd is Overwritten in [Q]
    call dgesvd(jobu, jobvt, m, n, Q, lda, s, u, 1, vt, 1, work, lwork, INFO)

    if (INFO .ne. 0) then
       print *, 'something is wrong in SVD in orth(...)! stop'
       stop
    end if

    ! clean ups
    deallocate(u, vt)
    deallocate(s, work)

    ! done here
  end subroutine orth

  ! adds the double precision val to the
  ! end of the array <a> 
  subroutine pushit(a, val)
    implicit none
    real*8, dimension(:), allocatable :: a
    real*8, intent(in) :: val

    ! local vals
    integer :: tsize, n
    real*8, dimension(:), allocatable :: tmp

    if ( .not. allocated(a) ) then 
       tsize = 0
    else 
       tsize = size(a)
    end if

    n = tsize + 1
    allocate(tmp(n))
    if (tsize >= 1) tmp(1:tsize) = a(1:tsize)
    tmp(n) = val

    call move_alloc(tmp, a)

    ! done here
  end subroutine pushit

  ! performs LAPACK general LU solve DGESV
  subroutine lin_solve(A, b, x)
    implicit none
    real*8, dimension(:, :), intent(in) :: A
    real*8, dimension(:), intent(in) :: b
    real*8, dimension(:,:), intent(out) :: x

    ! local vars
    integer :: N, NRHS, LDA, LDB, INFO
    real*8, dimension(size(A,1), size(A,1)) :: Atmp
    integer, dimension(size(A,1)) :: IPIV

    ! init
    N = size(A,1)
    NRHS = 1
    Atmp = A ! take copy; DONT modify the original one!
    LDA = N
    x(:,1) = b
    LDB = N
    INFO = 0
    call DGESV( N, NRHS, Atmp, LDA, IPIV, x, LDB, INFO )

    if ( INFO .ne. 0) then
       print *, 'something is wrong in lin_solve! stop'
       stop
    end if

    ! done here
  end subroutine lin_solve

end module mfree_st
