module ps
  use ortho_poly
  use globals
  use grid_opt
  use fem_utils
  use spem_2d
  use interp

  implicit none

  private
  ! ! DEBUG
  ! real*8, dimension(10,10) :: uu

  type point
     integer :: tag
     real*8  :: x, y
  end type point

  type connector

     ! tag of this face
     integer :: tag    ! for relating to FEM stuff
     logical :: is_bn_con ! true if it is a boundary connector 
     type(point), pointer :: pt1, pt2

     ! fem coupling data structure
     real*8, dimension(:), pointer :: ps_x,ps_y
     real*8, dimension(:,:), pointer :: ps_u !in
     real*8, dimension(:,:), pointer :: ps_q !out

     real*8, dimension(:), allocatable :: fem_x,fem_y
     real*8, dimension(:,:), allocatable :: fem_u !out
     real*8, dimension(:,:), allocatable :: fem_q !in

     integer, dimension(:), allocatable :: thetri
     real*8, dimension(:), allocatable :: fem_r,fem_s

  end type connector

  type blk

     integer :: tag
     integer :: fem_cons_tag

     ! this HEX block has 6 faces so here they are:
     type(connector), dimension(:), allocatable :: con

     integer :: neqs
     real*8 :: delta_t
     ! block size in each direction
     integer, dimension(2) :: ni

     real*8, dimension(:,:), allocatable :: x, y

     ! primary and secondary vars u, q and qq is the second derivative
     ! i.e. qx(ieq, x, y)
     real*8, dimension(:,:,:), allocatable :: u, qx, qy, qqx, qqy
     real*8, dimension(:,:,:), allocatable :: rhs, Ax

     ! diff matrices
     real*8, dimension(:,:), allocatable  :: D1, D2

     integer, dimension(2) ::  x_range, y_range

     ! the initialization flag for each block
     ! this will contain 1221360 after the data 
     ! is initialized.
     integer :: blk_initialized_flag

  end type blk

  type gmres_struct
     ! num = number of Krylov sub. and nrst = num restarts
     integer :: num, nrst 
     real(rk) :: epsil ! gmres error tolerance
     ! the following is with size "nrst"
     real(rk), dimension(:), allocatable :: res
     ! the following are arrays needed for method 
     ! of the givens to solve upper Hassenburg h and
     ! update the solution
     real(rk), dimension(:), allocatable :: g, s, c, alpha

     ! the following contains the upper Hassenburg matrix 
     real(rk), dimension(:,:), allocatable :: h

     ! The diagonal preconditioner
     real(rk) :: N = 1._rk ! in future N should be like following
     !                   i j k eq ns eq

     ! The following contains arrays with the shape of space-time
     real(rk), dimension(:, :, :), pointer :: b => null() &
          , Ax => null() &
          , sol => null() 

     ! the following contains the Krylov subspace
     !                  eq i j num
     real(rk), dimension(:,:,:, :), allocatable :: v
     real(rk), dimension(:,:,:), allocatable :: yk

  end type gmres_struct


  public :: point, connector, blk
  public :: init_ps_blk
  public :: gmres_struct, init_gmres_struct, ps_blk_solve
  public :: write_blk_tecplot
  public :: find_overlap, print_con, sync_u_fem2ps
  public :: sync_u_ps2fem

contains 

  subroutine init_ps_blk(this_blk, ni, neqs, alpha, beta, delta_t)
    implicit none
    type(blk), intent(inout), target :: this_blk
    integer, dimension(:), intent(in) :: ni
    integer, intent(in) :: neqs
    real*8, intent(in) :: alpha, beta, delta_t

    ! locals
    integer :: i
    real*8, dimension(:,:), allocatable :: DD1, DD2 ! temp second der
    real*8, dimension(:), allocatable :: xx, yy

    if (.not. allocated(this_blk%con) ) then
       print *, 'error : the block connectors must be associated before initialization of the block'
       stop
    end if

    this_blk%neqs = neqs
    this_blk%delta_t = delta_t
    this_blk%ni = ni

    allocate(this_blk%x(ni(1), ni(2)), this_blk%y(ni(1), ni(2)))
    allocate(DD1(ni(1), ni(1)), DD2(ni(2), ni(2)))
    allocate(xx(ni(1)), yy(ni(2)))


    ! primary and secondary vars u, q and qq is the second derivative
    allocate(this_blk%u(neqs,ni(1),ni(2)), this_blk%qx(neqs,ni(1),ni(2)), &
         this_blk%qy(neqs,ni(1),ni(2)), this_blk%qqx(neqs,ni(1),ni(2)), &
         this_blk%qqy(neqs,ni(1),ni(2)))
    allocate(this_blk%rhs(neqs,ni(1),ni(2)), this_blk%Ax(neqs,ni(1),ni(2)))

    ! diff matrices
    allocate(this_blk%D1(ni(1),ni(1)), this_blk%D2(ni(2),ni(2)) )

    call comp_jacobi_diff_matrix(alpha, beta, ni(1)-1, this_blk%D1, DD1, xx)


    do i = 1, ni(2)
       this_blk%x(:,i) = xx
    end do

    call comp_jacobi_diff_matrix(alpha, beta, ni(2)-1, this_blk%D2, DD2, yy)

    do i = 1, ni(1)
       this_blk%y(i,:) = yy
    end do
    ! setting the range of the structured block
    ! bottom 
    if ( this_blk%con(1)%is_bn_con .eqv. .true. ) then 
       this_blk%y_range(1) = 2
    else
       this_blk%y_range(1) = 1
    end if
    ! top
    if ( this_blk%con(3)%is_bn_con .eqv. .true. ) then 
       this_blk%y_range(2) = ni(2) - 1
    else
       this_blk%y_range(2) = ni(2)
    end if
    ! left
    if ( this_blk%con(4)%is_bn_con .eqv. .true. ) then 
       this_blk%x_range(1) = 2
    else
       this_blk%x_range(1) = 1
    end if
    ! right
    if ( this_blk%con(2)%is_bn_con .eqv. .true. ) then 
       this_blk%x_range(2) = ni(1) - 1
    else
       this_blk%x_range(2) = ni(1)
    end if

    ! setting up block connectors props
    ! con 1
    this_blk%con(1)%ps_x => this_blk%x(:,1)
    this_blk%con(1)%ps_y => this_blk%y(:,1)
    this_blk%con(1)%ps_u => this_blk%u(:,:,1) 
    allocate(this_blk%con(1)%thetri(ni(1)))
    allocate(this_blk%con(1)%fem_r(ni(1)))
    allocate(this_blk%con(1)%fem_s(ni(1)))
    ! con 2
    this_blk%con(2)%ps_x => this_blk%x(ni(1),:)
    this_blk%con(2)%ps_y => this_blk%y(ni(1),:)
    this_blk%con(2)%ps_u => this_blk%u(:,ni(1),:) 
    allocate(this_blk%con(2)%thetri(ni(2)))
    allocate(this_blk%con(2)%fem_r(ni(2)))
    allocate(this_blk%con(2)%fem_s(ni(2)))
    ! con 3
    this_blk%con(3)%ps_x => this_blk%x(:,ni(2))
    this_blk%con(3)%ps_y => this_blk%y(:,ni(2))
    this_blk%con(3)%ps_u => this_blk%u(:,:,ni(2)) 
    allocate(this_blk%con(3)%thetri(ni(1)))
    allocate(this_blk%con(3)%fem_r(ni(1)))
    allocate(this_blk%con(3)%fem_s(ni(1)))
    ! con 4
    this_blk%con(4)%ps_x => this_blk%x(1,:)
    this_blk%con(4)%ps_y => this_blk%y(1,:)
    this_blk%con(4)%ps_u => this_blk%u(:,1,:) 
    allocate(this_blk%con(4)%thetri(ni(2)))
    allocate(this_blk%con(4)%fem_r(ni(2)))
    allocate(this_blk%con(4)%fem_s(ni(2)))

    ! seal it !!!
    this_blk%blk_initialized_flag = 1221360

    ! little clean-ups
    deallocate(DD1, DD2, xx, yy)

    ! done here
  end subroutine init_ps_blk

  ! computes the matrix-vercor multiplication in GMRE

  subroutine comp_Ax(tblk, yk)
    implicit none
    type(blk), intent(inout) :: tblk
    real*8, dimension(:,:,:), intent(in) :: yk

    integer :: i, j, ieq

    ! compute the derivative in x
    do j = 1, tblk%ni(2)
       do ieq = 1, tblk%neqs
          tblk%qx(ieq, :, j) = matmul(tblk%D1, yk(ieq,:,j))
       end do
    end do

    ! compute the derivative in y
    do i = 1, tblk%ni(1)
       do ieq = 1, tblk%neqs
          tblk%qy(ieq, i, :) = matmul(tblk%D2, yk(ieq,i,:))
       end do
    end do

    ! compute the 2nd derivative in x
    do j = 1, tblk%ni(2)
       do ieq = 1, tblk%neqs
          tblk%qqx(ieq, :, j) = matmul(tblk%D1, tblk%qx(ieq,:,j))
       end do
    end do

    ! compute the 2nd derivative in y
    do i = 1, tblk%ni(1)
       do ieq = 1, tblk%neqs
          tblk%qqy(ieq, i, :) = matmul(tblk%D2, tblk%qy(ieq,i,:))
       end do
    end do

    ! ! DEBUG 
    ! uu = (tblk%x + tblk%y)*cos(tblk%x * tblk%y)
    ! ! print *, shape(uu), shape(tblk%qqx)
    ! print *, 'hey you div-exact ', maxval(abs(tblk%qx(1,:,:) + tblk%qy(1,:,:) - uu))

    ! uu = -(tblk%x*tblk%x + tblk%y*tblk%y)*sin(tblk%x * tblk%y)
    ! ! print *, shape(uu), shape(tblk%qqx)
    ! print *, 'hey you laplacian - exact ', maxval(abs(tblk%qqx(1,:,:) + tblk%qqy(1,:,:) - uu))
    ! stop

    ! finally compute matrix - vector product Ax
    tblk%Ax = yk - tblk%delta_t * (tblk%qqx + tblk%qqy)
  
    ! done here
  end subroutine comp_Ax

  ! The following initializes the gmres data
  ! structure just after the GMRES routine is called
  ! this acts like an adaptor in a way that it 
  ! connects the GMRES internal data structure to the
  ! arrays avalilable to the given block "blk" 
  ! also it allocates appropriate temp buffers.
  subroutine init_gmres_struct(tblk, gm, num, nrst, epsil)
    implicit none
    type(blk), target :: tblk
    type(gmres_struct), intent(inout) :: gm
    integer , intent(in) :: num, nrst
    real(rk), intent(in) :: epsil

    ! local vars
    ! range variables
    integer :: neqs
    ! i-range, j-range
    integer, dimension(:), pointer :: ir, jr
    integer :: i1, i2, j1, j2

    if(verbose_flag) print *, 'init. GMRES allocated' &
         , ' buffer and temp arrays ...'

    neqs = tblk%neqs ! numb. of equations in that block

    ! setting up the range arrays
    ir => tblk%x_range
    jr => tblk%y_range

    i1 = ir(1); i2 =ir(2)
    j1 = jr(1); j2 =jr(2)
 

    ! initializing the scalar values
    gm%num = num
    gm%nrst = nrst
    gm%epsil = epsil

    ! initializing / allocating array types
    allocate(gm%res(nrst)) 
    gm%res = 0._rk
    allocate(gm%g(num+1), gm%s(num), gm%c(num), gm%alpha(num))
    gm%g = 0._rk; gm%s = 0._rk; gm%c = 0._rk; gm%alpha = 0._rk

    ! allocate the matrix types
    allocate(gm%h(num+1,num))
    gm%h = 0._rk

    ! init N here when diagonal preconditioner is needed.!

    ! init pointer to the interior node section 
    ! of the the rhs vector tblk%rhs
    gm%b => tblk%rhs(1:neqs, i1:i2, j1:j2)
    gm%Ax => tblk%Ax(1:neqs, i1:i2, j1:j2)
    gm%sol => tblk%u(1:neqs, i1:i2, j1:j2)

    !
    allocate(gm%v(neqs, (i2- i1 + 1), (j2-j1+1) , num+1))
    allocate(gm%yk(size(tblk%u,1), size(tblk%u,2), size(tblk%u,3)))

    ! print *, 'shape(gm%Ax) = ', shape(gm%Ax)
    ! print *, 'lbound(gm%Ax) = ', lbound(gm%Ax,3)
    ! print *, 'shape(gm%v)   = ', shape(gm%v)
    ! stop

    if(verbose_flag) print *, 'DONE init. GMRES allocated' &
         , ' buffer and temp arrays ...'

    ! done here
  end subroutine init_gmres_struct


  ! pseudo-spectral matrix-free GMREs solver
  ! tblk = the current active block.
  ! num = number of Krylov subspace
  ! nrst = number of restarts
  ! epsil = GMRES tolerance 
  subroutine ps_gmres(tblk, gm)
    implicit none
    type(blk), target :: tblk
    type(gmres_struct) :: gm


    ! locals
    integer :: k, j, i
    integer :: finish_flag, jrst
    real(rk) :: norm_tmp, delta, gamma, resi

    ! i-range, j-range
    integer, dimension(:), pointer :: ir, jr
    integer :: i1, i2, j1, j2

    ! setting up the range arrays
    ir => tblk%x_range
    jr => tblk%y_range

    i1 = ir(1); i2 =ir(2)
    j1 = jr(1); j2 =jr(2)


    if(verbose_flag) print *, 'starting DPI-GMRES ...' 

    ! resetting the flag
    finish_flag = 0
    resi = 0._rk
    gm%yk = 0.0d0

    ! do j number of restarts
    do jrst = 1, gm%nrst    
       ! reset ------------
       gm%v = 0._rk; gm%h = 0._rk
       gm%c = 0._rk; gm%s = 0._rk
       gm%g = 0._rk
       gm%alpha = 0._rk
       ! ------------------
       ! compute Ax, result goes to => tblk%SAx
       ! where gm%Ax is pointing to it!
       call comp_Ax(tblk, tblk%u)
       ! r0 = gm%Ax to preserve the memory
       gm%Ax = gm%b - gm%Ax
       ! print *, 'dddd????%%%$$$', maxval(abs(gm%Ax(:,:,:,1:1,:)))
       ! !stop
       ! ! stop
       ! !print *,' --=-===>', maxval(abs( gm%b ))
       ! !stop
       norm_tmp = norm2(gm%Ax)
       gm%v(:,:,:,1) = gm%Ax/norm_tmp
       gm%g(1) = norm_tmp
       ! Arnoldi iterations    
       do k = 1,  gm%num
          gm%g(k+1) = 0._rk
          ! the following does yk = v(:,k) / N
          ! result is stored in  gm%sol =section> tblk%q 
          ! instead of yk
          gm%yk(:,i1:i2,j1:j2) =  gm%v(:,:,:,k) /  gm%N
          ! the following does uk = matmul(A, yk)
          ! the result is in gm%Ax =section> tblk%SAx 
          ! instead of uk
          call comp_Ax(tblk, gm%yk)
          do j = 1, k
             gm%h(j,k) = sum(gm%v(:,:,:,j) * gm%Ax)
             !grahm-schmidt orthogonalization
             gm%Ax = gm%Ax - gm%h(j,k) * gm%v(:,:,:,j) 
          end do
          gm%h(k+1,k) = norm2(gm%Ax)        
          gm%v(:,:,:,k+1) = gm%Ax/gm%h(k+1,k)
          ! applying givens
          do j = 1, (k-1)
             delta = gm%h(j,k)
             gm%h(j,k) = gm%c(j)*delta + gm%s(j)*gm%h(j+1,k)
             gm%h(j+1,k) = -gm%s(j)*delta + gm%c(j) * gm%h(j+1,k)
          end do
          gamma = sqrt(gm%h(k,k)**(2._rk) + gm%h(k+1,k)**(2._rk))
          gm%c(k) = gm%h(k,k) / gamma
          gm%s(k) = gm%h(k+1,k) / gamma
          gm%h(k,k) = gamma
          gm%h(k+1,k) = 0._rk
          delta = gm%g(k)
          gm%g(k) = gm%c(k) * delta + gm%s(k) * gm%g(k+1)
          gm%g(k+1) = -gm%s(k) * delta + gm%c(k) * gm%g(k+1)
          resi = abs(gm%g(k+1))
          print *,'>>> resi = ',  resi
          if( resi <= gm%epsil) then
             finish_flag = 1
             goto 100
          end if
          !if(verbose_flag) print *,'Arnoldi itr = <',k,'> completed!'
          !print *, resi
       end do
       k = k - 1
       ! solving backward for alpha
100    gm%alpha(k) = gm%g(k) / gm%h(k,k)
       do i = k-1, 1, -1
          gm%alpha(i) = gm%g(i) / gm%h(i,i)
          do j = i+1, k      
             gm%alpha(i) = gm%alpha(i) - gm%h(i,j)/gm%h(i,i) * gm%alpha(j)     
          end do
       end do
       ! compute the final directional vector using
       ! the combination of search vectors
       ! note: we used gm%Ax as zk symbol in original GMRES. 
       gm%Ax = 0._rk
       do j = 1, k
          gm%Ax = gm%Ax + gm%alpha(j)*gm%v(:,:,:,j)
       end do
       ! updating solution
       gm%sol = gm%b + gm%Ax / gm%N
       gm%res(jrst) = resi !records final residual
       if(finish_flag == 1) then 
          goto 200        
       end if
    end do
    ! subroutine final clean-ups
200 if(verbose_flag) print *, 'Done DPI-GMRES! Deallocating' &
         ,' DPI-GMRES temp arrays ...'
    print *, ' GMRES residual : ', gm%res
    !stop

    if(verbose_flag) print *, 'DONE deallocating!'

    return !done <here>

  end subroutine ps_gmres
  
  ! takes second norm 
  function norm2(x)
    implicit none
    real(rk), dimension(:,:,:), intent(in) :: x
    real(rk) :: norm2

    norm2 = sqrt(sum(x*x))

    ! done
  end function norm2

  ! writes blk to tecplot
  subroutine write_blk_tecplot(tblk, plt)
    implicit none
    type(blk), target, intent(in) :: tblk
    type(plt_spec), intent(in) :: plt

    !locals
    integer :: i, j, zz, imax, jmax
    integer :: zzmax, dimmax
    character ::  tipchar
    real(rk), dimension(:,:), pointer :: x, y
    real(rk), dimension(:,:,:), pointer :: u

    ! pointer association
    x => tblk%x
    y => tblk%y
    u => tblk%u


    imax = size(x,1); !dim 1
    jmax = size(x,2); !dim 2
    zzmax = size(u,1); !number of equations
    dimmax = 2; !number of dimension


    if (plt%append_flag) then
       ! append
       open(unit = 24, file = plt%filename, position = 'append')
    else !or just open a new file
       open(unit = 24, file = plt%filename)
    end if

    write(24,*) trim('TITLE = '), trim(plt%title)
100 FORMAT(A) !formats the character string after this
300 FORMAT(1D25.16) ! formats floating point numbers
    !write the VARIABLE block
    write(24, 100, advance = 'no') trim('VARIABLES = ')
    do i = 1, (zzmax+dimmax)
       !check if we have enough legends
       tipchar = trim(plt%legends(i))
       if ( iachar(tipchar) .ne. 34 ) then !34 is ASCII for "
          print *, 'WARNING, the 3d-plot legend #', i,'is unreadable or is blank!'
       end if
       write(24, 100, advance = 'no') plt%legends(i)
    end do

    write(24,*)
    !write the ZONE block
    write(24,*) trim('ZONE, I= ') , imax, ', J=', jmax, trim(', F=POINT')
    !write the data 
    do j = 1, jmax
       do i = 1, imax
          write(24, 300, advance = 'no') x(i,j), y(i,j)  
          do zz = 1, zzmax
             write(24, 300, advance = 'no') u(zz,i,j)
          end do
          write(24,*)
       end do
    end do
    !close the output file
    close(24)
    !done

  end subroutine write_blk_tecplot
  !----------------------------------------------
  
  ! the following solves one pseudo-spectral block
  ! with given dirichlet BCs on boundaries
  ! the results are stored in that block data
  ! structure
  subroutine ps_blk_solve(tblk, gm, itrmax , tol, echo)
    implicit none
    type(blk)  :: tblk
    type(gmres_struct) :: gm
    integer, intent(in) :: itrmax 
    real*8, intent(in) :: tol ! du <= tol is termination
    logical, intent(in) :: echo ! verbose mode is true

    ! local vars
    integer :: itr
    real*8 :: du

    if (echo) print *, 'ps-solving blk # ', tblk%tag

    do itr = 1, itrmax

       du = sqrt(sum((tblk%rhs - tblk%u)**2.0d0))
       if ( echo ) print *, 'itr = ', itr, 'du = ', du

       if ( du <= tol ) exit

       tblk%rhs = tblk%u
       call ps_gmres(tblk, gm)

    end do

    if (echo) print *, 'the ps-solve for blk', tblk%tag, ' converged to du = ', du, ' in ', itr, 'iterations'

    ! done here
  end subroutine ps_blk_solve

  ! find the overlap region of the given block
  ! with the given finite element grid and
  ! store it in the connector structure
  subroutine find_overlap(tblk, grd, tol)
    implicit none
    type(blk), target :: tblk
    type(grid), intent(in) :: grd
    real*8, intent(in) :: tol

    ! local vars
    integer :: l, i, seed, thetri
    type(connector), pointer :: tcon

    seed = 1

    do l = 1, 4 ! all connectors on this blk

       tcon => tblk%con(l)

       if (.not. tcon%is_bn_con) cycle

       do i = 1, size(tcon%ps_x)

          call search_tri(seed, tcon%ps_x(i), tcon%ps_y(i), grd, thetri)
          tcon%thetri(i) = thetri
          seed = thetri
          call xy2rs(tcon%ps_x(i), tcon%ps_y(i), thetri, grd, tol &
               , 10, tol, tcon%fem_r(i), tcon%fem_s(i))
       end do

    end do

    ! done here
  end subroutine find_overlap

  subroutine print_con(tcon)
    implicit none
    type(connector), intent(in) :: tcon

    print *, 'properties of connector : ', tcon%tag
    print *, ' '
    print *, ' '
    print *, 'is boundary connector ? ', tcon%is_bn_con
    print *, 'pt1 = ', tcon%pt1%tag, 'pt2 = ', tcon%pt2%tag
    print *, 'ps_x = ', tcon%ps_x
    print *, 'ps_y = ', tcon%ps_y
    print *, 'ps_u = ', tcon%ps_u
    print *, 'thetri = ', tcon%thetri
    print *, 'fem_r = ', tcon%fem_r 
    print *, 'fem_s = ', tcon%fem_s 

    ! done here
  end subroutine print_con

  subroutine sync_u_fem2ps(tblk, fem)
    implicit none
    type(blk), target :: tblk
    type(fem_struct), intent(in) :: fem

    ! local vars
    integer :: l, i
    type(connector), pointer :: tcon


    do l = 1, 4 ! all connectors on this blk

       tcon => tblk%con(l)

       if (.not. tcon%is_bn_con) cycle

       do i = 1, size(tcon%ps_x)
          call fem_interp_u(tcon%fem_r(i), tcon%fem_s(i), tcon%thetri(i), fem%grd, fem%u, tcon%ps_u(:,i))
       end do

    end do

    ! done here
  end subroutine sync_u_fem2ps
  
  subroutine sync_u_ps2fem(tblk, fem)
    implicit none
    type(blk), intent(in) :: tblk
    type(fem_struct), target :: fem

    ! local vars
    integer :: i, j, k, pt
    integer :: neqs, id
    real*8 :: x,y
    integer :: ieq
    type(grid), pointer :: grd
    real*8, dimension(:,:,:,:), pointer :: KK
    real*8, dimension(:,:), pointer :: rhs

    grd => fem%grd
    neqs = fem%neqs
    KK => fem%KK
    rhs => fem%rhs

    ! loop over ALL boundary edges
    do i =1, grd%nbedgeg 

       ! get id of that edge
       id = grd%ibedgeBC(i)
       if ( id .ne. tblk%fem_cons_tag) cycle

       ! we are good to go now ...
       do ieq = 1, neqs

          do k = 1, size(grd%ibedge,2) ! loop over nodes on that edge
             pt = grd%ibedge(i,k) ! global number of that node
             x = grd%x(pt)
             y = grd%y(pt)
             KK(ieq,:,pt,:) = 0.0d0 ! freeze that row 
             do j = ieq, ieq
                KK(j,j,pt,pt) = 1.0d0 ! ones on the diagonal
             end do
             call lag2d(tblk%x(:,1), tblk%y(1,:), tblk%u(ieq,:,:), x, y, rhs(ieq,pt))
          end do ! k - loop

       end do ! loop over equations

    end do ! done all boundary edges

    ! done here
  end subroutine sync_u_ps2fem

end module ps
