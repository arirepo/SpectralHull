module approx_fekete
  use quadri_elem
  implicit none
  private

  type fekete
     private

     ! private components
     integer, dimension(:), allocatable :: n
     integer :: nw, d, d_big, d_orig
     real*8, dimension(:), allocatable :: coords
     real*8, dimension(:), allocatable :: rhs
     real*8, dimension(:, :), allocatable :: x0
     integer ::  magnify, s
     character(len = 128) :: spacing
     logical :: initialized = .false., echo = .false.

     ! public components
     real*8, dimension(:), allocatable, public :: xy
     procedure(gen_xy), public, pointer, nopass  :: comp_xy, moment
     real*8, dimension(:), allocatable, public :: w
     real*8, dimension(:, :), allocatable, public :: fin_coords, xy0
     character(len = 128) :: name

   contains

     ! procedures
     procedure, public :: init => initialize
     procedure, public :: set_pt => set_point_coord
     procedure, private :: find => find_approx_fekete_pts
     procedure, public :: clean => clean_fekete_obj
     procedure, public :: export => write_pts_to_file

  end type fekete

  abstract interface
     subroutine gen_xy(this)
       import
       type(fekete), intent(inout) :: this
     end subroutine gen_xy
  end interface

  type fekete_table
     private

     type(fekete), dimension(:), allocatable :: feketes

   contains

     procedure, public :: lookup => lookup_in_fekete_table 

  end type fekete_table

  public :: fekete, fekete_table

contains

  ! find a given fekete element of any type and dimension,
  ! i.e. based on requested "p" and "name",
  ! in the current lookup table. If not found, then
  ! adds the newly requested element to the end of the lookup table
  ! and also, in the same time, returns that element  
  subroutine lookup_in_fekete_table(this, d, name, fekete_out &
       , magnify, s, spacing, echo)
    implicit none
    class(fekete_table), intent(inout), target :: this
    integer, intent(in) :: d
    character(len = *), intent(in) :: name
    type(fekete), intent(out) :: fekete_out ! return value 
    integer, intent(in), optional :: magnify, s
    character(len = *), intent(in), optional :: spacing
    logical, intent(in), optional :: echo

    ! local vars
    integer :: i, nf
    class(fekete), pointer :: tfekete => null()
    type(fekete), dimension(:), allocatable, target :: tmp

    ! if the table has not been created yet, then make it
    ! this happens in the first time call only 
    if( .not. allocated(this%feketes) ) then
       allocate(this%feketes(1))
       tfekete => this%feketes(1)
       call tfekete%init(d, name, magnify, s, spacing, echo)
       ! return the added element
       fekete_out = this%feketes(1)
       return
    end if

    ! we already have the table, then look up the fekete element 
    ! to see if it is there. if not, proceed to add it then  
    do i = 1, size(this%feketes)

       tfekete => this%feketes(i)

       ! compare 
       if ( (tfekete%d_orig .eq. d) .and. (tfekete%name .eq. name) ) then !found
          fekete_out = this%feketes(i)
          return
       end if

    end do

    ! OK, the requested element is not in the table, then add it!
    nf = size(this%feketes) !number of current fekete elems in the table
    nf = nf + 1 ! add one more
    allocate(tmp(nf))
    tmp(1:(nf-1)) = this%feketes !copy previous elems 
    tfekete => tmp(nf) !last one in the new one to be added as follows
    call tfekete%init(d, name, magnify, s, spacing, echo)

    ! move tmp to the main table
    call move_alloc(tmp, this%feketes)

    ! prepare the return value 
    fekete_out = this%feketes(nf)

    ! done here
  end subroutine lookup_in_fekete_table

  subroutine initialize(this, d, name, magnify, s, spacing, echo)
    implicit none
    class(fekete), intent(inout) :: this
    integer, intent(in) :: d
    character(len = *), intent(in) :: name
    integer, intent(in), optional :: magnify, s
    character(len = *), intent(in), optional :: spacing
    logical, intent(in), optional :: echo

    ! local vars
    integer :: dim
    real*8, dimension(:), allocatable :: x, y, z

    ! bullet proofing
    if (this%initialized) then
       print *, 'the fekete object is already been initialized! stop'
       stop
    end if

    ! setting up the element name and original degree "d_orig"
    this%name = name
    this%d_orig = d

    ! setting default values
    if ( .not. present(magnify)) then
       this%magnify = 4
    else
       this%magnify = magnify
    end if
    if ( .not. present(s)) then 
       this%s = 15
    else
       this%s = s
    end if
    if ( .not. present(spacing)) then 
       this%spacing = 'equal_space'
    else
       this%spacing = spacing
    end if
    if ( .not. present(echo)) then 
       this%echo = .false.
    else
       this%echo = echo
    end if

    ! choose the type of element
    select case (name)

    case ('quadri') ! 2D quadrilateral 

       ! init.
       dim = 2
       allocate(this%n(dim))
       this%n = d+1
       this%nw = this%n(1) * this%n(2)
       allocate(this%coords(dim))
       allocate(this%w(this%nw), this%rhs(this%nw), this%xy(this%nw))
       allocate(this%fin_coords(dim, this%nw), this%xy0(this%nw, this%nw))

       ! generate master element point distribution of the dense magnified grid
       call gen_master_quadri((this%magnify * this%n(1)) &
            , (this%magnify * this%n(2)), -1.0d0, 1.0d0, -1.0d0, 1.0d0 &
            , this%spacing, x, y)

       ! store it in the current fekete object 
       allocate(this%x0(dim, size(x)))
       this%x0(1, :) = x; this%x0(2, :) = y

       ! clean temp mem buffers
       deallocate(x, y)

       ! setting up pointers to appropriate xy eval and moment subroutines
       this%comp_xy => xy_at_pt_quadri
       this%moment => moment_quadri

    case ('1d') ! one-dimensional  

       ! init.
       dim = 1
       allocate(this%n(dim))
       this%n = d+1
       this%nw = this%n(1)
       allocate(this%coords(dim))
       allocate(this%w(this%nw), this%rhs(this%nw), this%xy(this%nw))
       allocate(this%fin_coords(dim, this%nw), this%xy0(this%nw, this%nw))

       ! generate master element point distribution of the dense magnified grid
       allocate(this%x0(1, (this%magnify * this%n(1))))
       call pt_dist_1d(n = (this%magnify * this%n(1)) &
            , xmin = -1.0d0, xmax = 1.0d0 &
            , method = this%spacing, x = this%x0(1, :))

       ! setting up pointers to appropriate xy eval and moment subroutines
       this%comp_xy => xy_at_pt_1d
       this%moment => moment_1d

    case ('tri') ! triangle 

       ! init.
       dim = 2
       this%d = d
       this%nw = (this%d + 1) * (this%d + 2) / 2
       this%d_big = magnify * this%d
       allocate(this%coords(dim))
       allocate(this%w(this%nw), this%rhs(this%nw), this%xy(this%nw))
       allocate(this%fin_coords(dim, this%nw), this%xy0(this%nw, this%nw))

       ! generate master element point distribution of the dense magnified grid
       call coord_tri(this%d_big, x, y)

       ! store it in the current fekete object 
       allocate(this%x0(dim, size(x)))
       this%x0(1, :) = x; this%x0(2, :) = y

       ! clean temp mem buffers
       deallocate(x, y)

       ! setting up pointers to appropriate xy eval and moment subroutines
       this%comp_xy => xy_at_pt_tri
       this%moment => moment_tri

    case ('tet') ! tet 

       ! init.
       dim = 3
       this%d = d
       this%nw = (this%d + 1) * (this%d + 2) * (this%d + 3) / 6
       this%d_big = magnify * this%d
       allocate(this%coords(dim))
       allocate(this%w(this%nw), this%rhs(this%nw), this%xy(this%nw))
       allocate(this%fin_coords(dim, this%nw), this%xy0(this%nw, this%nw))

       ! generate master element point distribution of the dense magnified grid
       call coord_tet(this%d_big, x, y, z)

       ! store it in the current fekete object 
       allocate(this%x0(dim, size(x)))
       this%x0(1, :) = x; this%x0(2, :) = y; this%x0(3, :) = z

       ! clean temp mem buffers
       deallocate(x, y, z)

       ! setting up pointers to appropriate xy eval and moment subroutines
       this%comp_xy => xy_at_pt_tet
       this%moment => moment_tet

    case ('hex') ! 3D hex 

       ! init.
       dim = 3
       allocate(this%n(dim))
       this%n = d+1
       this%nw = this%n(1) * this%n(2) * this%n(3)
       allocate(this%coords(dim))
       allocate(this%w(this%nw), this%rhs(this%nw), this%xy(this%nw))
       allocate(this%fin_coords(dim, this%nw), this%xy0(this%nw, this%nw))

       ! generate master element point distribution of the dense magnified grid
       call coord_hex((this%magnify * this%n), x, y, z)

       ! store it in the current fekete object 
       allocate(this%x0(dim, size(x)))
       this%x0(1, :) = x; this%x0(2, :) = y; this%x0(3, :) = z

       ! clean temp mem buffers
       deallocate(x, y, z)

       ! setting up pointers to appropriate xy eval and moment subroutines
       this%comp_xy => xy_at_pt_hex
       this%moment => moment_hex

    case ('prism') ! prism 

       ! init.
       dim = 3
       this%d = d
       this%nw = (this%d + 1) * (this%d + 2) / 2 * (this%d + 1)
       this%d_big = magnify * this%d
       allocate(this%coords(dim))
       allocate(this%w(this%nw), this%rhs(this%nw), this%xy(this%nw))
       allocate(this%fin_coords(dim, this%nw), this%xy0(this%nw, this%nw))

       ! generate master element point distribution of the dense magnified grid
       call coord_prism(this%d_big, x, y, z)

       ! store it in the current fekete object 
       allocate(this%x0(dim, size(x)))
       this%x0(1, :) = x; this%x0(2, :) = y; this%x0(3, :) = z

       ! clean temp mem buffers
       deallocate(x, y, z)

       ! setting up pointers to appropriate xy eval and moment subroutines
       this%comp_xy => xy_at_pt_prism
       this%moment => moment_prism

    case default
       print *, 'unknown name/type of element to generate' &
            , ' approx fekete points! stop'
       stop

    end select

    ! set the flag that init. is done!
    this%initialized = .true.

    ! calc approx fekete pts by computing
    ! quadrature points
    call this%find(this%x0, this%s)
 
    ! done here
  end subroutine initialize
  
  subroutine clean_fekete_obj(this)
    implicit none
    class(fekete), intent(inout) :: this

    ! bullet proofing
    if (.not. this%initialized) then
       print *, 'the fekete object must be initalized before cleaning! stop'
       stop
    end if

    ! private components
    if(allocated(this%n)) deallocate(this%n)
    if(allocated(this%coords)) deallocate(this%coords)
    if(allocated(this%w)) deallocate(this%w)
    if(allocated(this%rhs)) deallocate(this%rhs)
    if(allocated(this%fin_coords)) deallocate(this%fin_coords)
    if(allocated(this%xy0)) deallocate(this%xy0)
    if(allocated(this%x0)) deallocate(this%x0)
    if(allocated(this%xy)) deallocate(this%xy)

    ! final step
    this%initialized = .false.

    ! done here
  end subroutine clean_fekete_obj

  subroutine set_point_coord(this, coords)
    implicit none
    class(fekete), intent(inout) :: this
    real*8, dimension(:), intent(in) :: coords

    this%coords = coords

    ! done here
  end subroutine set_point_coord

  subroutine xy_at_pt_quadri(this)
    implicit none
    type(fekete), intent(inout) :: this

    ! local val
    integer :: i, j, jj

    jj = 1

    do i = 1, this%n(1)
       do j = 1, this%n(2)
          this%xy(jj) = (this%coords(1)**(i-1)) * (this%coords(2)**(j-1))
          jj = jj + 1
       end do
    end do

    !done here
  end subroutine xy_at_pt_quadri

  subroutine xy_at_pt_1d(this)
    implicit none
    type(fekete), intent(inout) :: this

    ! local val
    integer :: i

    do i = 1, this%n(1)
       this%xy(i) = this%coords(1)**(i-1)
    end do

    !done here
  end subroutine xy_at_pt_1d

  subroutine xy_at_pt_tri(this)
    implicit none
    type(fekete), intent(inout) :: this

    ! local val
    integer :: i, j, jj

    jj = 1

    do i = 0, this%d
       do j = 0, i
          this%xy(jj) = (this%coords(1)**(i-j)) * (this%coords(2)**j)
          jj = jj + 1
       end do
    end do

    !done here
  end subroutine xy_at_pt_tri

  subroutine xy_at_pt_tet(this)
    implicit none
    type(fekete), intent(inout) :: this

    ! local val
    integer :: i, j, k, jj

    jj = 1

    do i = 0, this%d
       do j = 0, i
          do k = 0, j
             this%xy(jj) = (this%coords(1)**(i-j)) * (this%coords(2)**(j-k)) &
                  * (this%coords(3)**k)
             jj = jj + 1
          end do
       end do
    end do

    !done here
  end subroutine xy_at_pt_tet

  subroutine xy_at_pt_hex(this)
    implicit none
    type(fekete), intent(inout) :: this

    ! local val
    integer :: i, j, k, jj

    jj = 1

    do i = 1, this%n(1)
       do j = 1, this%n(2)
          do k = 1, this%n(3)
             this%xy(jj) = (this%coords(1)**(i-1)) * (this%coords(2)**(j-1)) &
                  * (this%coords(3)**(k-1))
             jj = jj + 1
          end do
       end do
    end do

    !done here
  end subroutine xy_at_pt_hex

  subroutine xy_at_pt_prism(this)
    implicit none
    type(fekete), intent(inout) :: this

    ! local val
    integer :: i, j, k, jj

    jj = 1

    do i = 0, this%d
       do j = 0, i
          do k = 1, (this%d + 1)
             this%xy(jj) = (this%coords(1)**(i-j)) * (this%coords(2)**j) * (this%coords(3)**(k-1))
             jj = jj + 1
          end do
       end do
    end do

    !done here
  end subroutine xy_at_pt_prism

  subroutine moment_quadri(this)
    implicit none
    type(fekete), intent(inout) :: this

    ! local vars
    integer :: i, j, K 

    K = 1
    do i = 1, this%n(1)
       do j = 1, this%n(2)
          this%rhs(K) = 1.0d0/dble(i*j) * ( 1.0d0 - (-1.0d0)**(i)) &
               * (1.0d0 - (-1.0d0)**(j))
          K = K + 1
       end do
    end do

    ! done here
  end subroutine moment_quadri

  subroutine moment_1d(this)
    implicit none
    type(fekete), intent(inout) :: this

    ! local vars
    integer :: i 

    do i = 1, this%n(1)
       this%rhs(i) = 1.0d0/dble(i) * ( 1.0d0 - (-1.0d0)**i)
    end do

    ! done here
  end subroutine moment_1d

  subroutine moment_tri(this)
    implicit none
    type(fekete), intent(inout) :: this

    ! local vars
    integer :: i, j, K 

    K = 1
    do i = 0, this%d
       do j = 0, i
          this%rhs(K) = gamma(dble(j+2)) * gamma(dble(i - j + 1)) &
               / ( dble(j+1) * gamma(dble(i+3)) )
          K = K + 1
       end do
    end do

    ! done here
  end subroutine moment_tri

  subroutine moment_tet(this)
    implicit none
    type(fekete), intent(inout) :: this

    ! local vars
    integer :: i, j, k, jj 

    jj = 1
    do i = 0, this%d
       do j = 0, i
          do k = 0, j
             this%rhs(jj) = gamma(dble(k+1)) * gamma(dble(-k + j + 1)) &
                  * gamma(dble(i - j + 1)) / gamma(dble(i+4))
             jj = jj + 1
          end do
       end do
    end do

    ! done here
  end subroutine moment_tet

  subroutine moment_hex(this)
    implicit none
    type(fekete), intent(inout) :: this

    ! local vars
    integer :: i, j, k, jj 

    jj = 1
    do i = 1, this%n(1)
       do j = 1, this%n(2)
          do k = 1, this%n(3)
             this%rhs(jj) = 1.0d0/dble(i*j*k) * ( 1.0d0 - (-1.0d0)**(i)) &
                  * (1.0d0 - (-1.0d0)**(j)) * (1.0d0 - (-1.0d0)**(k))
             jj = jj + 1
          end do
       end do
    end do

    ! done here
  end subroutine moment_hex

  subroutine moment_prism(this)
    implicit none
    type(fekete), intent(inout) :: this

    ! local vars
    integer :: i, j, k, jj 

    jj = 1
    do i = 0, this%d
       do j = 0, i
          do k = 1, (this%d + 1)
             this%rhs(jj) = gamma(dble(j+2)) * gamma(dble(i - j + 1)) &
                  / ( dble(j+1) * gamma(dble(i+3)) ) * 1.0d0 / dble(k)
             jj = jj + 1
          end do
       end do
    end do

    ! done here
  end subroutine moment_prism

  ! Finds multidimensional approximate Fekete points
  ! by finding "nw" number of Gauss-Lobatto quadrature points among
  ! initially provided dense scattered points x0(dim, M) 
  ! The inappropriate points which have zero weights 
  ! are eliminated during iterative QR-factorization 
  ! and the final points are good approximate of Fekete points
  ! and have small growth rate of Lebesgue constant.
  ! the quadrature weights "w" is also reported for
  ! the case if these points are used for numerical integration.
  ! the xy0 matrix is a well conditioned Vandermonde matrix 
  ! which is very suitable for numerical interpolation (if needed). 
  subroutine find_approx_fekete_pts(this, x0, s)
    implicit none
    class(fekete), intent(inout) :: this
    real*8, dimension(:, :), intent(in) :: x0
    integer, intent(in) :: s

    ! local vars
    integer :: M, N, k
    real*8, dimension(:, :), allocatable :: Vk, Pk, Rk, mu, w, VkT
    integer, dimension(:), allocatable :: indx

    ! bullet proofing
    if ( .not. this%initialized ) then
       print *, 'the fekete object needs to be initialized first before find! stop'
       stop
    end if

    ! allocate
    allocate(Vk(size(x0, 2), this%nw), Pk(this%nw, this%nw), Rk(this%nw, this%nw))
    allocate(VkT(this%nw, size(x0, 2)))
    allocate(mu(this%nw, 1), w(size(x0, 2), 1), indx(size(x0, 2))) 

    ! init
    M = size(x0, 2)
    N = this%nw
    Vk = 0.0d0
    w = 0.0d0

    ! computing Vandermonde matrix
    if (this%echo) print *, 'computing Vandermonde matrix ...'
    do k = 1, M
       ! freez the point
       call this%set_pt(x0(:, k))

       ! comp. xy polynomials
       call this%comp_xy(this)

       ! store the complete row
       Vk(k, :) = this%xy
    end do

    ! computing moments and store in this%rhs
    if (this%echo) print *, 'computing moments ...'
    call this%moment(this)

    ! init QR-iteration matrices
    Pk = 0.0d0
    do k = 1, N
       Pk(k, k) = 1.0d0
    end do

    ! perform QR
    do k = 1, s

       if (this%echo) print *, 'QR-iteration # ', k , ' of total ', s
       if (this%echo) print *, '===================================>'
       if (this%echo) print *, 'computing QR ...'
       call wrapper_qr(Vk, Rk)

       if (this%echo) print *, 'solving Vk ...'
       call wrapper_solve(Rk, Vk, Vk, 'T', .true.)

       if (this%echo) print *, 'solving Pk ...'
       call wrapper_solve(Rk, Pk, Pk, 'T', .true.)

    end do

    ! computes wights over dense grid
    Rk = transpose(Pk)
    mu(:,1) = matmul(Rk, this%rhs)

    ! solving underdetermined system to find weights on coarse grid
    VkT = transpose(Vk)
    call solve_underdet(VkT, mu, w, indx)

    ! extract the significant weights (nonzero)
    do k = 1, N 
          this%w(k) = w(indx(k), 1)
          this%fin_coords(:, k) = x0(:, indx(k))
          this%xy0(k, :) = Vk(indx(k), :)
    end do

    ! clean ups
    deallocate(Vk, Pk, Rk, VkT)
    deallocate(mu, w, indx) 

    ! done here
  end subroutine find_approx_fekete_pts

  ! wrapper around LAPACK QR factorization
  subroutine wrapper_qr(A, R, indx, Q, dopiv)
    implicit none
    real*8, dimension(:, :), intent(in) :: A
    real*8, dimension(:, :), intent(out) :: R
    integer, dimension(size(A, 2)), optional :: indx
    real*8, dimension(size(A, 1), size(A, 1)), optional :: Q
    logical, intent(in), optional :: dopiv

    ! local vars
    integer :: M, N, LDA, LWORK, INFO
    real*8, dimension(:, :), allocatable :: Atmp, H, v
    integer, dimension(:), allocatable :: JPVT
    real*8, dimension(:), allocatable :: TAU, WORK
    integer :: i, j

    ! allocate
    allocate(Atmp(size(A, 1), size(A, 2)), H(size(A, 1), size(A, 1)))
    allocate(v(size(A, 1), 1))

    !HARD Reset
    R = 0.0d0

    ! take a copy first!
    Atmp = A

    ! init
    M = size(A, 1)
    N = size(A, 2)
    LDA = max(1, M)
    LWORK = 2 * (3*N+1)
    allocate(JPVT(N), TAU(min(M, N)), WORK(max(1, LWORK)))
    TAU = 0.0d0; WORK = 0.0d0
    JPVT = (/ ( i , i = 1, N) /)
    if ( present(dopiv) ) then
       if (dopiv) JPVT = 0
    end if
    ! subroutine dgeqp3(integer M,
    !   integer N,
    !   double precision, dimension( lda, * ) A,
    !   integer LDA,
    !   integer, dimension( * ) JPVT,
    !   double precision, dimension( * ) TAU,
    !   double precision, dimension( * ) WORK,
    !   integer LWORK,
    !   integer INFO 
    !   )

    call DGEQP3(M, N, Atmp, LDA, JPVT, TAU, WORK, LWORK, INFO)

    if ( INFO .ne. 0 ) then
       print *, 'something is wrong in QR factorization! INFO = ', INFO, ' stop'
       stop
    end if

    ! extract R
    do i = 1, min(M, N)
       do j = i, N
          R(i, j) = Atmp(i, j)
       end do
    end do

    ! extract pivot
    if (present(indx)) indx = JPVT

    ! Extract Q
    if ( present(Q) ) then

       Q = 0.0d0
       do j = 1, M
          Q(j, j) = 1.0d0
       end do

       do i = 1, min(M, N)

          if( i > 1 ) v(1:(i-1),1) = 0.0d0
          v(i,1) = 1.0d0
          v((i+1):M,1) = Atmp((i+1):M, i)
          H = -1.0d0 * TAU(i) * matmul(v, transpose(v))
          do j = 1, M
             H(j, j) = H(j, j) + 1.0d0
          end do

          Q = matmul(Q, H)
       end do

    end if

    ! clean ups
    deallocate(Atmp, H, v)
    deallocate(JPVT, TAU, WORK)

    ! done here
  end subroutine wrapper_qr

  subroutine wrapper_solve(A, X, B, TRANS, trans_rhs)
    implicit none
    real*8, dimension(:, :), intent(in) :: A
    real*8, dimension(:, :), intent(inout) :: X, B
    character*1, intent(in) :: TRANS
    logical, intent(in) :: trans_rhs

    ! local vars
    character :: FACT, EQUED
    integer :: N, NRHS, LDA, LDAF, LDB, LDX, INFO
    real*8, dimension( :, : ), allocatable :: AF, Atmp, Btmp, Xtmp
    integer, dimension(:), allocatable :: IPIV, IWORK
    real*8, dimension(:), allocatable :: R, C, FERR, BERR, WORK
    real*8 :: RCOND

    ! subroutine dgesvx(character FACT,
    !   character TRANS,
    !   integer N,
    !   integer NRHS,
    !   double precision, dimension( lda, * ) A,
    !   integer LDA,
    !   double precision, dimension( ldaf, * ) AF,
    !   integer LDAF,
    !   integer, dimension( * ) IPIV,
    !   character EQUED,
    !   double precision, dimension( * ) R,
    !   double precision, dimension( * ) C,
    !   double precision, dimension( ldb, * ) B,
    !   integer LDB,
    !   double precision, dimension( ldx, * ) X,
    !   integer LDX,
    !   double precision RCOND,
    !   double precision, dimension( * ) FERR,
    !   double precision, dimension( * ) BERR,
    !   double precision, dimension( * ) WORK,
    !   integer, dimension( * ) IWORK,
    !   integer INFO 
    !   )

    ! init
    if ( trans_rhs ) then
       allocate(Btmp(size(B,2), size(B,1)), Xtmp(size(X,2), size(X, 1)))
       Btmp = transpose(B); Xtmp = transpose(X)
    else
       allocate(Btmp(size(B,1), size(B,2)), Xtmp(size(X,1), size(X, 2)))
       Btmp = B; Xtmp = X
    end if

    FACT = 'E'
    N = size(A, 1)
    NRHS = size(Btmp, 2)
    LDA = max(1, N)
    LDAF = max(1, N)
    allocate(AF(LDAF, N), Atmp(size(A, 1), size(A, 2))) 
    Atmp = A
    allocate(IPIV(N))

    allocate(R(N), C(N))
    LDB = max(1, N)
    LDX = max(1, N)
    allocate(FERR(NRHS), BERR(NRHS), WORK(4*N), IWORK(N))

    ! expert solve
    call dgesvx(FACT,TRANS,N,NRHS, Atmp, LDA, AF, LDAF, IPIV &
         , EQUED, R, C, Btmp, LDB, Xtmp, LDX, RCOND, FERR, BERR &
         , WORK, IWORK, INFO) 

    ! bullet proofing
    if ( INFO .ne. 0 ) then 
       print *, 'something is wrong in expert linear' &
            ,' solver dgesvx! INFO = ', INFO, ' stop'
       stop
    end if

    ! extract final sol
    if (trans_rhs) then
       X = transpose(Xtmp)
    else
       X = Xtmp
    end if

    ! clean ups
    deallocate(AF, Atmp, Btmp, Xtmp, IPIV, IWORK, R, C, FERR, BERR, WORK)

    ! done here
  end subroutine wrapper_solve

  subroutine solve_underdet(A, b, x, indx)
    implicit none
    real*8, dimension(:, :), intent(in) :: A
    real*8, dimension(size(A, 1), 1), intent(in) :: b
    real*8, dimension(size(A, 2), 1), intent(out) :: x
    integer, dimension(size(A, 2)), intent(out) :: indx

    ! local vars
    integer :: M, N
    real*8, dimension(:, :), allocatable :: Q, R1, R, btmp, xx

    ! allocate
    allocate(Q(size(A, 1), size(A, 1)), R1(size(A, 1), size(A, 1)))
    allocate(R(size(A, 1), size(A, 2)))
    allocate(btmp(size(b,1), 1), xx(size(b,1), 1))

    ! init
    x = 0.0d0
    M = size(A, 1); N = size(A, 2)

    ! bullet proofing
    if ( M >= N ) then 
       print *, 'the given system is not under determined! stop'
       stop
    end if

    ! first obtain QR decomp. with pivoting
    call wrapper_qr(A, R, indx, Q, .true.)

    ! extracting strong(significant) part
    R1 = R(:, 1:M)

    ! preparing RHS
    btmp = matmul(transpose(Q), b)

    ! solve for significant solution
    call wrapper_solve(R1, xx, btmp, 'N', .false.)

    ! report
    x(indx(1:M), 1) = xx(:, 1)

    ! clean ups
    deallocate(Q, R1, R, btmp, xx)

    ! done here
  end subroutine solve_underdet
  
  subroutine write_pts_to_file(this, filename)
    implicit none
    class(fekete), intent(in) :: this
    character(len=*), intent(in) :: filename

    ! local vars
    integer :: i, j

    ! bullet proofing
    if ( .not. this%initialized ) then
       print *, 'the fekete object must be initialized before exporting! stop'
       stop
    end if

    ! open the output file first
    open(unit = 10, file = filename)

    do j = 1, size(this%fin_coords, 2)

       do i = 1, size(this%fin_coords, 1)
          write(10, '(F30.17, A)', advance = 'no') this%fin_coords(i, j), ' '
       end do

       write(10, '(F30.17)', advance = 'no') this%w(j)
       write(10,*) ! new line

    end do

    ! close the file
    close(10)

    ! done here
  end subroutine write_pts_to_file

  subroutine coord_tri(d, x, y)
    implicit none
    integer, intent(in) :: d
    real*8, dimension(:), allocatable :: x, y

    ! local vars
    integer :: npe, i, j, jj
    real*8 :: dx, dy, xloc, yloc

    npe = (d+1) * (d+2) / 2
    dx = 1.0d0 / dble(d)
    dy = 1.0d0 / dble(d)
    allocate(x(npe), y(npe))
    x = 0.0d0; y = 0.0d0

    jj = 1
    xloc = 1.0d0 
    do i = 0, d
       yloc = 0.0d0
       do j = 0, i
          x(jj) = xloc
          y(jj) = yloc
          yloc = yloc + dy 
          jj = jj + 1
       end do
       xloc = xloc - dx
    end do

    ! done here
  end subroutine coord_tri

  subroutine coord_tet(d, x, y, z)
    implicit none
    integer, intent(in) :: d
    real*8, dimension(:), allocatable :: x, y, z

    ! local vars
    integer :: npe, i, j, k, jj
    real*8 :: dx, dy, dz, xloc, yloc, zloc

    npe = (d+1) * (d+2) * (d+3) / 6
    dx = 1.0d0 / dble(d)
    dy = 1.0d0 / dble(d)
    dz = 1.0d0 / dble(d)

    allocate(x(npe), y(npe), z(npe))
    x = 0.0d0; y = 0.0d0; z = 0.0d0
    jj = 1
    xloc = 1.0d0 
    do i = 0, d
       yloc = 1.0d0 - xloc
       do j = 0, i
          zloc = 1.0d0 - xloc - yloc
          do k = 0, j
             x(jj) = xloc
             y(jj) = yloc
             z(jj) = zloc
             zloc = zloc - dz
             jj = jj + 1
          end do
          yloc = yloc - dy
       end do
       xloc = xloc - dx
    end do

    ! done here
  end subroutine coord_tet

  subroutine coord_hex(n, x, y, z)
    implicit none
    integer, dimension(3), intent(in) :: n
    real*8, dimension(:), allocatable :: x, y, z

    ! local vars
    integer :: npe, i, j, k, jj
    real*8 :: dx, dy, dz, xloc, yloc, zloc

    dx = 2.0d0 / dble(n(1) - 1)
    dy = 2.0d0 / dble(n(2) - 1)
    dz = 2.0d0 / dble(n(3) - 1)

    npe = n(1) * n(2) * n(3)
    allocate(x(npe), y(npe), z(npe))
    x = 0.0d0; y = 0.0d0; z = 0.0d0

    jj = 1
    xloc = -1.0d0 
    do i = 1, n(1) 
       yloc = -1.0d0
       do j = 1, n(2)
          zloc = -1.0d0
          do k = 1, n(3)
             x(jj) = xloc
             y(jj) = yloc
             z(jj) = zloc
             zloc = zloc + dz
             jj = jj + 1
          end do
          yloc = yloc + dy
       end do
       xloc = xloc + dx
    end do

    ! done here
  end subroutine coord_hex

  subroutine coord_prism(d, x, y, z)
    implicit none
    integer, intent(in) :: d
    real*8, dimension(:), allocatable :: x, y, z

    ! local vars
    integer :: npe, i, j, k, jj
    real*8 :: dx, dy, dz, xloc, yloc, zloc

    npe = (d+1) * (d+2) / 2 * (d+1)
    dx = 1.0d0 / dble(d)
    dy = 1.0d0 / dble(d)
    dz = 1.0d0 / dble(d)
    allocate(x(npe), y(npe), z(npe))
    x = 0.0d0; y = 0.0d0; z = 0.0d0

    jj = 1
    xloc = 1.0d0 
    do i = 0, d
       yloc = 0.0d0
       do j = 0, i
          zloc = 0.0d0
          do k = 1, (d+1)
             x(jj) = xloc
             y(jj) = yloc
             z(jj) = zloc
             zloc = zloc + dz 
             jj = jj + 1
          end do
          yloc = yloc + dy 
       end do
       xloc = xloc - dx
    end do

    ! done here
  end subroutine coord_prism

end module approx_fekete

! program tester
!   use approx_fekete
!   implicit none

!   ! local vars
!   type(fekete) :: tfekete

!   print *, '<<< coarse quadrilateral >>>'
!   call tfekete%init(d = 6, name = 'quadri', magnify = 10, s = 10 &
!        , spacing = 'equal_space', echo = .true.)
!   call tfekete%export('quad_d_6.dat')
!   call tfekete%clean()

!   print *, '<<< dense quadrilateral >>>'
!   call tfekete%init(d = 16, name = 'quadri', s = 6, echo = .true.)
!   call tfekete%export('quad_d_16.dat')
!   call tfekete%clean()

!   print *, '<<< one dimensional >>>'
!   call tfekete%init(d = 16, name = '1d', echo = .true.)
!   call tfekete%export('1d_d_16.dat')
!   call tfekete%clean()

!   print *, '<<< dense triangle >>>'
!   call tfekete%init(d = 25, name = 'tri', magnify = 2, s = 6, echo = .true.)
!   call tfekete%export('tri_d_25.dat')
!   call tfekete%clean()

!   print *, '<<< dense tetrahederal >>>'
!   call tfekete%init(d = 8, name = 'tet', magnify = 3, echo = .true.)
!   call tfekete%export('tet_d_8.dat')
!   call tfekete%clean()

!   print *, '<<< hex >>>'
!   call tfekete%init(d = 4, name = 'hex', magnify = 6, s = 3, echo = .true.)
!   call tfekete%export('hex_d_4.dat')
!   call tfekete%clean()

!   print *, '<<< prism >>>'
!   call tfekete%init(d = 8, name = 'prism', magnify = 4, s = 3, echo = .true.)
!   call tfekete%export('prism_d_4.dat')
!   call tfekete%clean()

!   ! done here
! end program tester
