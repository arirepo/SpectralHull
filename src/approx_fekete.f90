module approx_fekete
  use quadri_elem
  use globals
  ! use gen_basis
  implicit none
  private

  ! data type for handeling complex 2d edge geometies for master element
  type edg

     real*8 :: x1, y1
     real*8 :: x2, y2
     real*8 :: nx, ny
     integer :: tag

  end type edg

  type edgs

     type(edg), dimension(:), allocatable :: val
     integer :: maxtag
     real*8 :: tol

   contains

     procedure, public :: add => add_edg
     procedure, public :: sample => sample_edgs

  end type edgs

  type fekete
     private

     ! private components
     integer, dimension(:), allocatable :: n
     integer :: nw, d, d_big, d_orig, npedg, nedg
     real*8, dimension(:), allocatable :: coords
     real*8, dimension(:), allocatable :: rhs
     real*8, dimension(:, :), allocatable :: x0
     integer ::  magnify, s
     character(len = 128) :: spacing
     logical :: initialized = .false., echo = .false.
     character(len = 128) :: name
     type(edgs) :: tedgs

     ! public components
     real*8, dimension(:), allocatable, public :: xy
     procedure(gen_xy), public, pointer, nopass  :: comp_xy, moment
     real*8, dimension(:), allocatable, public :: w
     real*8, dimension(:, :), allocatable, public :: fin_coords, xy0
     class(fekete), pointer :: one_d => null()

   contains

     ! procedures
     procedure, public :: init => initialize
     procedure, public :: set_pt => set_point_coord
     procedure, private :: find => find_approx_fekete_pts
     procedure, public :: clean => clean_fekete_obj
     procedure, public :: export => write_pts_to_file
     procedure, public :: export_custom2D => write_custom2d_to_file
     ! procedure, public :: interp => sample_interpolate_export

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

  real*8, parameter :: HARD_TOL = 1.0d-12

  ! parameters for gravitational distribution
  type grav_param

     real*8 :: z ! gravitational law power
     real*8 :: dt ! pseudo time step size
     integer :: itrmax ! maximum number of pseudo time steps
     real*8 :: conv_epsil ! the successive convergence epsilon
     logical :: tightly_coupled = .false. ! set true for stretching 
                                          ! near boundaries in convex shapes
     logical :: no_fekete = .false. ! set true if only radial dist. is wanted 
  end type grav_param

  public :: fekete, fekete_table, edgs, read_segment_file
  public :: grav_param
  public :: grav_animation
  public :: coord_tri, edg

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

  recursive subroutine initialize(this, d, name, magnify, s, spacing, echo, tedgs, tgrav_param)
    implicit none
    class(fekete), intent(inout) :: this
    integer, intent(in) :: d
    character(len = *), intent(in) :: name
    integer, intent(in), optional :: magnify, s
    character(len = *), intent(in), optional :: spacing
    logical, intent(in), optional :: echo
    type(edgs), intent(in), optional :: tedgs
    type(grav_param), optional :: tgrav_param

    ! local vars
    integer :: dim
    real*8, dimension(:), allocatable :: x, y, z
    type(fekete), allocatable, target :: tmp_one_d
    integer :: nbn
    real*8 :: min_dist, max_dist, distance

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

    case ('custom2d') ! 2D general custom shape master element! 

       ! bullet proofing 
       if ( .not. present(tedgs) ) then
          print *, 'the tedgs structure containing edge list' &
               , ' for general custom2d element must be present' &
               , ' in the initializer subroutine argument list! stop'
          stop
       end if
  
       ! init.
       dim = 2
       allocate(this%n(dim))
       this%npedg = d+1
       this%nedg = size(tedgs%val)
       this%n(1) = ceiling(dble(this%nedg)/4.0D0 * dble(this%npedg))
       this%n(2) = this%n(1) 
       this%nw = this%n(1) * this%n(2)
       allocate(this%coords(dim))
       allocate(this%w(this%nw), this%rhs(this%nw), this%xy(this%nw))
       allocate(this%fin_coords(dim, this%nw), this%xy0(this%nw, this%nw))

       ! store the edge list
       this%tedgs = tedgs

       ! generate master element point distribution of the dense magnified grid
       ! call gen_master_custom2d_uniform(this, (this%magnify * this%npedg), this%x0)
       call gen_master_custom2d_uniform(this, (this%magnify * this%npedg), this%x0, nbn, min_dist, max_dist)
       ! call gen_master_custom2d(this, (this%magnify * this%npedg), this%x0, nbn, min_dist, max_dist)
 
       ! apply gravitational distribution to smoothen the spacing (if requested)
       if ( present(tgrav_param) ) then
          if ( tgrav_param%tightly_coupled ) then
             distance = max_dist
          else
             distance = min_dist
          end if
          ! subroutine gravitional_equilib(x, nbn, tedgs, z, distance, dt, itrmax, conv_epsil)
          call gravitional_equilib(this%x0, nbn, this%tedgs, tgrav_param%z, distance &
               , tgrav_param%dt, tgrav_param%itrmax, tgrav_param%conv_epsil, min_dist)

          if ( tgrav_param%no_fekete ) then ! finish initialization here
             deallocate(this%fin_coords, this%w)
             allocate(this%fin_coords(dim, size(this%x0,2)), this%w(size(this%x0,2)))
             this%fin_coords = this%x0
             this%w = 0.0d0
             this%initialized = .true.
             ! call this%export('custom.dat')
             ! print *, 'hey ...'
             return ! interupt
          end if
          ! stop

       end if

       ! initialize the one-dimensional Gauss-Lobatto rule
       ! for boundary integrals
       allocate(tmp_one_d)
       call tmp_one_d%init(d = (this%n(1) + this%n(2)+1), magnify = 10, s = 6, name = '1d', echo = .false.)
       this%one_d => tmp_one_d

       ! setting up pointers to appropriate xy eval and moment subroutines
       this%comp_xy => xy_at_pt_quadri
       this%moment => moment_custom2d


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
    if(allocated(this%tedgs%val)) deallocate(this%tedgs%val)
    ! if(associated(this%one_d)) call this%one_d%clean()
 
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

  subroutine moment_custom2d(this)
    implicit none
    type(fekete), intent(inout) :: this

    ! local vars
    integer :: i, j, K 

    K = 1
    do i = 1, this%n(1)
       do j = 1, this%n(2)
          this%rhs(K) = polygon_bn_integ(this, i, j)
          K = K + 1
       end do
    end do

    ! done here
  end subroutine moment_custom2d

  function polygon_bn_integ(this, i, j)
    implicit none
    type(fekete), intent(in), target :: this
    integer, intent(in) :: i, j
    real*8 :: polygon_bn_integ

    ! local vars
    integer :: k, l
    integer :: nedgs
    real*8 :: xkl, ykl, nuij
    real*8, dimension(:), pointer :: t => null(), w => null()
    type(edg), pointer :: ed => null()
    real*8 :: ds

    ! HARD reset
    nuij = 0.0d0

    nedgs = size(this%tedgs%val)
    t => this%one_d%fin_coords(1, :)
    w => this%one_d%w

    do k = 1, nedgs

       ed => this%tedgs%val(k)

       ds = sqrt((ed%x2 - ed%x1)**2 + (ed%y2 - ed%y1)**2)

       do l = 1, size(t)
          xkl = ed%x1 + 0.5d0 * (t(l) + 1.0d0) * ( ed%x2 - ed%x1 )
          ykl = ed%y1 + 0.5d0 * (t(l) + 1.0d0) * ( ed%y2 - ed%y1 )

          nuij = nuij + ( ( (xkl**i) * (ykl**(j-1)) ) / dble(i) * ed%nx &
               + ( (xkl**(i-1)) * (ykl**j) ) / dble(j) * ed%ny ) * 0.5d0 * ds * w(l) 
       end do

    end do

    polygon_bn_integ = 0.5d0 * nuij

    ! done here
  end function polygon_bn_integ

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

    if ( allocated(x) ) deallocate(x)
    if ( allocated(y) ) deallocate(y)
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

  subroutine gen_master_custom2d_uniform(this, npedg, x0, nbn, min_dist, max_dist)  
    implicit none
    class(fekete), intent(inout), target :: this
    integer, intent(in) :: npedg 
    real*8, dimension(:,:), allocatable :: x0 ! output
    integer, intent(out), optional :: nbn
    real*8, intent(out), optional :: min_dist, max_dist

    ! local vars
    integer :: m, k, tnp, nn
    integer :: i, j
    type(edg), pointer :: ed => null()
    real*8 :: length, xk, yk, minlen
    real*8 :: xmin, xmax, ymin, ymax, delta
    real*8, dimension(size(this%tedgs%val)) :: lengths

    ! compute unit edge normal vectors
    ! and min edge length 
    do m = 1, size(this%tedgs%val) ! all edges

       ! grab edge # m
       ed => this%tedgs%val(m)
       ! compute normals
       ed%nx = ed%y2 - ed%y1
       ed%ny = ed%x1 - ed%x2
       length = sqrt( ed%nx**2 + ed%ny**2 )
       ed%nx = ed%nx / length
       ed%ny = ed%ny / length 

       ! store edge length
       lengths(m) = length

       ! detect and save minimum edge length
       ! which will be used for point distribution
       if ( m .eq. 1 ) then
          minlen = length
       else
          if (length < minlen) then
             minlen = length
          end if
       end if

    end do

    ! first add edge vertices
    do m = 1, size(this%tedgs%val)

       ! grab edge # m
       ed => this%tedgs%val(m)

       ! add pt1
       xk = ed%x1  
       yk = ed%y1
       if ( is_pt_in_array(x0, (/ xk, yk/), HARD_TOL ) .eq. 0 ) then
          call add_double(x0, (/ xk, yk/) ) 
       end if

       ! add pt2
       xk = ed%x2  
       yk = ed%y2
       if ( is_pt_in_array(x0, (/ xk, yk/), HARD_TOL ) .eq. 0 ) then
          call add_double(x0, (/ xk, yk/) ) 
       end if

    end do

    ! then add interior points per edges
    do m = 1, size(this%tedgs%val)

       ! grab edge # m
       ed => this%tedgs%val(m)
       tnp = nint(lengths(m) / minlen * dble(npedg))

       do k = 2, (tnp-1)

          xk = ed%x1 + dble(k-1) / dble(tnp - 1) * (ed%x2 - ed%x1)  
          yk = ed%y1 + dble(k-1) / dble(tnp - 1) * (ed%y2 - ed%y1)

          if ( is_pt_in_array(x0, (/ xk, yk/), HARD_TOL ) .eq. 0 ) then
             call add_double(x0, (/ xk, yk/) ) 
          end if

       end do

    end do

    ! return the number of boundary points
    ! if requested to do so 
    if ( present( nbn ) ) nbn = size(x0, 2)

    ! finally, add interior (bubble) points
    ! find the box around the generic polygon
    xmin = minval(x0(1, :)); xmax = maxval(x0(1, :))
    ymin = minval(x0(2, :)); ymax = maxval(x0(2, :))
    delta = minlen / dble(npedg - 1)
    nn = ceiling( max( (xmax-xmin)/delta , (ymax-ymin)/delta ) ) + 1

    ! return min_dist and max_dist if requested
    if ( present(min_dist) ) min_dist = 2.0d0 * delta
    if ( present(max_dist) ) max_dist = 1.2d0 * sqrt((xmax - xmin)**2 + (ymax-ymin)**2)
  
    ! add uniform points inside the big box containing all edges
    do i = 1, nn
       xk = xmin + dble(i-1) / dble(nn-1) * (xmax - xmin)
       do j = 1, nn
          yk = ymin + dble(j-1) / dble(nn-1) * (ymax - ymin)

          if ( is_inside(this%tedgs, xk, yk) .and. &
               ( is_pt_in_array(x0, (/ xk, yk/), (delta/10.0d0)) .eq. 0) ) then
             call add_double(x0, (/ xk, yk/) )
          end if

       end do
    end do

    ! done here
  end subroutine gen_master_custom2d_uniform

  subroutine gen_master_custom2d(this, npedg, x0, nbn, min_dist, max_dist)  
    implicit none
    class(fekete), intent(inout), target :: this
    integer, intent(in) :: npedg 
    real*8, dimension(:,:), allocatable :: x0 ! output
    integer, intent(out), optional :: nbn
    real*8, intent(out), optional :: min_dist, max_dist

    ! local vars
    integer :: m, k, tnp, nn
    type(edg), pointer :: ed => null()
    real*8 :: length, xk, yk, minlen
    real*8 :: xmin, xmax, ymin, ymax, delta
    real*8, dimension(size(this%tedgs%val)) :: lengths
    integer :: ninter, jj
    real*8 :: rnd

    ! compute unit edge normal vectors
    ! and min edge length 
    do m = 1, size(this%tedgs%val) ! all edges

       ! grab edge # m
       ed => this%tedgs%val(m)
       ! compute normals
       ed%nx = ed%y2 - ed%y1
       ed%ny = ed%x1 - ed%x2
       length = sqrt( ed%nx**2 + ed%ny**2 )
       ed%nx = ed%nx / length
       ed%ny = ed%ny / length 

       ! store edge length
       lengths(m) = length

       ! detect and save minimum edge length
       ! which will be used for point distribution
       if ( m .eq. 1 ) then
          minlen = length
       else
          if (length < minlen) then
             minlen = length
          end if
       end if

    end do

    ! first add edge vertices
    do m = 1, size(this%tedgs%val)

       ! grab edge # m
       ed => this%tedgs%val(m)

       ! add pt1
       xk = ed%x1  
       yk = ed%y1
       if ( is_pt_in_array(x0, (/ xk, yk/), HARD_TOL ) .eq. 0 ) then
          call add_double(x0, (/ xk, yk/) ) 
       end if

       ! add pt2
       xk = ed%x2  
       yk = ed%y2
       if ( is_pt_in_array(x0, (/ xk, yk/), HARD_TOL ) .eq. 0 ) then
          call add_double(x0, (/ xk, yk/) ) 
       end if

    end do

    ! then add interior points per edges
    do m = 1, size(this%tedgs%val)

       ! grab edge # m
       ed => this%tedgs%val(m)
       tnp = nint(lengths(m) / minlen * dble(npedg))

       do k = 2, (tnp-1)

          xk = ed%x1 + dble(k-1) / dble(tnp - 1) * (ed%x2 - ed%x1)  
          yk = ed%y1 + dble(k-1) / dble(tnp - 1) * (ed%y2 - ed%y1)

          if ( is_pt_in_array(x0, (/ xk, yk/), HARD_TOL ) .eq. 0 ) then
             call add_double(x0, (/ xk, yk/) ) 
          end if

       end do

    end do

    ! return the number of boundary points
    ! if requested to do so 
    if ( present( nbn ) ) nbn = size(x0, 2)

    ! finally, add interior (bubble) points
    ! find the box around the generic polygon
    xmin = minval(x0(1, :)); xmax = maxval(x0(1, :))
    ymin = minval(x0(2, :)); ymax = maxval(x0(2, :))
    delta = minlen / dble(npedg - 1)
    nn = ceiling( max( (xmax-xmin)/delta , (ymax-ymin)/delta ) ) + 1

    ! return min_dist and max_dist if requested
    if ( present(min_dist) ) min_dist = 2.0d0 * delta
    if ( present(max_dist) ) max_dist = 1.2d0 * sqrt((xmax - xmin)**2 + (ymax-ymin)**2)
  
    ! determine the number of interior points
    ninter = nn * nn

    ! add random points
    jj = 0
    do

       ! random x
       call random_number(rnd)
       xk = xmin + rnd * (xmax - xmin)

       ! random y 
       call random_number(rnd)
       yk = ymin + rnd * (ymax - ymin)

       if ( is_inside(this%tedgs, xk, yk) .and. &
            ( is_pt_in_array(x0, (/ xk, yk/), this%tedgs%tol) .eq. 0 ) ) then
          call add_double(x0, (/ xk, yk/) )
          jj = jj + 1
       end if

       if (jj .eq. ninter) exit

    end do

    ! done here
  end subroutine gen_master_custom2d

  ! finds a tuple coords in a arrays of coords
  ! if found returns the location else returns 0 
  function is_pt_in_array(x0, x, tol)
    implicit none
    real*8, dimension(:, :), allocatable :: x0
    real*8, dimension(2), intent(in) :: x
    real*8, intent(in), optional :: tol
    integer :: is_pt_in_array

    ! local vars
    integer :: i

    is_pt_in_array = 0

    ! bullet proofing ...
    if( .not. allocated(x0) ) return
 
    do i = 1, size(x0, 2)

       ! standard check
       if ( all(x0(:, i) .eq. x) ) then
          is_pt_in_array = i
          return
       end if

       ! check with tolerance option if present
       if ( present(tol) ) then
          if (sqrt(sum((x0(:, i) - x)**2)) <= tol) then
             is_pt_in_array = i
             return
          end if
       end if

    end do

    ! done here
  end function is_pt_in_array

  ! add tuple xend to the end of x
  subroutine add_double(x, xend)
    implicit none
    real*8, dimension(:,:), allocatable :: x
    real*8, dimension(2), intent(in) :: xend

    integer :: n
    real*8, dimension(:,:), allocatable :: tmp

    if (.not. allocated(x)) then
       n = 0
    else
       n = size(x, 2)
    end if

    n = n + 1

    allocate(tmp(2, n))

    if (n > 1) tmp(:, 1:(n-1)) = x(:, 1:(n-1))

    tmp(:,n) = xend

    call move_alloc(tmp, x) 

    ! done here
  end subroutine add_double

  ! ! checks to see if point (xk, yk) is inside
  ! ! the given polygon. if yes, return .true.
  ! ! otherwise returns .false.
  ! function is_inside(tedgs, xk, yk)
  !   implicit none
  !   type(edgs), intent(in), target :: tedgs
  !   real*8, intent(in) :: xk, yk
  !   logical :: is_inside

  !   ! local vars
  !   integer :: m
  !   real*8 :: rx, ry, xc, yc, length, dot
  !   type(edg), pointer :: ed => null()
  !   integer, dimension(tedgs%maxtag, 2) :: signs ! +, -

  !   signs = 0

  !   do m = 1, size(tedgs%val) ! check all edges

  !      ed => tedgs%val(m)
  !      xc = 0.5d0 * (ed%x1 + ed%x2)
  !      yc = 0.5d0 * (ed%y1 + ed%y2)
  !      rx = xc - xk
  !      ry = yc - yk
  !      length = sqrt(rx * rx + ry * ry)
  !      rx = rx / length
  !      ry = ry / length
  !      dot = rx*ed%nx + ry*ed%ny

  !      ! quick exit
  !      if ( (ed%tag .eq. 1) .and. (dot <= 0.0d0) ) then
  !         is_inside = .false.
  !         return
  !      end if

  !      if (dot > 0.0d0) then
  !         signs(ed%tag, 1) = signs(ed%tag, 1) + 1 !add to plus 
  !      else
  !         signs(ed%tag, 2) = signs(ed%tag, 2) + 1 !add to minus
  !      end if

  !   end do

  !   ! decide
  !   if (tedgs%maxtag > 1 ) then
  !      if ( any(signs(2:,:) .eq. 0) ) then
  !         is_inside = .false.
  !         return
  !      end if
  !   end if

  !   is_inside = .true. ! otherwise not proven

  !   ! done here
  ! end function is_inside

  ! checks to see if point (x, y) is inside
  ! the given polygon. if yes, return .true.
  ! otherwise returns .false.
  function is_inside(tedgs, x, y)
    implicit none
    type(edgs), intent(in), target :: tedgs
    real*8, intent(in) :: x, y
    logical :: is_inside

    ! local vars
    integer :: m
    type(edg), pointer :: ed => null()
    real*8 :: xi, yi, xj, yj

    is_inside = .false.

    do m = 1, size(tedgs%val) ! check all edges

       ed => tedgs%val(m)
       xi = ed%x2; yi = ed%y2
       xj = ed%x1; yj = ed%y1
       if ( ((yi<y) .and. (yj>=y)) .or. ((yj<y) .and. (yi>=y)) ) then
          if ( (xi+(y-yi)/(yj-yi)*(xj-xi)) < x ) then
             is_inside = .not. is_inside 
          end if
       end if

    end do

    ! done here
  end function is_inside

  ! adds one edg to the end of edgs list
  subroutine add_edg(this, tedg)
    implicit none
    class(edgs), intent(inout) :: this
    type(edg), intent(in) :: tedg

    integer :: nedgs
    type(edg), dimension(:), allocatable :: tmp

    if ( .not. allocated(this%val) ) then
       nedgs = 0
    else 
       nedgs = size(this%val)
    end if

    nedgs = nedgs + 1

    allocate(tmp(nedgs))
    if (nedgs > 1) tmp(1:(nedgs-1)) = this%val(1:(nedgs-1))

    tmp(nedgs) = tedg

    call move_alloc(tmp, this%val)

    ! done here
  end subroutine add_edg

  ! creates a sample edg list
  subroutine sample_edgs(tedgs)
    implicit none
    class(edgs), intent(inout) :: tedgs

    ! local vars
    type(edg) :: tedg

    ! ==============================
    ! Big "T"
    ! ==============================

    ! edg1
    tedg%x1 = 0.0d0; tedg%y1 = 0.0d0
    tedg%x2 = 0.5d0; tedg%y2 = 0.0d0
    tedg%tag = 1
    call tedgs%add(tedg)

    ! edg2
    tedg%x1 = 0.5d0; tedg%y1 = 0.0d0
    tedg%x2 = 0.5d0; tedg%y2 = 2.0d0
    tedg%tag = 1
    call tedgs%add(tedg)

    ! edg3
    tedg%x1 = 0.5d0; tedg%y1 = 2.0d0
    tedg%x2 = 1.0d0; tedg%y2 = 2.0d0
    tedg%tag = 1
    call tedgs%add(tedg)

    ! edg4
    tedg%x1 = 1.0d0; tedg%y1 = 2.0d0
    tedg%x2 = 1.0d0; tedg%y2 = 2.5d0
    tedg%tag = 1
    call tedgs%add(tedg)

    ! edg5
    tedg%x1 = 1.0d0; tedg%y1 = 2.5d0
    tedg%x2 = -1.0d0; tedg%y2 = 2.5d0
    tedg%tag = 1
    call tedgs%add(tedg)

    ! edg6
    tedg%x1 = -1.0d0; tedg%y1 = 2.5d0
    tedg%x2 = -1.0d0; tedg%y2 = 2.0d0
    tedg%tag = 1
    call tedgs%add(tedg)

    ! edg7
    tedg%x1 = -1.0d0; tedg%y1 = 2.0d0
    tedg%x2 = -0.5d0; tedg%y2 = 2.0d0
    tedg%tag = 1
    call tedgs%add(tedg)

    ! edg8
    tedg%x1 = -0.5d0; tedg%y1 = 2.0d0
    tedg%x2 = -0.5d0; tedg%y2 = 0.0d0
    tedg%tag = 1
    call tedgs%add(tedg)

    ! edg9
    tedg%x1 = -0.5d0; tedg%y1 = 0.0d0
    tedg%x2 = 0.0d0; tedg%y2 = 0.0d0
    tedg%tag = 1
    call tedgs%add(tedg)

    ! ! ==============================
    ! ! simple triangle
    ! ! ==============================
    ! ! edg1
    ! tedg%x1 = 0.0d0; tedg%y1 = 0.0d0
    ! tedg%x2 = 1.0d0; tedg%y2 = 0.0d0
    ! tedg%tag = 1
    ! call tedgs%add(tedg)

    ! ! edg2
    ! tedg%x1 = 1.0d0; tedg%y1 = 0.0d0
    ! tedg%x2 = 0.0d0; tedg%y2 = 1.0d0
    ! tedg%tag = 1
    ! call tedgs%add(tedg)

    ! ! edg3
    ! tedg%x1 = 0.0d0; tedg%y1 = 1.0d0
    ! tedg%x2 = 0.0d0; tedg%y2 = 0.0d0
    ! tedg%tag = 1
    ! call tedgs%add(tedg)


    ! ! ==============================
    ! ! hollow rectangle
    ! ! ==============================
    ! ! edg1
    ! tedg%x1 = 0.0d0; tedg%y1 = 0.0d0
    ! tedg%x2 = 1.0d0; tedg%y2 = 0.0d0
    ! tedg%tag = 1
    ! call tedgs%add(tedg)

    ! ! edg2
    ! tedg%x1 = 1.0d0; tedg%y1 = 0.0d0
    ! tedg%x2 = 1.0d0; tedg%y2 = 1.0d0
    ! tedg%tag = 1
    ! call tedgs%add(tedg)

    ! ! edg3
    ! tedg%x1 = 1.0d0; tedg%y1 = 1.0d0
    ! tedg%x2 = 0.0d0; tedg%y2 = 1.0d0
    ! tedg%tag = 1
    ! call tedgs%add(tedg)

    ! ! edg4
    ! tedg%x1 = 0.0d0; tedg%y1 = 1.0d0
    ! tedg%x2 = 0.0d0; tedg%y2 = 0.0d0
    ! tedg%tag = 1
    ! call tedgs%add(tedg)

    ! ! edg5
    ! tedg%x1 = 0.25d0; tedg%y1 = 0.25d0
    ! tedg%x2 = 0.75d0; tedg%y2 = 0.25d0
    ! tedg%tag = 2
    ! call tedgs%add(tedg)

    ! ! edg6
    ! tedg%x1 = 0.75d0; tedg%y1 = 0.25d0
    ! tedg%x2 = 0.75d0; tedg%y2 = 0.75d0
    ! tedg%tag = 2
    ! call tedgs%add(tedg)

    ! ! edg7
    ! tedg%x1 = 0.75d0; tedg%y1 = 0.75d0
    ! tedg%x2 = 0.25d0; tedg%y2 = 0.75d0
    ! tedg%tag = 2
    ! call tedgs%add(tedg)

    ! ! edg8
    ! tedg%x1 = 0.25d0; tedg%y1 = 0.75d0
    ! tedg%x2 = 0.25d0; tedg%y2 = 0.25d0
    ! tedg%tag = 2
    ! call tedgs%add(tedg)
 
    !set maximum number of tags
    ! tedgs%maxtag = 2
    tedgs%maxtag = 1

    tedgs%tol = 2.6d-2

    ! call sample_nsides(tedgs, 16)

    ! done here
  end subroutine sample_edgs

  ! creates a sample edg list for
  ! n-sided element
  subroutine sample_nsides(tedgs, n)
    implicit none
    class(edgs), intent(inout) :: tedgs
    integer, intent(in) :: n ! num of sides

    ! local vars
    integer :: m
    type(edg) :: tedg
    real*8 :: theta, theta_m

    theta = dble(2) * PI / dble(n)
    tedg%x1 = cos(-theta); tedg%y1 = sin(-theta)
    tedg%tag = 1

    do m = 1, n

       theta_m = dble(m-1) * theta
       tedg%x2 = cos(theta_m); tedg%y2 = sin(theta_m)

       call tedgs%add(tedg)

       tedg%x1 = tedg%x2; tedg%y1 = tedg%y2 

    end do

    !set maximum number of tags
    tedgs%maxtag = 1
    tedgs%tol = 2.6d-2

    ! done here
  end subroutine sample_nsides

  ! add "integer" tuple xend to the end of x
  subroutine add_integer(x, xend)
    implicit none
    integer, dimension(:,:), allocatable :: x
    integer, dimension(2), intent(in) :: xend

    integer :: n
    integer, dimension(:,:), allocatable :: tmp

    if (.not. allocated(x)) then
       n = 0
    else
       n = size(x, 2)
    end if

    n = n + 1

    allocate(tmp(2, n))

    if (n > 1) tmp(:, 1:(n-1)) = x(:, 1:(n-1))

    tmp(:,n) = xend

    call move_alloc(tmp, x) 

    ! done here
  end subroutine add_integer

  subroutine add_bn_edg_pts_to_x0(this, x0, bn_edgs)
    implicit none
    class(fekete), intent(in), target :: this
    real*8, dimension(:,:), allocatable :: x0
    integer, dimension(:,:), allocatable :: bn_edgs

    ! local vars
    integer :: i, pt1, pt2
    type(edg), pointer :: ed => null()
    real*8 :: xi, yi

    do i = 1, size(this%tedgs%val) ! all bn edgs

       ed => this%tedgs%val(i)
       xi = ed%x1; yi = ed%y1
       pt1 = is_pt_in_array(x0, (/xi, yi/), 1.0d-14)
       if ( pt1 .eq. 0) then
          call add_double(x0, (/ xi, yi/) )
          pt1 = size(x0, 2)
       end if

       xi = ed%x2; yi = ed%y2
       pt2 = is_pt_in_array(x0, (/xi, yi/), 1.0d-14)
       if ( pt2 .eq. 0) then
          call add_double(x0, (/ xi, yi/) )
          pt2 = size(x0, 2)
       end if

       call add_integer(bn_edgs, (/pt1, pt2/))

    end do

    ! done here
  end subroutine add_bn_edg_pts_to_x0

  subroutine write_advanced_to_file(ptfile, edgfile, x0, val, bn_edgs)
    implicit none
    character(len=*), intent(in) :: ptfile, edgfile
    real*8, dimension(:,:), intent(in) :: x0
    real*8, dimension(:), intent(in) :: val
    integer, dimension(:,:), intent(in) :: bn_edgs

    ! local vars
    integer :: i, j

    ! first write the final points
    ! open the output file first
    open(unit = 10, file = ptfile)

    do j = 1, size(x0, 2)

       do i = 1, size(x0, 1)
          write(10, '(F30.17, A)', advance = 'no') x0(i, j), ' '
       end do

       write(10, '(F30.17)', advance = 'no') val(j)

       write(10,*) ! new line

    end do

    ! close the file
    close(10)

    ! then write the boundary edges
    ! open the output file first
    open(unit = 10, file = edgfile)

    do j = 1, size(bn_edgs, 2)

       do i = 1, size(bn_edgs, 1)
          write(10, '(I8, A)', advance = 'no') bn_edgs(i, j), ' '
       end do

       write(10,*) ! new line

    end do

    ! close the file
    close(10)

    ! done here
  end subroutine write_advanced_to_file

  subroutine write_custom2d_to_file(this, filename)
    implicit none
    class(fekete), intent(in) :: this
    character(len=*), intent(in) :: filename

    ! local vars
    real*8, dimension(:, :), allocatable :: x0
    real*8, dimension(:), allocatable :: val
    integer, dimension(:, :), allocatable :: bn_edgs

    ! bullet proofing
    if ( .not. this%initialized ) then
       print *, 'the fekete object must be initialized before exporting! stop'
       stop
    end if

    if ( this%name .ne. 'custom2d' ) then
       print *, 'write_custom2d_to_file(...) procedure is only for custom2d elements! stop'
       stop
    end if

    ! add boundary edges points to fekete points just
    ! for visuallization
    allocate(x0(size(this%fin_coords, 1), size(this%fin_coords, 2)))
    x0 = this%fin_coords
    call add_bn_edg_pts_to_x0(this, x0, bn_edgs)
    ! prepare the value for contouring
    ! here is the weight which is zero for
    ! non fekete points (bn edge points)
    allocate(val(size(x0, 2)))
    val(1:size(this%fin_coords, 2)) = this%w
    if (size(this%fin_coords, 2) < size(x0, 2)) then 
       val((size(this%fin_coords, 2)+1):) = 0.0d0
    end if

    ! then write advanced to file
    call  write_advanced_to_file(filename, (filename//'.edg.dat'), x0, val, bn_edgs)

    ! clean ups
    deallocate(x0, val, bn_edgs)

    ! done here
  end subroutine write_custom2d_to_file

  ! ! sample interpolation and export subroutine
  ! subroutine sample_interpolate_export(this, zoom, err, filename, p, use_RBF)
  !   implicit none
  !   class(fekete), intent(inout) :: this
  !   integer, intent(in) :: zoom
  !   real*8, optional :: err
  !   character(len=*), intent(in), optional :: filename
  !   integer, optional :: p
  !   logical, intent(in), optional :: use_RBF

  !   ! local vars
  !   real*8, dimension(:, :), allocatable :: x0
  !   integer, dimension(:, :), allocatable :: bn_edgs
  !   real*8, dimension(:), allocatable :: val, psi
  !   real*8, dimension(:), allocatable :: u0, famps
  !   type(basis) :: tbasis
  !   integer :: npts, k
  !   real*8 :: xk, yk, exact_val(1)

  !   ! bullet proofing
  !   if ( .not. this%initialized ) then
  !      print *, 'the fekete object must be initialized before exporting! stop'
  !      stop
  !   end if

  !   if ( this%name .ne. 'custom2d' ) then
  !      print *, 'write_custom2d_to_file(...) procedure is only for custom2d elements! stop'
  !      stop
  !   end if

  !   ! initialize the nodal values at Fekete points
  !   allocate(u0(size(this%fin_coords, 2)), famps(size(this%fin_coords, 2)))
  !   call fxy(this%fin_coords(1,:), this%fin_coords(2,:), u0)

  !   ! generate the set of nodal points for 
  !   ! doing interpolation
  !   call gen_master_custom2d_uniform(this, (this%npedg * zoom), x0)

  !   ! add boundary edges points to the selected interpolation points just
  !   ! for visuallization of interpolated values
  !   call add_bn_edg_pts_to_x0(this, x0, bn_edgs)

  !   ! prepare the value for contouring
  !   npts = size(x0, 2)
  !   allocate(val(npts))
  !   allocate(psi(size(this%fin_coords,2)))

  !   ! construct the set of basis functions using Fekete points
  !   call tbasis%init(x = this%fin_coords(1,:), y = this%fin_coords(2,:), elname = 2, do_radial = use_RBF)

  !   ! print Fourier amplitudes
  !   if (this%echo) then
  !      call tbasis%fourier_amps(u0, famps)
  !      print *, 'Fourier Amplitudes = ', famps
  !   end if

  !   ! evaluate the basis functions at all interpolation points
  !   if (present(err)) err = 0.0d0
  !   do k = 1, npts
  !      xk = x0(1, k); yk = x0(2, k)
  !      call tbasis%eval(xk, yk, 0, psi, resol = 1.0d0)
  !      val(k) = sum(psi * u0) ! the first basis function
  !      ! compute the error
  !      if ( present(err) ) then
  !         call fxy( (/ xk /), (/ yk /), exact_val)
  !         err = err + (exact_val(1) - val(k))**2
  !      end if
  !   end do

  !   ! if error calculation is asked for
  !   if ( present(err) ) err = sqrt(err)

  !   ! then write advanced to file (if asked for)
  !   if ( present(filename) ) then
  !      call  write_advanced_to_file(filename, (filename//'.edg.dat') &
  !           , x0, val, bn_edgs)
  !   end if

  !   ! return degree p if available
  !   if (present(p)) p = nint(sqrt(dble(size(this%fin_coords,2))))

  !   ! clean ups
  !   deallocate(x0, bn_edgs, val, psi, u0, famps)

  !   ! done here
  ! end subroutine sample_interpolate_export

  ! generates a sample f(x,y) distribution
  ! over the element which is suitable to test
  ! interpolation and filtering
  subroutine fxy(x, y, u)
    implicit none
    real*8, dimension(:), intent(in) :: x, y
    real*8, dimension(:), intent(out) :: u

    !u = BESSEL_JN( 0, 25.0d0 * sqrt(x**2 + y**2)) * (x+y) * (x-y)  !((x - 0.25d0) * (x**2 - 0.75d0) * (y**2 - 0.5d0) * (x + y - 0.1d0) ) * BESSEL_JN (0, 10.0d0 * sqrt(x**2 + y**2))!&
         !/ (3.0d0 * x**2 + 2.0d0 * y**4 + 1.0d0) !sin(3 * PI * x) * sin(3 * PI * y)

    u = sin(2.0d0 * PI * x) * sin(2.0d0 * PI * y)

    ! done here
  end subroutine fxy

  ! main subroutine for reading segment files
  ! in the form :
  !
  !
  ! number of points of connector 1
  ! x1 y1 z1
  ! x2 y2 z2
  ! .
  ! .
  ! number of points of connector 2
  ! x1 y1 z1
  ! x2 y2 z2
  ! .
  ! .
  ! where all points are SEQUNTIAL
  
  subroutine read_segment_file(infile, tedgs)
    implicit none
    character(len = *), intent(in) :: infile
    class(edgs), intent(inout) :: tedgs

    ! local vars
    integer :: i, j
    integer :: istat
    integer :: npt
    real*8 :: xx, yy, zz
    real*8 :: xx1, yy1, zz1
    real*8 :: xx2, yy2, zz2
    integer :: tot_curves, tot_seg
    type(edg) :: tedg

    ! opening the input file
    open ( unit=9, file=infile , status = 'old', &
         iostat = istat)
    if ( istat /= 0) then 
       print *, 'fatal: could not open <', infile, '> file! exit'
       stop
    end if

    ! first count the total number of boundary curves (connectors)
    ! and the total number of segments
    tot_curves = 0; tot_seg = 0

    ! read all lines of the input file
    do
       read(9,*,iostat = istat) npt 
       if(istat > 0) then
          print *, 'fatal error : file <', infile &
               ,'> is bad or corrupted. unable to read. exit.'
          stop
       end if
       if( istat < 0 ) exit ! EOF encountered. exit loop

       tot_curves = tot_curves + 1
       tot_seg = tot_seg + npt - 1

       do i = 1, npt
          read(9,*,iostat = istat) xx, yy, zz
       end do

    end do ! reading finishes after this completes

    ! print *, 'total boundary curves (connector) : ', tot_curves
    ! print *, 'total boundary segments : ', tot_seg

    rewind( unit = 9) ! go back to header

    ! read all lines of the input file (all connectors)
    do j = 1, tot_curves

       read(9,*,iostat = istat) npt 
       ! print *, 'for bn_curves(',j,')'

       !setting edge tag to the connector number of these edges
       tedg%tag = j 
       read(9,*,iostat = istat) xx1, yy1, zz1
       do i = 2, npt

          tedg%x1 = xx1; tedg%y1 = yy1 
          read(9,*,iostat = istat) xx2, yy2, zz2
          tedg%x2 = xx2; tedg%y2 = yy2 
          call tedgs%add(tedg)
          xx1 = xx2; yy1 = yy2

       end do

    end do ! reading finishes after this completes

    tedgs%maxtag = tot_curves
    tedgs%tol = 2.6d-2 !HARD coded

    ! final check for consistency
    if (size(tedgs%val) .ne. tot_seg ) then
       print *, 'something is wrong in counting the number of added segments! stop'
       stop
    end if

    ! closing the input file
    close(9)

    ! done here
  end subroutine read_segment_file

  ! perform pseudo temporal marching of
  ! a gravitional field about some scatterd nodes (masses)
  ! and finds the steady state field. 
  ! this can be used to smoothen point distribution
  ! in a complex polygon which is very helpful for
  ! approximation of Fekete points and radial basis interpolation.
  ! 
  subroutine gravitional_equilib(x, nbn, tedgs, z, distance, dt, itrmax, conv_epsil, min_dist)
    implicit none
    real*8, dimension(:, :), allocatable :: x !x(dim, npts) -> x(dim, final_npts)
    integer, intent(in) :: nbn
    type(edgs), intent(in) :: tedgs
    real*8, intent(in) :: z ! power of gravitional law
    real*8, intent(in) :: distance
    real*8, intent(inout) :: dt ! might be refined to impose stability
    integer, intent(in) :: itrmax
    real*8, intent(in) :: conv_epsil, min_dist

    ! local vars
    integer :: i, itr, max_size
    real*8, dimension(size(x,1)) :: fors !fx, fy, fz ...
    real*8, dimension(size(x,1)) :: xnew
    real*8 :: tot_displ, tot_displ0, tot_displ00, conv_ratio, dx

    do itr = 1, itrmax
       tot_displ = 0.0d0 ! HARD reset
       i = nbn + 1 ! just one point after last bn point 
       do
          ! compute gravitational force on each point
          call find_grav_force(i, x, distance, z, fors)

111       continue
          ! check see if large and unacceptable displacements might happen
          ! then start over immediately with a refined time step
          dx = sum((fors * dt)**2)
          if ( sqrt(dx) >= (min_dist/2.0d0) ) then
             print *, 'Warning: refining the pseudo time step in' &
                  , ' gravitational point distribution ...' 
             dt = dt / 2.0d0
             goto 111
          end if

          ! accumulate to the total displacement for convergence monit.
          tot_displ = tot_displ + dx

          ! compute possible update coordinates
          xnew = x(:, i) + fors * dt

          ! check if the possible coords are inside
          if ( is_inside(tedgs, xnew(1), xnew(2)) ) then ! is really inside!
             x(:, i) = xnew !OK update it then
          else
             call del_pt(x, i)
             i = i - 1 ! point doesn't count
          end if
          i = i + 1 ! go to the next point 
          max_size = size(x, 2)
          if ( i > max_size ) exit
       end do

       ! convergence monitor
       if ( itr .eq. 1 ) then
          tot_displ00 = tot_displ
       else
          conv_ratio = abs(sqrt(tot_displ) - sqrt(tot_displ0)) / sqrt(tot_displ00) 
          print *, conv_ratio
          if ( conv_ratio <= conv_epsil ) return
       end if
       tot_displ0 = tot_displ

    end do ! pseuo time iterations of grav. pot. eqs

    print *, 'Warning : iterations reached to maximum number' &
         , ' without convergence of gravitational field yielded!'

    ! done here

  contains

    ! computes the gravitional force at point x(:,i)
    ! in the field containing points with equal mass
    ! 
    ! distance = largest -> tightly coupled -> Fekete points
    !          = smallest-> lossly coupled  -> almost equally spaced
    !
    ! pow = 2         ==> Fekete gravity
    ! pow = 3..higher ==> more equally spaced
    !
    subroutine find_grav_force(i, x, distance, pow, fors)
      implicit none
      integer, intent(in) :: i
      real*8, dimension(:, :), intent(in) :: x
      real*8, intent(in) :: distance, pow
      real*8, dimension(:), intent(out) :: fors

      ! local vars
      integer :: j
      real*8 :: l
      real*8, dimension(size(fors)) :: nx

      fors = 0.0d0
      do j = 1, size(x, 2) ! loop over all masses

         if (j .eq. i) cycle

         ! find the distance
         l = sqrt(sum((x(:, j) - x(:, i))**2))

         if (l <= distance) then ! include the force
            nx = (x(:, j) - x(:, i)) / l
            fors = fors - 1.0d0 / (l**pow) * nx
         end if

      end do

      ! done here
    end subroutine find_grav_force

    ! deletes entry "i" in real*8 vector x(:,i)
    subroutine del_pt(x, i)
      implicit none
      real*8, dimension(:, :), allocatable :: x
      integer, intent(in) :: i

      ! local vars
      integer :: dim, leng, n
      real*8, dimension(:, :), allocatable :: tmp

      dim = size(x, 1)
      leng = size(x, 2)
      n =  leng - 1

      allocate(tmp(dim, n))

      if (i > 1) tmp(:, 1:(i-1)) = x(:, 1:(i-1))

      if (i < leng) tmp(:, i:n) = x(:, (i+1):leng) 

      call move_alloc(tmp, x)

      ! safe clean ups
      if (allocated(tmp)) deallocate(tmp)

      ! done here
    end subroutine del_pt

  end subroutine gravitional_equilib

  ! makes an animation to show how
  ! gravitational point distibution
  ! works on polygons 
  subroutine grav_animation(infile, d, unif, tgrav_param)
    implicit none
    character(len = *), intent(in) :: infile
    integer, intent(in) :: d
    logical, intent(in) :: unif
    type(grav_param) :: tgrav_param

    ! local vars
    type(fekete) :: this
    type(edgs) :: tedgs
    integer :: dim, max_anim_seq, itrs, nbn, stat
    character(len = 128) :: jpg_num, jpg_name
    real*8 :: min_dist, max_dist, distance

    ! first read the CAD input file containing edges
    print *, 'starting to make animation sequnce ...'
    call read_segment_file(infile, tedgs)

    ! init.
    dim = 2
    allocate(this%n(dim))
    this%npedg = d+1
    this%nedg = size(tedgs%val)
    this%n(1) = ceiling(dble(this%nedg)/4.0D0 * dble(this%npedg))
    this%n(2) = this%n(1) 
    this%nw = this%n(1) * this%n(2)
    allocate(this%coords(dim))
    allocate(this%w(this%nw))
    allocate(this%fin_coords(dim, this%nw))

    ! store the edge list
    this%tedgs = tedgs

    ! generate master element point distribution
    if ( unif ) then
       call gen_master_custom2d_uniform(this, this%npedg, this%x0, nbn, min_dist, max_dist)
    else ! random
       call gen_master_custom2d(this, this%npedg, this%x0, nbn, min_dist, max_dist)
    end if

    ! set coupling behavior
    if ( tgrav_param%tightly_coupled ) then
       distance = max_dist
    else
       distance = min_dist
    end if

    max_anim_seq = tgrav_param%itrmax
    tgrav_param%itrmax = 1 ! one time step only

    do itrs = 1, max_anim_seq

       ! subroutine gravitional_equilib(x, nbn, tedgs, z, distance, dt, itrmax, conv_epsil)
       if ( itrs .eq. 1 ) then
          call gravitional_equilib(this%x0, nbn, this%tedgs, tgrav_param%z, distance &
               , tgrav_param%dt, 0, tgrav_param%conv_epsil, min_dist)
       else

          call gravitional_equilib(this%x0, nbn, this%tedgs, tgrav_param%z, distance &
               , tgrav_param%dt, tgrav_param%itrmax, tgrav_param%conv_epsil, min_dist)

       end if

       deallocate(this%fin_coords, this%w)
       allocate(this%fin_coords(dim, size(this%x0,2)), this%w(size(this%x0,2)))
       this%fin_coords = this%x0
       this%w(1:nbn) = 1.0d0
       this%w((nbn+1)::1) = 0.0d0
       this%initialized = .true.

       ! choose the right name
       write(jpg_num, '(I5.5)') itrs
       jpg_name = 'seq'//trim(jpg_num)//'.dat'

       ! export coordinates
       call this%export(jpg_name)

       ! perform OS call to python script and create Jpeg image
       stat = system('python viselem.py -f '//jpg_name//' -c -e pdf DAEMON')

       ! show status
       print *, itrs, ' animation sequences were created!'

    end do

    ! 
    print *, 'converting image sequence to animation gif file ...' 
    stat = system('convert -delay 50 -loop 0 *.pdf '//infile//'.gif')

    ! clean up temporary files and arrays
    stat = system('rm seq*.pdf seq*.dat')

    ! OK now
    print *, 'animation was created successfully!'

    ! done here
  end subroutine grav_animation

end module approx_fekete

! program tester
!   use approx_fekete
!   implicit none

!   ! local vars
!   type(fekete) :: tfekete
!   type(edgs) :: tedgs0
!   real*8 :: err
!   integer :: p, i
!   type(grav_param) :: tgrav_param

!   ! set up the gravitational field parameters
!   tgrav_param%z = 2.0d0
!   tgrav_param%dt = 1.0d0
!   tgrav_param%itrmax = 100
!   tgrav_param%conv_epsil = 1.0d-9
! ! tgrav_param%tightly_coupled = .true.
! tgrav_param%no_fekete = .true.
!   ! call grav_animation(infile = './segs/heart.dat', d = 1, unif = .true., tgrav_param = tgrav_param)
!   ! stop

!   print *, '<<<custom2d-input file>>>'
!   call read_segment_file('./segs/cross.dat', tedgs0)
!   i = 0
!   call tfekete%init(d = i, name = 'custom2d', s = 3, echo = .true., tedgs = tedgs0, tgrav_param = tgrav_param)
!   call tfekete%export('custom2d_file.dat')
!   call tfekete%export_custom2D('custom2d_file.dat')
!   call tfekete%interp(15, err, 'interp_file.dat', p, use_RBF = .true.)
!   print *, p , err
!   call tfekete%clean()

!   stop

!   print *, '<<<custom2d>>>'
!   call tedgs0%sample()
!   print *, 'polynomial order <p>, L2 error'
!   do i = 0, 4
!      call tfekete%init(d = i, name = 'custom2d', s = 3, tedgs = tedgs0)
!      call tfekete%export('custom2d_d_3.dat')
!      call tfekete%export_custom2D('custom2d_all.dat')
!      call tfekete%interp(4, err, 'interp.dat', p)
!      print *, p , err
!      call tfekete%clean()
!   end do

!   stop

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
!   call tfekete%init(d = 20, name = 'tri', magnify = 5, s = 6, echo = .true.)
!   call tfekete%export('tri_d_20.dat')
!   call tfekete%clean()

!   print *, '<<< dense tetrahederal >>>'
!   call tfekete%init(d = 4, name = 'tet', magnify = 3, echo = .true.)
!   call tfekete%export('tet_d_4.dat')
!   call tfekete%clean()

!   print *, '<<< hex >>>'
!   call tfekete%init(d = 2, name = 'hex', magnify = 6, s = 3, echo = .true.)
!   call tfekete%export('hex_d_2.dat')
!   call tfekete%clean()

!   print *, '<<< prism >>>'
!   call tfekete%init(d = 4, name = 'prism', magnify = 4, s = 3, echo = .true.)
!   call tfekete%export('prism_d_4.dat')
!   call tfekete%clean()

!   ! done here

! end program tester
