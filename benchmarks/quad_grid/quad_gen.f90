module quad_gen
  use globals
  use spline
  use grid_opt, only : curve

  implicit none

  private

  public :: quad_grid_gen

  integer, parameter :: triangle = 3
  integer, parameter :: quadrilateral = 4
  real(rk) :: tol = 1.0e-10

  type mygrid
    integer :: nelem, nn, nq, nt, n_b_edges
    integer, dimension(:),    allocatable :: elem
    integer, dimension(:),    allocatable :: etype
    integer, dimension(:),    allocatable :: elemidx
    integer, dimension(:),    allocatable :: flag
    integer, dimension(:,:),  allocatable :: b_edges
    integer, dimension(:,:),  allocatable :: be_elem
    integer, dimension(:,:),  allocatable :: e2e
    integer, dimension(:),    allocatable :: b_vals
  end type mygrid

  contains
    subroutine quad_grid_gen(inputfile, elem, elemidx, etype, x, y, nq, nt)
      character(len=*),                      intent(in)     :: inputfile
      integer,  dimension(:),   allocatable, intent(in out) :: elem
      integer,  dimension(:),   allocatable, intent(in out) :: elemidx
      integer,  dimension(:),   allocatable, intent(in out) :: etype
      real(rk), dimension(:),   allocatable, intent(in out) :: x, y
      integer,                               intent(in out) :: nq, nt
      integer, dimension(:), allocatable :: nntoc, ntoc
      type(mygrid) :: g
      type(curve), dimension(:), allocatable :: bc
      integer :: nb, i, j, k

      call read_segment_file(inputfile, bc, nb, nurbs_i)

      call create_grid(g, bc, nb, x, y, g%b_edges, g%b_vals, g%n_b_edges, g%be_elem)

      call get_node_to_quad(g%nn, g%elem, g%etype, g%flag, g%elemidx, nntoc, ntoc)
      call get_e2e(g%elem, g%elemidx, g%etype, nntoc, ntoc, g%e2e)

      allocate(elemidx(size(g%elemidx)))
      allocate(elem(size(g%elem)))
      allocate(etype(size(g%etype)))
      elem(:) = g%elem(:)
      elemidx(:) = g%elemidx(:)
      etype(:) = g%etype(:)
      nq = g%nq
      nt = g%nt

      deallocate(nntoc, ntoc)
    end subroutine quad_grid_gen

    subroutine read_segment_file(inputfile, bc, nb, bt)
      character(len=*), intent(in) :: inputfile
      type(curve), dimension(:), allocatable, intent(in out) :: bc
      integer, intent(in out) :: nb
      integer, optional :: bt
      integer :: i, j, npt, tot_seg, tot_curves, istat, nc, btype
      real(rk) :: xx, yy, zz
      character(len=81) :: mode='natural'

      if(present(bt))then
        btype = bt
      else
        btype = spline_i
      end if
      ! opening the input file
      open ( unit=9, file=inputfile , status = 'old', &
           iostat = istat)
      if ( istat /= 0) then 
         print *, 'fatal: could not open <', inputfile, '> file! exit'
         stop
      end if

      ! first count the total number of boundary curves (connectors)
      ! and the total number of segments
      tot_curves = 0; tot_seg = 0
      do
        read(9,*,iostat=istat)npt
        if(istat > 0) then
           print *, 'fatal error : file <', inputfile &
                ,'> is bad or corrupted. unable to read. exit.'
           stop
        end if
        if( istat < 0 ) exit ! EOF encountered. exit loop
        tot_curves = tot_curves + 1
        tot_seg = tot_seg + npt - 1
        do i = 1, npt
          read(9,*,iostat=istat)xx,yy,zz
        end do
      end do

      print *, 'total boundary curves (connector) : ', tot_curves
      print *, 'total boundary segments : ', tot_seg
      allocate(bc(tot_curves))

      ! now fill grd%bn_curves
      rewind( unit = 9) ! go back to header

      ! read all lines of the input file
      do j = 1, tot_curves
         read(9,*,iostat=istat) npt 

         allocate(bc(j)%x(npt), bc(j)%y(npt), bc(j)%t(npt))

         print *, 'for bn_curves(',j,')'
         do i = 1, npt
            read(9,*,iostat=istat) bc(j)%x(i), bc(j)%y(i), zz
            print *, 'x = ', bc(j)%x(i), ' y = ', bc(j)%y(i)
         end do

      end do ! reading finishes after this completes

      ! closing the input file
      close(9)

      ! fit and store spline
      do j = 1, tot_curves

         nc = size(bc(j)%x)
         if(.not. allocated(bc(j)%Mx))then
           call spline_nurbs_alloc( bc(j)%Mx, bc(j)%My, bc(j)%a, bc(j)%b, bc(j)%c,&
                                      & bc(j)%d, bc(j)%cp, bc(j)%dp, nc, bc(j)%btype)
         end if
         call spline_nurbs_comp(bc(j)%x, bc(j)%y,bc(j)%a, bc(j)%b, bc(j)%c, bc(j)%d &
                            , bc(j)%cp, bc(j)%dp, bc(j)%Mx, bc(j)%My, bc(j)%t, &
                            & mode, bc(j)%btype)

      end do

      nb = tot_curves

      ! done here

    end subroutine read_segment_file

    subroutine create_grid(g, bc, nb, x, y, b_edges, b_vals, nbe, be_elem)
      type(mygrid),                             intent(in out) :: g
      type(curve), dimension(:),                intent(in out) :: bc
      integer,                                  intent(in)     :: nb
      real(rk),    dimension(:),   allocatable, intent(in out) :: x, y
      integer,     dimension(:,:), allocatable, intent(in out) :: b_edges
      integer,     dimension(:),   allocatable, intent(in out) :: b_vals
      integer,     dimension(:,:), allocatable, intent(in out) :: be_elem
      integer,                                  intent(in out) :: nbe
      integer,     dimension(:), allocatable :: loops, bpl
      real(rk) :: xn, xx, yn, yx
      integer :: nl, yes, i, k

      call get_extrema(xn, xx, yn, yx, bc, nb)

      !==========================================!
      !   join boundaries into connected loops
      call get_loops(bc, loops, bpl, nl, nb)
      write(*,*)'got loops'

      !==========================================!
      !  create overlay grid to cover entire domain
      call initialize_grid(bc, g, x, y, loops, bpl, nl, nb, xn, xx, yn, yx)
      write(*,*)'overlay grid initialized'

      !==========================================!
      !      cut out all quads with boundary
      ! intersection and rejoin boundaries to mesh
      call interface_boundary_grid(g, bc, x, y, loops, bpl, nl, nb, b_edges, b_vals, nbe)
      write(*,*)'interface between boundaries and overlay grid created'

      write(*,*)'enter 1 to subdivide bad quads, 0 otherwise'
      read(*,*)yes
      if(yes > 0)then
        call subdivide(x, y, g%elem, g%etype, g%nq, g%nt)
        g%nelem = g%nq + g%nt
      end if
      allocate(g%elemidx(g%nelem))
      k = 1
      do i = 1, g%nelem
        g%elemidx(i) = k
        k = k + g%etype(i)
      end do

      call get_elem_edge(g, b_edges, be_elem, nbe)

      call angle_smooth(g%flag, g%elem, g%etype, g%elemidx, x, y)
      write(*,*)'smoothing completed'

      deallocate(loops, bpl)
    end subroutine create_grid

    subroutine get_elem_edge(g, b_edges, be_elem, nbe)
      type(mygrid), intent(in) :: g
      integer, dimension(:,:), intent(in) :: b_edges
      integer, dimension(:,:), allocatable, intent(in out) :: be_elem
      integer, intent(in out) :: nbe
      integer :: i, j, k, n1, n2, n3, n4

      allocate(be_elem(nbe,2))
      be_elem(:,:) = 0

      do i = 1, g%nelem
        do j = 0, g%etype(i) - 1
          n1 = g%elem(g%elemidx(i) + j)
          n2 = g%elem(g%elemidx(i) + MOD(1 + j, g%etype(i)))

          do k = 1, nbe
            n3 = b_edges(k,1)
            n4 = b_edges(k,2)
            if(n1 .eq. n3 .and. n2 .eq. n4)then
              be_elem(k,1) = i
              be_elem(k,2) = j + 1
            end if
          end do
        end do
      end do

    end subroutine get_elem_edge

    subroutine get_e2e(elem, elemidx, etype, nntoc, ntoc, e2e)
      integer, dimension(:),                intent(in)     :: elem
      integer, dimension(:),                intent(in)     :: elemidx, etype, ntoc, nntoc
      integer, dimension(:,:), allocatable, intent(in out) :: e2e
      integer :: i, j, n1, n2, e1, e2, e11, e22, f

      allocate(e2e(size(etype), 4))

      do e1 = 1, size(etype)
        do f = 0, etype(e1) - 1
          n1 = elem(elemidx(e1) + f)
          n2 = elem(elemidx(e1) + MOD(1 + f, etype(e1)))

          e2 = 0
          do i = nntoc(n1), nntoc(n1 + 1) - 1
            e11 = ntoc(i)
            do j = nntoc(n2), nntoc(n2 + 1) - 1
              e22 = ntoc(j)
              if( abs(e11 - e22) < 1 .and. abs(e11 - e1) > 0)then
                e2 = e11
              end if
            end do
          end do

          e2e(e1, f + 1) = e2

        end do
      end do
    end subroutine get_e2e

    ! distribute np(i) points evenly along boundary(i) 
    subroutine distribute_points(bc, np, bx, by, x, y, bconn, bvals, bpl, loops, l)
      type(curve),  dimension(:), intent(in)     :: bc
      integer,      dimension(:), intent(in)     :: np
      real(rk),     dimension(:), intent(in out) :: bx
      real(rk),     dimension(:), intent(in out) :: by
      real(rk),     dimension(:), intent(in)     :: x, y
      integer,      dimension(:), intent(in)     :: bconn
      integer,      dimension(:), intent(in out) :: bvals
      integer,      dimension(:), intent(in)     :: bpl, loops
      integer, intent(in) :: l
      integer :: b, i, p, prime, nb, lc
      real(rk) :: d, xv, yv
      real(rk) :: t1 = 0.0_rk
      real(rk) :: t2 = 1.0_rk
      character(len=81) :: opt

      opt = 'interp'
      bx(:) = 0.0_rk
      by(:) = 0.0_rk
      nb = size(np)
      prime = 0
      i = 1

      do lc = bpl(l), bpl(l + 1) - 1
        b = loops(lc)
        do p = 1, np(b)
          d = (p - 1.0_rk) / (np(b))
          d = bc(b)%t(1) + d * (bc(b)%t(size(bc(b)%t)) - bc(b)%t(1))
          call spline_nurbs_eval2(d, bc(b)%x, bc(b)%y, bc(b)%a, bc(b)%b, &
                      & bc(b)%c, bc(b)%Mx, bc(b)%My, bc(b)%t, xv, yv, &
                      & opt, bc(b)%btype)
          bvals(i) = b
          bx(i) = xv
          by(i) = yv
          i = i + 1
        end do
      end do
    end subroutine distribute_points

    subroutine create_tmp_points(loops, bpl, l, bc, bx, by, npb)
      integer,       dimension(:), intent(in)     :: loops, bpl, npb
      type(curve), dimension(:), intent(in)     :: bc
      integer,                     intent(in)     :: l
      real(rk),      dimension(:), intent(in out) :: bx, by
      integer :: i, j, k, b
      real(rk) :: d, xv, yv
      character(len=81) :: opt

      opt = 'interp'
      k = 1
      do i = bpl(l), bpl(l + 1) - 1
        b = loops(i)
        do j = 1, npb(b)
          d = (j - 1.0_rk) / (npb(b) - 0.0_rk)
          d = bc(b)%t(1) + d * (bc(b)%t(size(bc(b)%t)) - bc(b)%t(1))
          call spline_nurbs_eval2(d, bc(b)%x, bc(b)%y, bc(b)%a, bc(b)%b, &
                      & bc(b)%c, bc(b)%Mx, bc(b)%My, bc(b)%t, xv, yv, &
                      & opt, bc(b)%btype)
          bx(k) = xv
          by(k) = yv
          k = k + 1
        end do
      end do
    end subroutine create_tmp_points

    subroutine get_closest_index(x, y, bconn, idx, xv, yv, b1, b2)
      real(rk), dimension(:), intent(in)     :: x, y
      integer,  dimension(:), intent(in)     :: bconn
      integer,                intent(in out) :: idx
      real(rk),               intent(in)     :: xv, yv
      integer,                intent(in)     :: b1, b2
      integer :: i, nn, n
      real(rk) :: diff, v

      idx = b1 + 1
      diff = (x(bconn(idx)) - xv)**2 + (y(bconn(idx)) - yv)**2
      nn = size(x)

      do i = b1 + 1, b2 - 1
        n = bconn(i)
        v = (x(n) - xv)**2 + (y(n) - yv)**2
        if(v < diff)then
          diff = v
          idx = i
        end if
      end do
    end subroutine get_closest_index

    subroutine get_extrema(xn, xx, yn, yx, bc, nb)
      real(rk), intent(in out) :: xn, xx, yn, yx
      type(curve), dimension(:), intent(in) :: bc
      integer, intent(in) :: nb
      integer :: b
      xn = MINVAL(bc(1)%x)
      xx = xn
      yn = MINVAL(bc(1)%y)
      yx = yn
      do b = 1, nb
        xn = MIN(xn, MINVAL(bc(b)%x))
        xx = MAX(xx, MAXVAL(bc(b)%x))
        yn = MIN(yn, MINVAL(bc(b)%y))
        yx = MAX(yx, MAXVAL(bc(b)%y))
      end do
    end subroutine get_extrema
  
    subroutine initialize_grid(bc, g, x, y, loops, bpl, nl, nb, gxn, gxx, gyn, gyx)
      type(curve), dimension(:),          intent(in)     :: bc
      type(mygrid),                         intent(in out) :: g
      real(rk),  dimension(:), allocatable, intent(in out) :: x
      real(rk),  dimension(:), allocatable, intent(in out) :: y
      integer,   dimension(:),              intent(in)     :: loops, bpl
      integer,                              intent(in)     :: nl, nb
      real(rk),                             intent(in out) :: gxn, gxx, gyn, gyx
      integer :: nn, nx, ny
      real(rk) :: dx, dy
      integer :: i, j, k, n1, n2, row, col
      real(rk), dimension(:), allocatable :: ddx,ddy
      integer :: check

      call get_rays(gxn, gxx, gyn, gyx, nx, ny, dx, dy)

      allocate(ddx(nx * ny), ddy(nx * ny))
      do i = 1, ny
        do j = 1, nx
          ddx( (i - 1) * nx + j ) = gxn + dx * (j - 1)
          ddy( (i - 1) * nx + j ) = gyn + dy * (i - 1)
        end do
      end do

      call get_spacing(g, loops, bpl, nl, nb, bc, ddx, ddy, nx, ny, gxn, gyn)

      check = 1
      do while(check > 0)
        check = -1
        call preprocess_quads(bc, ddx, ddy, nb, nx, ny, check)
      end do

      nn = nx * ny
      g%nelem = (nx - 1) * (ny - 1)
      allocate( g%elem(4 * g%nelem) )
      ! quads holds the 4 nodes attached to each quad 
      ! stored              4----3
      ! with the            |    |
      ! following           |    |
      ! winding             1----2
      ! when considering local node numbers
      ! e.g. quad i node 3 is accessed as g%quad(4 * (i - 1) + 3)
      k = 0
      do j = 0, ny - 2
        do i = 1, nx - 1
          n1 = i + j * nx
          n2 = i + (j + 1) * nx
          g%elem(4 * k + 1) = n1
          g%elem(4 * k + 2) = n1 + 1
          g%elem(4 * k + 3) = n2 + 1
          g%elem(4 * k + 4) = n2
          k = k + 1
        end do
      end do

      allocate(x(nn))
      allocate(y(nn))

      ! fill in coordinate points x, y
      x(:) = ddx(:)
      y(:) = ddy(:)
      g%nn = nn
    end subroutine initialize_grid

    subroutine get_rays(xmin, xmax, ymin, ymax, nx, ny, dx, dy)
      real(rk),  intent(in out) :: xmin, xmax, ymin, ymax
      integer,   intent(in out) :: nx, ny
      real(rk),  intent(in out) :: dx, dy

      write(*,*)'Enter average spacing in x and y directions <dx dy>'
      read(*,*)dx, dy

      xmin = xmin - 0.5_rk * dx
      ymin = ymin - 0.5_rk * dy
      xmax = xmax + 0.5_rk * dy
      ymax = ymax + 0.5_rk * dy

      nx = MAX(int( (xmax - xmin) / dx + 0.5) + 1, 6)
      ny = MAX(int( (ymax - ymin) / dy + 0.5) + 1, 6)

      dx = (xmax - xmin) / (nx - 1.0_rk)
      dy = (ymax - ymin) / (ny - 1.0_rk)
    end subroutine get_rays

    subroutine get_spacing(g, loops, bpl, nl, nb, bc, ddx, ddy, nx, ny, xmin, ymin)
      type(mygrid),                            intent(in)     :: g
      integer,      dimension(:),              intent(in)     :: loops, bpl
      integer,                                 intent(in)     :: nl, nb
      type(curve),dimension(:),              intent(in)     :: bc
      real(rk),     dimension(:), allocatable, intent(in out) :: ddx, ddy
      integer,                                 intent(in out) :: nx, ny
      real(rk),                                intent(in)     :: xmin, ymin
      integer, dimension(:), allocatable :: divide_x, divide_y
      integer :: i, j, k, n1, n2, outer_loop, nnx, nny
      real(rk) :: xv1, xv2, yv1, yv2
      real(rk) :: dx, dy, xn, yn, xx, yx, g_xx, g_yx, g_xn, g_yn, xv, yv, odx, ody
      real(rk), dimension(:), allocatable :: ndx, ndy

      odx = ddx(2) - ddx(1)
      ody = ddy(nx + 1) - ddy(1)

      allocate(divide_x(nx), divide_y(ny))

      divide_x(:) = 0
      divide_y(:) = 0

      g_xx = MINVAL(bc(1)%x) - 1
      g_yx = MINVAL(bc(1)%y) - 1
      g_xn = MAXVAL(bc(1)%x) + 1
      g_yn = MAXVAL(bc(1)%y) + 1

      outer_loop = -1
      do i = 1, nl
        xn = MINVAL(bc(loops(bpl(i)))%x)
        yn = MINVAL(bc(loops(bpl(i)))%y)
        xx = MAXVAL(bc(loops(bpl(i)))%x)
        yx = MAXVAL(bc(loops(bpl(i)))%y)
        do j = bpl(i), bpl(i + 1) - 1
          xn = MIN(xn, MINVAL(bc(loops(j))%x))
          yn = MIN(yn, MINVAL(bc(loops(j))%y))
          xx = MAX(xx, MAXVAL(bc(loops(j))%x))
          yx = MAX(yx, MAXVAL(bc(loops(j))%y))
        end do

        if( xn > g_xn .and. xx < g_xx .and. yn > g_yn .and. yx < g_yx)then
        else
          outer_loop = i
          g_xx = xx
          g_yx = yx
          g_xn = xn
          g_yn = yn
        end if
      end do

      if(outer_loop > 0)then
        write(*,*)'enter spacing around interior loops'
        read(*,*)dx, dy
        n1 = MAX( int(odx / dx), 1)
        n2 = MAX( int(ody / dy), 1)
        do i = 1, nl
          if(abs(outer_loop - i) > 0)then
            xn = MINVAL(bc(loops(bpl(i)))%x)
            xx = MAXVAL(bc(loops(bpl(i)))%x)
            yn = MINVAL(bc(loops(bpl(i)))%y)
            yx = MAXVAL(bc(loops(bpl(i)))%y)

            do j = bpl(i), bpl(i + 1) - 1
              xn = MIN(xn, MINVAL(bc(loops(j))%x))
              yn = MIN(yn, MINVAL(bc(loops(j))%y))
              xx = MAX(xx, MAXVAL(bc(loops(j))%x))
              yx = MAX(yx, MAXVAL(bc(loops(j))%y))
            end do

            do j = 1, nx - 1
              xv1 = xmin + j * odx
              xv2 = xmin + (j - 1) * odx
              if( xv1 >= (xn - 2 * dx) .and. xv2 <= (xx + 2 * dx) )then
                divide_x(j) = MAX(divide_x(j), n1 - 1)
              end if
            end do
            do j = 1, ny - 1
              yv1 = ymin + j * ody
              yv2 = ymin + (j - 1) * ody
              if( yv1 >= (yn - 2 * dy) .and. yv2 <= (yx + 2 * dy) )then
                divide_y(j) = MAX(divide_y(j), n2 - 1)
              end if
            end do
          end if
        end do
      end if

      call check_enclosed_loops(loops, bpl, nl, bc, xmin, ymin, odx, ody, dx, dy, divide_x, divide_y)

      nnx = nx + SUM(divide_x)
      nny = ny + SUM(divide_y)

      do i = 1, size(divide_x)
        if(divide_x(i) > 0)then
          nnx = nnx - 1
        end if
      end do
      do i = 1, size(divide_y)
        if(divide_y(i) > 0)then
          nny = nny - 1
        end if
      end do

      allocate(ndx(nnx))
      allocate(ndy(nny))

      ndx(1) = ddx(1)
      ndx(nnx) = ddx(nx)
      ndy(1) = ddy(1)
      ndy(nny) = ddy(ny)

      j = 2
      do i = 1, nx - 1
        if(divide_x(i) > 0)then
          do k = 1, divide_x(i)
            ndx(j) = ddx(i) + k * (odx / real(divide_x(i)))
            j = j + 1
          end do
        else
          ndx(j) = ddx(i + 1)
          j = j + 1
        end if
      end do

      j = 2
      do i = 1, ny - 1
        if(divide_y(i) > 0)then
          do k = 1, divide_y(i)
            ndy(j) = ddy(1 + (i - 1) * nx) + k * (ody / real(divide_y(i)))
            j = j + 1
          end do
        else
          ndy(j) = ddy(1 + i * nx)
          j = j + 1
        end if
      end do

      deallocate(ddx, ddy)
      allocate(ddx(nnx * nny))
      allocate(ddy(nnx * nny))

      do i = 0, nny - 1
        ddx( i * nnx + 1 : (i + 1) * nnx) = ndx( 1 : nnx)
      end do

      do i = 1, nny
        do j = 1, nnx
          ddy( (i - 1) * nnx + j) = ndy(i)
        end do
      end do

      nx = nnx
      ny = nny
      deallocate(ndx, ndy)
      deallocate(divide_x, divide_y)
    end subroutine get_spacing

    subroutine check_enclosed_loops(loops, bpl, nl, bc, xmin, ymin, odx, ody, dx, dy, nx, ny)
      integer,       dimension(:), intent(in)     :: loops, bpl
      integer,                     intent(in)     :: nl
      type(curve), dimension(:), intent(in)     :: bc
      real(rk),                    intent(in)     :: xmin, ymin
      real(rk),                    intent(in out) :: odx, ody, dx, dy
      integer,       dimension(:), intent(in out) :: nx, ny
      integer  :: l, lc, b, row, col, nrow, ncol, i, dnx, dny
      real(rk) :: xn, xx, yn, yx, x, y, exn, exx, eyn, eyx
      nrow = size(ny)
      ncol = size(nx)

      do l = 1, nl
        xn = MINVAL(bc(loops(bpl(l)))%x)
        yn = MINVAL(bc(loops(bpl(l)))%y)
        xx = MAXVAL(bc(loops(bpl(l)))%x)
        yx = MAXVAL(bc(loops(bpl(l)))%y)
        dnx = MAX(int((odx + 0.5_rk) / (xx - xn)), int((odx + 0.5_rk) / dx), 1)
        dny = MAX(int((ody + 0.5_rk) / (yx - yn)), int((ody + 0.5_rk) / dy), 1)

        do lc = bpl(l), bpl(l + 1) - 1
          b = loops(lc)
          xn = MIN(xn, MINVAL(bc(b)%x))
          yn = MIN(yn, MINVAL(bc(b)%y))
          xx = MAX(xx, MAXVAL(bc(b)%x))
          yx = MAX(yx, MAXVAL(bc(b)%y))
        end do ! lc : boundaries per loop l

        do i = 1, ncol
          x = xmin + (i - 1) * odx
          if(x < xn)then
            col = i
            exn = x
            exx = xmin + i * odx
          end if
        end do
        do i = 1, nrow
          y = ymin + (i - 1) * ody
          if(y < yn)then
            row = i
            eyn = y
            eyx = ymin + i * ody
          end if
        end do
        if(exn < xn .and. exx > xx .and. eyn < yn .and. eyx > yx)then
          nx(col) = MAX(nx(col), dnx)
          ny(row) = MAX(ny(row), dny)
        end if
      end do ! l : boundary loop count
    end subroutine check_enclosed_loops

    subroutine preprocess_quads(bc, x, y, nb, nx, ny, mycheck)
      type(curve), dimension(:),              intent(in)     :: bc
      real(rk),      dimension(:), allocatable, intent(in out) :: x, y
      integer,                                  intent(in)     :: nb
      integer,                                  intent(in out) :: nx, ny, mycheck
      integer :: n1, n2, p1, p2, c1, c2
      integer :: ect, ne, n_int, b, i, j, iftrue, row, col, nrow, ncol
      real(rk) :: x1, x2, y1, y2, ix1, ix2, ix3, ix4, iy1, iy2, iy3, iy4, xx, yx
      real(rk) :: bx1, bx2, by1, by2
      real(rk), dimension(:), allocatable :: me
      integer,  dimension(:), allocatable :: e

      mycheck = -1

      nx = nx - 1
      ny = ny - 1
      nrow = ny
      ncol = nx
      ne = nx * (ny + 1) + ny * (nx + 1)
      allocate(e(2 * ne))
      ect = 0
      do row = 1, nrow + 1
        do col = 1, ncol
          n1 = (row - 1) * (nx + 1) + col
          n2 = (row - 1) * (nx + 1) + col + 1
          e(2 * ect + 1) = n1
          e(2 * ect + 2) = n2
          ect = ect + 1
        end do
      end do
      do row = 1, nrow
        do col = 1, ncol + 1
          n1 = (row - 1) * (ncol + 1) + col
          n2 =  row      * (ncol + 1) + col
          e(2 * ect + 1) = n1
          e(2 * ect + 2) = n2
          ect = ect + 1
        end do
      end do

      xx = MAXVAL(x)
      yx = MAXVAL(y)
      allocate(me(2 * ne))
      me(:) = xx + 2

      do ect = 1, ne                     ! loop through all edges
        n1 = e(2 * ect - 1)              ! get points of edge
        n2 = e(2 * ect)

        x1 = x(n1); x2 = x(n2)           ! (x1,y1), (x2,y2) are endpoints of edge
        y1 = y(n1); y2 = y(n2)
        n_int = 0
        call find_intersection_nos(x1, y1, x2, y2, bc, nb, ix1, iy1, n_int)

        ! check that no more than 2 sides of any quad is intersected
        call check_max_quad_intersection(nrow, ncol, ne, me, x, y, e, xx + 2)

        if(n_int > 1)then
          mycheck = 1
          me(2 * ect - 1) = ix1
          me(2 * ect)     = iy1
        end if
      end do

      call add_new_vals(x, y, me, e, nrow, ncol, xx)

      nx = ncol
      ny = nrow

      deallocate(e)
      deallocate(me)
    end subroutine preprocess_quads

    subroutine find_intersection_nos(x1, y1, x2, y2, bc, nb, ix, iy, n_int)
      real(rk),                    intent(in)     :: x1, y1, x2, y2
      type(curve), dimension(:), intent(in)     :: bc
      integer,                     intent(in)     :: nb
      real(rk),                    intent(in out) :: ix, iy
      integer,                     intent(in out) :: n_int
      real(rk) :: x3, y3, x4, y4, ix1, iy1, ix2, iy2, ix3, iy3
      integer :: p1, p2, p3, p4, c1, c2, cc1, cc2
      integer :: i, j, k, b, check, npts

      ix1 = MAX(x1, x2) + 2.0_rk
      cc1 = 0
      cc2 = 0
      check = 1
      n_int = 0
      b = 0

      do while(b < nb .and. check > 0)
        b = b + 1
        npts = size(bc(b)%x)
        do i = 1, npts - 1
          if(check > 0)then
            x3 = bc(b)%x(i)
            y3 = bc(b)%y(i)
            x4 = bc(b)%x(i + 1)
            y4 = bc(b)%y(i + 1)

            call line_line_intersect(x1, y1, x2, y2, x3, y3, x4, y4, p1)

            if(p1 > 0)then
              call check_parallel_lines(x1, y1, x2, y2, x3, y3, x4, y4, p4)
              call line_point_intersect(x3, y3, x4, y4, x1, y1, c1)
              call line_point_intersect(x3, y3, x4, y4, x2, y2, c2)

              ! if intersection occurs in the middle of the segment(nonparallel)
              if(c1 < 1 .and. c2 < 1 .and. p4 < 1)then

                call get_intersection_point(ix3, iy3, x1, y1, x2, y2, x3, y3, x4, y4)

                if(ix1 < (x1 + 1) .and. ix2 < (x2 + 1))then
                  if( sqrt( (ix1 - ix3)**2 + (iy1 - iy3)**2) < tol)then
                  else
                    ix2 = ix3
                    iy2 = iy3
                    ix = 0.5_rk * (ix1 + ix2)
                    iy = 0.5_rk * (iy1 + iy2)
                    n_int = 2
                    check = -1
                  end if
                else
                  ix1 = ix3
                  iy1 = iy3
                  n_int = 1
                end if
              ! if segments are parallel and non-b segment fully overlapping
              else if(p4 > 0 .and. c1 < 1 .and. c2 < 1)then
                ix = 0.5_rk * (x3 + x4)
                iy = 0.5_rk * (y3 + y4)
                n_int = 2
                check = -1
              ! if boundary intersects at edge corner
              else if(c1 > 0)then
                cc1 = 1
              ! if boundary intersects at edge corner
              else if(c2 > 0)then
                cc2 = 1
              end if
            end if
          end if
        end do
      end do
    end subroutine find_intersection_nos

    subroutine check_max_quad_intersection(nrow, ncol, nedges, me, x, y, e, xx)
      integer,  intent(in) :: nrow, ncol, nedges
      real(rk), dimension(:), intent(in out) :: me
      integer,  dimension(:), intent(in) :: e
      real(rk), dimension(:), intent(in) :: x, y
      real(rk), intent(in) :: xx
      integer :: e1, e2, e3, e4, check, row, col, vct, hct
      integer :: c1, c2, c3, c4, n1, n2, n3, n4

      hct = 1
      vct = ncol * (nrow + 1) + 1

      do row = 1, nrow
        do col = 1, ncol
          e1 = hct 
          e2 = hct + ncol
          e3 = vct
          e4 = vct + 1

          c1 = 0
          c2 = 0
          c3 = 0
          c4 = 0

          if(me(2 * e1 - 1) < xx)then
            c1 = 1
          end if
          if(me(2 * e2 - 1) < xx)then
            c2 = 1
          end if
          if(me(2 * e3 - 1) < xx)then
            c3 = 1
          end if
          if(me(2 * e4 - 1) < xx)then
            c4 = 1
          end if

          if( (c1 + c2 + c3 + c4) > 2)then
            n1 = e(2 * e1 - 1)
            n2 = e(2 * e1)
            n3 = e(2 * e2 - 1)
            n4 = e(2 * e2)

            me(2 * e1 - 1) = 0.5_rk * (x(n1) + x(n2))
            me(2 * e1    ) = 0.5_rk * (y(n1) + y(n2))
            me(2 * e2 - 1) = 0.5_rk * (x(n2) + x(n3))
            me(2 * e2    ) = 0.5_rk * (y(n2) + y(n3))
            me(2 * e3 - 1) = 0.5_rk * (x(n3) + x(n4))
            me(2 * e3    ) = 0.5_rk * (y(n3) + y(n4))
            me(2 * e4 - 1) = 0.5_rk * (x(n4) + x(n1))
            me(2 * e4    ) = 0.5_rk * (y(n4) + y(n1))
          end if

          hct = hct + 1
          vct = vct + 1
        end do
        vct = vct + 1
      end do
    end subroutine check_max_quad_intersection

    subroutine check_parallel_lines(x1, y1, x2, y2, x3, y3, x4, y4, p)
      real(rk), intent(in) :: x1, y1, x2, y2, x3, y3, x4, y4
      integer, intent(in out) :: p
      real(rk) :: v1, v2, v3, v4, l1, l2, val
      p = -1

      if( (MIN(x1, x2) > MAX(x3, x4)) .or. &
        & (MIN(y1, y2) > MAX(y3, y4)) .or. &
        & (MAX(x1, x2) < MIN(x3, x4)) .or. & 
        & (MAX(y1, y2) < MIN(y3, y4)) )then
      else
        call line_point_intersect(x1, y1, x2, y2, x3, y3, p)
        if(p > 0)then
          call line_point_intersect(x1, y1, x2, y2, x4, y4, p)
        end if
      end if
    end subroutine check_parallel_lines

    subroutine line_line_intersect(bx1, by1, bx2, by2, x1, y1, x2, y2, isct)
      ! checks if two lines intersect and returns value of 0 if not, 1 if they
      ! do (stored in (int) isct)
      real(rk), intent(in) :: bx1, by1, bx2, by2, x1, y1, x2, y2
      integer, intent(in out) :: isct
      real(rk) :: nx1, nx2, ny1, ny2, bnx, bny, v1, v2
      real(rk) :: mxn, mxb, nxn, nxb, myn, myb, nyn, nyb
      real(rk) :: a, b, c, d, e, f, det

      mxn = MAX(x1,  x2);  nxn = MIN(x1,  x2)
      myn = MAX(y1,  y2);  nyn = MIN(y1,  y2)
      mxb = MAX(bx1, bx2); nxb = MIN(bx1, bx2)
      myb = MAX(by1, by2); nyb = MIN(by1, by2)

      isct = 0
      ! if the two lines can be contained in two nonoverlapping boxes
      if( ((mxn + tol) < nxb) .or. ((mxb + tol) < nxn) .or. &
        & ((myn + tol) < nyb) .or. ((myb + tol) < nyn) )then
        isct = 0

      ! otherwise check if they intersect
      else 
        a = x2 - x1
        b = bx1 - bx2
        c = y2 - y1
        d = by1 - by2
        e = x1 - bx1
        f = y1 - by1

        det = (a * d) - (b * c)
        if(abs(det) > 0)then
          v1 = ( (b * f) - (d * e) ) / det
          v2 = ( (c * e) - (a * f) ) / det

          if( (v1 + tol) > 0.0_rk .and. (v2 + tol) > 0.0_rk .and. &
            & v1 < (1.0_rk + tol) .and. v2 < (1.0_rk + tol) )then
            isct = 1
          end if
        else
          isct = 1
        end if
      end if
    end subroutine line_line_intersect

    subroutine line_point_intersection(x1,y1,x2,y2,x3,y3,p)
      real(rk), intent(in)     :: x1,y1,x2,y2,x3,y3
      integer,  intent(in out) :: p
      real(rk) :: xx, yx, xn, yn
      real(rk) :: d

      xx = MAX(x1, x2)
      xn = MIN(x1, x2)
      yx = MAX(y1, y2)
      yn = MIN(y1, y2)
      if( (x3 > (xx + tol)) .or. ((x3 + tol) < xn) .or. &
        & (y3 > (yx + tol)) .or. ((y3 + tol) < yn) )then
        p = -1
      else
        d = x2 * y3 - y2 * x3 + y1 * (x3 - x2) - x1 * (y3 - y2)
        if(abs(d) > tol)then
          p = -1
        else
          p = 1
        end if
      end if
    end subroutine line_point_intersection

    subroutine line_point_intersect(x1, y1, x2, y2, x3, y3, iftrue)
      real(rk), intent(in)     :: x1, y1, x2, y2, x3, y3
      integer,  intent(in out) :: iftrue
      real(rk) :: v1, v2, v3, v4, l1, l2
      iftrue = 0

      if( (x3 + tol) > MAX(x1, x2) .or. (x3 - tol) < MIN(x1, x2) .or. &
        & (y3 + tol) > MAX(y1, y2) .or. (y3 - tol) < MIN(y1, y2) )then
      else
        v1 = x1 - x3
        v2 = y1 - y3
        v3 = x2 - x3
        v4 = y2 - y3
        l1 = sqrt(v1**2 + v2**2)
        l2 = sqrt(v3**2 + v4**2)
        if(l1 > 0)then
          v1 = v1 / l1
          v2 = v2 / l1
        end if
        if(l2 > 0)then
          v3 = v3 / l2
          v4 = v4 / l2
        end if
        if( abs(v1 * v4 - v2 * v3) < tol)then
          iftrue = 1
        end if
      end if
    end subroutine line_point_intersect

    subroutine get_intersection_point(ix1,iy1,x1,y1,x2,y2,x3,y3,x4,y4)
      real(rk), intent(in out) :: ix1, iy1
      real(rk), intent(in) :: x1,y1,x2,y2,x3,y3,x4,y4
      real(rk) :: m1, m2, m3, m4, b1, b2
      real(rk) :: val

      if( abs((y2 - y1) * (x4 - x3) - (y4 - y3) * (x2 - x1)) > tol)then
        m1 = y2 - y1
        m2 = x2 - x1
        m3 = y4 - y3
        m4 = x4 - x3
        b1 = y1 * m2 - x1 * m1
        b2 = y3 * m4 - x3 * m3

        val = m1 * m4 - m2 * m3

        if(abs(val) > 0)then
          ix1 =-(b1 * m4 - b2 * m2) / val
          iy1 = (b2 * m1 - b1 * m3) / val
        end if
      else
        ix1 = 0.5_rk * (x2 + x1)
        iy1 = 0.5_rk * (y2 + y1)
      end if
    end subroutine get_intersection_point

    subroutine add_new_vals(x, y, me, e, nrow, ncol, xx)
      real(rk), dimension(:), allocatable, intent(in out) :: x, y
      real(rk), dimension(:),              intent(in)     :: me
      integer,  dimension(:),              intent(in)     :: e
      integer,                             intent(in out) :: nrow, ncol
      real(rk),                            intent(in)     :: xx
      integer, dimension(:), allocatable :: rtag, ctag
      real(rk), dimension(:), allocatable :: new_x, new_y
      integer :: new_nn, i, j, ect, rct, row, col, cct, nx, ny
      integer :: n1, n2, n3, n4, n5, nn1, nn2, nn3
      nx = ncol
      ny = nrow

      allocate(rtag(nrow + 1))
      allocate(ctag(ncol + 1))

      rtag(:) = 0
      ctag(:) = 0
      ect = 1
      do i = 1, nrow + 1
        do j = 1, ncol
          if(me(2 * ect - 1) < (xx + 1))then
            ctag(j) = 1
          end if
          ect = ect + 1
        end do
     end do
     do i = 1, nrow
        do j = 1, ncol + 1
          if(me(2 * ect - 1) < (xx + 1))then
            rtag(i) = 1
          end if
          ect = ect + 1
        end do
      end do

      new_nn = (nx + SUM(ctag) + 1) * (ny + SUM(rtag) + 1)
      allocate(new_x(new_nn))
      allocate(new_y(new_nn))
      new_x(:) = xx + 1

      nrow = ny
      ncol = nx
      ect = 0
      rct = 0
      do row = 1, nrow + 1
        cct = 0
        do col = 1, ncol
          n1 = e(2 * ect + 1)
          n2 = e(2 * ect + 2)
          nn1 = (row + rct - 1) * (ncol + 1 + SUM(ctag)) + col + cct
          nn2 = (row + rct - 1) * (ncol + 1 + SUM(ctag)) + col + 1 + ctag(col) + cct

          new_x(nn1) = x(n1)
          new_x(nn2) = x(n2)
          new_y(nn1) = y(n1)
          new_y(nn2) = y(n2)
          if(ctag(col) > 0)then
            nn3 = nn1 + 1
            if(me(2 * ect + 1) < (xx + 1))then
              new_x(nn3) = me(2 * ect + 1)
              new_y(nn3) = me(2 * ect + 2)
            else
              new_x(nn3) = ( new_x(nn1) + new_x(nn2) ) * 0.5_rk
              new_y(nn3) = ( new_y(nn1) + new_y(nn2) ) * 0.5_rk
            end if
          end if
          if(ctag(col) > 0)then
            cct = cct + 1
          end if
          ect = ect + 1
        end do
        if(rtag(row) > 0)then
          rct = rct + 1
        end if
      end do
      rct = 0
      do row = 1, nrow
        cct = 0
        do col = 1, ncol + 1
          n1 = e(2 * ect + 1)
          n2 = e(2 * ect + 2)
          nn1 = (row + rct - 1        ) * (ncol + 1 + SUM(ctag)) + col + cct
          nn2 = (row + rct + rtag(row)) * (ncol + 1 + SUM(ctag)) + col + cct
          new_x(nn1) = x(n1)
          new_x(nn2) = x(n2)
          new_y(nn1) = y(n1)
          new_y(nn2) = y(n2)
          if(rtag(row) > 0)then
            nn3 = nn1 + ncol + 1 + SUM(ctag)
            if(me(2 * ect + 1) < (xx + 1))then
              new_x(nn3) = me(2 * ect + 1)
              new_y(nn3) = me(2 * ect + 2)
            else
              new_x(nn3) = ( new_x(nn1) + new_x(nn2) ) * 0.5_rk
              new_y(nn3) = ( new_y(nn1) + new_y(nn2) ) * 0.5_rk
            end if
          end if
          if(ctag(col) > 0)then
            cct = cct + 1
          end if
          ect = ect + 1
        end do 
        if(rtag(row) > 0)then
          rct = rct + 1
        end if
      end do

      nx = ncol + SUM(ctag)
      ny = nrow + SUM(rtag)
      do row = 2, ny
        do col = 2, nx
          n1 = (row - 1) * (nx + 1) + col
          if(new_x(n1) > xx)then
            n2 = (row - 1) * (nx + 1) + col - 1
            n3 = (row - 1) * (nx + 1) + col + 1
            n4 = (row - 2) * (nx + 1) + col
            n5 = (row    ) * (nx + 1) + col
            new_x(n1) = 0.25_rk * (new_x(n2) + new_x(n3) + new_x(n4) + new_x(n5))
            new_y(n1) = 0.25_rk * (new_y(n2) + new_y(n3) + new_y(n4) + new_y(n5))
          end if
        end do
      end do

      nrow = ny + 1
      ncol = nx + 1

      deallocate(x)
      deallocate(y)
      allocate(x(new_nn))
      allocate(y(new_nn))
      x(:) = new_x(:)
      y(:) = new_y(:)
      deallocate(new_x)
      deallocate(new_y)
      deallocate(rtag)
      deallocate(ctag)
    end subroutine add_new_vals

    subroutine add_color(pcolors, p, nbpl, color)
      integer, dimension(:), intent(in out) :: pcolors
      integer,               intent(in)     :: p, nbpl
      integer,               intent(in)     :: color

      integer :: i, check
      check = 1
      if(p < 1)then
        write(*,*)'Error in routine add_color: array index out of bounds'
        stop
      else
        do i = 1, nbpl
          if( (abs(pcolors(nbpl * (p - 1) + 1)) < 1) .or. &
            & (abs(pcolors(nbpl * (p - 1) + i) - color) < 1) )then
            check = -1
          end if
          if( (pcolors(nbpl * (p - 1) + i) < 0) .and. (check > 0) )then
            pcolors(nbpl * (p - 1) + i) = color
            check = -1
          end if
        end do
      end if
    end subroutine add_color

    subroutine check_intersection(x, y, bx1, by1, bx2, by2, iftrue, rotate)
      ! checks the type (if any) of intersection between edge with endpoints
      ! (bx1, by1), (bx2, by2) and the quad whose corners are stored in x(4), y(4)
      real(rk), dimension(:), intent(in) :: x
      real(rk), dimension(:), intent(in) :: y
      real(rk),               intent(in) :: bx1, by1, bx2, by2
      integer,                intent(in out) :: iftrue
      integer,                intent(in out) :: rotate
      real(rk) :: bnx, bny
      real(rk) :: nx1, ny1
      real(rk) :: x1, x2, x3, x4, y1, y2, y3, y4
      integer :: i, j, isct1, isct2, isct3, isct4
      integer :: p1, p2, p3, p4, check, p11, p12, p13, p14

      iftrue = -1 ! return value: > 0 indicates intersection, < 0 not.

      x1 = x(1); x2 = x(2); x3 = x(3); x4 = x(4)
      y1 = y(1); y2 = y(2); y3 = y(3); y4 = y(4)

      bnx = by2 - by1
      bny = bx1 - bx2

      ! find out if the element either intersects with
      ! the boundary segment or else is entirely outside of the boundary

      call line_line_intersect(bx1, by1, bx2, by2, x1, y1, x2, y2, isct1)
      call line_line_intersect(bx1, by1, bx2, by2, x2, y2, x3, y3, isct2)
      call line_line_intersect(bx1, by1, bx2, by2, x3, y3, x4, y4, isct3)
      call line_line_intersect(bx1, by1, bx2, by2, x4, y4, x1, y1, isct4)

      if( (isct1 + isct2 + isct3 + isct4) > 0)then

        if( (isct1 + isct2 + isct3 + isct4) < 2)then
          iftrue = 4

          if(isct1 > 0)then
            nx1 = y2 - y1
            if(abs(nx1) > 0)then
              nx1 = nx1 / (x1 - x2)
            end if
            ny1 = bnx
            if(abs(bnx) > 0)then
              ny1 = ny1 / bny
            end if
            if(abs(nx1 - ny1) < tol)then
              iftrue = -1
            else
              nx1 = x1 - x2
              ny1 = y1 - y2
              if((bnx * nx1 + bny * ny1) > 0)then
                rotate = 1
              else
                rotate = 6
              end if
            end if
          else if(isct2 > 0)then
            nx1 = x2 - x3
            ny1 = y2 - y3
            if((bnx * nx1 + bny * ny1) > 0)then
              rotate = 2
            else
              rotate = 7
            end if
          else if(isct3 > 0)then
            nx1 = x3 - x4
            ny1 = y3 - y4
            if((bnx * nx1 + bny * ny1) > 0)then
              rotate = 3
            else
              rotate = 8
            end if
          else if(isct4 > 0)then
            nx1 = x4 - x1
            ny1 = y4 - y1
            if((bnx * nx1 + bny * ny1) > 0)then
              rotate = 4
            else
              rotate = 5
            end if
          end if
        else if( (isct1 + isct2 + isct3 + isct4) < 3 )then
          call line_point_intersection(bx1, by1, bx2, by2, x1, y1, p1)
          call line_point_intersection(bx1, by1, bx2, by2, x2, y2, p2)
          call line_point_intersection(bx1, by1, bx2, by2, x3, y3, p3)
          call line_point_intersection(bx1, by1, bx2, by2, x4, y4, p4)

          !if the boundary edge does not intersect at a corner
          if( (p1 < 1) .and. (p2 < 1) .and. (p3 < 1) .and. (p4 < 1) )then
            if( (isct1 + isct2) > 1 )then
                                         ! 4---3/
              nx1 = x2 - x4              ! |   /
              ny1 = y2 - y4              ! 1--/2
                                         !   /
              if( (bnx * nx1 + bny * ny1) > 0)then
                iftrue = 1
                rotate = 3
              else
                iftrue = 3
                rotate = 3
              end if
            else if( (isct2 + isct3) > 1 )then
                                         !   \
              nx1 = x3 - x1              ! 4--\3
              ny1 = y3 - y1              ! |   \
                                         ! 1---2\
              if( (bnx * nx1 + bny * ny1) > 0)then
                iftrue = 1
                rotate = 0
              else
                iftrue = 3
                rotate = 0
              end if
            else if( (isct3 + isct4) > 1 )then
                                         !   /
              nx1 = x4 - x2              ! 4/--3
              ny1 = y4 - y2              ! /   |
                                         !/1---2
              if( (bnx * nx1 + bny * ny1) > 0)then
                iftrue = 1
                rotate = 1
              else
                iftrue = 3
                rotate = 1
              end if
            else if( (isct4 + isct1) > 1 )then
                                         !\4---3
              nx1 = x1 - x3              ! \   |
              ny1 = y1 - y3              ! 1\--2
                                         !   \
              if( (bnx * nx1 + bny * ny1) > 0)then
                iftrue = 1
                rotate = 2
              else
                iftrue = 3
                rotate = 2
              end if
            else if( (isct1 + isct3) > 1)then
                                        ! 4-|-3 
              nx1 = x2 - x1             ! | | | 
              ny1 = y2 - y1             ! 1-|-2 
                                        !   |
              if( (bnx * nx1 + bny * ny1) > 0)then
                iftrue = 2
                rotate = 0
              else
                iftrue = 2
                rotate = 2
              end if

            else if( (isct2 + isct4) > 1)then
                                      ! 4---3
              nx1 = x3 - x2           !-|---|-
              ny1 = y3 - y2           ! 1---2
                                      !
              if( (bnx * nx1 + bny * ny1) > 0)then
                iftrue = 2
                rotate = 1
              else
                iftrue = 2
                rotate = 3
              end if
            end if
          else ! if it's one corner that the line has intersected
               ! (counted twice because each edge does intersect with it)
            p11 = 0
            if(p1 > 0)then
              nx1 = x2 - x1
              ny1 = y2 - y1
              if( (bnx * nx1 + bny * ny1) > 0)then
                nx1 = x4 - x1
                ny1 = y4 - y1
                if( (bnx * nx1 + bny * ny1) > 0)then
                  iftrue = -1
                else
                  iftrue = 4
                  rotate = 5 + p11 * 8
                end if
              else
                nx1 = x4 - x1
                ny1 = y4 - y1
                if( (bnx * nx1 + bny * ny1) > 0)then
                  iftrue = 4
                  rotate = 1 + p11 * 8
                else
                  iftrue = 1
                  rotate = 2
                end if
              end if
            else if(p2 > 0)then
              nx1 = x3 - x2
              ny1 = y3 - y2
              if( (bnx * nx1 + bny * ny1) > 0)then
                nx1 = x1 - x2
                ny1 = y1 - y2
                if( (bnx * nx1 + bny * ny1) > 0)then
                  iftrue = -1
                else
                  iftrue = 4
                  rotate = 6 + p11 * 8
                end if
              else
                nx1 = x1 - x2
                ny1 = y1 - y2
                if( (bnx * nx1 + bny * ny1) > 0)then
                  iftrue = 4
                  rotate = 2 + p11 * 8
                else
                  iftrue = 1
                  rotate = 3
                end if
              end if
            else if(p3 > 0)then
              nx1 = x4 - x3
              ny1 = y4 - y3
              if( (bnx * nx1 + bny * ny1) > 0)then
                nx1 = x2 - x3
                ny1 = y2 - y3
                if( (bnx * nx1 + bny * ny1) > 0)then
                  iftrue = -1
                else
                  iftrue = 4
                  rotate = 7 + p11 * 8
                end if
              else
                nx1 = x2 - x3
                ny1 = y2 - y3
                if( (bnx * nx1 + bny * ny1) > 0)then
                  iftrue = 4
                  rotate = 3 + p11 * 8
                else
                  iftrue = 1
                  rotate = 0
                end if
              end if
            else if(p4 > 0)then
              nx1 = x1 - x4
              ny1 = y1 - y4
              if( (bnx * nx1 + bny * ny1) > 0)then
                nx1 = x3 - x4
                ny1 = y3 - y4
                if( (bnx * nx1 + bny * ny1) > 0)then
                  iftrue = -1
                else
                  iftrue = 4
                  rotate = 8 + p11 * 8
                end if
              else
                nx1 = x3 - x4
                ny1 = y3 - y4
                if( (bnx * nx1 + bny * ny1) > 0)then
                  iftrue = 4
                  rotate = 4 + p11 * 8
                else
                  iftrue = 1
                  rotate = 1
                end if
              end if
            end if
          end if
        else if( (isct1 + isct2 + isct3 + isct4) < 4)then

          call line_point_intersection(bx1, by1, bx2, by2, x1, y1, p1)
          call line_point_intersection(bx1, by1, bx2, by2, x2, y2, p2)
          call line_point_intersection(bx1, by1, bx2, by2, x3, y3, p3)
          call line_point_intersection(bx1, by1, bx2, by2, x4, y4, p4)

          if( (p1 > 0) .and. (isct2 > 0) )then
            nx1 = x4 - x2
            ny1 = y4 - y2
            if( (bnx * nx1 + bny * ny1) > 0 )then
              iftrue = 3
              rotate = 3
            else
              iftrue = 2
              rotate = 3
            end if
          else if( (p1 > 0) .and. (isct3 > 0) )then
            nx1 = x4 - x2
            ny1 = y4 - y2
            if( (bnx * nx1 + bny * ny1) > 0 )then
              iftrue = 2
              rotate = 2
            else
              iftrue = 3
              rotate = 1
            end if
          else if( (p2 > 0) .and. (isct3 > 0) )then
            nx1 = x1 - x3
            ny1 = y1 - y3
            if( (bnx * nx1 + bny * ny1) > 0 )then
              iftrue = 3
              rotate = 0
            else
              iftrue = 2
              rotate = 0
            end if
          else if( (p2 > 0) .and. (isct4 > 0) )then
            nx1 = x1 - x3
            ny1 = y1 - y3
            if( (bnx * nx1 + bny * ny1) > 0 )then
              iftrue = 2
              rotate = 3
            else
              iftrue = 3
              rotate = 2
            end if
          else if( (p3 > 0) .and. (isct4 > 0) )then
            nx1 = x2 - x4
            ny1 = y2 - y4
            if( (bnx * nx1 + bny * ny1) > 0)then
              iftrue = 3
              rotate = 1
            else
              iftrue = 2
              rotate = 1
            end if
          else if( (p3 > 0) .and. (isct1 > 0) )then
            nx1 = x2 - x4
            ny1 = y2 - y4
            if( (bnx * nx1 + bny * ny1) > 0)then
              iftrue = 2
              rotate = 0
            else
              iftrue = 3
              rotate = 3
            end if
          else if( (p4 > 0) .and. (isct1 > 0) )then
            nx1 = x3 - x1
            ny1 = y3 - y1
            if( (bnx * nx1 + bny * ny1) > 0)then
              iftrue = 3
              rotate = 2
            else
              iftrue = 2
              rotate = 2
            end if
          else if( (p4 > 0) .and. (isct2 > 0) )then
            nx1 = x3 - x1
            ny1 = y3 - y1
            if( (bnx * nx1 + bny * ny1) > 0)then
              iftrue = 2
              rotate = 1
            else
              iftrue = 3
              rotate = 0
            end if
          else
            iftrue = -1
            write(*,*)'Error in routine check_intersection: ', &
                    & 'invalid number of edge intersections selected'
          end if
        else if( (isct1 + isct2 + isct3 + isct4) < 5)then

          call line_point_intersection(bx1, by1, bx2, by2, x1, y1, p1)
          call line_point_intersection(bx1, by1, bx2, by2, x2, y2, p2)
          call line_point_intersection(bx1, by1, bx2, by2, x3, y3, p3)
          call line_point_intersection(bx1, by1, bx2, by2, x4, y4, p4)

          if( (p1 > 0) .and. (p3 > 0) )then
            iftrue = 3
            nx1 = x2 - x4
            ny1 = y2 - y4
            if( (bnx * nx1 + bny * ny1) > 0)then
              rotate = 1
            else
              rotate = 3
            end if
          else if( (p2 > 0) .and. (p4 > 0) )then
            iftrue = 3
            nx1 = x1 - x3
            ny1 = y1 - y3
            if( (bnx * nx1 + bny * ny1) > 0)then
              rotate = 0
            else 
              rotate = 2
            end if
          else if( (p1 > 0) .and. (p2 > 0) )then
            nx1 = y2 - y1
            ny1 = x1 - x2
            if( (bnx * nx1 + bny * ny1) > 0)then
              iftrue = 2
              rotate = 3
            else
              iftrue = -1
            end if
          else if( (p1 > 0) .and. (p4 > 0) )then
            nx1 = y1 - y4
            ny1 = x4 - x1
            if( (bnx * nx1 + bny * ny1) > 0)then
              iftrue = 2
              rotate = 2
            else
              iftrue = -1
            end if
          else if( (p2 > 0) .and. (p3 > 0) )then
            nx1 = y3 - y2
            ny1 = x2 - x3
            if( (bnx * nx1 + bny * ny1) > 0)then
              iftrue = 2
              rotate = 0
            else
              iftrue = -1
            end if
          else if( (p3 > 0) .and. (p4 > 0) )then
            nx1 = y4 - y3
            ny1 = x3 - x4
            if( (bnx * nx1 + bny * ny1) > 0)then
              iftrue = 2
              rotate = 1
            else
              iftrue = -1
            end if
          else
            write(*,*)'Error in routine check_intersection: ', &
                    & 'invalid number of edge intersections selected'
          end if
        end if
      end if
    end subroutine check_intersection

    subroutine angle_smooth(flag, elem, etype, elemidx, x, y)
      integer,  dimension(:), intent(in)     :: flag, elem, etype, elemidx
      real(rk), dimension(:), intent(in out) :: x, y
      integer :: n1, n2, n3, n4, i, j, k, qn, nn, nelem
      real(rk) :: vx_1, vx_2, vx_3, vy_1, vy_2, vy_3, alpha1, alpha2, beta
      real(rk), dimension(:), allocatable :: newx, newy
      integer,  dimension(:), allocatable :: ntoc, nntoc
      integer,  dimension(:), allocatable :: nton, nnton
      integer :: iteration
      character(len=1024) :: fn

      nn = size(x)
      if( abs(size(x) - size(y) ) > 0 )then
        write(*,*)'Error in routine angle_smooth: x,y arrays of different size'
        stop
      end if
      nelem = size(etype)

      call get_node_to_quad(nn, elem, etype, flag, elemidx, nntoc, ntoc)

      call get_node_to_node(elem, etype, elemidx, nntoc, ntoc, nnton, nton, nn)

      allocate(newx(nn), newy(nn))

      do iteration = 1, 500
        newx(:) = 0.0_rk
        newy(:) = 0.0_rk

        do i = 1, nn
          if(flag(i) < 0)then
            newx(i) = 0.0_rk
            newy(i) = 0.0_rk
            do j = nnton(i), nnton(i + 1) - 1
              newx(i) = newx(i) + x(nton(j))
              newy(i) = newy(i) + y(nton(j))
            end do
            beta = real(nnton(i + 1) - nnton(i))
            x(i) = newx(i) / beta
            y(i) = newy(i) / beta
          end if
        end do
      end do
      deallocate(newx, newy)
    end subroutine angle_smooth

    subroutine get_node_to_node(elem, etype, elemidx, nntoc, ntoc, nnton, nton, nn)
      integer, dimension(:), intent(in) :: elem, etype, elemidx, nntoc, ntoc
      integer, dimension(:), allocatable, intent(in out) :: nnton, nton
      integer, intent(in) :: nn
      integer :: i, j, nelem, k, n1, n2, idx

      allocate(nnton(nn + 1))
      do i = 1, nn
        j = nntoc(i + 1) - nntoc(i)
        nnton(i + 1) = 2 * j
      end do
      nnton(1) = 1
      do i = 2, nn + 1
        nnton(i) = nnton(i) + nnton(i - 1)
      end do
      allocate(nton(nnton(nn + 1)))

      do i = 1, nn
        do j = nntoc(i), nntoc(i + 1) - 1
          nelem = ntoc(j)
          do k = 1, etype(nelem)
            if( abs( elem(elemidx(nelem) + k - 1) - i ) < 1 )then
              idx = k
            end if
          end do
          n1 = elem(elemidx(nelem) + 1 + MOD(idx, etype(nelem)) - 1)
          n2 = elem(elemidx(nelem) + 1 + MOD(idx + etype(nelem) - 2, etype(nelem)) - 1)
          k = j - nntoc(i)
          nton( nnton(i) + 2 * k    ) = n1
          nton( nnton(i) + 2 * k + 1) = n2
        end do
      end do
    end subroutine get_node_to_node

    subroutine get_node_to_quad(nn, elem, etype, flag, elemidx, nntoc, ntoc)
      integer,                            intent(in)     :: nn
      integer, dimension(:),              intent(in)     :: elem, etype, flag, elemidx
      integer, dimension(:), allocatable, intent(in out) :: nntoc
      integer, dimension(:), allocatable, intent(in out) :: ntoc
      integer :: i, j, k, e, e1, e2, idx, eidx
      integer :: check, n, n1, n2, n3, nelem, idx1, idx2

      nelem = size(etype)
      allocate(nntoc(nn + 2))
      nntoc(:) = 0
      k = 0
      do i = 1, nelem
        do j = 1, etype(i)
          n = elem(k + j)
          nntoc(n + 1) = nntoc(n + 1) + 1
        end do
        k = k + etype(i)
      end do
      nntoc(1) = 1
      do i = 2, nn + 2
        nntoc(i) = nntoc(i) + nntoc(i - 1)
      end do

      allocate(ntoc(nntoc(nn + 2)))
      k = 0
      do i = 1, nelem
        do j = 1, etype(i)
          n = elem(k + j)
          idx = nntoc(n)
          ntoc(idx) = i
          nntoc(n) = idx + 1
        end do
        k = k + etype(i)
      end do
      do i = nn + 1, 2, -1
        nntoc(i) = nntoc(i - 1)
      end do
      nntoc(1) = 1

      ! now order the elem connected to each node to be wound
      ! in a counterclockwise direction
      do i = 1, nn
        ! if node is connected to more than one element
        if( (nntoc(i + 1) - nntoc(i)) > 1)then

          ! if it's a boundary node(neighboring elems will be
          ! wound beginning with a specific elem)
          if(flag(i) > 0)then
            ! first element in winding will have at least
            ! two boundary nodes with the first one i
            do j = nntoc(i), nntoc(i + 1) - 1
              e = ntoc(j)
              do k = 0, etype(e) - 1
                n1 = elem(elemidx(e) + k)
                n2 = elem(elemidx(e) + MOD(k + 1, etype(e)))
                n3 = elem(elemidx(e) + MOD(k + etype(e) - 1, etype(e)))
                if(abs(n1 - i) < 1 .and. flag(n2) > 0 .and. flag(n3) < 0)then
                  eidx = j
                  idx = k
                end if
              end do
            end do

            e1 = ntoc(eidx)
            n1 = elem(elemidx(e1) + MOD(1 + idx + etype(e1) - 2, etype(e1)))
            n2 = elem(elemidx(e1) + idx)

            k = ntoc(nntoc(i))
            ntoc(nntoc(i)) = e1
            ntoc(eidx) = k
            idx = nntoc(i) + 1

            do e = nntoc(i), nntoc(i + 1) - 1
              check = 1
              idx1 = elemidx(e1)

              do j = idx, nntoc(i + 1) - 1
                e2 = ntoc(j)
                idx2 = elemidx(e2)
                do k = 0, etype(e2) - 1
                  if( abs(n2 - elem(idx2 + k)) < 1 .and. &
                    & abs(n1 - elem(idx2 + MOD(1 + k, etype(e2)))) < 1)then
                    eidx = j
                    check = -1
                    n = k
                  end if
                end do
              end do

              if(check < 0)then
                j = ntoc(eidx)
                ntoc(eidx) = ntoc(idx)
                ntoc(idx) = j
                e1 = ntoc(idx)
                idx = idx + 1
                n1 = elem(elemidx(e1) + MOD(1 + n + etype(e1) - 2, etype(e1)))
                n2 = elem(elemidx(e1) + n)
              end if
            end do
          else
            e = ntoc(nntoc(i))
            idx = nntoc(i) + 1
            do e2 = nntoc(i), nntoc(i + 1) - 1

              idx1 = elemidx(e)
              check = 1

              ! find las two nodes of element e(will be shared
              ! with the next element in winding)
              do j = 0, etype(e) - 1
                if(abs(elem(idx1 + j) - i) < 1)then
                  n1 = elem(idx1 + MOD(1 + j + etype(e) - 2, etype(e)))
                  n2 = elem(idx1 + j)
                end if
              end do

              ! loop through all elements surrounding node i
              ! to find the one next to element e
              do j = idx, nntoc(i + 1) - 1
                e1 = ntoc(j)
                idx2 = elemidx(e1)
                do k = 1, etype(e1)
                  if( abs(n2 - elem(idx2 + k - 1)) < 1 .and. &
                    & abs(n1 - elem(idx2 + MOD(k, etype(e1)))) < 1)then
                    eidx = j
                    check = -1
                  end if
                end do
              end do
              ! next element (stored at eidx) swap with 
              ! element stored in space idx and increment count
              if(check < 0)then
                j = ntoc(eidx)
                ntoc(eidx) = ntoc(idx)
                ntoc(idx) = j
                e = j
                idx = idx + 1
              end if
            end do
          end if
        else
          ! do nothing. Node only connected to one element
        end if
      end do
!!FOR DEBUGGING: prints out node together with all elements connected to it
!     do i = 1, nn
!       write(*,*)'node: ',i
!       do j = nntoc(i), nntoc(i + 1) - 1
!         write(*,*),'...',ntoc(j)
!       end do
!     end do
    end subroutine get_node_to_quad

    subroutine get_loops(bc, loops, bpl, nl, nb)
      type(curve),  dimension(:),              intent(in)     :: bc
      integer,        dimension(:), allocatable, intent(in out) :: loops
      integer,        dimension(:), allocatable, intent(in out) :: bpl
      integer,                                   intent(in out) :: nl
      integer,                                   intent(in)     :: nb
      integer, dimension(:), allocatable :: tag
      integer :: b, i, j, k, lcount, check, nlc
      real(rk) :: d, x1, y1, x2, y2, x3, y3

      lcount = 0
      nl = nb
      do i = 1, nb
        check = 1
        x1 = bc(i)%x(1)
        y1 = bc(i)%y(1)
        do j = 1, nb
          if( (check > 0) .and. abs(j - i) > 0)then
            x2 = bc(j)%x(size(bc(j)%x))
            y2 = bc(j)%y(size(bc(j)%x))
            d = sqrt( (x1 - x2)**2 + (y1 - y2)**2 )
            if(d < tol)then
              check = 0
              lcount = lcount + 1
            end if
          end if
        end do
      end do

      nl = nl - (lcount + 1) / 2
      allocate(bpl(nl + 1))
      allocate(loops(nb))
      allocate(tag(nb))
      tag(:) = -1
      loops(:) = -1
      bpl(:) = 1

      nlc = 1
      lcount = 1
      do i = 1, nl
        check = -1

        k = 1
        do j = 1, nb
          if(tag(j) < 0)then
            k = -1
          end if
        end do
        if(k > 0)then
        else
          nlc = i
          do j = 1, nb
            if(tag(j) > 0 .or. check > 0)then
            else if(check < 0)then
              check = 1
              b = j
            end if
          end do
          loops(lcount) = b

          lcount = lcount + 1
          bpl(i + 1) = lcount
          tag(b) = 1

          x1 = bc(b)%x(1)
          y1 = bc(b)%y(1)
          x2 = bc(b)%x(size(bc(b)%x))
          y2 = bc(b)%y(size(bc(b)%x))

          d = sqrt( (x1 - x2)**2 + (y1 - y2)**2 )
          if(d < tol)then
          else
            check = -1
            do while(check < 0)
              do j = 1, nb
              ! loop over boundaries

                if(tag(j) < 0 .and. check < 0)then 
                ! if boundary j has not been included in a loop yet, and 
                ! if the loop is not closed
                  x3 = bc(j)%x(1)
                  y3 = bc(j)%y(1)

                  d = sqrt( (x2 - x3)**2 + (y2 - y3)**2 )
                  if(d < tol)then
                  ! if this boundary connects to the last
                    loops(lcount) = j
                    lcount = lcount + 1
                    bpl(i + 1) = lcount
                    tag(j) = 1
                    x2 = bc(j)%x(size(bc(j)%x))
                    y2 = bc(j)%y(size(bc(j)%x))
                    x3 = bc(j)%x(size(bc(j)%x))
                    y3 = bc(j)%y(size(bc(j)%x))
                    d = sqrt( (x1 - x3)**2 + (y1 - y3)**2 )
                    if(d < tol)then
                      check = 1
                    end if
                  end if
                end if
              end do
            end do
          end if
        end if
      end do
      nl = nlc
      deallocate(tag)
    end subroutine get_loops

    ! remove all elements in vicinity of boundary loops from overlay mesh
    ! and create a mesh interface between boundary and overlay
    subroutine interface_boundary_grid(g, bc, x, y, loops, bpl, nl, nb, b_edges, be_vals, nbe)
      type(mygrid),                              intent(in out) :: g
      type(curve),  dimension(:),                intent(in out) :: bc
      real(rk),     dimension(:),   allocatable, intent(in out) :: x, y
      integer,      dimension(:),                intent(in)     :: loops, bpl
      integer,                                   intent(in)     :: nl, nb
      integer,      dimension(:,:), allocatable, intent(in out) :: b_edges
      integer,      dimension(:),   allocatable, intent(in out) :: be_vals
      integer,                                   intent(in out) :: nbe
      integer,      dimension(:), allocatable :: bconn
      integer,      dimension(:), allocatable :: npb
      integer,      dimension(:), allocatable :: nq
      real(rk),     dimension(:), allocatable :: bx, by
      integer :: i, j, k, l, n, nnq, nnp, nbpl
      integer :: n1, n2, n3, n4, q
      integer, dimension(:,:), allocatable :: tmp
      integer, dimension(:), allocatable :: pcolor, qcolor
      integer :: be_ct

      allocate(b_edges(1,2))
      allocate(tmp(1,2))
      allocate(be_vals(1))
      be_ct = 0

      nbpl = bpl(2) - bpl(1)
      allocate(npb(nb), pcolor(nbpl * g%nn), qcolor(g%nelem))

      do l = 1, nl ! loop through fully connected boundary loops
        nbpl = bpl(l + 1) - bpl(l)

        ! color elements and nodes based on their intersection with boundaries
        call color_elem(l, loops, bpl, bc, g, x, y, nbpl, pcolor, qcolor)

        ! remove the quads that have been colored for deleting
        call delete_external_elem(g, nbpl, pcolor, qcolor)

        ! npb(b) contains the number of points to be distributed along 
        ! boundary curve b (also the number of grid nodes to connect to)
        call get_npb(g%nn, pcolor, nbpl, npb)

        ! bconn contains an ordered list of the nodes surrounding the 
        ! boundary loop l wound clockwise
        nnp = SUM(npb)
        if(nnp < 1)then
          write(*,*)'Error in routine interface_boundary_grid. all nodes deleted leaving empty grid'
          stop
        end if
        allocate(bconn(nnp))
        call order_b_conn(bconn, npb, g, l, bpl, loops, nbpl, pcolor, qcolor)
        nnp = g%nn + SUM(npb)
        nnq = SUM(npb)

        deallocate(b_edges)
        allocate(b_edges(be_ct + SUM(npb), 2))
        b_edges(:be_ct,:2) = tmp(:be_ct,:2)

        ! renumber arrays after deleting external elem and nodes
        call renumber_arrays(g, bconn, x, y, nbpl, b_edges, be_ct, pcolor, qcolor)

        allocate(bx(nnp),  by(nnp))
        allocate(nq(4 * nnq + 4 * g%nelem))

        ! create the new elem that connect the grid with boundary loop l
        call connect_interior_boundary(g, bc, x, y, bx, by, npb, &
                            & bconn, nq, l, bpl, loops, nbpl, &
                            & b_edges, be_ct, be_vals, pcolor)

        deallocate(bx)
        deallocate(by)
        deallocate(bconn)
        deallocate(nq)
        deallocate(qcolor)
        deallocate(pcolor)
        allocate(qcolor(g%nelem))
        be_ct = be_ct + sum(npb)
        if(l < nl)then
          nbpl = bpl(l + 2) - bpl(l + 1)
          deallocate(tmp)
          allocate(tmp(be_ct, 2))
          tmp(:be_ct,:2) = b_edges(:be_ct,:2)
        end if
        allocate(pcolor(nbpl * g%nn))
      end do

      allocate(g%flag(g%nn), g%etype(g%nelem))
      g%flag(:) = -1
      g%etype(:) = quadrilateral

      do i = 1, be_ct
        g%flag(b_edges(i, 1)) = 1
        g%flag(b_edges(i, 2)) = 1
      end do
      g%nq = g%nelem
      g%nt = 0
      nbe = be_ct
      deallocate(tmp)
    end subroutine interface_boundary_grid

    ! flag the elements (and corresponding nodes) 
    !     that intersect boundary loop l
    subroutine color_elem(l, loops, bpl, bc, g, x, y, nbpl, pcolor, qcolor)
      integer,                    intent(in)     :: l, nbpl
      integer,      dimension(:), intent(in)     :: loops, bpl
      type(curve),dimension(:), intent(in)     :: bc
      type(mygrid),               intent(in out) :: g
      real(rk),     dimension(:), intent(in)     :: x, y
      integer,      dimension(:), intent(in out) :: pcolor, qcolor
      integer :: b, i, j, q, n1, n2, n3, n4, nn1, nn2, nn3, nn4
      integer :: side_1, side_2, side_3, side_4
      integer :: rotate, intersects
      integer :: b3, b4, b5, b6, dn1, dn2, dn3, dn4
      real(rk) :: bx1, by1, bx2, by2, xx, xn, yx, yn
      real(rk), dimension(4) :: qx, qy
      real(rk) :: lxn, lxx, lyn, lyx, qxn, qxx, qyn, qyx

      pcolor(:) = -1
      qcolor(:) = -1

      lxn = MINVAL(bc(loops(bpl(l)))%x)
      lyn = MAXVAL(bc(loops(bpl(l)))%y)
      lxx = MINVAL(bc(loops(bpl(l)))%x)
      lyx = MAXVAL(bc(loops(bpl(l)))%y)
      do i = bpl(l), bpl(l + 1) - 1
        lxn = MIN(lxn, MINVAL(bc(loops(i))%x))
        lyn = MIN(lyn, MAXVAL(bc(loops(i))%y))
        lxx = MAX(lxx, MINVAL(bc(loops(i))%x))
        lyx = MAX(lyx, MAXVAL(bc(loops(i))%y))
      end do

      do q = 0, g%nelem - 1 ! loop over all elem
        qx = (/ x( g%elem(4 * q + 1) ), x( g%elem(4 * q + 2) ), &
              & x( g%elem(4 * q + 3) ), x( g%elem(4 * q + 4) ) /)
        qy = (/ y( g%elem(4 * q + 1) ), y( g%elem(4 * q + 2) ), &
              & y( g%elem(4 * q + 3) ), y( g%elem(4 * q + 4) ) /)

        side_1 = -1
        side_2 = -1
        side_3 = -1
        side_4 = -1
        n1 = g%elem(4 * q + 1)
        n2 = g%elem(4 * q + 2)
        n3 = g%elem(4 * q + 3)
        n4 = g%elem(4 * q + 4)

        do i = bpl(l), bpl(l + 1) - 1! loop through boundaries in loop l
          b = loops(i)
          xn = MINVAL(bc(b)%x); xx = MAXVAL(bc(b)%x)
          yn = MINVAL(bc(b)%y); yx = MAXVAL(bc(b)%y)

          if( MINVAL(qx) <= xx .or. MAXVAL(qx) >= xn .or. & 
            & MINVAL(qy) <= yx .or. MAXVAL(qy) >= yn )then
            ! if quad is in range of boundary
            do j = 1, size(bc(b)%x) - 1
              intersects = -1
              bx1 = bc(b)%x(j); bx2 = bc(b)%x(j + 1)
              by1 = bc(b)%y(j); by2 = bc(b)%y(j + 1)
              call check_intersection(qx, qy, bx1, by1, bx2, by2, intersects, rotate)
              if(intersects > 0)then

                if(intersects < 2)then    ! 4---3/
                 ! one corner deleted     ! |   /-->
                 !                        ! 1--/*

                  qcolor(q + 1) = 1

                  nn1 = 1 + MOD(1 + rotate, 4)
                  nn2 = 1 + MOD(3 + rotate, 4)
                  nn3 = 1 + MOD(2 + rotate, 4)
                  nn4 = 1 + rotate

                  nn1 = g%elem(4 * q + nn1)
                  nn2 = g%elem(4 * q + nn2)
                  nn3 = g%elem(4 * q + nn3)
                  nn4 = g%elem(4 * q + nn4)

                  call add_color(pcolor, nn1, nbpl, b)
                  call add_color(pcolor, nn2, nbpl, b)
                  call add_color(pcolor, nn4, nbpl, b)

                  pcolor(nbpl * (nn3 - 1) + 1) = 0

                else if(intersects < 3)then     ! 4-|-*
                  ! half of the square deleted  ! | | |
                  !                             ! 1-|-*

                  qcolor(q + 1) = 1

                  nn1 = 1 + rotate
                  nn2 = 1 + MOD(3 + rotate, 4)
                  nn1 = g%elem(4 * q + nn1)
                  nn2 = g%elem(4 * q + nn2)

                  call add_color(pcolor, nn1, nbpl, b)
                  call add_color(pcolor, nn2, nbpl, b)

                  nn1 = 1 + MOD(1 + rotate, 4)
                  nn2 = 1 + MOD(2 + rotate, 4)
                  nn1 = g%elem(4 * q + nn1)
                  nn2 = g%elem(4 * q + nn2)

                  pcolor(nbpl * (nn1 - 1) + 1) = 0
                  pcolor(nbpl * (nn2 - 1) + 1) = 0

                else if(intersects < 4)then       ! *---*/ 
                  ! all but one corner to delete  ! | <-/  
                  !                               ! *--/2  

                  qcolor(q + 1) = 1

                  nn1 = 1 + MOD(2 + rotate, 4)
                  nn1 = g%elem(4 * q + nn1)

                  call add_color(pcolor, nn1, nbpl, b)

                  nn1 = 1 + MOD(1 + rotate, 4)
                  nn2 = 1 + MOD(3 + rotate, 4)
                  nn3 = 1 + rotate
                  nn1 = g%elem(4 * q + nn1)
                  nn2 = g%elem(4 * q + nn2)
                  nn3 = g%elem(4 * q + nn3)

                  pcolor(nbpl * (nn1 - 1) + 1) = 0
                  pcolor(nbpl * (nn2 - 1) + 1) = 0

                else if(intersects < 5)then
                  ! only one edge intersected or else
                  ! only one corner intersected
                  if(side_1 < 0)then
                    side_1 = rotate
                    b3 = b
                  else if(side_2 < 0)then
                    side_2 = rotate
                    b4 = b
                  else if(side_3 < 0)then
                    side_3 = rotate
                    b5 = b
                  else if(side_4 < 0)then
                    side_4 = rotate
                    b6 = b
                  end if
                else
                  ! do nothing other intersection types do not give enough
                  ! information to delete nodes
                end if
              else
                ! do nothing. all elem default marked -1 for interior
              end if
            end do
          end if
        end do ! end loop over boundaries in loop l

        ! if only one side of the quad was intersected by each boundary segment
        if(side_1 > 0)then
          if(side_4 > 0)then
            if( (side_4 - side_3) < 1)then
              side_4 = -1
            end if
            if( (side_4 - side_2) < 1)then
              side_4 = -1
            end if
          end if
          if(side_3 > 0)then
            if( (side_3 - side_2) < 1)then
              side_3 = side_4
              b5 = b6
              side_4 = -1
            end if
            if( (side_3 - side_1) < 1)then
              side_3 = -1
            end if
          end if
          if( (side_1 - side_2) < 1)then
            if(side_3 > 0)then
              side_2 = side_3
              b4 = b5
              side_3 = side_4
              b5 = b6
              side_4 = -1
            end if
          end if

          if(side_2 > 0)then
            qcolor(q + 1) = 1
            dn1 = -1
            dn2 = -1
            nn1 = -1
            nn2 = -1
            if(side_1 < 5)then
              dn1 = side_1
              nn1 = 1 + MOD(dn1, 4)
            else
              dn1 = side_1 - 4
              nn1 = 1 + MOD(dn1 + 2, 4)
            end if
            if(side_2 < 5)then
              dn2 = side_2
              nn2 = 1 + MOD(dn2, 4)
            else
              dn2 = side_2 - 4
              nn2 = 1 + MOD(dn2 + 2, 4)
            end if

            if( abs(side_1 - side_2) < 1)then
              if(side_1 < 5)then
                dn1 = side_1
              else
                dn1 = side_1 - 4
              end if
              qcolor(q + 1) = 1
              nn1 = 1 + MOD(dn1, 4)
              nn2 = 1 + MOD(dn1 + 1, 4)
              nn3 = 1 + MOD(dn1 + 2, 4)
              nn1 = g%elem(4 * q + nn1)
              nn2 = g%elem(4 * q + nn2)
              nn3 = g%elem(4 * q + nn3)
              nn4 = g%elem(4 * q + dn1)
              call add_color(pcolor, nn1, nbpl, b3)
              call add_color(pcolor, nn2, nbpl, b3)
              call add_color(pcolor, nn2, nbpl, b4)
              call add_color(pcolor, nn3, nbpl, b4)
              pcolor(nbpl * (nn4 - 1) + 1) = 0
            else if(side_1 > 8 .and. side_2 > 8)then
              ! do nothing: another boundary segment 
              ! intersects two of the edges of this quad

            else if( abs(dn1 - dn2) < 1)then
              qcolor(q + 1) = 1
              nn3 = 1 + MOD(dn1 + 1, 4)
              nn1 = g%elem(4 * q + nn1)
              nn2 = g%elem(4 * q + nn2)
              nn3 = g%elem(4 * q + nn3)
              nn4 = g%elem(4 * q + dn1)
              call add_color(pcolor, nn1, nbpl, b3)
              call add_color(pcolor, nn2, nbpl, b3)
              call add_color(pcolor, nn2, nbpl, b4)
              call add_color(pcolor, nn3, nbpl, b4)
              pcolor(nbpl * (nn4 - 1) + 1) = 0
            else
              qcolor(q + 1) = 1
              nn1 = g%elem(4 * q + nn1)
              nn2 = g%elem(4 * q + nn2)
              nn3 = g%elem(4 * q + dn1)
              nn4 = g%elem(4 * q + dn2)
              call add_color(pcolor, nn1, nbpl, b3)
              call add_color(pcolor, nn2, nbpl, b4)
              pcolor(nbpl * (nn3 - 1) + 1) = 0
              pcolor(nbpl * (nn4 - 1) + 1) = 0
            end if
          else
          end if
        end if
      end do ! end loop over elem

      do q = 0, g%nelem - 1
        n1 = g%elem(4 * q + 1)
        n2 = g%elem(4 * q + 2)
        n3 = g%elem(4 * q + 3)
        n4 = g%elem(4 * q + 4)
        nn1 = pcolor(nbpl * (n1 - 1) + 1)
        nn2 = pcolor(nbpl * (n2 - 1) + 1)
        nn3 = pcolor(nbpl * (n3 - 1) + 1)
        nn4 = pcolor(nbpl * (n4 - 1) + 1)
        if( abs(nn1) < 1 )then
          if(nn2 > 0 .and. nn4 > 0)then
            do i = 1, nbpl
              if(pcolor(nbpl * (n2 - 1) + i) > 0)then
                call add_color(pcolor, n3, nbpl, pcolor(nbpl * (n2 - 1) + i))
              end if
              if(pcolor(nbpl * (n4 - 1) + i) > 0)then
                call add_color(pcolor, n3, nbpl, pcolor(nbpl * (n4 - 1) + i))
              end if
            end do
          end if
        end if
        if( abs(nn2) < 1 )then
          if(nn3 > 0 .and. nn1 > 0)then
            do i = 1, nbpl
              if(pcolor(nbpl * (n3 - 1) + i) > 0)then
                call add_color(pcolor, n4, nbpl, pcolor(nbpl * (n3 - 1) + i))
              end if
              if(pcolor(nbpl * (n1 - 1) + i) > 0)then
                call add_color(pcolor, n4, nbpl, pcolor(nbpl * (n1 - 1) + i))
              end if
            end do
          end if
        end if
        if( abs(nn3) < 1 )then
          if(nn4 > 0 .and. nn2 > 0)then
            do i = 1, nbpl
              if(pcolor(nbpl * (n4 - 1) + i) > 0)then
                call add_color(pcolor, n1, nbpl, pcolor(nbpl * (n4 - 1) + i))
              end if
              if(pcolor(nbpl * (n2 - 1) + i) > 0)then
                call add_color(pcolor, n1, nbpl, pcolor(nbpl * (n2 - 1) + i))
              end if
            end do
          end if
        end if
        if( abs(nn4) < 1 )then
          if(nn1 > 0 .and. nn3 > 0)then
            do i = 1, nbpl
              if(pcolor(nbpl * (n1 - 1) + i) > 0)then
                call add_color(pcolor, n2, nbpl, pcolor(nbpl * (n1 - 1) + i))
              end if
              if(pcolor(nbpl * (n3 - 1) + i) > 0)then
                call add_color(pcolor, n2, nbpl, pcolor(nbpl * (n3 - 1) + i))
              end if
            end do
          end if
        end if
      end do ! end loop over elem
    end subroutine color_elem

    ! delete all quads that are completely outside of domain
    subroutine delete_external_elem(g, nbpl, pcolor, qcolor)
      type(mygrid), intent(in out) :: g
      integer,   intent(in) :: nbpl
      integer, dimension(:), intent(in out) :: pcolor, qcolor
      integer :: q, n1, n2, n3, n4, check
      integer :: c1, c2, c3, c4
      integer, dimension(:), allocatable :: tag

      check = -1
      do while(check < 0)

        check = 1

        do q = 1, g%nelem                             ! loop through elem
          if(qcolor(q) < 1)then
            n1 = g%elem(4 * (q - 1) + 1) - 1
            n2 = g%elem(4 * (q - 1) + 2) - 1
            n3 = g%elem(4 * (q - 1) + 3) - 1
            n4 = g%elem(4 * (q - 1) + 4) - 1

            c1 = pcolor(nbpl * n1 + 1)
            c2 = pcolor(nbpl * n2 + 1)
            c3 = pcolor(nbpl * n3 + 1)
            c4 = pcolor(nbpl * n4 + 1)

            if( & !if all of the nodes in quad q are either deleted or not marked
              & ( (c1 < 1) .and. (c2 < 1) .and. (c3 < 1) .and. (c4 < 1) &
              & ) .and. ( & ! and at least one is marked for deleting
              &   (abs(c1) < 1) .or. (abs(c2) < 1) .or. & 
              &   (abs(c3) < 1) .or. (abs(c4) < 1) ) )then

              qcolor(q) = 1                         ! color quad for delete

              check = -1
              if(pcolor(nbpl * n1 + 1) < 1)then     ! delete all other nodes
                pcolor(nbpl * n1 + 1) = 0           ! that are not colored
              end if                                  ! with a boundary color
              if(pcolor(nbpl * n2 + 1) < 1)then
                pcolor(nbpl * n2 + 1) = 0
              end if
              if(pcolor(nbpl * n3 + 1) < 1)then
                pcolor(nbpl * n3 + 1) = 0
              end if
              if(pcolor(nbpl * n4 + 1) < 1)then
                pcolor(nbpl * n4 + 1) = 0
              end if

            end if
          end if
        end do
      end do

      do q = 1, g%nelem
        n1 = g%elem(4 * (q - 1) + 1) - 1
        n2 = g%elem(4 * (q - 1) + 2) - 1
        n3 = g%elem(4 * (q - 1) + 3) - 1
        n4 = g%elem(4 * (q - 1) + 4) - 1

        c1 = abs(pcolor(nbpl * n1 + 1))
        c2 = abs(pcolor(nbpl * n2 + 1))
        c3 = abs(pcolor(nbpl * n3 + 1))
        c4 = abs(pcolor(nbpl * n4 + 1))

        if((c1 < 1) .or. (c2 < 1) .or. (c3 < 1) .or.(c4 < 1) )then
          qcolor(q) = 1
        end if
      end do

      ! delete any node that is only connected to deleted elem
      allocate(tag(g%nn))
      tag(:) = -1
      do q = 0, g%nelem - 1 
       n1 = g%elem(4 * q + 1)
       n2 = g%elem(4 * q + 2)
       n3 = g%elem(4 * q + 3)
       n4 = g%elem(4 * q + 4)
       if(abs(pcolor(nbpl * (n1 - 1) + 1)) > 0)then
         tag(n2) = 1
         tag(n4) = 1
       end if
       if(abs(pcolor(nbpl * (n2 - 1) + 1)) > 0)then
         tag(n1) = 1
         tag(n3) = 1
       end if
       if(abs(pcolor(nbpl * (n3 - 1) + 1)) > 0)then
         tag(n2) = 1
         tag(n4) = 1
       end if
       if(abs(pcolor(nbpl * (n4 - 1) + 1)) > 0)then
         tag(n1) = 1
         tag(n3) = 1
       end if

      end do
      do q = 1, g%nn
        if(tag(q) < 0)then
          pcolor(nbpl * (q - 1) + 1) = 0
        end if
      end do
      deallocate(tag)
    end subroutine delete_external_elem

    ! determine number of nodes to distribute about boundary curves
    subroutine get_npb(nnodes, pcolor, nbpl, b_nodes)
      integer,               intent(in)     :: nnodes, nbpl
      integer, dimension(:), intent(in)     :: pcolor
      integer, dimension(:), intent(in out) :: b_nodes
      integer :: i, j, k, m

      b_nodes(:) = 0
      do i = 0, nnodes - 1
        m = pcolor(nbpl * i + 1)
        do j = 1, nbpl
          k = pcolor(nbpl * i + j)
          if( (k > 0) .and. (m > 0) )then
            b_nodes(k) = b_nodes(k) + 1
          end if
        end do
      end do
    end subroutine get_npb

    ! order consecutive edges that form connected loop around boundary 
    subroutine order_b_conn(bconn, npb, g, l, bpl, loops, nbpl, pcolor, qcolor)
      integer,   dimension(:), intent(in out) :: bconn, npb
      type(mygrid),            intent(in out) :: g
      integer,                 intent(in)     :: l, nbpl
      integer,   dimension(:), intent(in)     :: bpl, loops, pcolor, qcolor
      integer :: i, j, k, m, n, n1, n2, n3, n4, c1, c2, check, tct, bct
      integer, dimension(:), allocatable :: tmp, mark

      allocate(tmp(4 * SUM(npb)))
      allocate(mark(g%nn))
      mark(:) = -1
      tmp(:) = -1

      tct = 0
      do i = 0, g%nelem - 1
        if(qcolor(i + 1) > 0)then
          do j = 1, 4
            n1 = g%elem(4 * i + j)
            n2 = g%elem(4 * i + 1 + MOD(j, 4))

            check = -1
            if(n1 < 0 .or. n2 < 0)then
            else if(abs(pcolor(nbpl * (n1 - 1) + 1)) < 1 .or. &
                  & abs(pcolor(nbpl * (n2 - 1) + 1)) < 1)then
            else
              do m = 1, nbpl
                do n = 1, nbpl
                  c1 = pcolor( nbpl * (n1 - 1) + m )
                  c2 = pcolor( nbpl * (n2 - 1) + n )

                  if(c1 > 0 .and. c2 > 0)then
                    check = 1
                  end if

                  if(check > 0)then
                    do k = 1, SUM(npb)
                      n3 = tmp(2 * k - 1)
                      n4 = tmp(2 * k)
                      if(abs(n1 - n3) < 1 .and. abs(n2 - n4) < 1)then
                        check = -1
                      end if
                    end do
                    if(check > 0)then
                      tmp(2 * tct + 1) = n1
                      tmp(2 * tct + 2) = n2
                      tct = tct + 1
                      check = -1
                    end if
                  end if
                end do ! m
              end do ! n
            end if  ! if n1 and n2 > 0 and neither colored for deleting
          end do  ! j
        end if ! if quad not deleted
      end do ! i (quad loop)
      tct = tct - 1

      bconn(:) = -1
      bconn(1) = tmp(2)
      bconn(2) = tmp(1)
      bct = 3
      tmp(1:2) = -1

      n1 = bconn(1)
      n2 = bconn(2)
      check = -1
      i = 0
      do while(i < 2 * tct .and. check < 0)
        i = i + 1
        do j = 1, SUM(npb)
          if(check < 0)then
            n4 = tmp(2 * j - 1)
            n3 = tmp(2 * j)

            if(n3 > 0 .and. n4 > 0)then
              if(abs(n2 - n3) < 1)then
                if(abs(n4 - n1) > 0)then
                  bconn(bct) = n4
                  bct = bct + 1
                  n2 = n4

                  tmp(2 * j - 1) = -1
                  tmp(2 * j) = -1
                else
                  check = 1
                end if
              end if ! n2 == n3
            end if ! n3, n4 > 0
          end if ! check
        end do ! j
      end do !do while
      deallocate(tmp)

      bct = bct - 1
      deallocate(mark)

      npb(:) = 0
      do i = 1, bct
        c1 = pcolor(nbpl * (bconn(i) - 1) + 1)
        npb(c1) = npb(c1) + 1
      end do

    end subroutine order_b_conn

    ! update quad and node numbering after deletion
    subroutine renumber_arrays(g, bconn, x, y, nbpl, b_edges, be_ct, pcolor, qcolor)
      type(mygrid),              intent(in out) :: g
      integer,   dimension(:),   intent(in out) :: bconn
      real(rk),  dimension(:),   intent(in out) :: x, y
      integer,                   intent(in)     :: nbpl
      integer,   dimension(:,:), intent(in out) :: b_edges
      integer,                   intent(in)     :: be_ct
      integer,   dimension(:),   intent(in out) :: pcolor, qcolor
      integer, dimension(:), allocatable :: qrenumber
      integer, dimension(:), allocatable :: nrenumber
      integer, dimension(:), allocatable :: inverse
      integer :: i, j, n, q, nnq, nnp, n1, n2, n3, n4

      allocate(qrenumber(g%nelem))
      allocate(nrenumber(g%nn))
      allocate(inverse(g%nn))

      !==========================================!
      !    renumber elem, deleting any with
      !             flagged nodes
      nnq = 1
      do q = 1, g%nelem        ! loop through elem
        if(qcolor(q) < 1)then  ! if all nodes are interior, keep
          qrenumber(nnq) = q
          nnq = nnq + 1
        else                   ! otherwise skip
          ! do nothing
        end if
      end do
      qrenumber(nnq:) = -1
      nnq = nnq - 1
      !==========================================!


      !==========================================!
      !      now update the quad numbering
      do i = 1, nnq
        n = qrenumber(i)
        qcolor(i) = qcolor(n)
        g%elem(4 * (i - 1) + 1 : 4 * (i - 1) + 4) = &
        & g%elem(4 * (n - 1) + 1 : 4 * (n - 1) + 4)
      end do
      g%nelem = nnq
      !==========================================!


      !==========================================!
      ! renumber nodes, deleting any marked '0'
      inverse(:) = -1
      nrenumber(:) = -1
      nnp = 1
      do n = 1, g%nn
        if(abs(pcolor(nbpl * (n - 1) + 1)) > 0)then
          nrenumber(nnp) = n
          inverse(n) = nnp
          nnp = nnp + 1
        else
          ! do nothing
        end if
      end do
      nnp = nnp - 1
      !==========================================!


      !==========================================!
      ! now update the quad-node connectivity
      do i = 0, nnq - 1
        do j = 1, 4
          n = g%elem(4 * i + j)
          g%elem(4 * i + j) = inverse(n)
        end do
      end do
      g%elem(4 * nnq + 1 : ) = -1
      qcolor(nnq + 1 : ) = -1
      j = nnq
      do i = nnq, 1, -1
        do n = 1, 4
          if(g%elem(4 * (i - 1) + n) < 1)then
            j = i - 1
          end if
        end do
      end do
      g%nelem = j
      !==========================================!


      !==========================================!
      !update x, y coordinates with new numbering
      do i = 1, nnp
        n = nrenumber(i)
        x(i) = x(n)
        y(i) = y(n)
        pcolor(nbpl * (i - 1) + 1 : nbpl * (i - 1) + nbpl) = &
        pcolor(nbpl * (n - 1) + 1 : nbpl * (n - 1) + nbpl)
      end do
      g%nn = nnp
      !==========================================!

      do i = 1, be_ct
        b_edges(i,1) = inverse(b_edges(i,1))
        b_edges(i,2) = inverse(b_edges(i,2))
      end do

      !==========================================!
      !  renumber bconn with new numbering
      do i = 1, size(bconn)
        if(bconn(i) > 0)then
          bconn(i) = inverse(bconn(i))
        end if
      end do
      n = 0
      do i = size(bconn), 1,-1
        if(n > 0)then
          if(bconn(i) > 0)then
          else
            n = 1
            j = i
          end if
        end if
      end do
      if(n > 0)then 
        bconn(j : ) = -1
      end if
      !==========================================!

      deallocate(qrenumber)
      deallocate(nrenumber)
      deallocate(inverse)
    end subroutine renumber_arrays

    subroutine connect_interior_boundary(g, bc, x, y, bx, by, npb, bconn, & 
                      & new_elem, l, bpl, loops, nbpl, b_edges, be_ct, be_vals, pcolor)
      type(mygrid),                              intent(in out) :: g
      type(curve),  dimension(:),                intent(in out)     :: bc
      real(rk),     dimension(:), allocatable,   intent(in out) :: x, y
      real(rk),     dimension(:),                intent(in out) :: bx, by
      integer,      dimension(:),                intent(in out) :: npb
      integer,      dimension(:),                intent(in out) :: bconn
      integer,      dimension(:),                intent(in out) :: new_elem
      integer,                                   intent(in)     :: l, nbpl
      integer,      dimension(:),                intent(in)     :: bpl, loops
      integer,      dimension(:,:),              intent(in out) :: b_edges
      integer,                                   intent(in)     :: be_ct
      integer,      dimension(:),   allocatable, intent(in out) :: be_vals
      integer,      dimension(:),                intent(in)     :: pcolor
      integer,      dimension(:), allocatable :: bidx
      integer,      dimension(:), allocatable :: bnum
      integer,      dimension(:), allocatable :: tbv
      integer :: nq, nnc, nnp, nnq
      integer :: i, j, n1, n2, n3, n4, q, k, bct, b

      nnc = SUM(npb)
      allocate(tbv(be_ct))
      tbv(:be_ct) = be_vals(:be_ct)
      deallocate(be_vals)
      allocate(be_vals(be_ct + nnc))
      be_vals(:be_ct) = tbv(:be_ct)
      deallocate(tbv)

      npb(:) = 0
      do i = 1, nnc
        n1 = bconn(i)
        n2 = pcolor(nbpl * (n1 - 1) + 1)
        npb(n2) = npb(n2) + 1
      end do

      nq = SUM(npb)
      allocate( bidx(SUM(npb)) )
      allocate( bnum(SUM(npb)) )
      bct = 1
      do i = bpl(l), bpl(l + 1) - 1
        do j = 1, npb(loops(i))
          bnum(bct) = g%nn + bct
          bct = bct + 1
        end do
      end do

      call distribute_points(bc, npb, bx, by, x, y, bconn, bidx, bpl, loops, l)
      call get_best_connections(bc, bpl, loops, l, nbpl, g%nn, x, y, bx, by, bconn, npb, bnum)
      call distribute_points(bc, npb, bx, by, x, y, bconn, bidx, bpl, loops, l)

      new_elem(:) = -1
      do q = 0, nnc - 1
        n1 = bnum(q + 1)
        n2 = bnum(1 + MOD(q + 1, SUM(npb)))
        n3 = bconn(1 + MOD(q + 1, SUM(npb)))
        n4 = bconn(q + 1)
        if(n1 > 0 .and. n2 > 0 .and. n3 > 0 .and. n4 > 0)then
          new_elem(4 * q + 1) = n1
          new_elem(4 * q + 2) = n2
          new_elem(4 * q + 3) = n3
          new_elem(4 * q + 4) = n4

          b_edges(be_ct + 1 + q, 1) = n1
          b_edges(be_ct + 1 + q, 2) = n2

          be_vals(be_ct + 1 + q) = bidx(n1 - g%nn)
          be_vals(be_ct + 1 + q) = bidx(n1 - g%nn)
          be_vals(be_ct + 1 + q) = bidx(n2 - g%nn)
        end if
      end do

      ! change x, y coordinates of bc to be the recently added values
      bct = 1
      do i = bpl(l), bpl(l + 1) - 1
        b = loops(i)
        j = npb(b)
        deallocate(bc(b)%x, bc(b)%y)
        allocate(bc(b)%x(j), bc(b)%y(j))
        bc(b)%x(:j) = bx(bct : bct + j - 1)
        bc(b)%y(:j) = by(bct : bct + j - 1)
        bct = bct + j
      end do

      deallocate( bidx )
      deallocate( bnum )

      nnp = g%nn + SUM(npb)

      bx(g%nn + 1 : nnp) = bx(1:SUM(npb))
      by(g%nn + 1 : nnp) = by(1:SUM(npb))
      bx(1:g%nn) = x(1:g%nn)
      by(1:g%nn) = y(1:g%nn)

      deallocate(x)
      deallocate(y)
      allocate(x(nnp))
      allocate(y(nnp))
      x(:nnp) = bx(:nnp)
      y(:nnp) = by(:nnp)
      g%nn = nnp

      nnq = SUM(npb)

      new_elem(4 * g%nelem + 1 : 4 * g%nelem + 4 * nnq) =  new_elem(1 : 4 * nnq)
      new_elem(1 : 4 * g%nelem) = g%elem(1 : 4 * g%nelem)

      call get_nq(new_elem, nnq)
      g%nelem = nnq

      deallocate(g%elem)
      allocate(g%elem(4 * g%nelem))
      g%elem(:4 * g%nelem) = new_elem(:4*g%nelem)
    end subroutine connect_interior_boundary

    subroutine get_best_connections(bc, bpl, loops, l, nbpl, nn, x, y, bx, by, bconn, npb, bnum)
      type(curve), dimension(:), intent(in)     :: bc
      integer,       dimension(:), intent(in)     :: bpl, loops
      integer,                     intent(in)     :: nbpl, l, nn
      real(rk),      dimension(:), intent(in)     :: x, y, bx, by
      integer,       dimension(:), intent(in out) :: bconn, npb, bnum

      ! first find closest point between boundary loop l and bx, by(bnum)
      call find_closest_point(bconn, bnum, x, y, bx, by, nn)

      ! then try walking over boundaries
      if(nbpl > 1)then
        call get_optimal_npb(bc, bconn, bnum, x, y, npb, bpl, loops, l, nn)
      end if

    end subroutine get_best_connections

    subroutine find_closest_point(edge_num, boundary_num, x, y, bx, by, nn)
      integer,  dimension(:), intent(in out) :: edge_num, boundary_num
      real(rk), dimension(:), intent(in)     :: x, y, bx, by
      integer,                intent(in)     :: nn
      integer, dimension(:), allocatable :: tmp      
      integer :: i, j, idx1, idx2, n1, n2, new_size
      real(rk) :: d, nd

      if(size(bx) > 0 .and. size(x) > 0)then

        idx1 = 1
        idx2 = 1
        d = (x(idx1) - bx(idx1))**2 + (y(idx1) - by(idx1))**2

        do i = 1, size(boundary_num)
          n1 = boundary_num(i) - nn
          if(n1 > 0)then
            do j = 1, size(edge_num)
              n2 = edge_num(j)
              if(n2 > 0)then
                nd = (x(n2) - bx(n1) )**2 + (y(n2) - by(n1))**2
                if(nd < d)then
                  idx1 = i
                  idx2 = j
                  d = nd
                end if
              end if
            end do
          end if
        end do

        new_size = size(boundary_num)
        do i = size(boundary_num), 2, -1
          if(boundary_num(i) < 1)then
            new_size = i - 1
          end if
        end do
        allocate(tmp(new_size))

        tmp(:new_size) = boundary_num(:new_size)
        do i = 0, new_size - 1
          boundary_num(i + 1) = tmp(1 + MOD(idx1 + i - 1, new_size))
        end do

        tmp(:new_size) = edge_num(:new_size)

        do i = 0, new_size - 1
          edge_num(i + 1) = tmp(1 + MOD(idx2 + i - 1 , new_size))
        end do

        deallocate(tmp)
      end if
    end subroutine find_closest_point

    subroutine get_optimal_npb(bc, bconn, bnum, x, y, npb, bpl, loops, l, ngn)
      type(curve), dimension(:), intent(in)     :: bc
      integer,       dimension(:), intent(in out) :: bconn, npb
      integer,       dimension(:), intent(in out) :: bnum
      integer,       dimension(:), intent(in)     :: bpl, loops
      real(rk),      dimension(:), intent(in)     :: x, y
      integer,                     intent(in)     :: l, ngn
      integer,  dimension(:), allocatable :: nnpb, tmp, tbc
      real(rk), dimension(:), allocatable :: bx, by
      integer :: nnp, i, j, k, b1, b2, n1, n2, ninv, npos, ninv1, npos1
      integer :: i1, i2, c1, c2, n_1, n_2
      nnp = SUM(npb)
      allocate(bx(nnp), by(nnp))
      allocate(nnpb(size(npb)))
      allocate(tmp(nnp), tbc(nnp))

      call create_tmp_points(loops, bpl, l, bc, bx, by, npb)
      call get_n_inverted(ninv, npos, bconn, bnum, x, y, bx, by, ngn)

      if(ninv < 1)then
      else
        do i = bpl(l), bpl(l + 1) - 1
          b1 = loops(i)
          if(i < (bpl(l + 1) - 1))then
            b2 = loops(i + 1)
          else
            b2 = loops(bpl(l))
          end if

          n_1 = npb(b1)
          n_2 = npb(b2)
          c1 = n_1 / 2
          tbc(:nnp) = bconn(:nnp)

          do j = 1, c1
            n1 = n_1 - j
            n2 = n_2 + j
            nnpb(:) = npb(:)
            nnpb(b1) = n1
            nnpb(b2) = n2

            if(i > (bpl(l + 1) - 2))then
              do k = 1, nnp
                tmp(k) = tbc(1 + MOD(k + nnp - (j + 1), nnp))
              end do
            else
              tmp(:nnp) = bconn(:nnp)
            end if
            call create_tmp_points(loops, bpl, l, bc, bx, by, nnpb)
            call get_n_inverted(ninv1, npos1, tmp, bnum, x, y, bx, by, ngn)
            if( (npos1 .ge. npos) .and. (ninv1 .le. ninv) )then
              npb(:) = nnpb(:)
              ninv = ninv1
              npos = npos1
              bconn(:nnp) = tmp(:nnp)
            else
            end if
          end do

          tbc(:nnp) = bconn(:nnp)
          n_1 = npb(b1)
          n_2 = npb(b2)
          c2 = n_2 / 2

          do j = 1, c2
            n1 = n_1 + j
            n2 = n_2 - j
            nnpb(:) = npb(:)
            nnpb(b1) = n1
            nnpb(b2) = n2
            if(i > (bpl(l + 1) - 2))then
              do k = 1, nnp
                tmp(k) = tbc(1 + MOD(k + j - 1, nnp))
              end do
            else
              tmp(:nnp) = bconn(:nnp)
            end if
            call create_tmp_points(loops, bpl, l, bc, bx, by, nnpb)
            call get_n_inverted(ninv1, npos1, tmp, bnum, x, y, bx, by, ngn)
            if( (npos1 .ge. npos) .and. (ninv1 .le. ninv) )then
              npb(:) = nnpb(:)
              ninv = ninv1
              npos = npos1
              bconn(:nnp) = tmp(:nnp)
            else
            end if
          end do
        end do
      end if

      deallocate(bx, by)
      deallocate(nnpb, tmp)
    end subroutine get_optimal_npb

    subroutine get_n_inverted(ninv, npos, bconn, bnum, x, y, bx, by, ngn)
      integer,                intent(in out) :: ninv, npos
      integer,  dimension(:), intent(in)     :: bconn, bnum
      real(rk), dimension(:), intent(in)     :: x, y, bx, by
      integer,                intent(in)     :: ngn
      integer :: i, n1, n2, n3, n4, check
      real(rk) :: x1, y1, x2, y2, x3, y3, x4, y4

      ninv = 0
      npos = 0

      do i = 1, size(bnum)
        n1 = bconn(1 + MOD(i, size(bnum)))
        n2 = bconn(i)
        n3 = bnum(i) - ngn
        n4 = bnum(1 + MOD(i, size(bnum))) - ngn

        x1 = x(n1);  x2 = x(n2)
        y1 = y(n1);  y2 = y(n2)
        x3 = bx(n3); x4 = bx(n4)
        y3 = by(n3); y4 = by(n4)

        check = 0
        if( ( (x2 - x1) * (y4 - y1) - (y2 - y1) * (x4 - x1)) < tol)then
          ninv = ninv + 1
          check = 1
        end if
        if( ( (x3 - x2) * (y1 - y2) - (y3 - y2) * (x1 - x2)) < tol)then
          ninv = ninv + 1
          check = 1
        end if
        if( ( (x4 - x3) * (y2 - y3) - (y4 - y3) * (x2 - x3)) < tol)then
          ninv = ninv + 1
          check = 1
        end if
        if( ( (x1 - x4) * (y3 - y4) - (y1 - y4) * (x3 - x4)) < tol)then
          ninv = ninv + 1
          check = 1
        end if
        if(check < 1)then
          npos = npos + 1
        end if
      end do
    end subroutine get_n_inverted

    subroutine get_nq(nq, nnq)
      integer, dimension(:), intent(in)     :: nq
      integer,               intent(in out) :: nnq
      integer :: i, j, tnq, n1, n2, n3, n4

      j = size(nq) / 4
      tnq = j
      do i = 1, j
        n1 = nq(4 * (i - 1) + 1)
        n2 = nq(4 * (i - 1) + 2)
        n3 = nq(4 * (i - 1) + 3)
        n4 = nq(4 * (i - 1) + 4)

        if(n1 < 0 .or. n2 < 0 .or. n3 < 0 .or. n4 < 0)then
          tnq = tnq - 1
        end if
      end do
      nnq = tnq
    end subroutine get_nq

    ! check if any quads have bad condition number. If so, then 
    ! subdivide into triangles
    subroutine subdivide(x, y, elem, etype, nq, nt)
      real(rk), dimension(:),              intent(in)     :: x, y
      integer,  dimension(:), allocatable, intent(in out) :: elem, etype
      integer, intent(in out) :: nq, nt
      integer :: n1, n2, n3, n4
      integer, dimension(:), allocatable :: mark
      integer, dimension(:), allocatable :: tri
      integer :: q, t, choice1
      real(rk) :: x1, x2, x3, x4, y1, y2, y3, y4
      real(rk) :: a1, a2, a3, a4

      allocate(mark(nq))
      nt = 0
      mark(:) = 0
      do q = 0, nq - 1
        n1 = elem(4 * q + 1)
        n2 = elem(4 * q + 2)
        n3 = elem(4 * q + 3)
        n4 = elem(4 * q + 4)

        x1 = x(n1)
        x2 = x(n2)
        x3 = x(n3)
        x4 = x(n4)

        y1 = y(n1)
        y2 = y(n2)
        y3 = y(n3)
        y4 = y(n4)

        call check_condition_number(1, x1, y1, x2, y2, x3, y3, x4, y4, a1)
        if(abs(a1 - 1) > 1.0)then
          mark(q + 1) = 2
        end if
      end do
      nt = SUM(mark)


      allocate(tri(3 * SUM(mark) + 4 * (nq - (nt / 2))))
      nt = 0
      do q = 1, nq
        if(mark(q) > 0)then
          n1 = elem(4 * (q - 1) + 1)
          n2 = elem(4 * (q - 1) + 2)
          n3 = elem(4 * (q - 1) + 3)
          n4 = elem(4 * (q - 1) + 4)

          x1 = x(n1)
          x2 = x(n2)
          x3 = x(n3)
          x4 = x(n4)

          y1 = y(n1)
          y2 = y(n2)
          y3 = y(n3)
          y4 = y(n4)

          ! check both triangle options: 
          !
          ! choice1:  choice2:
          !  4____3    4____3  
          !  |   /|    |\   |  
          !  |  / |    | \  |  
          !  | /  |    |  \ |  
          !  |/   |    |   \|  
          !  1----2    1----2  
          ! and choose the best option
          call get_triangle_choice(choice1, x1, x2, x3, x4, y1, y2, y3, y4)

          if(choice1 > 0)then
            tri(3 * nt + 1) = n1
            tri(3 * nt + 2) = n2
            tri(3 * nt + 3) = n3
            nt = nt + 1
            tri(3 * nt + 1) = n1
            tri(3 * nt + 2) = n3
            tri(3 * nt + 3) = n4
            nt = nt + 1
          else
            tri(3 * nt + 1) = n1
            tri(3 * nt + 2) = n2
            tri(3 * nt + 3) = n4
            nt = nt + 1
            tri(3 * nt + 1) = n2
            tri(3 * nt + 2) = n3
            tri(3 * nt + 3) = n4
            nt = nt + 1
          end if
        end if
      end do
      
      n2 = 3 * nt 
      do q = 0, nq - 1
        if(mark(q + 1) < 1)then
          tri(n2 + 1 : n2 + 4) = elem(4 * q + 1 : 4 * q + 4)
          n2 = n2 + 4
        end if
      end do
      nq = nq - ((nt + 1) / 2)

      deallocate(elem)
      allocate(elem(size(tri)))
      elem(:) = tri(:)
      deallocate(etype)
      allocate(etype(nq + nt))
      etype(1 : nt) = triangle
      etype(nt + 1 : nt + nq) = quadrilateral
      deallocate(tri)
      deallocate(mark)
    end subroutine subdivide

    ! n.b. this actually computes aspect ratio. works better than condition number
    subroutine check_condition_number(etype, x1, y1, x2, y2, x3, y3, x4, y4, cn)
      integer, intent(in) :: etype
      real(rk), intent(in) :: x1, y1, x2, y2, x3, y3, x4, y4
      real(rk), intent(in out) :: cn
      real(rk) :: ux, uy, vx, vy, wx, wy, zx, zy, tx, ty
      real(rk) :: a, b, c, d, e

      if(etype < 1)then
        ux = x2 - x1
        uy = y2 - y1
        vx = x3 - x2
        vy = y3 - y2
        wx = x3 - x1
        wy = y3 - y1
        a = sqrt(ux**2 + uy**2)
        b = sqrt(vx**2 + vy**2)
        c = sqrt(wx**2 + wy**2)
        d = 0.5_rk * (a + b + c)
        e = 8.0_rk * (d - a) * (d - b) * (d - c)
        if( abs(e) > 0.0_rk )then
          cn = a * b * c / e
        else
          cn = 1.0_rk
        end if
      else if(etype < 2)then
        a = 0.25_rk * (x3 + x4 - (x1 + x2))**2
        b = 0.25_rk * (y3 + y4 - (y1 + y2))**2
        c = 0.25_rk * (x2 + x3 - (x4 + x1))**2
        d = 0.25_rk * (y2 + y3 - (y4 + y1))**2
        
        a = sqrt(a + b)
        b = sqrt(c + d)
        if(a > 0.0_rk .and. b > 0.0_rk)then
          cn = MAX(a / b, b / a)
        else
          cn = MAX(a, b)
        end if

        ux = x2 - x1
        uy = y2 - y1
        vx = x3 - x1
        vy = y3 - y1
        wx = x4 - x1
        wy = y4 - y1
        zx = x4 - x2
        zy = y4 - y2
        tx = x3 - x2
        ty = y3 - y2
        ! now check that the element has no inverted corners.
        ! if it does, then set condition number high above tolerance
        if( (ux * vy - uy * vx) < 0 .or. &
          & (vx * wy - vy * wx) < 0 .or. &
          & (tx * zy - ty * zx) < 0.0_rk .or. &
          & (ux * zy - zx * uy) < 0.0_rk)then
          cn = 50.0_rk
        end if
      else
        write(*,*)'error: in routine check_condition_number: element type not recognized'
      end if
    end subroutine check_condition_number

    subroutine get_triangle_choice(choice, x1, x2, x3, x4, y1, y2, y3, y4)
      integer, intent(in out) :: choice
      real(rk), intent(in) :: x1, x2, x3, x4, y1, y2, y3, y4
      real(rk) :: cn1, cn2, cn3, cn4

      call check_condition_number(0, x1, y1, x2, y2, x3, y3, 0.0_rk, 0.0_rk, cn1)
      call check_condition_number(0, x3, y3, x4, y4, x1, y1, 0.0_rk, 0.0_rk, cn2)
      call check_condition_number(0, x1, y1, x2, y2, x4, y4, 0.0_rk, 0.0_rk, cn3)
      call check_condition_number(0, x3, y3, x4, y4, x2, y2, 0.0_rk, 0.0_rk, cn4)
      if( (cn1 + cn2) < (cn3 + cn4) )then
        choice = 1
      else
        choice = -1
      end if
      if( ((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)) < 0 .or. &
        & ((x3 - x1) * (y4 - y1) - (x4 - x1) * (y3 - y1)) < 0)then
        choice = -1
      else if( ((x2 - x1) * (y4 - y1) - (x4 - x1) * (y2 - y1)) < 0 .or. &
        & ((x3 - x2) * (y4 - y2) - (x4 - x2) * (y3 - y2)) < 0)then
        choice = 1
      end if
    end subroutine get_triangle_choice

end module quad_gen
