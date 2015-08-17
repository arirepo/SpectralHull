module polygauss
  implicit none

  private

  real*8, parameter :: PI = 4.D0*DATAN(1.D0)

  ! ! the following two lines are public only when the tester program 
  ! ! in this file is active and uncommented otherwise comment them
  ! public :: repmat, points2distances, cubature_rules_1D, auto_rotation
  ! public :: eval_bench_integ, write_quad_to_file

  ! the following is always public
  public :: polygon_gauss_leg

contains

  ! This is a Fortran implementation of the Gauss-Legendre-Like quadrature
  ! algorithm described in the following paper:
  ! 
  ! A. SOMMARIVA and M. VIANELLO
  ! "Gauss-like and triangulation-free cubature over polygons".
  ! 
  ! input:
  ! 
  ! N     : degree of the one-dimensional gauss-legendre rule.
  ! 
  ! polygon_sides: if the polygon has "l" sides, is a
  !          variable containing its vertices, ordered counterclockwise.
  !          as last row must have the components of the first vertex.
  !          in other words, the first row and last row are equal.
  !          Therefore "polygon_sides" is a "l+1 x 2" matrix.
  !
  !  P, Q: direction that fixes the rotation. only for case rotation = 2.
  !
  !
  !  output:
  ! 
  !  xyw     : the gauss like formula produces the nodes (xyw(:,1),xyw(:,2))
  !            which are the xy coordinates of the Gauss-Legendre quadrature points
  !            and the weights xyw(:,3) of a cubature rule on the polygon.
  !
  !             --------- optional arguments ---------
  ! 
  !  rotation: not present : no rotation.
  !            1: automatic.
  !            2: preferred direction rotation by P, Q.
  ! 
  
  !
  !
  ! Fortran2003 version with a few modifications by:
  !
  ! Arash Ghasemi
  ! PhD candidate, University of Tennessee at Chattanooga
  !
  ! ghasemi.arash@gmail.com
  ! 
  ! %--------------------------------------------------------------------------
  
  subroutine polygon_gauss_leg(N, polygon_sides, P, Q, xyw, rotation)
    implicit none
    integer, intent(in) :: N
    real*8, dimension(:, :), intent(inout) :: polygon_sides
    real*8, dimension(:), allocatable :: P, Q
    real*8, dimension(:, :), allocatable :: xyw
    integer, intent(in), optional :: rotation

    ! local vars
    real*8, dimension(size(polygon_sides, 1)) :: x_bd, y_bd
    real*8 :: x_min, x_max, y_min, y_max
    real*8, dimension(2, 2) :: rot_matrix 
    real*8, dimension(2) :: axis_abscissa
    real*8, dimension(size(polygon_sides, 1), size(polygon_sides, 2)) :: polygon_bd_rot
    real*8 :: rot_angle, nrm_vect, a, x1, x2, y1, y2, half_pt_x, half_pt_y
    real*8, dimension(:), allocatable, target :: s_N, w_N, s_M, w_M
    integer :: N_length, M, L, index_side, M_length
    real*8 :: half_length_x, half_length_y
    real*8, dimension(:), pointer :: s_M_loc => null(), w_M_loc => null()
    real*8, dimension(:), allocatable :: x_gauss_side, local_weights, xx, yy, x_rot0, y_rot0
    real*8, dimension(:, :), allocatable :: scaling_fact_plus, scaling_fact_minus, y_gauss_side 
    real*8, dimension(:, :), allocatable :: term_1, term_2, s_N_trans, rep_s_N,x,y,xtot,rot_gauss_pts
    real*8, dimension(:, :), allocatable :: x_rot, y_rot, nodes_x, nodes_y, weights, weights1d
    integer :: number_rows, number_cols, sum_num

    !----------------------------------------------------------------------
    ! boundary pts.
    !----------------------------------------------------------------------
    x_bd = polygon_sides(:,1)
    y_bd = polygon_sides(:,2)

    !----------------------------------------------------------------------
    ! "minimum" rectangle containing polygon.
    !----------------------------------------------------------------------
    x_min = minval(x_bd); x_max = maxval(x_bd)
    y_min = minval(y_bd); y_max = maxval(y_bd)



    !--------------------------------------------------------------------------
    ! polygon rotation (if necessary).
    !--------------------------------------------------------------------------
    if (.not. present(rotation) ) then ! no rotation!
       rot_matrix = reshape( (/ 1.0d0, 0.0d0, 0.0d0, 1.0d0 /), (/ 2, 2 /) )
       axis_abscissa = (/x_min, y_max /)-(/ x_min, y_min /)
    else

       select case (rotation)

       case (1) ! automatic

          ! bullet proofing ...
          if ( allocated(P) ) deallocate(P)
          if ( allocated(Q) ) deallocate(Q)

          call auto_rotation(polygon_sides, P, Q, polygon_bd_rot &
               , rot_matrix, rot_angle, axis_abscissa)
          polygon_sides = polygon_bd_rot ! updatinf the polygon frame (sides) with the rotated one

       case (2)
          ! bullet proofing ...
          if ( (.not. allocated(P)) .or. (.not. allocated(Q)) ) then
             print *, 'P and Q must be allocated and initialized in order to do manual rotation! stop'
             stop
          end if
          nrm_vect = norm(Q-P)
          if ( nrm_vect > 0.0d0 ) then
             ! direction_axis = (Q-P) / nrm_vect;
             call auto_rotation(polygon_sides, P, Q, polygon_bd_rot &
                  , rot_matrix, rot_angle, axis_abscissa)
             polygon_sides = polygon_bd_rot ! updatinf the polygon frame (sides) with the rotated one

          else
             print *, 'norm(Q-P) must be greater than zero in order to do manual rotation! stop'
             stop

          end if

       case default

          print *, 'unrecognized option in polygon_gauss_leg(...) subroutine! stop'
          stop

       end select

    end if

    ! use danielle funaro's subroutines to compute 1D gauss-legendre rule
    call cubature_rules_1D((N-1), s_N, w_N)
    N_length = size(s_N)

    M = N + 1
    call cubature_rules_1D((M-1), s_M, w_M)

    L = size(polygon_sides, 1) - 1

    a = axis_abscissa(1)

    !----------------------------------------------------------------------
    ! ready to compute 2d nodes (nodes_x,nodes_y) and weights
    !----------------------------------------------------------------------
    do index_side = 1, L
        x1 = polygon_sides(index_side, 1); x2 = polygon_sides((index_side+1), 1)
        y1 = polygon_sides(index_side, 2); y2 = polygon_sides((index_side+1), 2)
        if (.not. ((x1 .eq. a) .and. (x2 .eq. a)) ) then
            if ((y2-y1) .ne. 0.0d0) then

                if ((x2-x1) .ne. 0.0d0) then
                    s_M_loc => s_M
                    w_M_loc => w_M
                else
                    s_M_loc => s_N
                    w_M_loc => w_N
                 end if

                M_length = size(s_M_loc)

                half_pt_x=(x1+x2)/2.0d0; half_pt_y=(y1+y2)/2.0d0
                half_length_x=(x2-x1)/2.0d0; half_length_y=(y2-y1)/2.0d0

                ! gaussian points on the side.
                if ( allocated(x_gauss_side)) deallocate(x_gauss_side)
                allocate(x_gauss_side(M_length))
                if ( allocated(y_gauss_side)) deallocate(y_gauss_side)
                allocate(y_gauss_side(M_length, 1))

                x_gauss_side=half_pt_x+half_length_x*s_M_loc 
                y_gauss_side(:, 1)=half_pt_y+half_length_y*s_M_loc

                if ( allocated(scaling_fact_plus)) deallocate(scaling_fact_plus)
                allocate(scaling_fact_plus(M_length, 1))
                if ( allocated(scaling_fact_minus)) deallocate(scaling_fact_minus)
                allocate(scaling_fact_minus(M_length, 1))

                scaling_fact_plus(:,1)=(x_gauss_side+a)/2.0d0 
                scaling_fact_minus(:,1)=(x_gauss_side-a)/2.0d0

                if ( allocated(local_weights)) deallocate(local_weights)
                allocate(local_weights(M_length))
                local_weights=(half_length_y*scaling_fact_minus(:,1))*w_M_loc

                if ( allocated(term_1)) deallocate(term_1)
                allocate(term_1(M_length, N_length))

                term_1=repmat(scaling_fact_plus,1,N_length)

                if ( allocated(term_2)) deallocate(term_2)
                allocate(term_2(M_length, N_length))

                term_2=repmat(scaling_fact_minus,1,N_length)

                if ( allocated(s_N_trans) ) deallocate(s_N_trans)
                allocate(s_N_trans(1, size(s_N)))
                s_N_trans(1, :) = s_N

                if ( allocated(rep_s_N) ) deallocate(rep_s_N)
                allocate(rep_s_N(M_length, size(s_N)))

                rep_s_N=repmat(s_N_trans,M_length,1)

                if ( allocated(x)) deallocate(x)
                allocate(x(size(term_1, 1), size(term_1, 2)))
 
                x=term_1+term_2 * rep_s_N

                if ( allocated(y)) deallocate(y)
                allocate(y(size(term_1, 1), size(term_1, 2)))

                y=repmat(y_gauss_side, 1 , N_length)

                number_rows=size(x,1)
                number_cols=size(x,2)

                call  mat2crs(x, xx)
                call  mat2crs(y, yy)

                if (allocated(xtot) ) deallocate(xtot)
                allocate(xtot(2, size(xx)))

                xtot(1, :) = xx
                xtot(2, :) = yy

                if (allocated(rot_gauss_pts)) deallocate(rot_gauss_pts)
                allocate(rot_gauss_pts(2, size(xx)))

                rot_gauss_pts=matmul(transpose(rot_matrix), xtot) 

                if ( allocated(x_rot0)) deallocate(x_rot0)
                allocate(x_rot0(size(rot_gauss_pts, 2)))
                if ( allocated(y_rot0)) deallocate(y_rot0)
                allocate(y_rot0(size(rot_gauss_pts, 2)))

                x_rot0=rot_gauss_pts(1,:)
                y_rot0=rot_gauss_pts(2,:)

                if ( allocated(x_rot)) deallocate(x_rot)
                allocate(x_rot(number_rows,number_cols))
                if ( allocated(y_rot)) deallocate(y_rot)
                allocate(y_rot(number_rows,number_cols))


                x_rot=reshape(x_rot0, (/ number_rows,number_cols /) )
                y_rot=reshape(y_rot0, (/ number_rows,number_cols /) )

                call add_double(nodes_x, x_rot)
                call add_double(nodes_y, y_rot)
                call add_double(weights1d, reshape(local_weights, (/ size(local_weights), 1 /) ) )

             end if
          end if
       end do
       
       if(allocated(weights)) deallocate(weights)
       allocate(weights(size(weights1d,1), size(w_N)))
       weights=matmul(weights1d, reshape(w_N, (/ 1, size(w_N) /) ) )


       sum_num = size(weights)
       allocate(xyw(sum_num, 3))
       xyw(:, 1) = reshape(nodes_x, (/ sum_num /))
       xyw(:, 2) = reshape(nodes_y, (/ sum_num /))  
       xyw(:, 3) = reshape(weights, (/ sum_num /))  

       ! major cleanup
       !
       if(allocated(s_N)) deallocate(s_N)
       if(allocated(w_N)) deallocate(w_N)
       if(allocated(s_M)) deallocate(s_M)
       if(allocated(w_M)) deallocate(w_M)
       if(allocated(x_gauss_side)) deallocate(x_gauss_side)
       if(allocated(local_weights)) deallocate(local_weights)
       if(allocated(xx)) deallocate(xx)
       if(allocated(yy)) deallocate(yy)
       if(allocated(x_rot0)) deallocate(x_rot0)
       if(allocated(y_rot0)) deallocate(y_rot0)
       if(allocated(scaling_fact_plus)) deallocate(scaling_fact_plus)
       if(allocated(scaling_fact_minus)) deallocate(scaling_fact_minus)
       if(allocated(y_gauss_side)) deallocate(y_gauss_side) 
       if(allocated(term_1)) deallocate(term_1)
       if(allocated(term_2)) deallocate(term_2)
       if(allocated(s_N_trans)) deallocate(s_N_trans)
       if(allocated(rep_s_N)) deallocate(rep_s_N)
       if(allocated(x)) deallocate(x)
       if(allocated(y)) deallocate(y)
       if(allocated(xtot)) deallocate(xtot)
       if(allocated(rot_gauss_pts)) deallocate(rot_gauss_pts)
       if(allocated(x_rot)) deallocate(x_rot)
       if(allocated(y_rot)) deallocate(y_rot)
       if(allocated(nodes_x)) deallocate(nodes_x)
       if(allocated(nodes_y)) deallocate(nodes_y)
       if(allocated(weights)) deallocate(weights)
       if(allocated(weights1d)) deallocate(weights1d)

    ! done here
  end subroutine polygon_gauss_leg

  subroutine mat2crs(A, b)
    implicit none
    real*8, dimension(:, :), intent(in) :: A
    real*8, dimension(:), allocatable :: b

    ! local vars
    integer :: m, n, j, i1, i2

    m = size(A, 1); n = size(A, 2)

    if ( allocated(b) ) deallocate(b)

    allocate(b(m * n))

    do j = 1, n

       i1 = (j-1) * m + 1
       i2 = j * m

       b(i1:i2) = A(:, j)

    end do

    ! done here
  end subroutine mat2crs

  ! adds double matrix "x" to the bottom of matrix "a"
  ! they should be compatible (same number of columns)
  !
  subroutine add_double(a, x)
    implicit none
    real*8, dimension(:, :), allocatable :: a
    real*8, dimension(:, :), intent(in) :: x

    ! local vars
    integer :: m, n
    real*8, dimension(:, :), allocatable :: tmp

    if ( .not. allocated(a) ) then ! first timer!

       allocate(a(size(x,1), size(x,2)))
       a = x

    else

       m = size(a, 1) + size(x, 1)
       n = size(a, 2)

       ! bullet proofing
       if ( n .ne. size(x, 2) ) then
          print *, 'matrices are not compatible to be concatenated! stop'
          stop
       end if

       allocate(tmp(m , n))
       if ( size(a, 1) > 0 ) tmp(1:size(a, 1), 1:n) = a(1:size(a, 1), 1:n)
       tmp((size(a, 1)+1):m, 1:n) = x
       call move_alloc(tmp, a)

    end if

    ! cleanup for safety (compiler dependent!) 
    if (allocated(tmp)) deallocate(tmp)

    ! done here
  end subroutine add_double

  ! automatic rotation of a convex polygon so that "gaussian points",
  ! as in the paper they are all contained in the convex polygon.
  subroutine auto_rotation(polygon_bd, vertex_1, vertex_2, polygon_bd_rot &
       , rot_matrix, rot_angle, axis_abscissa)
    implicit none
    real*8, dimension(:, :), intent(in) :: polygon_bd
    real*8, dimension(:), allocatable :: vertex_1, vertex_2
    real*8, dimension(:, :), intent(out) :: polygon_bd_rot
    real*8, dimension(2, 2), intent(out) :: rot_matrix
    real*8 :: rot_angle
    real*8, dimension(2), intent(out) :: axis_abscissa

    ! local vars
    real*8, dimension(:, :), allocatable :: distances
    real*8, dimension(size(polygon_bd, 1)) :: max_distances
    integer, dimension(size(polygon_bd, 1)) :: max_col_comp
    real*8 :: max_distance
    integer :: max_row_comp
    real*8, dimension(size(polygon_bd, 2)) :: direction_axis
    real*8 :: rot_angle_x, rot_angle_y
    integer :: number_sides

    ! find direction and rotation angle.
    if (.not. allocated(vertex_1)) then

       ! computing all the distances between points.a little time consuming
       ! as procedure.
       allocate(distances(size(polygon_bd, 1), size(polygon_bd, 1)))
       distances = points2distances(polygon_bd)
       max_distances = maxval(distances, 2)
       max_col_comp = maxloc(distances, 2)
       max_distance = maxval(max_distances, 1)
       max_row_comp = maxloc(max_distances, 1)
       allocate(vertex_1(size(polygon_bd, 2)), vertex_2(size(polygon_bd, 2)) )

       vertex_1=polygon_bd(max_col_comp(max_row_comp),:)
       vertex_2=polygon_bd(max_row_comp,:)

       direction_axis=(vertex_2-vertex_1)/max_distance

    else
       direction_axis=(vertex_2-vertex_1)/norm(vertex_2-vertex_1)
    end if

    rot_angle_x=acos(direction_axis(1))
    rot_angle_y=acos(direction_axis(2))

    if (rot_angle_y <= PI/2.0d0) then
       if (rot_angle_x <= PI/2.0d0) then
          rot_angle = -rot_angle_y
       else
          rot_angle = rot_angle_y
       end if
    else
       if (rot_angle_x <= PI/2.0d0) then
          rot_angle = PI - rot_angle_y
       else
          rot_angle = rot_angle_y
       end if
    end if

    ! CLOCKWISE ROTATION.
    rot_matrix = reshape( (/ cos(rot_angle), -sin(rot_angle)  &
                           , sin(rot_angle),  cos(rot_angle) /), (/ 2, 2/) )

    number_sides = size(polygon_bd, 1) - 1

    polygon_bd_rot = transpose( matmul(rot_matrix, transpose(polygon_bd)) )

    axis_abscissa = matmul(rot_matrix, vertex_1)

    ! done here
  end subroutine auto_rotation


  ! is only limited to Gauss-Legendre quadrature rule
  subroutine cubature_rules_1D(n, nodes, weights)
    implicit none
    integer, intent(in) :: n
    real*8, dimension(:), allocatable :: nodes, weights

    ! local vars
    integer :: r
    real*8, parameter :: alpha = 0.0d0, beta  = 0.0d0
    real*8, dimension(:), allocatable :: derjac

    ! bullet proofing
    if (allocated(nodes)) deallocate(nodes)
    if (allocated(weights)) deallocate(weights)

    ! specific way!
    r = n + 1

    ! computing Legendre-Gauss-Jacobi points for integration
    ! and corresponding weight functions
    allocate(derjac(r), nodes(r), weights(r))
    call ZEJAGA(r, alpha, beta, nodes, derjac)
    call WEJAGA(r, alpha, beta, nodes, derjac, weights)

    ! clean ups
    deallocate(derjac)

    ! done here
  end subroutine cubature_rules_1D

  ! Create euclidean distance matrix from point matrix.
  function points2distances(points) result(distances)
    implicit none
    real*8, dimension(:, :), intent(in) :: points
    real*8, dimension(:, :), allocatable :: distances

    ! local vars
    integer :: numpoints,dim
    real*8, dimension(:, :), allocatable :: lsq

    ! Get dimensions.
    numpoints = size(points, 1)
    dim = size(points, 2)

    ! All inner products between points.
    allocate(distances(numpoints, numpoints))
    distances = matmul(points, transpose(points))

    ! Vector of squares of norms of points.
    allocate(lsq(numpoints, 1))
    lsq = diag(distances)

    ! Distance matrix.
    distances = sqrt(repmat(lsq,1,numpoints) + &
         transpose(repmat(lsq,1,numpoints)) - &
         2.0d0 * distances)

    ! clean ups
    if ( allocated(lsq) ) deallocate(lsq)

    ! done here
  end function points2distances

  ! similar to MATLAB's diag
  function diag(A)
    implicit none
    real*8, dimension(:, :), intent(in) :: A
    real*8, dimension(:, :), allocatable :: diag

    ! local vars
    integer :: i

    if (size(A, 1) .ne. size(A, 2) ) then
       print *, 'the matrix should be square in order to get the diagonal! stop'
       stop
    end if

    allocate(diag(size(A, 1), 1))

    do i = 1, size(A, 1)
       diag(i, 1) = A(i, i)
    end do

    ! done here
  end function diag

  ! similar functionality to MATLAB's repmat
  !
  ! makes an "mxn" prationed matrix where
  ! each partition is simply a copy of matrix "A" 
  ! 
  function repmat(A, m, n)
    implicit none
    real*8, dimension(:, :), intent(in) :: A
    integer, intent(in) :: m, n
    real*8, dimension(:, :), allocatable :: repmat

    ! local vars
    integer :: Am, An
    integer :: i, j, i1, i2, j1, j2

    Am = size(A, 1)
    An = size(A, 2)
    allocate( repmat(Am * m , An * n) )

    do i = 1, m
       do j = 1, n
          i1 = (i - 1) * Am + 1
          i2 = i * Am
          j1 = (j - 1) * An + 1
          j2 = j * An
          repmat(i1:i2, j1:j2) = A
       end do
    end do

    ! done here
  end function repmat

  function norm(A)
    implicit none
    real*8, dimension(:), intent(in) :: A
    real*8 :: norm

    norm = sqrt(sum(A**2))

    ! done here
  end function norm

  ! writes the quadrature points to the output file
  ! to be used for visualization with python script viselem.py
  ! 
  subroutine write_quad_to_file(xyw, filename)
    implicit none
    real*8, dimension(:, :), intent(in) :: xyw
    character(len=*), intent(in) :: filename

    ! local vars
    integer :: i, j

    ! open the output file first
    open(unit = 10, file = filename)


    do i = 1, size(xyw, 1)

       do j = 1, size(xyw, 2)
          write(10, '(F30.17, A)', advance = 'no') xyw(i, j), ' '
       end do

       write(10,*) ! new line

    end do

    ! close the file
    close(10)

    ! done here
  end subroutine write_quad_to_file

  ! generates some benchmark functions for
  ! validation of the gauss like quadrature
  ! for polygons
  function test_func(x, y, bench)
    implicit none
    real*8, intent(in) :: x, y
    integer, intent(in) :: bench
    real*8 :: test_func

    select case (bench)
    case (4)

       test_func = exp((x-0.5d0)**2+(y-0.5d0)**2)

    case (6)

       test_func = cos(30.0d0 * (x+y))

    case default

       print *, 'unknown benchmark test function! stop'
       stop

    end select

    ! done here
  end function test_func

  ! generates a benchmarking polygon numbered "polynum"
  subroutine get_polygon(polynum, polygon_bd)
    implicit none
    integer, intent(in) :: polynum
    real*8, dimension(:, :), allocatable :: polygon_bd

    ! bullet proofing 
    if(allocated(polygon_bd)) deallocate(polygon_bd)

    select case (polynum)

    case (1)

       allocate(polygon_bd(5, 2))
       polygon_bd = reshape( (/ 0.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0 &
            , 0.0d0, 0.0d0, 1.0d0, 1.0d0, 0.0d0 /), (/ 5, 2 /) )

    case (2)

       allocate(polygon_bd(7, 2))
       polygon_bd = reshape( (/ 0.1d0, 0.7d0, 1.0d0, 0.75d0, 0.5d0, 0.0d0, 0.1d0 &
            , 0.0d0, 0.2d0, 0.5d0, 0.85d0, 1.0d0, 0.25d0, 0.0d0 /), (/7, 2 /) )

    case (3)

       allocate(polygon_bd(10, 2))
       polygon_bd = reshape( (/ 0.25d0, 0.75d0, 0.75d0, 1.0d0, 0.75d0, 0.75d0, 0.5d0, 0.0d0, 0.25d0, 0.25d0 &
            , 0.0d0, 0.5d0, 0.0d0, 0.5d0, 0.75d0, 0.85d0, 1.0d0, 0.75d0, 0.5d0, 0.0d0 /), (/10, 2 /) )

    case default

       print *, 'unknown polygon number! stop'
       stop

    end select

    ! done here
  end subroutine get_polygon

  ! some exact definite integrals over the
  ! given benchmark polygons for the given
  ! benchmark test functions
  function exact_result(polynum, bench)
    implicit none
    integer, intent(in) :: polynum, bench
    real*8 :: exact_result

    select case(polynum)

    case (1)

       select case (bench)
       case (4)
          exact_result=1.188043774905800D0 
       case (6)
          exact_result=0.000289906533545D0
       case default
          print *, 'unknown becnhmarking function! stop'
          stop
       end select

    case (2)

       select case (bench)
       case (4)
          exact_result=0.593459365720563D0
       case (6)
          exact_result=0.008421180941490D0
       case default
          print *, 'unknown becnhmarking function! stop'
          stop
       end select

    case (3)

       select case (bench)
       case (4)
          exact_result=0.531841554503018D0
       case (6)
          exact_result=0.014222050981512D0
       case default
          print *, 'unknown becnhmarking function! stop'
          stop
       end select

    case default

       print *, 'unknown polygon number! stop'
       stop

    end select

    ! done here
  end function exact_result

  ! evaluates the numerical integral of order "N" over the 
  ! given polygon with number "polynum" for the given
  ! benchmark function "bench" and computes the absolute error "abs_err".
  ! An option can be set for automatic (or no) rotation of principle axis
  ! in calculating quadrature points and weights.  
  subroutine eval_bench_integ(N, polynum, bench, abs_err, option, ngauss, tofile)
    implicit none
    integer, intent(in) :: N, polynum, bench
    real*8, intent(out) :: abs_err
    integer, intent(in), optional :: option
    integer, intent(out), optional :: ngauss
    integer, intent(in), optional :: tofile

    ! local vars
    integer :: k
    real*8 :: integ, tmp
    real*8, dimension(:), allocatable :: P, Q
    real*8, dimension(:, :), allocatable :: polygon_bd, xyw

    ! get the requested polygon
    call get_polygon(polynum, polygon_bd)

    ! get the gauss-legendre like quadrature rule for the given polygon
    if ( present(option) ) then
       call polygon_gauss_leg(N, polygon_bd, P, Q, xyw, option)
    else
       call polygon_gauss_leg(N, polygon_bd, P, Q, xyw)
    end if

    ! return the number of generated gausslike points (if requested)
    if ( present(ngauss) ) ngauss = size(xyw, 1)

    ! further write to file (if requested)
    if ( present(tofile) ) then
       call write_quad_to_file(xyw, 'eval_bench_integ_output.dat')
    end if

    ! evaluate the numerical integral of the selected benchmark function
    integ = 0.0d0
    do k = 1, size(xyw, 1)
       tmp = test_func(xyw(k, 1), xyw(k, 2), bench)
       integ = integ + tmp * xyw(k, 3)
    end do

    ! finally, find the absolute value of the error
    abs_err = abs(integ - exact_result(polynum, bench))  


    ! clean ups
    if(allocated(P)) deallocate(P)
    if(allocated(Q)) deallocate(Q)
    if(allocated(polygon_bd)) deallocate(polygon_bd)
    if(allocated(xyw)) deallocate(xyw)

    ! done here
  end subroutine eval_bench_integ

end module polygauss

! program tester
!   use dispmodule
!   use polygauss
!   implicit none

!   integer :: n
!   real*8, dimension(:, :), allocatable :: A, full, B
!   real*8, dimension(:), allocatable :: nodes, weights
!   real*8, dimension(:, :), allocatable :: polygon_bd, polygon_bd_rot
!   real*8, dimension(:), allocatable :: vertex_1, vertex_2
!   real*8, dimension(2, 2) :: rot_matrix
!   real*8 :: rot_angle
!   real*8, dimension(2) :: axis_abscissa
!   real*8, dimension(:), allocatable :: P, Q
!   real*8, dimension(:, :), allocatable :: xyw

!   integer :: i, ngauss0
!   real*8, dimension(:), allocatable :: error_i

!   ! --------- test repmat ---------
!   allocate(A(5, 2), full(10, 4))
!   A = reshape( (/ 0.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0 &
!        , 0.0d0, 0.0d0, 1.0d0, 1.0d0, 0.0d0 /), (/ 5, 2 /) )
!   full = repmat(A, 2, 2)

!   ! show
!   call disp(A)
!   call disp(full)

!   ! --------- test points2distances ---------
!   allocate(B(5,5))
!   B = points2distances(A)

!   call disp('B = ', B)

!   ! --------- test cubature rule---------
!   n = 9
!   call cubature_rules_1D(n, nodes, weights)

!   call disp('nodes = ', nodes, DIGMAX = 15)
!   call disp('weights = ', weights, DIGMAX = 15)


!   ! --------- test auto rotation ---------
!   allocate( polygon_bd (5, 2), polygon_bd_rot(5, 2) )
!   polygon_bd = reshape( (/ 0.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0 &
!                          , 0.0d0, 0.0d0, 1.0d0, 1.0d0, 0.0d0 /), (/5, 2 /) )

!   call auto_rotation(polygon_bd, vertex_1, vertex_2, polygon_bd_rot &
!        , rot_matrix, rot_angle, axis_abscissa)

!   call disp('polygon_bd = ', polygon_bd)
!   call disp('vertex_1 = ', vertex_1)
!   call disp('vertex_2 = ', vertex_2)
!   call disp('polygon_bd_rot = ', polygon_bd_rot)
!   call disp('rot_matrix = ', rot_matrix)
!   print *, 'rot_angle = ', rot_angle 
!   print *, 'axis_abscissa = ', axis_abscissa 

!   !
!   ! --------- test final gauss-legendre quad rule ---------
!   !
!   ! simple polygon (quad elem) with automatic rotation
!   polygon_bd = reshape( (/ 0.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0 &
!                          , 0.0d0, 0.0d0, 1.0d0, 1.0d0, 0.0d0 /), (/5, 2 /) )

!   call polygon_gauss_leg(10, polygon_bd, P, Q, xyw, 1)
!   call disp('xyw_quad = ', xyw)
!   ! further write to file
!   call write_quad_to_file(xyw, 'xyw_quad.dat')
!   deallocate(polygon_bd, P, Q, xyw)

!   ! a benchmark convex polygon with automatic rotation
!   allocate(polygon_bd(7, 2))
!   polygon_bd = reshape( (/ 0.1d0, 0.7d0, 1.0d0, 0.75d0, 0.5d0, 0.0d0, 0.1d0 &
!        , 0.0d0, 0.2d0, 0.5d0, 0.85d0, 1.0d0, 0.25d0, 0.0d0 /), (/7, 2 /) )

!   call polygon_gauss_leg(50, polygon_bd, P, Q, xyw, 1)
!   call disp('xyw_convex = ', xyw)
!   ! further write to file
!   call write_quad_to_file(xyw, 'xyw_convex.dat')
!   deallocate(polygon_bd, P, Q, xyw)

!   ! compute the error
!   allocate(error_i(50))
!   do i = 1, size(error_i) 
!      call eval_bench_integ(N = i, polynum = 2, bench = 4 &
!           , abs_err = error_i(i), option = 1, ngauss = ngauss0, tofile = 1)
!      print *, i, ngauss0, error_i(i)
!   end do



!   ! clean ups
!   if ( allocated(A) ) deallocate(A)
!   if ( allocated(full) ) deallocate(full)
!   if ( allocated(B) ) deallocate(B)
!   if ( allocated(nodes) ) deallocate(nodes)
!   if ( allocated(weights) ) deallocate(weights)
!   if ( allocated(polygon_bd) ) deallocate(polygon_bd)
!   if ( allocated(polygon_bd_rot) ) deallocate(polygon_bd_rot)
!   if ( allocated(vertex_1) ) deallocate(vertex_1)
!   if ( allocated(vertex_2) ) deallocate(vertex_2)
!   if ( allocated(P)) deallocate(P)
!   if ( allocated(Q)) deallocate(Q)
!   if ( allocated(xyw)) deallocate(xyw)
!   if ( allocated(error_i)) deallocate(error_i)

!   ! done here
! end program tester
