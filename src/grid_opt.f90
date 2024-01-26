!> @ingroup Terrible
!> @author Arash Ghasemi and Mohammad Mahtabi
!> @brief
!> A Module to handle grid-related operations
!> @details
!> This modules includes the type of grid object as well as
!> subroutines and functions to
!> handle grid-related operations such as reading, writing the
!> output, etc.
!> @see elem_opt, feket, ncc_triangle, spline, fem_reordering
!> quadri_elem, and approx_fekete modules.

module grid_opt
  use fekete
  use ncc_triangle
  use spline
  use fem_reordering
  use quadri_elem
  use approx_fekete, only : fekete_table

  ! use ifport

  ! performs grid operations
  !    read
  !    write
  !    manipulation
  ! for unstructered grid in TETREX format
  implicit none

  private

  !> @brief
  !> Data types for master element
  !> @details
  !> Includes arrays of Jacobi polynomial constants \$f \xi and \eta \$f for
  !> the nodes on the master element.
  type master_elem
    real*8, dimension(:), allocatable :: xi  !< Jacobi polynomial constant \$f \xi \$f
    real*8, dimension(:), allocatable :: eta  !< Jacobi polynomial constant \$f \eta \$f
  end type master_elem

  !> @brief
  !> Data types for edges
  !> @details
  !> Includes arrays of edges which is the number of the edges.
  type edges
    integer, dimension(:) , allocatable :: edg1, edg2, edg3, edg4
  end type edges

  !> @brief
  !> Data types for curves
  !> @details
  !> Includes arrays of edges which is the number of the edges.
  type curve
    ! type of interpolation curve use
    integer :: btype = spline_i
    !    integer :: btype = nurbs_i
    ! point coords
    real*8, dimension(:), allocatable :: x, y
    ! spline stuff
    real*8, dimension(:), allocatable :: t, Mx, My, a, b, c, d, cp, dp
  end type curve

  !> @brief
  !> Data types for the grid object, i.e. the whole mesh.
  !> @details
  !> Includes various attributes for the grid object
  !> such as: number of total cells (elements), nodes, element connectivity, etc.
  !> neqs (# of equations considered for the element)
  !> \f$ \alpha \f$ and \f$ \beta \f$ (Jacobi polynomial coefficients)
  !> \f$ \xi \f$ and \f$ \eta \f$
  !> init_lock (showing whether and element is initiated or not)
  type grid
    integer :: nnodesg,ncellsg,nbedgeg,ntri,nquad4
    real(kind(0.d0)), allocatable :: x(:),y(:)
    integer, allocatable :: icon(:,:)  !< Element connectivity array. First dimension indicates the &
      !!element number and the second lists the related node numbers
    integer, allocatable :: ibedge(:,:)  !< Edge connectivity array. First dimension indicates the &
      !!edge number and the second lists the related node numbers
    integer, allocatable :: ibedgeBC(:)  !< Boundary condition values
    integer, allocatable :: ibedgeELEM(:)  !< Elements having boundary edges
    integer, allocatable :: ibedgeELEM_local_edg(:) !< Local edges of the element having boundary edges
    character*40 meshFile  !< Name of the mesh file.
    integer, allocatable :: e2e(:,:) !< @todo Add description here!

    ! hook-ups for higher order elements
    integer :: tot_bn_seg !< @todo Add description here!
    integer :: tot_repeated_bn_nodes !< @todo Add description here!
    integer, allocatable :: dup_nodes(:) !< @todo Add description here!
    type(curve), dimension(:), allocatable :: bn_curves !< @todo Add description here!
    ! el2bn(1:ncellsg, 1:2) = (/ tag , edgnum /)
    ! tag == 0 (interior element) : no need for BC
    ! tag <> 0 (boundary element) : contains ibedgeBC(edgnum)
    ! edgnum = the global edge number of the edge in this triangle which is
    !          is on the curved boundary.
    integer, dimension(:,:), allocatable :: el2bn !< Element to boundary map
    type(edges), dimension(:), allocatable :: el2edg !< Element to boundary edge

    integer, dimension(:), allocatable :: p, eltype, npe, elname, skleton_pts
    integer :: nnodesg0, ncellsg0, nbedgeg0 !initial config before hp-adapt

    type(master_elem), dimension(:), allocatable :: maselem

    character(len = 2) :: galtype ! Galerkin type, = DG .or. PG

    type(fekete_table) :: tfekete_table ! generic fekete points and quadratures

    ! (1:ncellsg, 1:4) includes tri and quad simultanously!
    integer, dimension(:, :), allocatable :: local_edg_bc

    logical :: linear_boundaries = .false.

  end type grid

  ! enumerate types for element names grd%elname
  !
  integer, parameter, public :: GEN_TRIANGLE = 1, GEN_QUADRI = 2, GEN_SPHULL = 3
  real*8, parameter, public :: piecewise_tol = 1.0d-9

  ! public data structure
  public :: grid, curve

  ! public subroutines
  public :: read_grid_TETREX, print_grid_props, write_u_tecplot
  public :: write_u_cell_centered_tecplot
  public :: search_tri
  public :: read_segment_file, visual_curved_bonds
  public :: p2node, add_more_points
  public :: scatter_tecplot
  public :: fill_ibedgeELEM_local_edg
  public :: find_dup_nodes, find_t
  public :: fill_local_edg_bc
  public :: comp_quad4_area
  public :: snap2curve

contains

  ! reads the input grid generated by CAD software (like Pointwise).
  ! the CAD grid written in TETREX format. The grid should be 2D
  ! and should be either triangle or quad.

  subroutine read_grid_TETREX(in_grd, grd, noe2e)
    implicit none
    character(len=*), intent(in)      :: in_grd
    type(grid), intent(inout) :: grd
    logical, intent(in), optional :: noe2e


    ! local vars
    integer :: i, ibedgecell, lface, npe, nzone
    integer :: tmp1,tmp2,cell_type
    integer :: btype,cell_indx

    !  =======       initialization         =========
    ! data initialization ---
    grd%meshFile= in_grd
    ! hard reset
    ! =-----------------=
    grd%nnodesg=0; grd%ncellsg=0; grd%nbedgeg=0; grd%ntri=0; grd%nquad4=0

    i = 0; ibedgecell = 0; lface = 0; npe = 0; nzone = 0
    tmp1 = 0; tmp2 = 0; cell_type = 0
    btype=0; cell_indx = 0
    ! =-----------------=

    ! formatting stuff
    1   format(4(/), 4x, 2(I9),(/))
    2   format(10x, 2(I10))

    write(*,'(A,/)') "READING/Checking the mesh file ..."

    open(10, file = grd%meshFile, status = 'old')

    read(10,1) nzone, grd%nnodesg

    do i=1,nzone
      read(10,2) tmp1, tmp2
      grd%ncellsg = grd%ncellsg + tmp1
      grd%nbedgeg = grd%nbedgeg + tmp2
    enddo

    read(10,*)  !locate the cursor to cell to point map

    do i = 1, grd%ncellsg
      read(10,*) cell_indx, cell_type
      if (cell_type==6) then
        grd%ntri = grd%ntri + 1         !Triangle cell
      else if (cell_type==5) then
        grd%nquad4 = grd%nquad4 + 1       !Quadrilateral cell
      else
        write(*,'(A)') "Error: undefined cell tpye at i=",i
        stop
      endif
    enddo
    if( (grd%ntri .ne. 0) .and. (grd%nquad4 .ne. 0) ) then
      stop 'can not handle mixed element mesh read in this version... stopping'
    endif

    do i = 1, (grd%nnodesg + 2) !locate the cursor to bface to point map
      read(10,*)
    enddo

    !Check the boundariy types: ---------------------------------
    do i = 1, grd%nbedgeg
      read(10,*) btype, cell_indx, lface

      !Check the local faces: ---------------------------------
      if ((lface < 1) .or. (lface>4) ) then
        write(*,'(A)') "Error: invalid local face number for a 2D mesh"
        stop
      endif
    enddo

    close(10)
    !----------------------------------------------------------------

    write(*,'(A)')      "Mesh Statistics:"
    write(*,'(A)')      "----------------"
    write(*,'(A,I7)')   "-Number of Zones:         ", nzone
    write(*,'(A,I7)')   "-Number of Points:        ", grd%nnodesg
    write(*,'(A,I7)')   "-Number of Cells:         ", grd%ncellsg
    write(*,'(A,I7)')   "-Number of Boundary faces:", grd%nbedgeg
    write(*,'(A,I7)')   "-Number of Tri-cells:     ", grd%ntri
    write(*,'(A,I7,/)') "-Number of Quad-cells:    ", grd%nquad4


    !...now allocate to read mesh/connectivity
    if(grd%ntri.ne.0) then
      npe=3
    elseif(grd%nquad4.ne.0) then
      npe=4
    endif
    ALLOCATE(grd%x(grd%nnodesg)) ; ALLOCATE(grd%y(grd%nnodesg))
    ALLOCATE(grd%icon(grd%ncellsg,npe))
    ALLOCATE(grd%ibedge(grd%nbedgeg,2))
    ALLOCATE(grd%ibedgeBC(grd%nbedgeg), grd%ibedgeELEM(grd%nbedgeg))  ! bc will be an integer

    open(10, file = grd%meshFile, status = 'old')

    do i = 1, nzone+7   !locate the cursor
      read(10,*)
    enddo

    do i = 1, grd%ncellsg
      if (grd%ntri.ne.0) then      !Triangle cell
        read(10,*) cell_indx, cell_type, grd%icon(i,1),grd%icon(i,2),grd%icon(i,3)
      else if (grd%nquad4.ne.0) then !Quadrilateral cell
        read(10,*) cell_indx, cell_type, grd%icon(i,1),grd%icon(i,2),grd%icon(i,3),grd%icon(i,4)
      endif
    enddo

    read(10,*)

    !Reading coordinates of points:------------------------------

    do i = 1, grd%nnodesg
      read(10,*) grd%x(i), grd%y(i)
    enddo

    read(10,*)

    do i = 1, grd%nbedgeg
      read(10,*) grd%ibedgeBC(i), ibedgeCell, lface
      grd%ibedgeELEM(i) = ibedgeCell ! store element number
      if (grd%ntri.ne.0) then
        select case(lface)
          case(1)
            grd%ibedge(i,1) = grd%icon(ibedgeCell,1)
            grd%ibedge(i,2) = grd%icon(ibedgeCell,2)
          case(2)
            grd%ibedge(i,1) = grd%icon(ibedgeCell,2)
            grd%ibedge(i,2) = grd%icon(ibedgeCell,3)
          case(3)
            grd%ibedge(i,1) = grd%icon(ibedgeCell,3)
            grd%ibedge(i,2) = grd%icon(ibedgeCell,1)
        endselect
      else if (grd%nquad4.ne.0) then
        select case(lface)
          case(1)
            grd%ibedge(i,1) = grd%icon(ibedgeCell,1)
            grd%ibedge(i,2) = grd%icon(ibedgeCell,2)
          case(2)
            grd%ibedge(i,1) = grd%icon(ibedgeCell,2)
            grd%ibedge(i,2) = grd%icon(ibedgeCell,3)
          case(3)
            grd%ibedge(i,1) = grd%icon(ibedgeCell,3)
            grd%ibedge(i,2) = grd%icon(ibedgeCell,4)
          case(4)
            grd%ibedge(i,1) = grd%icon(ibedgeCell,4)
            grd%ibedge(i,2) = grd%icon(ibedgeCell,1)
        endselect
      endif
    enddo
    !------------------------------------------------------------
    close(10)

    ! create e2e map
    if ( present(noe2e) ) then
      if ( noe2e ) then
        ! do nothing !
      else
        call create_e2e(grd)
      end if
    else
      call create_e2e(grd)
    end if

    ! done here
  end subroutine read_grid_TETREX

  ! prints grid properties and conectivity
  subroutine print_grid_props(grd)
    implicit none
    type(grid), intent(in) :: grd

    ! local vars
    integer :: i, j

    ! print statistics
    print *, 'total number of nodes: ', grd%nnodesg
    print *, 'total number of elements : ', grd%ncellsg
    print *, 'total number of boundary edges: ', grd%nbedgeg
    print *, 'total number of triangles: ', grd%ntri
    print *, 'total number of quads: ', grd%nquad4

    ! print coordinates
    do i = 1, size(grd%x)
      print *, 'point ', i, ' = (', grd%x(i),' , ', grd%y(i), ')'
    end do

    ! print connectivity matrix
    do i = 1, size(grd%icon,1)
      write(*, '(A,I5,A)', advance = 'no') 'element (', i, ') contains nodes: '
      do j = 1, size(grd%icon, 2)
        write(*, '(I5, A)', advance = 'no') grd%icon(i,j), ', '
      end do
      write(*,*)
    end do

    ! print boundary edges info
    do i = 1, size(grd%ibedge,1)
      print *, 'edge ', i, ' = (', grd%ibedge(i,1),',', grd%ibedge(i,2) &
        , ') with BC = ', grd%ibedgeBC(i), 'in the elem num. ' &
        , grd%ibedgeELEM(i) &
        , ' with local edg num. = ', grd%ibedgeELEM_local_edg(i)
    end do

    print *, 'the grid file is located in ', grd%meshFile

    ! print e2e map
    do i = 1, size(grd%e2e,1)
      print *, 'elems connected to triangle :', i,' are : ', grd%e2e(i, :)
    end do

    !
    print *, 'total boundary segments from CAD : ', grd%tot_bn_seg
    print *, 'total repeated boundary nodes : ', grd%tot_repeated_bn_nodes
    ! print *, 'dup_nodes = [', grd%dup_nodes , ']'

    ! print el2bn
    do i = 1, grd%ncellsg
      print *, 'grd%el2bn(', i, ':) = ', grd%el2bn(i, :)
    end do


    ! print el2edg (if allocated!)
    do i = 1, grd%ncellsg
      if ( allocated(grd%el2edg(i)%edg1) .and. &
          allocated(grd%el2edg(i)%edg2) .and. &
          allocated(grd%el2edg(i)%edg3) ) then

        print *, 'in element num.', i , 'we have :'
        print *, 'edg1 = [', grd%el2edg(i)%edg1, ']'
        print *, 'edg2 = [', grd%el2edg(i)%edg2, ']'
        print *, 'edg3 = [', grd%el2edg(i)%edg3, ']'
      else
        print *, 'warning: edg1..3 not allocated for element num. ', i ,'. OK'
      end if
    end do

    ! print master element interpolation points
    do i = 1, grd%ncellsg
      if (allocated(grd%maselem(i)%xi) .and. allocated(grd%maselem(i)%eta)) then
        print *, 'in element num.', i , 'master elem interpolation points are:'
        print *, 'xi  = [', grd%maselem(i)%xi , ']'
        print *, 'eta = [', grd%maselem(i)%eta, ']'
      else
        print *, 'warning: master element is not allocated' &
          , ' for element #', i, ' yet!'
      end if
    end do

    ! done here
  end subroutine print_grid_props

  !< @Detail
  !! Writes the unstructured grid + solution to
  !! Tecplot format.
  !! the format of 'u' is assumed to be:
  !>
  !! u(neqs, nnodesg)
  !!
  subroutine write_u_tecplot(outfile, grd, u, appendit)
    implicit none
    character(len=*), intent(in) :: outfile   !< Name of the outputfile
    type(grid), intent(in) :: grd     !< Grid under analysis
    real*8, dimension(:,:), intent(in) :: u !< Solution matrix
    logical, optional :: appendit   !< A logical variable indicating whether or not to append
    !! the results to the outfile!

    ! local vars
    integer :: i, j, k, neqs, nnodes

    ! init
    neqs = size(u,1)
    nnodes = size(u,2)
    if ( nnodes .ne. grd%nnodesg ) then
      print *, 'nnodes .ne. grd%nnodesg! something is wrong! stop'
      stop
    end if

    ! opening for rewrite
    if ( present( appendit ) ) then
      if ( appendit ) then
        open(10, file = outfile, status="old", position="append", action="write")
      else
        open(10, file = outfile, status = 'unknown')
      end if
    else
      open(10, file = outfile, status = 'unknown')
    end if

    ! write header
    !< @todo We may need to change this title.
    write(10, *) 'title = "spem2d solution"'

    write(10, '(A)', advance = 'no') 'variables = "x", "y"'
    do i = 1, neqs
      write(10, '(A, I1, A)', advance = 'no') ', "u', i,'"'
    end do

    write(10,*) ! new line!


    write(10, '(A, I7, A, I7, A, A)', advance = 'no') 'zone n = ' &
      , nnodes, ', e = ', grd%ncellsg, ', f = fepoint, ' &
      , 'et = quadrilateral'
    write(10,*) ! new line!

    ! write coordinates and values of the vector field [u]
    do k = 1, nnodes

      write(10, '(F30.17, A, F30.17)', advance = 'no') &
        grd%x(k), ' ',  grd%y(k)
      do j = 1, neqs
        write(10, '(A, F30.17)', advance = 'no') ' ',  u(j,k)
      end do

      write(10,*)

    end do

    ! writing the connectivity matrix
    write(10, *)

    if ( grd%nquad4 .eq. 0 ) then ! do old way
      do k = 1, grd%ntri
        write(10, *) ' ',  grd%icon(k,1) &
          , ' ',  grd%icon(k,2), ' ',  grd%icon(k,3) &
          , ' ',  grd%icon(k,3)
      end do

    else

      do k = 1, grd%ncellsg
        select case ( grd%elname(k) )
          case ( GEN_TRIANGLE )
            write(10, *) ' ',  grd%icon(k,1) &
              , ' ',  grd%icon(k,2), ' ',  grd%icon(k,3) &
              , ' ',  grd%icon(k,3)
          case ( GEN_QUADRI )
            write(10, *) ' ',  grd%icon(k,1) &
              , ' ',  grd%icon(k,2), ' ',  grd%icon(k,3) &
              , ' ',  grd%icon(k,4)
          case default
            print *, 'unknown name of element! stop'
            stop
        end select
      end do

    end if

    ! close the output file
    close(10)

    ! done here
  end subroutine write_u_tecplot


  ! writes cell-centered data to tecplot
  ! USE ONLY FOR ONE CELL CENTER (GAUSS POINT)
  ! THAT MEANS TRANSFORM u(neqs, ngauss, ncellsg) to
  !             u(:,:) = u(neqs, 1, ncellsg)
  ! BEFORE THIS!
  subroutine write_u_cell_centered_tecplot(outfile, grd, u, appendit)
    implicit none
    character(len=*), intent(in) :: outfile
    type(grid), intent(in) :: grd
    ! (1..neqs, 1..ncellsg)
    real*8, dimension(:,:), intent(in) :: u
    logical, optional :: appendit

    ! local vars
    integer i, j, jmax, neqs, ncellsg
    ! init
    neqs = size(u,1)
    ncellsg = size(u,2)

    ! error checking
    if (ncellsg .ne. grd%ncellsg) then
      print *,' error in write_u_cell_centered_tecplot(...) : ncellsg of '&
        ,' [u] is not equal to grd%ncellsg!. stop'
      stop
    end if

    ! opening for rewrite
    if ( present( appendit ) ) then
      if ( appendit ) then
        open(10, file = outfile, status="old", position="append", action="write")
      else
        open(10, file = outfile, status = 'unknown')
      end if
    else
      open(10, file = outfile, status = 'unknown')
    end if

    ! write header
    write(10, *) 'title = "spem2d solution CELL CENTERED"'

    write(10, '(A)', advance = 'no') 'variables = "x", "y"'
    do i = 1, neqs
      write(10, '(A, I1, A)', advance = 'no') ', "u', i,'"'
    end do
    write(10,*)

    write(10,*) 'zone T = "this zone"'

    write(10,*) 'STRANDID=1, SOLUTIONTIME=0'

    write(10, '(A,I7,A,I7,A)', advance = 'no')'Nodes=',grd%nnodesg,', Elements=', ncellsg,', ZONETYPE='

    if ( (grd%ntri .ne. 0) .and. (grd%nquad4 .eq. 0) ) then
      write(10, *) 'FETriangle'
      jmax = 3
    elseif ( (grd%ntri .eq. 0) .and. (grd%nquad4 .ne. 0) ) then
      write(10, *) 'FEQuadrilateral'
      jmax = 4
    else
      print *, 'can not export mixed element grid to ' &
        , 'tecplot format at this version. stop.'
      stop
    end if

    write(10,*)'DATAPACKING=BLOCK'
    write(10,'(A,I1,A)') 'VARLOCATION=([3-', (neqs + 2) ,']=CELLCENTERED)'

    ! ! write grid x,y
    ! write(10,*)(grd%x(i),i=1,grd%nnodesg)
    ! write(10,*)(grd%y(i),i=1,grd%nnodesg)

    ! write grid x,y
    do i = 1, grd%nnodesg
      write(10,'(F30.17,A)', advance = 'no') grd%x(i), ' '
      if ( mod(i, 50) == 0) then
        write(10,*)
      end if
    end do

    do i = 1, grd%nnodesg
      write(10,'(F30.17,A)', advance = 'no') grd%y(i), ' '
      if ( mod(i, 50) == 0) then
        write(10,*)
      end if
    end do

    ! write cell centered solution
    do j = 1, neqs
      do i = 1, ncellsg
        write(10,'(F30.17,A)', advance = 'no') u(j,i), ' '
        if ( mod(i, 50) == 0) then
          write(10,*)
        end if
      end do
    end do

    ! writing the connectivity matrix
    ! write(10, *)
    do i = 1, grd%ncellsg
      do j = 1, jmax
        write(10, '(A, I7)', advance = 'no') ' ',  grd%icon(i,j)
      end do
      write(10,*)
    end do

    close(10)

    ! done here
  end subroutine write_u_cell_centered_tecplot

  ! searches triangles in a grid for a target point using
  ! normal dot product method which is a fast method.
  subroutine search_tri(seed, xt, yt, grd, thetri)
    implicit none
    integer :: seed
    real*8, intent(in):: xt, yt
    type(grid), target, intent(in) :: grd
    integer, intent(out) :: thetri

    ! local vars
    integer ::  i,j
    real*8, dimension(2) :: n_hat, r_hat
    real*8, dimension(3) :: dot
    real*8 :: nx,ny, rx,ry, M
    real*8 :: xc, yc, max_dot
    integer :: indx_max_dot
    integer :: pt
    real*8, dimension(:), pointer :: x, y
    integer, dimension(:,:), pointer :: tri, nbr

    integer :: wallcount, rnd_max_dot
    integer, parameter :: default_wallcount = 300

    ! initialization
    x => grd%x(1:grd%nnodesg); y => grd%y(1:grd%nnodesg)
    tri => grd%icon
    nbr => grd%e2e
    wallcount = 0

    !NOTE : seed is the initial triangle
    do  !do it with the current seed
      !computing dot products and filling dot[] for three edges
      do pt = 1, 3 !looping over points in the element

        if (pt .eq. 3) then !selecting two points of an edge with winding
          i = pt
          j = 1
        else
          i = pt
          j = pt + 1
        end if

        !computing normals of that edge
        nx = y(tri(seed,j)) - y(tri(seed,i))
        ny = x(tri(seed,i)) - x(tri(seed,j))
        M = sqrt(nx*nx + ny*ny)
        if( M .eq. 0.0d0 ) then
          print *, 'Fatal error: one side of triangle' &
            , seed, ' has zero length! exit.'
          stop
        end if

        n_hat = (/ (nx / M), (ny / M) /)

        ! finding the center of the edge
        xc = 0.5d0 * (x(tri(seed,i)) + x(tri(seed,j)))
        yc = 0.5d0 * (y(tri(seed,i)) + y(tri(seed,j)))

        ! computing the distance vector for that edge
        rx = xt - xc
        ry = yt - yc
        M = sqrt(rx * rx + ry * ry)
        r_hat = (/ (rx / M), (ry / M) /)

        !removing singularity for the case when seed
        ! target location and center location are EXACTLY same
        if( M .ne. 0.0d0) then !storing the dot product
          dot(i) = n_hat(1)*r_hat(1) + n_hat(2)*r_hat(2)
        else
          dot(i) = -1.0d0
        end if

      end do
      !decision making section
      !the target is inside or on this triangle
      if( (dot(1) <= 0.0d0) .and. (dot(2) <= 0.0d0) .and. (dot(3) <= 0.0d0) ) then
        thetri = seed !found it!
        return
      else !find the maximum dot to define the direction of the marching

        max_dot = dot(1) !initializing max_dot
        indx_max_dot = 1
        do i = 2, 3
          if (dot(i) > max_dot) then
            max_dot = dot(i)
            indx_max_dot = i
          end if
        end do
      end if

      !correct the winding
      if     ( indx_max_dot .eq. 1 ) then
        indx_max_dot = 3
      elseif ( indx_max_dot .eq. 2 ) then
        indx_max_dot = 1
      else
        indx_max_dot = 2
      end if

      ! print *, 'indx_max_dot = ', indx_max_dot

      if (nbr(seed,indx_max_dot) .eq. -1) then
        ! print *, ' warning : we reached to a wall.'
        wallcount = wallcount + 1
        do
          rnd_max_dot = ceiling(rand()*3.0)
          if ( nbr(seed,rnd_max_dot) .ne. -1 ) then
            indx_max_dot = rnd_max_dot
            exit
          end if
        end do
      end if

      if( wallcount > default_wallcount ) then
        print *, ' wall encounter exceeds the default value ', default_wallcount &
          , ' for point (',xt, yt,'). starting using' &
          , ' brute force search ...'
        call search_tri_brute(xt, yt, grd, thetri)
        if ( thetri .eq. -10 ) then
          print *, ' the target point (', xt, yt,') not found' &
            , ' in any triangle even using a brute-force method!!! stop.'
          stop
        end if
        return
      end if

      seed = nbr(seed,indx_max_dot)
      ! print *, 'new seed =', seed

    end do !keep going, you're gonna find it soon!

  end subroutine search_tri

  ! searches triangles in a grid for a target point using
  ! brute-force method
  subroutine search_tri_brute(xt, yt, grd, thetri)
    implicit none
    real*8, intent(in):: xt, yt
    type(grid), target, intent(in) :: grd
    integer, intent(out) :: thetri

    ! local vars
    integer ::  i,j, seed
    real*8, dimension(2) :: n_hat, r_hat
    real*8, dimension(3) :: dot
    real*8 :: nx,ny, rx,ry, M
    real*8 :: xc, yc
    integer :: pt
    real*8, dimension(:), pointer :: x, y
    integer, dimension(:,:), pointer :: tri


    ! initialization
    x => grd%x(1:grd%nnodesg); y => grd%y(1:grd%nnodesg)
    tri => grd%icon

    !NOTE : seed is the current triangle and loops over all tris.
    do  seed = 1, grd%ntri !do it with the current seed
      !computing dot products and filling dot[] for three edges
      do pt = 1, 3 !looping over points in the element

        if (pt .eq. 3) then !selecting two points of an edge with winding
          i = pt
          j = 1
        else
          i = pt
          j = pt + 1
        end if

        !computing normals of that edge
        nx = y(tri(seed,j)) - y(tri(seed,i))
        ny = x(tri(seed,i)) - x(tri(seed,j))
        M = sqrt(nx*nx + ny*ny)
        if( M .eq. 0.0d0 ) then
          print *, 'Fatal error: one side of triangle' &
            , seed, ' has zero length! exit.'
          stop
        end if

        n_hat = (/ (nx / M), (ny / M) /)

        ! finding the center of the edge
        xc = 0.5d0 * (x(tri(seed,i)) + x(tri(seed,j)))
        yc = 0.5d0 * (y(tri(seed,i)) + y(tri(seed,j)))

        ! computing the distance vector for that edge
        rx = xt - xc
        ry = yt - yc
        M = sqrt(rx * rx + ry * ry)
        r_hat = (/ (rx / M), (ry / M) /)

        !removing singularity for the case when seed
        ! target location and center location are EXACTLY same
        if( M .ne. 0.0d0) then !storing the dot product
          dot(i) = n_hat(1)*r_hat(1) + n_hat(2)*r_hat(2)
        else
          dot(i) = -1.0d0
        end if

      end do
      !decision making section
      !the target is inside or on this triangle
      if( (dot(1) <= 0.0d0) .and. (dot(2) <= 0.0d0) .and. (dot(3) <= 0.0d0) ) then
        thetri = seed !found it!
        return
      end if

    end do !keep going, you're gonna find it soon!

    ! Not found yet!!!, oh that is bad, terrible infact :(
    thetri = -10 ! NULL it!

    !done here
  end subroutine search_tri_brute

  ! manually creates e2e map
  ! must be used when grid is read from CAD software
  !
  ! HINT:
  !
  ! there is a builtin e2e gen in Triangle library
  ! so DO NOT use this with that unless the physical element
  ! node placement is reordered in "transform_elem(...)".
  subroutine create_e2e(grd)
    implicit none
    type(grid), intent(inout) :: grd

    ! local vars
    integer :: i, pt1, pt2

    ! ! no mixed element, although it is very
    ! ! easy to extend it :)
    ! if (grd%nquad4 .ne. 0 ) then
    !    print *, 'sorry! create_e2e() is for a full trimesh grid in this version.' &
      !         , ' some quads were found in FEM grid. stop.'
    !    stop
    ! end if

    ! allocate
    if ( allocated(grd%e2e) ) then
      print *, 'warning : e2e in fem-grid is already initialized!'
      deallocate(grd%e2e)
    end if

    if (grd%nquad4 .eq. 0 ) then !only tri, old way!
      allocate(grd%e2e(grd%ncellsg, 3))
    else
      allocate(grd%e2e(grd%ncellsg, 4))
    end if

    ! create e2e map
    do i = 1, size(grd%icon,1)

select case (grd%elname(i))
case (GEN_TRIANGLE) 

       ! neigh 1
       pt1 = grd%icon(i,2); pt2 = grd%icon(i,3)
       grd%e2e(i,1) = find_elem_cont_pts(grd%icon, i, pt1, pt2)

       ! neigh 2
       pt1 = grd%icon(i,3); pt2 = grd%icon(i,1)
       grd%e2e(i,2) = find_elem_cont_pts(grd%icon, i, pt1, pt2)

       ! neigh 3
       pt1 = grd%icon(i,1); pt2 = grd%icon(i,2)
       grd%e2e(i,3) = find_elem_cont_pts(grd%icon, i, pt1, pt2)

case (GEN_QUADRI)

       ! neigh 1
       pt1 = grd%icon(i,2); pt2 = grd%icon(i,3)
       grd%e2e(i,1) = find_elem_cont_pts(grd%icon, i, pt1, pt2)

       ! neigh 2
       pt1 = grd%icon(i,3); pt2 = grd%icon(i,4)
       grd%e2e(i,2) = find_elem_cont_pts(grd%icon, i, pt1, pt2)

       ! neigh 3
       pt1 = grd%icon(i,4); pt2 = grd%icon(i,1)
       grd%e2e(i,3) = find_elem_cont_pts(grd%icon, i, pt1, pt2)

       ! neigh 4
       pt1 = grd%icon(i,1); pt2 = grd%icon(i,2)
       grd%e2e(i,4) = find_elem_cont_pts(grd%icon, i, pt1, pt2)

end select

    end do

    ! done here

  contains 

    function find_elem_cont_pts(icon, elem, pt1, pt2)
      implicit none
      integer, dimension(:,:), intent(in) :: icon
      integer, intent(in) :: elem, pt1, pt2
      integer :: find_elem_cont_pts

      ! local vars
      integer :: i, j
      logical :: pt1_found, pt2_found


      do i = 1, size(icon,1) 
         if ( i .eq. elem) cycle

         pt1_found = .false.; pt2_found = .false.
         do j = 1, size(icon,2)
            if (icon(i,j) .eq. pt1) pt1_found = .true.
            if (icon(i,j) .eq. pt2) pt2_found = .true.
         end do

         if ( pt1_found .and. pt2_found ) then
            find_elem_cont_pts = i
            return
         end if

      end do

      ! not found yet? that's wall man!
      find_elem_cont_pts = -1 ! wall :(

      ! done gere
    end function find_elem_cont_pts

  end subroutine create_e2e

  !> @brief
  !> Main subroutine for reading segment files.
  !> @details
  !> The segment files will be in the form: <br>
  !> number of points of connector 1 <br>
  !> x1 y1 z1 <br>
  !> x2 y2 z2 <br>
  !> .  .  .  <br>
  !> .  .  .  <br>
  !> number of points of connector 2  <br>
  !> x1 y1 z1 <br>
  !> x2 y2 z2 <br>
  !> .  .   .   <br>
  !> .  .   .   <br>
  !> where all points are SEQUNTIAL.
  !> It reads all point coordinates into corresponding
  !> boundary curve data struct in "grd" type.

  subroutine read_segment_file(infile, mode, grd)
    implicit none
    character(len = *), intent(in) :: infile  !<   Name of the input segment file.
    character(len = *), intent(in) :: mode    !<   Mode of the input segment file.
    type(grid), target, intent(inout) :: grd  !<   The related grid object.

    !> @todo Do we need definition of the local vars or how to exclude them from the documentation.
    integer :: i, j
    integer :: istat
    integer :: npt
    real*8 :: xx, yy, zz
    integer :: tot_curves, tot_seg
    real*8, dimension(:), pointer :: t, x, y
    integer :: nc ! n of curves
    type(curve), pointer :: tcurve 

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
    grd%tot_repeated_bn_nodes = 0

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
       grd%tot_repeated_bn_nodes = grd%tot_repeated_bn_nodes + npt

       do i = 1, npt
          read(9,*,iostat = istat) xx, yy, zz
       end do

    end do ! reading finishes after this completes

    print *, 'total boundary curves (connector) : ', tot_curves
    print *, 'total boundary segments : ', tot_seg

    grd%tot_bn_seg = tot_seg

    allocate(grd%bn_curves(tot_curves))

    ! now fill grd%bn_curves
    rewind( unit = 9) ! go back to header

    ! read all lines of the input file
    do j = 1, tot_curves
       read(9,*,iostat = istat) npt 

       allocate(grd%bn_curves(j)%x(npt), grd%bn_curves(j)%y(npt))

       print *, 'for bn_curves(',j,')'
       do i = 1, npt
          read(9,*,iostat = istat) grd%bn_curves(j)%x(i), grd%bn_curves(j)%y(i), zz
          print *, 'x = ', grd%bn_curves(j)%x(i), ' y = ', grd%bn_curves(j)%y(i)
       end do

    end do ! reading finishes after this completes

    ! closing the input file
    close(9)

    ! fit and store spline
    do j = 1, tot_curves

       tcurve => grd%bn_curves(j)
       x => tcurve%x
       y => tcurve%y
       nc = size(x)
       allocate(tcurve%t(nc))
       t => tcurve%t
       t = (/ (dble(i) , i = 1, nc) /)
       if (.not. allocated(tcurve%Mx) ) then
          call spline_nurbs_alloc(tcurve%Mx, tcurve%My, tcurve%a, tcurve%b, &
                & tcurve%c, tcurve%d, tcurve%cp, tcurve%dp, nc, tcurve%btype)
!         allocate(tcurve%Mx(nc), tcurve%My(nc), tcurve%a(nc) &
!                , tcurve%b(nc), tcurve%c(nc), tcurve%d(nc)    &
!                , tcurve%cp(nc), tcurve%dp(nc) )
       end if
       call spline_nurbs_comp(x, y, tcurve%a, tcurve%b, tcurve%c, tcurve%d, &
                           & tcurve%cp, tcurve%dp, tcurve%Mx, tcurve%My, t, mode, tcurve%btype)
!      ! parametrize the x-coords
!      call comp_spline_weights(x, t , tcurve%Mx, tcurve%a, tcurve%b, tcurve%c &
!                             , tcurve%d, tcurve%cp, tcurve%dp, mode)
!      ! parametrize the y-coords
!      call comp_spline_weights(y, t , tcurve%My, tcurve%a, tcurve%b, tcurve%c &
!                             , tcurve%d, tcurve%cp, tcurve%dp, mode)


    end do

    ! done here
  end subroutine read_segment_file

  ! scatter plot very handy for visualization of
  ! curves boundaries using additional scattered 
  ! points interpolated on boundary parametrization.  
  subroutine scatter_tecplot(outfile, append_flag, title, x, y, u)
    implicit none
    character(len = *), intent(in) :: outfile, title
    logical, intent(in) :: append_flag
    real*8, dimension(:), intent(in) :: x, y
    real*8, dimension(:,:), intent(in) :: u

    ! local vars
    integer i, j, neqs, nnodesg 
    ! init
    neqs = size(u,1)
    nnodesg = size(u,2)

    ! error checking
    if (nnodesg .ne. size(x)) then
       print *,' error in scatter_tecplot(...) : nnodesg of '&
            ,' [u] is not equal to the size of scatter x,y!. stop'
       stop
    end if

    ! opening the output file
    if (append_flag) then
       ! append
       open(unit = 10, file = outfile, position = 'append')
    else !or just open a new file
       open(unit = 10, file = outfile)
    end if

    ! write header
    write(10, *) 'title = "',title,'"'

    write(10, '(A)', advance = 'no') 'variables = "x", "y"'
    do i = 1, neqs
       write(10, '(A, I1, A)', advance = 'no') ', "u', i,'"'
    end do

    write(10,*)
    write(10, '(A, I7, A)', advance = 'no') 'zone I = ' &
         , nnodesg, ', datapacking = point'

    write(10,*)
    ! write coordinates and values of the vector field [u]
    do i = 1, nnodesg
       write(10, '(F30.17, A, F30.17)', advance = 'no') x(i), ' ',  y(i)
       do j = 1, neqs
          write(10, '(A, F30.17)', advance = 'no') ' ',  u(j,i)
       end do
       write(10,*) 
    end do

    ! shut down the output
    close(10)

    ! done here
  end subroutine scatter_tecplot

  ! visualization of curves boundaries using 
  ! additional scattered points interpolated 
  ! on boundary parametrization. 
  subroutine visual_curved_bonds(outfile, grd, factor)
    implicit none
    character(len = *), intent(in) :: outfile
    type(grid), target, intent(in) :: grd
    integer, intent(in) :: factor

    ! local vars
    integer :: i, j, nc, nscat
    type(curve), pointer :: tcurve
    real*8, dimension(:), pointer :: t, x, y
    real*8, dimension(:), allocatable :: tt, xx, yy
    real*8, dimension(:,:), allocatable :: u
    logical :: append_flag
    character(len=81) :: title, opt


    append_flag = .true.
    title = 'curve'
    opt = 'interp'

    do j = 1, size(grd%bn_curves)
       tcurve => grd%bn_curves(j)
       t => tcurve%t
       x => tcurve%x
       y => tcurve%y
       nc = size(x)
       nscat = factor * nc
       allocate(tt(nscat), xx(nscat), yy(nscat), u(1, nscat))
       u = 1.0d0
       tt = (/ (1.0d0 + dble(i - 1)/dble(nscat - 1) * dble(nc - 1) &
             , i = 1, nscat) /)
       do i = 1, nscat
          call spline_nurbs_eval2(tt(i), x, y, tcurve%a, tcurve%b, tcurve%c, &
                       & tcurve%Mx, tcurve%My, t, xx(i), yy(i), opt, tcurve%btype)
!         call spline_eval2(t, x, tt(i), xx(i), opt, tcurve%Mx)
!         call spline_eval2(t, y, tt(i), yy(i), opt, tcurve%My)
       end do
       ! write to the file
       call scatter_tecplot(outfile, append_flag, title, xx, yy, u)
       ! clean it before go to next curve
       deallocate(tt, xx, yy, u)
    end do

    ! done here
  end subroutine visual_curved_bonds

  
  ! eltypein = 0 : Lagrange elements
  !          = 1 : Chebyshev (fekete) elements
  subroutine add_more_points(grd, pin, eltypein, tolerance, galtype, n1, n2)
    implicit none
    type(grid), intent(inout) :: grd
    integer, dimension(:), intent(in) :: pin, eltypein
    real*8, intent(in) :: tolerance
    character(len = 2), intent(in), optional :: galtype ! Galerkin type = PG, DG
    integer, dimension(:), optional :: n1, n2

    ! local vars
    integer :: j, k, ntri
    integer :: pmax, npemax, npedmax
    ! integer :: nped_tmp
    integer, dimension(:,:), allocatable :: tmp
    integer :: rule
    real*8, dimension(:), allocatable :: w, xnew, ynew 
    real*8, dimension(:,:), allocatable :: xy, rtmp 
    integer :: li ! last index
    integer :: tn1, tn2

    ! init
    ntri = grd%ntri

    ! bullet proofing 
    if ( all( pin .eq. grd%p ) .and. &
         all( eltypein .eq. grd%eltype ) .and. (grd%galtype .eq. galtype) ) then 
       print *, 'warning : both "p" and element types are' &
            , ' equal to the current config. no midification is done!'
       return
    end if
    ! if(grd%nquad4 .ne. 0 ) then
    !    print *, 'current version of p-refinement is only for triangles. quads were' &
    !         , '  found in the gird! stop.'
    !    stop
    ! end if

    ! resize(reset) everything if previously 
    ! made high-order (nnodesg and ncellsg change as a flag)
    if ( (grd%nnodesg .ne. grd%nnodesg0) .or. &
         (grd%ncellsg .ne. grd%ncellsg0) ) then 
       print *, 'warning : resetting the grid to its original snapshot ...'
       call reset_grid(grd)
    end if


    ! set the order of accuracy and element type
    grd%p = pin
    grd%eltype = eltypein
    if (present(galtype) ) grd%galtype = galtype

    ! before anything, compute the maximum number of points 
    ! that would be inserted per element after this.
    ! then resize "grd" structure accordingly
    pmax = maxval(grd%p)
    if ( grd%nquad4 .eq. 0 ) then ! do it the old way
       ! compute the maximum number of points per element and edge
       call p2node(pmax, npemax, npedmax)
       ! resize grd%icon
       allocate(tmp(ntri, npemax))
       tmp = -1
       tmp(:,1:3) = grd%icon
    else ! mixed element
       npemax = (pmax + 1)**2
       ! resize grd%icon
       allocate(tmp(grd%ncellsg, npemax))
       tmp = -1

       do k = 1, grd%ncellsg
          select case (grd%elname(k))
          case (GEN_TRIANGLE)
             tmp(k,1:3) = grd%icon(k,1:3)
          case (GEN_QUADRI)
             tmp(k,1:4) = grd%icon(k,1:4)
          end select
       end do

    end if
    call move_alloc(tmp, grd%icon)

    ! snap boundary points to the boundary curves
    ! if they are not already on the curve
    call snap2curve(grd, tolerance)
    print *, 'done snapping!'

    ! now start to fill everything required
    do j = 1, grd%ncellsg   

select case ( grd%elname(j) )

case (GEN_TRIANGLE )
       ! first compute the master element coordinates 
       ! of high-order points 
       if ( grd%eltype(j) .eq. 0 ) then !lagrange triangles
          rule = grd%p(j) + 1
          call ncc_triangle_order_num( rule, grd%npe(j) )
       elseif ( grd%eltype(j) .eq. 1 ) then ! Chebyshev triangle
          call fekete_degree2rule( grd%p(j), rule)
          call fekete_order_num( rule, grd%npe(j) )
       else
          print *, 'unknown element type! stop'
          stop
       end if

       ! computing nped
       ! call npe2nped( grd%npe(j), nped_tmp )

       ! allocate space for that
       allocate( xy(2,grd%npe(j)), w(grd%npe(j)) &
            , xnew(grd%npe(j)), ynew(grd%npe(j)) )

       if ( grd%eltype(j) .eq. 0 ) then 
          ! compute the absicca and weights for Lagrange element
          call ncc_triangle_rule( rule, grd%npe(j), xy, w )

          ! odd rule needs more three corner points
          if ( mod(rule, 2) .ne. 0) then ! is odd!
             li = grd%npe(j)
             allocate( rtmp(2, (li + 3) ) )
             rtmp(:, 4:(li+3) ) = xy(:, 1:li )
 
             rtmp(1, 1 ) = 0.0d0; rtmp(2, 1 ) = 0.0d0
             rtmp(1, 2 ) = 1.0d0; rtmp(2, 2 ) = 0.0d0
             rtmp(1, 3 ) = 0.0d0; rtmp(2, 3 ) = 1.0d0

             call move_alloc(rtmp, xy)
             grd%npe(j) = grd%npe(j) + 3 ! more three points added
             deallocate(w, xnew, ynew)
             allocate( w(grd%npe(j)), xnew(grd%npe(j)), ynew(grd%npe(j)) ) 
          end if

       elseif ( grd%eltype(j) .eq. 1 ) then 

          ! compute the absicca and weights for fekete
          call fekete_rule( rule, grd%npe(j), xy, w )

       end if

       ! reorder master element points according to counterclockwise def.
       call fem_reorder(xy, tolerance)

       ! store master element (xi, eta) interpolation points
       ! for any future use
   if (allocated(grd%maselem(j)%xi)  ) deallocate(grd%maselem(j)%xi)
   if (allocated(grd%maselem(j)%eta) ) deallocate(grd%maselem(j)%eta)
       allocate(grd%maselem(j)%xi(grd%npe(j)), grd%maselem(j)%eta(grd%npe(j)) )
       grd%maselem(j)%xi  = xy(1,:)
       grd%maselem(j)%eta = xy(2,:)
 
       ! then transform to the physical element coordinates
       ! using exact transformation
       ! transform_elem(grd, elem, xi, eta, x, y, tol)
       do k = 1, grd%npe(j)
          call transform_elem(grd, j , xy(1,k), xy(2,k), xnew(k), ynew(k), tolerance)
       end do

case (GEN_QUADRI)

   ! compute directional point distribution for this quad elem 
   if( present(n1) .and. present(n2) ) then
      tn1 = n1(j); tn2 = n2(j)
      ! check see if n1 and n2 are consistent ...
      if ( (tn1 * tn2) .ne. (grd%p(j) + 1)**2 ) then
         print *, 'n1(',j, ') and n2(',j,',) are not consistent with given <p(j)>! stop'
         stop
      end if
   else
      tn1 = (grd%p(j) + 1)
      tn2 = tn1
   end if

   ! compute points per element
   grd%npe(j) = tn1 * tn2

   ! allocate space for tmp arrays
   allocate( xy(2,grd%npe(j)), w(grd%npe(j)) &
        , xnew(grd%npe(j)), ynew(grd%npe(j)) )

   ! compute the master element coordinates 
   ! of high-order points 
   if (allocated(grd%maselem(j)%xi)  ) deallocate(grd%maselem(j)%xi)
   if (allocated(grd%maselem(j)%eta) ) deallocate(grd%maselem(j)%eta)
   if ( grd%eltype(j) .eq. 0 ) then !lagrange quad
      call gen_master_quadri(tn1, tn2, -1.0d0, 1.0d0, -1.0d0, 1.0d0 &
           , 'equal_space', grd%maselem(j)%xi, grd%maselem(j)%eta)
   elseif ( grd%eltype(j) .eq. 1 ) then ! Chebyshev quad
      call gen_master_quadri(tn1, tn2, -1.0d0, 1.0d0, -1.0d0, 1.0d0 &
           , 'cheby', grd%maselem(j)%xi, grd%maselem(j)%eta)
   else
      print *, 'unknown element type for quad! stop'
      stop
   end if
   xy(1, :) = grd%maselem(j)%xi
   xy(2, :) = grd%maselem(j)%eta

   ! already reordered!!! so proceed to
   ! then transform to the physical element coordinates
   ! using exact transformation
   do k = 1, grd%npe(j)
      ! subroutine transform_elem_quadri(grd, elem, xi, eta, x, y, tol)
      call transform_elem_quadri(grd, j, xy(1,k), xy(2,k) &
           , xnew(k), ynew(k), tolerance)
   end do

case default
print *, 'unknown element name in add_more_points(...)! stop'
stop

end select

       ! now add the new points to the "grd" data structure if accepted
       call add_xnew(grd, j , xy, xnew, ynew, tolerance)


       ! clean-ups
       deallocate(xy, w, xnew, ynew)

end do

! update e2e map for both PG and DG methods 
! based on coordinate searching
call create_e2e_based_physical_coords(grd)
! refill ibedgeELEM_local_edg
call fill_ibedgeELEM_local_edg_coords(grd)

! gather local edge info in one place suitable for 
! DG grid processor and solvers
call fill_local_edg_bc(grd)

! if (grd%galtype .eq. 'PG') then
!     ! fix the e2e map if any changes is made 
!     ! in "transform_elem" in the winding of the
!     ! physical elements
!     ! deallocate(grd%e2e) 
!     ! call create_e2e(grd)

!     ! refill ibedgeELEM_local_edg
!     ! call fill_ibedgeELEM_local_edg(grd)

! end if


    ! done here
  end subroutine add_more_points

  ! removes all higher-order elements and
  ! reset all elements back to the original
  ! linear first-order elements
  subroutine reset_grid(grd)
    implicit none
    type(grid), intent(inout) :: grd

    ! local vars
    integer :: i 
    real*8, dimension(:), allocatable :: rtmp !real tmp
    integer, dimension(:,:), allocatable :: tmp

    ! ! bullet proofing
    ! if ( grd%nquad4 .ne. 0 ) then
    !    print *, 'current version of reset_grid(...) is only' &
    !           , ' for triangle. some quad elems were found! stop'
    !    stop
    ! end if

    ! resize x, y
    allocate(rtmp(grd%nnodesg0))
    ! first "x"
    rtmp = grd%x(1:grd%nnodesg0)
    deallocate(grd%x)
    allocate(grd%x(grd%nnodesg0))
    grd%x = rtmp
    ! second "y"
    rtmp = grd%y(1:grd%nnodesg0)
    deallocate(grd%y)
    allocate(grd%y(grd%nnodesg0))
    grd%y = rtmp
    ! done
    deallocate(rtmp)

    ! resize grd%icon
    if (grd%nquad4 .eq. 0) then ! old way

       allocate(tmp(grd%ncellsg0, 3))
       tmp = grd%icon(:,1:3)

    else

       allocate(tmp(grd%ncellsg0, 4))
       tmp = 0
       do i = 1, grd%ncellsg0

          select case ( grd%elname(i) )
          case (GEN_TRIANGLE)
             tmp(i, 1:3) = grd%icon(i,1:3)
          case (GEN_QUADRI) 
             tmp(i, 1:4) = grd%icon(i,1:4)
          case default
             print *, 'unknown element name in reset_grid(...)! stop'
             stop
          end select

       end do

    end if

    call move_alloc(tmp, grd%icon)


    ! reset everything else back to its original 
    ! status
    grd%nnodesg = grd%nnodesg0
    grd%ncellsg = grd%ncellsg0
    grd%nbedgeg = grd%nbedgeg0 
    grd%p = 1
    grd%eltype = 0

    ! done here
  end subroutine reset_grid

  ! converts order "p" to the number of required points;
  ! npe : number of points per element
  ! nped : number of points per edge
  subroutine p2node(p, npe, nped)
    implicit none
    integer, intent(in)  :: p
    integer, intent(out) :: npe, nped

    if( p < 1) then
       print *,' minimum order should be p = 1. stop.'
       stop
    end if

    nped = p+1
    npe = (p+1) * (p+2) / 2

    ! done here
  end subroutine p2node

  ! find points per edge of the element
  ! given the number of points inside
  subroutine npe2nped(npe, nped)
    implicit none
    integer, intent(in) :: npe
    integer, intent(out) :: nped

    ! local vars
    real*8 :: delta

    delta = sqrt(1.0d0 + 8.0d0 * dble(npe))
    nped = maxval( (/ nint(-1.5d0 + 0.5d0 * delta) , nint(-1.5d0 - 0.5d0 * delta) /) )  
    nped = nped + 1

    if ( nped <= 1 ) then
       print *, 'error: points per edge are <=1! should be >= 2. stop'
       stop
    end if

    ! done here
  end subroutine npe2nped

  ! finsds the parameter t on the boundaries
  ! based on physical coords (x, y)
  !
  !
  ! if found then 0.0d0 <= t <= 1.0d0
  ! if not found t = -1.0d0  
  !
  subroutine find_t(grd, tag, x, y, tol, t)
    implicit none
    type(grid), target, intent(in) :: grd
    integer, intent(in) :: tag
    real*8, intent(in) :: x, y, tol
    real*8, intent(out) :: t

    ! local vars
    integer :: i
    real*8 :: t_star, abs_dx, abs_dy
    type(curve), pointer :: tcurve
    real*8, dimension(:), pointer :: xs, ys

    ! init
    tcurve => grd%bn_curves(tag)
    xs => tcurve%x
    ys => tcurve%y
    t = 1.0d0

! print *, 'tag = ', tag
! print *, 'xs = ', xs
! print *, 'ys = ', ys

    do i = 1, (size(xs) - 1) ! search all static points on that curve

       ! bullet proofing ...
       if ( sqrt( (xs(i+1) - xs(i)) * (xs(i+1) - xs(i)) &
                + (ys(i+1) - ys(i)) * (ys(i+1) - ys(i)) ) <= (sqrt(2.0d0) * tol) ) then
          print *, 'fatal : the distance between two static curve points' &
               , ' i = ', i, ', i+1 = ', (i+1), ' on curve = ', tag, ' is smaller than' & 
               , ' the tolerance! we better to quit to prevent wierd bugs! stop.'
          stop
       end if

       ! determine t_star in a robust way !   
       abs_dx = abs(xs(i+1) - xs(i))
       abs_dy = abs(ys(i+1) - ys(i))
       if ( (abs_dy >= abs_dx) .and.  (abs_dy >= tol) ) then
          t_star = (y - ys(i)) / (ys(i+1) - ys(i))
       elseif ( (abs_dx > abs_dy) .and.  (abs_dx >= tol) ) then
          t_star = (x - xs(i)) / (xs(i+1) - xs(i))
       else
          print *, 'fatal : inconsistency in computing t_star!' &
                 , ' none of dx or dy of boundary edge is greater' &
                 , ' than the tolerance!!! stop.'
          stop 
       end if

       ! print *, 'i = ', i
       ! print *, 't_star = ', t_star

       if( (abs( (y - ys(i)) * (xs(i+1) - xs(i)) - (x - xs(i)) * (ys(i+1) - ys(i)) ) <= tol) &
            .and. (t_star >= (0.0d0 - tol) ) .and. (t_star <= (1.0d0 + tol) ) ) then
          t = t + t_star
          ! print *, 'final t = ', t 
          return
       end if
       t = t + 1.0d0
    end do
    ! print *, 'x = ', x
    ! print *, 'y = ', y
    ! print *, 't = ', t
    ! print *, 'tol = ', tol
    ! print *, 't_star = ', t_star

    t = -1.0d0 ! not found

    ! done here
  end subroutine find_t

  ! snaps the BOUNDARY grid points
  ! to the boundary curves. 
  subroutine snap2curve(grd, tol)
    implicit none
    type(grid), target, intent(inout) :: grd
    real*8, intent(in) :: tol

    ! local vars
    integer :: elem, tag, edgnum, pt1, pt2
    real*8 :: tt, xx, yy
    type(curve), pointer :: tcurve
    real*8, dimension(:), pointer :: t, xs, ys
    logical, dimension(:), allocatable :: snapped

    ! init
    allocate(snapped(grd%nnodesg))
    snapped = .false. ! nothing snapped yet!

    do elem = 1, grd%ncellsg

       ! first determine points 1, 2 
       ! on the boundary edge of this element
       ! if this is a boundary element
       tag = grd%el2bn(elem, 1)

       if (tag .eq. 0) cycle ! interior element no need to snap
       ! then it is a boundary element

       ! get that curve
       tcurve => grd%bn_curves(tag)
       t => tcurve%t
       xs => tcurve%x
       ys => tcurve%y

       ! compute that boundary edge and points
       edgnum  = grd%el2bn(elem, 2)
       pt1 = grd%ibedge(edgnum,1)
       pt2 = grd%ibedge(edgnum,2)

       ! snapping point1
       if ( .not. snapped(pt1) ) then 
          call find_t(grd, tag, grd%x(pt1), grd%y(pt1), tol, tt)
          if ( tt .eq. -1.0d0) then
             print *, 'error : the parameter t can not be computed' &
                  , ' for point1 # ',pt1,' of edg # ', edgnum,'. stop.'
             stop
          end if

          ! computing snapped coordinate
          call spline_nurbs_eval2(tt, xs, ys, tcurve%a, tcurve%b, tcurve%c &
                       , tcurve%Mx, tcurve%My, t, xx, yy, 'interp', tcurve%btype)
!         call spline_eval2(t, xs, tt, xx, 'interp', tcurve%Mx)
!         ! updating x
          grd%x(pt1) = xx
!         ! computing snapped coordinate
!         call spline_eval2(t, ys, tt, yy, 'interp', tcurve%My)
!         ! updating y
          grd%y(pt1) = yy

          ! set the flag
          snapped(pt1) = .true.

       end if

       ! snapping point2
       if ( .not. snapped(pt2) ) then 
          call find_t(grd, tag, grd%x(pt2), grd%y(pt2), tol, tt)
          if ( tt .eq. -1.0d0) then
             print *, 'error : the parameter t can not be computed' &
                  , ' for point2 # ',pt2,' of edg # ', edgnum,'. stop.'
             stop
          end if

          ! computing snapped coordinate
          call spline_nurbs_eval2(tt, xs, ys, tcurve%a, tcurve%b, tcurve%c &
                       , tcurve%Mx, tcurve%My, t, xx, yy, 'interp', tcurve%btype)
!         call spline_eval2(t, xs, tt, xx, 'interp', tcurve%Mx)
!         ! updating x
          grd%x(pt2) = xx
!         ! computing snapped coordinate
!         call spline_eval2(t, ys, tt, yy, 'interp', tcurve%My)
!         ! updating y
          grd%y(pt2) = yy

          ! set the flag
          snapped(pt2) = .true.

       end if

    end do ! for all elements

    ! clean ups
    deallocate(snapped)

    ! done here
  end subroutine snap2curve

  ! general master element r : 0->1, s : 0->1
  ! to physical coordinate transformation
  ! for linear and curved (nonlinear) elements
  
  subroutine transform_elem(grd, elem, xi, eta, x, y, tol)
    implicit none
    type(grid), target, intent(inout) :: grd
    integer, intent(in) :: elem
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
    character(len=*), parameter :: opt = 'interp'
    integer :: i, edgnum
    real*8 :: tangx, tangy, norm_tang, tmpx, tmpy

    ! first determine points 1, 2, 3
    tag = grd%el2bn(elem, 1)
    if (tag .eq. 0) then ! interior
       pt1 = grd%icon(elem, 1)
       pt2 = grd%icon(elem, 2)
       pt3 = grd%icon(elem, 3)
    else ! boundary
       edgnum  = grd%el2bn(elem, 2)
       pt1 = grd%ibedge(edgnum,1)
       pt2 = grd%ibedge(edgnum,2)
       do i = 1, 3
          pt3 = grd%icon(elem, i)
          if ( (pt3 .eq. pt1) .or. (pt3 .eq. pt2) ) then
             cycle
          else
             exit
          end if
       end do

       ! update local edge number in that element
       ! NOW IT IS ALWAYS 1
       grd%ibedgeELEM_local_edg(edgnum) = 1
       ! changing physical elem winding to match with master elem
       grd%icon(elem, 1) = pt1
       grd%icon(elem, 2) = pt2
       grd%icon(elem, 3) = pt3
    end if

    ! find coordinates
    x1_x = grd%x(pt1)
    x1_y = grd%y(pt1)
    x2_x = grd%x(pt2)
    x2_y = grd%y(pt2)
    x3_x = grd%x(pt3)
    x3_y = grd%y(pt3)

    ! find tangential vector to side y1
norm_tang = sqrt( (grd%x(pt2) - grd%x(pt1))**2.0d0 + (grd%y(pt2) - grd%y(pt1))**2.0d0 )
tangx = (grd%x(pt2) - grd%x(pt1)) / norm_tang
tangy = (grd%y(pt2) - grd%y(pt1)) / norm_tang

    ! compute Cx, Cy
    if( tag .eq. 0 ) then ! no curve bn
       Cx = x1_x + xi * (x2_x - x1_x)
       Cy = x1_y + xi * (x2_y - x1_y)
    else ! do boundary curve interpolation

       tcurve => grd%bn_curves(tag)
       xs => tcurve%x
       ys => tcurve%y
       t  => tcurve%t
       nc = size(xs)

       ! find_t(grd, tag, x, y, tol, t)
tmpx = grd%x(pt1) + piecewise_tol * tangx
tmpy = grd%y(pt1) + piecewise_tol * tangy 
       call find_t(grd, tag, tmpx, tmpy, tol, t1)
tmpx = grd%x(pt2) - piecewise_tol * tangx
tmpy = grd%y(pt2) - piecewise_tol * tangy 
       call find_t(grd, tag, tmpx, tmpy, tol, t2)
t1 = dble(nint(t1))
t2 = dble(nint(t2))

       tt = t1 + xi * (t2 - t1)

       ! print *, 'finding Cx for curve tag = ', tag
       ! print *, 'x = ', xs
       ! print *, 'y = ', ys
       ! print *, 't = ', t
       ! print *, 't1 = ', t1, 't2 = ', t2, 'tt = ', tt

       call spline_nurbs_eval2(tt, xs, ys, tcurve%a, tcurve%b, tcurve%c &
                   , tcurve%Mx, tcurve%My, t, Cx, Cy, 'interp', tcurve%btype)
!      call spline_eval2(t, xs, tt, Cx, opt, tcurve%Mx)
!      call spline_eval2(t, ys, tt, Cy, opt, tcurve%My)
       ! print *, 'done computing C for curve ', tag

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
  end subroutine transform_elem

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
  ! line definition. Refer to subroutine transform_elem(...)
  ! to see how it is done for a triangle.   
  
  subroutine transform_elem_quadri(grd, elem, xi, eta, x, y, tol)
    implicit none
    type(grid), target, intent(inout) :: grd
    integer, intent(in) :: elem
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
    character(len=*), parameter :: opt = 'interp'
    integer :: i, edgnum, j, k
    logical :: duplicate
    real*8 :: y1x, y1y, y3x, y3y
    real*8 :: area
    integer :: twist, pt_tmp
    real*8 :: tangx, tangy, norm_tang, tmpx, tmpy

    ! first determine points 1, 2, 3, 4
    tag = grd%el2bn(elem, 1)
    if (tag .eq. 0) then ! interior
       pt = grd%icon(elem, 1:4)

    else ! boundary

       edgnum  = grd%el2bn(elem, 2)
       pt(1) = grd%ibedge(edgnum,1)
       pt(2) = grd%ibedge(edgnum,2)

       ! now find pt(3) and pt(4)
       do k = 3, 4
          do j = 1, 4
             pt(k) = grd%icon(elem, j)
             duplicate = .false. 
             do i = 1, (k-1)          
                if ( pt(k) .eq. pt(i) ) then
                   duplicate = .true.
                   exit             
                end if
             end do
             if ( duplicate ) then
                cycle
             else
                exit
             end if
          end do
       end do

    end if

    ! see if the chosen points lead to
    ! a good or twisted element, if
    ! twisted then swap points 3, 4
    ! to make it good element
    !
    ! subroutine comp_quad4_area(x, y, area, twist)
    !
    call comp_quad4_area(grd%x(pt), grd%y(pt), area, twist)
    if ( twist .eq. 1 ) then ! twisted!

       !swapping ...
       pt_tmp = pt(4)
       pt(4) = pt(3)
       pt(3) = pt_tmp
       ! double check
       call comp_quad4_area(grd%x(pt), grd%y(pt), area, twist)
       if (twist .eq. 1) then ! still twisted wow !!!
          print *, 'the quad element #', elem, 'is twisted and cant be fixed! stop'
          stop
       end if

    end if

    ! updating physical elem winding to match with master elem
    if( any(pt .ne. grd%icon(elem, 1:4)) ) then
       grd%icon(elem, 1:4) = pt

       ! update local edge number in that element
       ! NOW IT IS ALWAYS 1
       grd%ibedgeELEM_local_edg(edgnum) = 1

    end if

    ! find coordinates
    xx = grd%x(pt)
    yy = grd%y(pt)

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
  end subroutine transform_elem_quadri

  ! computes the area of a general quad4 and indicates
  ! if it is twisted or not!
  subroutine comp_quad4_area(x, y, area, twist)
    implicit none
    real*8, dimension(:), intent(in) :: x, y
    real*8, intent(out) :: area
    integer, intent(out) :: twist

    ! local vars
    real*8 :: area1, area2

    area1 = tri_area(x( (/1, 2, 4/) ), y( (/1, 2, 4/) ) )
    area2 = tri_area(x( (/2, 3, 4/) ), y( (/2, 3, 4/) ) )

    area = area1 + area2

    if ( (area1 > 0.0d0) .and. (area2 > 0.0d0) ) then
       twist = 0
    else
       twist = 1
    end if

    ! done here
  contains

    function tri_area(xx, yy)
      implicit none
      real*8, dimension(3), intent(in) :: xx, yy
      real*8 :: tri_area

      ! local vars
      real*8 :: ax, ay, bx, by

      ax = xx(2) - xx(1)
      ay = yy(2) - yy(1)

      bx = xx(3) - xx(1)
      by = yy(3) - yy(1)

      tri_area = 0.5d0 * (ax * by - bx * ay)

      ! done here
    end function tri_area

  end subroutine comp_quad4_area

  ! add new points "xnew", "ynew" in physical coordinates
  ! to the current element "elem"
  subroutine add_xnew(grd, elem , xy, xnew, ynew, tolerance)
    implicit none
    type(grid), intent(inout) :: grd
    integer, intent(in) :: elem
    real*8, dimension(:,:), intent(in) :: xy
    real*8, dimension(:), intent(in) :: xnew, ynew
    real*8, intent(in) :: tolerance

    ! local vars
    integer :: i, j, neigh, npts, indx
    real*8 :: x, y
    integer :: exists
    real*8, dimension(:), allocatable :: rtmp
    integer :: max_neigh

    ! init
    select case ( grd%elname(elem) )

    case ( GEN_TRIANGLE)
       indx = 4
       max_neigh = 3
 
    case ( GEN_QUADRI)
       indx = 5
       max_neigh = 4

    end select

    if ( grd%galtype .eq. 'DG') then
       indx = 1
       ! grd%icon(elem, :) = -1
    end if

    do j = 1, size(xnew)
       x = xnew(j); y = ynew(j)

       ! only for PG, check if the point exists in the elem and/or neighbors
       ! if exist, don't update x, y and simply just add it to connectivities
       if ( grd%galtype .eq. 'PG') then

          ! check to see if the point already exist 
          ! in the element itself 
          exists = does_point_exist(grd, elem, x, y, tolerance)
          if ( exists .ne. -1 ) go to 200 ! oh! exists and is one of skleton vertices

          ! neighbors ... 
          do i = 1, max_neigh
             neigh = grd%e2e(elem, i)
             if ( neigh .eq. -1) cycle ! wall! 
             exists = does_point_exist(grd, neigh, x, y, tolerance)
             if ( exists .ne. -1 ) go to 100
          end do

       end if
       !
       ! point (x,y) not in element and neighbors. well :) add it then!
       !
       ! 1- add to grd%x grd%y
       npts = grd%nnodesg + 1
       ! add x
       allocate(rtmp(npts))
       rtmp(1:(npts-1)) = grd%x(1:(npts-1))
       rtmp(npts) = x
       call move_alloc(rtmp, grd%x)
       ! add y
       allocate(rtmp(npts))
       rtmp(1:(npts-1)) = grd%y(1:(npts-1))
       rtmp(npts) = y
       call move_alloc(rtmp, grd%y)
       exists = npts
       ! update nodes number
       grd%nnodesg = npts

100    continue ! add to icon
       
       grd%icon(elem, indx) = exists
       indx = indx + 1

200    continue ! add to edges

       select case ( grd%elname(elem) )

       case (GEN_TRIANGLE)       

          if  (  (abs(xy(1, j) - 0.0d0) <= tolerance)  .and. & !xi  = 0
               (abs(xy(2, j) - 0.0d0) <= tolerance)  ) then  !eta = 0

          elseif (  (abs(xy(1, j) - 1.0d0) <= tolerance)  .and. & !xi  = 1
               (abs(xy(2, j) - 0.0d0) <= tolerance)  ) then  !eta = 0

          elseif (  (abs(xy(1, j) - 0.0d0) <= tolerance)  .and. & !xi  = 0
               (abs(xy(2, j) - 1.0d0) <= tolerance)  ) then  !eta = 1

          elseif (  (abs(xy(1, j) - 0.0d0) <= tolerance)  ) then !xi  = 0

             call add_pt2edg(edg = grd%el2edg(elem)%edg3, pt = exists)

          elseif (  (abs(xy(2, j) - 0.0d0) <= tolerance)  ) then !eta  = 0

             call add_pt2edg(edg = grd%el2edg(elem)%edg1, pt = exists)

          elseif (  (abs(xy(1, j) + xy(2, j) - 1.0d0) <= tolerance)  ) then 
             !xi + eta  = 1
             call add_pt2edg(edg = grd%el2edg(elem)%edg2, pt = exists)

          else

             ! do nothing in this version

          end if

       case (GEN_QUADRI)

          if  (  (abs(xy(1, j) + 1.0d0) <= tolerance)  .and. & !xi  = -1
               (abs(xy(2, j) + 1.0d0) <= tolerance)  ) then  !eta = -1

          elseif (  (abs(xy(1, j) - 1.0d0) <= tolerance)  .and. & !xi  = 1
               (abs(xy(2, j) + 1.0d0) <= tolerance)  ) then  !eta = -1

          elseif (  (abs(xy(1, j) - 1.0d0) <= tolerance)  .and. & !xi  = 1
               (abs(xy(2, j) - 1.0d0) <= tolerance)  ) then  !eta = 1

          elseif (  (abs(xy(1, j) + 1.0d0) <= tolerance)  .and. & !xi  = -1
               (abs(xy(2, j) - 1.0d0) <= tolerance)  ) then  !eta = 1

          elseif (  abs(xy(1, j) + 1.0d0) <= tolerance  ) then !xi  = -1

             call add_pt2edg(edg = grd%el2edg(elem)%edg4, pt = exists)

          elseif (  abs(xy(1, j) - 1.0d0) <= tolerance  ) then !xi  = 1

             call add_pt2edg(edg = grd%el2edg(elem)%edg2, pt = exists)

          elseif (  abs(xy(2, j) + 1.0d0) <= tolerance  ) then !eta  = -1

             call add_pt2edg(edg = grd%el2edg(elem)%edg1, pt = exists)

          elseif (  abs(xy(2, j) - 1.0d0) <= tolerance  ) then !eta  = 1

             call add_pt2edg(edg = grd%el2edg(elem)%edg3, pt = exists)

          else

             ! do nothing in this version

          end if

       end select

    end do ! next point in xnew

    ! done here
  end subroutine add_xnew

  ! searches just one element "elem" for point (x, y)
  ! within the tolerance.
  ! if already exists, then return the global number 
  ! of that point (node) which is unique >= 1
  ! if (x,y) not found, then returns -1.
  function does_point_exist(grd, elem, x, y, tolerance) 
    implicit none
    type(grid), intent(in) :: grd
    integer, intent(in) :: elem
    real*8, intent(in) :: x, y, tolerance
    integer :: does_point_exist

    ! local vars
    integer :: j, pt

    does_point_exist = -1 ! like NULL 

    do j = 1, grd%npe(elem)

       pt = grd%icon(elem, j) 
       if (pt .eq. -1) exit !end of points
       if ( (abs( x - grd%x(pt)) <= tolerance) .and. &
            (abs( y - grd%y(pt)) <= tolerance) ) then ! found it!
          does_point_exist = pt
          return
       end if

    end do ! not found after this completes

    ! done here
  end function does_point_exist
  
  ! adds global point (node) number "pt >= 1" 
  ! to the end of array edg
  subroutine add_pt2edg(edg, pt)
    implicit none
    integer, dimension(:), allocatable, intent(inout) :: edg
    integer, intent(in) :: pt

    ! local vars
    integer :: leng
    integer, dimension(:), allocatable :: tmp

    ! bulletproofing
    if (.not. allocated(edg) ) then
       allocate(edg(1))
       edg(1) = pt 
       return
    end if
    ! else the size of "edg" is at least 1
    leng = size(edg) + 1
    allocate(tmp(leng))
    tmp(1:(leng - 1)) = edg(1:(leng - 1))
    tmp(leng) = pt
    call move_alloc(tmp, edg)

    ! done here
  end subroutine add_pt2edg

  ! 
  subroutine fill_ibedgeELEM_local_edg(grd)
    implicit none
    type(grid), intent(inout) :: grd

    ! local vars
    integer :: i, pt1, pt2, ielem
    integer :: ps1, ps2, ps3, ps4

    ! bulletproofing
    if ( .not. allocated(grd%ibedgeELEM_local_edg) ) then
       print *, 'fatal: grd%ibedgeELEM_local_edg must be first allocated! stop'
       stop
    end if

    ! search all the edges to find local edg number
    do i = 1, grd%nbedgeg

       pt1 = grd%ibedge(i, 1)
       pt2 = grd%ibedge(i, 2)
       ielem = grd%ibedgeELEM(i)
       ps1 = grd%icon(ielem, 1)
       ps2 = grd%icon(ielem, 2)
       ps3 = grd%icon(ielem, 3)
       if ( grd%elname(ielem) .eq. GEN_QUADRI) ps4 = grd%icon(ielem, 4) 

       select case ( grd%elname(ielem) )

       case ( GEN_TRIANGLE )
          if     (  (pt1 .eq. ps1) .and. &
               (pt2 .eq. ps2) ) then
             grd%ibedgeELEM_local_edg(i) = 1

          elseif (  (pt1 .eq. ps2) .and. &
               (pt2 .eq. ps3) ) then
             grd%ibedgeELEM_local_edg(i) = 2

          elseif (  (pt1 .eq. ps3) .and. &
               (pt2 .eq. ps1) ) then
             grd%ibedgeELEM_local_edg(i) = 3

          else
             print *, 'fatal : the edge with nodes ibedge(i,:) = [' &
                  , grd%ibedge(i, :), '] does not match any edges' &
                  , ' in tri element #', ielem,'. stop'
             stop

          end if
       case ( GEN_QUADRI )
          if     (  (pt1 .eq. ps1) .and. &
               (pt2 .eq. ps2) ) then
             grd%ibedgeELEM_local_edg(i) = 1

          elseif (  (pt1 .eq. ps2) .and. &
               (pt2 .eq. ps3) ) then
             grd%ibedgeELEM_local_edg(i) = 2

          elseif (  (pt1 .eq. ps3) .and. &
               (pt2 .eq. ps4) ) then
             grd%ibedgeELEM_local_edg(i) = 3

          elseif (  (pt1 .eq. ps4) .and. &
               (pt2 .eq. ps1) ) then
             grd%ibedgeELEM_local_edg(i) = 4


          else
             print *, 'fatal : the edge with nodes ibedge(i,:) = [' &
                  , grd%ibedge(i, :), '] does not match any edges' &
                  , ' in quad element #', ielem,'. stop'
             stop

          end if
       end select

    end do

    ! done here
  end subroutine fill_ibedgeELEM_local_edg

  ! finds the duplicated nodes with the same (x,y) tuples
  ! with the tolerance "tol" and put them sequentially in 
  ! integer array "dup_nodes" like:
  ! 
  ! for
  ! 
  ! x = (/ 1.0d0,  2.0d0,  2.0d0, 1.0d0, -1.0d0 /)
  ! y = (/ 3.0d0, -6.0d0, -6.0d0, 3.0d0,  7.0d0 /)
  !
  ! we get :
  !
  ! dup_nodes = (/ 1, 4, 2, 3 /)
  ! so node 1 is duplicated by node 4 and node 2 
  ! is also duplicated by node 3 and so on.
  ! 
  ! NOTE : this also checks if the node is duplicated 
  !        more than one time it stops the program and 
  !        prevents nasty unexpected bugs.
  !
  
  subroutine find_dup_nodes(x, y, tol, dup_nodes)
    implicit none
    real*8, dimension(:), intent(in) :: x, y
    real*8, intent(in) :: tol
    integer, dimension(:), allocatable, intent(out) :: dup_nodes

    ! local vars
    integer :: n, i, j, n_found, n_prev
    logical :: found
    integer, dimension(:), allocatable :: tmp

    ! bulletproofing
    if ( allocated(dup_nodes) ) then
       print *, 'warning : the duplicate node list is already allocated!' &
            , ' deallocating it now ...'
       deallocate(dup_nodes)
    end if

    if( size(x) .ne. size(y) ) then
       print *, 'fatal : the size of input vectors <x> and <y> does' &
            ,' not match in find_dup_nodes! stop.'
       stop
    end if

    n = size(x)
    n_found = 0 ! not found any duplicated node yet!

    do i = 1, (n-1) ! brute force for all nodes

       found = .false. ! node <i> is not yet found!
       do j = i+1, n

          if ( (abs(x(j) - x(i)) <= tol) .and. & 
               (abs(y(j) - y(i)) <= tol) ) then ! duplicate just found!

             if ( found ) then ! check if it is more than one time

                print *, 'fatal : the node #', i ,' is duplicated' &
                       , ' more than one time! stop'
                stop
 
             else

                ! add this duplicate node to the list
                n_found = n_found + 1
                allocate( tmp( 2 * n_found) )
                if (n_found > 1) then
                   n_prev = n_found - 1
                   tmp(1:(2*n_prev)) = dup_nodes(1:(2*n_prev))
                end if
                tmp( 2*n_found - 1 ) = i
                tmp(   2*n_found   ) = j
                call move_alloc(tmp, dup_nodes)
                found = .true. 

             end if

          end if

       end do

    end do


    ! done here
  end subroutine find_dup_nodes

  ! manually creates e2e map based on physical 
  ! coordinates of the points on the edges
  ! 
  ! must be used when grid is read from CAD software
  ! specially in DG grids where points are NOT repeated
  ! in the neighboring elements
  !
  subroutine create_e2e_based_physical_coords(grd)
    implicit none
    type(grid), intent(inout) :: grd

    ! local vars
    integer :: i, pt1, pt2

    ! allocate
    if ( allocated(grd%e2e) ) then
       print *, 'warning : e2e in fem-grid is already initialized!'
       deallocate(grd%e2e)
    end if

    if (grd%nquad4 .eq. 0 ) then !only tri, old way!
       allocate(grd%e2e(grd%ncellsg, 3))
    else
       allocate(grd%e2e(grd%ncellsg, 4))
    end if

    ! create e2e map
    do i = 1, size(grd%icon,1) ! loop over all cells

select case (grd%elname(i))
case (GEN_TRIANGLE) 

       ! neigh 1
       pt1 = grd%icon(i,2); pt2 = grd%icon(i,3)
       grd%e2e(i,1) = find_elem_cont_pts_coords(grd, i, pt1, pt2)

       ! neigh 2
       pt1 = grd%icon(i,3); pt2 = grd%icon(i,1)
       grd%e2e(i,2) = find_elem_cont_pts_coords(grd, i, pt1, pt2)

       ! neigh 3
       pt1 = grd%icon(i,1); pt2 = grd%icon(i,2)
       grd%e2e(i,3) = find_elem_cont_pts_coords(grd, i, pt1, pt2)

case (GEN_QUADRI)

       ! neigh 1
       pt1 = grd%icon(i,2); pt2 = grd%icon(i,3)
       grd%e2e(i,1) = find_elem_cont_pts_coords(grd, i, pt1, pt2)

       ! neigh 2
       pt1 = grd%icon(i,3); pt2 = grd%icon(i,4)
       grd%e2e(i,2) = find_elem_cont_pts_coords(grd, i, pt1, pt2)

       ! neigh 3
       pt1 = grd%icon(i,4); pt2 = grd%icon(i,1)
       grd%e2e(i,3) = find_elem_cont_pts_coords(grd, i, pt1, pt2)

       ! neigh 4
       pt1 = grd%icon(i,1); pt2 = grd%icon(i,2)
       grd%e2e(i,4) = find_elem_cont_pts_coords(grd, i, pt1, pt2)

end select

    end do

    ! done here

  contains 

    function find_elem_cont_pts_coords(grd, elem, pt1, pt2)
      implicit none
      type(grid), intent(in), target :: grd
      integer, intent(in) :: elem, pt1, pt2
      integer :: find_elem_cont_pts_coords

      ! local vars
      integer :: i, j, ptj
      logical :: pt1_found, pt2_found
      integer, dimension(:,:), pointer :: icon => null()
      real*8 :: x1, y1, x2, y2, xj, yj

      ! init
      icon => grd%icon
x1 = grd%x(pt1)
y1 = grd%y(pt1)
x2 = grd%x(pt2)
y2 = grd%y(pt2)
 
      do i = 1, size(icon,1) 
         if ( i .eq. elem) cycle

         pt1_found = .false.; pt2_found = .false.
         do j = 1, size(icon,2)
            ptj = icon(i,j)
            if ( ptj <= 0 ) exit ! tail of connectivity!
xj = grd%x(ptj)
yj = grd%y(ptj)
 
            if ( sqrt( (x1 - xj)**2 + (y1 - yj)**2 ) <= 1.0d-14 ) pt1_found = .true.
            if ( sqrt( (x2 - xj)**2 + (y2 - yj)**2 ) <= 1.0d-14 ) pt2_found = .true.
         end do

         if ( pt1_found .and. pt2_found ) then
            find_elem_cont_pts_coords = i
            return
         end if

      end do

      ! not found yet? that's wall man!
      find_elem_cont_pts_coords = -1 ! wall :(

      ! done gere
    end function find_elem_cont_pts_coords

  end subroutine create_e2e_based_physical_coords

  ! 
  subroutine fill_ibedgeELEM_local_edg_coords(grd)
    implicit none
    type(grid), intent(inout) :: grd

    ! local vars
    integer :: i, pt1, pt2, ielem
    integer :: ps1, ps2, ps3, ps4
    real*8 :: x1(2), x2(2), xs1(2), xs2(2), xs3(2), xs4(2)

    ! bulletproofing
    if ( .not. allocated(grd%ibedgeELEM_local_edg) ) then
       print *, 'fatal: grd%ibedgeELEM_local_edg must be first allocated! stop'
       stop
    end if

    ! search all the edges to find local edg number
    do i = 1, grd%nbedgeg

       pt1 = grd%ibedge(i, 1)
       x1 = (/ grd%x(pt1), grd%y(pt1) /)
 
       pt2 = grd%ibedge(i, 2)
       x2 = (/ grd%x(pt2), grd%y(pt2) /)

       ielem = grd%ibedgeELEM(i)

       ps1 = grd%icon(ielem, 1)
       xs1 = (/ grd%x(ps1), grd%y(ps1) /)

       ps2 = grd%icon(ielem, 2)
       xs2 = (/ grd%x(ps2), grd%y(ps2) /)

       ps3 = grd%icon(ielem, 3)
       xs3 = (/ grd%x(ps3), grd%y(ps3) /)

       if ( grd%elname(ielem) .eq. GEN_QUADRI) then
          ps4 = grd%icon(ielem, 4)
          xs4 = (/ grd%x(ps4), grd%y(ps4) /)
       end if

       select case ( grd%elname(ielem) )

       case ( GEN_TRIANGLE )
          if     (  is_close(x1, xs1) .and. &
               is_close(x2, xs2) ) then
             grd%ibedgeELEM_local_edg(i) = 1

          elseif (  is_close(x1, xs2) .and. &
               is_close(x2, xs3) ) then
             grd%ibedgeELEM_local_edg(i) = 2

          elseif (  is_close(x1, xs3) .and. &
               is_close(x2, xs1) ) then
             grd%ibedgeELEM_local_edg(i) = 3

          else
             print *, 'fatal : the edge with nodes x1= [' &
                  , x1, '] and x2 = [', x2 , '] does not match any edges' &
                  , ' in tri element #', ielem,'. stop'
             stop

          end if
       case ( GEN_QUADRI )
          if     (  is_close(x1, xs1) .and. &
               is_close(x2, xs2) ) then
             grd%ibedgeELEM_local_edg(i) = 1

          elseif (  is_close(x1, xs2) .and. &
               is_close(x2, xs3) ) then
             grd%ibedgeELEM_local_edg(i) = 2

          elseif (  is_close(x1, xs3) .and. &
               is_close(x2, xs4) ) then
             grd%ibedgeELEM_local_edg(i) = 3

          elseif (  is_close(x1, xs4) .and. &
               is_close(x2, xs1) ) then
             grd%ibedgeELEM_local_edg(i) = 4

          else
             print *, 'fatal : the edge with nodes x1= [' &
                  , x1, '] and x2 = [', x2 , '] does not match any edges' &
                  , ' in quad element #', ielem,'. stop'
             stop

          end if

       end select

    end do

    ! done here
  contains

    function is_close(x1, x2)
      implicit none
      real*8, dimension(:), intent(in) :: x1, x2
      logical :: is_close

      is_close = ( sqrt( sum( (x1 - x2)**2 ) ) <= 1.0d-14 )

      ! done here
    end function is_close

  end subroutine fill_ibedgeELEM_local_edg_coords

  subroutine fill_local_edg_bc(grd)
    implicit none
    type(grid), intent(inout) :: grd

    ! local vars
    integer :: i, ielem, loc_num, bctag

    ! bullet proofing 
    if ( allocated(grd%local_edg_bc) ) deallocate(grd%local_edg_bc)
    allocate(grd%local_edg_bc(grd%ncellsg, 4))
    grd%local_edg_bc = 0 ! initially all edges are interior 

    do i = 1, size(grd%ibedge, 1) ! all boundary edges available

       ! ibedge(:,:),ibedgeBC(:), ibedgeELEM(:), ibedgeELEM_local_edg(:)
       ielem = grd%ibedgeELEM(i)
       loc_num = grd%ibedgeELEM_local_edg(i)
       bctag = grd%ibedgeBC(i)

       ! store info
       grd%local_edg_bc(ielem, loc_num) = bctag  

    end do

    ! done here
  end subroutine fill_local_edg_bc

end module grid_opt
