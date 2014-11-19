module trimesher
  use iso_c_binding, only : c_int, c_char, c_double, c_null_char
  use grid_opt
  use ps
  implicit none

  private

  interface
     subroutine trimesh(cmd_arg_string, numberofpoints, pointlist, &
          numberofsegments, segmentlist, &
          segmentmarkerlist, numberofholes, &
          holelist, report_before, report_after, &
          default_numberofpoints, default_numberoftriangles, default_numberofsegments, &
          new_numberofpoints, new_pointlist, &
          numberoftriangles, numberofcorners, trianglelist, neighborlist, &
          new_numberofsegments, new_segmentlist, &
          new_segmentmarkerlist) bind( C, name = "custom_Tri_wrapper")
       use iso_c_binding, only : c_int, c_char, c_double
       implicit none
       ! inputs
       character (c_char) :: cmd_arg_string(*)
       integer (c_int), value :: numberofpoints
       real (c_double) :: pointlist(2*numberofpoints)
       integer (c_int), value :: numberofsegments 
       integer (c_int) :: segmentlist(2*numberofsegments), segmentmarkerlist(numberofsegments)
       integer (c_int), value :: numberofholes
       real (c_double) :: holelist(2*numberofholes)
       integer (c_int), value :: report_before, report_after
       integer(c_int), value :: default_numberofpoints
       integer(c_int), value :: default_numberoftriangles
       integer(c_int), value :: default_numberofsegments

       ! outputs
       integer(c_int) :: new_numberofpoints
       real(c_double) :: new_pointlist(2 * default_numberofpoints)
       integer(c_int) :: numberoftriangles
       integer(c_int) :: numberofcorners
       integer(c_int) :: trianglelist(default_numberoftriangles*3)
       integer(c_int) :: neighborlist(default_numberoftriangles*3)
       integer(c_int) :: new_numberofsegments
       integer(c_int) :: new_segmentlist(2*default_numberofsegments)
       integer(c_int) :: new_segmentmarkerlist(default_numberofsegments)

     end subroutine trimesh

  end interface

  ! a type for holding the edges of a triangle
  type tri_edges
     ! edg(edge number  = 1:3, edge pt1, edge pt2)
     integer, dimension(3,2) :: edg
  end type tri_edges


  public :: meshit, trigen, mesh_a_triangle

contains

  ! meshes a region given by the set of connectors 
  ! and returns the results in grid object "grd"
  subroutine meshit(cmd_arg_string_in, pts, cons, holes, report_before_in, report_after_in, grd)
    implicit none
    character (len = *), intent(in) :: cmd_arg_string_in
    type(point), dimension(:), intent(in) :: pts
    type(connector), dimension(:), intent(in) :: cons
    type(point), dimension(:), intent(in) :: holes
    integer, intent(in) :: report_before_in, report_after_in
    type(grid), intent(inout) :: grd


    ! local vars
    character (kind = c_char, len = 200) :: cmd_arg_string
    integer (c_int) :: report_before, report_after
    integer (c_int) :: numberofpoints
    real (c_double), dimension(:), allocatable :: pointlist
    integer (c_int) :: numberofsegments 
    integer (c_int), dimension(:), allocatable  :: segmentlist, segmentmarkerlist
    integer (c_int) :: numberofholes
    real (c_double), dimension(:), allocatable  :: holelist

    integer(c_int), parameter :: default_numberofpoints = 500000
    integer(c_int), parameter :: default_numberoftriangles = 500000
    integer(c_int), parameter :: default_numberofsegments = 500000

    ! outputs
    integer(c_int) :: new_numberofpoints
    real(c_double) :: new_pointlist(2 * default_numberofpoints)
    integer(c_int) :: numberoftriangles
    integer(c_int) :: numberofcorners
    integer(c_int) :: trianglelist(default_numberoftriangles*3)
    integer(c_int) :: neighborlist(default_numberoftriangles*3)
    integer(c_int) :: new_numberofsegments
    integer(c_int) :: new_segmentlist(2*default_numberofsegments)
    integer(c_int) :: new_segmentmarkerlist(default_numberofsegments)

    ! indices
    integer :: ii, telem

    ! initialization
    cmd_arg_string = cmd_arg_string_in//C_NULL_CHAR 
    report_before = report_before_in; report_after = report_after_in

    numberofpoints = size(pts)
    allocate(pointlist(2*numberofpoints))
    numberofsegments = size(cons) 
    allocate(segmentlist(2*numberofsegments), segmentmarkerlist(numberofsegments))
    numberofholes = size(holes)
    allocate(holelist(2*numberofholes))


    ! cmd_arg_string = 'punq35.0'
    !cmd_arg_string = 'pn'

    ! filling out point list
    do ii = 1, numberofpoints
       pointlist(2 * ii - 1) = pts(ii)%x
       pointlist(2 * ii    ) = pts(ii)%y
    end do

    ! filling out segment list
    do ii = 1, numberofsegments
       segmentlist(2 * ii - 1) = cons(ii)%pt1%tag
       segmentlist(2 * ii    ) = cons(ii)%pt2%tag
       ! filling out segment marker list
       segmentmarkerlist(ii) = cons(ii)%tag
    end do

    do ii = 1, size(holes)
       holelist( 2 * ii - 1 ) = holes(ii)%x
       holelist( 2 * ii     ) = holes(ii)%y
    end do

    call trimesh(cmd_arg_string, numberofpoints, pointlist, &
         numberofsegments, segmentlist, &
         segmentmarkerlist, numberofholes, &
         holelist, report_before, report_after, &
         default_numberofpoints, default_numberoftriangles, default_numberofsegments, &
         new_numberofpoints, new_pointlist, &
         numberoftriangles, numberofcorners, trianglelist, neighborlist, &
         new_numberofsegments, new_segmentlist, &
         new_segmentmarkerlist)


    ! print *, 'new_numberofpoints = ', new_numberofpoints
    ! print *, 'new_pointlist = ', new_pointlist(1:(2*new_numberofpoints))
    ! print *, 'numberofcorners = ', numberofcorners  
    ! print *, 'trianglelist = ', trianglelist(1:(numberofcorners*numberoftriangles))
    ! print *, 'neighborlist = ', neighborlist(1:(3*numberoftriangles))
    ! print *, 'new_numberofsegments = ', new_numberofsegments
    ! print *, 'new_segmentlist = ', new_segmentlist(1:(2*new_numberofsegments))
    ! print *, 'new_segmentmarkerlist = ', new_segmentmarkerlist(1:new_numberofsegments)

    ! bug checking
    if (numberoftriangles .eq. 0) then
       print *,'error : no triangle is generated! stopped!'
       stop
    end if

    if (numberofcorners > 3) then
       print *, 'error : higher order triangular elements not implemented in this version'
       print *, 'check meshit(...) subroutine for more info. stopped'
       stop
    end if

    ! lets put everything back into the grid object
    grd%nnodesg = new_numberofpoints
    grd%ncellsg = numberoftriangles
    grd%nbedgeg = new_numberofsegments
    grd%ntri = grd%ncellsg
    grd%nquad4 = 0
    allocate(grd%x(new_numberofpoints), grd%y(new_numberofpoints))
    do ii = 1, new_numberofpoints
       grd%x(ii) = new_pointlist(2 * ii - 1)
       grd%y(ii) = new_pointlist(2 * ii    )
    end do

    allocate(grd%icon(numberoftriangles, 3))
    do ii = 1, numberoftriangles
       grd%icon(ii, 1) = trianglelist(3*ii-2)
       grd%icon(ii, 2) = trianglelist(3*ii-1)
       grd%icon(ii, 3) = trianglelist(3*ii  )
    end do

    allocate(grd%ibedge(new_numberofsegments,2))
    do ii = 1, new_numberofsegments
       grd%ibedge(ii,1) = new_segmentlist( 2 * ii - 1)
       grd%ibedge(ii,2) = new_segmentlist( 2 * ii    )
    end do
    allocate(grd%ibedgeBC(new_numberofsegments))
    do ii = 1, new_numberofsegments
       grd%ibedgeBC(ii) = new_segmentmarkerlist(ii)
    end do

    ! implement ibedgeELEM ??? here
    call bnedge2elem(grd)
    grd%meshFile = 'nofile'

    ! add e2e
    allocate(grd%e2e(numberoftriangles, 3))
    do ii = 1, numberoftriangles
       grd%e2e(ii, 1) = neighborlist(3*ii-2)
       grd%e2e(ii, 2) = neighborlist(3*ii-1)
       grd%e2e(ii, 3) = neighborlist(3*ii  )
    end do

    ! initializing el2bn
    if ( .not. allocated(grd%el2bn) ) allocate(grd%el2bn(grd%ncellsg,2))
    grd%el2bn = 0 ! all assumed interior elements initially!
    ! and then we add el2bn maps for boundary elements accordingly
    do ii = 1, grd%nbedgeg
       telem = grd%ibedgeELEM(ii)
       grd%el2bn(telem, 1) = grd%ibedgeBC(ii)
       grd%el2bn(telem, 2) = ii 
    end do

    ! find duplicated nodes and store them.
    ! this will be used to impose BCs
    call find_dup_nodes(grd%x, grd%y, 1.0d-15, grd%dup_nodes)
    if (allocated(grd%dup_nodes)) then
       grd%tot_repeated_bn_nodes &
            = size(grd%dup_nodes) / 2
    else
       grd%tot_repeated_bn_nodes = 0
    end if

    ! done here
  end subroutine meshit

  ! 1- finds the boundary edge to element map
  ! 2- fix winding of the bnedges according
  !    to the winding of the elements.  
  subroutine bnedge2elem(grd)
    implicit none
    type(grid), target, intent(inout) :: grd

    ! local vars
    integer :: i, j, l, k, pt1, pt2
    integer, dimension(:,:), pointer :: icon, ibedge
    logical :: found1, found2
    integer :: l1, l2 

    type(tri_edges), dimension(:), target, allocatable :: tris_edges
    type(tri_edges), pointer :: this_one

    if (.not. allocated(grd%ibedgeELEM) ) then
       allocate(grd%ibedgeELEM(grd%nbedgeg), grd%ibedgeELEM_local_edg(grd%nbedgeg))
    end if

    icon => grd%icon
    ibedge => grd%ibedge

    print *, 'creating bnedge2elem map ...'

    ! error checking
    if ( grd%ntri .eq. 0 ) then
       print *, 'the FEM grid does not have triangle!!! stopped'
       stop
    end if

    ! first collecting the edges from triangles for real correct winding
    allocate(tris_edges(grd%ntri))

    do j = 1, grd%ntri
       this_one => tris_edges(j) 

       this_one%edg(1,1) = icon(j,1); this_one%edg(1,2) = icon(j,2)
       this_one%edg(2,1) = icon(j,2); this_one%edg(2,2) = icon(j,3)
       this_one%edg(3,1) = icon(j,3); this_one%edg(3,2) = icon(j,1)
    end do

    ! STARTING THE BRUTE FORCE METHOD
    ! 
    do i = 1, grd%nbedgeg ! over all boundary edges
       pt1 = ibedge(i,1); pt2 = ibedge(i,2)

       do j = 1, grd%ntri ! over all tri

          do k = 1, 3 ! over edges in that triangle
             found1 = .false.; found2 = .false.

             do l = 1, 2 ! over point of that edge

                if (tris_edges(j)%edg(k,l) .eq. pt1) then
                   found1 = .true.
                   l1 = l
                end if

                if (tris_edges(j)%edg(k,l) .eq. pt2) then
                   found2 = .true.
                   l2 = l
                end if

             end do ! l

             if ( found1 .and. found2) then ! both points were found!
                grd%ibedgeELEM(i) = j ! saving the element number
                if(l1 > l2 ) then ! fix the winding 
                   print *, 'warning : winding of the boundary edge ', i, 'was wrong. fixed successfully!'
                   ! print*, grd%ibedge(i,1) , '<=', tris_edges(j)%edg(k,l2)
                   grd%ibedge(i,1) = tris_edges(j)%edg(k,l2)
                   ! print*, grd%ibedge(i,2) , '<=', tris_edges(j)%edg(k,l1)
                   grd%ibedge(i,2) = tris_edges(j)%edg(k,l1)
                end if
                exit
             end if

          end do ! k - edges of a single tri

          ! don't search for more triangles
          ! already found and fixed that boundary edge
          if ( found1 .and. found2) exit

       end do ! j - all tris

    end do ! i - all bn edges

    ! fill ibedgeELEM_local_edg
    call fill_ibedgeELEM_local_edg(grd)

    ! clean up
    deallocate( tris_edges )

    ! done here
  end subroutine bnedge2elem

  ! creates triangular mesh for a geometry
  ! read from segment file.
  !
  subroutine trigen(method, grd)
    implicit none
    character(len = *), intent(in) :: method
    type(grid), target, intent(inout) :: grd

    ! local vars
    integer :: i, j, indx, npts, ncons, ncurves, indx2
    type(point), dimension(:), allocatable, target :: pts
    type(connector), dimension(:), allocatable :: cons
    type(point) :: holes(1)
    real*8, dimension(:), pointer :: x, y
    logical :: begin
    integer :: npe, nped

    if (.not. allocated(grd%bn_curves)) then
       print *, 'first read segments before trigen! call '&
            ,' read_segment_file subroutine before this. stop.'
       stop
    end if

    ncurves = size(grd%bn_curves)
    npts = grd%tot_repeated_bn_nodes
    ncons = grd%tot_bn_seg

    allocate(pts(npts), cons(ncons))

    indx = 1
    indx2 = 1
    do i = 1, ncurves 
       x => grd%bn_curves(i)%x
       y => grd%bn_curves(i)%y
       begin = .true.

       do j = 1, size(x)

          pts(indx)%x = x(j)
          pts(indx)%y = y(j)
          pts(indx)%tag = indx

          if (.not. begin) then
             cons(indx2)%tag = i
             cons(indx2)%pt1 => pts(indx-1)
             cons(indx2)%pt2 => pts(indx)   
             cons(indx2)%is_bn_con = .true.
             indx2 = indx2 + 1
          end if

          indx = indx + 1
          begin = .false.

       end do

    end do

    ! triangulate the fem region
    holes(1)%tag = -1
    holes(1)%x = 0.0d0
    holes(1)%y = 0.0d0
    ! holes(2)%tag = -2
    ! holes(2)%x =  2.41d0
    ! holes(2)%y = -0.94d0

    ! holes(1)%tag = -1
    ! holes(1)%x = 0.921496d0
    ! holes(1)%y = -0.006996d0
    ! holes(2)%tag = -2
    ! holes(2)%x =  0.3335193d0
    ! holes(2)%y = -0.001835d0
    ! holes(3)%tag = -3
    ! holes(3)%x =  -0.06288048d0
    ! holes(3)%y = -0.0883964d0

    ! holes(1)%tag = -1
    ! holes(1)%x = 500.0d0
    ! holes(1)%y = 0.0d0

    ! holes(1)%tag = -1
    ! holes(1)%x = 0.01011d0
    ! holes(1)%y = 0.00406d0

    ! holes(2)%tag = -2
    ! holes(2)%x = 0.006244d0
    ! holes(2)%y = -0.0006355d0

    ! holes(3)%tag = -3
    ! holes(3)%x = 0.014008d0
    ! holes(3)%y = -0.00374d0

    ! holes(4)%tag = -4
    ! holes(4)%x = 0.01905d0
    ! holes(4)%y = -0.00983d0

    ! holes(5)%tag = -5
    ! holes(5)%x = 0.02311d0
    ! holes(5)%y = -0.01646d0

    ! holes(6)%tag = -6
    ! holes(6)%x = 0.02626d0
    ! holes(6)%y = -0.02255d0

    call meshit(method, pts, cons, holes, 1, 1, grd)

    ! bullet proofing 
    if ( grd%ncellsg .eq. 0 ) then
       print *, 'no triangle was generated by trimesher.' &
            , ' error happend in trigen subroutine! stop.'
       stop
    end if

    ! setup the order of the elements
    if( .not. allocated(grd%p) ) allocate(grd%p(grd%ncellsg))

    grd%p = 1  ! initially first-order
               ! this is just the skeleton
               ! of the element.

    ! setup the type of the elements
    if( .not. allocated(grd%eltype) ) allocate(grd%eltype(grd%ncellsg))
    grd%eltype = 0 ! initially all elements are Lagrange 0
                   ! then after this we can add more Chebyshev or fekete points
                   ! to have different element types 
    ! setup the number of points per elem and edges
    ! default is npe = 3 nped = 2 for skeleton (linear triangle)
    call p2node(grd%p(1), npe, nped) ! lets compute it!

    if ( .not. allocated(grd%npe)  ) allocate(grd%npe(grd%ncellsg))
    grd%npe = npe

    if ( .not. allocated(grd%el2edg)  ) allocate(grd%el2edg(grd%ncellsg))

    ! taking a snapshot of initial grid config.
    ! so we can retrieve original unadapted grid back
    grd%nnodesg0 = grd%nnodesg
    grd%ncellsg0 = grd%ncellsg
    grd%nbedgeg0 = grd%nbedgeg 

    ! initializing master element information
    ! for the default 3 points linear element
    if ( .not. allocated(grd%maselem) ) then
       allocate(grd%maselem(grd%ncellsg))

       do i = 1, grd%ncellsg
          allocate(grd%maselem(i)%xi(3), grd%maselem(i)%eta(3))
          grd%maselem(i)%xi  = (/ 0.0d0, 1.0d0, 0.0d0 /)
          grd%maselem(i)%eta = (/ 0.0d0, 0.0d0, 1.0d0 /) 
       end do

    end if

    ! done here
  end subroutine trigen

  ! meshes a single high-order curvelinear triangle
  ! with element number "ielem" in input grid "grd_in"
  ! and writes the resulting grid in "grd_out" 
  subroutine mesh_a_triangle(grd_in, ielem, grd_out)
    implicit none
    type(grid), intent(in) :: grd_in
    integer, intent(in) :: ielem
    type(grid), intent(out) :: grd_out

    ! local vars
    integer :: i, i1, i2, indx, npts, pts_per_edg, cons_per_edg, dpt
    type(point), dimension(:), allocatable :: pts
    type(connector), dimension(:), allocatable :: cons 
    type(point) :: holes(1)
    integer, dimension(:), allocatable :: pts_list, int_list
    integer :: report_before_in, report_after_in

    ! init.
    report_before_in = 0 
    report_after_in = 0

    ! ! refresh grd_out structure
    ! if ( allocated(grd_out) ) deallocate(grd_out)


    ! first fill the points
    npts = grd_in%npe(ielem)
    allocate(pts(npts))

    do i = 1, npts
       pts(i)%x = grd_in%x(grd_in%icon(ielem, i))
       pts(i)%y = grd_in%y(grd_in%icon(ielem, i))
       pts(i)%tag = i
    end do

    ! then dimension all edges
    pts_per_edg = size(grd_in%el2edg(ielem)%edg1) + 2 ! two for end points
    cons_per_edg = pts_per_edg - 1
    dpt = size(grd_in%el2edg(ielem)%edg1)

    ! proceeding to fill connectors array
    allocate(cons(3 * cons_per_edg))
    indx = 1

    ! -----------------------------------
    !                EDG1
    ! -----------------------------------
    ! adding connectors on edg1 of the element
    if (pts_per_edg .eq. 2) then
       pts_list = (/ 1, 2 /)
    else
       i1 = 4
       i2 = 4 + dpt - 1
       int_list = (/ (i, i = i1, i2) /)
       pts_list = (/ 1, int_list, 2 /)
    end if
    call add_to_cons(pts, pts_list, indx, cons)


    ! -----------------------------------
    !                EDG2
    ! -----------------------------------
    ! adding connectors on edg2 of the element
    if (pts_per_edg .eq. 2) then
       pts_list = (/ 2, 3 /)
    else
       i1 = 4 + dpt
       i2 = 4 + 2* dpt - 1
       int_list = (/ (i, i = i1, i2) /)
       pts_list = (/ 2, int_list, 3 /)
    end if
    call add_to_cons(pts, pts_list, indx, cons)

    ! -----------------------------------
    !                EDG3
    ! -----------------------------------
    ! adding connectors on edg3 of the element
    if (pts_per_edg .eq. 2) then
       pts_list = (/ 3, 1 /)
    else
       i1 = 4 + 2* dpt
       i2 = 4 + 3* dpt - 1
       int_list = (/ (i, i = i1, i2) /)
       pts_list = (/ 3, int_list, 1 /)
    end if
    call add_to_cons(pts, pts_list, indx, cons)


    ! setting up the holes which there is no hole
    ! here in this case
    holes(1)%tag = -1
    holes(1)%x = 1.0d20 ! some point outside the domain
    holes(1)%y = 1.0d20 ! some point outside the domain

    !
    call meshit('pnjY', pts, cons, holes, report_before_in &
                      , report_after_in, grd_out)

    ! clean ups
    deallocate(pts, cons)

    ! done here

  contains

    subroutine add_to_cons(pts, pts_list, indx, cons)
      implicit none
      type(point), dimension(:), target, intent(in) :: pts
      integer, dimension(:), intent(in) :: pts_list
      integer, intent(inout) :: indx 
      type(connector), dimension(:), intent(inout):: cons 

      ! local vars
      integer :: i , leng, pt1, pt2 

      leng = size(pts_list)

      do i = 1, (leng-1)

         pt1 = pts_list( i )
         pt2 = pts_list(i+1)

         cons(indx)%tag = 1
         cons(indx)%pt1 => pts(pt1)
         cons(indx)%pt2 => pts(pt2)   
         cons(indx)%is_bn_con = .true.
         indx = indx + 1

      end do

      ! done here
    end subroutine add_to_cons

  end subroutine mesh_a_triangle

end module trimesher

! program testit
! use trimesher
! use iso_c_binding, only : c_int, c_char, c_double
! implicit none

!   ! local vars
!   character (kind = c_char, len = 20) :: cmd_arg_string
!   integer (c_int), parameter :: numberofpoints = 8
!   real (c_double) :: pointlist(2*numberofpoints)
!   integer (c_int), parameter :: numberofsegments = 8 
!   integer (c_int) :: segmentlist(2*numberofsegments), segmentmarkerlist(numberofsegments)
!   integer (c_int), parameter :: numberofholes = 1
!   real (c_double) :: holelist(2*numberofholes)
!   integer (c_int) :: report_before, report_after
!   integer(c_int), parameter :: default_numberofpoints = 500000
!   integer(c_int), parameter :: default_numberoftriangles = 500000
!   integer(c_int), parameter :: default_numberofsegments = 500000

!   ! outputs
!   integer(c_int) :: new_numberofpoints
!   real(c_double) :: new_pointlist(2 * default_numberofpoints)
!   integer(c_int) :: numberoftriangles
!   integer(c_int) :: numberofcorners
!   integer(c_int) :: trianglelist(default_numberoftriangles*3)
!   integer(c_int) :: neighborlist(default_numberoftriangles*3)
!   integer(c_int) :: new_numberofsegments
!   integer(c_int) :: new_segmentlist(2*default_numberofsegments)
!   integer(c_int) :: new_segmentmarkerlist(default_numberofsegments)

!   cmd_arg_string = 'punq36.1'
!   cmd_arg_string = 'pn'
!   pointlist = (/ 0.0d0, 0.0d0, 1.0d0, 0.0d0, 1.0d0, 1.0d0, 0.0d0, 1.0d0, &
!        0.25d0, 0.25d0, 0.75d0, 0.25d0, 0.75d0, 0.75d0, 0.25d0, 0.75d0 /)
!   segmentlist = (/ 1,2,2,3,3,4,4,1,5,6,6,7,7,8,8,5 /)
!   segmentmarkerlist = (/ 1,1,1,1,2,2,2,2 /)
!   holelist = (/ 0.5d0, 0.5d0 /)
!   report_before = 0; report_after = 1

!   print *, default_numberofpoints, default_numberoftriangles, default_numberofsegments

! call trimesh(cmd_arg_string, numberofpoints, pointlist, &
!                        numberofsegments, segmentlist, &
!                        segmentmarkerlist, numberofholes, &
!                        holelist, report_before, report_after, &
!                        default_numberofpoints, default_numberoftriangles, default_numberofsegments, &
!                        new_numberofpoints, new_pointlist, &
!                        numberoftriangles, numberofcorners, trianglelist, neighborlist, &
!                        new_numberofsegments, new_segmentlist, &
!                        new_segmentmarkerlist)

!   ! call trimesh(cmd_arg_string, numberofpoints, pointlist, &
!   !      numberofsegments, segmentlist, &
!   !      segmentmarkerlist, numberofholes, &
!   !      holelist, report_before, report_after)
! print *, 'new_numberofpoints = ', new_numberofpoints
! print *, 'new_pointlist = ', new_pointlist(1:(2*new_numberofpoints))
! print *, 'numberofcorners = ', numberofcorners  
! print *, 'trianglelist = ', trianglelist(1:(numberofcorners*numberoftriangles))
! print *, 'neighborlist = ', neighborlist(1:(3*numberoftriangles))
! print *, 'new_numberofsegments = ', new_numberofsegments
! print *, 'new_segmentlist = ', new_segmentlist(1:(2*new_numberofsegments))
! print *, 'new_segmentmarkerlist = ', new_segmentmarkerlist(1:new_numberofsegments)


!   print *, 'done!'
!   stop


! end program testit
