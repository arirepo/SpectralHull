module grd2hull
  ! version 2.00; %x is now local to each edge
  use grid_opt
  use approx_fekete
  implicit none

  private

  type ej

     ! (1=x;2=y , 1=pt1;2=pt2) 
     real*8, dimension(2,2) :: x
     integer :: bc, neigh

  end type ej


  type hull

     logical :: is_active, is_bn
     type(ej), dimension(:), allocatable :: ejs


  end type hull

  type hulls

     type(hull), dimension(:), allocatable :: hl

  end type hulls

  real*8, parameter :: PI_NUM = 4.D0*DATAN(1.D0)

  public :: hull, hulls
  public :: convert_grd_to_hull, write_hulls_gnuplot
  public :: agglomerate_hulls, hull2array, hull2edgs
  public :: print_hulls

contains

  ! converts a standard type(grid) object to hull-based representation
  ! implemented in type(hulls) class.
  ! only call trigen or trimesh or quadgen before this
  ! do not call add_more_points in either PG or DG form before
  ! this. This is only based on the skleton points (main vertices)
  ! so it is better to never call add_more_points before this!
  !
  subroutine convert_grd_to_hull(grd, hls)
    implicit none
    type(grid), intent(in) :: grd
    type(hulls), intent(out), target :: hls

    ! local vars
    integer :: i, j, pt1, pt2, max_ejs
    type(hull), pointer :: thull => null()

    ! bullet proofing ...
    if ( grd%nnodesg <= 0 ) then
       print *, 'before converting to hull, the grid object first needs to be initialized! stop'
       stop
    end if

    ! init all hulls
    allocate(hls%hl(size(grd%icon, 1)))

    ! fill out each hull 
    do i = 1, size(hls%hl)

       thull => hls%hl(i) !this hull

       ! select maximum number of the edges
       ! for the given element in the input grid
       select case (grd%elname(i))
       case ( GEN_TRIANGLE )
          max_ejs = 3
       case ( GEN_QUADRI )
          max_ejs = 4
       case default
          print *, 'unknown name of the element in converting grd to hull! stop'
          stop
       end select

       ! set it as an active hull if it is not a boundary element
       if ( any(grd%local_edg_bc(i, 1:max_ejs) .ne. 0) ) then !boundary hull
          thull%is_active = .false.
          thull%is_bn = .true.
       else
          thull%is_active = .true.
          thull%is_bn = .false.
       end if

       ! allocate and fill ejes of each hull
       allocate(thull%ejs(max_ejs))
       do j = 1, size(thull%ejs)

          ! add point numbers to ejs
          pt1 = grd%icon(i, j)
          if (j .eq. size(thull%ejs) ) then
             pt2 = grd%icon(i, 1)
          else
             pt2 = grd%icon(i, (j+1))
          end if
          thull%ejs(j)%x(:,1) = (/ grd%x(pt1), grd%y(pt1) /)
          thull%ejs(j)%x(:,2) = (/ grd%x(pt2), grd%y(pt2) /)

          ! detect neighbor per that edge 
          if ( j .eq. 1 ) then
             thull%ejs(j)%neigh = grd%e2e(i, max_ejs)
          else
             thull%ejs(j)%neigh = grd%e2e(i, (j-1))
          end if

          ! specify the bc of this edge
          thull%ejs(j)%bc = grd%local_edg_bc(i, j)

       end do !over edges of each hull

    end do ! over all hulls

    ! done here
  end subroutine convert_grd_to_hull

  ! find local edge number of the given edge 
  ! in the given hull using a brute force search method
  ! and comparing the points number of the two ends of 
  ! the edge.
  ! 
  ! HINT : winding of the edge is not important
  ! 
  function find_local_ej_in_hull(tej, thull)
    implicit none
    type(ej), intent(in) :: tej
    type(hull), intent(in), target :: thull
    integer :: find_local_ej_in_hull 

    ! local vars
    integer :: i
    type(ej), pointer :: lej => null() ! local ej
    real*8, dimension(2,2) :: pt_rev

    ! init
    find_local_ej_in_hull = -1 ! initially no local ej!
    pt_rev(:,1) = tej%x(:,2)
    pt_rev(:,2) = tej%x(:,1)

    do i = 1, size(thull%ejs)

       lej => thull%ejs(i)

       if ( all(lej%x .eq. tej%x) .or. all(lej%x .eq. pt_rev)  ) then
          find_local_ej_in_hull = i
          return
       end if

    end do

    ! edj not found! well this is bad. stop it
    print *, 'the given local edge was not found in the hull! stop'
    stop  

    ! done here
  end function find_local_ej_in_hull

  ! in this hull (thull), the following returns
  ! the local ej number of a shift of either +1 or -1
  ! with respect to the current local ej number, i.e. ejnum
  ! 
  function cyclic_ej_num(thull, ejnum, shift)
    implicit none
    type(hull), intent(in) :: thull
    integer, intent(in) :: ejnum, shift
    integer :: cyclic_ej_num

    ! local vars
    integer :: max_ejs

    max_ejs = size(thull%ejs)

    ! bullet proofing ...
    if ( (ejnum < 1) .or. ( ejnum > max_ejs) ) then
       print *, 'the given ej number is not within the range ' &
            , ' of the given hull! cannot find the shift! stop'
       stop
    end if
    select case (shift)
    case (1, -1)
       ! do nothing, this is normal
    case default
       print *, 'wrong shift value in cyclic_ej_num! stop'
       stop
    end select

    ! now proceed to compute the local edge number of
    ! the result of applying shift to the current ej
    !
    if ( (ejnum .eq. 1) .and. (shift .eq. -1) ) then
       cyclic_ej_num = max_ejs
    elseif ( (ejnum .eq. max_ejs) .and. (shift .eq. 1)) then
       cyclic_ej_num = 1
    else
       cyclic_ej_num = ejnum + shift
    end if

    ! done here
  end function cyclic_ej_num

  ! computes the cross product of
  ! two edges in two dimensions. 
  ! can be positive and negative
  ! this is used to assess the convexity 
  !
  ! if (z = ej1 <cross> ej2) >= 0 then convex
  ! else concave!
  ! 
  function is_convex(ej1, ej2)
    implicit none
    type(ej), intent(in) :: ej1, ej2
    logical :: is_convex

    ! local vars
    real*8 :: x1, y1, x2, y2
    real*8 :: a, b, c, d, z

    ! compute the components of vector ej1 and ej2
    ! NOTE : orientation is also accounted 
    ! ej1 ...
    x1 = ej1%x(1,1); y1 = ej1%x(2,1)
    x2 = ej1%x(1,2); y2 = ej1%x(2,2)

    a = x2 - x1
    b = y2 - y1

    ! ej2 ...
    x1 = ej2%x(1,1); y1 = ej2%x(2,1)
    x2 = ej2%x(1,2); y2 = ej2%x(2,2)

    c = x2 - x1
    d = y2 - y1
    ! compute cross product
    z = a * d - c * b 

    ! final decision
    if ( z > 0.0d0) then
       is_convex = .true.
    else
       is_convex = .false.
    end if

    ! done here
  end function is_convex

  ! given an edge number in a hull
  ! the following starts from next ej 
  ! and gather all ejs sequentially
  ! until it reachs back to one ej 
  ! before the given ej.
  ! the result is packed in "other_ejs"
  subroutine pack_other_ejs(thull, ejnum, other_ejs)
    implicit none
    type(hull), intent(in) :: thull
    integer, intent(in) :: ejnum
    type(ej), dimension(:), allocatable :: other_ejs 

    ! local vars
    integer :: i, max_ejs, tejnum, next_ej

    if ( allocated(other_ejs) ) deallocate(other_ejs)

    max_ejs = size(thull%ejs)
    allocate(other_ejs(max_ejs-1)) ! one will be removed

    tejnum = ejnum
    do i = 1, size(other_ejs)

       next_ej = cyclic_ej_num(thull, tejnum, 1) !counter clockwise +1 shift
       other_ejs(i) = thull%ejs(next_ej) ! full copy

       ! update this current ej
       tejnum = next_ej 

    end do

    ! done here
  end subroutine pack_other_ejs

  ! inserts a pack of "other_ejs" at a location "loc"
  ! of the contigous type array ejs 
  !
  subroutine insert_ejs(ejs, other_ejs, loc)
    implicit none
    type(ej), dimension(:), allocatable :: ejs, other_ejs
    integer, intent(in) :: loc

    ! local vars
    integer :: tot_ejs, i1, i2
    type(ej), dimension(:), allocatable :: tmp

    tot_ejs = size(ejs) - 1 + size(other_ejs)

    allocate(tmp(tot_ejs))
    tmp(1:(loc-1)) = ejs(1:(loc-1))
    i1 = loc
    i2 = loc + size(other_ejs) - 1  
    tmp(i1:i2) = other_ejs
    i1 = i2+1
    i2 = tot_ejs  
    if (loc < size(ejs)) tmp(i1:i2) = ejs((loc+1):size(ejs))

    call move_alloc(tmp, ejs)

    ! clean ups
    if( allocated(tmp)) deallocate(tmp)

    ! done here
  end subroutine insert_ejs

  ! try to agglomerate hull numbered 1
  ! with the neighboring hull to local edge loc_ej
  ! of the hull 1.
  ! if successfully agglomerated then returns .true.
  ! lse return .false.
  !  
  function try_agglomerate(thulls, hullnum1, loc_ej)
    implicit none
    type(hulls), intent(inout), target :: thulls
    integer, intent(in) :: hullnum1, loc_ej
    logical :: try_agglomerate


    ! local vars
    integer :: i, ej1num, ej2num, hullnum2, loc2, loc2_bak
    type(ej), pointer :: ej1 => null(), ej2 => null(), tej => null()
    type(hull), pointer :: thull1 => null(), neigh_hull => null(), second_neigh => null()
    type(ej), dimension(:), allocatable :: other_ejs

    ! init
    thull1 => thulls%hl(hullnum1)
    tej => thull1%ejs(loc_ej)
    hullnum2 = tej%neigh

    ! bullet proofing
    if ( hullnum2 .eq. -1 ) then !wall
       try_agglomerate = .false.
       return
    end if

    if ( .not. thulls%hl(hullnum2)%is_active) then !inactive hull
       try_agglomerate = .false.
       return
    end if

    if ( (size(thull1%ejs) >= 4) .or. (size(thulls%hl(hullnum2)%ejs) >= 4) ) then
       try_agglomerate = .false.
       return
    end if

    do i = 1, size(thull1%ejs)
       if ( any(thull1%ejs(i)%x(1,:) < .5d0) ) then
             try_agglomerate = .false.
             return
       end if
    end do

    do i = 1, size(thulls%hl(hullnum2)%ejs)
       if ( any(thulls%hl(hullnum2)%ejs(i)%x(1,:) < .5d0) ) then
             try_agglomerate = .false.
             return
       end if
    end do

    ej1num = cyclic_ej_num(thull1, loc_ej, -1)
    ej1 => thull1%ejs(ej1num)

    neigh_hull => thulls%hl(hullnum2)
    loc2 = find_local_ej_in_hull(tej, neigh_hull)
    loc2_bak = loc2
    ej2num = cyclic_ej_num(neigh_hull, loc2, 1)
    ej2 => neigh_hull%ejs(ej2num)

    if ( (.not. is_convex(ej1, ej2)) &
         .or. (.not. is_angle_good(ej1, ej2, 10.0d0)) ) then
       try_agglomerate = .false.
       return
    end if

    ej1num = cyclic_ej_num(neigh_hull, loc2, -1)
    ej1 => neigh_hull%ejs(ej1num)

    ej2num = cyclic_ej_num(thull1, loc_ej, 1)
    ej2 => thull1%ejs(ej2num)

    if ( (.not. is_convex(ej1, ej2)) &
         .or. (.not. is_angle_good(ej1, ej2, 10.0d0)) ) then
       try_agglomerate = .false.
       return
    end if

    ! allright, we proceed to agglomerate

    ! update the neighbor info in the neighbor hulls 
    ! of the hull to me eliminated, i.e. hull2 
    do i = 1, size(neigh_hull%ejs)
       tej => neigh_hull%ejs(i)
       if ( tej%neigh .eq. -1 ) cycle ! wall
       second_neigh => thulls%hl(tej%neigh)

       loc2 = find_local_ej_in_hull(tej, second_neigh)
       second_neigh%ejs(loc2)%neigh = hullnum1 
    end do

    ! then pack other ejs
    call pack_other_ejs(neigh_hull, loc2_bak, other_ejs)

    ! then insert
    call insert_ejs(thull1%ejs, other_ejs, loc_ej)

    ! finally make the neighbor hull deactive 
    neigh_hull%is_active = .false.

    ! little clean up
    if ( allocated(other_ejs) ) deallocate(other_ejs)

    ! final stamp
    try_agglomerate = .true.

    ! done here
  end function try_agglomerate

  ! connects neighboring hulls if they are 
  ! satisfying a condition like convexity
  !
  subroutine agglomerate_hulls(thulls, maxtry)
    implicit none
    type(hulls), intent(inout), target :: thulls
    integer, intent(in) :: maxtry

    ! local vars
    integer :: i, j, randej, max_ejs
    real*8 :: rnd
    type(hull), pointer :: thull => null()
    logical :: success

    do i = 1, size(thulls%hl)

       thull => thulls%hl(i)
       if ( .not. thull%is_active ) cycle

       success = .false.
       do j = 1, maxtry
          max_ejs = size(thull%ejs)
          call random_number(rnd)
          randej = ceiling(rnd * dble(max_ejs))
          success = try_agglomerate(thulls, i, randej)
       end do !trys

    end do ! agglomeration of all hulls

    ! only select active or boundary hulls 
    ! and update numbering
    call elim_inact_hulls(thulls)

    ! done here
  end subroutine agglomerate_hulls

  subroutine write_hulls_gnuplot(fname, thulls)
    implicit none
    character(len=*), intent(in) :: fname
    type(hulls), intent(in), target :: thulls

    ! local vars
    integer :: i, j
    type(ej), dimension(:), pointer :: ejs => null()


    ! init
    open(unit = 10, file = fname, status = 'unknown')

    do i = 1, size(thulls%hl)

       if ( (.not. thulls%hl(i)%is_active) .and. (.not. thulls%hl(i)%is_bn) ) cycle

       ejs => thulls%hl(i)%ejs

       do j = 1, size(ejs)

          write(10, *) ejs(j)%x(1, 1), ejs(j)%x(2, 1)

          if ( j .eq. size(ejs) ) then
             write(10, *) ejs(j)%x(1, 2), ejs(j)%x(2, 2)
          end if

       end do

       write(10, *)

    end do ! end of all hulls

    !
    close(10)

    ! done here
  end subroutine write_hulls_gnuplot

  ! converts a given hull to an array
  ! containing the nodes connected sequntially
  ! in conterclockwise direction and the last node
  ! is repeated.
  !
  ! for a polygon(hull) that has "l" sides, "ar" is an
  ! array containing its vertices, ordered counterclockwise.
  ! as last row must have the components of the first vertex.
  ! in other words, the first row and last row are equal.
  ! Therefore "ar" is a "l+1 x 2" matrix.
  !
  subroutine hull2array(hl, ar)
    implicit none
    class(hull), intent(in) :: hl
    real*8, dimension(:, :), allocatable :: ar

    ! local vars
    integer :: i, nrow

    ! determine the number of rows 
    ! of the output array "ar"
    nrow = size(hl%ejs) + 1 ! one repeated vertex

    ! allocate the output array
    if ( allocated(ar) ) deallocate(ar)
    allocate(ar(nrow, 2))

    ! packing the vertices of hull into array
    do i = 1, size(hl%ejs)
       ar(i, :) = hl%ejs(i)%x(:, 1)
    end do

    ! inserting the last repeated node
    ar(nrow, :) = hl%ejs(size(hl%ejs))%x(:, 2)

    ! final check
    if ( any( ar(1, :) .ne. ar(nrow, :) ) ) then
       print *, 'the given hull could not be converted to array! stop'
       stop
    end if

    ! done here
  end subroutine hull2array

  ! converts an "hull" class to "edgs" class 
  ! defined in approximate Fekete module
  subroutine hull2edgs(hl, tedgs)
    implicit none
    class(hull), intent(in), target :: hl
    class(edgs), intent(out) :: tedgs

    ! local vars
    integer :: i
    type(edg) :: tedg
    type(ej), pointer :: tej => null()

    do i = 1, size(hl%ejs)
       tej => hl%ejs(i)

       tedg%x1 = tej%x(1,1); tedg%y1 = tej%x(2,1)
       tedg%x2 = tej%x(1,2); tedg%y2 = tej%x(2,2)
       tedg%tag = tej%bc

       call tedgs%add(tedg)

    end do

    ! final settings
    tedgs%maxtag = 1
    tedgs%tol = 2.6d-2

    ! done here
  end subroutine hull2edgs

  ! eliminates inactive interior 
  ! hulls and renumber the elements
  ! and connectivities based on
  ! a contigous (all-active) numbering
  !
  subroutine elim_inact_hulls(hls)
    implicit none
    type(hulls), intent(inout), target :: hls

    ! local vars
    integer :: i, hnum, j, k
    class(hull), pointer :: thull => null()
    type(hull), dimension(:), allocatable :: tmp
    class(ej), pointer :: tej => null()

    ! counting the number of hulls
    hnum = 0
    do i = 1, size(hls%hl)
       thull => hls%hl(i)
       if ( (.not. thull%is_active) &
            .and. (.not. thull%is_bn) ) cycle
       hnum = hnum + 1
    end do

    ! allocate final hulls
    allocate(tmp(hnum))

    ! renumber accordingly ...
    hnum = 1
    do i = 1, size(hls%hl)

       ! skip inactive interior hulls
       thull => hls%hl(i)
       if ( (.not. thull%is_active) &
            .and. (.not. thull%is_bn) ) cycle

       ! update this hull number in all other hulls
       do j = 1, size(hls%hl)

          thull => hls%hl(j)
          if ( (.not. thull%is_active) &
               .and. (.not. thull%is_bn) ) cycle

          do k = 1, size(thull%ejs)
             tej => thull%ejs(k)

             if (tej%neigh .eq. i) tej%neigh = hnum

          end do

       end do

       ! increase the hull number
       hnum = hnum + 1 

    end do


    ! copy the active or boundary hulls
    ! into the new temporary hull buffer
    hnum = 1
    do i = 1, size(hls%hl)

       ! skip inactive interior hulls
       thull => hls%hl(i)
       if ( (.not. thull%is_active) &
            .and. (.not. thull%is_bn) ) cycle

       tmp(hnum) = thull

       ! increase the hull number
       hnum = hnum + 1 

    end do

    ! move the temp array into the final result
    call move_alloc(tmp , hls%hl)

    ! final clean ups
    if ( allocated(tmp) ) deallocate(tmp)
    if ( associated(thull)) nullify(thull)
    if ( associated(tej)) nullify(tej)

    ! done here
  end subroutine elim_inact_hulls

  ! prints info about the hulls
  !
  subroutine print_hulls(hls)
    implicit none
    type(hulls), intent(in), target :: hls

    ! local vars
    integer :: i, j
    class(hull), pointer :: thull => null()
    class(ej), pointer :: tej => null()

    do i = 1, size(hls%hl)

       thull => hls%hl(i)

       print *, '-----------------------------------------------------'
       print *, 'hull # ', i , ' has ', size(thull%ejs) , ' edges! '
       print *, 'the activation status : ', thull%is_active

       if (thull%is_bn) then
          print *, 'this is a boundary hull!'
       else
          print *, 'this is an interior hull!'
       end if

       ! print edges
       do j = 1, size(thull%ejs)
          tej => thull%ejs(j)
          print *, 'edge # ', j, ' start : (', tej%x(:, 1) &
               , ') and end : (', tej%x(:, 2), ').'
          print *, 'ege%bc = ', tej%bc, 'ege%neigh = ', tej%neigh     

       end do

       print *, '****************************************************'

    end do

    ! cleanups
    if ( associated(thull)) nullify(thull)
    if ( associated(tej)) nullify(tej)

    ! done here
  end subroutine print_hulls

  ! check the angle between two edges using
  ! dot product and if satisfies the given angle
  ! in degrees, then return .true. otherwise
  ! returns .false.
  !
  function is_angle_good(ej1, ej2, ang)
    implicit none
    type(ej), intent(in) :: ej1, ej2
    real*8, intent(in) :: ang
    logical :: is_angle_good

    ! local vars
    real*8 :: ax, ay, bx, by
    real*8 :: a_dot_b, norm_a, norm_b, ang0

    ! "a" is alias for edges 1
    ax = ej1%x(1, 2) - ej1%x(1, 1)
    ay = ej1%x(2, 2) - ej1%x(2, 1)

    ! "b" is alias for edges 2
    bx = ej2%x(1, 2) - ej2%x(1, 1)
    by = ej2%x(2, 2) - ej2%x(2, 1)

    ! compute dot product
    a_dot_b = ax*bx + ay*by
    norm_a = sqrt(ax*ax + ay*ay)
    norm_b = sqrt(bx*bx + by*by)

    ! compute the angle between edges
    ang0 = (180.0d0 / PI_NUM) * dacos(a_dot_b / (norm_a * norm_b))

    ! decide on angle
    if ( (ang0 < ang) .or. ( (180.0d0 - ang0) < ang ) ) then
       is_angle_good = .false.
    else
       is_angle_good = .true.
    end if

    ! done here
  end function is_angle_good

end module grd2hull
