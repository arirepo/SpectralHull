module plot_curved_elems_dg
  use grid_opt
  use ps
  use element_opt_dg2d
  use trimesher
  use dg_workspace
  use renka_trimesh
  implicit none

  private


  public :: vis_curved_grid_dg, vis_curved_grid_dg_renka
  public :: vis_curved_grid_dg2

contains


  ! meshes a single high-order curvelinear 
  ! Discontinious Galerkin element
  ! and writes the resulting grid in "grd_out" 
  subroutine mesh_a_dg_elem(elem, grd_out)
    implicit none
    type(element_dg2d), intent(in) :: elem
    type(grid), intent(out) :: grd_out

    ! local vars
    integer :: ii, jj, indx, cons_per_edg, int1, int2, last_pt, pt1, pt2
    type(point), dimension(:), allocatable :: pts
    type(connector), dimension(:), allocatable :: cons 
    type(point) :: holes(1)
    integer, dimension(:), allocatable :: pts_list
    integer :: report_before_in, report_after_in

    ! init.
    report_before_in = 0 
    report_after_in = 0

    ! first fill the points
    allocate(pts(elem%npe))

    do ii = 1, elem%npe
       pts(ii)%x = elem%x(1, ii)
       pts(ii)%y = elem%x(2, ii)
       pts(ii)%tag = ii
    end do

    ! then dimension all edges
    cons_per_edg = elem%npedg - 1
    allocate(cons(elem%nedgs * cons_per_edg))
    allocate(pts_list(elem%npedg))
    indx = 1
    last_pt = elem%nedgs + 1

    ! proceeding to fill connectors array
    do ii = 1, size(elem%edgs) ! loop over all edges
       pt1 = ii
       pt2 = ii + 1
       if ( ii .eq. size(elem%edgs) ) pt2 = 1 !loop closed

       if ( elem%p >= 2 ) then !higher order elements
          int1 = last_pt
          int2 = last_pt + elem%npedg - 3
          pts_list = (/ pt1, (/ (jj, jj = int1, int2) /) , pt2 /)
          last_pt = int2 + 1
       else
          pts_list = (/ pt1, pt2 /)
       end if

       call add_to_cons(pts, pts_list, indx, cons)

    end do

    ! setting up the holes which there is no hole
    ! here in this case
    holes(1)%tag = -1
    holes(1)%x = 1.0d20 ! some point outside the domain
    holes(1)%y = 1.0d20 ! some point outside the domain

    !
    call meshit('pnjYY', pts, cons, holes, report_before_in &
         , report_after_in, grd_out)

    ! clean ups
    deallocate(pts_list, pts, cons)

    ! done here
  end subroutine mesh_a_dg_elem

  !
  ! visualizes the curved element DG grid
  !
  !
  subroutine vis_curved_grid_dg(wspace, outfile)
    implicit none
    type(dg_wspace), intent(in), target :: wspace
    character(len = *), intent(in) :: outfile

    ! local vars
    integer :: ii, neqs, jj
    type(grid) :: grd2
    real*8, dimension(:,:), allocatable :: utmp
    type(element_dg2d), pointer :: telem => null()

    ! init
    neqs = wspace%elems(1)%neqs

    do ii = 1, size(wspace%elems) ! loop over all elems and export each

       telem => wspace%elems(ii)

       ! produce sub-triangle for a single curved SEM DG element
       ! just to use in the visuallization procedure
       call mesh_a_dg_elem(telem, grd2)

       allocate(utmp(neqs, grd2%nnodesg))

       ! project elem%U on sub-element triangulation
       do jj = 1, grd2%nnodesg
          ! simple one to one
          utmp(:, jj) = telem%U(:, jj)
       end do

       ! append the sub-element triangulation as a 
       ! new zone to already generated zonal grid 
       if ( ii .eq. 1 ) then
          call write_u_tecplot(outfile, grd2, utmp)
       else
          call write_u_tecplot(outfile, grd2, utmp, .true.)
       end if

       ! little clean-up
       deallocate(utmp)

    end do ! next DG element

    ! done here
  end subroutine vis_curved_grid_dg

  ! meshes a single high-order curvelinear 
  ! Discontinious Galerkin element using
  ! Renka triangulation (independent of bn edges)
  ! and writes the resulting grid in "grd_out" 
  subroutine mesh_a_dg_elem_renka(elem, grd_out)
    implicit none
    class(element_dg2d), intent(in) :: elem
    type(grid), intent(out) :: grd_out


    ! 
    call rtrimesh(xin = elem%x, icon = grd_out%icon &
         , xout = grd_out%x, yout = grd_out%y)

    ! lets put everything back into the grid object
    grd_out%nnodesg = size(grd_out%x)
    grd_out%ncellsg = size(grd_out%icon, 1)
    grd_out%ntri = grd_out%ncellsg
    grd_out%nquad4 = 0

    ! done here
  end subroutine mesh_a_dg_elem_renka

  !
  ! visualizes the curved element DG grid
  ! using a more generalized algorithm independent
  ! of the boundary limitation of the element.
  !
  subroutine vis_curved_grid_dg_renka(wspace, outfile)
    implicit none
    class(dg_wspace), intent(in), target :: wspace
    character(len = *), intent(in) :: outfile

    ! local vars
    integer :: ii, neqs, jj
    type(grid) :: grd2
    real*8, dimension(:,:), allocatable :: utmp
    class(element_dg2d), pointer :: telem => null()
    real*8 :: r, s

    ! init
    neqs = wspace%elems(1)%neqs

    do ii = 1, size(wspace%elems) ! loop over all elems and export each

       telem => wspace%elems(ii)

       ! produce sub-triangle for a single curved SEM DG element
       ! just to use in the visuallization procedure
       call mesh_a_dg_elem_renka(elem = telem, grd_out = grd2)

       allocate(utmp(neqs, grd2%nnodesg))

       ! project elem%U on sub-element triangulation
       do jj = 1, grd2%nnodesg
          ! find (r,s) of this (x,y) of subtriangulation
          call telem%xy2rs(x = grd2%x(jj), y = grd2%y(jj) &
               , maxitr = 40, tolrs = 1.0d-6, r = r, s = s) 

          ! project elem%u at that point
          call telem%comp_u(r, s, utmp(:, jj))

       end do

       ! append the sub-element triangulation as a 
       ! new zone to already generated zonal grid 
       if ( ii .eq. 1 ) then
          call write_u_tecplot(outfile, grd2, utmp)
       else
          call write_u_tecplot(outfile, grd2, utmp, .true.)
       end if

       ! little clean-up
       deallocate(utmp)

    end do ! next DG element

    ! done here
  end subroutine vis_curved_grid_dg_renka

  ! meshes a single high-order curvelinear 
  ! Discontinious Galerkin element
  ! and writes the resulting grid in "grd_out"
  !
  ! NOTE : this is version 2; more generalized than
  ! previous version "mesh_a_dg_elem"  
  !
  subroutine mesh_a_dg_elem2(elem, grd_out)
    implicit none
    type(element_dg2d), intent(in), target :: elem
    type(grid), intent(out) :: grd_out

    ! local vars
    integer :: ii, jj, pt1, pt2, indx
    type(point), dimension(:), allocatable :: pts
    type(connector), dimension(:), allocatable :: cons 
    type(point) :: holes(1)
    integer :: report_before_in, report_after_in

    class(edg_dg), pointer :: tedg => null()
    class(neigh_dg), pointer :: tneigh => null()
    integer :: npt, total_points
    integer, dimension(2) :: pts_list

    ! init.
    report_before_in = 0 
    report_after_in = 14

    ! first count the points on the edges and sub-edges (neighbors)
    ! this forms the skeleton (boundary) of high-order element
    ! 
    npt = 0
    do ii = 1, size(elem%edgs) ! loop over all edges
       tedg => elem%edgs(ii)
       do jj = 1, size(tedg%neighs)
          ! tneigh => tedg%neighs(jj)
          npt = npt + 1
       end do
    end do

    ! allocate space for all points
    total_points = npt + size(elem%x, 2)
    allocate(pts(total_points)) 

    ! now fill the x/y of thoses skeleton points
    npt = 1
    do ii = 1, size(elem%edgs) ! loop over all edges
       tedg => elem%edgs(ii)
       do jj = 1, size(tedg%neighs)
          tneigh => tedg%neighs(jj)
          ! add the point of the start of this segment
          pts(npt)%x = tneigh%xs(1)
          pts(npt)%y = tneigh%xs(2)
          pts(npt)%tag = npt
          ! next point
          npt = npt + 1
       end do
    end do

    ! then dimension all edges
    npt = npt - 1
    allocate(cons(npt))
    indx = 1

    ! proceeding to fill connectors array
    do ii = 1, npt ! loop over all edges
       pt1 = ii
       pt2 = ii + 1
       if ( ii .eq. npt ) pt2 = 1 !loop closed

       pts_list = (/ pt1, pt2 /)

       call add_to_cons(pts, pts_list, indx, cons)

    end do

    ! then add interpolation points to the end
    ! of the points array
    indx = 1
    do ii = (npt+1), total_points
       pts(ii)%x = elem%x(1, indx)
       pts(ii)%y = elem%x(2, indx)
       pts(ii)%tag = ii
       indx = indx + 1
    end do

    ! setting up the holes which there is no hole
    ! here in this case
    holes(1)%tag = -1
    holes(1)%x = 1.0d20 ! some point outside the domain
    holes(1)%y = 1.0d20 ! some point outside the domain

    !
    call meshit_visual('j', pts, cons, holes, report_before_in &
         , report_after_in, grd_out)

    ! clean ups
    deallocate(pts, cons)

    ! done here
  end subroutine mesh_a_dg_elem2

  !
  ! visualizes the curved element DG grid
  ! using a more generalized algorithm independent
  ! of the boundary limitation of the element.
  ! NOTE : this is still based on Triangle package!
  !
  subroutine vis_curved_grid_dg2(wspace, outfile)
    implicit none
    class(dg_wspace), intent(in), target :: wspace
    character(len = *), intent(in) :: outfile

    ! local vars
    integer :: ii, neqs, jj
    type(grid) :: grd2
    real*8, dimension(:,:), allocatable :: utmp
    class(element_dg2d), pointer :: telem => null()
    real*8 :: r, s

    ! init
    neqs = wspace%elems(1)%neqs

    do ii = 1, size(wspace%elems) ! loop over all elems and export each

       telem => wspace%elems(ii)

       ! produce sub-triangle for a single curved SEM DG element
       ! just to use in the visuallization procedure
       call mesh_a_dg_elem2(elem = telem, grd_out = grd2)
       ! call mesh_a_dg_elem_renka(elem = telem, grd_out = grd2)

       allocate(utmp(neqs, grd2%nnodesg))

       ! project elem%U on sub-element triangulation
       do jj = 1, grd2%nnodesg
          ! find (r,s) of this (x,y) of subtriangulation
          call telem%xy2rs(x = grd2%x(jj), y = grd2%y(jj) &
               , maxitr = 40, tolrs = 1.0d-6, r = r, s = s) 

          ! project elem%u at that point
          call telem%comp_u(r, s, utmp(:, jj))

       end do

       ! append the sub-element triangulation as a 
       ! new zone to already generated zonal grid 
       if ( ii .eq. 1 ) then
          call write_u_tecplot(outfile, grd2, utmp)
       else
          call write_u_tecplot(outfile, grd2, utmp, .true.)
       end if

       ! little clean-up
       deallocate(utmp)

    end do ! next DG element

    ! done here
  end subroutine vis_curved_grid_dg2

end module plot_curved_elems_dg
