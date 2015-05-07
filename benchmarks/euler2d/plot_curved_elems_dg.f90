module plot_curved_elems_dg
  use grid_opt
  use ps
  use element_opt_dg2d
  use trimesher
  use dg_workspace
  implicit none

  private


  public :: vis_curved_grid_dg

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

end module plot_curved_elems_dg
