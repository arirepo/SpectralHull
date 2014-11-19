module bcs
  ! implements boundary conditions
  use grid_opt
  use element_opt
  use bn_integral
  implicit none

  private

  integer, parameter :: nneqs = 1 ! num. of eqs.

  interface

     subroutine ff(x,y,fout)
       implicit none
       real*8, dimension(:), intent(out) :: fout
       real*8, intent(in) :: x,y
     end subroutine ff

  end interface

  ! a generic pointer-based function 
  ! association for imposing BCs
  ! so each BC can be defined as 
  ! an analytic function
  type func_pointer
     procedure(ff), pointer, nopass :: f => null()
     integer, dimension(nneqs) :: bc_type ! 1 = Dirichlet
                                          ! 2 = Numann
  end type func_pointer

  ! contains pointers to BC functions 
  type(func_pointer), dimension(5), save :: f_bc

  public :: impose_bcs, comp_bn_length, fill_Q, print_cp
 
contains 

  ! subroutine defined for boundary tag 1
  subroutine f1(x,y,fout)
    implicit none
    real*8, dimension(:),intent(out) :: fout
    real*8, intent(in) :: x,y ! global coords
    ! f1 = 50.0d0 !0.2d0 * (1.0d0 + x * (x-1.0d0))
    ! if ( (x .eq. 0.0d0) .or. (x .eq. 1.0d0) ) then
    !    f1 = 0.2d0 
    ! else
    !    f1 = 0.15d0
    ! end if
    !f1 = 0.0d0
    ! f1 = 50.0d0
    ! f1 = cos(2.0d0 * x) * sinh(2.0d0 * y)
    ! f1 = 323.15D0 
    ! fout = (/ 0.0d0, 0.0d0 /)
    ! fout = (/ 1.0d0/8.0d0 * (x**2.0d0 + y**2.0d0) /)
    ! fout = (/ cos(2.0d0 * x) * sinh(2.0d0 * y) /)
    fout = (/ 150.0d0 /)
  end subroutine f1

  ! subroutine defined for boundary tag 2
  subroutine f2(x,y,fout)
    implicit none
    real*8, dimension(:), intent(out) :: fout
    real*8, intent(in) :: x,y ! global coords
    ! f2 = 50.0d0 !0.2d0 !* y
    ! f2 = 0.2d0
    ! f2 = 50.0d0
    ! f2 = 0.0d0
    ! f2 = cos(2.0d0 * x) * sinh(2.0d0 * y)
    ! f2 = 323.15D0
    ! fout = (/ 0.0d0, 0.0d0 /)
    ! fout = (/ cos(2.0d0 * x) * sinh(2.0d0 * y) /)
    fout = (/ 150.0d0 /)
    ! fout = (/ -1.0d0 * x /)
  end subroutine f2

  ! subroutine defined for boundary tag 3
  subroutine f3(x,y,fout)
    implicit none
    real*8, dimension(:), intent(out) :: fout
    real*8, intent(in) :: x,y ! global coords
    !f3 = 50.0d0 * (1.0d0 + x * (x-1.0d0)) !150.0d0 !0.2d0 !* y
    ! f3 = 1.0d0 * y
    ! f3 = 0.2d0
    ! if ( (abs(x - 0.0d0) < 10.0D-9) .or. (abs(x - 2.0d0) < 10.0D-9) ) then
    !    f3 = 323.15D0
    ! else
    !    f3 = 423.15D0
    ! end if
    ! f3 = 0.0d0
    ! f3 = cos(2.0d0 * x) * sinh(2.0d0 * y)
    ! if ( (abs(x - 24.0d0) < 10.0D-9) .and. (abs(y - 0.0d0) < 10.0D-9) ) then
    !    fout = (/ 0.0d0, 0.0d0 /)
    ! else
    !    fout = (/ 5000.0d0, 0.0d0 /)
    ! end if

    ! fout = (/ 5000.0d0, 0.0d0 /)
    ! fout = (/ cos(2.0d0 * x) * sinh(2.0d0 * y) /)
    ! fout = (/ 3.0d0 /)
    fout = (/ 150.0d0 /)
  end subroutine f3

  ! subroutine defined for boundary tag 4
  subroutine f4(x,y,fout)
    implicit none
    real*8, dimension(:), intent(out) :: fout
    real*8, intent(in) :: x,y ! global coords
    !f4 = 50.0d0 ! 0.2d0 !* y
    ! f4 = 1.0d0 * y
    ! f4 = 0.2d0
    ! f4 = 50.0d0

    ! f4 = 1.0d0 
    ! f4 = cos(2.0d0 * x) * sinh(2.0d0 * y)
    ! f4 = 323.15D0
    ! fout = (/ 0.0d0, 0.0d0 /)
    ! fout = (/ cos(2.0d0 * x) * sinh(2.0d0 * y) /)
    ! fout = (/ 0.0d0 /)
    fout = (/ 150.0d0 /)
  end subroutine f4

  ! subroutine defined for boundary tag 5
  subroutine f5(x,y,fout)
    implicit none
    real*8, dimension(:), intent(out) :: fout
    real*8, intent(in) :: x,y ! global coords

    ! fout = (/ 0.0d0, 0.0d0 /)
    fout = (/ cos(2.0d0 * x) * sinh(2.0d0 * y) /)
    ! fout = (/ 0.0d0 /)
  end subroutine f5

  !      init BCs function array
  ! this subroutine initializes the array of
  ! data type "func_pointer", i.e., f_bc, where each entry 
  ! points to the corresponding analytical definition of
  ! the boundary condition and bc type, say:
  !
  ! f_bc(tag)%f => corresponfing f(x,y) for that boundary curve
  ! f_bc(tag)%bc_type = 1 or 2.
  subroutine init_bcs_func_array()
    implicit none

    f_bc(1)%f => f1
    ! f_bc(1)%bc_type = (/ 2, 2 /) 
    f_bc(1)%bc_type = (/ 1 /) 
    f_bc(2)%f => f2
    ! f_bc(2)%bc_type = (/ 2, 1 /)
    f_bc(2)%bc_type = (/ 1 /)
    f_bc(3)%f => f3
    ! f_bc(3)%bc_type = (/ 2, 2 /)
    f_bc(3)%bc_type = (/ 1 /)
    f_bc(4)%f => f4
    ! f_bc(4)%bc_type = (/ 2, 2 /)
    f_bc(4)%bc_type = (/ 1 /)
    f_bc(5)%f => f5
    ! f_bc(5)%bc_type = (/ 1, 2 /)
    f_bc(5)%bc_type = (/ 2 /)

    ! done here
  end subroutine init_bcs_func_array

  ! this subroutine imposes the boundary condition 
  ! in a general sense. The idea is to use analytical 
  ! functions for BC (either nuimann or dirichlet)
  ! and when higher order elements are used this approach
  ! is more staright and accurate because the value 
  ! of given BCs (either natural or essential) are
  ! always exact.
  
  subroutine impose_bcs(grd, KK, rhs)
    implicit none
    type(grid), target, intent(in) :: grd
    real*8, dimension(:,:,:,:), intent(inout) :: KK
    real*8, dimension(:,:), intent(inout) :: rhs

    ! local vars
    integer :: i, j, k, pt, ielem, iedg
    integer :: neqs, id
    real*8 :: x, y
    ! real*8 :: x1, x2, y1, y2, ll
    integer :: pt1, pt2
    integer :: ieq
    real*8, dimension(nneqs) :: tmp
    integer, dimension(:), pointer :: inter_pts => null(), dup => null()

    ! 
    neqs = size(KK,1)

    ! check fatal erros first
    if ( neqs .ne. nneqs ) then
       print *, 'the number of eqs in the system (i.e. neqs) does not match with' &
            , 'the number of eqs defined in this bcs module (i.e. nneqs)! stop'
       stop
    end if

    ! first initialize the boundary condition function array
    if( .not. associated(f_bc(1)%f) ) then
       call init_bcs_func_array()
    end if

    ! loop over ALL boundary edges
    do i =1, grd%nbedgeg 

       ! get id of that edge
       id = grd%ibedgeBC(i)
       select case (id)

       ! case(1:7)
       case(1:2)
          id = 1
          cycle
       ! case(8:9)
       case(3:4)

          id = 2
       case default 
          print *, 'wrong id! stop'
          stop

       end select
       ! find the elem containing that edge
       ! and also local edge number in that elem
       ielem = grd%ibedgeELEM(i)
       iedg = grd%ibedgeELEM_local_edg(i)

       ! ! should be between 1 and the number of fi defined here
       ! if ( (id < 1) .or. (id > size(f_bc)) ) then
       !    ! id = 3
       !    print *, 'boundary tag is outside of the range of defined' &
       !         , ' bcs. Please fix the number in CAD or add new' &
       !         , ' function definitions' &
       !         , ' to "bcs" module. stop.'  
       !    stop
       ! end if

       ! we are good to go now ...
       do ieq = 1, neqs
          select case (f_bc(id)%bc_type(ieq)) ! decide on BC type

          case(1) ! general dirichlet BCs on that edge

             ! loop over 2 start and end nodes on that edge
             do k = 1, 2 
                pt = grd%ibedge(i,k) ! global number of that node
                x = grd%x(pt)
                y = grd%y(pt)
                KK(ieq,:,pt,:) = 0.0d0 ! freeze that row 
                do j = ieq, ieq
                   KK(j,j,pt,pt) = 1.0d0 ! ones on the diagonal
                end do
                call f_bc(id)%f(x,y,tmp)
                rhs(ieq,pt) = tmp(ieq)
             end do ! k - loop

             ! start looping over nodes inside that edge
             if ( allocated(grd%el2edg(ielem)%edg1) ) then
                ! associate a universal pointer to points on that edge.
                select case (iedg)
                case (1)
                   inter_pts => grd%el2edg(ielem)%edg1
                case (2)
                   inter_pts => grd%el2edg(ielem)%edg2
                case (3)
                   inter_pts => grd%el2edg(ielem)%edg3
                case default
                   print *, 'something is wring in applying BCs' &
                        , ' on points inside a boundary edge! stop'
                   stop
                end select
                ! if at least one point exists on that edge
                if ( size(inter_pts) >= 1 ) then
                   ! loop over all interior pts on that edge
                   do k = 1, size(inter_pts) 
                      pt = inter_pts(k) ! global number of that node
                      x = grd%x(pt)
                      y = grd%y(pt)
                      KK(ieq,:,pt,:) = 0.0d0 ! freeze that row 
                      do j = ieq, ieq
                         KK(j,j,pt,pt) = 1.0d0 ! ones on the diagonal
                      end do
                      call f_bc(id)%f(x,y,tmp)
                      rhs(ieq,pt) = tmp(ieq)
                   end do ! k - loop
                end if

             end if ! there exists bn pts on that edge


          case (2) ! general nuimann (not implemented yet!)

             ! !                 NNOOTTEE
             ! !         CAN CAUSE FATAL PROBLEMS 
             ! ! THE FOLLOWING IS JUST FOR LINEAR ELEMENTS
             ! if ( size(grd%ibedge,2) > 2 ) then
             !    print * , 'The current implementation of nuimann BCs is only' &
             !         , ' valid for linear elements! higher order boundary' &
             !         , ' element detected! stop.'
             !    stop
             ! end if

             ! pt1 = grd%ibedge(i,1) ! global number of node1
             ! pt2 = grd%ibedge(i,2) ! global number of node2
             ! x1 = grd%x(pt1); x2 = grd%x(pt2)
             ! y1 = grd%y(pt1); y2 = grd%y(pt2)
             ! ll = sqrt((x2 - x1)**2.0d0 + (y2 - y1)**2.0d0)
             ! call f_bc(id)%f(x1,y1,tmp)
             ! rhs(ieq,pt1) = rhs(ieq,pt1) + tmp(ieq) * ll / 2.0d0
             ! rhs(ieq,pt2) = rhs(ieq,pt2) + tmp(ieq) * ll / 2.0d0

          case default ! bc type vague

             print *,'unknown type of boundary condition. stop!'
             stop

          end select
       end do ! loop over equations
    end do ! done all boundary edges

    ! handling duplicate nodes (if any)
    if ( allocated(grd%dup_nodes) .and. ( size(grd%dup_nodes) >= 1) ) then
       dup => grd%dup_nodes
       do k = 1, grd%tot_repeated_bn_nodes
          pt1 = dup(2 * k - 1)
          pt2 = dup(  2 * k  )
          if(maxval(abs(KK(:,:,pt1,:))) .eq. 0.0d0) call swap(pt1, pt2)
          KK(:,:,pt2,:) = 0.0d0 ! freeze that row 
          rhs(:,pt2) = 0.0d0
          do j = 1, neqs
             KK(j,j,pt2,pt1) = 1.0d0
             KK(j,j,pt2,pt2) = -1.0d0
          end do
       end do
    end if

    ! done here                                         
  end subroutine impose_bcs
  
  subroutine swap(a, b)
    implicit none
    integer, intent(inout) :: a, b

    ! local vars
    integer :: tmp

    tmp = a
    a = b
    b = tmp

    ! done here
  end subroutine swap

  subroutine comp_bn_length(grd, elems)
    implicit none
    type(grid), intent(in) :: grd
    type(element), dimension(:), intent(inout) :: elems

    ! local vars
    integer :: i, loc_edg, ielem, pt1, pt2
    real*8 :: tol
    real*8, dimension(1) :: val , tmp

    tol = 1.0d-14
    val = 0.0d0; tmp = 0.0d0

    do i = 1, grd%nbedgeg
       pt1 = grd%ibedge(i,1)
       pt2 = grd%ibedge(i,2)
       ielem = grd%ibedgeELEM(i)
       loc_edg = grd%ibedgeELEM_local_edg(i)

       call comp_elem_bn_integral(grd, elems(ielem), ielem, loc_edg, tol, func, tmp)
       print *, 'tmp = ', tmp, '(x1,y1) = ', grd%x(pt1), grd%y(pt1) &
            , '(x2,y2) = ', grd%x(pt2), grd%y(pt2) 
       val = val + tmp

    end do

    print *, 'val = ', val

    ! done here
  end subroutine comp_bn_length

  ! ! custom subroutine to compute the integrand over the edge
  ! subroutine func(elem, x0, y0, x, y, val)
  !   implicit none
  !   type(element), intent(inout) :: elem
  !   real*8, intent(in) :: x0, y0, x, y
  !   real*8, dimension(:), intent(out) :: val

  !   ! val = 1.0d0
  !   val = x**2.0d0
  !   print *, 'func = ', val

  !   ! done here
  ! end subroutine func

  ! custom subroutine to compute the integrand over the edge
  subroutine func(elem, x0, y0, x, y, val)
    implicit none
    type(element), intent(inout) :: elem
    real*8, intent(in) :: x0, y0, x, y
    real*8, dimension(:), intent(out) :: val


    call elem%tbasis%eval(x0, y0, 0,  val)

    print *, 'size(val) inside = ', size(val)
    ! stop

    ! done here
  end subroutine func
  
  subroutine fill_Q(grd, elem, ielem)
    implicit none
    type(grid), intent(in) :: grd
    type(element), intent(inout) :: elem
    integer, intent(in) :: ielem

    ! local vars
    integer :: loc_edg
    real *8 :: tol 
    real*8, dimension(grd%npe(ielem)) :: tmp

    tol = 1.0d-10
    tmp = 0.0d0
    elem%Q = 0.0d0

    do loc_edg = 1, 3
print *, 'size(tmp) = ', size(tmp)
! stop
       call comp_elem_bn_integral(grd, elem, ielem, loc_edg, tol, func, tmp)
       elem%Q(1,:) = elem%Q(1,:) + tmp
    end do

    ! done here
  end subroutine fill_Q

  !
  subroutine print_cp(u, grd, elems)
    implicit none
    real*8, dimension(:,:), intent(in) :: u
    type(grid), target, intent(in) :: grd
    type(element), dimension(:), target, intent(in) :: elems

    ! local vars
    integer :: i, k, id, ielem, iedg, pt
    integer, dimension(:), pointer :: inter_pts => null()
    real*8, dimension(1) :: dudx, dudy
    real*8 :: cp, r, s, x, y
    type(element), pointer :: elem => null()

    ! loop over ALL boundary edges
    do i =1, grd%nbedgeg 

       ! get id of that edge
       id = grd%ibedgeBC(i)
       select case (id)

       case(1:2)
       ! case(1:7)
          ! cycle
       case default 
cycle
          ! print *, 'wrong id! stop'
          ! stop

       end select
       ! find the elem containing that edge
       ! and also local edge number in that elem
       ielem = grd%ibedgeELEM(i)
       iedg = grd%ibedgeELEM_local_edg(i)
       elem => elems(ielem)

       ! start looping over nodes inside that edge
       if ( allocated(grd%el2edg(ielem)%edg1) ) then
          ! associate a universal pointer to points on that edge.
          select case (iedg)
          case (1)
             inter_pts => grd%el2edg(ielem)%edg1
          case (2)
             inter_pts => grd%el2edg(ielem)%edg2
          case (3)
             inter_pts => grd%el2edg(ielem)%edg3
          case default
             print *, 'something is wring in applying BCs' &
                  , ' on points inside a boundary edge! stop'
             stop
          end select
          ! if at least one point exists on that edge
          if ( size(inter_pts) >= 1 ) then
             ! loop over all interior pts on that edge
             ! and print cp
             do k = 1, size(inter_pts) 
                pt = inter_pts(k) ! global number of that node
                x = grd%x(pt)
                y = grd%y(pt)
                call xy2rs(x, y, elem, ielem, grd, 40, 1.0d-12, r, s)
                call comp_grad_u_point(u, grd, elem, ielem, &
                                       r, s, dudx, dudy)
                cp = 1.0d0 - (dudx(1)**2.0d0 + dudy(1)**2.0d0)
                print *, '+1221360', x, y, cp

             end do ! k - loop
          else
             print *, ' fatal : size(inter_pts) is' &
                  , ' not >= 1 for CP computation! stop'
             stop
          end if
       else
          print *, ' fatal : at least one point on boundary' &
               , '  edge (or equivalently p2) is required' &
               , ' to compute Cp. stop'
          stop
       end if ! there exists bn pts on that edge

    end do ! done all boundary edges

    ! done here
  end subroutine print_cp

end module bcs
