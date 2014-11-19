module mfree
  use grid_opt
  use fem_reordering
  use element_opt
  implicit none

  private

  type mf_struct
     type(grid), pointer :: grd => null()
     type(element), dimension(:), pointer :: elems => null()
     integer, dimension(:), allocatable :: tags
     real*8, dimension(:), allocatable :: values
     logical :: ordered, unique
     integer, dimension(:), allocatable :: fixed, free


   contains
     procedure :: init => init_mf_struct
     procedure :: Ax_gen => A_x
     procedure :: Ax => full_A_x

  end type mf_struct

  public :: mf_struct
  public :: init_cond
  public :: gmres_orig
  public :: norm

contains

  subroutine comp_fixed_free(grd, neqs, dirch_tags, mf)
    implicit none
    type(grid), intent(in), target :: grd
    integer, intent(in) :: neqs
    integer, dimension(:), intent(in) :: dirch_tags
    type(mf_struct), intent(inout) :: mf

    ! local vars
    integer :: i, id, k, pt, nfixed, nfree
    integer :: ielem, iedg
    integer, dimension(:), allocatable :: itmp, itmp2
    integer, dimension(:), pointer :: inter_pts => null()

    ! first compute and store the node number of fixed bcs 
    nfixed = 0
    do i = 1, grd%nbedgeg ! loop over boundary edges 
       id = grd%ibedgeBC(i)
       ! skip if this edge's tag is not a dirichlet tag
       if ( .not. any( dirch_tags .eq. id ) ) cycle

       ! add 2 start and end nodes on that edge
       do k = 1, 2 
          pt = grd%ibedge(i,k) ! global number of that node
          call add_nodes(neqs, nfixed, itmp, pt, mf%fixed)
       end do ! k - loop

       ! find the element containing that edge and its locallity
       ielem = grd%ibedgeELEM(i)
       iedg = grd%ibedgeELEM_local_edg(i)

       ! associate a universal pointer to the points on that edge.
       select case (iedg)
       case (1)
          inter_pts => grd%el2edg(ielem)%edg1
       case (2)
          inter_pts => grd%el2edg(ielem)%edg2
       case (3)
          inter_pts => grd%el2edg(ielem)%edg3
       case default
          print *, 'something is wrong in computing fixed' &
               , ' points between two end points of a boundary edge. stop'
          stop
       end select

       ! if at least one point exists on that edge
       if ( size(inter_pts) >= 1 ) then

          ! loop over all interior pts on that edge
          do k = 1, size(inter_pts) 
             pt = inter_pts(k) ! global number of that node
             call add_nodes(neqs, nfixed, itmp, pt, mf%fixed)
          end do ! k - loop

       end if



    end do ! done adding fixed nodes

    ! compute free nodes
    nfree = 0
    do i = 1, grd%nnodesg
       if ( any(i .eq. mf%fixed) ) cycle
       call add_nodes(neqs, nfree, itmp, i, mf%free)
    end do

    ! make lists unique if requested
    if ( mf%unique ) then
       call unique(mf%fixed, itmp)
       mf%fixed = itmp

       call unique(mf%free, itmp)
       mf%free = itmp
    end if

    ! reorder the list if requested
    if ( mf%ordered ) then

       ! reorder mf%fixed
       if (allocated(itmp)) deallocate(itmp)
       if (allocated(itmp2)) deallocate(itmp2)
       allocate(itmp(size(mf%fixed)))
       allocate(itmp2(size(mf%fixed)))
       call arash_sort(dble(mf%fixed), itmp, '+')
       itmp2 = mf%fixed(itmp)
       mf%fixed = itmp2

       ! reorder mf%free
       if (allocated(itmp)) deallocate(itmp)
       if (allocated(itmp2)) deallocate(itmp2)
       allocate(itmp(size(mf%free)))
       allocate(itmp2(size(mf%free)))
       call arash_sort(dble(mf%free), itmp, '+')
       itmp2 = mf%free(itmp)
       mf%free = itmp2

    end if

    ! done here
  end subroutine comp_fixed_free

  ! adds sequntial node indices after given 
  ! node number, ie. "pt" to the end of out.
  ! the number of added nodes are 1..neqs.
  subroutine add_nodes(neqs, n, itmp, pt, out)
    implicit none
    integer, intent(in) :: neqs, pt
    integer, intent(inout) :: n
    integer, dimension(:), allocatable :: itmp, out

    ! local vars
    integer :: j

    do j = 1, neqs
       n = n + 1
       allocate( itmp(n) )
       if ( n > 1 ) itmp(1:(n-1)) = out(1:(n-1))
       itmp(n) = pt + (j-1)
       call move_alloc(itmp, out)
    end do

    ! done here
  end subroutine add_nodes

  ! converts the possibly non-unique array "a"
  ! to an array with unique integer elements.
  subroutine unique(a, b)
    implicit none
    integer, dimension(:), intent(in) :: a
    integer, dimension(:), allocatable :: b

    ! local vars
    integer :: i, j, val, count, n
    integer, dimension(:), allocatable :: itmp

    if ( allocated(b) ) deallocate(b)

    n = 0
    do i = 1, size(a)

       val = a(i)
       count = 0
       do j = i, size(a)
          if ( val .eq. a(j) ) count = count + 1
       end do
       if (count > 1) cycle ! do not add point

       n = n + 1
       allocate(itmp(n))
       if (n > 1) itmp(1:(n-1)) = b(1:(n-1))
       itmp(n) = val

       call move_alloc(itmp, b)

    end do

    ! done here
  end subroutine unique

  ! NOTE : USE THIS IN NORMAL SITUATIONS
  ! this is the fast and scalable version of
  ! previous subroutine "comp_fixed_free"
  !
  subroutine comp_fixed_free_fast(grd, neqs, dirch_tags, mf)
    implicit none
    type(grid), intent(in), target :: grd
    integer, intent(in) :: neqs
    integer, dimension(:), intent(in) :: dirch_tags
    type(mf_struct), intent(inout) :: mf

    ! local vars
    integer :: i, id, k, pt, nfixed, nfree
    integer :: ielem, iedg
    integer, dimension(:), allocatable :: itmp
    integer, dimension(:), pointer :: inter_pts => null()
    logical, dimension(grd%nnodesg) :: nailed

    nailed = .false. ! all nodes are free except specified 
    ! first compute the node number of fixed bcs 
    do i = 1, grd%nbedgeg ! loop over boundary edges 
       id = grd%ibedgeBC(i)
       ! skip if this edge's tag is not a dirichlet tag
       if ( .not. any( dirch_tags .eq. id ) ) cycle

       ! add 2 start and end nodes on that edge
       do k = 1, 2 
          pt = grd%ibedge(i,k) ! global number of that node
          nailed(pt) = .true.
       end do ! k - loop

       ! find the element containing that edge and its locallity
       ielem = grd%ibedgeELEM(i)
       iedg = grd%ibedgeELEM_local_edg(i)

       ! associate a universal pointer to the points on that edge.
       select case (iedg)
       case (1)
          inter_pts => grd%el2edg(ielem)%edg1
       case (2)
          inter_pts => grd%el2edg(ielem)%edg2
       case (3)
          inter_pts => grd%el2edg(ielem)%edg3
       case default
          print *, 'something is wrong in computing fixed' &
               , ' points between two end points of a boundary edge. stop'
          stop
       end select

       ! if at least one point exists on that edge
       if ( size(inter_pts) >= 1 ) then

          ! loop over all interior pts on that edge
          do k = 1, size(inter_pts) 
             pt = inter_pts(k) ! global number of that node
             nailed(pt) = .true.
          end do ! k - loop

       end if

    end do ! done adding fixed nodes

    ! store fixed and free nodes
    nfixed = 0; nfree = 0
    do i = 1, grd%nnodesg
       if ( nailed(i) ) then
          call add_nodes(neqs, nfixed, itmp, i, mf%fixed)
       else
          call add_nodes(neqs, nfree, itmp, i, mf%free)
       end if
    end do

    ! done here
  end subroutine comp_fixed_free_fast

  subroutine Ae_x(grd, ielem, neqs, mf, A, x, xout)
    implicit none
    type(grid), intent(in) :: grd
    integer, intent(in) :: ielem, neqs
    type(mf_struct), intent(inout) :: mf
    real*8, dimension(:,:,:,:), intent(in) :: A
    real*8, dimension(:), intent(in) :: x
    real*8, dimension(:), intent(inout) :: xout

    ! local vars
    integer :: i, j, pt_i, pt_j
    integer :: i1, i2, j1, j2
    real*8, dimension(neqs) :: tmp

    do i = 1, size(A, 3)

       pt_i = grd%icon(ielem, i)
       if (allocated(mf%fixed) .and. (size(mf%fixed) .ge. 1)) then 
          if ( any ( pt_i .eq. mf%fixed ) ) then
             cycle
          end if
       end if
       tmp = 0.0d0

       do j = 1, size(A, 4)

          pt_j = grd%icon(ielem, j)
          j1 = pt_j
          j2 = pt_j + (neqs - 1) 
          tmp = tmp + matmul(A(:,:, i, j) , x(j1:j2))   

       end do

       i1 = pt_i
       i2 = pt_i + (neqs - 1) 
       xout(i1:i2) = xout(i1:i2) + tmp 

    end do

    ! done here
  end subroutine Ae_x

  subroutine A_x(mf, x, xout)
    implicit none
    class(mf_struct), intent(inout) :: mf
    real*8, dimension(:), intent(in) :: x
    real*8, dimension(:), intent(inout) :: xout

    ! local vars
    integer :: i

    do i = 1, size(mf%elems)
       call Ae_x(mf%grd, i, mf%elems(i)%neqs, mf, mf%elems(i)%A, x, xout)
    end do

    ! at the end, set the values of fix nodes (dirichlet)
    ! back to their original values
    if ( allocated(mf%fixed) ) then
       if (size(mf%fixed) .ge. 1 ) xout(mf%fixed) = x(mf%fixed)
    end if

    ! done here
  end subroutine A_x

  subroutine full_A_x(mf, x, xout)
    implicit none
    class(mf_struct), intent(inout) :: mf
    real*8, dimension(:), intent(in) :: x
    real*8, dimension(:), intent(out) :: xout

    ! local vars
    integer :: i

    ! HARD reset
    xout = 0.0d0

    ! ---------------------------------------------------------
    ! repeat the following part between dasht lines for different 
    !                     elemental matrices
    ! if other form of matrix vector product is desired.
    ! default is mass matrix but you can set it to any elemental 
    ! matrix like mass matrix or damping matrix if you want.
    ! then repeat matrix-vector product A_x. The result would be
    ! accumulated. Basically copy and paste and change the area between
    ! the dasht lines.
    !
    ! setting up the pointers to appropriate elemental matrices
    do i = 1, size(mf%elems)
       mf%elems(i)%A => mf%elems(i)%K 
    end do
    ! perform a matrix-vector product
    call A_x(mf, x, xout)
    ! ---------------------------------------------------------

    ! done here
  end subroutine full_A_x

  ! initialization based on given tag values
  subroutine init_cond(grd, tags, vals, method, x)
    implicit none
    type(grid), intent(in), target :: grd
    integer, dimension(:), intent(in) :: tags
    real*8, dimension(size(tags)), intent(in) :: vals
    character(len = *), intent(in) :: method
    real*8, dimension(:), intent(inout) :: x

    ! local vars
    integer :: i, j, k, loc, pt, ielem, iedg
    integer, dimension(:), pointer :: inter_pts => null()
    real*8 :: ave

    ! choose initialization method
    select case (method)
    case ('ave')
       ave = sum(vals) / dble(size(vals))
       x = ave
    case ('zero')
       x = 0.0d0
    case('none')
       ! do nothing 
    case default
       print *, 'initilization method is not specified correctly! stop.'
       stop
    end select

    ! proceeding to set boundary values based on 
    ! given boundary tags
    !
    do i = 1, grd%nbedgeg

       loc = minloc(array=tags, dim = 1, mask = tags .eq. grd%ibedgeBC(i))
       if (loc .eq. 0 ) cycle ! this bn tag is not in the given list

       ! fixing ICs at two end nodes of bedges
       do j = 1, 2
          pt = grd%ibedge(i, j)
          x(pt) = vals(loc) 
       end do

       ! find the element containing that edge and its locallity
       ielem = grd%ibedgeELEM(i)
       iedg = grd%ibedgeELEM_local_edg(i)

       ! associate a universal pointer to the points on that edge.
       select case (iedg)
       case (1)
          inter_pts => grd%el2edg(ielem)%edg1
       case (2)
          inter_pts => grd%el2edg(ielem)%edg2
       case (3)
          inter_pts => grd%el2edg(ielem)%edg3
       case default
          print *, 'something is wrong in computing fixed' &
               , ' points between two end points of a boundary edge. stop'
          stop
       end select

       ! if at least one point exists on that edge
       if ( size(inter_pts) >= 1 ) then

          ! loop over all interior pts on that edge
          do k = 1, size(inter_pts) 
             pt = inter_pts(k) ! global number of that node
             x(pt) = vals(loc) 
          end do ! k - loop

       end if

    end do

    ! done here
  end subroutine init_cond


  subroutine init_mf_struct(this, grd, elems, tags, values)
    implicit none
    class(mf_struct), intent(inout) :: this
    type(grid), intent(in), target :: grd
    type(element), dimension(:), intent(in), target :: elems
    integer, dimension(:), intent(in) :: tags
    real*8, dimension(:), intent(in) :: values


    this%grd => grd
    this%elems => elems

    ! init this%tags
    if ( allocated(this%tags) ) deallocate(this%tags)
    allocate(this%tags(size(tags)))
    this%tags = tags

    ! init this%values
    if ( allocated(this%values) ) deallocate(this%values)
    allocate(this%values(size(values)))
    this%values = values

    this%ordered = .true.
    this%unique = .true.

    ! determine fixed and free nodes
    call comp_fixed_free_fast(this%grd, this%elems(1)%neqs, this%tags, this)


    ! done here
  end subroutine init_mf_struct

  subroutine gmres_orig(mf, b, N, epsil, num, nrst, sol, res)
    implicit none
    ! inputs
    class(mf_struct), intent(inout) :: mf
    real*8, dimension(:), intent(in) :: b
    real*8, dimension(size(mf%free)), intent(in) :: N
    real*8, intent(in) :: epsil
    integer, intent(in) :: num, nrst

    !outputs
    real*8, dimension(:), intent(inout) :: sol
    real*8, dimension(:), allocatable :: res

    ! locals
    integer :: k, j, i
    integer :: finish_flag, jrst
    real*8, dimension(:,:), allocatable :: v, h
    real*8, dimension(:), allocatable :: r0, g, yk, uk, s, c, zk, alpha, tmp
    real*8 :: norm_tmp, delta, gamma, resi
    integer :: nrows, nvars, len_res


    ! resetting the flag
    finish_flag = 0
    resi = 0.0d0
    nrows = size(b)
    nvars = size(mf%free)
    len_res = 0

    ! we don't use dynamic mem allocation in the gmres inner-outer loops 
    ! to increase spped. this might not be memory efficient though.
    allocate(v(nvars,num+1), r0(nrows), g(num+1), s(num), c(num))
    allocate(yk(nrows), uk(nvars), zk(nvars) )
    allocate(h(num+1,num), alpha(num))

    ! do j number of restarts
    do jrst = 1, nrst    
       ! reset ------------
       v = 0.0d0; h = 0.0d0
       c = 0.0d0; s = 0.0d0
       zk = 0.0d0; g = 0.0d0
       alpha = 0.0d0
       yk = 0.0d0
       ! ------------------
       ! ***********************
       call mf%Ax(sol, r0)
       r0 = b - r0
       ! r0 = b - matmul(A,sol)
       ! ***********************
       v(:,1) = r0(mf%free)
       norm_tmp = norm(v(:,1))
       v(:,1) = v(:,1)/norm_tmp
       g(1) = norm_tmp
       ! Arnoldi iterations    
       do k = 1, num
          g(k+1) = 0.0d0
          ! ***********************
          yk = 0.0d0
          yk(mf%free) = v(:,k) / N
          call mf%Ax(yk, r0)
          uk = r0(mf%free)
          ! yk = v(:,k) / N
          ! uk = matmul(A, yk)
          ! ***********************
          do j = 1, k
             h(j,k) = dot_product(v(:,j),uk)
             !grahm-schmidt orthogonalization
             uk = uk - h(j,k) * v(:,j) 
          end do
          h(k+1,k) = norm(uk)        
          v(:,k+1) = uk/h(k+1,k)
          ! applying givens
          do j = 1, (k-1)
             delta = h(j,k)
             h(j,k) = c(j)*delta + s(j)*h(j+1,k)
             h(j+1,k) = -s(j)*delta + c(j) * h(j+1,k)
          end do
          gamma = sqrt(h(k,k)**(2.0d0) + h(k+1,k)**(2.0d0))
          c(k) = h(k,k) / gamma
          s(k) = h(k+1,k) / gamma
          h(k,k) = gamma
          h(k+1,k) = 0.0d0
          delta = g(k)
          g(k) = c(k) * delta + s(k) * g(k+1)
          g(k+1) = -s(k) * delta + c(k) * g(k+1)
          resi = abs(g(k+1))
          ! add to residual history
          len_res = len_res + 1
          allocate(tmp(len_res))
          if(len_res > 1) tmp(1:(len_res-1)) = res(1:(len_res-1)) 
          tmp(len_res) = resi
          call move_alloc(tmp, res)
          !
          if( resi <= epsil) then
             finish_flag = 1
             goto 100
          end if

       end do
       k = k - 1
       ! solving backward for alpha
100    alpha(k) = g(k) / h(k,k)
       do i = k-1, 1, -1
          alpha(i) = g(i) / h(i,i)
          do j = i+1, k      
             alpha(i) = alpha(i) - h(i,j)/h(i,i) * alpha(j)     
          end do
       end do
       ! compute the final directional vector using
       ! the combination of search vectors 
       do j = 1, k
          zk = zk + alpha(j)*v(:,j)
       end do
       ! updating solution
       ! ***********************
       zk = zk/N
       r0 = 0.0d0
       r0(mf%free) = zk
       sol = sol + r0     
       ! sol = sol + zk/N    
       ! ***********************
       if(finish_flag == 1) then 
          goto 200        
       end if
    end do
    ! subroutine final clean-ups
200 deallocate(v, r0, g, s, c)
    deallocate(yk, uk, zk, h, alpha)

    return !done <here>

  end subroutine gmres_orig

  function norm(x)
    implicit none
    real*8, dimension(:), intent(in) :: x
    real*8 :: norm

    ! locals
    integer :: i

    norm = 0.0d0 !reset
    do i = 1, size(x)
       norm = norm + x(i)* x(i)
    end do

    norm = sqrt(norm)

    ! done
  end function norm

end module mfree
