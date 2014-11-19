module time_integration
  use spem_2d
  use fem_utils
  use grid_opt
  implicit none

  private

  type temporal

     integer :: tag
     real*8 :: dt
     real*8, dimension(:,:), allocatable :: Keff, LU, BB
     real*8, dimension(:), allocatable :: Reff
     logical :: lu_started
     integer, dimension(:), allocatable :: IPIV

     ! only for Newmark method
     real*8 :: alpha, delta

  end type temporal

  public :: temporal, init_time, march_newmark

contains

  ! initializes vars related to time-dependent
  ! nature of the equations
  ! un, dun, d2un at t = 0
  subroutine init_time(ttemporal, tfem, u0, du0)
    implicit none
    type(temporal), intent(inout) :: ttemporal
    type(fem_struct), intent(inout) :: tfem
    real*8, dimension(:), intent(in) :: u0, du0

    ! local vars
    real*8, dimension(size(tfem%KK_mat, 1)) :: b

    ! init.
    tfem%un  = u0
    tfem%dun = du0

    b = -1.0d0 * matmul(tfem%KK_mat, u0)

    call linear_solve(tfem%MM_mat, b, tfem%d2un, tfem%param)

    ! 
    allocate( ttemporal%Keff(size(tfem%KK_mat, 1), size(tfem%KK_mat, 2)) )
    allocate( ttemporal%LU(size(tfem%KK_mat, 1), size(tfem%KK_mat, 2)) )
    allocate( ttemporal%BB(size(tfem%KK_mat, 1), 1) )
    allocate( ttemporal%Reff(size(tfem%KK_mat, 1)) )
    ttemporal%Keff = 0.0d0; ttemporal%LU = 0.0d0; ttemporal%BB = 0.0d0 
    ttemporal%Reff = 0.0d0
    ttemporal%lu_started = .false.
    allocate(ttemporal%IPIV(size(tfem%KK_mat, 1)))
    ttemporal%IPIV = 0


    ! done here
  end subroutine init_time

  ! implicit Newmark method
  subroutine march_newmark(tfem, ttemporal)
    implicit none
    type(fem_struct), intent(inout) :: tfem
    type(temporal), intent(inout) :: ttemporal

    ! local vars
    real*8 :: a0, a1, a2, a3, a4, a5, a6, a7
    real*8, dimension(size(tfem%KK_mat,1)) :: sum_u_1
    ! real*8, dimension(size(tfem%KK_mat,1)) :: sum_u_2
    integer :: INFO

    associate( M => tfem%MM_mat, K => tfem%KK_mat, Keff => ttemporal%Keff, &
         Reff => ttemporal%Reff, dt => ttemporal%dt, &
         alpha => ttemporal%alpha, delta => ttemporal%delta, &
         un => tfem%un, dun => tfem%dun, d2un => tfem%d2un, &
         unp1 => tfem%unp1, dunp1 => tfem%dunp1, d2unp1 => tfem%d2unp1, &
         LU => ttemporal%LU, IPIV => ttemporal%IPIV, BB => ttemporal%BB )

      ! comp rhs 
    a0 = 1.0d0 / (alpha * dt*dt)
    a1 = delta / (alpha * dt)
    a2 = 1.0d0 / (alpha * dt)
    a3 = 1.0d0 / (2.0d0 * alpha) - 1.0d0
    a4 = delta/ alpha - 1.0d0
    a5 = 0.5d0 * dt * (delta / alpha - 2.0d0)
    a6 = dt * (1.0d0 - delta)
    a7 = delta * dt
    sum_u_1 = a0 * un + a2 * dun + a3 * d2un
    ! sum_u_2 = a1 * un + a4 * dun + a5 * d2un
    Reff = matmul(M, sum_u_1)
    ! Reff = Rext + matmul(M, sum_u_1) + matmul(C, sum_u_2)

    ! comp lhs
    Keff = a0 * M + K
    ! Keff = a0 * M + a1 * C + K

    ! apply boundary conditions
    ! call sys%bcs()
    call freeze_bcs(tfem%grd, Keff, Reff, 1, (/1,2,3, 4,5,6,7/) )

    ! solve
    ! call linear_solve(Keff, Reff, unp1, tfem%param)

    if ( .not. ttemporal%lu_started ) then
       ! first time only do LU and store it
       print *, 'computing LU ...'
       LU = Keff
       call DGETRF( size(LU,1), size(LU,2), LU, size(LU,1), IPIV, INFO )
       if ( INFO .ne. 0 ) then 
          print *, 'something is wrong in LU factorization! stop'
          stop
       end if
       ttemporal%lu_started = .true.
    end if

    ! use LU to solve it
    BB(:,1) = Reff
    call DGETRS( 'No transpose', size(LU,1), 1, LU, size(LU,1) &
         , IPIV, BB, size(LU,1), INFO )
    if ( INFO .ne. 0 ) then 
       print *, 'something is wrong in LU-solve! stop'
       stop
    end if

    ! update
    unp1 = BB(:, 1)

    ! update
    d2unp1 = a0 * (unp1 - un) - a2 * dun - a3 * d2un
    dunp1 = dun + a6 * d2un + a7 * d2unp1

    end associate

    ! done here
  end subroutine march_newmark

  !
  subroutine freeze_bcs(grd, KK, rhs, neqs, tag_range)
    implicit none
    type(grid), target, intent(in) :: grd
    real*8, dimension(:,:), intent(inout) :: KK
    real*8, dimension(:), intent(inout) :: rhs
    integer, intent(in) :: neqs
    integer, dimension(:), intent(in) :: tag_range

    ! local vars
    integer :: i, j, k, pt, ielem, iedg
    integer :: id
    integer :: i1, i2
    real*8 :: ave
    integer, dimension(:), pointer :: inter_pts => null()
    ! integer, dimension(:), pointer :: dup => null()

    ! find the average of the main diagonal
    ave = 0.0d0
    do i = 1, size(KK, 1)
       ave = ave + abs(KK(i, i))
    end do
    ave = ave / dble(size(KK, 1))

    ! loop over ALL boundary edges
    do i =1, grd%nbedgeg 

       ! ! get id of that edge
       id = grd%ibedgeBC(i)
       ! decide if it should be freezed
       if (any( id .eq. tag_range) ) cycle

       ! find the elem containing that edge
       ! and also local edge number in that elem
       ielem = grd%ibedgeELEM(i)
       iedg = grd%ibedgeELEM_local_edg(i)


       ! loop over 2 start and end nodes on that edge
       do k = 1, 2 
          pt = grd%ibedge(i,k) ! global number of that node
          i1 = neqs * (pt - 1) + 1
          i2 = neqs * pt

          KK(i1:i2,:) = 0.0d0 ! freeze that row 
          do j = 1, neqs
             KK((i1+j-1),(i1+j-1)) = ave ! ones on the diagonal
             rhs(i1+j-1) = 0.0d0
          end do

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
                i1 = neqs * (pt - 1) + 1
                i2 = neqs * pt

                KK(i1:i2,:) = 0.0d0 ! freeze that row 
                do j = 1, neqs
                   KK((i1+j-1),(i1+j-1)) = ave ! ones on the diagonal
                   rhs(i1+j-1) = 0.0d0
                end do
             end do ! k - loop
          end if

       end if ! there exists bn pts on that edge

    end do ! done all boundary edges

    ! ! handling duplicate nodes (if any)
    ! if ( allocated(grd%dup_nodes) .and. ( size(grd%dup_nodes) >= 1) ) then
    !    dup => grd%dup_nodes
    !    do k = 1, grd%tot_repeated_bn_nodes
    !       pt1 = dup(2 * k - 1)
    !       pt2 = dup(  2 * k  )
    !       if(maxval(abs(KK(:,:,pt1,:))) .eq. 0.0d0) call swap(pt1, pt2)
    !       KK(:,:,pt2,:) = 0.0d0 ! freeze that row 
    !       rhs(:,pt2) = 0.0d0
    !       do j = 1, neqs
    !          KK(j,j,pt2,pt1) = 1.0d0
    !          KK(j,j,pt2,pt2) = -1.0d0
    !       end do
    !    end do
    ! end if

    ! done here                                         
  end subroutine freeze_bcs

end module time_integration
