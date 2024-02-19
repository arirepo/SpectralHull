module spem_2d
  ! SPectral Element in 2D
  !
  use grid_opt ! for grid operations
  use element_opt
  use fem_utils
  use bcs
  use condition
  use elasticity2D

  implicit none

  private

  type fem_struct
     integer :: tag
     real*8, dimension(:,:) , allocatable :: u, dudx, dudy
     type(grid) :: grd
     integer :: neqs
     type(element), dimension(:) , allocatable :: elems
     real*8, dimension(:,:,:,:), allocatable :: KK
     real*8, dimension(:,:), allocatable :: rhs
     real*8, dimension(:,:), allocatable :: KK_mat
     real*8, dimension(:), allocatable :: rhs_mat, xx
     type(lin_solv_param) :: param
	 integer :: fem_type	!< Indicates the type of problem: 1 for 2D plane stress
							!! linear elasticity, 2: acoustic wave

  end type fem_struct

  public :: fem_struct, fem_init, fem_solve

contains

!> @brief Initializes the FEM object and its cell objects by allocating
!! the solution variables, reseting them, and writing them into Tecplot format.
!> @details Initializes an FEM object by setting the \n
!! neqs: # of equations \n
!! lin_solve_method: method of solving Ax=b. Foe example 'LAPACK_LU', etc. \n
!! tag: is a specific tag of this FEM object in case multiple existed. \n
!!


  subroutine fem_init(fem, neqs, lin_solve_method, tag, fem_type, mater)
    implicit none
    type(fem_struct), target :: fem   !< The FEM object to be initialized.
    integer, intent(in) :: neqs       !< # of equations of this FEM problem.
    character(len = *), intent(in) :: lin_solve_method  !< Method of solving Ax=b.
    integer , intent(in) :: tag     !< Specific tag for this FEM object.
    integer , intent(in) :: fem_type     !< Type of problem
    type(isoMat2D), intent(in) :: mater   !< Material properties for solid mechanics



    ! local vars
    integer :: i
    type(grid), pointer :: grd

    if(fem%grd%nnodesg .eq. 0) then ! fatal
       print *, 'error : fem%grd%nnodesg = 0! grid must be assigned'&
            , ' to fem_struct before initializing it! stopped.'
       stop
    end if

    !
    grd => fem%grd

    ! number of equations
    fem%neqs = neqs

    ! tag of this FEM object.
    fem%tag = tag

    ! type of the analysis of this FEM object.
    fem%fem_type = fem_type


    ! allocate solution vars
    allocate(fem%u(neqs,grd%nnodesg), fem%dudx(neqs,grd%nnodesg) &
         , fem%dudy(neqs,grd%nnodesg))
    ! hard RESET !
    fem%u = 0.0d0; fem%dudx = 0.0d0; fem%dudy = 0.0d0

    call write_u_tecplot('fem.dat', fem%grd, fem%u)
    ! ! initialize the elements
    ! call init_material_props(nu = .33d0, E = 10600000.0d0)

    allocate(fem%elems(grd%ncellsg))
    if(fem_type .eq. 1) then !Plane stress linear elasticity
      do i = 1, grd%ncellsg
         call init_elem(fem%elems(i), i, grd, neqs)
         call comp_2D_elem_K(fem%elems(i), i , grd, mater)
         ! call fill_Q(grd, fem%elems(i), i)
      end do
    elseif(fem_type .eq. 2) then !Acoustic wave analysi
      do i = 1, grd%ncellsg
         call init_elem(fem%elems(i), i, grd, neqs)
         call comp_elem_KQf(fem%elems(i), i , grd)
         ! call fill_Q(grd, fem%elems(i), i)
      end do
    else
      print *, "Unidentified type of analysis in fem_init. Stop!"
      stop
    end if

    ! compute grid area
    call comp_grid_area(fem)
    ! stop

    ! proceeding to assemble all
    allocate(fem%KK(neqs, neqs, grd%nnodesg, grd%nnodesg), fem%rhs(neqs, grd%nnodesg))
    call assemble_all(fem%elems, grd, fem%KK, fem%rhs)

    ! do i = 1, grd%nnodesg
    !    print *, 'node i = ', i, 'rhs = ', fem%rhs(1, i)
    ! end do

    ! imposing boundary conditions for fem grid
    ! NOTE : additional BCs coming from synchronization
    ! between fem and ps blocks are imposed dynamically
    call impose_bcs(grd, fem%KK, fem%rhs)

    ! call comp_bn_length(grd, fem%elems)

    ! allocate matrix vars before filling
    allocate(fem%KK_mat(neqs * grd%nnodesg, neqs * grd%nnodesg) &
         ,  fem%rhs_mat(neqs * grd%nnodesg), fem%xx(neqs * grd%nnodesg))

    ! setting up the linear solver algorithm
    fem%param%method = lin_solve_method

    ! extra addition
    call convert_to_matrix(fem%KK, fem%rhs, fem%KK_mat, fem%rhs_mat)

    ! done here
  end subroutine fem_init

  subroutine comp_grid_area(fem)
    implicit none
    type(fem_struct), intent(in) :: fem

    ! local vars
    integer :: i
    real*8 :: area

    ! init
    area = 0.0d0

    do i = 1, fem%grd%ncellsg

       select case (fem%grd%elname(i))

       case ( GEN_QUADRI)
          area = area + sum(fem%elems(i)%W * fem%elems(i)%JJ)
       case ( GEN_TRIANGLE)
          area = area + 0.5d0 * sum(fem%elems(i)%W * fem%elems(i)%JJ)
       case default
          print *, 'Could not recognize element name when computing' &
               , ' total grid area! stop'
          stop

       end select

    end do

    ! report
    print *, 'Total grid area is : ', area

    ! done here
  end subroutine comp_grid_area

  ! solves the fem region
  ! the result is stored in fem itself (see fem%u)
  subroutine fem_solve(fem, echo)
    implicit none
    type(fem_struct) :: fem
    logical, intent(in) :: echo

    ! local vars
    real*8 :: cond_num

    ! first convert to full matrix
    ! can be substituted by a custom GMRES solve
    ! which solves in block form so fem%KK, fem%rhs are
    ! ready to go in that case
    call convert_to_matrix(fem%KK, fem%rhs, fem%KK_mat, fem%rhs_mat)

    ! solve full matrix
    print *, 'fem%rhs_mat = ', fem%rhs_mat
    call condition_hager( size(fem%KK_mat,1), fem%KK_mat, cond_num )
    print *, 'cond_num = ', cond_num, 'npe = ', fem%grd%npe(1)
    ! print *, 'npe = ', fem%grd%npe(1)

    ! stop
    call linear_solve(fem%KK_mat, fem%rhs_mat, fem%xx, fem%param)

    ! convert the solution back to block structure
    call convert_to_block_struct(fem%xx, fem%u)

    if (echo) print *, 'fem_solve completed for fem region', fem%tag, ' successfully!'

    call print_cp(fem%u, fem%grd, fem%elems)

    ! done here
  end subroutine fem_solve

end module spem_2d

! program spem_2d
!   ! SPectral Element in 2D
!   !
!   use grid_opt ! for grid operations
!   use element_opt
!   use fem_utils
!   use bcs

!   implicit none

!   ! local vars
!   integer :: i
!   real*8, dimension(:,:) , allocatable :: u, dudx, dudy, u_all
!   type(grid) :: grd
!   integer :: neqs
!   type(element), dimension(:) , allocatable :: elems
!   real*8, dimension(:,:,:,:), allocatable :: KK
!   real*8, dimension(:,:), allocatable :: rhs
!   real*8, dimension(:,:), allocatable :: KK_mat
!   real*8, dimension(:), allocatable :: rhs_mat, xx
!   type(lin_solv_param) :: param
!   real*8, dimension(:,:,:) , allocatable :: cell_grad_x, cell_grad_y
!   real*8, dimension(:,:) , allocatable :: grad_u_all

!   ! number of equations
!   neqs = 1

!   ! read the grid from TETREX format
!   call read_grid_TETREX('../grids/arash.grd', grd)
!   ! show grid properties for double check
!   !call print_grid_props(grd)
!   ! allocate solution vars
!   allocate(u(neqs,grd%nnodesg), dudx(neqs,grd%nnodesg) &
!          , dudy(neqs,grd%nnodesg), u_all(3*neqs,grd%nnodesg))
!   u = 0.0d0; dudx = 0.0d0; dudy = 0.0d0; u_all = 0.0d0

!   ! initialize the elements
!   call init_material_props(nu = .33d0, E = 10600000.0d0)
!   allocate(elems(grd%ncellsg))
!   do i = 1, grd%ncellsg
!      call init_elem(elems(i), i, grd, neqs)
!      call comp_elem_KQf(elems(i), i , grd)
!   end do
!   ! proceeding to assemble all
!   allocate(KK(neqs, neqs, grd%nnodesg, grd%nnodesg), rhs(neqs, grd%nnodesg))
!   call assemble_all(elems, grd, KK, rhs)

!   ! imposing boundary conditions
!   call impose_bcs(grd, KK, rhs)

!   ! allocate matrix vars before filling
!   allocate(KK_mat(neqs* grd%nnodesg, neqs* grd%nnodesg), rhs_mat(neqs* grd%nnodesg), xx(neqs* grd%nnodesg) )
!   call convert_to_matrix(KK, rhs, KK_mat, rhs_mat)

!   ! solve the system
!   param%method = 'LAPACK_LU'
!   call linear_solve(KK_mat, rhs_mat, xx, param)

!   ! convert the solution back to block structure
!   call convert_to_block_struct(xx, u)

!   ! write down the solution itself
!   call write_u_tecplot('../runs/final.dat', grd, u)

!   ! print *, C11, C12, C21, C22, C66

!   ! ! comp the exact solution (if required)
!   ! call comp_T_plate_exact(2.0d0, 1.0d0, 323.15D0, 423.15D0, 100, grd, dudy)

!   ! ! store in u_all form
!   ! u_all(1::3, :) = u
!   ! u_all(2::3, :) = dudy
!   ! u_all(3::3, :) = abs(dudy - u)
!   ! call write_u_tecplot('../runs/u_all.dat', grd, u_all)

!   ! ! show error arithmatic average
!   ! print *, 'error arith. average : ', sum(abs(dudy - u))/dble(size(abs(dudy - u)))

!   ! ! sample section for cell-based gradients
!   ! allocate(cell_grad_x(neqs, elems(1)%ngauss, grd%ncellsg) &
!   !      , cell_grad_y(neqs, elems(1)%ngauss, grd%ncellsg) )
!   ! allocate( grad_u_all(neqs, grd%ncellsg) )
!   ! call comp_grad_u(u, grd, elems, cell_grad_x, cell_grad_y)
!   ! grad_u_all(1, :) = C11 * cell_grad_x(1,1,:) + C12 * cell_grad_y(2,1,:)
!   ! ! grad_u_all(2, :) = 1 - (cell_grad_x(1,1,:)**2.0d0 + cell_grad_y(1,1,:)**2.0d0) !-cell_grad_x(1,1,:)
!   ! grad_u_all(2, :) = C21 * cell_grad_x(1,1,:) + C22 * cell_grad_y(2,1,:)
!   ! call write_u_cell_centered_tecplot('../runs/grad_test.dat', grd, grad_u_all)
!   ! !
!   ! call write_cell_center_to_bedge(grd, grad_u_all, 5, '../runs/boundary_elasticity_5.dat')
!   ! call write_cell_center_to_bedge(grd, grad_u_all, 2, '../runs/boundary_elasticity_2.dat')

!   print *, 'Program terminated successfully!'
!   ! done here

! end program spem_2d
