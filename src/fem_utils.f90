module fem_utils
! contains some utility subroutines for
! block structure to matrix and matrix to block
! structure conversion, linear solve and etc.
use grid_opt
use lag_basis

implicit none

private

type lin_solv_param
   character(len = 30) :: method
   ! you can add anything you want
   ! like arrays for residual etc.
   ! when GMRES or PCG is added.
   ! ..
   ! ..
   ! ..

end type lin_solv_param

public :: lin_solv_param, convert_to_matrix, convert_to_block_struct
public :: linear_solve, comp_T_plate_exact, write_cell_center_to_bedge 
public :: xy2rs, fem_interp_u

contains

  ! converts multidimensional blocked KK and rhs into
  ! matrix and vector form suitable for linear solve
  subroutine convert_to_matrix(KK, rhs, KK_mat, rhs_mat)
    implicit none
    real*8, dimension(:,:,:,:), intent(in) :: KK
    real*8, dimension(:,:), intent(in) :: rhs
    real*8, dimension(:,:), intent(out) :: KK_mat
    real*8, dimension(:), intent(out) :: rhs_mat

    integer :: i, j, imax, jmax, dk
    integer :: i1, i2, j1, j2

    ! hard reset
    KK_mat = 0.0d0;  rhs_mat = 0.0d0

    dk = size(KK,1)
    imax = size(KK,3)
    jmax = size(KK,4)

    ! convert stiffness matrix
    do i = 1, imax
       i1 = (i-1)* dk + 1
       i2 = i* dk
       do j = 1, jmax
          j1 = (j-1)* dk + 1
          j2 = j* dk
          KK_mat(i1:i2, j1:j2) = KK(:,:,i,j) 
       end do
    end do

    ! convert rhs
    do i = 1, imax
       i1 = (i-1)* dk + 1
       i2 = i* dk
       rhs_mat(i1:i2) = rhs(:,i) 
    end do

    ! done here
  end subroutine convert_to_matrix

  ! converts a vector to block structure
  ! after linear solve, we do this to 
  ! be able to plot in tecplot
  subroutine convert_to_block_struct(u_mat, u)
    implicit none
    real*8, dimension(:), intent(in)  :: u_mat
    real*8, dimension(:, :), intent(out)   :: u

    ! local vars
    integer :: i, i1, i2, imax, z

    z = size(u,1)
    imax = size(u,2)

    do i = 1, imax

       i1 = (i-1) * z + 1
       i2 = i * z

       u(:,i) = u_mat(i1:i2)

    end do

    ! done here
  end subroutine convert_to_block_struct

  ! implements a generic linear solver
  ! can be LU = using Lapack
  ! can be GMRES = my implementation
  ! or anything else
  subroutine linear_solve(A, b, x, param)
    implicit none
    real*8 , dimension(:,:), intent(in) :: A
    real*8 , dimension(:), intent(in) :: b
    real*8 , dimension(:), intent(out) :: x
    type(lin_solv_param), intent(inout) :: param

    ! local vars
    ! LAPACK variables
    integer :: INFO, LDA, LDB, N, NRHS
    real*8, dimension(size(A,1), size(A,2)) :: Atmp
    integer, dimension(size(A,1))  :: IPIV
    real*8, dimension(size(A,1),1) :: BB

    select case (param%method)

    case('LAPACK_LU')

       ! solve using LAPACK
       ! LAPACK variables
       N = size(A,1)
       NRHS = 1
       LDA = N
       LDB = N
       Atmp = A
       BB(:,1) = b
       call DGESV( N, NRHS, Atmp, LDA, IPIV, BB, LDB, INFO)
       if( INFO .ne. 0 ) then
          print *, ' something is wrong in the linear solver! stop'
          stop
       end if
       x = BB(:,1) 

    case default

       print *, 'no method was found to solve Ax=b'
       print *, 'error happened in linear_solve(..). stop!'
       stop

    end select

    ! done here
  end subroutine linear_solve

  ! exact solution for heat transfer problem
  subroutine comp_T_plate_exact(L, W, T1, T2, n_max, grd, T)
    implicit none
    real*8, intent(in) :: L, W, T1, T2
    integer, intent(in) :: n_max
    type(grid), intent(in) :: grd
    real*8, dimension(:,:), intent(out) :: T ! solution

    ! local vars
    integer :: i, n
    real*8, parameter :: PI = 4.D0*DATAN(1.0D0)
    real*8 :: x , y, tmp

    ! hard reset
    T = 0.0d0

    do i = 1, grd%nnodesg !comp exact sol for all nodes

       x = grd%x(i)
       y = grd%y(i)
       tmp = 0.0d0

       do n = 1, n_max
          tmp = tmp + dble((-1)**(n+1) + 1) / dble(n) * sin( dble(n) * PI * x/L) &
               * sinh( dble(n) * PI * y/L) / sinh( dble(n) * PI * W/L)
       end do

       T(:,i) = T1 + (T2 - T1) * 2.0d0 / PI * tmp 

       ! T(:,i) = cos(2.0d0 * x) * sinh(2.0d0 * y)

    end do

    ! done here
  end subroutine comp_T_plate_exact
  
  ! given a boundary tag, prints the cell centered quantity 
  ! at the middle of those boundary edges.
  ! the output is written to the selected file.
  !
  ! Note : this is valid only for constant elements
  !        where cell centerd values are constant
  !        across the elements. 
  !
  ! Note: can be used to write the cell-centerd gradients 
  ! of linear elements since gradients are all constant
  ! across the element.
  
  subroutine write_cell_center_to_bedge(grd, u, btag, outfile)
    implicit none
    type(grid), intent(in) :: grd
    ! u(z,ncellsg) : cell-based u
    real*8, dimension(:,:), intent(in) :: u
    integer, intent(in) :: btag
    character(len = *), intent(in) :: outfile

    ! local vars
    integer :: i, j, pt1, pt2, z, the_cell
    real *8 :: xc, yc ! the edge center coords

    z = size(u,1) ! degree of freedom

    open(10, file = outfile, status = 'unknown')

    ! start writing boundary edge values
    do i = 1, grd%nbedgeg

       if (grd%ibedgeBC(i) .eq. btag) then 
          pt1 = grd%ibedge(i,1); pt2 = grd%ibedge(i,2)
          xc = (grd%x(pt1) + grd%x(pt2)) / 2.0d0
          yc = (grd%y(pt1) + grd%y(pt2)) / 2.0d0 
          the_cell = grd%ibedgeELEM(i) 

          write (10, '(F30.17, F30.17)', advance = 'no') xc, yc

          do j = 1, z
             write (10, '(F30.17, A)', advance = 'no') u(j,the_cell), ' ' 
          end do

          write (10, *)

       end if ! done that bedge

    end do ! done all bedges

    close(10) ! close the output

    ! done here
  end subroutine write_cell_center_to_bedge
  
  ! computes local coords (r,s) given the
  ! global coordinates (x,y) using a robust
  ! Newton method.
  ! 
  ! This should be applied consistently to
  ! higher-order elements also.
  
  subroutine xy2rs(x, y, elem, grd, tol, maxitr, tolrs, r, s)
    implicit none
    real*8, intent(in) :: x, y
    integer, intent(in) :: elem
    type(grid), intent(in) :: grd
    real*8, intent(in) :: tol
    integer, intent(in) :: maxitr
    real*8, intent(in) :: tolrs
    real*8, intent(out) :: r, s

    ! local vars
    integer :: i, j, npt, pt ! points per triangle
    real*8 :: val
    real*8, dimension(2), save :: der, F, delta
    real*8, dimension(2,2), save:: JJ, Jinv


    ! initial guess
    r = 0.25d0; s =0.25d0
    npt = size(grd%icon,2)

    do j = 1, maxitr

       F = (/ x, y /)
       JJ = 0.0d0; Jinv = 0.0d0
       val = 0.0d0; der = 0.0d0

       do i = 1, npt
          pt = grd%icon(elem, i)
          call psi(1, i, r, s, val, der)
          F = F - (/ (grd%x(pt) * val) , (grd%y(pt) * val) /)
          JJ(1,1) =  JJ(1,1) + grd%x(pt) * der(1); JJ(1,2) =  JJ(1,2) + grd%x(pt) * der(2)
          JJ(2,1) =  JJ(2,1) + grd%y(pt) * der(1); JJ(2,2) =  JJ(2,2) + grd%y(pt) * der(2)
       end do

       Jinv(1,1) = JJ(2,2)
       Jinv(2,2) = JJ(1,1)
       Jinv(2,1) = -JJ(2,1)
       Jinv(1,2) = -JJ(1,2)
       Jinv = 1.0d0 /(JJ(1,1) * JJ(2,2) - JJ(2,1) * JJ(1,2)) * Jinv
       delta = matmul(Jinv, F)

       r = r + delta(1)
       s = s + delta(2)

       if ( sqrt(sum(delta*delta)) .le. tol ) then
          ! robust bug checking before returning
          if ( (npt .eq. 3) .and. (j > 2) ) then
             print *, 'number of itr is ', j, ' and is greater than 2 ' &
                    , 'for linear element with number', elem, '! stop.'
             stop
          else if (  (r > (1.0d0 + tolrs)) .or. (s > (1.0d0 + tolrs)) &
               .or.  (r < (0.0d0 - tolrs)) .or. (s < (0.0d0 - tolrs)) ) then
             print *, '(r,s) = ', r, s ,' out of range! stop.'
             stop 
          end if
          ! print *, 'j = ', j , 'error = ', sqrt(sum(delta*delta))
          return
       end if

    end do

    print *, 'the xy2rs(...) subroutine did not converge' &
         , ' to desired tolerance in ', maxitr, 'iterations' &
         , ' at point (', x, y, ')! stop.'
    stop

    ! done here
  end subroutine xy2rs

  ! evaluates the value of the primary var "u"
  ! at somewhere (r,s) inside the element "elem"
  ! in grid "grd". The result is stored in "uout" 
  subroutine fem_interp_u(r, s, elem, grd, u, uout)
    implicit none
    real*8, intent(in) :: r, s
    integer, intent(in) :: elem
    type(grid), intent(in) :: grd
    real*8, dimension(:,:), intent(in) :: u
    real*8, dimension(:), intent(out) :: uout

    ! local vars
    integer :: i, pt, pts
    real*8 :: val
    real*8, dimension(2), save :: der

    pts = size(grd%icon,2)

    uout = 0.0d0

    do i = 1, pts
       pt = grd%icon(elem, i)
       call psi(1, i, r, s, val, der)
       uout = uout + val * u(:,pt)
    end do

    ! done here
  end subroutine fem_interp_u
 
end module fem_utils
