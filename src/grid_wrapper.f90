module grid_wrapper
! testing git push 
  use spline
  use nurbs
  implicit none
  private

  integer, parameter :: rk = selected_real_kind(12)
  integer, parameter :: spline_i = 1
  integer, parameter :: nurbs_i = 2

  public :: boundary_curve_alloc
  public :: boundary_comp
  public :: boundary_eval
  public :: spline_i, nurbs_i

  contains
    subroutine boundary_curve_alloc(ct, cmx, cmy, ca, cb, cc, cd, cp, dp, x, y, np, bt)
      real(rk), dimension(:), allocatable, intent(in out) :: ct
      real(rk), dimension(:), allocatable, intent(in out) :: cmx
      real(rk), dimension(:), allocatable, intent(in out) :: cmy
      real(rk), dimension(:), allocatable, intent(in out) :: ca
      real(rk), dimension(:), allocatable, intent(in out) :: cb
      real(rk), dimension(:), allocatable, intent(in out) :: cc
      real(rk), dimension(:), allocatable, intent(in out) :: cd
      real(rk), dimension(:), allocatable, intent(in out) :: cp
      real(rk), dimension(:), allocatable, intent(in out) :: dp
      real(rk), dimension(:), allocatable, intent(in out) :: x
      real(rk), dimension(:), allocatable, intent(in out) :: y
      integer, intent(in) :: np
      integer, optional :: bt
      integer :: btype, i, p
      if(present(bt))then
        btype = bt
      else
        btype = spline_i
      end if

      if(.not. allocated(cmx))then
        if(btype == spline_i)then
          allocate(ct(np))
          ct = (/ (dble(i), i = 1, np) /)
          if(.not. allocated(x))then
            allocate(x(np))
          end if
          if(.not. allocated(y))then
            allocate(y(np))
          end if
          allocate(cmx(np), cmy(np), ca(np), cb(np), &
          & cc(np), cd(np), cp(np), dp(np))
        else if(btype == nurbs_i)then
          p = 3
          if(np < 3)then
            write(*,*)'Error in subroutine boundary_curve_alloc:'&
                       ,' cannot use nurbs for curve with fewer than 3 points'
            stop
          end if
          do while((np - p) < 1)
            p = p - 1
          end do

          allocate(ct(np))
          ct = (/ (dble(i), i = 1, np) /) ! x, y parametrization in range [1, np]
          if(.not. allocated(x))then
            allocate(x(np))
          end if
          if(.not. allocated(y))then
            allocate(y(np))
          end if
          allocate(cmx(np), cmy(np))      ! control points
          allocate(ca(np))                ! weights
          allocate(cb(np+p+1))            ! knot vector
          allocate(cc(np))                ! x, y parametrization in range [0, 1]
          allocate(cd(1), cp(1), dp(1))
        else
          write(*,*)'error: boundary type neither spline nor nurbs chosen'
        end if
      end if
    end subroutine boundary_curve_alloc

    subroutine boundary_comp(x, y, a, b, c, d, cp, dp, Mx, My, t, mode, bt)
      real(rk), dimension(:), intent(in out) :: x
      real(rk), dimension(:), intent(in out) :: y
      real(rk), dimension(:), intent(in out) :: a
      real(rk), dimension(:), intent(in out) :: b
      real(rk), dimension(:), intent(in out) :: c
      real(rk), dimension(:), intent(in out) :: d
      real(rk), dimension(:), intent(in out) :: cp
      real(rk), dimension(:), intent(in out) :: dp
      real(rk), dimension(:), intent(in out) :: Mx
      real(rk), dimension(:), intent(in out) :: My
      real(rk), dimension(:), intent(in out) :: t
      character(len = *), intent(in) :: mode
      integer, optional :: bt
      integer :: btype
      if(present(bt))then
        btype = bt
      else
        btype = spline_i
      end if

      if(btype == spline_i)then
        call comp_spline_weights(x, t, Mx, a, b, c, d, cp, dp, mode)
        call comp_spline_weights(y, t, My, a, b, c, d, cp, dp, mode)
      else if(btype == nurbs_i)then
        call parametrize_points(0.5_rk, x, y, c)
        call comp_nurbs(x, y, c, Mx, My, a, b)
      end if
    end subroutine boundary_comp

    subroutine boundary_eval(val, x, y, a, b, c, Mx, My, t, xv, yv, opt, bt)
      real(rk), intent(in) :: val
      real(rk), dimension(:), intent(in) :: x
      real(rk), dimension(:), intent(in) :: y
      real(rk), dimension(:), intent(in) :: a
      real(rk), dimension(:), intent(in) :: b
      real(rk), dimension(:), intent(in) :: c
      real(rk), dimension(:), intent(in) :: Mx
      real(rk), dimension(:), intent(in) :: My
      real(rk), dimension(:), intent(in) :: t
      real(rk), intent(in out) :: xv
      real(rk), intent(in out) :: yv
      character(len = *), intent(in) :: opt
      integer, optional :: bt
      integer :: np, idx, i, btype
      real(rk) :: inv

      if(present(bt))then
        btype = bt
      else
        btype = spline_i
      end if

      if(btype == spline_i)then
        call spline_eval2(t, x, val, xv, opt, Mx)
        call spline_eval2(t, y, val, yv, opt, My)
      else if(btype == nurbs_i)then
        inv = (val - t(1)) / (t(size(t)) - t(1))
        call nurbs_eval(Mx, My, np, inv, xv, yv, a, b, 0)
      end if
    end subroutine boundary_eval
end module grid_wrapper
