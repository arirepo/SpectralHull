module grid_vis
  implicit none
  private
  integer, parameter :: rk = selected_real_kind(12)

  public :: printout
  contains
    subroutine printout(elem, x, y, nq, nt)
      integer,  dimension(:), intent(in) :: elem
      real(rk), dimension(:), intent(in) :: x, y
      integer, intent(in) :: nq, nt

      integer :: i, fn, ct

      ct = 0
      fn = 12
      open(fn,file='elements')
      do i = 1, nt + nq
        if(i .le. nt)then
          write(fn,*)x(elem(ct + 1)),y(elem(ct + 1))
          write(fn,*)x(elem(ct + 2)),y(elem(ct + 2))
          write(fn,*)x(elem(ct + 3)),y(elem(ct + 3))
          write(fn,*)x(elem(ct + 1)),y(elem(ct + 1))
          write(fn,*)
          ct = ct + 3
        else
          write(fn,*)x(elem(ct + 1)),y(elem(ct + 1))
          write(fn,*)x(elem(ct + 2)),y(elem(ct + 2))
          write(fn,*)x(elem(ct + 3)),y(elem(ct + 3))
          write(fn,*)x(elem(ct + 4)),y(elem(ct + 4))
          write(fn,*)x(elem(ct + 1)),y(elem(ct + 1))
          write(fn,*)
          ct = ct + 4
        end if
      end do
      close(fn)
    end subroutine printout
end module grid_vis
