module globals
  implicit none
  ! accounts for the precision of real kind
  ! rk = 12 is double precision
  integer, parameter  :: rk = selected_real_kind(12)
  real(rk), parameter :: PI = 4.D0*DATAN(1.D0)
  type plt_spec
     character (len = 120) :: filename
     character (len = 250) :: title
     !500 data with 20 char limit for name
     character :: legends(500)*20 
     character (len = 120) :: data_format
     logical append_flag
  end type plt_spec
  logical, parameter :: verbose_flag = .true.
end module globals


