!> @ingroup Terrible
!> @author Arash Ghasemi and Mohammad Mahtabi
!> @brief
!> A Module to have the global variables
!> @details
!> This modules includes the definition of global variables
!> such as \f$ \pi \f$, rank of real variables, filenames, title
!> type the plot specifics, format of the data, etc.

module globals
  implicit none
  ! accounts for the precision of real kind
  ! rk = 12 is double precision
  integer, parameter  :: rk = selected_real_kind(12)  !< rank of the real variables for high precision
  real(rk), parameter :: PI = 4.D0*DATAN(1.D0)  !< Constant value of \f$ \pi = 4atan(1) \f$
  
  !> @brief
  !> Data types plotting results
  !> @details
  !> Includes various attributes for plotting
  !> such as: filename, title, legend, data format, append flag
  !> \f$ \alpha \f$ and \f$ \beta \f$ (Jacobi polynomial coefficients)
  !> \f$ \xi \f$ and \f$ \eta \f$
  !> init_lock (showing whether and element is initiated or not)
  type plt_spec
     character (len = 120) :: filename  !< Name of the file
     character (len = 250) :: title  !< Title of the plot
     !500 data with 20 char limit for name
     character :: legends(500)*20  !< Legends strings
     character (len = 120) :: data_format  !< A string indicating the format of the data
     logical append_flag  !< A flag indicating appending status.
  end type plt_spec
  logical, parameter :: verbose_flag = .true.
end module globals


