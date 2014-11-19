module ncc_triangle
  implicit none

  private

  public :: ncc_triangle_order_num, ncc_triangle_rule

contains

  subroutine file_name_inc ( file_name )

    !*****************************************************************************80
    !
    !! FILE_NAME_INC increments a partially numeric filename.
    !
    !  Discussion:
    !
    !    It is assumed that the digits in the name, whether scattered or
    !    connected, represent a number that is to be increased by 1 on
    !    each call.  If this number is all 9's on input, the output number
    !    is all 0's.  Non-numeric letters of the name are unaffected.
    !
    !    If the name is empty, then the routine stops.
    !
    !    If the name contains no digits, the empty string is returned.
    !
    !  Example:
    !
    !      Input            Output
    !      -----            ------
    !      'a7to11.txt'     'a7to12.txt'
    !      'a7to99.txt'     'a8to00.txt'
    !      'a9to99.txt'     'a0to00.txt'
    !      'cat.txt'        ' '
    !      ' '              STOP!
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    14 September 2005
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input/output, character ( len = * ) FILE_NAME.
    !    On input, a character string to be incremented.
    !    On output, the incremented string.
    !
    implicit none

    character c
    integer ( kind = 4 ) change
    integer ( kind = 4 ) digit
    character ( len = * ) file_name
    integer ( kind = 4 ) i
    integer ( kind = 4 ) lens

    lens = len_trim ( file_name )

    if ( lens <= 0 ) then
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'FILE_NAME_INC - Fatal error!'
       write ( *, '(a)' ) '  The input string is empty.'
       stop
    end if

    change = 0

    do i = lens, 1, -1

       c = file_name(i:i)

       if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

          change = change + 1

          digit = ichar ( c ) - 48
          digit = digit + 1

          if ( digit == 10 ) then
             digit = 0
          end if

          c = char ( digit + 48 )

          file_name(i:i) = c

          if ( c /= '0' ) then
             return
          end if

       end if

    end do

    if ( change == 0 ) then
       file_name = ' '
       return
    end if

    return
  end subroutine file_name_inc
  subroutine get_unit ( iunit )

    !*****************************************************************************80
    !
    !! GET_UNIT returns a free FORTRAN unit number.
    !
    !  Discussion:
    !
    !    A "free" FORTRAN unit number is a value between 1 and 99 which
    !    is not currently associated with an I/O device.  A free FORTRAN unit
    !    number is needed in order to open a file with the OPEN command.
    !
    !    If IUNIT = 0, then no free FORTRAN unit could be found, although
    !    all 99 units were checked (except for units 5, 6 and 9, which
    !    are commonly reserved for console I/O).
    !
    !    Otherwise, IUNIT is a value between 1 and 99, representing a
    !    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
    !    are special, and will never return those values.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    18 September 2005
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Output, integer ( kind = 4 ) IUNIT, the free unit number.
    !
    implicit none

    integer ( kind = 4 ) i
    integer ( kind = 4 ) ios
    integer ( kind = 4 ) iunit
    logical lopen

    iunit = 0

    do i = 1, 99

       if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

          inquire ( unit = i, opened = lopen, iostat = ios )

          if ( ios == 0 ) then
             if ( .not. lopen ) then
                iunit = i
                return
             end if
          end if

       end if

    end do

    return
  end subroutine get_unit
  function i4_modp ( i, j )

    !*****************************************************************************80
    !
    !! I4_MODP returns the nonnegative remainder of I4 division.
    !
    !  Discussion:
    !
    !    If
    !      NREM = I4_MODP ( I, J )
    !      NMULT = ( I - NREM ) / J
    !    then
    !      I = J * NMULT + NREM
    !    where NREM is always nonnegative.
    !
    !    The MOD function computes a result with the same sign as the
    !    quantity being divided.  Thus, suppose you had an angle A,
    !    and you wanted to ensure that it was between 0 and 360.
    !    Then mod(A,360) would do, if A was positive, but if A
    !    was negative, your result would be between -360 and 0.
    !
    !    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
    !
    !  Example:
    !
    !        I     J     MOD I4_MODP    Factorization
    !
    !      107    50       7       7    107 =  2 *  50 + 7
    !      107   -50       7       7    107 = -2 * -50 + 7
    !     -107    50      -7      43   -107 = -3 *  50 + 43
    !     -107   -50      -7      43   -107 =  3 * -50 + 43
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    02 March 1999
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) I, the number to be divided.
    !
    !    Input, integer ( kind = 4 ) J, the number that divides I.
    !
    !    Output, integer ( kind = 4 ) I4_MODP, the nonnegative remainder when I is
    !    divided by J.
    !
    implicit none

    integer ( kind = 4 ) i
    integer ( kind = 4 ) i4_modp
    integer ( kind = 4 ) j
    integer ( kind = 4 ) value

    if ( j == 0 ) then
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'I4_MODP - Fatal error!'
       write ( *, '(a,i8)' ) '  Illegal divisor J = ', j
       stop
    end if

    value = mod ( i, j )

    if ( value < 0 ) then
       value = value + abs ( j )
    end if

    i4_modp = value

    return
  end function i4_modp
  function i4_wrap ( ival, ilo, ihi )

    !*****************************************************************************80
    !
    !! I4_WRAP forces an I4 to lie between given limits by wrapping.
    !
    !  Example:
    !
    !    ILO = 4, IHI = 8
    !
    !    I  Value
    !
    !    -2     8
    !    -1     4
    !     0     5
    !     1     6
    !     2     7
    !     3     8
    !     4     4
    !     5     5
    !     6     6
    !     7     7
    !     8     8
    !     9     4
    !    10     5
    !    11     6
    !    12     7
    !    13     8
    !    14     4
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    19 August 2003
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) IVAL, an integer value.
    !
    !    Input, integer ( kind = 4 ) ILO, IHI, the desired bounds for the
    !    integer value.
    !
    !    Output, integer ( kind = 4 ) I4_WRAP, a "wrapped" version of IVAL.
    !
    implicit none

    ! integer ( kind = 4 ) i4_modp
    integer ( kind = 4 ) i4_wrap
    integer ( kind = 4 ) ihi
    integer ( kind = 4 ) ilo
    integer ( kind = 4 ) ival
    integer ( kind = 4 ) jhi
    integer ( kind = 4 ) jlo
    integer ( kind = 4 ) value
    integer ( kind = 4 ) wide

    jlo = min ( ilo, ihi )
    jhi = max ( ilo, ihi )

    wide = jhi - jlo + 1

    if ( wide == 1 ) then
       value = jlo
    else
       value = jlo + i4_modp ( ival - jlo, wide )
    end if

    i4_wrap = value

    return
  end function i4_wrap
  subroutine ncc_triangle_degree ( rule, degree )

    !*****************************************************************************80
    !
    !! NCC_TRIANGLE_DEGREE returns the degree of an NCC rule for the triangle.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    29 January 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Peter Silvester,
    !    Symmetric Quadrature Formulae for Simplexes,
    !    Mathematics of Computation,
    !    Volume 24, Number 109, January 1970, pages 95-100.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) RULE, the index of the rule.
    !
    !    Output, integer ( kind = 4 ) DEGREE, the polynomial degree of exactness of
    !    the rule.
    !
    implicit none

    integer ( kind = 4 ) degree
    integer ( kind = 4 ) rule

    if ( 1 <= rule .and. rule <= 9 ) then

       degree = rule - 1

    else

       degree = -1
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'NCC_TRIANGLE_DEGREE - Fatal error!'
       write ( *, '(a,i8)' ) '  Illegal RULE = ', rule
       stop

    end if

    return
  end subroutine ncc_triangle_degree
  subroutine ncc_triangle_order_num ( rule, order_num )

    !*****************************************************************************80
    !
    !! NCC_TRIANGLE_ORDER_NUM returns the order of an NCC rule for the triangle.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    29 January 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Peter Silvester,
    !    Symmetric Quadrature Formulae for Simplexes,
    !    Mathematics of Computation,
    !    Volume 24, Number 109, January 1970, pages 95-100.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) RULE, the index of the rule.
    !
    !    Output, integer ( kind = 4 ) ORDER_NUM, the order (number of points)
    !    of the rule.
    !
    implicit none

    integer ( kind = 4 ) order_num
    integer ( kind = 4 ) rule
    integer ( kind = 4 ), allocatable, dimension ( : ) :: suborder
    integer ( kind = 4 ) suborder_num

    call ncc_triangle_suborder_num ( rule, suborder_num )

    allocate ( suborder(1:suborder_num) )

    call ncc_triangle_suborder ( rule, suborder_num, suborder )

    order_num = sum ( suborder(1:suborder_num) )

    deallocate ( suborder )

    return
  end subroutine ncc_triangle_order_num
  subroutine ncc_triangle_rule ( rule, order_num, xy, w )

    !*****************************************************************************80
    !
    !! NCC_TRIANGLE_RULE returns the points and weights of an NCC rule.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    09 December 2006
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Peter Silvester,
    !    Symmetric Quadrature Formulae for Simplexes,
    !    Mathematics of Computation,
    !    Volume 24, Number 109, January 1970, pages 95-100.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) RULE, the index of the rule.
    !
    !    Input, integer ( kind = 4 ) ORDER_NUM, the order (number of points)
    !    of the rule.
    !
    !    Output, real ( kind = 8 ) XY(2,ORDER_NUM), the points of the rule.
    !
    !    Output, real ( kind = 8 ) W(ORDER_NUM), the weights of the rule.
    !
    implicit none

    integer ( kind = 4 ) order_num

    ! integer ( kind = 4 ) i4_wrap
    integer ( kind = 4 ) k
    integer ( kind = 4 ) o
    integer ( kind = 4 ) rule
    integer ( kind = 4 ) s
    integer ( kind = 4 ), allocatable, dimension ( : ) :: suborder
    integer ( kind = 4 ) suborder_num
    real    ( kind = 8 ), allocatable, dimension ( : ) :: suborder_w
    real    ( kind = 8 ), allocatable, dimension ( :, : ) :: suborder_xyz
    real    ( kind = 8 ) w(order_num)
    real    ( kind = 8 ) xy(2,order_num)
    !
    !  Get the suborder information.
    !
    call ncc_triangle_suborder_num ( rule, suborder_num )

    allocate ( suborder(suborder_num) )
    allocate ( suborder_xyz(3,suborder_num) )
    allocate ( suborder_w(suborder_num) )

    call ncc_triangle_suborder ( rule, suborder_num, suborder )

    call ncc_triangle_subrule ( rule, suborder_num, suborder_xyz, suborder_w )
    !
    !  Expand the suborder information to a full order rule.
    !
    o = 0

    do s = 1, suborder_num

       if ( suborder(s) == 1 ) then

          o = o + 1
          xy(1:2,o) = suborder_xyz(1:2,s)
          w(o) = suborder_w(s)

       else if ( suborder(s) == 3 ) then

          do k = 1, 3
             o = o + 1
             xy(1,o) = suborder_xyz ( i4_wrap(k,  1,3), s )
             xy(2,o) = suborder_xyz ( i4_wrap(k+1,1,3), s )
             w(o) = suborder_w(s)
          end do

       else if ( suborder(s) == 6 ) then

          do k = 1, 3
             o = o + 1
             xy(1,o) = suborder_xyz ( i4_wrap(k,  1,3), s )
             xy(2,o) = suborder_xyz ( i4_wrap(k+1,1,3), s )
             w(o) = suborder_w(s)
          end do

          do k = 1, 3
             o = o + 1
             xy(1,o) = suborder_xyz ( i4_wrap(k+1,1,3), s )
             xy(2,o) = suborder_xyz ( i4_wrap(k,  1,3), s )
             w(o) = suborder_w(s)
          end do

       else

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'NCC_TRIANGLE_RULE - Fatal error!'
          write ( *, '(a,i8,a,i8)' ) '  Illegal SUBORDER(', s, ') = ', suborder(s)
          write ( *, '(a,i8)' ) '  RULE =    ', rule
          write ( *, '(a,i8)' ) '  ORDER_NUM = ', order_num
          stop

       end if

    end do

    deallocate ( suborder )
    deallocate ( suborder_xyz )
    deallocate ( suborder_w )

    return
  end subroutine ncc_triangle_rule
  subroutine ncc_triangle_rule_num ( rule_num )

    !*****************************************************************************80
    !
    !! NCC_TRIANGLE_RULE_NUM returns the number of NCC rules available.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    29 January 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Peter Silvester,
    !    Symmetric Quadrature Formulae for Simplexes,
    !    Mathematics of Computation,
    !    Volume 24, Number 109, January 1970, pages 95-100.
    !
    !  Parameters:
    !
    !    Output, integer ( kind = 4 ) RULE_NUM, the number of rules available.
    !
    implicit none

    integer ( kind = 4 ) rule_num

    rule_num = 9

    return
  end subroutine ncc_triangle_rule_num
  subroutine ncc_triangle_suborder ( rule, suborder_num, suborder )

    !*****************************************************************************80
    !
    !! NCC_TRIANGLE_SUBORDER returns the suborders for an NCC rule.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    29 January 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Peter Silvester,
    !    Symmetric Quadrature Formulae for Simplexes,
    !    Mathematics of Computation,
    !    Volume 24, Number 109, January 1970, pages 95-100.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) RULE, the index of the rule.
    !
    !    Input, integer ( kind = 4 ) SUBORDER_NUM, the number of suborders
    !    of the rule.
    !
    !    Output, integer ( kind = 4 ) SUBORDER(SUBORDER_NUM), the suborders
    !    of the rule.
    !
    implicit none

    integer ( kind = 4 ) suborder_num

    integer ( kind = 4 ) rule
    integer ( kind = 4 ) suborder(suborder_num)

    if ( rule == 1 ) then
       suborder(1:suborder_num) = (/ &
            1 /)
    else if ( rule == 2 ) then
       suborder(1:suborder_num) = (/ &
            3 /)
    else if ( rule == 3 ) then
       suborder(1:suborder_num) = (/ &
            3 /)
    else if ( rule == 4 ) then
       suborder(1:suborder_num) = (/ &
            3, 6, 1 /)
    else if ( rule == 5 ) then
       suborder(1:suborder_num) = (/ &
            6, 3, 3 /)
    else if ( rule == 6 ) then
       suborder(1:suborder_num) = (/ &
            3, 6, 6, 3, 3 /)
    else if ( rule == 7 ) then
       suborder(1:suborder_num) = (/ &
            6, 6, 3, 3, 6, 1 /)
    else if ( rule == 8 ) then
       suborder(1:suborder_num) = (/ &
            3, 6, 6, 3, 6, 6, 3, 3 /)
    else if ( rule == 9 ) then
       suborder(1:suborder_num) = (/ &
            6, 6, 3, 6, 6, 3, 6, 3, 3 /)
    else

       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'NCC_TRIANGLE_SUBORDER - Fatal error!'
       write ( *, '(a,i8)' ) '  Illegal RULE = ', rule
       stop

    end if

    return
  end subroutine ncc_triangle_suborder
  subroutine ncc_triangle_suborder_num ( rule, suborder_num )

    !*****************************************************************************80
    !
    !! NCC_TRIANGLE_SUBORDER_NUM returns the number of suborders for an NCC rule.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    29 January 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Peter Silvester,
    !    Symmetric Quadrature Formulae for Simplexes,
    !    Mathematics of Computation,
    !    Volume 24, Number 109, January 1970, pages 95-100.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) RULE, the index of the rule.
    !
    !    Output, integer ( kind = 4 ) SUBORDER_NUM, the number of suborders
    !    of the rule.
    !
    implicit none

    integer ( kind = 4 ) rule
    integer ( kind = 4 ), dimension(1:9) :: suborder = (/ &
         1, 1, 1, 3, 3, 5, 6, 8, 9 /)

    integer ( kind = 4 ) suborder_num

    if ( 1 <= rule .and. rule <= 9 ) then
       suborder_num = suborder(rule)
    else
       suborder_num = -1
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'NCC_TRIANGLE_SUBORDER_NUM - Fatal error!'
       write ( *, '(a,i8)' ) '  Illegal RULE = ', rule
       stop

    end if

    return
  end subroutine ncc_triangle_suborder_num
  subroutine ncc_triangle_subrule ( rule, suborder_num, suborder_xyz, suborder_w )

    !*****************************************************************************80
    !
    !! NCC_TRIANGLE_SUBRULE returns a compressed NCC rule.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    29 January 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Peter Silvester,
    !    Symmetric Quadrature Formulae for Simplexes,
    !    Mathematics of Computation,
    !    Volume 24, Number 109, January 1970, pages 95-100.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) RULE, the index of the rule.
    !
    !    Input, integer ( kind = 4 ) SUBORDER_NUM, the number of suborders
    !    of the rule.
    !
    !    Output, real ( kind = 8 ) SUBORDER_XYZ(3,SUBORDER_NUM),
    !    the barycentric coordinates of the abscissas.
    !
    !    Output, real ( kind = 8 ) SUBORDER_W(SUBORDER_NUM), the
    !    suborder weights.
    !
    implicit none

    integer ( kind = 4 ) suborder_num

    integer ( kind = 4 ) rule
    real ( kind = 8 ) suborder_w(suborder_num)
    integer ( kind = 4 ) suborder_w_n(suborder_num)
    integer ( kind = 4 ) suborder_w_d
    real ( kind = 8 ) suborder_xyz(3,suborder_num)
    integer ( kind = 4 ) suborder_xyz_n(3,suborder_num)
    integer ( kind = 4 ) suborder_xyz_d

    if ( rule == 1 ) then

       suborder_xyz_n(1:3,1:suborder_num) = reshape ( (/ &
            1,  1, 1  &
            /), (/ 3, suborder_num /) )

       suborder_xyz_d = 3

       suborder_w_n(1:suborder_num) = (/ &
            1 /)

       suborder_w_d = 1

    else if ( rule == 2 ) then

       suborder_xyz_n(1:3,1:suborder_num) = reshape ( (/ &
            1, 0, 0  &
            /), (/ 3, suborder_num /) )

       suborder_xyz_d = 1

       suborder_w_n(1:suborder_num) = (/ &
            1 /)

       suborder_w_d = 3

    else if ( rule == 3 ) then

       suborder_xyz_n(1:3,1:suborder_num) = reshape ( (/ &
            1, 1, 0  &
            /), (/ 3, suborder_num /) )

       suborder_xyz_d = 2

       suborder_w_n(1:suborder_num) = (/ &
            1 /)

       suborder_w_d = 3

    else if ( rule == 4 ) then

       suborder_xyz_n(1:3,1:suborder_num) = reshape ( (/ &
            3, 0, 0,  &
            2, 1, 0,  &
            1, 1, 1   &
            /), (/ 3, suborder_num /) )

       suborder_xyz_d = 3

       suborder_w_n(1:suborder_num) = (/ &
            4, 9, 54 /)

       suborder_w_d = 120

    else if ( rule == 5 ) then

       suborder_xyz_n(1:3,1:suborder_num) = reshape ( (/ &
            3, 1, 0,  &
            2, 2, 0,  &
            2, 1, 1   &
            /), (/ 3, suborder_num /) )

       suborder_xyz_d = 4

       suborder_w_n(1:suborder_num) = (/ &
            4, -1, 8 /)

       suborder_w_d = 45

    else if ( rule == 6 ) then

       suborder_xyz_n(1:3,1:suborder_num) = reshape ( (/ &
            5, 0, 0,  &
            4, 1, 0,  &
            3, 2, 0,  &
            3, 1, 1,  &
            2, 2, 1  &
            /), (/ 3, suborder_num /) )

       suborder_xyz_d = 5

       suborder_w_n(1:suborder_num) = (/ &
            11, 25, 25, 200, 25 /)

       suborder_w_d = 1008

    else if ( rule == 7 ) then

       suborder_xyz_n(1:3,1:suborder_num) = reshape ( (/ &
            5, 1, 0,  &
            4, 2, 0,  &
            4, 1, 1,  &
            3, 3, 0,  &
            3, 2, 1,  &
            2, 2, 2   &
            /), (/ 3, suborder_num /) )

       suborder_xyz_d = 6

       suborder_w_n(1:suborder_num) = (/ &
            36, -27, 72, 64, 72, -54 /)

       suborder_w_d = 840

    else if ( rule == 8 ) then

       suborder_xyz_n(1:3,1:suborder_num) = reshape ( (/ &
            7, 0, 0,  &
            6, 1, 0,  &
            5, 2, 0,  &
            5, 1, 1,  &
            4, 3, 0,  &
            4, 2, 1,  &
            3, 3, 1,  &
            3, 2, 2   &
            /), (/ 3, suborder_num /) )

       suborder_xyz_d = 7

       suborder_w_n(1:suborder_num) = (/ &
            1336, 2989, 3577, 32242, 2695, -6860, 44590, 3430 /)

       suborder_w_d = 259200

    else if ( rule == 9 ) then

       suborder_xyz_n(1:3,1:suborder_num) = reshape ( (/ &
            7, 1, 0,  &
            6, 2, 0,  &
            6, 1, 1,  &
            5, 3, 0,  &
            5, 2, 1,  &
            4, 4, 0,  &
            4, 3, 1,  &
            4, 2, 2,  &
            3, 3, 2   &
            /), (/ 3, suborder_num /) )

       suborder_xyz_d = 8

       suborder_w_n(1:suborder_num) = (/ &
            368, -468, 704, 1136, 832, -1083, 672, -1448, 1472 /)

       suborder_w_d = 14175

    else

       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'NCC_TRIANGLE_SUBRULE - Fatal error!'
       write ( *, '(a,i8)' ) '  Illegal RULE = ', rule
       stop

    end if

    suborder_xyz(1:3,1:suborder_num) = &
         real ( suborder_xyz_n(1:3,1:suborder_num), kind = 8 ) &
         / real ( suborder_xyz_d,                     kind = 8 )

    suborder_w(1:suborder_num) = &
         real ( suborder_w_n(1:suborder_num), kind = 8 ) &
         / real ( suborder_w_d,                 kind = 8 )

    return
  end subroutine ncc_triangle_subrule

end module ncc_triangle
