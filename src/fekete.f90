module fekete
  implicit none

  private

  public :: fekete_order_num, fekete_rule
  public :: fekete_degree2rule

contains

  subroutine fekete_degree ( rule, degree )

    !*****************************************************************************80
    !
    !! FEKETE_DEGREE returns the degree of a Fekete rule for the triangle.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    02 October 2006
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Mark Taylor, Beth Wingate, Rachel Vincent,
    !    An Algorithm for Computing Fekete Points in the Triangle,
    !    SIAM Journal on Numerical Analysis,
    !    Volume 38, Number 5, 2000, pages 1707-1720.
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

    if ( rule == 1 ) then
       degree = 3
    else if ( rule == 2 ) then
       degree = 6
    else if ( rule == 3 ) then
       degree = 9
    else if ( rule == 4 ) then
       degree = 12
    else if ( rule == 5 ) then
       degree = 12
    else if ( rule == 6 ) then
       degree = 15
    else if ( rule == 7 ) then
       degree = 18
    else

       degree = -1
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'FEKETE_DEGREE - Fatal error!'
       write ( *, '(a,i8)' ) '  Illegal RULE = ', rule
       stop

    end if

    return
  end subroutine fekete_degree
  subroutine fekete_order_num ( rule, order_num )

    !*****************************************************************************80
    !
    !! FEKETE_ORDER_NUM returns the order of a Fekete rule for the triangle.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    02 October 2006
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Mark Taylor, Beth Wingate, Rachel Vincent,
    !    An Algorithm for Computing Fekete Points in the Triangle,
    !    SIAM Journal on Numerical Analysis,
    !    Volume 38, Number 5, 2000, pages 1707-1720.
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

    call fekete_suborder_num ( rule, suborder_num )

    allocate ( suborder(1:suborder_num) )

    call fekete_suborder ( rule, suborder_num, suborder )

    order_num = sum ( suborder(1:suborder_num) )

    deallocate ( suborder )

    return
  end subroutine fekete_order_num
  subroutine fekete_rule ( rule, order_num, xy, w )

    !*****************************************************************************80
    !
    !! FEKETE_RULE returns the points and weights of a Fekete rule.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    30 January 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Mark Taylor, Beth Wingate, Rachel Vincent,
    !    An Algorithm for Computing Fekete Points in the Triangle,
    !    SIAM Journal on Numerical Analysis,
    !    Volume 38, Number 5, 2000, pages 1707-1720.
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
    real ( kind = 8 ), allocatable, dimension ( : ) :: suborder_w
    real ( kind = 8 ), allocatable, dimension ( :, : ) :: suborder_xyz
    real ( kind = 8 ) w(order_num)
    real ( kind = 8 ) xy(2,order_num)
    !
    !  Get the suborder information.
    !
    call fekete_suborder_num ( rule, suborder_num )

    allocate ( suborder(suborder_num) )
    allocate ( suborder_xyz(3,suborder_num) )
    allocate ( suborder_w(suborder_num) )

    call fekete_suborder ( rule, suborder_num, suborder )

    call fekete_subrule ( rule, suborder_num, suborder_xyz, suborder_w )
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
          write ( *, '(a)' ) 'FEKETE_RULE - Fatal error!'
          write ( *, '(a,i8,a,i8)' ) '  Illegal SUBORDER(', s, ') = ', suborder(s)
          stop

       end if

    end do

    deallocate ( suborder )
    deallocate ( suborder_xyz )
    deallocate ( suborder_w )

    return
  end subroutine fekete_rule
  subroutine fekete_rule_num ( rule_num )

    !*****************************************************************************80
    !
    !! FEKETE_RULE_NUM returns the number of Fekete rules available.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    02 October 2006
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Mark Taylor, Beth Wingate, Rachel Vincent,
    !    An Algorithm for Computing Fekete Points in the Triangle,
    !    SIAM Journal on Numerical Analysis,
    !    Volume 38, Number 5, 2000, pages 1707-1720.
    !
    !  Parameters:
    !
    !    Output, integer ( kind = 4 ) RULE_NUM, the number of rules available.
    !
    implicit none

    integer ( kind = 4 ) rule_num

    rule_num = 7

    return
  end subroutine fekete_rule_num
  subroutine fekete_suborder ( rule, suborder_num, suborder )

    !*****************************************************************************80
    !
    !! FEKETE_SUBORDER returns the suborders for a Fekete rule.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    02 October 2006
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Mark Taylor, Beth Wingate, Rachel Vincent,
    !    An Algorithm for Computing Fekete Points in the Triangle,
    !    SIAM Journal on Numerical Analysis,
    !    Volume 38, Number 5, 2000, pages 1707-1720.
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
            1, 3, 6 /)
    else if ( rule == 2 ) then
       suborder(1:suborder_num) = (/ &
            1, 3, 3, 3, 6, 6, 6 /)
    else if ( rule == 3 ) then
       suborder(1:suborder_num) = (/ &
            1, 3, 3, 3, 3, 6, 6, 6, 6, 6, &
            6, 6 /)
    else if ( rule == 4 ) then
       suborder(1:suborder_num) = (/ &
            1, 3, 3, 3, 3, 3, 3, 6, 6, 6, &
            6, 6, 6, 6, 6, 6, 6, 6, 6 /)
    else if ( rule == 5 ) then
       suborder(1:suborder_num) = (/ &
            1, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
            3, 6, 6, 6, 6, 6, 6, 6, 6, 6, &
            6  /)
    else if ( rule == 6 ) then
       suborder(1:suborder_num) = (/ &
            1, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
            6, 6, 6, 6, 6, 6, 6, 6, 6, 6, &
            6, 6, 6, 6, 6, 6, 6, 6  /)
    else if ( rule == 7 ) then
       suborder(1:suborder_num) = (/ &
            1, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
            3, 3, 6, 6, 6, 6, 6, 6, 6, 6, &
            6, 6, 6, 6, 6, 6, 6, 6, 6, 6, &
            6, 6, 6, 6, 6, 6, 6, 6  /)
    else

       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'FEKETE_SUBORDER - Fatal error!'
       write ( *, '(a,i8)' ) '  Illegal RULE = ', rule
       stop

    end if

    return
  end subroutine fekete_suborder
  subroutine fekete_suborder_num ( rule, suborder_num )

    !*****************************************************************************80
    !
    !! FEKETE_SUBORDER_NUM returns the number of suborders for a Fekete rule.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    02 October 2006
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Mark Taylor, Beth Wingate, Rachel Vincent,
    !    An Algorithm for Computing Fekete Points in the Triangle,
    !    SIAM Journal on Numerical Analysis,
    !    Volume 38, Number 5, 2000, pages 1707-1720.
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
    integer ( kind = 4 ) suborder_num

    if ( rule == 1 ) then
       suborder_num = 3
    else if ( rule == 2 ) then
       suborder_num = 7
    else if ( rule == 3 ) then
       suborder_num = 12
    else if ( rule == 4 ) then
       suborder_num = 19
    else if ( rule == 5 ) then
       suborder_num = 21
    else if ( rule == 6 ) then
       suborder_num = 28
    else if ( rule == 7 ) then
       suborder_num = 38
    else

       suborder_num = -1
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'FEKETE_SUBORDER_NUM - Fatal error!'
       write ( *, '(a,i8)' ) '  Illegal RULE = ', rule
       stop

    end if

    return
  end subroutine fekete_suborder_num
  subroutine fekete_subrule ( rule, suborder_num, suborder_xyz, suborder_w )

    !*****************************************************************************80
    !
    !! FEKETE_SUBRULE returns a compressed Fekete rule.
    !
    !  Discussion:
    !
    !    The listed weights are twice what we want...since we want them
    !    to sum to 1/2, reflecting the area of a unit triangle.  So we
    !    simple halve the values before exiting this routine.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    30 January 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Mark Taylor, Beth Wingate, Rachel Vincent,
    !    An Algorithm for Computing Fekete Points in the Triangle,
    !    SIAM Journal on Numerical Analysis,
    !    Volume 38, Number 5, 2000, pages 1707-1720.
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
    integer ( kind = 4 ) s
    real ( kind = 8 ) suborder_w(suborder_num)
    real ( kind = 8 ) suborder_xyz(3,suborder_num)

    if ( rule == 1 ) then

       suborder_xyz(1:3,1:suborder_num) = reshape ( (/ &
            0.3333333333D+00,  0.3333333333D+00, 0.3333333334D+00, &
            1.0000000000D+00,  0.0000000000D+00, 0.0000000000D+00, &
            0.0000000000D+00,  0.2763932023D+00, 0.7236067977D+00  &
            /), (/ 3, suborder_num /) )

       suborder_w(1:suborder_num) = (/ &
            0.9000000000D+00, &
            0.0333333333D+00, &
            0.1666666667D+00 /)

    else if ( rule == 2 ) then

       suborder_xyz(1:3,1:suborder_num) = reshape ( (/ &
            0.3333333333D+00,  0.3333333333D+00, 0.3333333334D+00, &
            0.1063354684D+00,  0.1063354684D+00, 0.7873290632D+00, &
            0.5000000000D+00,  0.5000000000D+00, 0.0000000000D+00, &
            1.0000000000D+00,  0.0000000000D+00, 0.0000000000D+00, &
            0.1171809171D+00,  0.3162697959D+00, 0.5665492870D+00, &
            0.0000000000D+00,  0.2655651402D+00, 0.7344348598D+00, &
            0.0000000000D+00,  0.0848854223D+00, 0.9151145777D+00  &
            /), (/ 3, suborder_num /) )

       suborder_w(1:suborder_num) = (/ &
            0.2178563571D+00, &
            0.1104193374D+00, &
            0.0358939762D+00, &
            0.0004021278D+00, &
            0.1771348660D+00, &
            0.0272344079D+00, &
            0.0192969460D+00 /)

    else if ( rule == 3 ) then

       suborder_xyz(1:3,1:suborder_num) = reshape ( (/ &
            0.3333333333D+00,  0.3333333333D+00, 0.3333333334D+00, &
            0.1704318201D+00,  0.1704318201D+00, 0.6591363598D+00, &
            0.0600824712D+00,  0.4699587644D+00, 0.4699587644D+00, &
            0.0489345696D+00,  0.0489345696D+00, 0.9021308608D+00, &
            0.0000000000D+00,  0.0000000000D+00, 1.0000000000D+00, &
            0.1784337588D+00,  0.3252434900D+00, 0.4963227512D+00, &
            0.0588564879D+00,  0.3010242110D+00, 0.6401193011D+00, &
            0.0551758079D+00,  0.1543901944D+00, 0.7904339977D+00, &
            0.0000000000D+00,  0.4173602935D+00, 0.5826397065D+00, &
            0.0000000000D+00,  0.2610371960D+00, 0.7389628040D+00, &
            0.0000000000D+00,  0.1306129092D+00, 0.8693870908D+00, &
            0.0000000000D+00,  0.0402330070D+00, 0.9597669930D+00  &
            /), (/ 3, suborder_num /) )

       suborder_w(1:suborder_num) = (/ &
            0.1096011288D+00, &
            0.0767491008D+00, &
            0.0646677819D+00, &
            0.0276211659D+00, &
            0.0013925011D+00, &
            0.0933486453D+00, &
            0.0619010169D+00, &
            0.0437466450D+00, &
            0.0114553907D+00, &
            0.0093115568D+00, &
            0.0078421987D+00, &
            0.0022457501D+00 /)

    else if ( rule == 4 ) then

       suborder_xyz(1:3,1:suborder_num) = reshape ( (/ &
            0.3333333333D+00,  0.3333333333D+00, 0.3333333334D+00, &
            0.1988883477D+00,  0.4005558262D+00, 0.4005558261D+00, &
            0.2618405201D+00,  0.2618405201D+00, 0.4763189598D+00, &
            0.0807386775D+00,  0.0807386775D+00, 0.8385226450D+00, &
            0.0336975736D+00,  0.0336975736D+00, 0.9326048528D+00, &
            0.0000000000D+00,  0.5000000000D+00, 0.5000000000D+00, &
            0.0000000000D+00,  0.0000000000D+00, 1.0000000000D+00, &
            0.1089969290D+00,  0.3837518758D+00, 0.5072511952D+00, &
            0.1590834479D+00,  0.2454317980D+00, 0.5954847541D+00, &
            0.0887037176D+00,  0.1697134458D+00, 0.7415828366D+00, &
            0.0302317829D+00,  0.4071849276D+00, 0.5625832895D+00, &
            0.0748751152D+00,  0.2874821712D+00, 0.6376427136D+00, &
            0.0250122615D+00,  0.2489279690D+00, 0.7260597695D+00, &
            0.0262645218D+00,  0.1206826354D+00, 0.8530528428D+00, &
            0.0000000000D+00,  0.3753565349D+00, 0.6246434651D+00, &
            0.0000000000D+00,  0.2585450895D+00, 0.7414549105D+00, &
            0.0000000000D+00,  0.1569057655D+00, 0.8430942345D+00, &
            0.0000000000D+00,  0.0768262177D+00, 0.9231737823D+00, &
            0.0000000000D+00,  0.0233450767D+00, 0.9766549233D+00  &
            /), (/ 3, suborder_num /) )

       suborder_w(1:suborder_num) = (/ &
            0.0626245179D+00, &
            0.0571359417D+00, &
            0.0545982307D+00, &
            0.0172630326D+00, &
            0.0142519606D+00, &
            0.0030868485D+00, &
            0.0004270742D+00, &
            0.0455876390D+00, &
            0.0496701966D+00, &
            0.0387998322D+00, &
            0.0335323983D+00, &
            0.0268431561D+00, &
            0.0237377452D+00, &
            0.0177255972D+00, &
            0.0043097313D+00, &
            0.0028258057D+00, &
            0.0030994935D+00, &
            0.0023829062D+00, &
            0.0009998683D+00 /)

    else if ( rule == 5 ) then

       suborder_xyz(1:3,1:suborder_num) = reshape ( (/ &
            0.3333333333D+00,  0.3333333333D+00, 0.3333333334D+00, &
            0.2201371125D+00,  0.3169406831D+00, 0.4629222044D+00, &
            0.2201371125D+00,  0.4629222044D+00, 0.3169406831D+00, &
            0.1877171129D+00,  0.1877171129D+00, 0.6245657742D+00, &
            0.1403402144D+00,  0.4298298928D+00, 0.4298298928D+00, &
            0.0833252778D+00,  0.0833252778D+00, 0.8333494444D+00, &
            0.0664674598D+00,  0.0252297247D+00, 0.9083028155D+00, &
            0.0218884020D+00,  0.4890557990D+00, 0.4890557990D+00, &
            0.0252297247D+00,  0.0664674598D+00, 0.9083028155D+00, &
            0.0000000000D+00,  0.5000000000D+00, 0.5000000000D+00, &
            0.0000000000D+00,  0.0000000000D+00, 1.0000000000D+00, &
            0.1157463404D+00,  0.2842319093D+00, 0.6000217503D+00, &
            0.0672850606D+00,  0.3971764400D+00, 0.5355384994D+00, &
            0.0909839531D+00,  0.1779000668D+00, 0.7311159801D+00, &
            0.0318311633D+00,  0.3025963402D+00, 0.6655724965D+00, &
            0.0273518579D+00,  0.1733665506D+00, 0.7992815915D+00, &
            0.0000000000D+00,  0.3753565349D+00, 0.6246434651D+00, &
            0.0000000000D+00,  0.2585450895D+00, 0.7414549105D+00, &
            0.0000000000D+00,  0.1569057655D+00, 0.8430942345D+00, &
            0.0000000000D+00,  0.0768262177D+00, 0.9231737823D+00, &
            0.0000000000D+00,  0.0233450767D+00, 0.9766549233D+00  &
            /), (/ 3, suborder_num /) )

       suborder_w(1:suborder_num) = (/ &
            0.0485965670D+00, &
            0.0602711576D+00, &
            0.0602711576D+00, &
            0.0476929767D+00, &
            0.0453940802D+00, &
            0.0258019417D+00, &
            0.0122004614D+00, &
            0.0230003812D+00, &
            0.0122004614D+00, &
            0.0018106475D+00, &
            -0.0006601747D+00, &
            0.0455413513D+00, &
            0.0334182802D+00, &
            0.0324896773D+00, &
            0.0299402736D+00, &
            0.0233477738D+00, &
            0.0065962854D+00, &
            0.0021485117D+00, &
            0.0034785755D+00, &
            0.0013990566D+00, &
            0.0028825748D+00 /)

    else if ( rule == 6 ) then

       suborder_xyz(1:3,1:suborder_num) = reshape ( (/ &
            0.3333333333D+00,  0.3333333333D+00, 0.3333333334D+00, &
            0.2379370518D+00,  0.3270403780D+00, 0.4350225702D+00, &
            0.3270403780D+00,  0.2379370518D+00, 0.4350225702D+00, &
            0.1586078048D+00,  0.4206960976D+00, 0.4206960976D+00, &
            0.2260541354D+00,  0.2260541354D+00, 0.5478917292D+00, &
            0.1186657611D+00,  0.1186657611D+00, 0.7626684778D+00, &
            0.0477095725D+00,  0.4761452137D+00, 0.4761452138D+00, &
            0.0531173538D+00,  0.0531173538D+00, 0.8937652924D+00, &
            0.0219495841D+00,  0.0219495841D+00, 0.9561008318D+00, &
            0.0000000000D+00,  0.0000000000D+00, 1.0000000000D+00, &
            0.1585345951D+00,  0.3013819154D+00, 0.5400834895D+00, &
            0.0972525649D+00,  0.3853507643D+00, 0.5173966708D+00, &
            0.0875150140D+00,  0.2749910734D+00, 0.6374939126D+00, &
            0.1339547708D+00,  0.1975591066D+00, 0.6684861226D+00, &
            0.0475622627D+00,  0.3524012205D+00, 0.6000365168D+00, &
            0.0596194677D+00,  0.1978887556D+00, 0.7424917767D+00, &
            0.0534939782D+00,  0.1162464503D+00, 0.8302595715D+00, &
            0.0157189888D+00,  0.4176001732D+00, 0.5666808380D+00, &
            0.0196887324D+00,  0.2844332752D+00, 0.6958779924D+00, &
            0.0180698489D+00,  0.1759511193D+00, 0.8059790318D+00, &
            0.0171941515D+00,  0.0816639421D+00, 0.9011419064D+00, &
            0.0000000000D+00,  0.4493368632D+00, 0.5506631368D+00, &
            0.0000000000D+00,  0.3500847655D+00, 0.6499152345D+00, &
            0.0000000000D+00,  0.2569702891D+00, 0.7430297109D+00, &
            0.0000000000D+00,  0.1738056486D+00, 0.8261943514D+00, &
            0.0000000000D+00,  0.1039958541D+00, 0.8960041459D+00, &
            0.0000000000D+00,  0.0503997335D+00, 0.9496002665D+00, &
            0.0000000000D+00,  0.0152159769D+00, 0.9847840231D+00  &
            /), (/ 3, suborder_num /) )

       suborder_w(1:suborder_num) = (/ &
            0.0459710878D+00, &
            0.0346650571D+00, &
            0.0346650571D+00, &
            0.0384470625D+00, &
            0.0386013566D+00, &
            0.0224308157D+00, &
            0.0243531004D+00, &
            0.0094392654D+00, &
            0.0061105652D+00, &
            0.0001283162D+00, &
            0.0305412307D+00, &
            0.0262101254D+00, &
            0.0265367617D+00, &
            0.0269859772D+00, &
            0.0172635676D+00, &
            0.0188795851D+00, &
            0.0158224870D+00, &
            0.0127170850D+00, &
            0.0164489660D+00, &
            0.0120018620D+00, &
            0.0072268907D+00, &
            0.0023599161D+00, &
            0.0017624674D+00, &
            0.0018648017D+00, &
            0.0012975716D+00, &
            0.0018506035D+00, &
            0.0009919379D+00, &
            0.0004893506D+00 /)

    else if ( rule == 7 ) then

       suborder_xyz(1:3,1:suborder_num) = reshape ( (/ &
            0.3333333333D+00,  0.3333333333D+00, 0.3333333334D+00, &
            0.2515553103D+00,  0.3292984162D+00, 0.4191462735D+00, &
            0.3292984162D+00,  0.2515553103D+00, 0.4191462735D+00, &
            0.1801930996D+00,  0.4099034502D+00, 0.4099034502D+00, &
            0.2438647767D+00,  0.2438647767D+00, 0.5122704466D+00, &
            0.1512564554D+00,  0.1512564554D+00, 0.6974870892D+00, &
            0.0810689493D+00,  0.4594655253D+00, 0.4594655254D+00, &
            0.0832757649D+00,  0.0832757649D+00, 0.8334484702D+00, &
            0.0369065587D+00,  0.0369065587D+00, 0.9261868826D+00, &
            0.0149574850D+00,  0.0149574850D+00, 0.9700850300D+00, &
            0.0000000000D+00,  0.5000000000D+00, 0.5000000000D+00, &
            0.0000000000D+00,  0.0000000000D+00, 1.0000000000D+00, &
            0.1821465920D+00,  0.3095465041D+00, 0.5083069039D+00, &
            0.1246901255D+00,  0.3789288931D+00, 0.4963809814D+00, &
            0.1179441386D+00,  0.2868915642D+00, 0.5951642972D+00, &
            0.1639418454D+00,  0.2204868669D+00, 0.6155712877D+00, &
            0.0742549663D+00,  0.3532533654D+00, 0.5724916683D+00, &
            0.0937816771D+00,  0.2191980979D+00, 0.6870202250D+00, &
            0.0890951387D+00,  0.1446273457D+00, 0.7662775156D+00, &
            0.0409065243D+00,  0.4360543636D+00, 0.5230391121D+00, &
            0.0488675890D+00,  0.2795984854D+00, 0.6715339256D+00, &
            0.0460342127D+00,  0.2034211147D+00, 0.7505446726D+00, &
            0.0420687187D+00,  0.1359040280D+00, 0.8220272533D+00, &
            0.0116377940D+00,  0.4336892286D+00, 0.5546729774D+00, &
            0.0299062187D+00,  0.3585587824D+00, 0.6115349989D+00, &
            0.0132313129D+00,  0.2968103667D+00, 0.6899583204D+00, &
            0.0136098469D+00,  0.2050279257D+00, 0.7813622274D+00, &
            0.0124869684D+00,  0.1232146223D+00, 0.8642984093D+00, &
            0.0365197797D+00,  0.0805854893D+00, 0.8828947310D+00, &
            0.0118637765D+00,  0.0554881302D+00, 0.9326480933D+00, &
            0.0000000000D+00,  0.4154069883D+00, 0.5845930117D+00, &
            0.0000000000D+00,  0.3332475761D+00, 0.6667524239D+00, &
            0.0000000000D+00,  0.2558853572D+00, 0.7441146428D+00, &
            0.0000000000D+00,  0.1855459314D+00, 0.8144540686D+00, &
            0.0000000000D+00,  0.1242528987D+00, 0.8757471013D+00, &
            0.0000000000D+00,  0.0737697111D+00, 0.9262302889D+00, &
            0.0000000000D+00,  0.0355492359D+00, 0.9644507641D+00, &
            0.0000000000D+00,  0.0106941169D+00, 0.9893058831D+00  &
            /), (/ 3, suborder_num /) )

       suborder_w(1:suborder_num) = (/ &
            0.0326079297D+00, &
            0.0255331366D+00, &
            0.0255331366D+00, &
            0.0288093886D+00, &
            0.0279490452D+00, &
            0.0174438045D+00, &
            0.0203594338D+00, &
            0.0113349170D+00, &
            0.0046614185D+00, &
            0.0030346239D+00, &
            0.0012508731D+00, &
            0.0000782945D+00, &
            0.0235716330D+00, &
            0.0206304700D+00, &
            0.0204028340D+00, &
            0.0215105697D+00, &
            0.0183482070D+00, &
            0.0174161032D+00, &
            0.0155972434D+00, &
            0.0119269616D+00, &
            0.0147074804D+00, &
            0.0116182830D+00, &
            0.0087639138D+00, &
            0.0098563528D+00, &
            0.0096342355D+00, &
            0.0086477936D+00, &
            0.0083868302D+00, &
            0.0062576643D+00, &
            0.0077839825D+00, &
            0.0031415239D+00, &
            0.0006513246D+00, &
            0.0021137942D+00, &
            0.0004393452D+00, &
            0.0013662119D+00, &
            0.0003331251D+00, &
            0.0011613225D+00, &
            0.0004342867D+00, &
            0.0002031499D+00 /)

    else

       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'FEKETE_SUBRULE - Fatal error!'
       write ( *, '(a,i8)' ) '  Illegal RULE = ', rule
       stop

    end if
    !
    !  The listed weights are twice what we want!.
    !
    do s = 1, suborder_num
       suborder_w(s) = 0.5D+00 * suborder_w(s)
    end do

    return
  end subroutine fekete_subrule
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
    !    A "free" FORTRAN unit number is an integer between 1 and 99 which
    !    is not currently associated with an I/O device.  A free FORTRAN unit
    !    number is needed in order to open a file with the OPEN command.
    !
    !    If IUNIT = 0, then no free FORTRAN unit could be found, although
    !    all 99 units were checked (except for units 5, 6 and 9, which
    !    are commonly reserved for console I/O).
    !
    !    Otherwise, IUNIT is an integer between 1 and 99, representing a
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
    !    Input, integer ( kind = 4 ) ILO, IHI, the desired bounds.
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

  subroutine fekete_degree2rule( degree, rule)
    implicit none
    integer, intent(in) :: degree
    integer, intent(out) :: rule

    select case ( degree )
    case ( 3 )
       rule = 1
    case ( 6 )
       rule = 2
    case ( 9 )
       rule = 3
    case ( 12 )
       rule = 4
    case ( 15 )
       rule = 6
    case ( 18 )
       rule = 7
    case default
       print *, 'no fekete rule for the given degree(order)' &
            ,' of accuracy of the element! stop'
       stop
    end select

    ! done here
  end subroutine fekete_degree2rule

end module fekete
