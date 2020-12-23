subroutine qag ( nth, nph, d1, d2, f, tf, jj, a, b, epsabs, epsrel, key, result, abserr, neval, ier )

    !*****************************************************************************80
    !
    !! QAG approximates an integral over a finite interval.
    !
    !  Discussion:
    !
    !    The routine calculates an approximation RESULT to a definite integral   
    !      I = integral of F over (A,B),
    !    hopefully satisfying
    !      || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
    !
    !    QAG is a simple globally adaptive integrator using the strategy of 
    !    Aind (Piessens, 1973).  It is possible to choose between 6 pairs of
    !    Gauss-Kronrod quadrature formulae for the rule evaluation component. 
    !    The pairs of high degree of precision are suitable for handling
    !    integration difficulties due to a strongly oscillating integrand.
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, external real ( kind = 4 ) F, the name of the function routine, of the form
    !      function f ( x )
    !      real ( kind = 4 ) f
    !      real ( kind = 4 ) x
    !    which evaluates the integrand function.
    !
    !    Input, real ( kind = 4 ) A, B, the limits of integration.
    !
    !    Input, real ( kind = 4 ) EPSABS, EPSREL, the absolute and relative accuracy requested.
    !
    !    Input, integer ( kind = 4 ) KEY, chooses the order of the local integration rule:
    !    1,  7 Gauss points, 15 Gauss-Kronrod points,
    !    2, 10 Gauss points, 21 Gauss-Kronrod points,
    !    3, 15 Gauss points, 31 Gauss-Kronrod points,
    !    4, 20 Gauss points, 41 Gauss-Kronrod points,
    !    5, 25 Gauss points, 51 Gauss-Kronrod points,
    !    6, 30 Gauss points, 61 Gauss-Kronrod points.
    !
    !    Output, real ( kind = 4 ) RESULT, the estimated value of the integral.
    !
    !    Output, real ( kind = 4 ) ABSERR, an estimate of || I - RESULT ||.
    !
    !    Output, integer ( kind = 4 ) NEVAL, the number of times the integral was evaluated.
    !
    !    Output, integer ( kind = 4 ) IER, return code.
    !    0, normal and reliable termination of the routine.  It is assumed that the 
    !      requested accuracy has been achieved.
    !    1, maximum number of subdivisions allowed has been achieved.  One can 
    !      allow more subdivisions by increasing the value of LIMIT in QAG. 
    !      However, if this yields no improvement it is advised to analyze the
    !      integrand to determine the integration difficulties.  If the position
    !      of a local difficulty can be determined, such as a singularity or
    !      discontinuity within the interval) one will probably gain from 
    !      splitting up the interval at this point and calling the integrator 
    !      on the subranges.  If possible, an appropriate special-purpose 
    !      integrator should be used which is designed for handling the type 
    !      of difficulty involved.
    !    2, the occurrence of roundoff error is detected, which prevents the
    !      requested tolerance from being achieved.
    !    3, extremely bad integrand behavior occurs at some points of the
    !      integration interval.
    !    6, the input is invalid, because EPSABS < 0 and EPSREL < 0.
    !
    !  Local parameters:
    !
    !    LIMIT is the maximum number of subintervals allowed in
    !    the subdivision process of QAGE.
    !
    implicit none
  
    integer nth,nph,d1,d2,jj
    real(kind=8), dimension(d1,d2) :: f
    character*2 tf
    
    integer ( kind = 4 ), parameter :: limit = 500

    real ( kind = 4 ) a
    real ( kind = 4 ) abserr
    real ( kind = 4 ) alist(limit)
    real ( kind = 4 ) b
    real ( kind = 4 ) blist(limit)
    real ( kind = 4 ) elist(limit)
    real ( kind = 4 ) epsabs
    real ( kind = 4 ) epsrel
    integer ( kind = 4 ) ier
    integer ( kind = 4 ) iord(limit)
    integer ( kind = 4 ) key
    integer ( kind = 4 ) last
    integer ( kind = 4 ) neval
    real ( kind = 4 ) result
    real ( kind = 4 ) rlist(limit)

    call qage ( nth, nph, d1, d2, f, tf, jj, a, b, epsabs, epsrel, key, limit, result, abserr, neval, &
    ier, alist, blist, rlist, elist, iord, last )

    return
end
    
subroutine qage ( nth, nph, d1, d2, f, tf, jj, a, b, epsabs, epsrel, key, limit, result, abserr, neval, &
  ier, alist, blist, rlist, elist, iord, last )

!*****************************************************************************80
!
!! QAGE estimates a definite integral.
!
!  Discussion:
!
!    The routine calculates an approximation RESULT to a definite integral   
!      I = integral of F over (A,B),
!    hopefully satisfying
!      || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real ( kind = 4 ) F, the name of the function routine, of the form
!      function f ( x )
!      real ( kind = 4 ) f
!      real ( kind = 4 ) x
!    which evaluates the integrand function.
!
!    Input, real ( kind = 4 ) A, B, the limits of integration.
!
!    Input, real ( kind = 4 ) EPSABS, EPSREL, the absolute and relative accuracy requested.
!
!    Input, integer ( kind = 4 ) KEY, chooses the order of the local integration rule:
!    1,  7 Gauss points, 15 Gauss-Kronrod points,
!    2, 10 Gauss points, 21 Gauss-Kronrod points,
!    3, 15 Gauss points, 31 Gauss-Kronrod points,
!    4, 20 Gauss points, 41 Gauss-Kronrod points,
!    5, 25 Gauss points, 51 Gauss-Kronrod points,
!    6, 30 Gauss points, 61 Gauss-Kronrod points.
!
!    Input, integer ( kind = 4 ) LIMIT, the maximum number of subintervals that
!    can be used.
!
!    Output, real ( kind = 4 ) RESULT, the estimated value of the integral.
!
!    Output, real ( kind = 4 ) ABSERR, an estimate of || I - RESULT ||.
!
!    Output, integer ( kind = 4 ) NEVAL, the number of times the integral was evaluated.
!
!    Output, integer ( kind = 4 ) IER, return code.
!    0, normal and reliable termination of the routine.  It is assumed that the 
!      requested accuracy has been achieved.
!    1, maximum number of subdivisions allowed has been achieved.  One can 
!      allow more subdivisions by increasing the value of LIMIT in QAG. 
!      However, if this yields no improvement it is advised to analyze the
!      integrand to determine the integration difficulties.  If the position
!      of a local difficulty can be determined, such as a singularity or
!      discontinuity within the interval) one will probably gain from 
!      splitting up the interval at this point and calling the integrator 
!      on the subranges.  If possible, an appropriate special-purpose 
!      integrator should be used which is designed for handling the type 
!      of difficulty involved.
!    2, the occurrence of roundoff error is detected, which prevents the
!      requested tolerance from being achieved.
!    3, extremely bad integrand behavior occurs at some points of the
!      integration interval.
!    6, the input is invalid, because EPSABS < 0 and EPSREL < 0.
!
!    Workspace, real ( kind = 4 ) ALIST(LIMIT), BLIST(LIMIT), contains in entries 1 
!    through LAST the left and right ends of the partition subintervals.
!
!    Workspace, real ( kind = 4 ) RLIST(LIMIT), contains in entries 1 through LAST
!    the integral approximations on the subintervals.
!
!    Workspace, real ( kind = 4 ) ELIST(LIMIT), contains in entries 1 through LAST
!    the absolute error estimates on the subintervals.
!
!    Output, integer ( kind = 4 ) IORD(LIMIT), the first K elements of which are pointers 
!    to the error estimates over the subintervals, such that
!    elist(iord(1)), ..., elist(iord(k)) form a decreasing sequence, with
!    k = last if last <= (limit/2+2), and k = limit+1-last otherwise.
!
!    Output, integer ( kind = 4 ) LAST, the number of subintervals actually produced 
!    in the subdivision process.
!
!  Local parameters:
!
!    alist     - list of left end points of all subintervals
!                       considered up to now
!    blist     - list of right end points of all subintervals
!                       considered up to now
!    elist(i)  - error estimate applying to rlist(i)
!    maxerr    - pointer to the interval with largest error estimate
!    errmax    - elist(maxerr)
!    area      - sum of the integrals over the subintervals
!    errsum    - sum of the errors over the subintervals
!    errbnd    - requested accuracy max(epsabs,epsrel*abs(result))
!    *****1    - variable for the left subinterval
!    *****2    - variable for the right subinterval
!    last      - index for subdivision
!
    implicit none
  
    integer nth,nph,d1,d2,jj
    real(kind=8), dimension(d1,d2) :: f
    character*2 tf

    integer ( kind = 4 ) limit

    real ( kind = 4 ) a
    real ( kind = 4 ) abserr
    real ( kind = 4 ) alist(limit)
    real ( kind = 4 ) area
    real ( kind = 4 ) area1
    real ( kind = 4 ) area12
    real ( kind = 4 ) area2
    real ( kind = 4 ) a1
    real ( kind = 4 ) a2
    real ( kind = 4 ) b
    real ( kind = 4 ) blist(limit)
    real ( kind = 4 ) b1
    real ( kind = 4 ) b2
    real ( kind = 4 ) c
    real ( kind = 4 ) defabs
    real ( kind = 4 ) defab1
    real ( kind = 4 ) defab2
    real ( kind = 4 ) elist(limit)
    real ( kind = 4 ) epsabs
    real ( kind = 4 ) epsrel
    real ( kind = 4 ) errbnd
    real ( kind = 4 ) errmax
    real ( kind = 4 ) error1
    real ( kind = 4 ) error2
    real ( kind = 4 ) erro12
    real ( kind = 4 ) errsum
    integer ( kind = 4 ) ier
    integer ( kind = 4 ) iord(limit)
    integer ( kind = 4 ) iroff1
    integer ( kind = 4 ) iroff2
    integer ( kind = 4 ) key
    integer ( kind = 4 ) keyf
    integer ( kind = 4 ) last
    integer ( kind = 4 ) maxerr
    integer ( kind = 4 ) neval
    integer ( kind = 4 ) nrmax
    real ( kind = 4 ) resabs
    real ( kind = 4 ) result
    real ( kind = 4 ) rlist(limit)
!
!  Test on validity of parameters.
!
    ier = 0
    neval = 0
    last = 0
    result = 0.0E+00
    abserr = 0.0E+00
    alist(1) = a
    blist(1) = b
    rlist(1) = 0.0E+00
    elist(1) = 0.0E+00
    iord(1) = 0

    if ( epsabs < 0.0E+00 .and. epsrel < 0.0E+00 ) then
    ier = 6
    return
    end if
!
!  First approximation to the integral.
!
    keyf = key
    keyf = max ( keyf, 1 )
    keyf = min ( keyf, 6 )

    c = keyf
    neval = 0

    !if ( keyf == 1 ) then
    !  call qk15 ( f, a, b, result, abserr, defabs, resabs )
    !else if ( keyf == 2 ) then
    !  call qk21 ( f, a, b, result, abserr, defabs, resabs )
    !else if ( keyf == 3 ) then
    !  call qk31 ( f, a, b, result, abserr, defabs, resabs )
    !else if ( keyf == 4 ) then
    !  call qk41 ( f, a, b, result, abserr, defabs, resabs )
    !else if ( keyf == 5 ) then
    !  call qk51 ( f, a, b, result, abserr, defabs, resabs )
    !else if ( keyf == 6 ) then
    if ( keyf == 6 ) then
        call qk61 ( nth, nph, d1, d2, f, tf, jj, a, b, a, b, result, abserr, defabs, resabs )
    end if

    last = 1
    rlist(1) = result
    elist(1) = abserr
    iord(1) = 1
!
!  Test on accuracy.
!
    errbnd = max ( epsabs, epsrel * abs ( result ) )

    if ( abserr <= 5.0E+01 * epsilon ( defabs ) * defabs .and. &
    errbnd < abserr ) then
    ier = 2
    end if

    if ( limit == 1 ) then
    ier = 1
    end if

    if ( ier /= 0 .or. &
    ( abserr <= errbnd .and. abserr /= resabs ) .or. &
    abserr == 0.0E+00 ) then

    if ( keyf /= 1 ) then
        neval = (10*keyf+1) * (2*neval+1)
    else
        neval = 30 * neval + 15
    end if

    return

    end if
!
!  Initialization.
!
    errmax = abserr
    maxerr = 1
    area = result
    errsum = abserr
    nrmax = 1
    iroff1 = 0
    iroff2 = 0

    do last = 2, limit
!
!  Bisect the subinterval with the largest error estimate.
!
    a1 = alist(maxerr)
    b1 = 0.5E+00 * ( alist(maxerr) + blist(maxerr) )
    a2 = b1
    b2 = blist(maxerr)

    !if ( keyf == 1 ) then
    !  call qk15 ( f, a1, b1, area1, error1, resabs, defab1 )
    !else if ( keyf == 2 ) then
    !  call qk21 ( f, a1, b1, area1, error1, resabs, defab1 )
    !else if ( keyf == 3 ) then
    !  call qk31 ( f, a1, b1, area1, error1, resabs, defab1 )
    !else if ( keyf == 4 ) then
    !  call qk41 ( f, a1, b1, area1, error1, resabs, defab1)
    !else if ( keyf == 5 ) then
    !  call qk51 ( f, a1, b1, area1, error1, resabs, defab1 )
    !else if ( keyf == 6 ) then
    if ( keyf == 6 ) then
        call qk61 ( nth, nph, d1, d2, f, tf, jj, a, b, a1, b1, area1, error1, resabs, defab1 )
    end if

    !if ( keyf == 1 ) then
    !  call qk15 ( f, a2, b2, area2, error2, resabs, defab2 )
    !else if ( keyf == 2 ) then
    !  call qk21 ( f, a2, b2, area2, error2, resabs, defab2 )
    !else if ( keyf == 3 ) then
    !  call qk31 ( f, a2, b2, area2, error2, resabs, defab2 )
    !else if ( keyf == 4 ) then
    !  call qk41 ( f, a2, b2, area2, error2, resabs, defab2 )
    !else if ( keyf == 5 ) then
    !  call qk51 ( f, a2, b2, area2, error2, resabs, defab2 )
    !else if ( keyf == 6 ) then
    if ( keyf == 6 ) then
        call qk61 ( nth, nph, d1, d2, f, tf, jj, a, b, a2, b2, area2, error2, resabs, defab2 )
    end if
!
!  Improve previous approximations to integral and error and
!  test for accuracy.
!
    neval = neval + 1
    area12 = area1 + area2
    erro12 = error1 + error2
    errsum = errsum + erro12 - errmax
    area = area + area12 - rlist(maxerr)

    if ( defab1 /= error1 .and. defab2 /= error2 ) then

        if ( abs ( rlist(maxerr) - area12 ) <= 1.0E-05 * abs ( area12 ) &
        .and. 9.9E-01 * errmax <= erro12 ) then
        iroff1 = iroff1 + 1
        end if

        if ( 10 < last .and. errmax < erro12 ) then
        iroff2 = iroff2 + 1
        end if

    end if

    rlist(maxerr) = area1
    rlist(last) = area2
    errbnd = max ( epsabs, epsrel * abs ( area ) )
!
!  Test for roundoff error and eventually set error flag.
!
    if ( errbnd < errsum ) then

        if ( 6 <= iroff1 .or. 20 <= iroff2 ) then
        ier = 2
        end if
!
!  Set error flag in the case that the number of subintervals
!  equals limit.
!
        if ( last == limit ) then
        ier = 1
        end if
!
!  Set error flag in the case of bad integrand behavior
!  at a point of the integration range.
!
        if ( max ( abs ( a1 ), abs ( b2 ) ) <= ( 1.0E+00 + c * 1.0E+03 * &
        epsilon ( a1 ) ) * ( abs ( a2 ) + 1.0E+04 * tiny ( a2 ) ) ) then
        ier = 3
        end if

    end if
!
!  Append the newly-created intervals to the list.
!
    if ( error2 <= error1 ) then
        alist(last) = a2
        blist(maxerr) = b1
        blist(last) = b2
        elist(maxerr) = error1
        elist(last) = error2
    else
        alist(maxerr) = a2
        alist(last) = a1
        blist(last) = b1
        rlist(maxerr) = area2
        rlist(last) = area1
        elist(maxerr) = error2
        elist(last) = error1
    end if
!
!  Call QSORT to maintain the descending ordering
!  in the list of error estimates and select the subinterval
!  with the largest error estimate (to be bisected next).
!
    call qsort ( limit, last, maxerr, errmax, elist, iord, nrmax )
 
    if ( ier /= 0 .or. errsum <= errbnd ) then
        exit
    end if

    end do
!
!  Compute final result.
!
    result = sum ( rlist(1:last) )

    abserr = errsum

    if ( keyf /= 1 ) then
    neval = ( 10 * keyf + 1 ) * ( 2 * neval + 1 )
    else
    neval = 30 * neval + 15
    end if

    return
end    
    
subroutine qk61 ( nth, nph, d1, d2, f, tf, jj, a0, b0, a, b, result, abserr, resabs, resasc ) 

!*****************************************************************************80
!
!! QK61 carries out a 61 point Gauss-Kronrod quadrature rule.
!
!  Discussion:
!
!    This routine approximates
!      I = integral ( A <= X <= B ) F(X) dx
!    with an error estimate, and
!      J = integral ( A <= X <= B ) | F(X) | dx
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real ( kind = 4 ) F, the name of the function routine, of the form
!      function f ( x )
!      real ( kind = 4 ) f
!      real ( kind = 4 ) x
!    which evaluates the integrand function.
!
!    Input, real ( kind = 4 ) A, B, the limits of integration.
!
!    Output, real ( kind = 4 ) RESULT, the estimated value of the integral.
!                    result is computed by applying the 61-point
!                    Kronrod rule (resk) obtained by optimal addition of
!                    abscissae to the 30-point Gauss rule (resg).
!
!    Output, real ( kind = 4 ) ABSERR, an estimate of | I - RESULT |.
!
!    Output, real ( kind = 4 ) RESABS, approximation to the integral of the absolute
!    value of F.
!
!    Output, real ( kind = 4 ) RESASC, approximation to the integral | F-I/(B-A) | 
!    over [A,B].
!
!  Local Parameters:
!
!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc   - abscissa
!           fval*  - function value
!           resg   - result of the 30-point Gauss rule
!           resk   - result of the 61-point Kronrod rule
!           reskh  - approximation to the mean value of f
!                    over (a,b), i.e. to i/(b-a)
!
    USE Initialization
    implicit none
  
    integer nth,nph,d1,d2,jj
    real(kind=8), dimension(d1,d2) :: f
    character*2 tf
    integer iic,iicam,iicap

    real ( kind = 4 ) a,a0
    real ( kind = 4 ) absc
    real ( kind = 4 ) abserr
    real ( kind = 4 ) b,b0
    real ( kind = 4 ) centr
    real ( kind = 4 ) dhlgth
    real ( kind = 4 ) fc
    real ( kind = 4 ) fsum
    real ( kind = 4 ) fval1
    real ( kind = 4 ) fval2
    real ( kind = 4 ) fv1(30)
    real ( kind = 4 ) fv2(30)
    real ( kind = 4 ) hlgth
    integer ( kind = 4 ) jq   ! here just J is replaced by jq to avoid conflict with Initialization
    integer ( kind = 4 ) jtw
    integer ( kind = 4 ) jtwm1
    real ( kind = 4 ) resabs
    real ( kind = 4 ) resasc
    real ( kind = 4 ) resg
    real ( kind = 4 ) resk
    real ( kind = 4 ) reskh
    real ( kind = 4 ) result
    real ( kind = 4 ) wg(15)
    real ( kind = 4 ) wgk(31)
    real ( kind = 4 ) xgk(31)
!
!           the abscissae and weights are given for the
!           interval (-1,1). because of symmetry only the positive
!           abscissae and their corresponding weights are given.
!
!           xgk   - abscissae of the 61-point Kronrod rule
!                   xgk(2), xgk(4)  ... abscissae of the 30-point
!                   Gauss rule
!                   xgk(1), xgk(3)  ... optimally added abscissae
!                   to the 30-point Gauss rule
!
!           wgk   - weights of the 61-point Kronrod rule
!
!           wg    - weigths of the 30-point Gauss rule
!
    data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8), &
        xgk(9),xgk(10)/ &
        9.994844100504906E-01,     9.968934840746495E-01, &
        9.916309968704046E-01,     9.836681232797472E-01, &
        9.731163225011263E-01,     9.600218649683075E-01, &
        9.443744447485600E-01,     9.262000474292743E-01, &
        9.055733076999078E-01,     8.825605357920527E-01/
    data xgk(11),xgk(12),xgk(13),xgk(14),xgk(15),xgk(16),xgk(17), &
    xgk(18),xgk(19),xgk(20)/ &
        8.572052335460611E-01,     8.295657623827684E-01, &
        7.997278358218391E-01,     7.677774321048262E-01, &
        7.337900624532268E-01,     6.978504947933158E-01, &
        6.600610641266270E-01,     6.205261829892429E-01, &
        5.793452358263617E-01,     5.366241481420199E-01/
    data xgk(21),xgk(22),xgk(23),xgk(24),xgk(25),xgk(26),xgk(27), &
    xgk(28),xgk(29),xgk(30),xgk(31)/ &
        4.924804678617786E-01,     4.470337695380892E-01, &
        4.004012548303944E-01,     3.527047255308781E-01, &
        3.040732022736251E-01,     2.546369261678898E-01, &
        2.045251166823099E-01,     1.538699136085835E-01, &
        1.028069379667370E-01,     5.147184255531770E-02, &
        0.0E+00                   /
    data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8), &
    wgk(9),wgk(10)/ &
        1.389013698677008E-03,     3.890461127099884E-03, &
        6.630703915931292E-03,     9.273279659517763E-03, &
        1.182301525349634E-02,     1.436972950704580E-02, &
        1.692088918905327E-02,     1.941414119394238E-02, &
        2.182803582160919E-02,     2.419116207808060E-02/
    data wgk(11),wgk(12),wgk(13),wgk(14),wgk(15),wgk(16),wgk(17), &
    wgk(18),wgk(19),wgk(20)/ &
        2.650995488233310E-02,     2.875404876504129E-02, &
        3.090725756238776E-02,     3.298144705748373E-02, &
        3.497933802806002E-02,     3.688236465182123E-02, &
        3.867894562472759E-02,     4.037453895153596E-02, &
        4.196981021516425E-02,     4.345253970135607E-02/
    data wgk(21),wgk(22),wgk(23),wgk(24),wgk(25),wgk(26),wgk(27), &
    wgk(28),wgk(29),wgk(30),wgk(31)/ &
        4.481480013316266E-02,     4.605923827100699E-02, &
        4.718554656929915E-02,     4.818586175708713E-02, &
        4.905543455502978E-02,     4.979568342707421E-02, &
        5.040592140278235E-02,     5.088179589874961E-02, &
        5.122154784925877E-02,     5.142612853745903E-02, &
        5.149472942945157E-02/
    data wg(1),wg(2),wg(3),wg(4),wg(5),wg(6),wg(7),wg(8)/ &
        7.968192496166606E-03,     1.846646831109096E-02, &
        2.878470788332337E-02,     3.879919256962705E-02, &
        4.840267283059405E-02,     5.749315621761907E-02, &
        6.597422988218050E-02,     7.375597473770521E-02/
    data wg(9),wg(10),wg(11),wg(12),wg(13),wg(14),wg(15)/ &
        8.075589522942022E-02,     8.689978720108298E-02, &
        9.212252223778613E-02,     9.636873717464426E-02, &
        9.959342058679527E-02,     1.017623897484055E-01, &
        1.028526528935588E-01/

    centr = 5.0E-01*(b+a)
    hlgth = 5.0E-01*(b-a)
    dhlgth = abs(hlgth)
!
!  Compute the 61-point Kronrod approximation to the integral,
!  and estimate the absolute error.
!
    resg = 0.0E+00
    if (tf == "th") then 
        iic = nint(centr*(nth-1)/(b0-a0))+1;  ! ici je peux calculer iic sans avoir recours a center      
    else
        iic = nint(centr*(nph-1)/(b0-a0))+1;
    endif
    fc = f(iic, jj);
    !*****************************
    resk = wgk(31)*fc
    resabs = abs(resk)

    do jq = 1, 15
        jtw = jq*2
        absc = hlgth*xgk(jtw)
        if (tf == "th") then 
            iicam = nint((centr-absc)*(nth-1)/(b0-a0))+1;  
            iicap = nint((centr+absc)*(nth-1)/(b0-a0))+1;  
        else
            iicam = nint((centr-absc)*(nph-1)/(b0-a0))+1;
            iicap = nint((centr+absc)*(nph-1)/(b0-a0))+1;
        endif
        fval1 = f(iicam, jj);
        fval2 = f(iicap, jj)
        !**********************************************
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(jq)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
    end do

    do jq = 1, 15
        jtwm1 = jq*2-1
        absc = hlgth*xgk(jtwm1)
        if (tf == "th") then 
            iicam = nint((centr-absc)*(nth-1)/(b0-a0))+1;  
            iicap = nint((centr+absc)*(nth-1)/(b0-a0))+1;  
        else
            iicam = nint((centr-absc)*(nph-1)/(b0-a0))+1;
            iicap = nint((centr+absc)*(nph-1)/(b0-a0))+1;
        endif
        fval1 = f(iicam, jj);
        fval2 = f(iicap, jj)
        !**********************************************
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
    end do

    reskh = resk * 5.0E-01
    resasc = wgk(31)*abs(fc-reskh)

    do jq = 1, 30
    resasc = resasc+wgk(jq)*(abs(fv1(jq)-reskh)+abs(fv2(jq)-reskh))
    end do

    result = resk*hlgth
    resabs = resabs*dhlgth
    resasc = resasc*dhlgth
    abserr = abs((resk-resg)*hlgth)

    if ( resasc /= 0.0E+00 .and. abserr /= 0.0E+00) then
    abserr = resasc*min ( 1.0E+00,(2.0E+02*abserr/resasc)**1.5E+00)
    end if

    if ( resabs > tiny ( resabs ) / (5.0E+01* epsilon ( resabs ) )) then
    abserr = max ( ( epsilon ( resabs ) *5.0E+01)*resabs, abserr )
    end if

    return
end
    
subroutine qsort ( limit, last, maxerr, ermax, elist, iord, nrmax )

    !*****************************************************************************80
    !
    !! QSORT maintains the order of a list of local error estimates.
    !
    !  Discussion:
    !
    !    This routine maintains the descending ordering in the list of the 
    !    local error estimates resulting from the interval subdivision process. 
    !    At each call two error estimates are inserted using the sequential 
    !    search top-down for the largest error estimate and bottom-up for the
    !    smallest error estimate.
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) LIMIT, the maximum number of error estimates the list can
    !    contain.
    !
    !    Input, integer ( kind = 4 ) LAST, the current number of error estimates.
    !
    !    Input/output, integer ( kind = 4 ) MAXERR, the index in the list of the NRMAX-th 
    !    largest error.
    !
    !    Output, real ( kind = 4 ) ERMAX, the NRMAX-th largest error = ELIST(MAXERR).
    !
    !    Input, real ( kind = 4 ) ELIST(LIMIT), contains the error estimates.
    !
    !    Input/output, integer ( kind = 4 ) IORD(LAST).  The first K elements contain 
    !    pointers to the error estimates such that ELIST(IORD(1)) through
    !    ELIST(IORD(K)) form a decreasing sequence, with
    !      K = LAST 
    !    if 
    !      LAST <= (LIMIT/2+2), 
    !    and otherwise
    !      K = LIMIT+1-LAST.
    !
    !    Input/output, integer ( kind = 4 ) NRMAX.
    !
    implicit none

    integer ( kind = 4 ) last

    real ( kind = 4 ) elist(last)
    real ( kind = 4 ) ermax
    real ( kind = 4 ) errmax
    real ( kind = 4 ) errmin
    integer ( kind = 4 ) i
    integer ( kind = 4 ) ibeg
    integer ( kind = 4 ) iord(last)
    integer ( kind = 4 ) isucc
    integer ( kind = 4 ) j
    integer ( kind = 4 ) jbnd
    integer ( kind = 4 ) jupbn
    integer ( kind = 4 ) k
    integer ( kind = 4 ) limit
    integer ( kind = 4 ) maxerr
    integer ( kind = 4 ) nrmax
!
!  Check whether the list contains more than two error estimates.
!
    if ( last <= 2 ) then
    iord(1) = 1
    iord(2) = 2
    go to 90
    end if
!
!  This part of the routine is only executed if, due to a
!  difficult integrand, subdivision increased the error
!  estimate. in the normal case the insert procedure should
!  start after the nrmax-th largest error estimate.
!
    errmax = elist(maxerr)

    do i = 1, nrmax-1

    isucc = iord(nrmax-1)

    if ( errmax <= elist(isucc) ) then
        exit
    end if

    iord(nrmax) = isucc
    nrmax = nrmax-1

    end do
!
!  Compute the number of elements in the list to be maintained
!  in descending order.  This number depends on the number of
!  subdivisions still allowed.
!
    jupbn = last

    if ( (limit/2+2) < last ) then
    jupbn = limit+3-last
    end if

    errmin = elist(last)
!
!  Insert errmax by traversing the list top-down, starting
!  comparison from the element elist(iord(nrmax+1)).
!
    jbnd = jupbn-1
    ibeg = nrmax+1

    do i = ibeg, jbnd
    isucc = iord(i)
    if ( elist(isucc) <= errmax ) then
        go to 60
    end if
    iord(i-1) = isucc
    end do

    iord(jbnd) = maxerr
    iord(jupbn) = last
    go to 90
!
!  Insert errmin by traversing the list bottom-up.
!
60 continue

    iord(i-1) = maxerr
    k = jbnd

    do j = i, jbnd
    isucc = iord(k)
    if ( errmin < elist(isucc) ) then
        go to 80
    end if
    iord(k+1) = isucc
    k = k-1
    end do

    iord(i) = last
    go to 90

80 continue

    iord(k+1) = last
!
!  Set maxerr and ermax.
!
90 continue

    maxerr = iord(nrmax)
    ermax = elist(maxerr)

    return
end    
    
    