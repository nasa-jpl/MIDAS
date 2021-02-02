subroutine getGaussLegendreQuadPts1(Np,a,b,x,w)
    
        implicit none 
        
        ! in/out
        integer, intent(in) :: Np
        real(8), intent(in) :: a,b
        real(8), dimension(Np) :: x,w
        
        ! local
        integer :: i,p,k,N,N1,N2
        real(8) :: rp
        real(8), dimension(:,:), allocatable ::xu,y,L,Lp,t,y0
        real(8), parameter :: Pi = 4*atan(1.0_8);
        
        
        N = Np -1 ;
        N1=N+1; N2=N+2;
        
        allocate(xu(N1,1));
        do p=1,N1
            xu(p,1) = -1. + 2.*(p-1)/(N1-1);
        end do
        
        ! initial guess
        allocate(t(N1,1));t(:,1) = (/ (i,i=0,N) /);
        allocate(y(N1,1));
        y=cos((2*t+1)*Pi/(2*N+2))+(0.27/N1)*sin(Pi*xu*N/N2);
        
        ! legendre-Gauss Vandermonde Matrix
        Allocate(L(N1,N2)); L= 0; 
        ! Derivative of LGVM
        Allocate(Lp(N1,N2)); Lp= 0;
        
        ! compute the zeros of the N+1 Legendre Polynomial 
        ! using the recursion relation and the Newton-Raphson method
        allocate(y0(N1,1));
        y0 =2;
        
        ! Iterate until new points are uniformly within epsilon of old points
        do while ( maxval(abs(y-y0)) > EPSILON(y0))
            L(:,1) = 1.;
            Lp(:,1) = 0.;
            
            L(:,2)= y(:,1);
            Lp(:,2)=1.;
            
            do k=2,N1
                L(:,k+1)=( (2*k-1)*y(:,1)*L(:,k)-(k-1)*L(:,k-1) )/k;
            end do
            
            Lp(:,1)=N2*(L(:,N1)-y(:,1)*L(:,N2) )/(1-y(:,1)**2);   
            
            y0=y;
            y(:,1)=y0(:,1)-L(:,N2)/Lp(:,1);                 
        end do 
        
        !Linear map from[-1,1] to [a,b]
        x=(a*(1-y(:,1))+b*(1+y(:,1)))/2;     
        
        !Compute the weights
        rp = (real(N2)/real(N1))**2.;
        w=(b-a)/((1-y(:,1)**2.)*Lp(:,1)**2.)*rp;      
                
    end subroutine    
    
    subroutine getGaussLegendreQuadPts2(n,a,b,x,w)
        implicit none                                                                                
        
        ! in/out                                                                                     
        integer, intent(in) :: n                                                                     
        real(8), intent(in) :: a, b                                                                  
        real(8), dimension(n) :: xint,x, w                                                                
        ! local 
        real(8) :: alpha,beta,r;
        character*1 meth        
        meth = 'l';
        
        alpha=0.;beta=0;
        call cgqf(n,1,alpha,beta,a,b,xint,w)   
        x = xint(n:1:-1); 
        if (meth=='l') then
            r=1.1291;
            w = w*r;
        endif
    end subroutine                                                                                   
    
    subroutine getGaussLegendreQuadPts3(n,a,b,x,w)

        !*****************************************************************************80
        !
        !! MAIN (now getGaussLegendreQuadPts3) is the main program for LEGENDRE_RULE_FAST.
        !
        !  Discussion:
        !
        !    This program computes a standard Gauss-Legendre quadrature rule
        !    and writes it to a file.
        !
        !  Usage:
        !
        !    legendre_rule_fast ( n, a, b )
        !
        !    where
        !
        !    * n is the number of points in the rule;
        !    * a is the left endpoint;
        !    * b is the right endpoint.
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license. 
        !
        !  Modified:
        !
        !    28 June 2009
        !
        !  Author:
        !
        !    John Burkardt
        !
          implicit none
          
          ! in/out
          integer ( kind = 4 ), intent(in) :: n
          real ( kind = 8 ), intent(in) :: a,b
          real ( kind = 8 ), dimension (n), intent(out) :: x,w
          
          integer ( kind = 4 ) arg_num          
          integer ( kind = 4 ) iarg
          integer ( kind = 4 ) iargc
          integer ( kind = 4 ) ierror
          integer ( kind = 4 ) last
          character ( len = 255 ) string
          real ( kind = 8 ), dimension (n) :: xint                 
        !
        !  Construct the rule and output it.
        !
          call legendre_handle ( n, a, b , xint, w)
          x = xint(n:1:-1); 
        
    end subroutine
    
    subroutine getAdaptiveQuadPts1(Nth,Nph,dim1,dim2,fQscat,type_func,jj,a0,b0,I) 
    
        !based on "Adpative Quadrature" by Jim Lambers 2009
    
        USE stack_structure
        USE Initialization
    
        implicit none 
        
        ! IN/OUT
        Integer, INTENT(IN) :: Nth,Nph,dim1,dim2,jj
        Real(kind=8), Dimension(dim1,dim2), INTENT(IN) :: fQscat
        Character*2, INTENT(IN) :: type_func
        Real(kind=8), INTENT(IN) :: a0,b0 
        Real(kind=8), INTENT(OUT) :: I 
    
        ! LOCAL
        Real(kind=8) :: tol =1e-6;
        Integer :: N_div_min =2;
        Integer :: N_div_max =6;
        Integer iia, iib, iim,ctmp,ii;
        Integer :: niterdiv,cond_numb,nelts_max
        Real(kind=8) a,b,m,c,d,n, I1, I2
        
        niterdiv = 0; cond_numb = 0;
        nelts_max = 2**(N_div_max + 1);
        
        call push(a0); call push(b0); !! as each time, we're adding (push) ou deleting (pop) an interval, we always need to use twice push or pop ON THE OTHER WAY
                
        Do while (stack_top > 0)
            call pop(b); call pop(a); ! the interval [a,b] on top of S is recovered and removed from S 
            m = (a+b)/2.;
            If (type_func == "th") Then 
                iia = nint(a*(Nth-1)/Pi)+1;
                iib = nint(b*(Nth-1)/Pi)+1;
                iim = nint(m*(Nth-1)/Pi)+1;
            Else
                iia = nint(a*(Nph-1)/(2*Pi))+1;
                iib = nint(b*(Nph-1)/(2*Pi))+1;
                iim = nint(m*(Nph-1)/(2*Pi))+1;
            EndIf
            
            !! trapezoidal Rule : I1 = ((b - a) / 2) * (f(a) + f(b))
            I1 = ((b-a)/2.) * (fQscat(iia, jj) + fQscat(iib, jj));   
            !! composite trapezoidal rule with 2 subintervals ((b - a) / 4) * (f(a) + 2 * f(m) + f(b));
            I2 = ((b-a)/4.) * (fQscat(iia, jj) + 2 * fQscat(iim, jj) + fQscat(iib, jj));
        
            If ((cond_numb == 0) .AND. (stack_top >= nelts_max)) Then !once nelts_max achived we stop dividing and just compute the rest of intervals remaining in S
                cond_numb = 1;
            EndIf
        
            If ((cond_numb == 1) .OR. ((abs(I1-I2) < 3*(b-a)*tol) .AND. (niterdiv >= N_div_min))) Then ! from error term in trapezoidal Rule, |I(f) - I2| ~~ 1/3 |I1-I2| 
                I = I + I2; 
            Else
                ctmp = stack_top
                Do ii=1,ctmp,2
                    call pop(d); call pop(c)
                    n = (c+d)/2.
                    call push(n); call push(d); ! add the interval [m,b]
                    call push(c); call push(n);                    
                EndDo
                
                call push(m); call push(b); ! add the interval [m,b]
                call push(a); call push(m); ! add the interval [a,m] 
                niterdiv = niterdiv + 1;                             
            EndIf            
        EndDo
        
    end subroutine
    
    subroutine getAdaptiveQuadPts2(Nth,Nph,dim1,dim2,fQscat,type_func,jj,a0,b0,I)   
        
        !https://people.sc.fsu.edu/~jburkardt/f_src/quadpack/quadpack.html
        ! We use the subroutine qagp from QUADPACK 
        !(after modification : imput fxy instead of external function f)
        
        ! In/Out
        implicit none 
        
        ! IN/OUT
        Integer, INTENT(IN) :: Nth,Nph,dim1,dim2,jj
        Real(kind=8), Dimension(dim1,dim2), INTENT(IN) :: fQscat
        Character*2, INTENT(IN) :: type_func
        Real(kind=8), INTENT(IN) :: a0,b0 
        Real(kind=8), INTENT(OUT) :: I 
        
        ! Local 
        real ( kind = 4 ) abserr,a,b,rr
        real ( kind = 4 ), parameter :: epsabs = 0.0E+00
        real ( kind = 4 ), parameter :: epsrel = 0.001E+00
        real ( kind = 4 ), external :: f01
        integer ( kind = 4 ) ier,neval
        integer ( kind = 4), parameter :: key = 6 ! Key chooses the order of the integration, should be equal to 6 in 
                                                   ! our code since we only modified the quadpack subroutine qk61
         
        a=a0; b= b0;
        ! the quadpack qag is modified to fit our case, since f01 depends on jj, 
        ! we replace f01 by the 7 imputs [Nth,Nph,dim1,dim2,fQscat,type_func,jj]
        ! The change are mainly made in the quadpack subroutine qk61
        call qag ( Nth, Nph, dim1, dim2, fQscat, type_func, jj, a, b, epsabs, epsrel, key, rr, abserr, neval, ier)
        I =rr;
        
    end subroutine
    
    !! Integration of a function of two variables f(x,y)
    !! Integration of a function f(x,y) using cubature trapezoid rule (trap_2Dc.f90)
    !! http://ww2.odu.edu/~agodunov/computing/programs/
    
    subroutine trap_2Dc(f,a,b,c,d,integral,nx,ny)
    !==========================================================
    ! Integration of f(x,y) on [a,b] for x and [c,d] for y
    ! Method: cubature trapezoid rule for nx*ny points  
    ! written by: Alex Godunov (October 2009)
    !----------------------------------------------------------
    ! IN:
    ! f   - Function to integrate (supplied by a user)
    ! a	  - Lower limit of integration for x
    ! b	  - Upper limit of integration for x
    ! c	  - Lower limit of integration for y
    ! d	  - Upper limit of integration for y
    ! nx  - number of points along x
    ! ny  - number of points along y
    ! OUT:
    ! integral - Result of integration
    !==========================================================
    implicit none
    double precision f(nx,ny), a, b, c, d, integral, sum
    double precision hx, hy, x, y
    integer nx, ny
    integer i, j

    hx = (b-a)/dfloat(nx-1)
    hy = (d-c)/dfloat(ny-1)
   
    ! calculate the corner's terms
    !sum = f(a,c)+f(a,d)+f(b,c)+f(b,d)
    sum = f(1,1)+f(1,ny)+f(nx,1)+f(nx,ny)

    ! calculate single sums
    do i=2,nx-1
       !x = a + hx*(i-1)
       !sum = sum + 2.0*(f(x,c)+f(x,d))
        sum = sum + 2.0*(f(i,1)+f(i,ny))
    end do

    do j=2, ny-1
       !y = c + hy*(j-1)
       !sum = sum + 2.0*(f(a,y)+f(b,y))        
       sum = sum + 2.0*(f(1,j)+f(nx,j))
    end do

    ! calculate the double sum
    do i = 2,nx-1
       !x = a + hx*(i-1)
       do j = 2,ny-1
          !y = c + hy*(j-1)
          !sum = sum + 4.0*f(x,y)
          sum = sum + 4.0*f(i,j)
       end do
    end do
    integral = 0.25*hx*hy*sum
    return
    end subroutine trap_2Dc
    
    !! Automatic adaptive Integration of a function f(x,y) using Simpson rule (simpson2D.f90)
    !! http://ww2.odu.edu/~agodunov/computing/programs/
    recursive function simpson2D(f,a,b,eps)
    !==========================================================
    ! Integration of f(x) on [a,b]
    ! Method: Simpson rule with doubling number of intervals  
    !         till  error = coeff*|I_n - I_2n| < eps
    ! written by: Alex Godunov (October 2009)
    !----------------------------------------------------------
    ! IN:
    ! f   - Function to integrate (supplied by a user)
    ! a	  - Lower limit of integration
    ! b	  - Upper limit of integration
    ! eps - tolerance
    ! OUT:
    ! integral - Result of integration
    !==========================================================
    implicit none
    double precision f, a, b, eps, simpson2D
    double precision sn, s2n, h, x
    double precision, parameter :: coeff = 1.0/15.0 ! error estimate coeff
    integer, parameter :: nmax=1048576              ! max number of intervals
    integer n, i

    ! evaluate integral for 2 intervals (three points)
    h = (b-a)/2.0
    sn = (1.0/3.0)*h*(f(a)+4.0*f(a+h)+f(b))

    ! loop over number of intervals (starting from 4 intervals)
    n=4
    do while (n <= nmax)
       s2n = 0.0   
       h = (b-a)/dfloat(n)
       do i=2, n-2, 2
          x   = a+dfloat(i)*h
          s2n = s2n + 2.0*f(x) + 4.0*f(x+h)
       end do
       s2n = (s2n + f(a) + f(b) + 4.0*f(a+h))*h/3.0
       if(coeff*abs(s2n-sn) <= eps) then
          simpson2D = s2n
          exit
       end if
       sn = s2n
       n = n*2
    end do
    return
    end function simpson2D
    
    
        subroutine cgqf ( nt, kind, alpha, beta, a, b, t, wts )

        !*****************************************************************************80
        !
        !! CGQF computes knots and weights of a Gauss quadrature formula.
        !
        !  Discussion:
        !
        !    The user may specify the interval (A,B).
        !
        !    Only simple knots are produced.
        !
        !    Use routine EIQFS to evaluate this quadrature formula.
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license. 
        !
        !  Modified:
        !
        !    16 February 2010
        !
        !  Author:
        !
        !    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
        !    FORTRAN90 version by John Burkardt.
        !
        !  Reference:
        !
        !    Sylvan Elhay, Jaroslav Kautsky,
        !    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
        !    Interpolatory Quadrature,
        !    ACM Transactions on Mathematical Software,
        !    Volume 13, Number 4, December 1987, pages 399-415.
        !
        !  Parameters:
        !
        !    Input, integer ( kind = 4 ) NT, the number of knots.
        !
        !    Input, integer ( kind = 4 ) KIND, the rule.
        !    1, Legendre,             (a,b)       1.0
        !    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^-0.5)
        !    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
        !    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
        !    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
        !    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
        !    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
        !    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
        !    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
        !
        !    Input, real ( kind = 8 ) ALPHA, the value of Alpha, if needed.
        !
        !    Input, real ( kind = 8 ) BETA, the value of Beta, if needed.
        !
        !    Input, real ( kind = 8 ) A, B, the interval endpoints, or
        !    other parameters.
        !
        !    Output, real ( kind = 8 ) T(NT), the knots.
        !
        !    Output, real ( kind = 8 ) WTS(NT), the weights.
        !
          implicit none

          integer ( kind = 4 ) nt

          real ( kind = 8 ) a
          real ( kind = 8 ) alpha
          real ( kind = 8 ) b
          real ( kind = 8 ) beta
          integer ( kind = 4 ) i
          integer ( kind = 4 ) kind
          integer ( kind = 4 ), allocatable :: mlt(:)
          integer ( kind = 4 ), allocatable :: ndx(:)
          real ( kind = 8 ) t(nt)
          real ( kind = 8 ) wts(nt)
        !
        !  Compute the Gauss quadrature formula for default values of A and B.
        !
          call cdgqf ( nt, kind, alpha, beta, t, wts )
        !
        !  Prepare to scale the quadrature formula to other weight function with 
        !  valid A and B.
        !
          allocate ( mlt(1:nt) )

          mlt(1:nt) = 1

          allocate ( ndx(1:nt) )

          do i = 1, nt 
            ndx(i) = i
          end do

          call scqf ( nt, t, mlt, wts, nt, ndx, wts, t, kind, alpha, beta, a, b )

          deallocate ( mlt )
          deallocate ( ndx )

          return
    end       
    
    subroutine cdgqf ( nt, kind, alpha, beta, t, wts )

        !*****************************************************************************80
        !
        !! CDGQF computes a Gauss quadrature formula with default A, B and simple knots.
        !
        !  Discussion:
        !
        !    This routine computes all the knots and weights of a Gauss quadrature
        !    formula with a classical weight function with default values for A and B,
        !    and only simple knots.
        !
        !    There are no moments checks and no printing is done.
        !
        !    Use routine EIQFS to evaluate a quadrature computed by CGQFS.
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license. 
        !
        !  Modified:
        !
        !    04 January 2010
        !
        !  Author:
        !
        !    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
        !    FORTRAN90 version by John Burkardt.
        !
        !  Reference:
        !
        !    Sylvan Elhay, Jaroslav Kautsky,
        !    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
        !    Interpolatory Quadrature,
        !    ACM Transactions on Mathematical Software,
        !    Volume 13, Number 4, December 1987, pages 399-415.
        !
        !  Parameters:
        !
        !    Input, integer ( kind = 4 ) NT, the number of knots.
        !
        !    Input, integer ( kind = 4 ) KIND, the rule.
        !    1, Legendre,             (a,b)       1.0
        !    2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
        !    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
        !    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
        !    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
        !    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
        !    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
        !    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta
        !
        !    Input, real ( kind = 8 ) ALPHA, the value of Alpha, if needed.
        !
        !    Input, real ( kind = 8 ) BETA, the value of Beta, if needed.
        !
        !    Output, real ( kind = 8 ) T(NT), the knots.
        !
        !    Output, real ( kind = 8 ) WTS(NT), the weights.
        !
          implicit none

          integer ( kind = 4 ) nt

          real ( kind = 8 ) aj(nt)
          real ( kind = 8 ) alpha
          real ( kind = 8 ) beta
          real ( kind = 8 ) bj(nt)
          integer ( kind = 4 ) kind
          real ( kind = 8 ) t(nt)
          real ( kind = 8 ) wts(nt)
          real ( kind = 8 ) zemu

        !
        !  Get the Jacobi matrix and zero-th moment.
        !
          call class_matrix ( kind, nt, alpha, beta, aj, bj, zemu )
        !
        !  Compute the knots and weights.
        !
          call sgqf ( nt, aj, bj, zemu, t, wts )

          return
    end

    subroutine class_matrix ( kind, m, alpha, beta, aj, bj, zemu )

        !*****************************************************************************80
        !
        !! CLASS_MATRIX computes the Jacobi matrix for a quadrature rule.
        !
        !  Discussion:
        !
        !    This routine computes the diagonal AJ and sub-diagonal BJ
        !    elements of the order M tridiagonal symmetric Jacobi matrix
        !    associated with the polynomials orthogonal with respect to
        !    the weight function specified by KIND.
        !
        !    For weight functions 1-7, M elements are defined in BJ even
        !    though only M-1 are needed.  For weight function 8, BJ(M) is
        !    set to zero.
        !
        !    The zero-th moment of the weight function is returned in ZEMU.
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license. 
        !
        !  Modified:
        !
        !    27 December 2009
        !
        !  Author:
        !
        !    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
        !    FORTRAN90 version by John Burkardt.
        !
        !  Reference:
        !
        !    Sylvan Elhay, Jaroslav Kautsky,
        !    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
        !    Interpolatory Quadrature,
        !    ACM Transactions on Mathematical Software,
        !    Volume 13, Number 4, December 1987, pages 399-415.
        !
        !  Parameters:
        !
        !    Input, integer ( kind = 4 ) KIND, the rule.
        !    1, Legendre,             (a,b)       1.0
        !    2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
        !    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
        !    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
        !    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
        !    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
        !    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
        !    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta
        !
        !    Input, integer ( kind = 4 ) M, the order of the Jacobi matrix.
        !
        !    Input, real ( kind = 8 ) ALPHA, the value of Alpha, if needed.
        !
        !    Input, real ( kind = 8 ) BETA, the value of Beta, if needed.
        !
        !    Output, real ( kind = 8 ) AJ(M), BJ(M), the diagonal and subdiagonal
        !    of the Jacobi matrix.
        !
        !    Output, real ( kind = 8 ) ZEMU, the zero-th moment.
        !
          implicit none

          integer ( kind = 4 ) m

          real ( kind = 8 ) a2b2
          real ( kind = 8 ) ab
          real ( kind = 8 ) aba
          real ( kind = 8 ) abi
          real ( kind = 8 ) abj
          real ( kind = 8 ) abti
          real ( kind = 8 ) aj(m)
          real ( kind = 8 ) alpha
          real ( kind = 8 ) apone
          real ( kind = 8 ) beta
          real ( kind = 8 ) bj(m)
          integer ( kind = 4 ) i
          integer ( kind = 4 ) kind
          real ( kind = 8 ), parameter :: pi = 3.14159265358979323846264338327950D+00
          real ( kind = 8 ) r8_gamma
          real ( kind = 8 ) temp
          real ( kind = 8 ) temp2
          real ( kind = 8 ) zemu

          temp = epsilon ( temp )

          call parchk ( kind, 2 * m - 1, alpha, beta )

          temp2 = 0.5D+00

          !if ( 500.0D+00 * temp < abs ( ( r8_gamma ( temp2 ) )**2 - pi ) ) then
          !  write ( *, '(a)' ) ' '
          !  write ( *, '(a)' ) 'CLASS_MATRIX - Fatal error!'
          !  write ( *, '(a)' ) '  Gamma function does not match machine parameters.'
          !  stop
          !end if

          if ( kind == 1 ) then

            ab = 0.0D+00

            zemu = 2.0D+00 / ( ab + 1.0D+00 )

            aj(1:m) = 0.0D+00

            do i = 1, m
              abi = i + ab * mod ( i, 2 )
              abj = 2 * i + ab
              bj(i) = abi * abi / ( abj * abj - 1.0D+00 )
            end do
            bj(1:m) =  sqrt ( bj(1:m) )

          else if ( kind == 2 ) then

            zemu = pi

            aj(1:m) = 0.0D+00

            bj(1) =  sqrt ( 0.5D+00 )
            bj(2:m) = 0.5D+00

          else if ( kind == 3 ) then

            ab = alpha * 2.0D+00
            zemu = 2.0D+00**( ab + 1.0D+00 ) * r8_gamma ( alpha + 1.0D+00 )**2 &
              / r8_gamma ( ab + 2.0D+00 )

            aj(1:m) = 0.0D+00
            bj(1) = 1.0D+00 / ( 2.0D+00 * alpha + 3.0D+00 )
            do i = 2, m
              bj(i) = i * ( i + ab ) / ( 4.0D+00 * ( i + alpha )**2 - 1.0D+00 )
            end do
            bj(1:m) =  sqrt ( bj(1:m) )

          else if ( kind == 4 ) then

            ab = alpha + beta
            abi = 2.0D+00 + ab
            zemu = 2.0D+00**( ab + 1.0D+00 ) * r8_gamma ( alpha + 1.0D+00 ) &
              * r8_gamma ( beta + 1.0D+00 ) / r8_gamma ( abi )
            aj(1) = ( beta - alpha ) / abi
            bj(1) = 4.0D+00 * ( 1.0 + alpha ) * ( 1.0D+00 + beta ) &
              / ( ( abi + 1.0D+00 ) * abi * abi )
            a2b2 = beta * beta - alpha * alpha

            do i = 2, m
              abi = 2.0D+00 * i + ab
              aj(i) = a2b2 / ( ( abi - 2.0D+00 ) * abi )
              abi = abi**2
              bj(i) = 4.0D+00 * i * ( i + alpha ) * ( i + beta ) * ( i + ab ) &
                / ( ( abi - 1.0D+00 ) * abi )
            end do
            bj(1:m) =  sqrt ( bj(1:m) )

          else if ( kind == 5 ) then

            zemu = r8_gamma ( alpha + 1.0D+00 )

            do i = 1, m
              aj(i) = 2.0D+00 * i - 1.0D+00 + alpha
              bj(i) = i * ( i + alpha )
            end do
            bj(1:m) =  sqrt ( bj(1:m) )

          else if ( kind == 6 ) then

            zemu = r8_gamma ( ( alpha + 1.0D+00 ) / 2.0D+00 )

            aj(1:m) = 0.0D+00

            do i = 1, m
              bj(i) = ( i + alpha * mod ( i, 2 ) ) / 2.0D+00
            end do
            bj(1:m) =  sqrt ( bj(1:m) )

          else if ( kind == 7 ) then

            ab = alpha
            zemu = 2.0D+00 / ( ab + 1.0D+00 )

            aj(1:m) = 0.0D+00

            do i = 1, m
              abi = i + ab * mod ( i, 2 )
              abj = 2 * i + ab
              bj(i) = abi * abi / ( abj * abj - 1.0D+00 )
            end do
            bj(1:m) =  sqrt ( bj(1:m) )

          else if ( kind == 8 ) then

            ab = alpha + beta
            zemu = r8_gamma ( alpha + 1.0D+00 ) * r8_gamma ( - ( ab + 1.0D+00 ) ) &
              / r8_gamma ( - beta )
            apone = alpha + 1.0D+00
            aba = ab * apone
            aj(1) = - apone / ( ab + 2.0D+00 )
            bj(1) = - aj(1) * ( beta + 1.0D+00 ) / ( ab + 2.0D+00 ) / ( ab + 3.0D+00 )
            do i = 2, m
              abti = ab + 2.0D+00 * i
              aj(i) = aba + 2.0D+00 * ( ab + i ) * ( i - 1 )
              aj(i) = - aj(i) / abti / ( abti - 2.0D+00 )
            end do

            do i = 2, m - 1
              abti = ab + 2.0D+00 * i
              bj(i) = i * ( alpha + i ) / ( abti - 1.0D+00 ) * ( beta + i ) &
                / ( abti * abti ) * ( ab + i ) / ( abti + 1.0D+00 )
            end do

            bj(m) = 0.0D+00
            bj(1:m) =  sqrt ( bj(1:m) )

          end if

          return
    end
    
    subroutine sgqf ( nt, aj, bj, zemu, t, wts )

        !*****************************************************************************80
        !
        !! SGQF computes knots and weights of a Gauss Quadrature formula.
        !
        !  Discussion:
        !
        !    This routine computes all the knots and weights of a Gauss quadrature
        !    formula with simple knots from the Jacobi matrix and the zero-th
        !    moment of the weight function, using the Golub-Welsch technique.
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license. 
        !
        !  Modified:
        !
        !    04 January 2010
        !
        !  Author:
        !
        !    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
        !    FORTRAN90 version by John Burkardt.
        !
        !  Reference:
        !
        !    Sylvan Elhay, Jaroslav Kautsky,
        !    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
        !    Interpolatory Quadrature,
        !    ACM Transactions on Mathematical Software,
        !    Volume 13, Number 4, December 1987, pages 399-415.
        !
        !  Parameters:
        !
        !    Input, integer ( kind = 4 ) NT, the number of knots.
        !
        !    Input, real ( kind = 8 ) AJ(NT), the diagonal of the Jacobi matrix.
        !
        !    Input/output, real ( kind = 8 ) BJ(NT), the subdiagonal of the Jacobi 
        !    matrix, in entries 1 through NT-1.  On output, BJ has been overwritten.
        !
        !    Input, real ( kind = 8 ) ZEMU, the zero-th moment of the weight function.
        !
        !    Output, real ( kind = 8 ) T(NT), the knots.
        !
        !    Output, real ( kind = 8 ) WTS(NT), the weights.
        !
          implicit none

          integer ( kind = 4 ) nt,info

          real ( kind = 8 ) aj(nt)
          real ( kind = 8 ) bj(nt)
          integer ( kind = 4 ) i
          real ( kind = 8 ) t(nt)
          real ( kind = 8 ) wts(nt),z(nt,nt)
          real ( kind = 8 ) zemu
          double precision work(max(1,2*nt-2))
          real(8), parameter :: Pi = 4*atan(1.0_8);
        !
        !  Exit if the zero-th moment is not positive.
        !
          if ( zemu <= 0.0D+00 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'SGQF - Fatal error!'
            write ( *, '(a)' ) '  ZEMU <= 0.'
            stop
          end if
        !
        !  Set up vectors for IMTQLX.
        !
          t(1:nt) = aj(1:nt)                    
        !
        !  Diagonalize the Jacobi matrix. Change this
        !
          ! imtqlx
          !wts(1) = sqrt ( zemu )
          !wts(2:nt) = 0.0D+00
          !call imtqlx ( nt, t, bj, wts )
          !wts(1:nt) = wts(1:nt)**2
          
          ! OR lapack
          z=0.;
          z(1,1) = sqrt (zemu)
          forall (i=2:nt) z(i,i)=0.;
          call dsteqr('I',nt,t,bj(1:nt-1),z,nt,work,info);           
          !wts(1:nt) = z(1,:)**2.;
          wts(1:nt) = sqrt(Pi)*z(1,:)**2.;

          return
    end
    
    subroutine scqf ( nt, t, mlt, wts, nwts, ndx, swts, st, kind, alpha, beta, a, &
  b )

        !*****************************************************************************80
        !
        !! SCQF scales a quadrature formula to a nonstandard interval.
        !
        !  Discussion:
        !
        !    The arrays WTS and SWTS may coincide.
        !
        !    The arrays T and ST may coincide.
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license. 
        !
        !  Modified:
        !
        !    27 December 2009
        !
        !  Author:
        !
        !    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
        !    FORTRAN90 version by John Burkardt.
        !
        !  Reference:
        !
        !    Sylvan Elhay, Jaroslav Kautsky,
        !    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
        !    Interpolatory Quadrature,
        !    ACM Transactions on Mathematical Software,
        !    Volume 13, Number 4, December 1987, pages 399-415.
        !
        !  Parameters:
        !
        !    Input, integer ( kind = 4 ) NT, the number of knots.
        !
        !    Input, real ( kind = 8 ) T(NT), the original knots.
        !
        !    Input, integer ( kind = 4 ) MLT(NT), the multiplicity of the knots.
        !
        !    Input, real ( kind = 8 ) WTS(NWTS), the weights.
        !
        !    Input, integer ( kind = 4 ) NWTS, the number of weights.
        !
        !    Input, integer ( kind = 4 ) NDX(NT), used to index the array WTS.  
        !    For more details see the comments in CAWIQ.
        !
        !    Output, real ( kind = 8 ) SWTS(NWTS), the scaled weights.
        !
        !    Output, real ( kind = 8 ) ST(NT), the scaled knots.
        !
        !    Input, integer ( kind = 4 ) KIND, the rule.
        !    1, Legendre,             (a,b)       1.0
        !    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
        !    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
        !    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
        !    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
        !    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
        !    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
        !    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
        !    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
        !
        !    Input, real ( kind = 8 ) ALPHA, the value of Alpha, if needed.
        !
        !    Input, real ( kind = 8 ) BETA, the value of Beta, if needed.
        !
        !    Input, real ( kind = 8 ) A, B, the interval endpoints.
        !
          implicit none

          integer ( kind = 4 ) nt
          integer ( kind = 4 ) nwts

          real ( kind = 8 ) a
          real ( kind = 8 ) al
          real ( kind = 8 ) alpha
          real ( kind = 8 ) b
          real ( kind = 8 ) be
          real ( kind = 8 ) beta
          integer ( kind = 4 ) i
          integer ( kind = 4 ) k
          integer ( kind = 4 ) kind
          integer ( kind = 4 ) l
          integer ( kind = 4 ) mlt(nt)
          integer ( kind = 4 ) ndx(nt)
          real ( kind = 8 ) p
          real ( kind = 8 ) shft
          real ( kind = 8 ) slp
          real ( kind = 8 ) st(nt)
          real ( kind = 8 ) swts(nwts)
          real ( kind = 8 ) t(nt)
          real ( kind = 8 ) temp
          real ( kind = 8 ) tmp
          real ( kind = 8 ) wts(nwts)

          temp = epsilon ( temp )

          call parchk ( kind, 1, alpha, beta )

          if ( kind == 1 ) then

            al = 0.0D+00
            be = 0.0D+00

            if ( abs ( b - a ) <= temp ) then
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) 'SCQF - Fatal error!'
              write ( *, '(a)' ) '  |B - A| too small.'
              stop
            end if

            shft = ( a + b ) / 2.0D+00
            slp = ( b - a ) / 2.0D+00

          else if ( kind == 2 ) then

            al = -0.5D+00
            be = -0.5D+00

            if ( abs ( b - a ) <= temp ) then
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) 'SCQF - Fatal error!'
              write ( *, '(a)' ) '  |B - A| too small.'
              stop
            end if

            shft = ( a + b ) / 2.0D+00
            slp = ( b - a ) / 2.0D+00

          else if ( kind == 3 ) then

            al = alpha
            be = alpha

            if ( abs ( b - a ) <= temp ) then
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) 'SCQF - Fatal error!'
              write ( *, '(a)' ) '  |B - A| too small.'
              stop
            end if

            shft = ( a + b ) / 2.0D+00
            slp = ( b - a ) / 2.0D+00

          else if ( kind == 4 ) then

            al = alpha
            be = beta

            if ( abs ( b - a ) <= temp ) then
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) 'SCQF - Fatal error!'
              write ( *, '(a)' ) '  |B - A| too small.'
              stop
            end if

            shft = ( a + b ) / 2.0D+00
            slp = ( b - a ) / 2.0D+00

          else if ( kind == 5 ) then

            if ( b <= 0.0D+00 ) then
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) 'SCQF - Fatal error!'
              write ( *, '(a)' ) '  B <= 0'
              stop
            end if

            shft = a
            slp = 1.0D+00 / b
            al = alpha
            be = 0.0D+00

          else if ( kind == 6 ) then

            if ( b <= 0.0D+00 ) then
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) 'SCQF - Fatal error!'
              write ( *, '(a)' ) '  B <= 0.'
              stop
            end if

            shft = a
            slp = 1.0D+00 / sqrt ( b )
            al = alpha
            be = 0.0D+00

          else if ( kind == 7 ) then

            al = alpha
            be = 0.0D+00

            if ( abs ( b - a ) <= temp ) then
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) 'SCQF - Fatal error!'
              write ( *, '(a)' ) '  |B - A| too small.'
              stop
            end if

            shft = ( a + b ) / 2.0D+00
            slp = ( b - a ) / 2.0D+00

          else if ( kind == 8 ) then

            if ( a + b <= 0.0D+00 ) then
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) 'SCQF - Fatal error!'
              write ( *, '(a)' ) '  A + B <= 0.'
              stop
            end if

            shft = a
            slp = a + b
            al = alpha
            be = beta

          else if ( kind == 9 ) then

            al = 0.5D+00
            be = 0.5D+00

            if ( abs ( b - a ) <= temp ) then
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) 'SCQF - Fatal error!'
              write ( *, '(a)' ) '  |B - A| too small.'
              stop
            end if

            shft = ( a + b ) / 2.0D+00
            slp = ( b - a ) / 2.0D+00

          end if

          p = slp**( al + be + 1.0D+00 )

          do k = 1, nt

            st(k) = shft + slp * t(k)
            l = abs ( ndx(k) )

            if ( l /= 0 ) then
              tmp = p
              do i = l, l + mlt(k) - 1
                swts(i) = wts(i) * tmp
                tmp = tmp * slp
              end do
            end if

          end do

          return
    end
    
    subroutine imtqlx ( n, d, e, z )

    !*****************************************************************************80
    !
    !! IMTQLX diagonalizes a symmetric tridiagonal matrix.
    !
    !  Discussion:
    !
    !    This routine is a slightly modified version of the EISPACK routine to
    !    perform the implicit QL algorithm on a symmetric tridiagonal matrix.
    !
    !    The authors thank the authors of EISPACK for permission to use this
    !    routine.
    !
    !    It has been modified to produce the product Q' * Z, where Z is an input
    !    vector and Q is the orthogonal matrix diagonalizing the input matrix.
    !    The changes consist (essentially) of applying the orthogonal 
    !    transformations directly to Z as they are generated.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    27 December 2009
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Sylvan Elhay, Jaroslav Kautsky,
    !    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
    !    Interpolatory Quadrature,
    !    ACM Transactions on Mathematical Software,
    !    Volume 13, Number 4, December 1987, pages 399-415.
    !
    !    Roger Martin, James Wilkinson,
    !    The Implicit QL Algorithm,
    !    Numerische Mathematik,
    !    Volume 12, Number 5, December 1968, pages 377-383.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) N, the order of the matrix.
    !
    !    Input/output, real ( kind = 8 ) D(N), the diagonal entries of the matrix.
    !    On output, the information in D has been overwritten.
    !
    !    Input/output, real ( kind = 8 ) E(N), the subdiagonal entries of the
    !    matrix, in entries E(1) through E(N-1).  On output, the information in
    !    E has been overwritten.
    !
    !    Input/output, real ( kind = 8 ) Z(N).  On input, a vector.  On output,
    !    the value of Q' * Z, where Q is the matrix that diagonalizes the
    !    input symmetric tridiagonal matrix.
    !
      implicit none

      integer ( kind = 4 ) n

      real ( kind = 8 ) b
      real ( kind = 8 ) c
      real ( kind = 8 ) d(n)
      real ( kind = 8 ) e(n)
      real ( kind = 8 ) f
      real ( kind = 8 ) g
      integer ( kind = 4 ) i
      integer ( kind = 4 ) ii
      integer ( kind = 4 ), parameter :: itn = 30
      integer ( kind = 4 ) j
      integer ( kind = 4 ) k
      integer ( kind = 4 ) l
      integer ( kind = 4 ) m
      integer ( kind = 4 ) mml
      real ( kind = 8 ) p
      real ( kind = 8 ) prec
      real ( kind = 8 ) r
      real ( kind = 8 ) s
      real ( kind = 8 ) z(n)

      prec = epsilon ( prec )

      if ( n == 1 ) then
        return
      end if

      e(n) = 0.0D+00

      do l = 1, n

        j = 0

        do

          do m = l, n

            if ( m == n ) then
              exit
            end if

            if ( abs ( e(m) ) <= prec * ( abs ( d(m) ) + abs ( d(m+1) ) ) ) then
              exit
            end if

          end do

          p = d(l)

          if ( m == l ) then
            exit
          end if

          if ( itn <= j ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'IMTQLX - Fatal error!'
            write ( *, '(a)' ) '  Iteration limit exceeded.'
            write ( *, '(a,i8)' ) '  J = ', j
            write ( *, '(a,i8)' ) '  L = ', l
            write ( *, '(a,i8)' ) '  M = ', m
            write ( *, '(a,i8)' ) '  N = ', n
            stop
          end if

          j = j + 1
          g = ( d(l+1) - p ) / ( 2.0D+00 * e(l) )
          r =  sqrt ( g * g + 1.0D+00 )
          g = d(m) - p + e(l) / ( g + sign ( r, g ) )
          s = 1.0D+00
          c = 1.0D+00
          p = 0.0D+00
          mml = m - l

          do ii = 1, mml

            i = m - ii
            f = s * e(i)
            b = c * e(i)

            if ( abs ( g ) <= abs ( f ) ) then
              c = g / f
              r =  sqrt ( c * c + 1.0D+00 )
              e(i+1) = f * r
              s = 1.0D+00 / r
              c = c * s
            else
              s = f / g
              r =  sqrt ( s * s + 1.0D+00 )
              e(i+1) = g * r
              c = 1.0D+00 / r
              s = s * c
            end if

            g = d(i+1) - p
            r = ( d(i) - g ) * s + 2.0D+00 * c * b
            p = s * r
            d(i+1) = g + p
            g = c * r - b
            f = z(i+1)
            z(i+1) = s * z(i) + c * f
            z(i) = c * z(i) - s * f

          end do

          d(l) = d(l) - p
          e(l) = g
          e(m) = 0.0D+00

        end do

      end do
    !
    !  Sorting.
    !
      do ii = 2, n

        i = ii - 1
        k = i
        p = d(i)

        do j = ii, n
          if ( d(j) < p ) then
            k = j
            p = d(j)
          end if
        end do

        if ( k /= i ) then
          d(k) = d(i)
          d(i) = p
          p = z(i)
          z(i) = z(k)
          z(k) = p
        end if

      end do

      return
    end subroutine
    
    subroutine parchk ( kind, m, alpha, beta )
    
        !*****************************************************************************80
        !
        !! PARCHK checks parameters ALPHA and BETA for classical weight functions. 
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license. 
        !
        !  Modified:
        !
        !    27 December 2009
        !
        !  Author:
        !
        !    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
        !    FORTRAN90 version by John Burkardt.
        !
        !  Reference:
        !
        !    Sylvan Elhay, Jaroslav Kautsky,
        !    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
        !    Interpolatory Quadrature,
        !    ACM Transactions on Mathematical Software,
        !    Volume 13, Number 4, December 1987, pages 399-415.
        !
        !  Parameters:
        !
        !    Input, integer ( kind = 4 ) KIND, the rule.
        !    1, Legendre,             (a,b)       1.0
        !    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
        !    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
        !    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
        !    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
        !    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
        !    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
        !    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta
        !    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
        !
        !    Input, integer ( kind = 4 ) M, the order of the highest moment to
        !    be calculated.  This value is only needed when KIND = 8.
        !
        !    Input, real ( kind = 8 ) ALPHA, BETA, the parameters, if required
        !    by the value of KIND.
        !
          implicit none
    
          real ( kind = 8 ) alpha
          real ( kind = 8 ) beta
          integer ( kind = 4 ) kind
          integer ( kind = 4 ) m
          real ( kind = 8 ) tmp
    
          if ( kind <= 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'PARCHK - Fatal error!'
            write ( *, '(a)' ) '  KIND <= 0.'
            stop
          end if
        !
        !  Check ALPHA for Gegenbauer, Jacobi, Laguerre, Hermite, Exponential.
        !
          if ( 3 <= kind .and. kind <= 8 .and. alpha <= -1.0D+00 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'PARCHK - Fatal error!'
            write ( *, '(a)' ) '  3 <= KIND and ALPHA <= -1.'
            stop
          end if
        !
        !  Check BETA for Jacobi.
        !
          if ( kind == 4 .and. kind <= 8 .and. beta <= -1.0D+00 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'PARCHK - Fatal error!'
            write ( *, '(a)' ) '  KIND == 4 and BETA <= -1.0.'
            stop
          end if
        !
        !  Check ALPHA and BETA for rational.
        !
          if ( kind == 8 ) then
            tmp = alpha + beta + m + 1.0D+00
            if ( 0.0D+00 <= tmp .or. tmp <= beta ) then
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) 'PARCHK - Fatal error!'
              write ( *, '(a)' ) '  KIND == 8 but condition on ALPHA and BETA fails.'
              stop
            end if
          end if
    
          return
    end
    
    function r8_epsilon ( )

    !*****************************************************************************80
    !
    !! R8_EPSILON returns the R8 roundoff unit.
    !
    !  Discussion:
    !
    !    The roundoff unit is a number R which is a power of 2 with the
    !    property that, to the precision of the computer's arithmetic,
    !      1 < 1 + R
    !    but
    !      1 = ( 1 + R / 2 )
    !
    !    FORTRAN90 provides the superior library routine
    !
    !      EPSILON ( X )
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    01 September 2012
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Output, real ( kind = 8 ) R8_EPSILON, the round-off unit.
    !
      implicit none

      real ( kind = 8 ) r8_epsilon

      r8_epsilon = 2.220446049250313D-016

      return
    end
    
    function r8_gamma ( x )

        !*****************************************************************************80
        !
        !! R8_GAMMA evaluates Gamma(X) for a real argument.
        !
        !  Discussion:
        !
        !    This routine calculates the gamma function for a real argument X.
        !
        !    Computation is based on an algorithm outlined in reference 1.
        !    The program uses rational functions that approximate the gamma
        !    function to at least 20 significant decimal digits.  Coefficients
        !    for the approximation over the interval (1,2) are unpublished.
        !    Those for the approximation for 12 <= X are from reference 2.
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license. 
        !
        !  Modified:
        !
        !    15 April 2013
        !
        !  Author:
        !
        !    Original FORTRAN77 version by William Cody, Laura Stoltz.
        !    FORTRAN90 version by John Burkardt.
        !
        !  Reference:
        !
        !    William Cody,
        !    An Overview of Software Development for Special Functions,
        !    in Numerical Analysis Dundee, 1975,
        !    edited by GA Watson,
        !    Lecture Notes in Mathematics 506,
        !    Springer, 1976.
        !
        !    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
        !    Charles Mesztenyi, John Rice, Henry Thatcher,
        !    Christoph Witzgall,
        !    Computer Approximations,
        !    Wiley, 1968,
        !    LC: QA297.C64.
        !
        !  Parameters:
        !
        !    Input, real ( kind = 8 ) X, the argument of the function.
        !
        !    Output, real ( kind = 8 ) R8_GAMMA, the value of the function.
        !
          implicit none

          real ( kind = 8 ), dimension ( 7 ) :: c = (/ &
           -1.910444077728D-03, &
            8.4171387781295D-04, &
           -5.952379913043012D-04, &
            7.93650793500350248D-04, &
           -2.777777777777681622553D-03, &
            8.333333333333333331554247D-02, &
            5.7083835261D-03 /)
          real ( kind = 8 ) fact
          integer ( kind = 4 ) i
          integer ( kind = 4 ) n
          real ( kind = 8 ), dimension ( 8 ) :: p = (/ &
            -1.71618513886549492533811D+00, &
             2.47656508055759199108314D+01, &
            -3.79804256470945635097577D+02, &
             6.29331155312818442661052D+02, &
             8.66966202790413211295064D+02, &
            -3.14512729688483675254357D+04, &
            -3.61444134186911729807069D+04, &
             6.64561438202405440627855D+04 /)
          logical parity
          real ( kind = 8 ), parameter :: pi = 3.1415926535897932384626434D+00
          real ( kind = 8 ), dimension ( 8 ) :: q = (/ &
            -3.08402300119738975254353D+01, &
             3.15350626979604161529144D+02, &
            -1.01515636749021914166146D+03, &
            -3.10777167157231109440444D+03, &
             2.25381184209801510330112D+04, &
             4.75584627752788110767815D+03, &
            -1.34659959864969306392456D+05, &
            -1.15132259675553483497211D+05 /)
          real ( kind = 8 ) r8_epsilon
          real ( kind = 8 ) r8_gamma
          real ( kind = 8 ) res
          real ( kind = 8 ), parameter :: sqrtpi = 0.9189385332046727417803297D+00
          real ( kind = 8 ) sum
          real ( kind = 8 ) x
          real ( kind = 8 ), parameter :: xbig = 171.624D+00
          real ( kind = 8 ) xden
          real ( kind = 8 ), parameter :: xinf = 1.79D+308
          real ( kind = 8 ), parameter :: xminin = 2.23D-308
          real ( kind = 8 ) xnum
          real ( kind = 8 ) y
          real ( kind = 8 ) y1
          real ( kind = 8 ) ysq
          real ( kind = 8 ) z

          parity = .false.
          fact = 1.0D+00
          n = 0
          y = x
        !
        !  Argument is negative.
        !
          if ( y <= 0.0D+00 ) then

            y = - x
            y1 = aint ( y )
            res = y - y1

            if ( res /= 0.0D+00 ) then

              if ( y1 /= aint ( y1 * 0.5D+00 ) * 2.0D+00 ) then
                parity = .true.
              end if

              fact = - pi / sin ( pi * res )
              y = y + 1.0D+00

            else

              res = xinf
              r8_gamma = res
              return

            end if

          end if
        !
        !  Argument is positive.
        !
          if ( y < r8_epsilon ( ) ) then
        !
        !  Argument < EPS.
        !
            if ( xminin <= y ) then
              res = 1.0D+00 / y
            else
              res = xinf
              r8_gamma = res
              return
            end if

          else if ( y < 12.0D+00 ) then

            y1 = y
        !
        !  0.0 < argument < 1.0.
        !
            if ( y < 1.0D+00 ) then

              z = y
              y = y + 1.0D+00
        !
        !  1.0 < argument < 12.0.
        !  Reduce argument if necessary.
        !
            else

              n = int ( y ) - 1
              y = y - real ( n, kind = 8 )
              z = y - 1.0D+00

            end if
        !
        !  Evaluate approximation for 1.0 < argument < 2.0.
        !
            xnum = 0.0D+00
            xden = 1.0D+00
            do i = 1, 8
              xnum = ( xnum + p(i) ) * z
              xden = xden * z + q(i)
            end do

            res = xnum / xden + 1.0D+00
        !
        !  Adjust result for case  0.0 < argument < 1.0.
        !
            if ( y1 < y ) then

              res = res / y1
        !
        !  Adjust result for case 2.0 < argument < 12.0.
        !
            else if ( y < y1 ) then

              do i = 1, n
                res = res * y
                y = y + 1.0D+00
              end do

            end if

          else
        !
        !  Evaluate for 12.0 <= argument.
        !
            if ( y <= xbig ) then

              ysq = y * y
              sum = c(7)
              do i = 1, 6
                sum = sum / ysq + c(i)
              end do
              sum = sum / y - y + sqrtpi
              sum = sum + ( y - 0.5D+00 ) * log ( y )
              res = exp ( sum )

            else

              res = xinf
              r8_gamma = res
              return

            end if

          end if
        !
        !  Final adjustments and return.
        !
          if ( parity ) then
            res = - res
          end if

          if ( fact /= 1.0D+00 ) then
            res = fact / res
          end if

          r8_gamma = res

          return
    end function
    

    
