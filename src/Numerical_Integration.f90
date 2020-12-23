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