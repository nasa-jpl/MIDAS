MODULE DiverseUtil
    
    implicit none
    
CONTAINS
    SUBROUTINE unique_sort(N,val,val_final)
        implicit none
        
        integer, INTENT(IN) :: N
        integer, dimension(N), INTENT(IN) :: val
        integer, dimension(:), allocatable, INTENT(OUT):: val_final
        
        ! local
        integer :: i = 0, min_val, max_val
        integer, dimension(N) :: unique
        
        min_val = minval(val)-1
        max_val = maxval(val)
        
        do while (min_val<max_val)
            i = i+1
            min_val = minval(val, mask=val>min_val)
            unique(i) = min_val
        enddo
        allocate(val_final(i), source=unique(1:i))   !<-- Or, just use unique(1:i) 
        
    end subroutine unique_sort
    
    SUBROUTINE count_unique_vals(N,val,count)
        implicit none
        
        integer, INTENT(IN) :: N
        integer, dimension(N), INTENT(IN) :: val
        integer, INTENT(OUT):: count
        
        ! local
        integer :: ii, prev_val
        integer, dimension(N) :: unique
        
        prev_val = val(1);
        count = 1;
        do ii=2,N 
            if (val(ii) .ne. prev_val) then 
                count = count + 1;
            endif            
        enddo       
    end subroutine count_unique_vals
    
    SUBROUTINE RandomReal(real_a,real_b, real_al)
 
        !! Generate a random real(kind=8) between real_a and real_b
        USE Initialization
        Implicit none
    
        !IN/OUT
        Real(kind=8), INTENT(IN) :: real_a
        Real(kind=8), INTENT(IN) :: real_b
        Real(kind=8), INTENT(OUT) :: real_al

        !Local
        Real(kind=8) :: a
    
        
        call DATE_AND_TIME
        call random_number(harvest=a)
        real_al = (a * (real_b - real_a)) + real_a

        real_al = aint(10**(Round_D+1)*real_al)/10**(Round_D+1)

    End Subroutine RandomReal

    SUBROUTINE RandomInteger(int_a,int_b, int_al)
 
        !! Generate a random number between the two integers int_a and int_b
        USE Initialization
        Implicit none

        integer, INTENT(IN) :: int_a
        integer, INTENT(IN) :: int_b
        integer, INTENT(OUT) :: int_al
        Real(kind=8) :: a


    
        ! DATE_AND_TIME ensure more randomness 
        call DATE_AND_TIME
        call random_number(harvest=a)
        int_al = INT(a * (int_b + 1 - int_a)) + int_a

    End Subroutine RandomInteger
    
    Subroutine rand_normal(mean,stdev,c)
    
        USE Initialization
        Implicit none
    
        real(kind=8) :: r,theta
        real(kind=8), Dimension(2) :: temp
    
        real(kind=8), INTENT(IN) :: mean,stdev
        real(kind=8), INTENT(OUT) :: c
    
        IF(stdev <= 0.0d0) THEN
            WRITE(*,*) "Standard Deviation must be +ve"
        ELSE
            CALL RANDOM_NUMBER(temp)
            r=(-2.0d0*log(temp(1)))**0.5
            theta = 2.0d0*PI*temp(2)
            c= mean+stdev*r*sin(theta)
        END IF    
    
    End Subroutine rand_normal
    
    ! --------------------------------------------------------------------
    ! SUBROUTINE  Sort_desc():
    !    This subroutine receives an array x() and sorts it into descending
    ! order.
    ! --------------------------------------------------------------------

   SUBROUTINE  Sort_desc(x, x_size,order)
      IMPLICIT  NONE
      
      INTEGER, INTENT(IN) :: x_size
      INTEGER, DIMENSION(x_size), INTENT(INOUT) :: x
      INTEGER, DIMENSION(x_size), INTENT(INOUT) :: order
      
      INTEGER :: i
      INTEGER :: Location

      DO i = 1, x_size-1			! except for the last
         CALL FindMaximum(i, x_size, x, Location)	! find max from this to last
         CALL  Swap(x(i), x(Location))	! swap this and the maximum
         CALL  Swap(order(i), order(Location)); ! backup the swap action        
      END DO
   END SUBROUTINE  Sort_desc
   
    ! --------------------------------------------------------------------
    ! SUBROUTINE  Sort_asc():
    !    This subroutine receives an array x() and sorts it into ascending
    ! order.
    ! --------------------------------------------------------------------

   SUBROUTINE  Sort_asc(x, x_size,order)
      IMPLICIT  NONE
      
      INTEGER, INTENT(IN) :: x_size
      INTEGER, DIMENSION(x_size), INTENT(INOUT) :: x
      INTEGER, DIMENSION(x_size), INTENT(INOUT) :: order
      
      INTEGER :: i
      INTEGER :: Location

      DO i = 1, x_size-1			! except for the last
         CALL FindMinimum(i, x_size, x, Location)	! find min from this to last
         CALL  Swap(x(i), x(Location))	! swap this and the min
         CALL  Swap(order(i), order(Location)); ! backup the swap action        
      END DO
   END SUBROUTINE  Sort_asc



    ! --------------------------------------------------------------------
    ! INTEGER FUNCTION  FindMinimum():
    !    This function returns the location of the minimum in the section
    ! between Start and End.
    ! --------------------------------------------------------------------

   SUBROUTINE FindMinimum(xStart, xEnd, x, MinLoc)
      IMPLICIT  NONE
      
      INTEGER, INTENT(IN) :: xStart, xEnd
      INTEGER, DIMENSION(xEnd), INTENT(IN) :: x
      INTEGER, INTENT(OUT) :: MinLoc
      INTEGER :: Minimum
      INTEGER :: Location
      INTEGER :: i

      Minimum  = x(xStart)		! assume the first is the min
      Location = xStart			! record its position
      DO i = xStart+1, xEnd		! start with next elements
         IF (x(i) < Minimum) THEN	!   if x(i) less than the min?
            Minimum  = x(i)		!      Yes, a new minimum found
            Location = i                !      record its position
         END IF
      END DO
      MinLoc = Location        	! return the position
   END SUBROUTINE FindMinimum
   
    ! --------------------------------------------------------------------
    ! INTEGER FUNCTION  FindMaximum():
    !    This function returns the location of the maximum in the section
    ! between Start and End.
    ! --------------------------------------------------------------------

   SUBROUTINE FindMaximum(xStart, xEnd, x, MaxLoc)
      IMPLICIT  NONE
      
      INTEGER, INTENT(IN) :: xStart, xEnd
      INTEGER, DIMENSION(xEnd), INTENT(IN) :: x
      INTEGER, INTENT(OUT) :: MaxLoc
      INTEGER :: Maximum
      INTEGER :: Location
      INTEGER :: i

      Maximum  = x(xStart)		! assume the first is the max
      Location = xStart			! record its position
      DO i = xStart+1, xEnd		! start with next elements
         IF (x(i) > Maximum) THEN	!   if x(i) higher than the max?
            Maximum  = x(i)		!      Yes, a new maximum found
            Location = i                !      record its position
         END IF
      END DO
      MaxLoc = Location        	! return the position
   END SUBROUTINE FindMaximum
   
   
    ! --------------------------------------------------------------------
    ! SUBROUTINE  Swap():
    !    This subroutine swaps the values of its two formal arguments.
    ! --------------------------------------------------------------------

   SUBROUTINE  Swap(a, b)
      IMPLICIT  NONE
      INTEGER, INTENT(INOUT) :: a, b
      INTEGER                :: Temp

      Temp = a
      a    = b
      b    = Temp
   END SUBROUTINE  Swap

END MODULE DiverseUtil