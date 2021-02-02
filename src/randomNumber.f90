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