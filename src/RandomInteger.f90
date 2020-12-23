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