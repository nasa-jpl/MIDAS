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

    real_al = aint(10**(Round_Dp+1)*real_al)/10**(Round_Dp+1)

End Subroutine RandomReal