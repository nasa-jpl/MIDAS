SUBROUTINE TrouverDiviseurCommun(int_a,int_div,intDiv)
 
!! Trouve un diviseur superieur a 12 
USE Initialization

Implicit none

integer, INTENT(IN) :: int_a, int_div
integer :: bon, diviseur, q
integer, INTENT(OUT) :: intDiv

diviseur = int_div
bon = 0
Do while (bon==0)
    q = int_a/diviseur
    If (q*diviseur == int_a) then
        intDiv = diviseur 
        bon = 1
    else 
        diviseur = diviseur + 1
    Endif
Enddo

End Subroutine TrouverDiviseurCommun