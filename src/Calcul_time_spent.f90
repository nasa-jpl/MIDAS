SUBROUTINE Calcul_time_spent (instant_init,instant_final, time_s)

USE Initialization

Implicit none

integer, dimension(8), INTENT(IN) :: instant_init
integer, dimension(8), INTENT(IN) :: instant_final
integer, dimension(4), INTENT(OUT) :: time_s

time_s = 0

!! Si on n'est pas au meme mois 
If (instant_init(2)/=instant_final(2)) then 
    !! Mois a 31 jours
    If ((instant_init(2) == 1) .or. (instant_init(2) == 3) .or. (instant_init(2) == 5) .or. (instant_init(2) == 7)&
       .or. (instant_init(2) == 8) .or. (instant_init(2) == 10).or. (instant_init(2) == 12)) Then
       
       time_s(1) = 31 - instant_init(3) + instant_final(3)
    
    !! Mois a 30 jours   
    Elseif ((instant_init(2) == 4) .or. (instant_init(2) == 6) .or. (instant_init(2) == 9) .or. (instant_init(2) == 11)) Then
    
       time_s(1) = 30 - instant_init(3) + instant_final(3)
       
    !! Sinon alors Fevrier   
    Else
        If (((mod(instant_init(1),4) == 0).and. (mod(instant_init(1),100) /= 0)) .or. (mod(instant_init(1),400) == 0) ) Then
            time_s(1) = 29 - instant_init(3) + instant_final(3)
        Else
            time_s(1) = 28 - instant_init(3) + instant_final(3)
        Endif    
    Endif 
Endif

!! Maintenant  on passe a la comparaison des journees
If (instant_init(3)/=instant_final(3)) then 
    time_s(2) = 24 - instant_init(5) + instant_final(5) 
    if (time_s(2) >= 24) Then
        time_s(2) = time_s(2) - 24
    Else
        time_s(1) = time_s(1)-1
    Endif
!! Maintenant on passe a la comparaison des heures
Else
    time_s(2) = instant_final(5) - instant_init(5)     
Endif   

!!Ici on compte la difference en minutes 
If (instant_init(5)/=instant_final(5)) Then
    time_s(3) = 60 - instant_init(6) + instant_final(6)
    if (time_s(3) >= 60) Then
        time_s(3) = time_s(3) - 60         
    Else
        time_s(2) = time_s(2) - 1 
    Endif
Else 
    time_s(3) = instant_final(6) - instant_init(6)
Endif 

If (instant_init(6)/=instant_final(6)) Then
    time_s(4) = 60 - instant_init(7) + instant_final(7)
    if (time_s(4) >= 60) Then
        time_s(4) = time_s(4) - 60         
    Else
        time_s(3) = time_s(3) -1 
    Endif
Else 
    time_s(4) = instant_final(7) - instant_init(7)
Endif  

!! Finalement on compte la difference en secondes 
END SUBROUTINE Calcul_time_spent