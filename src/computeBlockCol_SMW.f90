SUBROUTINE computeBlockCol_SMW(nb1,nb2,CellsB1,CellsB2,icol,Matrix_Green_Col) 
 
    USE Initialization
    USE common_variables
  
    IMPLICIT NONE
    !! IN/OUT ******************************************************************
    INTEGER, INTENT(IN) :: nb1,nb2,icol
    type (Cell), Dimension(nb1), INTENT(IN) :: CellsB1
    type (Cell), Dimension(nb2), INTENT(IN) :: CellsB2
    COMPLEX(KIND=8), DIMENSION(3*nb1,1), INTENT(OUT) :: Matrix_Green_Col
    
    !local 
    INTEGER ii,Is, Isc,Io, Ioc
    INTEGER Matrix_Green_IndexRow_sup, Index_lig,nb_cels_i,Io_init
    COMPLEX :: Term1, Term2, Term3
    COMPLEX :: Fxx, Fyy, Fzz, Fxy, Fyz,Fxz,Fyx,Fzy,Fzx,Gx,Gy,Gz
    Complex :: K11x, K11z, K22z, K12z, K21z
    Real(kind=8) :: RE, RM, Ang_inc, P1, P2, norme, Px, Py
    REAL(kind=8) :: Distx, Disty, Distz, Dist, Distxy
    
    !! commencons par connaitre la cellule concernee par la colonne Matrix_Green_ICol
    Isc = mod(icol,3)
    If (Isc .NE. 0) Then
        Is = icol/3 + 1
    Else 
        Is = icol/3
    Endif  
  
    Index_lig = 0
    
    If (Isc == 1) Then !!.X
        Do Io=1,nb1    
            
            Distx = CellsB1(Io)%Xc - CellsB2(Is)%Xc
            Disty = CellsB1(Io)%Yc - CellsB2(Is)%Yc
            Distz = CellsB1(Io)%Zc - CellsB2(Is)%Zc
            Dist = sqrt(Distx*Distx+Disty*Disty+Distz*Distz)
      
            Term1 = exp(J*K_air*Dist)/Dist**2.
            Term2 = J*K_air-1/Dist
            Term3 = CellsB2(Is)%parameter_Const*CellsB2(Is)%parameter_Ce

            Fxx = (Term2 + K_air**2.*(Dist - Distx**2./Dist)-3.*Term2*Distx**2./Dist**2.)*Term1
            Fxy = (Distx*Disty*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
            Fxz = (Distx*Distz*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
      
            Matrix_Green_Col(Index_lig+1,1)= - Fxx*Term3 !!XX
            Matrix_Green_Col(Index_lig+2,1)= - Fxy*Term3 !!YX
            Matrix_Green_Col(Index_lig+3,1)= - Fxz*Term3 !!ZX

            Index_lig = Index_lig + 3
        Enddo
    ElseIf (Isc==2) Then !! .Y
        Do Io=1,nb1       
            
            Distx = CellsB1(Io)%Xc - CellsB2(Is)%Xc
            Disty = CellsB1(Io)%Yc - CellsB2(Is)%Yc
            Distz = CellsB1(Io)%Zc - CellsB2(Is)%Zc

            Dist = sqrt(Distx*Distx+Disty*Disty+Distz*Distz)

            Term1 = exp(J*K_air*Dist)/Dist**2.
            Term2 = J*K_air-1/Dist
            Term3 = CellsB2(Is)%parameter_Const*CellsB2(Is)%parameter_Ce

            Fxy = (Distx*Disty*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
            Fyy = (Term2 + K_air**2.*(Dist - Disty**2./Dist)-3.*Term2*Disty**2./Dist**2.)*Term1
            Fyz = (Disty*Distz*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1            
          
            Matrix_Green_Col(Index_lig+1,1)= - Fxy*Term3 !!XY
            Matrix_Green_Col(Index_lig+2,1)= - Fyy*Term3 !!YY
            Matrix_Green_Col(Index_lig+3,1)= - Fyz*Term3 !!ZY           
          
            Index_lig = Index_lig + 3 
        Enddo
    Else  !! .Z
        Do Io=1,nb1      
            
            Distx = CellsB1(Io)%Xc - CellsB2(Is)%Xc
            Disty = CellsB1(Io)%Yc - CellsB2(Is)%Yc
            Distz = CellsB1(Io)%Zc - CellsB2(Is)%Zc

            Dist = sqrt(Distx*Distx+Disty*Disty+Distz*Distz)

            Term1 = exp(J*K_air*Dist)/Dist**2.
            Term2 = J*K_air-1/Dist
            Term3 = CellsB2(Is)%parameter_Const*CellsB2(Is)%parameter_Ce
            
            Fxz = (Distx*Distz*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
            Fyz = (Disty*Distz*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
            Fzz = (Term2 + K_air**2.*(Dist - Distz**2./Dist)-3.*Term2*Distz**2./Dist**2.)*Term1
            
            Matrix_Green_Col(Index_lig+1,1)= -Fxz *Term3 !!XZ
            Matrix_Green_Col(Index_lig+2,1)= -Fyz*Term3  !!YZ
            Matrix_Green_Col(Index_lig+3,1)= -Fzz*Term3  !!ZZ           
          
            Index_lig = Index_lig + 3 
        Enddo
    EndIf

 END SUBROUTINE computeBlockCol_SMW