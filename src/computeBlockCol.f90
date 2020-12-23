SUBROUTINE computeBlockCol(Cells,Matrix_Green_IndexRow_inf,nbRows,Matrix_Green_ICol,Matrix_Green_Col) 
 
    USE Initialization
    USE common_variables
  
    IMPLICIT NONE
    !! IN/OUT ******************************************************************
    type (Cell), Dimension(Nbc), INTENT(IN) :: Cells
    INTEGER, INTENT(IN) :: Matrix_Green_IndexRow_inf,nbRows,Matrix_Green_ICol
    COMPLEX(KIND=8), DIMENSION(nbRows,1),INTENT(OUT) :: Matrix_Green_Col
    
    ! local
    INTEGER ii,Is, Isc,Io, Ioc
    INTEGER Matrix_Green_IndexRow_sup, curs_lig,nb_cels_i,Io_init
    COMPLEX :: Term1, Term2, Term3
    COMPLEX :: Fxx, Fyy, Fzz, Fxy, Fyz,Fxz,Fyx,Fzy,Fzx,Gx,Gy,Gz
    Complex :: K11x, K11z, K22z, K12z, K21z
    Real(kind=8) :: RE, RM, Ang_inc, P1, P2, norme, Px, Py
    REAL(kind=8) :: Distx, Disty, Distz, Dist, Distxy
    
    !! commencons par connaitre la cellule concernee par la colonne Matrix_Green_ICol
    Isc = mod(Matrix_Green_ICol,3)
    If (Isc .NE. 0) Then
        Is = Matrix_Green_ICol/3 + 1
    Else 
        Is = Matrix_Green_ICol/3
    Endif  
  
    Curs_lig = 0
    Io_init = Matrix_Green_IndexRow_inf/3 + 1
    nb_cels_i = nbRows/3
  
    Io = Io_init  
    If (Isc == 1) Then !!.X
        Do ii=1,nb_cels_i    
            
            Distx = Cells(Io)%Xc - Cells(Is)%Xc
            Disty = Cells(Io)%Yc - Cells(Is)%Yc
            Distz = Cells(Io)%Zc - Cells(Is)%Zc
            Dist = sqrt(Distx*Distx+Disty*Disty+Distz*Distz)
      
            Term1 = exp(J*K_air*Dist)/Dist**2.
            Term2 = J*K_air-1/Dist
            Term3 = Cells(Is)%parameter_Const*Cells(Is)%parameter_Ce

            Fxx = (Term2 + K_air**2.*(Dist - Distx**2./Dist)-3.*Term2*Distx**2./Dist**2.)*Term1
            Fxy = (Distx*Disty*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
            Fxz = (Distx*Distz*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
      
            Matrix_Green_Col(curs_lig+1,1)= - Fxx*Term3 !!XX
            Matrix_Green_Col(curs_lig+2,1)= - Fxy*Term3 !!YX
            Matrix_Green_Col(curs_lig+3,1)= - Fxz*Term3 !!ZX

            Curs_lig = Curs_lig + 3 
            Io = Io + 1
        Enddo
    ElseIf (Isc==2) Then !! .Y
        Do ii=1,nb_cels_i       
            
            Distx = Cells(Io)%Xc - Cells(Is)%Xc
            Disty = Cells(Io)%Yc - Cells(Is)%Yc
            Distz = Cells(Io)%Zc - Cells(Is)%Zc

            Dist = sqrt(Distx*Distx+Disty*Disty+Distz*Distz)

            Term1 = exp(J*K_air*Dist)/Dist**2.
            Term2 = J*K_air-1/Dist
            Term3 = Cells(Is)%parameter_Const*Cells(Is)%parameter_Ce

            Fxy = (Distx*Disty*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
            Fyy = (Term2 + K_air**2.*(Dist - Disty**2./Dist)-3.*Term2*Disty**2./Dist**2.)*Term1
            Fyz = (Disty*Distz*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
            
          
            Matrix_Green_Col(curs_lig+1,1)= - Fxy*Term3 !!XY
            Matrix_Green_Col(curs_lig+2,1)= - Fyy*Term3 !!YY
            Matrix_Green_Col(curs_lig+3,1)= - Fyz*Term3 !!ZY           
          
            Curs_lig = Curs_lig + 3 
            Io = Io + 1
        Enddo
    Else  !! .Z
        Do ii=1,nb_cels_i      
            
            Distx = Cells(Io)%Xc - Cells(Is)%Xc
            Disty = Cells(Io)%Yc - Cells(Is)%Yc
            Distz = Cells(Io)%Zc - Cells(Is)%Zc

            Dist = sqrt(Distx*Distx+Disty*Disty+Distz*Distz)

            Term1 = exp(J*K_air*Dist)/Dist**2.
            Term2 = J*K_air-1/Dist
            Term3 = Cells(Is)%parameter_Const*Cells(Is)%parameter_Ce
            
            Fxz = (Distx*Distz*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
            Fyz = (Disty*Distz*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
            Fzz = (Term2 + K_air**2.*(Dist - Distz**2./Dist)-3.*Term2*Distz**2./Dist**2.)*Term1
            
            Matrix_Green_Col(curs_lig+1,1)= -Fxz *Term3 !!XZ
            Matrix_Green_Col(curs_lig+2,1)= -Fyz*Term3  !!YZ
            Matrix_Green_Col(curs_lig+3,1)= -Fzz*Term3  !!ZZ           
          
            Curs_lig = Curs_lig + 3 
            Io = Io + 1
        Enddo
    EndIf

 END SUBROUTINE computeBlockCol