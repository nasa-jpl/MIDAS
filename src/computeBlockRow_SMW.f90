SUBROUTINE computeBlockRow_SMW(nb1,nb2,CellsB1,CellsB2,irow,Matrix_Green_row) 
  
    USE Initialization
    USE common_variables
    
    IMPLICIT NONE

    !! IN/OUT ******************************************************************
    !type (Cell), Dimension(Nbc), INTENT(IN) :: Cells
    !INTEGER, INTENT(IN) :: Matrix_Green_irow, Matrix_Green_IndexCol_inf,nbCols
    !COMPLEX(KIND=8), DIMENSION(1,nbCols), INTENT(OUT) :: Matrix_Green_row
    INTEGER, INTENT(IN) :: nb1,nb2,irow
    type (Cell), Dimension(nb1), INTENT(IN) :: CellsB1
    type (Cell), Dimension(nb2), INTENT(IN) :: CellsB2
    
    COMPLEX(KIND=8), DIMENSION(1,3*nb2), INTENT(OUT) :: Matrix_Green_row
    
    ! local 
    INTEGER jj,Is,Isc,Io,Ioc
    INTEGER Index_Col, Is_init, Nb_cels_j
    COMPLEX :: Term1, Term2, Term3
    COMPLEX :: Fxx, Fyy, Fzz, Fxy, Fyz,Fxz,Fyx,Fzy,Fzx,Gx,Gy,Gz
    Complex :: K11x, K11z, K22z, K12z, K21z
    Real(kind=8) :: RE, RM, Ang_inc, P1, P2, norme, Px, Py
    Real(kind=8) :: Distx, Disty, Distz, Dist, Distxy
  
    
    Ioc= mod(irow,3)
    If (Ioc .NE. 0) Then
        Io = irow/3 + 1
    Else
        Io = irow/3
    Endif
    Index_Col = 0  
    
    If (Ioc ==1) Then !! X.
        Do Is=1,nb2      
            
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

            Matrix_Green_row(1,Index_Col+1)= - Fxx*Term3 !XX
            Matrix_Green_row(1,Index_Col+2)= - Fxy*Term3 !XY
            Matrix_Green_row(1,Index_Col+3)= - Fxz*Term3 !XZ             
                    
            Index_Col = Index_Col +3
        Enddo
    ElseIf (Ioc == 2) Then  !!Y.
        Do Is=1,nb2      
            
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
            
            Matrix_Green_row(1,Index_Col+1)= - Fxy*Term3   !YX
            Matrix_Green_row(1,Index_Col+2)= - Fyy*Term3   !YY
            Matrix_Green_row(1,Index_Col+3)= - Fyz*Term3   !YZ
          
            Index_Col = Index_Col +3
        Enddo
    Else !!Z.
        Do Is=1,nb2      
            
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

            Matrix_Green_row(1,Index_Col+1)= - Fxz*Term3 !ZX
            Matrix_Green_row(1,Index_Col+2)= - Fyz*Term3 !ZY
            Matrix_Green_row(1,Index_Col+3)= - Fzz*Term3 !ZZ
          
            Index_Col = Index_Col +3
        Enddo
    EndIf  
 END SUBROUTINE computeBlockRow_SMW