SUBROUTINE Green_s_tr_total(Cells,Green_s_tr)

    USE Initialization
    USE common_variables
    Implicit none

    ! IN/OUT
    type (Cell), Dimension(Nbc), INTENT(IN) :: Cells
    COMPLEX(real64), Dimension(3*Nbc,3*Nbc), INTENT(OUT):: Green_s_tr
    
    ! Local 
    Complex :: Term1, Term2, Term3, Fxx, Fyy, Fzz, Fxy, Fyz, Fxz
    Real(kind=8) :: Distx, Disty, Distz, Dist
    Integer :: Is, Io, Iox, Ioy, Ioz, Isx, Isy, Isz

    Green_s_tr = 0

    Do Is=1, Nbc
      Do Io=1, Nbc  
      
        Iox = 3*(Io-1)+1
        Ioy = 3*(Io-1)+2
        Ioz = 3*(Io-1)+3

        Isx = 3*(Is-1)+1
        Isy = 3*(Is-1)+2
        Isz = 3*(Is-1)+3
    
        if (Io==Is)then
          Green_s_tr(Iox,Isx) = 1 - Cells(Is)%parameter_Sing
          Green_s_tr(Ioy,Isy) = 1 - Cells(Is)%parameter_Sing
          Green_s_tr(Ioz,Isz) = 1 - Cells(Is)%parameter_Sing
      
          Green_s_tr(Iox,Isy) = 0
          Green_s_tr(Iox,Isz) = 0
          Green_s_tr(Ioy,Isx) = 0
          Green_s_tr(Ioz,Isx) = 0
          Green_s_tr(Ioy,Isz) = 0
          Green_s_tr(Ioz,Isy) = 0    
	    else
          Distx = Cells(Io)%Xc - Cells(Is)%Xc
          Disty = Cells(Io)%Yc - Cells(Is)%Yc
          Distz = Cells(Io)%Zc - Cells(Is)%Zc

          Dist = sqrt(Distx*Distx+Disty*Disty+Distz*Distz)

          Term1 = exp(J*K_air*Dist)/Dist**2.
          Term2 = J*K_air-1/Dist
          Term3 = Cells(Is)%parameter_Const*Cells(Is)%parameter_Ce

          Fxx = (Term2 + K_air**2.*(Dist - Distx**2./Dist)-3.*Term2*Distx**2./Dist**2.)*Term1
          Fyy = (Term2 + K_air**2.*(Dist - Disty**2./Dist)-3.*Term2*Disty**2./Dist**2.)*Term1
          Fzz = (Term2 + K_air**2.*(Dist - Distz**2./Dist)-3.*Term2*Distz**2./Dist**2.)*Term1

          Fxy = (Distx*Disty*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
          Fyz = (Disty*Distz*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
          Fxz = (Distx*Distz*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
      
          Green_s_tr(Iox,Isx)= -Fxx*Term3 !XX
          Green_s_tr(Iox,Isy)= -Fxy*Term3 !XY
          Green_s_tr(Iox,Isz)= -Fxz*Term3 !XZ
          Green_s_tr(Ioy,Isx)= -Fxy*Term3 !YX
          Green_s_tr(Ioy,Isy)= -Fyy*Term3 !YY
          Green_s_tr(Ioy,Isz)= -Fyz*Term3 !YZ
          Green_s_tr(Ioz,Isx)= -Fxz*Term3 !ZX
          Green_s_tr(Ioz,Isy)= -Fyz*Term3 !ZY
          Green_s_tr(Ioz,Isz)= -Fzz*Term3 !ZZ
                 
        endIf    
      endDo
    endDo 
    
End Subroutine Green_s_tr_total