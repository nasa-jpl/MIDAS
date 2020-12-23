SUBROUTINE Green_s_tr_partial(sizeBlock1,CellsBlock1,sizeBlock2,CellsBlock2,Green_s_tr)

    USE Initialization
    USE common_variables
    Implicit none
    
    !IN/OUT
    Integer, INTENT(IN) :: sizeBlock1,sizeBlock2
    type (Cell), Dimension(sizeBlock1), INTENT(IN) :: CellsBlock1
    type (Cell), Dimension(sizeBlock2), INTENT(IN) :: CellsBlock2
    COMPLEX(real64), Dimension(3*sizeBlock1,3*sizeBlock2),INTENT(OUT)::Green_s_tr

    ! Local
    Complex	:: Term1, Term2, Term3, Fxx, Fyy, Fzz, Fxy, Fyz, Fxz
    Real (kind=8) :: Distx, Disty, Distz, Dist
    Integer :: Is, Io, Isg, Iog, Iox, Ioy, Ioz, Isx, Isy, Isz
        
    

    Do Is=1, sizeBlock2
      Do Io=1, sizeBlock1  
      
        Iox = 3*(Io-1)+1
        Ioy = 3*(Io-1)+2
        Ioz = 3*(Io-1)+3

        Isx = 3*(Is-1)+1
        Isy = 3*(Is-1)+2
        Isz = 3*(Is-1)+3
    
        Iog = CellsBlock1(Io)%num_cell
        Isg = CellsBlock2(Is)%num_cell
        if (Iog == Isg) then
          Green_s_tr(Iox,Isx) = 1 - CellsBlock2(Is)%parameter_Sing
          Green_s_tr(Ioy,Isy) = 1 - CellsBlock2(Is)%parameter_Sing
          Green_s_tr(Ioz,Isz) = 1 - CellsBlock2(Is)%parameter_Sing
      
          Green_s_tr(Iox,Isy) = 0
          Green_s_tr(Iox,Isz) = 0
          Green_s_tr(Ioy,Isx) = 0
          Green_s_tr(Ioz,Isx) = 0
          Green_s_tr(Ioy,Isz) = 0
          Green_s_tr(Ioz,Isy) = 0    
	    else
          Distx = CellsBlock1(Io)%Xc - CellsBlock2(Is)%Xc
          Disty = CellsBlock1(Io)%Yc - CellsBlock2(Is)%Yc
          Distz = CellsBlock1(Io)%Zc - CellsBlock2(Is)%Zc

          Dist = sqrt(Distx*Distx+Disty*Disty+Distz*Distz)

          Term1 = exp(J*K_air*Dist)/Dist**2.
          Term2 = J*K_air-1/Dist
          Term3 = CellsBlock2(Is)%parameter_Const*CellsBlock2(Is)%parameter_Ce

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


    End Subroutine Green_s_tr_partial
    
    
    SUBROUTINE Green_s_tr_partial_FN(sizeBlock1,CellsBlock1,sizeBlock2,CellsBlock2,FN_Green_s_tr)

    ! this subroutine will simply enable us to determine the Frobenius norm without storing Z patch (memory efficient)
    USE Initialization
    USE common_variables
    USE iso_fortran_env
    
    Implicit none
    
    !IN/OUT
    Integer, INTENT(IN) :: sizeBlock1,sizeBlock2
    type (Cell), Dimension(sizeBlock1), INTENT(IN) :: CellsBlock1
    type (Cell), Dimension(sizeBlock2), INTENT(IN) :: CellsBlock2
    Real(kind=8), INTENT(OUT) :: FN_Green_s_tr

    ! Local
    Complex	:: Term1, Term2, Term3, Fxx, Fyy, Fzz, Fxy, Fyz, Fxz
    Real (kind=8) :: Distx, Disty, Distz, Dist,val_fn
    Integer :: Is, Io, Isg, Iog, Iox, Ioy, Ioz, Isx, Isy, Isz,Io_beg
        
    
    val_fn = 0D0;
  
    if (homogs .eq. 0) then 
        Do Is=1, sizeBlock2
          Do Io=1, sizeBlock1  
      
            Iox = 3*(Io-1)+1
            Ioy = 3*(Io-1)+2
            Ioz = 3*(Io-1)+3

            Isx = 3*(Is-1)+1
            Isy = 3*(Is-1)+2
            Isz = 3*(Is-1)+3
    
            Iog = CellsBlock1(Io)%num_cell
            Isg = CellsBlock2(Is)%num_cell
            if (Iog == Isg) then
              val_fn = val_fn + 3.*(abs(1 - CellsBlock2(Is)%parameter_Sing))**2          
	        else
              Distx = CellsBlock1(Io)%Xc - CellsBlock2(Is)%Xc
              Disty = CellsBlock1(Io)%Yc - CellsBlock2(Is)%Yc
              Distz = CellsBlock1(Io)%Zc - CellsBlock2(Is)%Zc

              Dist = sqrt(Distx*Distx+Disty*Disty+Distz*Distz)

              Term1 = exp(J*K_air*Dist)/Dist**2.
              Term2 = J*K_air-1/Dist
              Term3 = CellsBlock2(Is)%parameter_Const*CellsBlock2(Is)%parameter_Ce

              Fxx = (Term2 + K_air**2.*(Dist - Distx**2./Dist)-3.*Term2*Distx**2./Dist**2.)*Term1
              Fyy = (Term2 + K_air**2.*(Dist - Disty**2./Dist)-3.*Term2*Disty**2./Dist**2.)*Term1
              Fzz = (Term2 + K_air**2.*(Dist - Distz**2./Dist)-3.*Term2*Distz**2./Dist**2.)*Term1

              Fxy = (Distx*Disty*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
              Fyz = (Disty*Distz*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
              Fxz = (Distx*Distz*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
          
              val_fn = val_fn + ((abs(Fxx))**2+2.*(abs(Fxy))**2+2.*(abs(Fxz))**2 &
                                +(abs(Fyy))**2+2.*(abs(Fyz))**2+(abs(Fzz))**2)*(abs(Term3))**2. 
            endIf    
          endDo
        endDo 
    else
        Do Is=1, sizeBlock2
          Do Io=Is, sizeBlock1  
      
            Iox = 3*(Io-1)+1
            Ioy = 3*(Io-1)+2
            Ioz = 3*(Io-1)+3

            Isx = 3*(Is-1)+1
            Isy = 3*(Is-1)+2
            Isz = 3*(Is-1)+3
    
            Iog = CellsBlock1(Io)%num_cell
            Isg = CellsBlock2(Is)%num_cell
            if (Iog == Isg) then
              val_fn = val_fn + 3.*(abs(1 - CellsBlock2(Is)%parameter_Sing))**2          
	        else
              Distx = CellsBlock1(Io)%Xc - CellsBlock2(Is)%Xc
              Disty = CellsBlock1(Io)%Yc - CellsBlock2(Is)%Yc
              Distz = CellsBlock1(Io)%Zc - CellsBlock2(Is)%Zc

              Dist = sqrt(Distx*Distx+Disty*Disty+Distz*Distz)

              Term1 = exp(J*K_air*Dist)/Dist**2.
              Term2 = J*K_air-1/Dist
              Term3 = CellsBlock2(Is)%parameter_Const*CellsBlock2(Is)%parameter_Ce

              Fxx = (Term2 + K_air**2.*(Dist - Distx**2./Dist)-3.*Term2*Distx**2./Dist**2.)*Term1
              Fyy = (Term2 + K_air**2.*(Dist - Disty**2./Dist)-3.*Term2*Disty**2./Dist**2.)*Term1
              Fzz = (Term2 + K_air**2.*(Dist - Distz**2./Dist)-3.*Term2*Distz**2./Dist**2.)*Term1

              Fxy = (Distx*Disty*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
              Fyz = (Disty*Distz*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
              Fxz = (Distx*Distz*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
          
              val_fn = val_fn + ((abs(Fxx))**2+2.*(abs(Fxy))**2+2.*(abs(Fxz))**2 &
                                +(abs(Fyy))**2+2.*(abs(Fyz))**2+(abs(Fzz))**2)*(abs(Term3))**2. 
            endIf    
          endDo
        endDo 
    endif
    

    FN_Green_s_tr = sqrt(val_fn);

End Subroutine Green_s_tr_partial_FN