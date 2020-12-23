SUBROUTINE DR_Green_s_tr_partial(sizeBlock,CellsBlock,klu_cel,klu,Green_s_tr)

    USE Initialization
    USE common_variables
    Implicit none
    
    !IN/OUT
    Integer, INTENT(IN) :: sizeBlock
    type (Cell), Dimension(sizeBlock), INTENT(IN) :: CellsBlock
    Integer, INTENT(IN):: klu_cel,klu
    COMPLEX(real64),Dimension(2*klu+1,3*sizeBlock),INTENT(OUT)::Green_s_tr

    ! Local
    Complex	:: Term1, Term2, Term3, Fxx, Fyy, Fzz, Fxy, Fyz, Fxz
    Real (kind=8) :: Distx, Disty, Distz, Dist
    Integer :: Is, Io, Isg, Iog, Iox, Ioy, Ioz, Isx, Isy, Isz  
    Integer :: Index_col_Inf,Index_col_Sup,Iox_band,Ioy_band,Ioz_band
    
    Green_s_tr(:,:) = 0.D0;

    Do Io = 1, sizeBlock
    
        Index_col_Inf = max(1,Io-klu_cel)
        Index_col_Sup = min(Io+klu_cel,sizeBlock);
    
        Do Is = Index_col_Inf, Index_col_Sup  
                
            !! indices ligne 
            Iox = 3*(Io-1)+1
            Ioy = 3*(Io-1)+2
            Ioz = 3*(Io-1)+3
        
            !! indices colonne 
            Isx = 3*(Is-1)+1
            Isy = 3*(Is-1)+2
            Isz = 3*(Is-1)+3

            if (Io==Is) then           
                !! BAND STORAGE : a(i,j) is stored in ab(ku+1+i-j,j)
                Iox_band = klu+1+Iox-Isx
                Green_s_tr(Iox_band,Isx) = 1 - CellsBlock(Is)%parameter_Sing
                Ioy_band = klu+1+Ioy-Isy
                Green_s_tr(Ioy_band,Isy) = 1 - CellsBlock(Is)%parameter_Sing
                Ioz_band = klu+1+Ioz-Isz
                Green_s_tr(Ioz_band,Isz) = 1 - CellsBlock(Is)%parameter_Sing
            
                !! les autres (xy, yx, xz ...) restent a 0 pour ce cas ()
	        else
              Distx = CellsBlock(Io)%Xc - CellsBlock(Is)%Xc
              Disty = CellsBlock(Io)%Yc - CellsBlock(Is)%Yc
              Distz = CellsBlock(Io)%Zc - CellsBlock(Is)%Zc

              Dist = sqrt(Distx*Distx+Disty*Disty+Distz*Distz)
      
              Term1 = exp(J*K_air*Dist)/Dist**2.
              Term2 = J*K_air-1/Dist
              Term3 = CellsBlock(Is)%parameter_Const*CellsBlock(Is)%parameter_Ce

              Fxx = (Term2 + K_air**2.*(Dist - Distx**2./Dist)-3.*Term2*Distx**2./Dist**2.)*Term1
              Fyy = (Term2 + K_air**2.*(Dist - Disty**2./Dist)-3.*Term2*Disty**2./Dist**2.)*Term1
              Fzz = (Term2 + K_air**2.*(Dist - Distz**2./Dist)-3.*Term2*Distz**2./Dist**2.)*Term1

              Fxy = (Distx*Disty*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
              Fyz = (Disty*Distz*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
              Fxz = (Distx*Distz*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1      
      
              Iox_band = klu+1+Iox-Isx
              Green_s_tr(Iox_band,Isx)= -Fxx*Term3 !XX
              Iox_band = klu+1+Iox-Isy
              Green_s_tr(Iox_band,Isy)= -Fxy*Term3 !XY
              Iox_band = klu+1+Iox-Isz
              Green_s_tr(Iox_band,Isz)= -Fxz*Term3 !XZ
          
              Ioy_band = klu+1+Ioy-Isx
              Green_s_tr(Ioy_band,Isx)= -Fxy*Term3 !YX
              Ioy_band = klu+1+Ioy-Isy
              Green_s_tr(Ioy_band,Isy)= -Fyy*Term3 !YY
              Ioy_band = klu+1+Ioy-Isz
              Green_s_tr(Ioy_band,Isz)= -Fyz*Term3 !YZ
          
              Ioz_band = klu+1+Ioz-Isx
              Green_s_tr(Ioz_band,Isx)= -Fxz*Term3 !ZX
              Ioz_band = klu+1+Ioz-Isy
              Green_s_tr(Ioz_band,Isy)= -Fyz*Term3 !ZY
              Ioz_band = klu+1+Ioz-Isz
              Green_s_tr(Ioz_band,Isz)= -Fzz*Term3 !ZZ  
            endIf    
          endDo
    endDo 


End Subroutine DR_Green_s_tr_partial