SUBROUTINE SR_Green_s_tr_partial(sizeBlock,CellsBlock,localfSR,nnz,Green_s_tr,irow,icol)           

    USE Initialization
    USE common_variables
    Implicit none
    
    !IN/OUT
    Integer, INTENT(IN) :: sizeBlock,nnz
    real(kind=8), INTENT(IN) :: localfSR
    Integer, Dimension(3*sizeBlock+1) :: irow
    Integer, Dimension(nnz) :: icol
    type (Cell), Dimension(sizeBlock),INTENT(IN) :: CellsBlock
    COMPLEX(real64),Dimension(nnz) ::Green_s_tr

    ! Local    
    Integer :: curs_nnz,ii,jj,Is, Io, Isg, Iog, Iox, Ioy, Ioz, Isx, Isy, Isz
    Integer :: Ioo,num_dir
    Real (kind=8) :: Distx, Disty, Distz, Dist
    Real (kind=8) :: v_max,Green_s_tr_max,v_irow,threshold
    Real (kind=8),dimension(3) :: abs_g
    Complex	:: Term1, Term2, Term3, Fxx, Fyy, Fzz, Fxy, Fyz, Fxz
    complex,dimension(3) :: g        
    
    Green_s_tr_max =0D0
    curs_nnz = 0;
    v_irow = 1;
    
    ! start by computing Green_s_tr(1,1), it will be used to elliminate 
    ! what we will consider as weak/non-significant interactions
    v_max = abs(1 - CellsBlock(1)%parameter_Sing)
    threshold = v_max/localfSR; !fct_SR;  ! the global fct_SR can be used if not all the blocks tested   
    
    If (homogs == 1) Then     
      Do Io = 1, sizeBlock   ! loop on row
        
        Do num_dir=1,3            
          ! Is = Io 
          Ioo = 3*(Io-1) + num_dir
          irow(Ioo) = v_irow ! start by storing the index (%nnz elts) 
                             ! of the first nnz elt in this row  
          
          !IoIo dir-dir (Io == Is)
          curs_nnz = curs_nnz + 1; icol(curs_nnz) = Ioo; 
          Green_s_tr(curs_nnz) = 1 - CellsBlock(Io)%parameter_Sing ;                
          v_irow = v_irow + 1 
          
          ! Now Is > Io
          Do Is=Io+1, sizeBlock     !loop on col  ! if homogs=0 we will scan the entire matrix    
            Isx = 3*(Is-1)+1
            Isy = 3*(Is-1)+2
            Isz = 3*(Is-1)+3
            
            Distx = CellsBlock(Io)%Xc - CellsBlock(Is)%Xc
            Disty = CellsBlock(Io)%Yc - CellsBlock(Is)%Yc
            Distz = CellsBlock(Io)%Zc - CellsBlock(Is)%Zc
    
            Dist = sqrt(Distx*Distx+Disty*Disty+Distz*Distz)
    
            Term1 = exp(J*K_air*Dist)/Dist**2.
            Term2 = J*K_air-1/Dist
            Term3 = CellsBlock(Is)%parameter_Const*CellsBlock(Is)%parameter_Ce
            
            if (num_dir == 1) then  !Iox
              !f
              Fxx = (Term2 + K_air**2.*(Dist - Distx**2./Dist)-3.*Term2*Distx**2./Dist**2.)*Term1
              Fxy = (Distx*Disty*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
              Fxz = (Distx*Distz*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
              !g
              g(1) = -Fxx*Term3 !XX :Green_s_tr(Iox,Isx)
              g(2) = -Fxy*Term3 !XY : Green_s_tr(Iox,Isy)
              g(3) = -Fxz*Term3 !XZ : Green_s_tr(Iox,Isz)
            elseif (num_dir == 2) then !Ioy
              !f
              Fxy = (Distx*Disty*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
              Fyy = (Term2 + K_air**2.*(Dist - Disty**2./Dist)-3.*Term2*Disty**2./Dist**2.)*Term1
              Fyz = (Disty*Distz*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
              !g
              g(1) = -Fxy*Term3 !YX :  Green_s_tr(Ioy,Isx)
              g(2) = -Fyy*Term3 !YY : Green_s_tr(Ioy,Isy)
              g(3) = -Fyz*Term3 !YZ : Green_s_tr(Ioy,Isz)
            else !Ioz
              !f
              Fxz = (Distx*Distz*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
              Fyz = (Disty*Distz*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
              Fzz = (Term2 + K_air**2.*(Dist - Distz**2./Dist)-3.*Term2*Distz**2./Dist**2.)*Term1
              !g
              g(1) = -Fxz*Term3 !ZX : Green_s_tr(Ioz,Isx)
              g(2) = -Fyz*Term3 !ZY : Green_s_tr(Ioz,Isy)
              g(3) = -Fzz*Term3 !ZZ : Green_s_tr(Ioz,Isz)
            endif
    
            abs_g = abs(g); 
            
            If (abs_g(1) .gt. threshold)  then
              curs_nnz = curs_nnz + 1; icol(curs_nnz) = Isx; 
              v_irow = v_irow + 1;
              Green_s_tr(curs_nnz)= g(1);
            endif 
            If (abs_g(2) .gt. threshold)  then
              curs_nnz = curs_nnz + 1; icol(curs_nnz) = Isy; 
              v_irow = v_irow + 1;
              Green_s_tr(curs_nnz)= g(2);  
            endif          
            If (abs_g(3) .gt. threshold)  then 
              curs_nnz = curs_nnz + 1; icol(curs_nnz) = Isz; 
              v_irow = v_irow + 1;
              Green_s_tr(curs_nnz)= g(3); 
            endif 
            !Write(*,*) 'nnz = ',nnz
            !Write(*,*) 'Green_s_tr(nnz) =', Green_s_tr(nnz)    
          EndDo
        Enddo
      EndDo 
    Else    
        Do Io = 1, sizeBlock   ! loop on row        
            Do num_dir=1,3  
          
                Ioo = 3*(Io-1) + num_dir
                irow(Ioo) = v_irow ! start by storing the index (%nnz elts) 
                                    ! of the first nnz elt in this row
                               
                ! here Is < Io
                Do Is=1, Io-1       
                Isx = 3*(Is-1)+1
                Isy = 3*(Is-1)+2
                Isz = 3*(Is-1)+3
            
                Distx = CellsBlock(Io)%Xc - CellsBlock(Is)%Xc
                Disty = CellsBlock(Io)%Yc - CellsBlock(Is)%Yc
                Distz = CellsBlock(Io)%Zc - CellsBlock(Is)%Zc
    
                Dist = sqrt(Distx*Distx+Disty*Disty+Distz*Distz)
    
                Term1 = exp(J*K_air*Dist)/Dist**2.
                Term2 = J*K_air-1/Dist
                Term3 = CellsBlock(Is)%parameter_Const*CellsBlock(Is)%parameter_Ce
            
                if (num_dir == 1) then  !Iox
                    !f
                    Fxx = (Term2 + K_air**2.*(Dist - Distx**2./Dist)-3.*Term2*Distx**2./Dist**2.)*Term1
                    Fxy = (Distx*Disty*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
                    Fxz = (Distx*Distz*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
                    !g
                    g(1) = -Fxx*Term3 !XX :Green_s_tr(Iox,Isx)
                    g(2) = -Fxy*Term3 !XY : Green_s_tr(Iox,Isy)
                    g(3) = -Fxz*Term3 !XZ : Green_s_tr(Iox,Isz)
                elseif (num_dir == 2) then !Ioy
                    !f
                    Fxy = (Distx*Disty*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
                    Fyy = (Term2 + K_air**2.*(Dist - Disty**2./Dist)-3.*Term2*Disty**2./Dist**2.)*Term1
                    Fyz = (Disty*Distz*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
                    !g
                    g(1) = -Fxy*Term3 !YX :  Green_s_tr(Ioy,Isx)
                    g(2) = -Fyy*Term3 !YY : Green_s_tr(Ioy,Isy)
                    g(3) = -Fyz*Term3 !YZ : Green_s_tr(Ioy,Isz)
                else !Ioz
                    !f
                    Fxz = (Distx*Distz*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
                    Fyz = (Disty*Distz*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
                    Fzz = (Term2 + K_air**2.*(Dist - Distz**2./Dist)-3.*Term2*Distz**2./Dist**2.)*Term1
                    !g
                    g(1) = -Fxz*Term3 !ZX : Green_s_tr(Ioz,Isx)
                    g(2) = -Fyz*Term3 !ZY : Green_s_tr(Ioz,Isy)
                    g(3) = -Fzz*Term3 !ZZ : Green_s_tr(Ioz,Isz)
                endif
    
                abs_g = abs(g); 
            
                If (abs_g(1) .gt. threshold)  then
                    curs_nnz = curs_nnz + 1; icol(curs_nnz) = Isx; 
                    v_irow = v_irow + 1;
                    Green_s_tr(curs_nnz)= g(1);
                endif 
                If (abs_g(2) .gt. threshold)  then
                    curs_nnz = curs_nnz + 1; icol(curs_nnz) = Isy; 
                    v_irow = v_irow + 1;
                    Green_s_tr(curs_nnz)= g(2);  
                endif          
                If (abs_g(3) .gt. threshold)  then 
                    curs_nnz = curs_nnz + 1; icol(curs_nnz) = Isz; 
                    v_irow = v_irow + 1;
                    Green_s_tr(curs_nnz)= g(3); 
                endif      
                EndDo
        
                ! Is = Io
                !IoIo dir-dir (Io == Is)
                curs_nnz = curs_nnz + 1; icol(curs_nnz) = Ioo; 
                Green_s_tr(curs_nnz) = 1 - CellsBlock(Io)%parameter_Sing ;                
                v_irow = v_irow + 1 
          
                ! Now Is > Io
                Do Is=Io+1, sizeBlock     !loop on col  ! if homogs=0 we will scan the entire matrix    
                Isx = 3*(Is-1)+1
                Isy = 3*(Is-1)+2
                Isz = 3*(Is-1)+3
            
                Distx = CellsBlock(Io)%Xc - CellsBlock(Is)%Xc
                Disty = CellsBlock(Io)%Yc - CellsBlock(Is)%Yc
                Distz = CellsBlock(Io)%Zc - CellsBlock(Is)%Zc
    
                Dist = sqrt(Distx*Distx+Disty*Disty+Distz*Distz)
    
                Term1 = exp(J*K_air*Dist)/Dist**2.
                Term2 = J*K_air-1/Dist
                Term3 = CellsBlock(Is)%parameter_Const*CellsBlock(Is)%parameter_Ce
            
                if (num_dir == 1) then  !Iox
                    !f
                    Fxx = (Term2 + K_air**2.*(Dist - Distx**2./Dist)-3.*Term2*Distx**2./Dist**2.)*Term1
                    Fxy = (Distx*Disty*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
                    Fxz = (Distx*Distz*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
                    !g
                    g(1) = -Fxx*Term3 !XX :Green_s_tr(Iox,Isx)
                    g(2) = -Fxy*Term3 !XY : Green_s_tr(Iox,Isy)
                    g(3) = -Fxz*Term3 !XZ : Green_s_tr(Iox,Isz)
                elseif (num_dir == 2) then !Ioy
                    !f
                    Fxy = (Distx*Disty*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
                    Fyy = (Term2 + K_air**2.*(Dist - Disty**2./Dist)-3.*Term2*Disty**2./Dist**2.)*Term1
                    Fyz = (Disty*Distz*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
                    !g
                    g(1) = -Fxy*Term3 !YX :  Green_s_tr(Ioy,Isx)
                    g(2) = -Fyy*Term3 !YY : Green_s_tr(Ioy,Isy)
                    g(3) = -Fyz*Term3 !YZ : Green_s_tr(Ioy,Isz)
                else !Ioz
                    !f
                    Fxz = (Distx*Distz*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
                    Fyz = (Disty*Distz*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
                    Fzz = (Term2 + K_air**2.*(Dist - Distz**2./Dist)-3.*Term2*Distz**2./Dist**2.)*Term1
                    !g
                    g(1) = -Fxz*Term3 !ZX : Green_s_tr(Ioz,Isx)
                    g(2) = -Fyz*Term3 !ZY : Green_s_tr(Ioz,Isy)
                    g(3) = -Fzz*Term3 !ZZ : Green_s_tr(Ioz,Isz)
                endif
    
                abs_g = abs(g); 
            
                If (abs_g(1) .gt. threshold)  then
                    curs_nnz = curs_nnz + 1; icol(curs_nnz) = Isx; 
                    v_irow = v_irow + 1;
                    Green_s_tr(curs_nnz)= g(1);
                endif 
                If (abs_g(2) .gt. threshold)  then
                    curs_nnz = curs_nnz + 1; icol(curs_nnz) = Isy; 
                    v_irow = v_irow + 1;
                    Green_s_tr(curs_nnz)= g(2);  
                endif          
                If (abs_g(3) .gt. threshold)  then 
                    curs_nnz = curs_nnz + 1; icol(curs_nnz) = Isz; 
                    v_irow = v_irow + 1;
                    Green_s_tr(curs_nnz)= g(3); 
                endif      
                EndDo
            Enddo
        EndDo  
    EndIf
        
    irow(3*sizeBlock+1) = v_irow;
    
    if (nnz .ne. 0) then ! nnz is ne 0 when SR_Green_s_tr_partial is called to calculate the CBFs, 
                         ! it is equal to 0 when the subroutine 
                         ! is called to determine fSR 
      if (curs_nnz .ne. nnz) then 
          Write(*,'(a)') 'Something went wrong when filling ZGreen_s_tr: curs_nnz .ne. nnz!'
          stop 1;
      endif
    endif
    
        
    End Subroutine SR_Green_s_tr_partial
    
SUBROUTINE SR_Green_s_tr_partial_FN(sizeBlock,CellsBlock,localfSR,nnz,FN_SR_Green_s_tr)

    USE Initialization
    USE common_variables
    Implicit none
    
    !IN/OUT
    Integer, INTENT(IN) :: sizeBlock
    real(kind=8), INTENT(IN) :: localfSR
    type (Cell), Dimension(sizeBlock),INTENT(IN) :: CellsBlock
    Integer, INTENT(OUT) :: nnz
    real(kind=8), INTENT(OUT) :: FN_SR_Green_s_tr
    
    ! Local    
    Integer :: curs_nnz,ii,jj,Is, Io, Isg, Iog, Iox, Ioy, Ioz, Isx, Isy, Isz
    Integer :: Ioo,num_dir
    Real (kind=8) :: Distx, Disty, Distz, Dist
    Real (kind=8) :: v_max,Green_s_tr_max,v_irow,threshold
    Real (kind=8),dimension(3) :: abs_g
    Complex	:: Term1, Term2, Term3, Fxx, Fyy, Fzz, Fxy, Fyz, Fxz
    Real(kind=8) :: G_elts
    complex,dimension(3) :: g        
    
    nnz = 0;
    G_elts = 0;
    
    ! start by computing Green_s_tr(1,1), it will be used to elliminate 
    ! what we will consider as weak/non-significant interactions
    v_max = abs(1 - CellsBlock(1)%parameter_Sing)
    threshold = v_max/localfSR; !fct_SR;  ! the global fct_SR can be used if not all the blocks tested   
    
    If (homogs == 1) Then     
      Do Io = 1, sizeBlock   ! loop on row
        
        Do num_dir=1,3            
          ! Is = Io 
          Ioo = 3*(Io-1) + num_dir
                    
          !IoIo dir-dir (Io == Is)
          nnz = nnz + 1; 
          G_elts = G_elts + (abs(1 - CellsBlock(Io)%parameter_Sing))**2.;  ;     ! + (abs(Zpatch_e_spr(cc)))**2.           
                    
          ! Now Is > Io
          Do Is=Io+1, sizeBlock     !loop on col  ! if homogs=0 we will scan the entire matrix    
            Isx = 3*(Is-1)+1
            Isy = 3*(Is-1)+2
            Isz = 3*(Is-1)+3
            
            Distx = CellsBlock(Io)%Xc - CellsBlock(Is)%Xc
            Disty = CellsBlock(Io)%Yc - CellsBlock(Is)%Yc
            Distz = CellsBlock(Io)%Zc - CellsBlock(Is)%Zc
    
            Dist = sqrt(Distx*Distx+Disty*Disty+Distz*Distz)
    
            Term1 = exp(J*K_air*Dist)/Dist**2.
            Term2 = J*K_air-1/Dist
            Term3 = CellsBlock(Is)%parameter_Const*CellsBlock(Is)%parameter_Ce
            
            if (num_dir == 1) then  !Iox
              !f
              Fxx = (Term2 + K_air**2.*(Dist - Distx**2./Dist)-3.*Term2*Distx**2./Dist**2.)*Term1
              Fxy = (Distx*Disty*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
              Fxz = (Distx*Distz*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
              !g
              g(1) = -Fxx*Term3 !XX :Green_s_tr(Iox,Isx)
              g(2) = -Fxy*Term3 !XY : Green_s_tr(Iox,Isy)
              g(3) = -Fxz*Term3 !XZ : Green_s_tr(Iox,Isz)
            elseif (num_dir == 2) then !Ioy
              !f
              Fxy = (Distx*Disty*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
              Fyy = (Term2 + K_air**2.*(Dist - Disty**2./Dist)-3.*Term2*Disty**2./Dist**2.)*Term1
              Fyz = (Disty*Distz*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
              !g
              g(1) = -Fxy*Term3 !YX :  Green_s_tr(Ioy,Isx)
              g(2) = -Fyy*Term3 !YY : Green_s_tr(Ioy,Isy)
              g(3) = -Fyz*Term3 !YZ : Green_s_tr(Ioy,Isz)
            else !Ioz
              !f
              Fxz = (Distx*Distz*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
              Fyz = (Disty*Distz*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
              Fzz = (Term2 + K_air**2.*(Dist - Distz**2./Dist)-3.*Term2*Distz**2./Dist**2.)*Term1
              !g
              g(1) = -Fxz*Term3 !ZX : Green_s_tr(Ioz,Isx)
              g(2) = -Fyz*Term3 !ZY : Green_s_tr(Ioz,Isy)
              g(3) = -Fzz*Term3 !ZZ : Green_s_tr(Ioz,Isz)
            endif
    
            abs_g = abs(g); 
            
            If (abs_g(1) .gt. threshold)  then
              nnz = nnz + 1; 
              G_elts= G_elts + abs(g(1))**2;
            endif 
            If (abs_g(2) .gt. threshold)  then
              nnz = nnz + 1;
              G_elts = G_elts + abs(g(2))**2;  
            endif          
            If (abs_g(3) .gt. threshold)  then 
              nnz = nnz + 1; 
              G_elts= G_elts + abs(g(3))**2; 
            endif   
          EndDo
        Enddo
      EndDo 
    Else    
        Do Io = 1, sizeBlock   ! loop on row        
            Do num_dir=1,3  
          
                Ioo = 3*(Io-1) + num_dir
                
                ! here Is < Io
                Do Is=1, Io-1       
                Isx = 3*(Is-1)+1
                Isy = 3*(Is-1)+2
                Isz = 3*(Is-1)+3
            
                Distx = CellsBlock(Io)%Xc - CellsBlock(Is)%Xc
                Disty = CellsBlock(Io)%Yc - CellsBlock(Is)%Yc
                Distz = CellsBlock(Io)%Zc - CellsBlock(Is)%Zc
    
                Dist = sqrt(Distx*Distx+Disty*Disty+Distz*Distz)
    
                Term1 = exp(J*K_air*Dist)/Dist**2.
                Term2 = J*K_air-1/Dist
                Term3 = CellsBlock(Is)%parameter_Const*CellsBlock(Is)%parameter_Ce
            
                if (num_dir == 1) then  !Iox
                    !f
                    Fxx = (Term2 + K_air**2.*(Dist - Distx**2./Dist)-3.*Term2*Distx**2./Dist**2.)*Term1
                    Fxy = (Distx*Disty*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
                    Fxz = (Distx*Distz*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
                    !g
                    g(1) = -Fxx*Term3 !XX :Green_s_tr(Iox,Isx)
                    g(2) = -Fxy*Term3 !XY : Green_s_tr(Iox,Isy)
                    g(3) = -Fxz*Term3 !XZ : Green_s_tr(Iox,Isz)
                elseif (num_dir == 2) then !Ioy
                    !f
                    Fxy = (Distx*Disty*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
                    Fyy = (Term2 + K_air**2.*(Dist - Disty**2./Dist)-3.*Term2*Disty**2./Dist**2.)*Term1
                    Fyz = (Disty*Distz*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
                    !g
                    g(1) = -Fxy*Term3 !YX :  Green_s_tr(Ioy,Isx)
                    g(2) = -Fyy*Term3 !YY : Green_s_tr(Ioy,Isy)
                    g(3) = -Fyz*Term3 !YZ : Green_s_tr(Ioy,Isz)
                else !Ioz
                    !f
                    Fxz = (Distx*Distz*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
                    Fyz = (Disty*Distz*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
                    Fzz = (Term2 + K_air**2.*(Dist - Distz**2./Dist)-3.*Term2*Distz**2./Dist**2.)*Term1
                    !g
                    g(1) = -Fxz*Term3 !ZX : Green_s_tr(Ioz,Isx)
                    g(2) = -Fyz*Term3 !ZY : Green_s_tr(Ioz,Isy)
                    g(3) = -Fzz*Term3 !ZZ : Green_s_tr(Ioz,Isz)
                endif
    
                abs_g = abs(g); 
            
                If (abs_g(1) .gt. threshold)  then
                    nnz = nnz + 1; 
                    G_elts = G_elts + abs(g(1))**2;
                endif 
                If (abs_g(2) .gt. threshold)  then
                    nnz = nnz + 1; 
                    G_elts = G_elts + abs(g(2))**2;  
                endif          
                If (abs_g(3) .gt. threshold)  then 
                    nnz = nnz + 1; 
                    G_elts = G_elts + abs(g(3))**2; 
                endif      
                EndDo
        
                ! Is = Io
                !IoIo dir-dir (Io == Is)
                nnz = nnz + 1; 
                G_elts = G_elts + abs((1 - CellsBlock(Io)%parameter_Sing))**2 ; 
          
                ! Now Is > Io
                Do Is=Io+1, sizeBlock     !loop on col  ! if homogs=0 we will scan the entire matrix    
                Isx = 3*(Is-1)+1
                Isy = 3*(Is-1)+2
                Isz = 3*(Is-1)+3
            
                Distx = CellsBlock(Io)%Xc - CellsBlock(Is)%Xc
                Disty = CellsBlock(Io)%Yc - CellsBlock(Is)%Yc
                Distz = CellsBlock(Io)%Zc - CellsBlock(Is)%Zc
    
                Dist = sqrt(Distx*Distx+Disty*Disty+Distz*Distz)
    
                Term1 = exp(J*K_air*Dist)/Dist**2.
                Term2 = J*K_air-1/Dist
                Term3 = CellsBlock(Is)%parameter_Const*CellsBlock(Is)%parameter_Ce
            
                if (num_dir == 1) then  !Iox
                    !f
                    Fxx = (Term2 + K_air**2.*(Dist - Distx**2./Dist)-3.*Term2*Distx**2./Dist**2.)*Term1
                    Fxy = (Distx*Disty*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
                    Fxz = (Distx*Distz*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
                    !g
                    g(1) = -Fxx*Term3 !XX :Green_s_tr(Iox,Isx)
                    g(2) = -Fxy*Term3 !XY : Green_s_tr(Iox,Isy)
                    g(3) = -Fxz*Term3 !XZ : Green_s_tr(Iox,Isz)
                elseif (num_dir == 2) then !Ioy
                    !f
                    Fxy = (Distx*Disty*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
                    Fyy = (Term2 + K_air**2.*(Dist - Disty**2./Dist)-3.*Term2*Disty**2./Dist**2.)*Term1
                    Fyz = (Disty*Distz*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
                    !g
                    g(1) = -Fxy*Term3 !YX :  Green_s_tr(Ioy,Isx)
                    g(2) = -Fyy*Term3 !YY : Green_s_tr(Ioy,Isy)
                    g(3) = -Fyz*Term3 !YZ : Green_s_tr(Ioy,Isz)
                else !Ioz
                    !f
                    Fxz = (Distx*Distz*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
                    Fyz = (Disty*Distz*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
                    Fzz = (Term2 + K_air**2.*(Dist - Distz**2./Dist)-3.*Term2*Distz**2./Dist**2.)*Term1
                    !g
                    g(1) = -Fxz*Term3 !ZX : Green_s_tr(Ioz,Isx)
                    g(2) = -Fyz*Term3 !ZY : Green_s_tr(Ioz,Isy)
                    g(3) = -Fzz*Term3 !ZZ : Green_s_tr(Ioz,Isz)
                endif
    
                abs_g = abs(g); 
            
                If (abs_g(1) .gt. threshold)  then
                    nnz = nnz + 1; 
                    G_elts = G_elts + abs(g(1))**2;
                endif 
                If (abs_g(2) .gt. threshold)  then
                    nnz = nnz + 1; 
                    G_elts = G_elts + abs(g(2))**2;  
                endif          
                If (abs_g(3) .gt. threshold)  then 
                    nnz = nnz + 1; 
                    G_elts= G_elts + abs(g(3))**2; 
                endif      
                EndDo
            Enddo
        EndDo  
    EndIf  
    
    FN_SR_Green_s_tr = sqrt(G_elts)
        
    End Subroutine SR_Green_s_tr_partial_FN
    
    



    


