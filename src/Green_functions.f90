!! SUBROUTINES : 
!! - Green_s_tr_total
!! - Green_s_tr_partial
!! - Green_s_tr_partial_FN
!! - SR_Green_s_tr_partial
!! - SR_Green_s_tr_partial_FN
!! - DR_Green_s_tr_partial 
!! - Green_ff_dt
!! - computeBlockCol.f90
!! - computeBlockCol_SMW.f90
!! - computeBlockRow.f90
!! - computeBlockRow_SMW.f90


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


SUBROUTINE Green_ff_dt(Nc,Cells_in,Receivers,num_capteur,theta_capteur,phi_capteur,Green_dt)

    USE Initialization
    USE common_variables
    IMPLICIT NONE
    
    Integer, INTENT(IN) :: Nc
    type (Cell), Dimension(Nc), INTENT(IN) :: Cells_in
    type (Dipole), Dimension(NRx_tot), INTENT(IN) :: Receivers
    Integer, INTENT(IN) :: num_capteur
    COMPLEX(real64), Dimension(3,3*Nc), INTENT(OUT) :: Green_dt
    Real(kind=8), INTENT(IN) :: theta_capteur,phi_capteur
    
    !Local
    Integer Is, Isx, Isy, Isz
    Real(kind=8) :: xc,yc,zc,x_cap,y_cap,z_cap
    COMPLEX(real64) :: Term1, Term2, Term3, Fxx, Fyy, Fzz, Fxy, Fyz, Fxz

    DO Is=1, Nc
        Isx=3*(Is-1)+1
        Isy=3*(Is-1)+2
        Isz=3*(Is-1)+3
        
        xc = Cells_in(Is)%Xc;
        yc = Cells_in(Is)%Yc;
        zc = Cells_in(Is)%Zc;
        
        !! BSA (Back Scattering Alignment) (see equations 1.65 for FSA and 
        !!transformation to BSA in Eq 1.68 Phd Bellez : ksb = -ksf)
        !!x_cap = -sin(theta_capteur*Pi/180.)*cos(phi_capteur*Pi/180.)
        !!y_cap = -sin(theta_capteur*Pi/180.)*sin(phi_capteur*Pi/180.)
        !!z_cap = -cos(theta_capteur*Pi/180.)
        
        !! since 3/13/2019 !
        x_cap = cos(theta_capteur*Pi/180.) 
        y_cap = sin(theta_capteur*Pi/180.)*cos(phi_capteur*Pi/180.) 
        z_cap = sin(theta_capteur*Pi/180.)*sin(phi_capteur*Pi/180.)
        
        ! translation theorem E2= E1*exp(-ik*delta.u) 
        Term1= Cells_in(Is)%parameter_Const*exp(-J*K_air*(x_cap*xc+y_cap*yc+z_cap*zc))
        Term3= Cells_in(Is)%parameter_Ce
        
        ! if simplified Fxx=Fyy=Fzz=Cells(Is)%parameter_Const*K_air**2.*exp(J*K_air*Dist)/Dist 
        Fxx= Term1*K_air**2. * K_air; 
        Fyy= Term1*K_air**2. * K_air;
        Fzz= Term1*K_air**2. * K_air;
        
        Fxy=0
        Fxz=0
        Fyz=0
        
        Green_dt(1,Isx)=Fxx*Term3 !!XX
        Green_dt(1,Isy)=Fxy*Term3 !!XY
        Green_dt(1,Isz)=Fxz*Term3 !!XZ
        
        Green_dt(2,Isx)=Fxy*Term3 !!YX
        Green_dt(2,Isy)=Fyy*Term3 !!YY
        Green_dt(2,Isz)=Fyz*Term3 !!YZ

        Green_dt(3,Isx)=Fxz*Term3 !!ZX
        Green_dt(3,Isy)=Fyz*Term3 !!ZY
        Green_dt(3,Isz)=Fzz*Term3 !!ZZ        
        
    ENDDO

END SUBROUTINE Green_ff_dt


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
 
 
 SUBROUTINE computeBlockRow(Cells,Matrix_Green_irow,Matrix_Green_IndexCol_inf,nbCols,Matrix_Green_row) 
  
    USE Initialization
    USE common_variables
    
    IMPLICIT NONE

    !! IN/OUT ******************************************************************
    type (Cell), Dimension(Nbc), INTENT(IN) :: Cells
    INTEGER, INTENT(IN) :: Matrix_Green_irow, Matrix_Green_IndexCol_inf,nbCols
    COMPLEX(KIND=8), DIMENSION(1,nbCols), INTENT(OUT) :: Matrix_Green_row
  
    !local
    INTEGER jj,Is,Isc,Io,Ioc
    INTEGER Index_Col, Is_init, Nb_cels_j
    COMPLEX :: Term1, Term2, Term3
    COMPLEX :: Fxx, Fyy, Fzz, Fxy, Fyz,Fxz,Fyx,Fzy,Fzx,Gx,Gy,Gz
    Complex :: K11x, K11z, K22z, K12z, K21z
    Real(kind=8) :: RE, RM, Ang_inc, P1, P2, norme, Px, Py
    Real(kind=8) :: Distx, Disty, Distz, Dist, Distxy
  
    
    Ioc= mod(Matrix_Green_irow,3)
    If (Ioc .NE. 0) Then
        Io = Matrix_Green_irow/3 + 1
    Else
        Io = Matrix_Green_irow/3
    Endif
    Index_Col = 0
    Is_init = Matrix_Green_IndexCol_inf/3 + 1
    nb_cels_j = nbCols/3  
    Is = Is_init
  
    If (Ioc ==1) Then !! X.
        Do jj=1,nb_cels_j      
            
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

            Matrix_Green_row(1,Index_Col+1)= - Fxx*Term3 !XX
            Matrix_Green_row(1,Index_Col+2)= - Fxy*Term3 !XY
            Matrix_Green_row(1,Index_Col+3)= - Fxz*Term3 !XZ             
                    
            Index_Col = Index_Col +3
            Is = Is + 1
        Enddo
    ElseIf (Ioc == 2) Then  !!Y.
        Do jj=1,nb_cels_j      
            
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
            
            Matrix_Green_row(1,Index_Col+1)= - Fxy*Term3   !YX
            Matrix_Green_row(1,Index_Col+2)= - Fyy*Term3   !YY
            Matrix_Green_row(1,Index_Col+3)= - Fyz*Term3   !YZ
          
            Index_Col = Index_Col +3
            Is = Is + 1
        Enddo
    Else !!Z.
        Do jj=1,nb_cels_j      
            
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

            Matrix_Green_row(1,Index_Col+1)= - Fxz*Term3 !ZX
            Matrix_Green_row(1,Index_Col+2)= - Fyz*Term3 !ZY
            Matrix_Green_row(1,Index_Col+3)= - Fzz*Term3 !ZZ
          
            Index_Col = Index_Col +3
            Is = Is + 1
        Enddo
    EndIf  
 END SUBROUTINE computeBlockRow
 
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