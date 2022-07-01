!! SUBROUTINES : 
!! 1) To compute ZMoM matrix elements : 
!! - Green_s_tr_total
!! - Green_s_tr_partial
!! - Green_s_tr_partial_FN
!! - SR_Green_s_tr_partial
!! - SR_Green_s_tr_partial_FN
!! - DR_Green_s_tr_partial

!! 2) To compute ZMoM Col/Row elements (for Adaptive Cross Approximation and Sherman-Morrison Woodebery)
!! - computeBlockCol.f90
!! - computeBlockCol_SMW.f90
!! - computeBlockRow.f90
!! - computeBlockRow_SMW.f90
 
!! 3) To compute Scattered fields/ Scattering matrices 
!! - Green_s_dt
!! - Green_ff_dt


SUBROUTINE Green_s_tr_total(Cells,Green_s_tr)

    USE Initialization
    USE common_variables
    Implicit none

    ! IN/OUT
    type (Cell), Dimension(Nbc), INTENT(IN) :: Cells
    COMPLEX(real64), Dimension(3*Nbc,3*Nbc), INTENT(OUT):: Green_s_tr
    
    ! Local 
    Complex :: Gr_mn, Tau_mn, f_kapChe, Fxx, Fyy, Fzz, Fxy, Fyz, Fxz
    Real(kind=8) :: rx, ry, rz, r_mn
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
          Green_s_tr(Iox,Isx) = 1 - Cells(Is)%Znnpp
          Green_s_tr(Ioy,Isy) = 1 - Cells(Is)%Znnpp
          Green_s_tr(Ioz,Isz) = 1 - Cells(Is)%Znnpp
      
          Green_s_tr(Iox,Isy) = 0
          Green_s_tr(Iox,Isz) = 0
          Green_s_tr(Ioy,Isx) = 0
          Green_s_tr(Ioz,Isx) = 0
          Green_s_tr(Ioy,Isz) = 0
          Green_s_tr(Ioz,Isy) = 0    
	    else
          rx = Cells(Io)%Xc - Cells(Is)%Xc
          ry = Cells(Io)%Yc - Cells(Is)%Yc
          rz = Cells(Io)%Zc - Cells(Is)%Zc

          r_mn = sqrt(rx**2.+ry**2.+rz**2.)   ! r_mn

          Gr_mn = exp(J*k_0*r_mn)/(4*Pi*r_mn**2.)      ! Green_mn
          Tau_mn = J*k_0 - 1/r_mn                      ! Tau_mn 
          f_kapChe = Cells(Is)%Kappa_n*Cells(Is)%Che_n ! factor : Kappa_n * Che_n
          
          ! if p .eq.  q
          Fxx = Gr_mn * (Tau_mn + k_0**2.*(r_mn - rx**2./r_mn)-3.*rx**2./r_mn**2. * Tau_mn)
          Fyy = Gr_mn * (Tau_mn + k_0**2.*(r_mn - ry**2./r_mn)-3.*ry**2./r_mn**2. * Tau_mn)
          Fzz = Gr_mn * (Tau_mn + k_0**2.*(r_mn - rz**2./r_mn)-3.*rz**2./r_mn**2. * Tau_mn)

          ! if p .ne. q
          Fxy = Gr_mn/r_mn * (rx*ry) * (-k_0**2. -3.*Tau_mn/r_mn)
          Fyz = Gr_mn/r_mn * (ry*rz) * (-k_0**2. -3.*Tau_mn/r_mn)
          Fxz = Gr_mn/r_mn * (rx*rz) * (-k_0**2. -3.*Tau_mn/r_mn)
          
          Green_s_tr(Iox,Isx)= -Fxx*f_kapChe !XX
          Green_s_tr(Iox,Isy)= -Fxy*f_kapChe !XY
          Green_s_tr(Iox,Isz)= -Fxz*f_kapChe !XZ
          Green_s_tr(Ioy,Isx)= -Fxy*f_kapChe !YX
          Green_s_tr(Ioy,Isy)= -Fyy*f_kapChe !YY
          Green_s_tr(Ioy,Isz)= -Fyz*f_kapChe !YZ
          Green_s_tr(Ioz,Isx)= -Fxz*f_kapChe !ZX
          Green_s_tr(Ioz,Isy)= -Fyz*f_kapChe !ZY
          Green_s_tr(Ioz,Isz)= -Fzz*f_kapChe !ZZ                 
        endIf    
      endDo
    endDo 
    
End Subroutine Green_s_tr_total

SUBROUTINE Green_s_tr_partial(sizeB1,CellsB1,sizeB2,CellsB2,Green_s_tr)

    USE Initialization
    USE common_variables
    Implicit none
    
    !IN/OUT
    Integer, INTENT(IN) :: sizeB1,sizeB2
    type (Cell), Dimension(sizeB1), INTENT(IN) :: CellsB1
    type (Cell), Dimension(sizeB2), INTENT(IN) :: CellsB2
    COMPLEX(real64), Dimension(3*sizeB1,3*sizeB2),INTENT(OUT)::Green_s_tr

    ! Local
    Complex	:: Gr_mn, Tau_mn, f_kapChe, Fxx, Fyy, Fzz, Fxy, Fyz, Fxz
    Real (kind=8) :: rx, ry, rz, r_mn
    Integer :: Is, Io, Isg, Iog, Iox, Ioy, Ioz, Isx, Isy, Isz           
    

    Do Is=1, sizeB2
      Do Io=1, sizeB1  
      
        Iox = 3*(Io-1)+1
        Ioy = 3*(Io-1)+2
        Ioz = 3*(Io-1)+3

        Isx = 3*(Is-1)+1
        Isy = 3*(Is-1)+2
        Isz = 3*(Is-1)+3
    
        Iog = CellsB1(Io)%n_cell
        Isg = CellsB2(Is)%n_cell
        if (Iog == Isg) then
          Green_s_tr(Iox,Isx) = 1 - CellsB2(Is)%Znnpp
          Green_s_tr(Ioy,Isy) = 1 - CellsB2(Is)%Znnpp
          Green_s_tr(Ioz,Isz) = 1 - CellsB2(Is)%Znnpp
      
          Green_s_tr(Iox,Isy) = 0
          Green_s_tr(Iox,Isz) = 0
          Green_s_tr(Ioy,Isx) = 0
          Green_s_tr(Ioz,Isx) = 0
          Green_s_tr(Ioy,Isz) = 0
          Green_s_tr(Ioz,Isy) = 0    
	    else
          rx = CellsB1(Io)%Xc - CellsB2(Is)%Xc
          ry = CellsB1(Io)%Yc - CellsB2(Is)%Yc
          rz = CellsB1(Io)%Zc - CellsB2(Is)%Zc

          r_mn = sqrt(rx**2.+ry**2.+rz**2.)   ! r_mn 

          Gr_mn = exp(J*k_0*r_mn)/(4*Pi*r_mn**2.)  ! Green_mn
          Tau_mn = J*k_0 - 1/r_mn                    ! Tau_mn 
          f_kapChe = CellsB2(Is)%Kappa_n*CellsB2(Is)%Che_n  ! factor : Kappa_n * Che_n

          ! if p .eq.  q
          Fxx = Gr_mn * (Tau_mn + k_0**2.*(r_mn - rx**2./r_mn)-3.*rx**2./r_mn**2. * Tau_mn)
          Fyy = Gr_mn * (Tau_mn + k_0**2.*(r_mn - ry**2./r_mn)-3.*ry**2./r_mn**2. * Tau_mn)
          Fzz = Gr_mn * (Tau_mn + k_0**2.*(r_mn - rz**2./r_mn)-3.*rz**2./r_mn**2. * Tau_mn)

          ! if p .ne. q
          Fxy = Gr_mn/r_mn * (rx*ry) * (-k_0**2. -3.*Tau_mn/r_mn)
          Fyz = Gr_mn/r_mn * (ry*rz) * (-k_0**2. -3.*Tau_mn/r_mn)
          Fxz = Gr_mn/r_mn * (rx*rz) * (-k_0**2. -3.*Tau_mn/r_mn)
      
          Green_s_tr(Iox,Isx)= -Fxx*f_kapChe !XX
          Green_s_tr(Iox,Isy)= -Fxy*f_kapChe !XY
          Green_s_tr(Iox,Isz)= -Fxz*f_kapChe !XZ
          Green_s_tr(Ioy,Isx)= -Fxy*f_kapChe !YX
          Green_s_tr(Ioy,Isy)= -Fyy*f_kapChe !YY
          Green_s_tr(Ioy,Isz)= -Fyz*f_kapChe !YZ
          Green_s_tr(Ioz,Isx)= -Fxz*f_kapChe !ZX
          Green_s_tr(Ioz,Isy)= -Fyz*f_kapChe !ZY
          Green_s_tr(Ioz,Isz)= -Fzz*f_kapChe !ZZ
                 
        endIf    
      endDo
    endDo  
    
    End Subroutine Green_s_tr_partial
    
    
    SUBROUTINE Green_s_tr_partial_FN(sizeB1,CellsB1,sizeB2,CellsB2,FN_Green_s_tr)

    ! this subroutine will simply enable us to determine the Frobenius norm without storing Z patch (memory efficient)
    USE Initialization
    USE common_variables
    USE iso_fortran_env
    
    Implicit none
    
    !IN/OUT
    Integer, INTENT(IN) :: sizeB1,sizeB2
    type (Cell), Dimension(sizeB1), INTENT(IN) :: CellsB1
    type (Cell), Dimension(sizeB2), INTENT(IN) :: CellsB2
    Real(kind=8), INTENT(OUT) :: FN_Green_s_tr

    ! Local
    Complex	:: Gr_mn, Tau_mn, f_kapChe, Fxx, Fyy, Fzz, Fxy, Fyz, Fxz
    Real (kind=8) :: rx, ry, rz, r_mn,val_fn
    Integer :: Is, Io, Isg, Iog, Iox, Ioy, Ioz, Isx, Isy, Isz,Io_beg
        
    
    val_fn = 0D0;
  
    if (homogs .eq. 0) then 
        Do Is=1, sizeB2
          Do Io=1, sizeB1  
      
            Iox = 3*(Io-1)+1
            Ioy = 3*(Io-1)+2
            Ioz = 3*(Io-1)+3

            Isx = 3*(Is-1)+1
            Isy = 3*(Is-1)+2
            Isz = 3*(Is-1)+3
    
            Iog = CellsB1(Io)%n_cell
            Isg = CellsB2(Is)%n_cell
            if (Iog == Isg) then
              val_fn = val_fn + 3.*(abs(1 - CellsB2(Is)%Znnpp))**2          
	        else
              rx = CellsB1(Io)%Xc - CellsB2(Is)%Xc
              ry = CellsB1(Io)%Yc - CellsB2(Is)%Yc
              rz = CellsB1(Io)%Zc - CellsB2(Is)%Zc

              r_mn = sqrt(rx**2.+ry**2.+rz**2.)

              Gr_mn = exp(J*k_0*r_mn)/(4*Pi*r_mn**2.)
              Tau_mn = J*k_0 - 1/r_mn
              f_kapChe = CellsB2(Is)%Kappa_n*CellsB2(Is)%Che_n

              ! if p .eq.  q
              Fxx = Gr_mn * (Tau_mn + k_0**2.*(r_mn - rx**2./r_mn)-3.*rx**2./r_mn**2. * Tau_mn)
              Fyy = Gr_mn * (Tau_mn + k_0**2.*(r_mn - ry**2./r_mn)-3.*ry**2./r_mn**2. * Tau_mn)
              Fzz = Gr_mn * (Tau_mn + k_0**2.*(r_mn - rz**2./r_mn)-3.*rz**2./r_mn**2. * Tau_mn)
    
              ! if p .ne. q
              Fxy = Gr_mn/r_mn * (rx*ry) * (-k_0**2. -3.*Tau_mn/r_mn)
              Fyz = Gr_mn/r_mn * (ry*rz) * (-k_0**2. -3.*Tau_mn/r_mn)
              Fxz = Gr_mn/r_mn * (rx*rz) * (-k_0**2. -3.*Tau_mn/r_mn)
          
              val_fn = val_fn + ((abs(Fxx))**2+2.*(abs(Fxy))**2+2.*(abs(Fxz))**2 &
                                +(abs(Fyy))**2+2.*(abs(Fyz))**2+(abs(Fzz))**2)*(abs(f_kapChe))**2. 
            endIf    
          endDo
        endDo 
    else
        Do Is=1, sizeB2
          Do Io=Is, sizeB1  
      
            Iox = 3*(Io-1)+1
            Ioy = 3*(Io-1)+2
            Ioz = 3*(Io-1)+3

            Isx = 3*(Is-1)+1
            Isy = 3*(Is-1)+2
            Isz = 3*(Is-1)+3
    
            Iog = CellsB1(Io)%n_cell
            Isg = CellsB2(Is)%n_cell
            if (Iog == Isg) then
              val_fn = val_fn + 3.*(abs(1 - CellsB2(Is)%Znnpp))**2          
	        else
              rx = CellsB1(Io)%Xc - CellsB2(Is)%Xc
              ry = CellsB1(Io)%Yc - CellsB2(Is)%Yc
              rz = CellsB1(Io)%Zc - CellsB2(Is)%Zc

              r_mn = sqrt(rx**2.+ry**2.+rz**2.)

              Gr_mn = exp(J*k_0*r_mn)/(4*Pi*r_mn**2.)
              Tau_mn = J*k_0 - 1/r_mn
              f_kapChe = CellsB2(Is)%Kappa_n*CellsB2(Is)%Che_n

              ! if p .eq.  q
              Fxx = Gr_mn * (Tau_mn + k_0**2.*(r_mn - rx**2./r_mn)-3.*rx**2./r_mn**2. * Tau_mn)
              Fyy = Gr_mn * (Tau_mn + k_0**2.*(r_mn - ry**2./r_mn)-3.*ry**2./r_mn**2. * Tau_mn)
              Fzz = Gr_mn * (Tau_mn + k_0**2.*(r_mn - rz**2./r_mn)-3.*rz**2./r_mn**2. * Tau_mn)
    
              ! if p .ne. q
              Fxy = Gr_mn/r_mn * (rx*ry) * (-k_0**2. -3.*Tau_mn/r_mn)
              Fyz = Gr_mn/r_mn * (ry*rz) * (-k_0**2. -3.*Tau_mn/r_mn)
              Fxz = Gr_mn/r_mn * (rx*rz) * (-k_0**2. -3.*Tau_mn/r_mn)
          
              val_fn = val_fn + ((abs(Fxx))**2+2.*(abs(Fxy))**2+2.*(abs(Fxz))**2 &
                                +(abs(Fyy))**2+2.*(abs(Fyz))**2+(abs(Fzz))**2)*(abs(f_kapChe))**2. 
            endIf    
          endDo
        endDo 
    endif      

    FN_Green_s_tr = sqrt(val_fn);

End Subroutine Green_s_tr_partial_FN


SUBROUTINE SR_Green_s_tr_partial(sizeB,CellsB,localfSR,nnz,Green_s_tr,irow,icol)           

    USE Initialization
    USE common_variables
    Implicit none
    
    !IN/OUT
    Integer, INTENT(IN) :: sizeB,nnz
    real(kind=8), INTENT(IN) :: localfSR
    Integer, Dimension(3*sizeB+1) :: irow
    Integer, Dimension(nnz) :: icol
    type (Cell), Dimension(sizeB),INTENT(IN) :: CellsB
    COMPLEX(real64),Dimension(nnz) ::Green_s_tr

    ! Local    
    Integer :: curs_nnz,ii,jj,Is, Io, Isg, Iog, Iox, Ioy, Ioz, Isx, Isy, Isz
    Integer :: Ioo,num_dir
    Real (kind=8) :: rx, ry, rz, r_mn
    Real (kind=8) :: v_max,Green_s_tr_max,v_irow,threshold
    Real (kind=8),dimension(3) :: abs_g
    Complex	:: Gr_mn, Tau_mn, f_kapChe, Fxx, Fyy, Fzz, Fxy, Fyz, Fxz
    complex,dimension(3) :: g        
    
    Green_s_tr_max =0D0
    curs_nnz = 0;
    v_irow = 1;
    
    ! start by computing Green_s_tr(1,1), it will be used to elliminate 
    ! what we will consider as weak/non-significant interactions
    v_max = abs(1 - CellsB(1)%Znnpp)
    threshold = v_max/localfSR; !fct_SR;  ! the global fct_SR can be used if not all the blocks tested   
    
    If (homogs == 1) Then     
      Do Io = 1, sizeB   ! loop on row
        
        Do num_dir=1,3            
          ! Is = Io 
          Ioo = 3*(Io-1) + num_dir
          irow(Ioo) = v_irow ! start by storing the index (%nnz elts) 
                             ! of the first nnz elt in this row  
          
          !IoIo dir-dir (Io == Is)
          curs_nnz = curs_nnz + 1; icol(curs_nnz) = Ioo; 
          Green_s_tr(curs_nnz) = 1 - CellsB(Io)%Znnpp ;                
          v_irow = v_irow + 1 
          
          ! Now Is > Io
          Do Is=Io+1, sizeB     !loop on col  ! if homogs=0 we will scan the entire matrix    
            Isx = 3*(Is-1)+1
            Isy = 3*(Is-1)+2
            Isz = 3*(Is-1)+3
            
            rx = CellsB(Io)%Xc - CellsB(Is)%Xc
            ry = CellsB(Io)%Yc - CellsB(Is)%Yc
            rz = CellsB(Io)%Zc - CellsB(Is)%Zc
    
            r_mn = sqrt(rx**2.+ry**2.+rz**2.)
    
            Gr_mn = exp(J*k_0*r_mn)/(4*Pi*r_mn**2.)
            Tau_mn = J*k_0 - 1/r_mn
            f_kapChe = CellsB(Is)%Kappa_n*CellsB(Is)%Che_n
            
            if (num_dir == 1) then  !Iox
              !f
              Fxx = Gr_mn * (Tau_mn + k_0**2.*(r_mn - rx**2./r_mn) -3.*rx**2./r_mn**2. *Tau_mn)
              Fxy = Gr_mn/r_mn * (rx*ry) *(-k_0**2. -3.*Tau_mn/r_mn)
              Fxz = Gr_mn/r_mn * (rx*rz) *(-k_0**2. -3.*Tau_mn/r_mn)
              !g
              g(1) = -Fxx*f_kapChe !XX :Green_s_tr(Iox,Isx)
              g(2) = -Fxy*f_kapChe !XY : Green_s_tr(Iox,Isy)
              g(3) = -Fxz*f_kapChe !XZ : Green_s_tr(Iox,Isz)
            elseif (num_dir == 2) then !Ioy
              !f
              Fxy = Gr_mn/r_mn * (rx*ry) *(-k_0**2. -3.*Tau_mn/r_mn)
              Fyy = Gr_mn * (Tau_mn + k_0**2.*(r_mn - ry**2./r_mn) -3.*ry**2./r_mn**2. *Tau_mn)
              Fyz = Gr_mn/r_mn * (ry*rz) *(-k_0**2. -3.*Tau_mn/r_mn)
              !g
              g(1) = -Fxy*f_kapChe !YX :  Green_s_tr(Ioy,Isx)
              g(2) = -Fyy*f_kapChe !YY : Green_s_tr(Ioy,Isy)
              g(3) = -Fyz*f_kapChe !YZ : Green_s_tr(Ioy,Isz)
            else !Ioz
              !f
              Fxz = Gr_mn/r_mn * (rx*rz) *(-k_0**2. -3.*Tau_mn/r_mn)
              Fyz = Gr_mn/r_mn * (ry*rz) *(-k_0**2. -3.*Tau_mn/r_mn)
              Fzz = Gr_mn * (Tau_mn + k_0**2.*(r_mn - rz**2./r_mn) -3.*rz**2./r_mn**2. *Tau_mn)
              !g
              g(1) = -Fxz*f_kapChe !ZX : Green_s_tr(Ioz,Isx)
              g(2) = -Fyz*f_kapChe !ZY : Green_s_tr(Ioz,Isy)
              g(3) = -Fzz*f_kapChe !ZZ : Green_s_tr(Ioz,Isz)
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
        Do Io = 1, sizeB   ! loop on row        
            Do num_dir=1,3  
          
                Ioo = 3*(Io-1) + num_dir
                irow(Ioo) = v_irow ! start by storing the index (%nnz elts) 
                                    ! of the first nnz elt in this row
                               
                ! here Is < Io
                Do Is=1, Io-1       
                Isx = 3*(Is-1)+1
                Isy = 3*(Is-1)+2
                Isz = 3*(Is-1)+3
            
                rx = CellsB(Io)%Xc - CellsB(Is)%Xc
                ry = CellsB(Io)%Yc - CellsB(Is)%Yc
                rz = CellsB(Io)%Zc - CellsB(Is)%Zc
    
                r_mn = sqrt(rx**2.+ry**2.+rz**2.)
    
                Gr_mn = exp(J*k_0*r_mn)/(4*Pi*r_mn**2.)
                Tau_mn = J*k_0 - 1/r_mn
                f_kapChe = CellsB(Is)%Kappa_n*CellsB(Is)%Che_n
            
                if (num_dir == 1) then  !Iox
                    !f
                    Fxx = Gr_mn * (Tau_mn + k_0**2.*(r_mn - rx**2./r_mn) -3.*rx**2./r_mn**2. *Tau_mn)
                    Fxy = Gr_mn/r_mn * (rx*ry) *(-k_0**2. -3.*Tau_mn/r_mn)
                    Fxz = Gr_mn/r_mn * (rx*rz) *(-k_0**2. -3.*Tau_mn/r_mn)
                    !g
                    g(1) = -Fxx*f_kapChe !XX :Green_s_tr(Iox,Isx)
                    g(2) = -Fxy*f_kapChe !XY : Green_s_tr(Iox,Isy)
                    g(3) = -Fxz*f_kapChe !XZ : Green_s_tr(Iox,Isz)
                elseif (num_dir == 2) then !Ioy
                    !f
                    Fxy = Gr_mn/r_mn * (rx*ry) * (-k_0**2.-3.*Tau_mn/r_mn)
                    Fyy = Gr_mn * (Tau_mn + k_0**2.*(r_mn - ry**2./r_mn) -3.*ry**2./r_mn**2. *Tau_mn)
                    Fyz = Gr_mn/r_mn * (ry*rz) * (-k_0**2.-3.*Tau_mn/r_mn)
                    !g
                    g(1) = -Fxy*f_kapChe !YX :  Green_s_tr(Ioy,Isx)
                    g(2) = -Fyy*f_kapChe !YY : Green_s_tr(Ioy,Isy)
                    g(3) = -Fyz*f_kapChe !YZ : Green_s_tr(Ioy,Isz)
                else !Ioz
                    !f
                    Fxz = Gr_mn/r_mn * (rx*rz) * (-k_0**2. -3.*Tau_mn/r_mn)
                    Fyz = Gr_mn/r_mn * (ry*rz) * (-k_0**2. -3.*Tau_mn/r_mn)
                    Fzz = Gr_mn * (Tau_mn + k_0**2.*(r_mn - rz**2./r_mn) -3.*rz**2./r_mn**2. *Tau_mn)
                    !g
                    g(1) = -Fxz*f_kapChe !ZX : Green_s_tr(Ioz,Isx)
                    g(2) = -Fyz*f_kapChe !ZY : Green_s_tr(Ioz,Isy)
                    g(3) = -Fzz*f_kapChe !ZZ : Green_s_tr(Ioz,Isz)
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
                Green_s_tr(curs_nnz) = 1 - CellsB(Io)%Znnpp ;                
                v_irow = v_irow + 1 
          
                ! Now Is > Io
                Do Is=Io+1, sizeB     !loop on col  ! if homogs=0 we will scan the entire matrix    
                Isx = 3*(Is-1)+1
                Isy = 3*(Is-1)+2
                Isz = 3*(Is-1)+3
            
                rx = CellsB(Io)%Xc - CellsB(Is)%Xc
                ry = CellsB(Io)%Yc - CellsB(Is)%Yc
                rz = CellsB(Io)%Zc - CellsB(Is)%Zc
    
                r_mn = sqrt(rx**2.+ry**2.+rz**2.)
    
                Gr_mn = exp(J*k_0*r_mn)/(4*Pi*r_mn**2.)
                Tau_mn = J*k_0-1/r_mn
                f_kapChe = CellsB(Is)%Kappa_n*CellsB(Is)%Che_n
            
                if (num_dir == 1) then  !Iox
                    !f
                    Fxx = Gr_mn * (Tau_mn + k_0**2.*(r_mn - rx**2./r_mn) -3.*rx**2./r_mn**2. *Tau_mn)
                    Fxy = Gr_mn/r_mn * (rx*ry) *(-k_0**2.-3.*Tau_mn/r_mn)
                    Fxz = Gr_mn/r_mn * (rx*rz) *(-k_0**2.-3.*Tau_mn/r_mn)
                    !g
                    g(1) = -Fxx*f_kapChe !XX :Green_s_tr(Iox,Isx)
                    g(2) = -Fxy*f_kapChe !XY : Green_s_tr(Iox,Isy)
                    g(3) = -Fxz*f_kapChe !XZ : Green_s_tr(Iox,Isz)
                elseif (num_dir == 2) then !Ioy
                    !f
                    Fxy = Gr_mn/r_mn * (rx*ry) *(-k_0**2. -3.*Tau_mn/r_mn)
                    Fyy = Gr_mn * (Tau_mn + k_0**2.*(r_mn - ry**2./r_mn) -3.*ry**2./r_mn**2. *Tau_mn)
                    Fyz = Gr_mn/r_mn * (ry*rz) *(-k_0**2. -3.*Tau_mn/r_mn)
                    !g
                    g(1) = -Fxy*f_kapChe !YX :  Green_s_tr(Ioy,Isx)
                    g(2) = -Fyy*f_kapChe !YY : Green_s_tr(Ioy,Isy)
                    g(3) = -Fyz*f_kapChe !YZ : Green_s_tr(Ioy,Isz)
                else !Ioz
                    !f
                    Fxz = Gr_mn/r_mn * (rx*rz) *(-k_0**2. -3.*Tau_mn/r_mn)
                    Fyz = Gr_mn/r_mn * (ry*rz) *(-k_0**2. -3.*Tau_mn/r_mn)
                    Fzz = Gr_mn * (Tau_mn + k_0**2.*(r_mn - rz**2./r_mn) -3.*rz**2./r_mn**2. *Tau_mn)
                    !g
                    g(1) = -Fxz*f_kapChe !ZX : Green_s_tr(Ioz,Isx)
                    g(2) = -Fyz*f_kapChe !ZY : Green_s_tr(Ioz,Isy)
                    g(3) = -Fzz*f_kapChe !ZZ : Green_s_tr(Ioz,Isz)
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
        
    irow(3*sizeB+1) = v_irow;
    
    if (nnz .ne. 0) then ! nnz is ne 0 when SR_Green_s_tr_partial is called to calculate the CBFs, 
                         ! it is equal to 0 when the subroutine 
                         ! is called to determine fSR 
      if (curs_nnz .ne. nnz) then 
          Write(*,'(a)') 'Something went wrong when filling ZGreen_s_tr: curs_nnz .ne. nnz!'
          stop 1;
      endif
    endif
    
        
    End Subroutine SR_Green_s_tr_partial
    
SUBROUTINE SR_Green_s_tr_partial_FN(sizeB,CellsB,localfSR,nnz,FN_SR_Green_s_tr)

    USE Initialization
    USE common_variables
    Implicit none
    
    !IN/OUT
    Integer, INTENT(IN) :: sizeB
    real(kind=8), INTENT(IN) :: localfSR
    type (Cell), Dimension(sizeB),INTENT(IN) :: CellsB
    Integer, INTENT(OUT) :: nnz
    real(kind=8), INTENT(OUT) :: FN_SR_Green_s_tr
    
    ! Local    
    Integer :: curs_nnz,ii,jj,Is, Io, Isg, Iog, Iox, Ioy, Ioz, Isx, Isy, Isz
    Integer :: Ioo,num_dir
    Real (kind=8) :: rx, ry, rz, r_mn
    Real (kind=8) :: v_max,Green_s_tr_max,v_irow,threshold
    Real (kind=8),dimension(3) :: abs_g
    Complex	:: Gr_mn, Tau_mn, f_kapChe, Fxx, Fyy, Fzz, Fxy, Fyz, Fxz
    Real(kind=8) :: G_elts
    complex,dimension(3) :: g        
    
    nnz = 0;
    G_elts = 0;
    
    ! start by computing Green_s_tr(1,1), it will be used to elliminate 
    ! what we will consider as weak/non-significant interactions
    v_max = abs(1 - CellsB(1)%Znnpp)
    threshold = v_max/localfSR; !fct_SR;  ! the global fct_SR can be used if not all the blocks tested   
    
    If (homogs == 1) Then     
      Do Io = 1, sizeB   ! loop on row
        
        Do num_dir=1,3            
          ! Is = Io 
          Ioo = 3*(Io-1) + num_dir
                    
          !IoIo dir-dir (Io == Is)
          nnz = nnz + 1; 
          G_elts = G_elts + (abs(1 - CellsB(Io)%Znnpp))**2.;  ;     ! + (abs(Zpatch_e_spr(cc)))**2.           
                    
          ! Now Is > Io
          Do Is=Io+1, sizeB     !loop on col  ! if homogs=0 we will scan the entire matrix    
            Isx = 3*(Is-1)+1
            Isy = 3*(Is-1)+2
            Isz = 3*(Is-1)+3
            
            rx = CellsB(Io)%Xc - CellsB(Is)%Xc
            ry = CellsB(Io)%Yc - CellsB(Is)%Yc
            rz = CellsB(Io)%Zc - CellsB(Is)%Zc
    
            r_mn = sqrt(rx**2.+ry**2.+rz**2.)
    
            Gr_mn = exp(J*k_0*r_mn)/(4*Pi*r_mn**2.)
            Tau_mn = J*k_0 - 1/r_mn
            f_kapChe = CellsB(Is)%Kappa_n*CellsB(Is)%Che_n
            
            if (num_dir == 1) then  !Iox
              !f
              Fxx = Gr_mn * (Tau_mn + k_0**2.*(r_mn - rx**2./r_mn) -3.*rx**2./r_mn**2. *Tau_mn)
              Fxy = Gr_mn/r_mn * (rx*ry) *(-k_0**2.-3.*Tau_mn/r_mn)
              Fxz = Gr_mn/r_mn * (rx*rz) *(-k_0**2.-3.*Tau_mn/r_mn)
              !g
              g(1) = -Fxx*f_kapChe !XX :Green_s_tr(Iox,Isx)
              g(2) = -Fxy*f_kapChe !XY : Green_s_tr(Iox,Isy)
              g(3) = -Fxz*f_kapChe !XZ : Green_s_tr(Iox,Isz)
            elseif (num_dir == 2) then !Ioy
              !f
              Fxy = Gr_mn/r_mn * (rx*ry) *(-k_0**2. -3.*Tau_mn/r_mn)
              Fyy = Gr_mn * (Tau_mn + k_0**2.*(r_mn - ry**2./r_mn) -3.*ry**2./r_mn**2. *Tau_mn)
              Fyz = Gr_mn/r_mn * (ry*rz) *(-k_0**2. -3.*Tau_mn/r_mn)
              !g
              g(1) = -Fxy*f_kapChe !YX :  Green_s_tr(Ioy,Isx)
              g(2) = -Fyy*f_kapChe !YY : Green_s_tr(Ioy,Isy)
              g(3) = -Fyz*f_kapChe !YZ : Green_s_tr(Ioy,Isz)
            else !Ioz
              !f
              Fxz = Gr_mn/r_mn * (rx*rz) *(-k_0**2. -3.*Tau_mn/r_mn)
              Fyz = Gr_mn/r_mn * (ry*rz) *(-k_0**2. -3.*Tau_mn/r_mn)
              Fzz = Gr_mn * (Tau_mn + k_0**2.*(r_mn - rz**2./r_mn) -3.*rz**2./r_mn**2. *Tau_mn)
              !g
              g(1) = -Fxz*f_kapChe !ZX : Green_s_tr(Ioz,Isx)
              g(2) = -Fyz*f_kapChe !ZY : Green_s_tr(Ioz,Isy)
              g(3) = -Fzz*f_kapChe !ZZ : Green_s_tr(Ioz,Isz)
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
        Do Io = 1, sizeB   ! loop on row        
            Do num_dir=1,3  
          
                Ioo = 3*(Io-1) + num_dir
                
                ! here Is < Io
                Do Is=1, Io-1       
                Isx = 3*(Is-1)+1
                Isy = 3*(Is-1)+2
                Isz = 3*(Is-1)+3
            
                rx = CellsB(Io)%Xc - CellsB(Is)%Xc
                ry = CellsB(Io)%Yc - CellsB(Is)%Yc
                rz = CellsB(Io)%Zc - CellsB(Is)%Zc
    
                r_mn = sqrt(rx**2.+ry**2.+rz**2.)
    
                Gr_mn = exp(J*k_0*r_mn)/(4*Pi*r_mn**2.)
                Tau_mn = J*k_0 - 1/r_mn
                f_kapChe = CellsB(Is)%Kappa_n*CellsB(Is)%Che_n
            
                if (num_dir == 1) then  !Iox
                    !f
                    Fxx = Gr_mn * (Tau_mn + k_0**2.*(r_mn - rx**2./r_mn) -3.*rx**2./r_mn**2. *Tau_mn)
                    Fxy = Gr_mn/r_mn * (rx*ry) * (-k_0**2. -3.*Tau_mn/r_mn)
                    Fxz = Gr_mn/r_mn * (rx*rz) * (-k_0**2. -3.*Tau_mn/r_mn)
                    !g
                    g(1) = -Fxx*f_kapChe !XX :Green_s_tr(Iox,Isx)
                    g(2) = -Fxy*f_kapChe !XY : Green_s_tr(Iox,Isy)
                    g(3) = -Fxz*f_kapChe !XZ : Green_s_tr(Iox,Isz)
                elseif (num_dir == 2) then !Ioy
                    !f
                    Fxy = Gr_mn/r_mn * (rx*ry) * (-k_0**2. -3.*Tau_mn/r_mn)
                    Fyy = Gr_mn * (Tau_mn + k_0**2.*(r_mn - ry**2./r_mn) -3.*ry**2./r_mn**2. *Tau_mn)
                    Fyz = Gr_mn/r_mn * (ry*rz) * (-k_0**2. -3.*Tau_mn/r_mn)
                    !g
                    g(1) = -Fxy*f_kapChe !YX :  Green_s_tr(Ioy,Isx)
                    g(2) = -Fyy*f_kapChe !YY : Green_s_tr(Ioy,Isy)
                    g(3) = -Fyz*f_kapChe !YZ : Green_s_tr(Ioy,Isz)
                else !Ioz
                    !f
                    Fxz = Gr_mn/r_mn * (rx*rz) * (-k_0**2. -3.*Tau_mn/r_mn)
                    Fyz = Gr_mn/r_mn * (ry*rz) * (-k_0**2. -3.*Tau_mn/r_mn)
                    Fzz = Gr_mn * (Tau_mn + k_0**2.*(r_mn - rz**2./r_mn) -3.*rz**2./r_mn**2. *Tau_mn)
                    !g
                    g(1) = -Fxz*f_kapChe !ZX : Green_s_tr(Ioz,Isx)
                    g(2) = -Fyz*f_kapChe !ZY : Green_s_tr(Ioz,Isy)
                    g(3) = -Fzz*f_kapChe !ZZ : Green_s_tr(Ioz,Isz)
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
                G_elts = G_elts + abs((1 - CellsB(Io)%Znnpp))**2 ; 
          
                ! Now Is > Io
                Do Is=Io+1, sizeB     !loop on col  ! if homogs=0 we will scan the entire matrix    
                Isx = 3*(Is-1)+1
                Isy = 3*(Is-1)+2
                Isz = 3*(Is-1)+3
            
                rx = CellsB(Io)%Xc - CellsB(Is)%Xc
                ry = CellsB(Io)%Yc - CellsB(Is)%Yc
                rz = CellsB(Io)%Zc - CellsB(Is)%Zc
    
                r_mn = sqrt(rx**2.+ry**2.+rz**2.)
    
                Gr_mn = exp(J*k_0*r_mn)/(4*Pi*r_mn**2.)
                Tau_mn = J*k_0 - 1/r_mn
                f_kapChe = CellsB(Is)%Kappa_n*CellsB(Is)%Che_n
            
                if (num_dir == 1) then  !Iox
                    !f
                    Fxx = Gr_mn * (Tau_mn + k_0**2.*(r_mn - rx**2./r_mn) -3.*rx**2./r_mn**2. *Tau_mn)
                    Fxy = Gr_mn/r_mn * (rx*ry) *(-k_0**2. -3.*Tau_mn/r_mn)
                    Fxz = Gr_mn/r_mn * (rx*rz) *(-k_0**2. -3.*Tau_mn/r_mn)
                    !g
                    g(1) = -Fxx*f_kapChe !XX :Green_s_tr(Iox,Isx)
                    g(2) = -Fxy*f_kapChe !XY : Green_s_tr(Iox,Isy)
                    g(3) = -Fxz*f_kapChe !XZ : Green_s_tr(Iox,Isz)
                elseif (num_dir == 2) then !Ioy
                    !f
                    Fxy = Gr_mn/r_mn * (rx*ry) *(-k_0**2. -3.*Tau_mn/r_mn)
                    Fyy = Gr_mn * (Tau_mn + k_0**2.*(r_mn - ry**2./r_mn) -3.*ry**2./r_mn**2. *Tau_mn)
                    Fyz = Gr_mn/r_mn * (ry*rz) *(-k_0**2. -3.*Tau_mn/r_mn)
                    !g
                    g(1) = -Fxy*f_kapChe !YX :  Green_s_tr(Ioy,Isx)
                    g(2) = -Fyy*f_kapChe !YY : Green_s_tr(Ioy,Isy)
                    g(3) = -Fyz*f_kapChe !YZ : Green_s_tr(Ioy,Isz)
                else !Ioz
                    !f
                    Fxz = Gr_mn/r_mn * (rx*rz) *(-k_0**2. -3.*Tau_mn/r_mn)
                    Fyz = Gr_mn/r_mn * (ry*rz) *(-k_0**2. -3.*Tau_mn/r_mn)
                    Fzz = Gr_mn * (Tau_mn + k_0**2.*(r_mn - rz**2./r_mn) -3.*rz**2./r_mn**2. *Tau_mn)
                    !g
                    g(1) = -Fxz*f_kapChe !ZX : Green_s_tr(Ioz,Isx)
                    g(2) = -Fyz*f_kapChe !ZY : Green_s_tr(Ioz,Isy)
                    g(3) = -Fzz*f_kapChe !ZZ : Green_s_tr(Ioz,Isz)
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
    

SUBROUTINE DR_Green_s_tr_partial(sizeB,CellsB,klu_cel,klu,Green_s_tr)

    USE Initialization
    USE common_variables
    Implicit none
    
    !IN/OUT
    Integer, INTENT(IN) :: sizeB
    type (Cell), Dimension(sizeB), INTENT(IN) :: CellsB
    Integer, INTENT(IN):: klu_cel,klu
    COMPLEX(real64),Dimension(2*klu+1,3*sizeB),INTENT(OUT)::Green_s_tr

    ! Local
    Complex	:: Gr_mn, Tau_mn, f_kapChe, Fxx, Fyy, Fzz, Fxy, Fyz, Fxz
    Real (kind=8) :: rx, ry, rz, r_mn
    Integer :: Is, Io, Isg, Iog, Iox, Ioy, Ioz, Isx, Isy, Isz  
    Integer :: Index_col_Inf,Index_col_Sup,Iox_band,Ioy_band,Ioz_band
    
    Green_s_tr(:,:) = 0.D0;

    Do Io = 1, sizeB
    
        Index_col_Inf = max(1,Io-klu_cel)
        Index_col_Sup = min(Io+klu_cel,sizeB);
    
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
                Green_s_tr(Iox_band,Isx) = 1 - CellsB(Is)%Znnpp
                Ioy_band = klu+1+Ioy-Isy
                Green_s_tr(Ioy_band,Isy) = 1 - CellsB(Is)%Znnpp
                Ioz_band = klu+1+Ioz-Isz
                Green_s_tr(Ioz_band,Isz) = 1 - CellsB(Is)%Znnpp
            
                !! les autres (xy, yx, xz ...) restent a 0 pour ce cas ()
	        else
              rx = CellsB(Io)%Xc - CellsB(Is)%Xc
              ry = CellsB(Io)%Yc - CellsB(Is)%Yc
              rz = CellsB(Io)%Zc - CellsB(Is)%Zc

              r_mn = sqrt(rx**2.+ry**2.+rz**2.)
      
              Gr_mn = exp(J*k_0*r_mn)/(4*Pi*r_mn**2.)
              Tau_mn = J*k_0-1/r_mn
              f_kapChe = CellsB(Is)%Kappa_n*CellsB(Is)%Che_n

              ! if p .eq.  q
              Fxx = Gr_mn * (Tau_mn + k_0**2.*(r_mn - rx**2./r_mn)-3.*rx**2./r_mn**2. * Tau_mn)
              Fyy = Gr_mn * (Tau_mn + k_0**2.*(r_mn - ry**2./r_mn)-3.*ry**2./r_mn**2. * Tau_mn)
              Fzz = Gr_mn * (Tau_mn + k_0**2.*(r_mn - rz**2./r_mn)-3.*rz**2./r_mn**2. * Tau_mn)
    
              ! if p .ne. q
              Fxy = Gr_mn/r_mn * (rx*ry) * (-k_0**2. -3.*Tau_mn/r_mn)
              Fyz = Gr_mn/r_mn * (ry*rz) * (-k_0**2. -3.*Tau_mn/r_mn)
              Fxz = Gr_mn/r_mn * (rx*rz) * (-k_0**2. -3.*Tau_mn/r_mn)      
      
              Iox_band = klu+1+Iox-Isx
              Green_s_tr(Iox_band,Isx)= -Fxx*f_kapChe !XX
              Iox_band = klu+1+Iox-Isy
              Green_s_tr(Iox_band,Isy)= -Fxy*f_kapChe !XY
              Iox_band = klu+1+Iox-Isz
              Green_s_tr(Iox_band,Isz)= -Fxz*f_kapChe !XZ
          
              Ioy_band = klu+1+Ioy-Isx
              Green_s_tr(Ioy_band,Isx)= -Fxy*f_kapChe !YX
              Ioy_band = klu+1+Ioy-Isy
              Green_s_tr(Ioy_band,Isy)= -Fyy*f_kapChe !YY
              Ioy_band = klu+1+Ioy-Isz
              Green_s_tr(Ioy_band,Isz)= -Fyz*f_kapChe !YZ
          
              Ioz_band = klu+1+Ioz-Isx
              Green_s_tr(Ioz_band,Isx)= -Fxz*f_kapChe !ZX
              Ioz_band = klu+1+Ioz-Isy
              Green_s_tr(Ioz_band,Isy)= -Fyz*f_kapChe !ZY
              Ioz_band = klu+1+Ioz-Isz
              Green_s_tr(Ioz_band,Isz)= -Fzz*f_kapChe !ZZ  
            endIf    
          endDo
    endDo    

End Subroutine DR_Green_s_tr_partial



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
    COMPLEX :: Gr_mn, Tau_mn, f_kapChe
    COMPLEX :: Fxx, Fyy, Fzz, Fxy, Fyz,Fxz,Fyx,Fzy,Fzx,Gx,Gy,Gz
    Complex :: K11x, K11z, K22z, K12z, K21z
    Real(kind=8) :: RE, RM, Ang_inc, P1, P2, norme, Px, Py
    REAL(kind=8) :: rx, ry, rz, r_mn, rxy
    
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
            
            rx = Cells(Io)%Xc - Cells(Is)%Xc
            ry = Cells(Io)%Yc - Cells(Is)%Yc
            rz = Cells(Io)%Zc - Cells(Is)%Zc
            r_mn = sqrt(rx**2.+ry**2.+rz**2.)
      
            Gr_mn = exp(J*k_0*r_mn)/(4*Pi*r_mn**2.)
            Tau_mn = J*k_0 - 1/r_mn
            f_kapChe = Cells(Is)%Kappa_n*Cells(Is)%Che_n

            Fxx = Gr_mn * (Tau_mn + k_0**2.*(r_mn - rx**2./r_mn) -3.*rx**2./r_mn**2. *Tau_mn)
            Fxy = Gr_mn/r_mn * (rx*ry) *(-k_0**2. -3.*Tau_mn/r_mn)
            Fxz = Gr_mn/r_mn * (rx*rz) *(-k_0**2. -3.*Tau_mn/r_mn)
      
            Matrix_Green_Col(curs_lig+1,1)= - Fxx*f_kapChe !!XX
            Matrix_Green_Col(curs_lig+2,1)= - Fxy*f_kapChe !!YX
            Matrix_Green_Col(curs_lig+3,1)= - Fxz*f_kapChe !!ZX

            Curs_lig = Curs_lig + 3 
            Io = Io + 1
        Enddo
    ElseIf (Isc==2) Then !! .Y
        Do ii=1,nb_cels_i       
            
            rx = Cells(Io)%Xc - Cells(Is)%Xc
            ry = Cells(Io)%Yc - Cells(Is)%Yc
            rz = Cells(Io)%Zc - Cells(Is)%Zc

            r_mn = sqrt(rx**2.+ry**2.+rz**2.)

            Gr_mn = exp(J*k_0*r_mn)/(4*Pi*r_mn**2.)
            Tau_mn = J*k_0-1/r_mn
            f_kapChe = Cells(Is)%Kappa_n*Cells(Is)%Che_n

            Fxy = Gr_mn/r_mn * (rx*ry) *(-k_0**2. -3.*Tau_mn/r_mn)
            Fyy = Gr_mn * (Tau_mn + k_0**2.*(r_mn - ry**2./r_mn) -3.*ry**2./r_mn**2. *Tau_mn)
            Fyz = Gr_mn/r_mn * (ry*rz) *(-k_0**2. -3.*Tau_mn/r_mn)
            
          
            Matrix_Green_Col(curs_lig+1,1)= - Fxy*f_kapChe !!XY
            Matrix_Green_Col(curs_lig+2,1)= - Fyy*f_kapChe !!YY
            Matrix_Green_Col(curs_lig+3,1)= - Fyz*f_kapChe !!ZY           
          
            Curs_lig = Curs_lig + 3 
            Io = Io + 1
        Enddo
    Else  !! .Z
        Do ii=1,nb_cels_i      
            
            rx = Cells(Io)%Xc - Cells(Is)%Xc
            ry = Cells(Io)%Yc - Cells(Is)%Yc
            rz = Cells(Io)%Zc - Cells(Is)%Zc

            r_mn = sqrt(rx**2.+ry**2.+rz**2.)

            Gr_mn = exp(J*k_0*r_mn)/(4*Pi*r_mn**2.)
            Tau_mn = J*k_0 - 1/r_mn
            f_kapChe = Cells(Is)%Kappa_n*Cells(Is)%Che_n
            
            Fxz = Gr_mn/r_mn * (rx*rz) *(-k_0**2. -3.*Tau_mn/r_mn)
            Fyz = Gr_mn/r_mn * (ry*rz) *(-k_0**2. -3.*Tau_mn/r_mn)
            Fzz = Gr_mn * (Tau_mn + k_0**2.*(r_mn - rz**2./r_mn) -3.*rz**2./r_mn**2. *Tau_mn)
            
            Matrix_Green_Col(curs_lig+1,1)= -Fxz *f_kapChe !!XZ
            Matrix_Green_Col(curs_lig+2,1)= -Fyz*f_kapChe  !!YZ
            Matrix_Green_Col(curs_lig+3,1)= -Fzz*f_kapChe  !!ZZ           
          
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
    COMPLEX :: Gr_mn, Tau_mn, f_kapChe
    COMPLEX :: Fxx, Fyy, Fzz, Fxy, Fyz,Fxz,Fyx,Fzy,Fzx,Gx,Gy,Gz
    Complex :: K11x, K11z, K22z, K12z, K21z
    Real(kind=8) :: RE, RM, Ang_inc, P1, P2, norme, Px, Py
    REAL(kind=8) :: rx, ry, rz, r_mn, rxy
    
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
            
            rx = CellsB1(Io)%Xc - CellsB2(Is)%Xc
            ry = CellsB1(Io)%Yc - CellsB2(Is)%Yc
            rz = CellsB1(Io)%Zc - CellsB2(Is)%Zc
            
            r_mn = sqrt(rx**2.+ry**2.+rz**2.)
      
            Gr_mn = exp(J*k_0*r_mn)/(4*Pi*r_mn**2.)
            Tau_mn = J*k_0 - 1/r_mn
            f_kapChe = CellsB2(Is)%Kappa_n*CellsB2(Is)%Che_n

            Fxx = Gr_mn * (Tau_mn + k_0**2.*(r_mn - rx**2./r_mn) -3.*rx**2./r_mn**2. *Tau_mn)
            Fxy = Gr_mn/r_mn * (rx*ry) *(-k_0**2.-3.*Tau_mn/r_mn)
            Fxz = Gr_mn/r_mn * (rx*rz) *(-k_0**2.-3.*Tau_mn/r_mn)
      
            Matrix_Green_Col(Index_lig+1,1)= - Fxx*f_kapChe !!XX
            Matrix_Green_Col(Index_lig+2,1)= - Fxy*f_kapChe !!YX
            Matrix_Green_Col(Index_lig+3,1)= - Fxz*f_kapChe !!ZX

            Index_lig = Index_lig + 3
        Enddo
    ElseIf (Isc==2) Then !! .Y
        Do Io=1,nb1       
            
            rx = CellsB1(Io)%Xc - CellsB2(Is)%Xc
            ry = CellsB1(Io)%Yc - CellsB2(Is)%Yc
            rz = CellsB1(Io)%Zc - CellsB2(Is)%Zc

            r_mn = sqrt(rx**2.+ry**2.+rz**2.)

            Gr_mn = exp(J*k_0*r_mn)/(4*Pi*r_mn**2.)
            Tau_mn = J*k_0 - 1/r_mn
            f_kapChe = CellsB2(Is)%Kappa_n*CellsB2(Is)%Che_n

            Fxy = Gr_mn/r_mn * (rx*ry) *(-k_0**2.-3.*Tau_mn/r_mn)
            Fyy = Gr_mn * (Tau_mn + k_0**2.*(r_mn - ry**2./r_mn) -3.*ry**2./r_mn**2. *Tau_mn)
            Fyz = Gr_mn/r_mn * (ry*rz) *(-k_0**2.-3.*Tau_mn/r_mn)            
          
            Matrix_Green_Col(Index_lig+1,1)= - Fxy*f_kapChe !!XY
            Matrix_Green_Col(Index_lig+2,1)= - Fyy*f_kapChe !!YY
            Matrix_Green_Col(Index_lig+3,1)= - Fyz*f_kapChe !!ZY           
          
            Index_lig = Index_lig + 3 
        Enddo
    Else  !! .Z
        Do Io=1,nb1      
            
            rx = CellsB1(Io)%Xc - CellsB2(Is)%Xc
            ry = CellsB1(Io)%Yc - CellsB2(Is)%Yc
            rz = CellsB1(Io)%Zc - CellsB2(Is)%Zc

            r_mn = sqrt(rx**2.+ry**2.+rz**2.)

            Gr_mn = exp(J*k_0*r_mn)/(4*Pi*r_mn**2.)
            Tau_mn = J*k_0 - 1/r_mn
            f_kapChe = CellsB2(Is)%Kappa_n*CellsB2(Is)%Che_n
            
            Fxz = Gr_mn/r_mn * (rx*rz) *(-k_0**2. -3.*Tau_mn/r_mn)
            Fyz = Gr_mn/r_mn * (ry*rz) *(-k_0**2. -3.*Tau_mn/r_mn)
            Fzz = Gr_mn * (Tau_mn + k_0**2.*(r_mn - rz**2./r_mn) -3.*rz**2./r_mn**2. *Tau_mn)
            
            Matrix_Green_Col(Index_lig+1,1)= -Fxz *f_kapChe !!XZ
            Matrix_Green_Col(Index_lig+2,1)= -Fyz*f_kapChe  !!YZ
            Matrix_Green_Col(Index_lig+3,1)= -Fzz*f_kapChe  !!ZZ           
          
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
    COMPLEX :: Gr_mn, Tau_mn, f_kapChe
    COMPLEX :: Fxx, Fyy, Fzz, Fxy, Fyz,Fxz,Fyx,Fzy,Fzx,Gx,Gy,Gz
    Complex :: K11x, K11z, K22z, K12z, K21z
    Real(kind=8) :: RE, RM, Ang_inc, P1, P2, norme, Px, Py
    Real(kind=8) :: rx, ry, rz, r_mn, rxy
  
    
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
            
            rx = Cells(Io)%Xc - Cells(Is)%Xc
            ry = Cells(Io)%Yc - Cells(Is)%Yc
            rz = Cells(Io)%Zc - Cells(Is)%Zc

            r_mn = sqrt(rx**2.+ry**2.+rz**2.)

            Gr_mn = exp(J*k_0*r_mn)/(4*Pi*r_mn**2.)
            Tau_mn = J*k_0 - 1/r_mn
            f_kapChe = Cells(Is)%Kappa_n*Cells(Is)%Che_n

            Fxx = Gr_mn * (Tau_mn + k_0**2.*(r_mn - rx**2./r_mn) -3.*rx**2./r_mn**2. *Tau_mn) 
            Fxy = Gr_mn/r_mn * (rx*ry) *(-k_0**2. -3.*Tau_mn/r_mn)
            Fxz = Gr_mn/r_mn * (rx*rz) *(-k_0**2. -3.*Tau_mn/r_mn)

            Matrix_Green_row(1,Index_Col+1)= - Fxx*f_kapChe !XX
            Matrix_Green_row(1,Index_Col+2)= - Fxy*f_kapChe !XY
            Matrix_Green_row(1,Index_Col+3)= - Fxz*f_kapChe !XZ             
                    
            Index_Col = Index_Col +3
            Is = Is + 1
        Enddo
    ElseIf (Ioc == 2) Then  !!Y.
        Do jj=1,nb_cels_j      
            
            rx = Cells(Io)%Xc - Cells(Is)%Xc
            ry = Cells(Io)%Yc - Cells(Is)%Yc
            rz = Cells(Io)%Zc - Cells(Is)%Zc

            r_mn = sqrt(rx**2.+ry**2.+rz**2.)

            Gr_mn = exp(J*k_0*r_mn)/(4*Pi*r_mn**2.)
            Tau_mn = J*k_0 - 1/r_mn
            f_kapChe = Cells(Is)%Kappa_n*Cells(Is)%Che_n

            Fxy = Gr_mn/r_mn * (rx*ry) *(-k_0**2. -3.*Tau_mn/r_mn)
            Fyy = Gr_mn * (Tau_mn + k_0**2.*(r_mn - ry**2./r_mn) -3.*ry**2./r_mn**2. *Tau_mn)
            Fyz = Gr_mn/r_mn * (ry*rz) *(-k_0**2. -3.*Tau_mn/r_mn)
            
            Matrix_Green_row(1,Index_Col+1)= - Fxy*f_kapChe   !YX
            Matrix_Green_row(1,Index_Col+2)= - Fyy*f_kapChe   !YY
            Matrix_Green_row(1,Index_Col+3)= - Fyz*f_kapChe   !YZ
          
            Index_Col = Index_Col +3
            Is = Is + 1
        Enddo
    Else !!Z.
        Do jj=1,nb_cels_j      
            
            rx = Cells(Io)%Xc - Cells(Is)%Xc
            ry = Cells(Io)%Yc - Cells(Is)%Yc
            rz = Cells(Io)%Zc - Cells(Is)%Zc

            r_mn = sqrt(rx**2.+ry**2.+rz**2.)

            Gr_mn = exp(J*k_0*r_mn)/(4*Pi*r_mn**2.)
            Tau_mn = J*k_0-1/r_mn
            f_kapChe = Cells(Is)%Kappa_n*Cells(Is)%Che_n
            
            Fxz = Gr_mn/r_mn * (rx*rz) *(-k_0**2.-3.*Tau_mn/r_mn)
            Fyz = Gr_mn/r_mn * (ry*rz) *(-k_0**2.-3.*Tau_mn/r_mn)
            Fzz = Gr_mn * (Tau_mn + k_0**2.*(r_mn - rz**2./r_mn) -3.*rz**2./r_mn**2. *Tau_mn)

            Matrix_Green_row(1,Index_Col+1)= - Fxz*f_kapChe !ZX
            Matrix_Green_row(1,Index_Col+2)= - Fyz*f_kapChe !ZY
            Matrix_Green_row(1,Index_Col+3)= - Fzz*f_kapChe !ZZ
          
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
    COMPLEX :: Gr_mn, Tau_mn, f_kapChe
    COMPLEX :: Fxx, Fyy, Fzz, Fxy, Fyz,Fxz,Fyx,Fzy,Fzx,Gx,Gy,Gz
    Complex :: K11x, K11z, K22z, K12z, K21z
    Real(kind=8) :: RE, RM, Ang_inc, P1, P2, norme, Px, Py
    Real(kind=8) :: rx, ry, rz, r_mn, rxy
  
    
    Ioc= mod(irow,3)
    If (Ioc .NE. 0) Then
        Io = irow/3 + 1
    Else
        Io = irow/3
    Endif
    Index_Col = 0  
    
    If (Ioc ==1) Then !! X.
        Do Is=1,nb2      
            
            rx = CellsB1(Io)%Xc - CellsB2(Is)%Xc
            ry = CellsB1(Io)%Yc - CellsB2(Is)%Yc
            rz = CellsB1(Io)%Zc - CellsB2(Is)%Zc

            r_mn = sqrt(rx**2.+ry**2.+rz**2.)

            Gr_mn = exp(J*k_0*r_mn)/(4*Pi*r_mn**2.)
            Tau_mn = J*k_0-1/r_mn
            f_kapChe = CellsB2(Is)%Kappa_n*CellsB2(Is)%Che_n

            Fxx = Gr_mn * (Tau_mn + k_0**2.*(r_mn - rx**2./r_mn) -3.*rx**2./r_mn**2. *Tau_mn)
            Fxy = Gr_mn/r_mn * (rx*ry) *(-k_0**2.-3.*Tau_mn/r_mn)
            Fxz = Gr_mn/r_mn * (rx*rz) *(-k_0**2.-3.*Tau_mn/r_mn)

            Matrix_Green_row(1,Index_Col+1)= - Fxx*f_kapChe !XX
            Matrix_Green_row(1,Index_Col+2)= - Fxy*f_kapChe !XY
            Matrix_Green_row(1,Index_Col+3)= - Fxz*f_kapChe !XZ             
                    
            Index_Col = Index_Col +3
        Enddo
    ElseIf (Ioc == 2) Then  !!Y.
        Do Is=1,nb2      
            
            rx = CellsB1(Io)%Xc - CellsB2(Is)%Xc
            ry = CellsB1(Io)%Yc - CellsB2(Is)%Yc
            rz = CellsB1(Io)%Zc - CellsB2(Is)%Zc

            r_mn = sqrt(rx**2.+ry**2.+rz**2.)

            Gr_mn = exp(J*k_0*r_mn)/(4*Pi*r_mn**2.)
            Tau_mn = J*k_0 - 1/r_mn
            f_kapChe = CellsB2(Is)%Kappa_n*CellsB2(Is)%Che_n

            Fxy = Gr_mn/r_mn * (rx*ry) *(-k_0**2.-3.*Tau_mn/r_mn)
            Fyy = Gr_mn * (Tau_mn + k_0**2.*(r_mn - ry**2./r_mn) -3.*ry**2./r_mn**2. *Tau_mn)
            Fyz = Gr_mn/r_mn * (ry*rz) *(-k_0**2.-3.*Tau_mn/r_mn)
            
            Matrix_Green_row(1,Index_Col+1)= - Fxy*f_kapChe   !YX
            Matrix_Green_row(1,Index_Col+2)= - Fyy*f_kapChe   !YY
            Matrix_Green_row(1,Index_Col+3)= - Fyz*f_kapChe   !YZ
          
            Index_Col = Index_Col +3
        Enddo
    Else !!Z.
        Do Is=1,nb2      
            
            rx = CellsB1(Io)%Xc - CellsB2(Is)%Xc
            ry = CellsB1(Io)%Yc - CellsB2(Is)%Yc
            rz = CellsB1(Io)%Zc - CellsB2(Is)%Zc

            r_mn = sqrt(rx**2.+ry**2.+rz**2.)

            Gr_mn = exp(J*k_0*r_mn)/(4*Pi*r_mn**2.)
            Tau_mn = J*k_0-1/r_mn
            f_kapChe = CellsB2(Is)%Kappa_n*CellsB2(Is)%Che_n
            
            Fxz = Gr_mn/r_mn * (rx*rz) *(-k_0**2. -3.*Tau_mn/r_mn)
            Fyz = Gr_mn/r_mn * (ry*rz) *(-k_0**2. -3.*Tau_mn/r_mn)
            Fzz = Gr_mn * (Tau_mn + k_0**2.*(r_mn - rz**2./r_mn) -3.*rz**2./r_mn**2. *Tau_mn)

            Matrix_Green_row(1,Index_Col+1)= - Fxz*f_kapChe !ZX
            Matrix_Green_row(1,Index_Col+2)= - Fyz*f_kapChe !ZY
            Matrix_Green_row(1,Index_Col+3)= - Fzz*f_kapChe !ZZ
          
            Index_Col = Index_Col +3
        Enddo
    EndIf  
 END SUBROUTINE computeBlockRow_SMW
 

SUBROUTINE Green_s_dt(Nc,Cells_in,theta_capteur,phi_capteur,Green_dt) 

    USE Initialization
    USE common_variables
    IMPLICIT NONE
    
    Integer, INTENT(IN) :: Nc
    type (Cell), Dimension(Nc), INTENT(IN) :: Cells_in
    Real(kind=8), INTENT(IN) :: theta_capteur,phi_capteur
    COMPLEX(real64), Dimension(3,3*Nc), INTENT(OUT) :: Green_dt
    
    !Local
    Integer Is, Isx, Isy, Isz
    Real(kind=8) :: xc,yc,zc,x_cap,y_cap,z_cap, rx, ry, rz, r_mn
    COMPLEX(real64) :: Gr_mn, Tau_mn, f_kapChe, Fxx, Fyy, Fzz, Fxy, Fyz, Fxz    
    
    DO Is=1, Nc
        
        Isx=3*(Is-1)+1
        Isy=3*(Is-1)+2
        Isz=3*(Is-1)+3
        
        ! Cell in scatterer x, y & z
        xc = Cells_in(Is)%Xc
        yc = Cells_in(Is)%Yc
        zc = Cells_in(Is)%Zc
        
        ! Receiver x, y & z
        x_cap = Rso*cos(theta_capteur*Pi/180.) 
        y_cap = Rso*sin(theta_capteur*Pi/180.)*cos(phi_capteur*Pi/180.) 
        z_cap = Rso*sin(theta_capteur*Pi/180.)*sin(phi_capteur*Pi/180.)            
        
        rx = x_cap - xc
        ry = y_cap - yc
        rz = z_cap - zc
        r_mn = sqrt(rx**2.+ry**2.+rz**2.)   ! r_mn
    
        Gr_mn = exp(J*k_0*r_mn)/(4*Pi*r_mn**2.)      ! Green_mn
        Tau_mn = J*k_0 - 1/r_mn                      ! Tau_mn 
        f_kapChe = Cells_in(Is)%Kappa_n*Cells_in(Is)%Che_n ! factor :  Kappa_n * Che_n   ! Check from equations !
    
        Fxx = Gr_mn * (Tau_mn + k_0**2. *(r_mn -rx**2./ r_mn)-3.*rx**2./r_mn**2. * Tau_mn)
        Fyy = Gr_mn * (Tau_mn + k_0**2. *(r_mn-ry**2. / r_mn)-3.*ry**2./r_mn**2. * Tau_mn)
        Fzz = Gr_mn * (Tau_mn + k_0**2. *(r_mn-rz**2. / r_mn)-3.*rz**2./r_mn**2. * Tau_mn)
                
        Fxy = Gr_mn/r_mn * (rx * ry) * (-k_0**2. -3.*Tau_mn/r_mn) 
        Fxz = Gr_mn/r_mn * (rx * rz) * (-k_0**2. -3.*Tau_mn/r_mn)
        Fyz = Gr_mn/r_mn * (ry * rz) * (-k_0**2. -3.*Tau_mn/r_mn) 
    
        Green_dt(1,Isx)=Fxx*f_kapChe !!XX
        Green_dt(1,Isy)=Fxy*f_kapChe !!XY
        Green_dt(1,Isz)=Fxz*f_kapChe !!XZ
    
        Green_dt(2,Isx)=Fxy*f_kapChe !!YX
        Green_dt(2,Isy)=Fyy*f_kapChe !!YY
        Green_dt(2,Isz)=Fyz*f_kapChe !!YZ
    
        Green_dt(3,Isx)=Fxz*f_kapChe !!ZX
        Green_dt(3,Isy)=Fyz*f_kapChe !!ZY
        Green_dt(3,Isz)=Fzz*f_kapChe !!ZZ
    ENDDO

END SUBROUTINE Green_s_dt