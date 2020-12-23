SUBROUTINE UpdateCellsParameters(SimParticle,Cells,Upd_Cells)!,CBFM_Blocks,CBFM_Blocks_Ext,Upd_CBFM_Blocks_Ext)
    
    ! Modifs 8/29/2019 : Here I started implementing adaptive mesh depending on the dielectric properties of each cell ! if Adapt_mesh == 1, each cell 
    ! which is not respecting the validity crieteari is divided to smaller cells, then the cells and the blocks parameters are updated accordingly
    ! if Adapt_mesh == 0, we only update cells properties (Lambda_s and D_lambda) depending on the simlated wavelength.
    
    USE Initialization
    USE common_variables
    USE iso_fortran_env
    Implicit NONE
    
    ! IN/OUT
    type (Particle), INTENT(INOUT) :: SimParticle
    type (Cell), Dimension(Nbc), INTENT(INOUT):: Cells
    type (Cell), Dimension(:), allocatable, INTENT(OUT):: Upd_Cells
    !type (CBFM_Block), Dimension(NBlocks), INTENT(INOUT):: CBFM_Blocks
    !Integer, Dimension(NBlocks,Nbc_ext), INTENT(IN):: CBFM_Blocks_Ext
    !Integer, Dimension(:,:),allocatable, INTENT(OUT):: Upd_CBFM_Blocks_Ext

    ! Local
    Integer :: ii, bb,curs_new, new_Nbc,fmr,Dlambda_min,ix,iy,iz
    Integer :: num_B, found,pos_cellii_here,N
    real(kind=8) :: Sc,A,B,C,Sc_o,Sc_n,x_o,y_o,z_o,x_n,y_n,z_n
    COMPLEX(real64) :: Eps_c,Ac,Ai, Bi, Ci
    
    ! let us start with a simple constant fmr at 2 (every cell not respecting the validity criteria is divided into 8 cells)
    fmr = 2;
    Dlambda_min = 10;
    
    ! This subroutine is specific to the Multi frequency code. It enables to update the parameters 
    ! parameter_Ce, parameter_Sing and parameter_Const  of the cells because they depend on k_air and Eps_p        
    
    if (Adapt_mesh .eq. 0) then 
        Allocate(Upd_Cells(Nbc)); 
        Upd_Cells = Cells;  
    Else
        new_Nbc = 0;
        Do ii=1,Nbc
            !if (Cells(ii)%Dlamb_cell .lt. Dlambda_min) then 
                new_Nbc = new_Nbc + fmr**3;
            !else
            !    new_Nbc = new_Nbc + 1;                
            !endif
        EndDo
        Allocate(Upd_Cells(new_Nbc));
        curs_new = 1;
        Do ii=1,Nbc
            !if (Cells(ii)%Dlamb_cell .lt. Dlambda_min) then 
                x_o = Cells(ii)%Xc;
                y_o = Cells(ii)%Yc;
                z_o = Cells(ii)%Zc;
                Sc_o = Cells(ii)%Sc;
                Sc_n = Sc_o/fmr;
                Do ix=1,fmr
                    Do iy = 1,fmr
                        Do iz = 1, fmr
                            x_n = x_o-0.5*Sc_o + (ix-0.5)*Sc_n;
                            y_n = y_o-0.5*Sc_o + (iy-0.5)*Sc_n;
                            z_n = z_o-0.5*Sc_o + (iz-0.5)*Sc_n;
                        
                            Upd_Cells(curs_new)%num_cell = curs_new
                            Upd_Cells(curs_new)%Xc = x_n
                            Upd_Cells(curs_new)%Yc = y_n
                            Upd_Cells(curs_new)%Zc = z_n
                            Upd_Cells(curs_new)%Sc = Sc_n
                            
                            Upd_Cells(curs_new)%m_cell = Cells(ii)%m_cell;
                            Upd_Cells(curs_new)%Eps_cell = Cells(ii)%Eps_cell;
                            Upd_Cells(curs_new)%lambda_cell = Cells(ii)%lambda_cell ;
                            ! Update cell Dlam
                            Upd_Cells(curs_new)%Dlamb_cell = Cells(ii)%lambda_cell/Sc_n;
                            
                            Upd_Cells(curs_new)%num_block = Cells(ii)%num_block
                            Upd_Cells(curs_new)%num_diel = Cells(ii)%num_diel
                            
                            curs_new = curs_new + 1;
                        EndDo
                    EndDo
                EndDo
            !else
            !    Upd_Cells(curs_new) = Cells(ii);
            !    curs_new = curs_new + 1;            
            !endif
        EndDo
        Nbc = curs_new - 1;
        SimParticle%Nbc_p = Nbc;
        !Nbc_ext = maxval(CBFM_Blocks%Nbc_ext);
    EndIf
    
    
    Do ii =1,Nbc                    
        Sc = Upd_Cells(ii)%Sc 
        Eps_c = Upd_Cells(ii)%Eps_cell
        
        ! update parameter_Ce, parameter_Sing and parameter_Const 
        ! with regard to the current k0=2pi/lambda and Eps_cell (depends also on the wavelength)
        ! parameter_Ce
        Upd_Cells(ii)%parameter_Rad = Sc*(0.75/Pi)**(1./3.)
        Upd_Cells(ii)%parameter_Ce = Eps_c - 1.
        ! parameter_Sing
        Ac = J*K_air*Upd_Cells(ii)%parameter_Rad
        Ai = (2./3.)*exp(Ac)
        Bi = 1. -J*K_air*Upd_Cells(ii)%parameter_Rad
        Ci = Upd_Cells(ii)%parameter_Ce
        Upd_Cells(ii)%parameter_Sing = (Ai*Bi-1)*Ci
        ! parameter_Const 
        A = sin(K_air*Upd_Cells(ii)%parameter_Rad)
        B = K_air*Upd_Cells(ii)%parameter_Rad*cos(K_air*Upd_Cells(ii)%parameter_Rad)
        C = K_air**3.        
        Upd_Cells(ii)%parameter_Const = (A-B)/C   
    EndDo
    
        

END SUBROUTINE UpdateCellsParameters