SUBROUTINE SetCellsParams(SimScatterer,Cells,UCells)
    
    ! Modifs 8/29/2019 : Here I started implementing adaptive mesh depending on the dielectric properties of each cell ! if Adapt_mesh == 1, each cell 
    ! which is not respecting the validity crieteari is divided to smaller cells, then the cells and the blocks parameters are updated accordingly
    ! if Adapt_mesh == 0, we only update cells properties (Lambda_s and D_lambda) depending on the simlated wavelength.
    
    USE Initialization
    USE common_variables
    USE iso_fortran_env
    Implicit NONE
    
    ! IN/OUT
    type (Scatterer), INTENT(INOUT) :: SimScatterer
    type (Cell), Dimension(Nbc), INTENT(INOUT):: Cells
    type (Cell), Dimension(:), allocatable, INTENT(OUT):: UCells ! Updated Cells 

    ! Local
    Integer :: nn, bb,cn, new_Nbc,fmr,Dlambda_min,ix,iy,iz
    Integer :: num_B, found,pos_cellii_here,N
    real(kind=8) :: Sc,a_n,Sc_o,Sc_n,x_o,y_o,z_o,x_n,y_n,z_n
    COMPLEX(real64) :: Eps_c,Che_n
    
    ! let us start with a simple constant fmr at 2 (every cell not respecting the validity criteria is divided into 8 cells)
    fmr = 2;
    Dlambda_min = 10;
    
    ! This subroutine is specific to the Multi frequency code. It enables to update the parameters 
    ! parameter_Ce, parameter_Sing and parameter_Const  of the cells because they depend on k_air and Eps_p        
    
    if (Adapt_mesh .eq. 0) then 
        Allocate(UCells(Nbc)); 
        UCells = Cells;  
    Else
        new_Nbc = 0;
        Do nn=1,Nbc
            !if (Cells(ii)%Dlamb_cell .lt. Dlambda_min) then 
                new_Nbc = new_Nbc + fmr**3;
            !else
            !    new_Nbc = new_Nbc + 1;                
            !endif
        EndDo
        Allocate(UCells(new_Nbc));
        cn = 1;
        Do nn=1,Nbc
            !if (Cells(nn)%Dlamb_n .lt. Dlambda_min) then 
                x_o = Cells(nn)%Xc;
                y_o = Cells(nn)%Yc;
                z_o = Cells(nn)%Zc;
                Sc_o = Cells(nn)%Sc;
                Sc_n = Sc_o/fmr;
                Do ix=1,fmr
                    Do iy = 1,fmr
                        Do iz = 1, fmr
                            x_n = x_o-0.5*Sc_o + (ix-0.5)*Sc_n;
                            y_n = y_o-0.5*Sc_o + (iy-0.5)*Sc_n;
                            z_n = z_o-0.5*Sc_o + (iz-0.5)*Sc_n;
                        
                            UCells(cn)%n_cell = cn
                            UCells(cn)%Xc = x_n
                            UCells(cn)%Yc = y_n
                            UCells(cn)%Zc = z_n
                            UCells(cn)%Sc = Sc_n
                            
                            UCells(cn)%m_n = Cells(nn)%m_n;
                            UCells(cn)%Eps_n = Cells(nn)%Eps_n;
                            UCells(cn)%lambda_n = Cells(nn)%lambda_n ;
                            ! Update cell Dlam
                            UCells(cn)%Dlamb_n = Cells(nn)%lambda_n/Sc_n;
                            
                            UCells(cn)%n_block = Cells(nn)%n_block
                            UCells(cn)%n_diel = Cells(nn)%n_diel
                            
                            cn = cn + 1;
                        EndDo
                    EndDo
                EndDo
            !else
            !    Upd_Cells(curs_new) = Cells(ii);
            !    curs_new = curs_new + 1;            
            !endif
        EndDo
        Nbc = cn - 1;
    EndIf
    
    
    Do nn = 1, Nbc                    
        Sc = UCells(nn)%Sc 
        Eps_c = UCells(nn)%Eps_n
        
        ! update parameters with regard to the current k0=2pi/lambda and Eps_cell (depends also on the wavelength)
        ! radius of equivalent sphere a_n, and dielectric constrast Che for each cell n
        a_n = Sc*(0.75/Pi)**(1./3.)
        UCells(nn)%a_n = a_n
        
        Che_n = Eps_c - 1.
        UCells(nn)%Che_n = Che_n
        
        ! Zmn,pq with m=n & p=q
        UCells(nn)%Znnpp = ((2./3.)*exp(J*k_0*a_n) * (1. - J*k_0*a_n) - 1) * Che_n
        
        ! kappa = coefficient for the approximation of the integration of the green's function term over
        ! a cubical cell of size Sc by an integral over a sphere of equivalent radius a_n       
        UCells(nn)%Kappa_n = 4*Pi*(sin(k_0*a_n) - k_0*a_n*cos(k_0*a_n) )/k_0**3.  
    EndDo      

END SUBROUTINE SetCellsParams