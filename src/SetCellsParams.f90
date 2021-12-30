SUBROUTINE SetCellsParams(SimScatterer,Cells,UCells)
    
    ! Modifs 12/6/2021 : Here we update cells properties (Lambda_s and D_lambda) depending on the simlated wavelength.
    
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
    
    Allocate(UCells(Nbc)); 
    UCells = Cells;
    
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