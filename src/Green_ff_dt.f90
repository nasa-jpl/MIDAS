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