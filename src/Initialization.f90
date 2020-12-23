MODULE Initialization
    
    USE iso_fortran_env
    Implicit none 
    
    !! Definition of the object dipole for the transmitter/receivers  
    type Dipole 
	    Real(kind=8) :: theta,phi !,x,y,z
    end type Dipole 

    !! Definition of the object particule
    type Particle
        Integer :: type_p
        Character*9 :: info_p
        Real(kind=8) :: Dp,ap
        COMPLEX(real64) :: m_p_min,Eps_p_min,m_p_max, Eps_p_max
        Real(kind=8) :: lambda_p_min,lambda_p_max, Sc_p ! Sc_p will be the original, or largest uniform cell size : needed for some division into blocks measurments
        Integer :: Dlamb,Nbc_p ! Nbc_p will be updated if adaptive_mesh==1
        Real(kind=8) :: pr_dx, pr_dy, pr_dz
        Real(kind=8) :: pr_xmin, pr_xmax,pr_ymin, pr_ymax 
        Real(kind=8) :: pr_zmin, pr_zmax
        Integer :: NbBl
    end type Particle 

    !! Definition of the object cell
    type Cell 
        Integer :: num_cell, num_block, num_diel,Dlamb_cell
        COMPLEX(real64) :: m_cell, Eps_cell
        COMPLEX(real64) ::  parameter_Ce, parameter_Sing
        Real(kind=8) :: parameter_Rad, parameter_Const
        Real(kind=8) :: lambda_cell,Sc, Xc, Yc, Zc
    end type Cell 
    
    type CBFM_Block
        Integer :: num_block,Nbc_b,Nbc_ext,NbBadj        
        Integer, Dimension(:), allocatable :: Badj
        Integer, Dimension(3):: gridPos
        Real(kind=8), Dimension(3,3) :: BCubCont ! for each direction x,y,z (columns) : val_min, val_max and h (rows)
        Real(kind=8), Dimension(4) :: BSphCont ! circumscribed sphere : center : (xc,yc,zc) and Radius
    end type CBFM_Block 

    !! CONSTANTS 
    Real(kind=8),Parameter ::Jr=0.;
    Real(kind=8),Parameter ::Ji=1.;
    Real(kind=8),Parameter :: Pi=Acos(-1.)
    Real(kind=8),Parameter :: offset=0.015
    Real(kind=8),Parameter :: ErrorEps = 1e-6
    COMPLEX(real64),Parameter :: J=(Jr,Ji)

    Real(kind=8),Parameter :: Eps0=8.8542E-12
    Real(kind=8),Parameter :: Rmu0=4.0E-7*Pi		!! Permeability of vaccum
    Real(kind=8),Parameter :: Eps_air=1            !! relative Permittivity of air
    Real(kind=8),Parameter :: Ro= 1E6              !(g/m3) Density of water
    Real(kind=8),Parameter :: C0= 3E8; !299792458;
    
    Integer, Parameter :: Round_Dp = 4   ! while generating the diameter of the particles 
                                         ! we keep 'Round_Dp' digit of precision (when Dp is expressed in mm) !
    Integer, Parameter :: Round_Sp = 6  ! same for Sc, always expressed in mm, we keep 3 digit of precison
    
    ! Parmeters used when discretizing the particles 
    Integer, Parameter :: NBc_max_per_part = 40000;
    
    ! Parmeters used when dividing the particles into blocks
    Integer, Parameter :: NbBlock_p_max = 10000;
    Integer, Parameter :: Fact_Nbext_max=10;  ! Nbext_max = Fact_Nbext_max*maxval(CBFM_Blocks(:)%Nbc_b);
    Integer,Parameter :: Nbcel_Blk_max = 2000;
    Real(kind=8),Parameter :: hB_test_step = (1./2.); !! this parameter is used when dividing into blocks; the next hB to test is equal to hB_test_step*current hB  
    integer,parameter,dimension(8) :: Nipws_f_rlamb = [91,190,190,231,325,496,703,861]   ! Number of Nipws (with uniform step on theta and Phi) required to ensure the accuracy 
                                                                                            ! of the CBFs depending on the ratio r/lambda_s (Fenni et al 2014) 
                                                                                            ! here we provide Nipws for 1< r/lambda_s <8
    integer,parameter,dimension(8) :: SphDes_Nipws_f_rlamb = [94,108,108,120,144,156,180,204]; ! this is a prelimineray table for the Nipws with spherical design
    integer,parameter,dimension(8) :: LebQuad_Nipws_f_rlamb = [86,110,110,146,170,194,230,266]; ! depending on hB
    
    ! Parameters used for the MPI parallelization 
    Integer, parameter :: tag = 1000
    Integer, Parameter :: M_B = 32; 
    Integer, Parameter :: N_B = 32;
    
End MODULE Initialization