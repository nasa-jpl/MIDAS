MODULE common_variables
        
    ! EM wave
    Real(kind=8) :: Freq_w,Omega_w,lambda_w,k_0
    Integer :: Nfreq,num_freq,freq_mag,lamb_mag
    Character(3) :: freq_unit
    Character(2) :: lamb_unit
    SAVE Nfreq, num_freq, Freq_w, freq_unit, freq_mag 
    SAVE Omega_w, lambda_w, lamb_unit, lamb_mag, k_0
    
    ! Transmitters/Receivers
    Integer :: NTrTheta, NTrPhi, NTr, NRxTheta, NRxPhi, NRx, NRx_tot
    Integer :: NTrTheta_CBFM, NTrPhi_CBFM, NTr_CBFM
    Integer :: sd_type, NPolBeta
    Real(kind=8) :: theta_init_trans_comp,theta_final_trans_comp
    Real(kind=8) :: phi_init_trans_comp,phi_final_trans_comp
    Real(kind=8) :: theta_init_Recei,theta_final_Recei,phi_init_Recei,phi_final_Recei
    Real(kind=8) :: beta_init_Pol, beta_final_Pol
    SAVE sd_type, NTrTheta, NTrPhi, NTr
    SAVE NRxTheta, NRxPhi, NRx, NRx_tot
    SAVE NTrTheta_CBFM, NTrPhi_CBFM, NTr_CBFM   
    SAVE theta_init_trans_comp,theta_final_trans_comp
    SAVE phi_init_trans_comp,phi_final_trans_comp
    SAVE theta_init_Recei,theta_final_Recei,phi_init_Recei,phi_final_Recei
    SAVE beta_init_Pol, beta_final_Pol, NPolBeta
    
    ! Far Field Approximation 
    Integer :: FFA
    Real(kind=8) :: Rso
    SAVE FFA, Rso
    
    ! Scatterer 
    Integer :: Nbc,homogs,NbintBl
    Integer :: Round_D,Round_S     
    CHARACTER(13) :: ap_str,lc_str, ac_str
    SAVE Nbc,homogs,NbintBl
    SAVE Round_D, Round_S       ! while generating the diameter of the particles  ENHANCEMENT needed here !!!! 
    SAVE ap_str,lc_str,ac_str   ! we keep 'Round_Dp' digit of precision (when Dp is expressed in mm) !
                                    ! same for Sc, always expressed in mm, we keep 3 digits of precison
     
    ! dielectric properties 
    CHARACTER(25) :: dielcomp_option
    Integer :: Ndiel
    SAVE dielcomp_option, Ndiel   
    
    ! Numerical methods parameters 
    Integer :: Nber_methods,leng_meth,CBFM, MLCBFM, MoM,RGE,ML_Extension
    Real(kind=8) :: Dlambda,hBlock,Epsilon_ACA,fct_SR,res_SR,Eps_SR_Zc 
    Integer :: Nblocks,Nccp_max,Navg_cells
    Integer :: Nipws,set_Nipws,distr_ipws,Nbc_ext,Nc_extended
    Integer :: NberLevels,NbBlksL2,Nber_IPWs_MLCBFM
    Integer :: Use_ACA,Nb_it_max,vrb_ACA,DR,SR,SR_Zc,SR_Zc_type
    Integer :: define_use_Copies,NcalBlks
    CHARACTER(2) :: NumIntType_t, NumIntType_r
    CHARACTER(3) :: div_type
    SAVE Nber_methods,leng_meth,CBFM, MLCBFM, MoM,RGE
    SAVE Dlambda,Nc_extended,Nbc_ext,Nipws,set_Nipws,distr_ipws
    SAVE div_type,hBlock,Nblocks,Nccp_max,Navg_cells
    SAVE DR,SR,fct_SR,res_SR,SR_Zc,SR_Zc_type,Eps_SR_Zc       
    SAVE NberLevels,NbBlksL2,Nber_IPWs_MLCBFM
    SAVE Use_ACA,Nb_it_max,Epsilon_ACA,vrb_ACA
    SAVE define_use_Copies,NcalBlks
    SAVE NumIntType_t, NumIntType_r
    
    ! Scattering Quantitites from Internal field
    Integer :: CextIntFields
    SAVE CextIntFields
    
    ! In/Output files
    Integer :: wr_Sij,wr_Qij,EqSph,shape_list
    Integer :: save_Zc,save_Eint,save_Eint_Nmax,save_Einc
    CHARACTER(100) :: Outfld_name
    CHARACTER(100) :: SimOutfld_name
    CHARACTER(240) :: ShapeFilePath
    SAVE Outfld_name, SimOutfld_name,EqSph,shape_list,ShapeFilePath
    SAVE wr_Sij,wr_Qij
    SAVE save_Zc,save_Eint,save_Eint_Nmax,save_Einc
    
    ! Environmemt
    CHARACTER(4) :: Env_type
    CHARACTER(1) :: Env_sep
    SAVE Env_type,Env_sep
    
    ! MPI Paralellization 
    Integer :: rank,nber_procs,INFO,code,Nbc_proc,Nblk_proc_max,MyNBlocks
    Integer :: Mlocal, Nlocal, NRHSlocal, Myrow, Mycol, NPROW, NPCOL
    SAVE rank,nber_procs,INFO,code,Nbc_proc,Nblk_proc_max 
    SAVE Mlocal, Nlocal, NRHSlocal, Myrow, Mycol, NPROW, NPCOL
    
    ! track memory performance
    Integer :: track_memory, Njob_max
    SAVE track_memory, Njob_max
       
    
END MODULE common_variables