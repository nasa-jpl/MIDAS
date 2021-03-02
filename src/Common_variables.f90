MODULE common_variables
        
    ! EM wave
    Real(kind=8) :: Freq_w,Omega_w,lambda_w,K_air
    Integer :: Nfreq,num_freq,freq_mag,lamb_mag
    Character(3) :: freq_unit
    Character(2) :: lamb_unit
    SAVE Nfreq, num_freq, Freq_w, freq_unit, freq_mag 
    SAVE Omega_w, lambda_w, lamb_unit, lamb_mag, K_air
    
    ! Transmitters/Receivers
    Integer :: NTrTheta, NTrPhi, NTr, NRxTheta, NRxPhi, NRx, NRx_tot
    Integer :: NTrTheta_CBFM, NTrPhi_CBFM, NTr_CBFM
    save NTrTheta, NTrPhi, NTr
    SAVE NRxTheta, NRxPhi, NRx, NRx_tot
    SAVE NTrTheta_CBFM, NTrPhi_CBFM, NTr_CBFM    
    
    ! Scatterer 
    Integer :: Nbc,homogs,Adapt_mesh,NbintBl
    SAVE Nbc,homogs,Adapt_mesh,NbintBl
     
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
    SAVE Nber_methods,leng_meth,CBFM, MLCBFM, MoM,RGE
    SAVE Dlambda,Nc_extended,Nbc_ext,Nipws,set_Nipws,distr_ipws
    SAVE hBlock,Nblocks,Nccp_max,Navg_cells
    SAVE DR,SR,fct_SR,res_SR,SR_Zc,SR_Zc_type,Eps_SR_Zc       
    SAVE NberLevels,NbBlksL2,Nber_IPWs_MLCBFM
    SAVE Use_ACA,Nb_it_max,Epsilon_ACA,vrb_ACA
    SAVE define_use_Copies,NcalBlks
    SAVE NumIntType_t, NumIntType_r
    
    ! Scattering Quantitites from Internal field
    Integer :: QextIFDisp
    SAVE QextIFDisp
    
    ! In/Output files
    Integer :: wr_Sij,wr_Qij,EqSph,shape_list
    Integer :: save_Zc,save_Eint,save_Eint_Nmax
    CHARACTER(100) :: Outfld_name
    CHARACTER(100) :: SimOutfld_name
    CHARACTER(240) :: ShapeFilePath
    SAVE Outfld_name, SimOutfld_name,EqSph,shape_list,ShapeFilePath
    SAVE wr_Sij,wr_Qij
    SAVE save_Zc,save_Eint,save_Eint_Nmax
    
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