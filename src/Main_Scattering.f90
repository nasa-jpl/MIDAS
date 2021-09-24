Program Main_Scattering
    
    USE Initialization
    USE common_variables
    USE iso_fortran_env
    USE MPI
    USE DiverseUtil
    
    Implicit none
    
    !! LOCAL *****************************************************************************************************************
    !! ***********************************************************************************************************************
    !! Transmitters/Receivers 
    Integer :: Ninc_sugg, Nscat_sugg, sd_type
    real(kind=8) :: theta_init_trans_comp,theta_final_trans_comp,Step_theta_trans_comp
    real(kind=8) :: phi_init_trans_comp, phi_final_trans_comp, step_phi_trans_comp
    real(kind=8) :: theta_init_Recei,theta_final_Recei,step_theta_Recei
    real(kind=8) :: phi_init_Recei,phi_final_Recei, step_phi_Recei   
    
    ! In/Output files
    character(200) :: file_name
    character(200) ::shape_folder_path
    
    ! Informations on CPU time 
    character(8)  :: date
    character(10) :: time
    character(5)  :: zone
    integer,dimension(8) :: values
    
    ! specific code Multi-frequency
    Integer :: Nwave, dd1,dd2
    logical :: dirExists
    Real(kind=8) :: freq_min, freq_max, Wave_min, Wave_max, step_wave, step_freq
    Real(kind=8), Dimension(:), allocatable :: Wavesle
    Complex, Dimension(:,:), allocatable :: m_lambdas
    
    ! Database code
    Integer :: eastat,numlines,NbSimulations,SimShape
    Real(kind=8) :: MaxDim, dim_ref, h_LargePart
    CHARACTER(60) ::  shapefile_path
    CHARACTER(:), allocatable ::  shapeSim_folder
    CHARACTER(300) ,allocatable::ShapesDirNames(:)
    CHARACTER(13) :: ap_str, lc_str, ac_str
    character(240) :: inputline
    logical :: fileExists    
    
    ! MPI 
    Integer :: countU
    Integer, Dimension(:,:), allocatable :: MPI_CBFM_Blocks, all_NbcBlocks
    Integer, Dimension(:), allocatable :: all_NBlocks
            
    ! others   
    Integer :: a,ii,jj,rr,Ind,I,K,m,ios,N_vals_m,Sim,old_Nbc
    Integer :: ind_line,tdistr_sca,Calc_EqSph,Nval_eps_r,Nval_eps_i
    Integer :: N,NBlks_exp,m_read_opt,err,Type_Par,pr_d,d,selected,num_bin,Nbins
    Real(kind=8) :: Volume,q, rp, ip,mrp , mip, p, Sc,Dp,h,ap,theta_dipole, phi_dipole
    Real(kind=8) :: x_l, y_l, z_l, xmax,xeq,xmax_m,xeq_m
    Real(kind=8) :: r_min,r_max,i_min,i_max, r_cyl, l_cyl,rp_min,rp_max
    Real(kind=8) :: xp,yp,zp,v_freq,wv,r_lambda,dmin,dmax,hselect,Deq_s,Dmax_s
    
    Character(9) :: SR_Zc_type_ch 
    Character(1) :: ch_tmp
    Character(11) :: freq_unit_tmp
    Character(10) :: lamb_unit_tmp
    Character(4) :: DataType
    Character(9) :: wave_descr,info_p_fl
    CHARACTER(:) ,allocatable::methods_names(:),stFreq
    CHARACTER(6) :: ty,tdata  
    character(250) :: m_file_name_0
    character(250),allocatable :: m_file_name(:)
    CHARACTER(300) analysis_fold_name
    CHARACTER(250) :: Sfold_name,Qfold_name,Solfold_name
    character(7) :: st_th_Tx,st_th_Rx,st_ph_Tx,st_ph_Rx
    CHARACTER(LEN=3) :: path
    
    type(Scatterer) :: SimScatterer
    type(Cell), Dimension(:),allocatable :: Cells, Upd_Cells
    type(Dipole), Dimension(:),allocatable :: Transmitters_Comp
    type(Dipole), Dimension(:),allocatable :: Receivers
    type(CBFM_Block), Dimension(:), allocatable :: CBFM_Blocks
    Integer, Dimension(7) :: Ncells_SphDomains
    Integer, Dimension(:,:), allocatable :: CBFM_Blocks_Ext,Upd_CBFM_Blocks_Ext, MLCBFM_BlDistr
    Integer, Dimension(:), allocatable :: Nbc_blocks,Diff_avg,NSims_bin,diel_comp_perc
    Real(kind=8), Dimension(:), allocatable :: hB_test,all_eps_r,all_eps_i,vals
    
    ! time 
    ! Time performances
    character(8)  :: date_init, date_final
    character(10) :: time_init, time_final
    character(5)  :: zone_init, zone_final    
    Integer,dimension(8) :: values_init, values_final
    Integer, dimension(4):: Comp_time_disc, Comp_time_write, Comp_time_div, Comp_time_ext
        
    !! END LOCAL *****************************************************************************************************************
    !! ***************************************************************************************************************************
    
    INTERFACE        
        SUBROUTINE Discretization(SimScatterer,Cells,Ncells_SphDomains)
            USE Initialization   
            USE common_variables
            USE MPI
            implicit none
            
            type (Scatterer), INTENT(INOUT) :: SimScatterer
            type (Cell), Dimension(:), allocatable, INTENT(OUT):: Cells
            Integer, Dimension(7), INTENT(OUT), OPTIONAL:: Ncells_SphDomains            
        END SUBROUTINE Discretization
        
        SUBROUTINE get_diel_values_lambdas(m_file_name,m_lambdas)    
            USE Initialization
            USE common_variables
            USE iso_fortran_env
            USE strings
    
            Implicit NONE
    
            ! IN/OUT 
            character(250), dimension(Ndiel), INTENT(IN) :: m_file_name
            Complex, Dimension(:,:), allocatable, INTENT(OUT) :: m_lambdas
            END SUBROUTINE get_diel_values_lambdas
        
        SUBROUTINE DielComposition(m_lambdas,Cells)    
            USE Initialization
            USE common_variables
            USE iso_fortran_env
        
            Implicit NONE
    
            ! IN/OUT 
            Complex, Dimension(Ndiel,Nfreq), INTENT(IN) :: m_lambdas
            type (Cell), Dimension(Nbc), INTENT(INOUT):: Cells            
        END SUBROUTINE DielComposition
        
        SUBROUTINE SetCellsParams(SimScatterer,Cells,Upd_Cells)    
            USE Initialization
            USE common_variables
            USE iso_fortran_env
            Implicit NONE
    
            ! IN/OUT
            type (Scatterer), INTENT(INOUT) :: SimScatterer
            type (Cell), Dimension(Nbc), INTENT(INOUT):: Cells
            type (Cell), Dimension(:), allocatable, INTENT(OUT):: Upd_Cells
        END SUBROUTINE SetCellsParams
        
        SUBROUTINE Division_blocks(SimScatterer,Cells,Ncells_SphDomains,CBFM_Blocks,MLCBFM_BlDistr)
    
            USE Initialization
            USE common_variables
            USE MPI
    
            !IN/OUT 
            type (Scatterer), INTENT(INOUT) :: SimScatterer
            type (Cell), Dimension(Nbc), INTENT(INOUT) :: Cells
            Integer, Dimension(7), INTENT(IN) :: Ncells_SphDomains
            type (CBFM_Block), Dimension(:), allocatable, INTENT(OUT) :: CBFM_Blocks
            Integer, Dimension(:,:), allocatable, INTENT(OUT):: MLCBFM_BlDistr
        END SUBROUTINE Division_blocks
        
        SUBROUTINE Extend_blocks(SimScatterer,Cells,CBFM_Blocks,CBFM_Blocks_Ext)

            USE Initialization
            USE common_variables
            !USE f95_precision
            USE MPI
        
            Implicit NONE
        
            !IN/OUT 
            type(Scatterer), INTENT(IN) :: SimScatterer
            type(Cell), Dimension(Nbc), INTENT(IN) :: Cells
            type(CBFM_Block), Dimension(Nblocks), INTENT(INOUT) :: CBFM_Blocks
            Integer, Dimension(:,:), allocatable, INTENT(OUT) :: CBFM_Blocks_Ext           
        END SUBROUTINE Extend_blocks 
        
        SUBROUTINE get_trans_Receiv(Ninc_in,Nscat_in,sd_type,thitrans,thftrans,phitrans,phftrans,thiRecei,thfRecei,phiRecei,phfRecei,Transmitters_Comp,Receivers);
            USE Initialization
            USE common_variables
    
            IMPLICIT NONE
    
            !! IN/OUT ******************************************************************
    
            Integer, INTENT(IN) :: Ninc_in,Nscat_in
            Integer, INTENT(OUT) :: sd_type
            Real(kind=8), INTENT(IN) :: thftrans,thitrans,phftrans,phitrans 
            Real(kind=8), INTENT(IN) :: thfRecei,thiRecei,phfRecei,phiRecei
            type (Dipole), Dimension(:), allocatable, INTENT(OUT) :: Transmitters_Comp
            type (Dipole), Dimension(:), allocatable, INTENT(OUT) :: Receivers                   
        END SUBROUTINE get_trans_Receiv
    
        SUBROUTINE MPI_distribution_blocks(CBFM_Blocks,MPI_CBFM_Blocks)

            ! HERE Distribution of the CBFM blocks among the available MPI jobs        
            USE Initialization
            USE common_variables
            USE iso_fortran_env
            USE DiverseUtil
            USE MPI
    
            Implicit NONE    
            !IN/OUT 
            type (CBFM_Block), Dimension(Nblocks), INTENT(IN) :: CBFM_Blocks
            Integer, Dimension(:,:), allocatable, INTENT(INOUT) :: MPI_CBFM_Blocks
        END SUBROUTINE MPI_distribution_blocks
    END INTERFACE
    
    !! ---------------------------------------------------------------------------------------------------------------------------!!
    !! ---------------------------------------------------- MAIN PROGRAM ---------------------------------------------------------!!
    !! ---------------------------------------------------------------------------------------------------------------------------!!
    
    ! MPI INITIALIZATION
    Call MPI_INIT(code);
    Call MPI_COMM_SIZE (MPI_COMM_WORLD,nber_procs,code)
    Call MPI_COMM_RANK (MPI_COMM_WORLD,rank,code)        
    
    !! ---------------------------------------------------------------------------------------------------------------------------!!
    !! ---------------------------------------------------- MAIN PROGRAM ---------------------------------------------------------!!
    !! ---------------------------------------------------------------------------------------------------------------------------!!
    
    !! for now for the MPI code I will not use the spherical shape option 'define & use Copies', I will see if interesting later
    define_use_Copies = 0;
    ! If IFDisp == 1, we will keep and display the Extenction cross section
    ! calculated from the Internal Field 
    QextIFDisp = 0; 
    
    ! track memory use 
    track_memory = 1; 
    Njob_max = 100;
    
    !! HERE ALL THE PROCS WILL READ THE SAME SIMULATION INPUT FILES :
    !! Reading the data file *****************************************************************************************************
    !! ***************************************************************************************************************************    
    Call GET_ENVIRONMENT_VARIABLE('PATH',path);
    env_sep = path(3:3);
    if (Env_sep .eq. '\') Then 
        Env_type = 'WIND';
    else
        Env_type = 'LINU';
        Env_sep = '/';
    endif
    
    !! Open the dat file and read the simulation parameters.    
    file_name = 'inputs'//Env_sep//'Simulation_data.dat'
    Open(11,File = file_name) 
    Read(11,*);Read(11,*)
    
    
    !! Parameters of the EM wave
    read(11,*)
    read(11,*), wave_descr !! wave description : wavelength or frequency     
     
    if (wave_descr == 'NFreq') Then 
        read(11,'(i3)'), Nfreq
        read(11,'(a6)'), tdata
        read(11,'(a)'), freq_unit_tmp
        if (Len(trim(freq_unit_tmp)) .eq. 11) then
            freq_unit = freq_unit_tmp(8:10);                
        else 
            If (rank .eq. 0) Then
                Write(*,'(a)') 'Error : Enable to read frequency unit !'
            endif
            go to 30;
        endif
        if (freq_unit == 'MHz') then
            lamb_unit = 'm'; 
            freq_mag = 6; lamb_mag = 0;                         
        elseif (freq_unit == 'GHz') then
            lamb_unit = 'mm';
            freq_mag = 9; lamb_mag = 3;
        elseif (freq_unit == 'THz') then
            lamb_unit = 'um';
            freq_mag = 12; lamb_mag = 6;
        else
            If (rank .eq. 0) Then
                Write(*,'(a)') 'Error : Invalid frequency unit !'
            endif
            go to 30;
        endif           
        
        Nwave = Nfreq; Allocate(Wavesle(Nwave));            
        if (tdata == 'minmax') then
            read(11,'(f7.3,a,f7.3)'), freq_min ,ch_tmp, freq_max
            !wave_min = C0/(freq_max*1E6); wave_max = C0/(freq_min*1E6);
            !wave_min = C0/(freq_max*1E6); wave_max = C0/(freq_min*1E6);  
            If (Nfreq .gt. 1) then           
              step_freq = (freq_max-freq_min)/(Nfreq-1)  
            Else
              step_freq = 0;
            EndIf                 
            Do ii =1,Nfreq
                Wavesle(ii) = C0/((10.**freq_mag)*(freq_min + (ii-1)*step_freq))   
            EndDo
        else
            If (Nfreq == 1) Then 
                read(11,'(f7.3,a)'), v_freq 
                Wavesle(1) = C0/(v_freq*(10.**freq_mag))*(10.**lamb_mag);    ! Wavesle is expressed inside the code in m if MHz, mm if GHz, um if THz    
            Else                
                Do ii=1, Nfreq-1
                  read(11,'(f7.3,a)',advance='no'), v_freq ,ch_tmp
                  Wavesle(ii) = C0/(v_freq*(10.**freq_mag))*(10.**lamb_mag);    
                EndDo
                read(11,'(f7.3)',advance='no'), v_freq 
                Wavesle(ii) = C0/(v_freq*(10.**freq_mag))*(10.**lamb_mag);      
                read(11,*)
            EndIf            
        endif       
    Else
        read(11,'(i3)'), Nwave
        read(11,'(a6)'), tdata      
        read(11,'(a)'), lamb_unit_tmp
        if (Len(trim(lamb_unit_tmp)) .eq. 9) then
            lamb_unit = trim(lamb_unit_tmp(8:8));             
        elseif (Len(trim(lamb_unit_tmp)) .eq. 10) then
            lamb_unit = trim(lamb_unit_tmp(8:9));
        else
            If (rank .eq. 0) Then
                Write(*,'(a)') 'Error : Enable to read wavelength unit !'
            endif
            go to 30;
        endif
        if (lamb_unit == 'm') then
            freq_unit = 'MHz'; 
            freq_mag = 6; lamb_mag = 0;                         
        elseif (lamb_unit == 'mm') then
            freq_unit = 'GHz';
            freq_mag = 9; lamb_mag = 3;
        elseif (lamb_unit == 'um') then
            freq_unit = 'THz';
            freq_mag = 12; lamb_mag = 6;
        else
            If (rank .eq. 0) Then
                Write(*,'(a)') 'Error : Invalid wavelength unit !'
            endif
            go to 30;
        endif  
        
        Nfreq = Nwave; Allocate(Wavesle(Nwave));        
        if (tdata == 'minmax') then
            read(11,*), wave_min ;read(11,*), wave_max 
            freq_min = C0/(wave_max*10.**(-lamb_mag)); 
            freq_max = C0/(wave_min*10.**(-lamb_mag));
            If (Nwave .gt. 1) then           
              step_wave = (wave_max-wave_min)/(Nwave-1)  
            Else
              step_wave = 0;
            EndIf                     
            Do ii =1,Nwave
                Wavesle(ii) = wave_max - (ii-1)*step_wave
            EndDo   
        else
            if (Nwave == 1) Then 
                read(11,'(f7.3)'), Wavesle(1);
            else
                Do ii=1, Nwave-1
                  read(11,'(f7.3,a)',advance='no'), Wavesle(Nwave-ii+1) ,ch_tmp
                EndDo 
                read(11,'(f7.3)',advance='no'), Wavesle(Nwave-ii+1)
                read(11,*)
            EndIf                        
        endif       
    EndIf
    
    read(11,*)    
    !! Parameters of the scatterer 
    read(11,*)
    read(11,'(i1)');read(11,*), shape_list
    read(11,'(a)'), shape_folder_path
    read(11,'(a)'), Outfld_name
    
    if (shape_list == 0) then 
      ! type_p
      read(11,*);
      !read(11,'(i1,a,a)'), SimScatterer%type_s,SimScatterer%info_s, type_shape_in
      read(11,'(i1,a)'), SimScatterer%type_s, inputline
      ind_line = INDEX(inputline,' ');
      SimScatterer%info_s = inputline(1:ind_line-1); 
      ind_line = INDEX(inputline,'cells');
      if (ind_line .ne. 0) then 
        SimScatterer%ty_shape_in = 'cells'
      else
        SimScatterer%ty_shape_in = 'shape'
      endif
            
      read(11,*);
      if (SimScatterer%type_s .eq. 3) then ! for the moment the only different type in reading param is the cylinder : we read a and L
          read(11,*), ac_str, lc_str
          read(ac_str,*), r_cyl
          read(lc_str,*),l_cyl ; ! (m or mm or um)
          SimScatterer%a = r_cyl/(10**lamb_mag)
          SimScatterer%dm = 2*r_cyl/(10**lamb_mag)
          SimScatterer%dy = SimScatterer%dm
          SimScatterer%dz = SimScatterer%dm
          SimScatterer%dx = l_cyl/(10**lamb_mag)
          SimScatterer%info_s ='Cylin';
      else          
          ! ap
          read(11,'(a)'), ap_str
          read(ap_str,*), ap
          SimScatterer%a = ap/(10**lamb_mag)
          SimScatterer%dm = 2*ap/(10**lamb_mag)
      endif
    else
      read(11,*); read(11,*);
      read(11,*); read(11,*);     
    EndIf
    
    ! Eps_p
    read(11,*);
    read(11,*), dielcomp_option, Ndiel 
    if (trim(dielcomp_option) == 'fromshapefile') then 
        Allocate(m_file_name(Ndiel),diel_comp_perc(Ndiel)); ! Ndiel is releavant for this option, diel_comp_perc can be caluclated once dielc composition read from shape.dat
        DO ii=1,Ndiel
            read(11,*), m_file_name_0
            m_file_name(ii) = trim(m_file_name_0);           
        Enddo
    elseif (trim(dielcomp_option) == 'fromonlymfile') then 
        Ndiel = 1;
        Allocate(m_file_name(Ndiel),diel_comp_perc(Ndiel)); ! Ndiel is simply equal to 1 here. one m per frequency !
        DO ii=1,Ndiel
            read(11,*), m_file_name_0
            m_file_name(ii) = trim(m_file_name_0);           
        Enddo
    elseif (trim(dielcomp_option) == 'random1') then
        Ndiel = 2;
        Allocate(m_file_name(Ndiel),diel_comp_perc(Ndiel)); ! because of the totally random process Ndiel is in theory =Nbc and diel_comp_perc is not releavant here  
        Do ii=1,2                                           
            read(11,*), m_file_name_0
            m_file_name(1) = trim(m_file_name_0);
        Enddo        
    elseif (trim(dielcomp_option) == 'random2') then
        Allocate(m_file_name(Ndiel),diel_comp_perc(Ndiel));
        DO ii=1,Ndiel
            read(11,*), m_file_name_0,p
            m_file_name(ii) = trim(m_file_name_0);  
            diel_comp_perc(ii) = p; 
        Enddo
    elseif (trim(dielcomp_option) == 'fromdielcompositionfile') then
        Allocate(m_file_name(1),diel_comp_perc(1)); ! Since each of the Nbc cell has a different refractive index here, diel_comp_perc is not releavant here 
        m_file_name(1) = 'inputs/dielcomposition.dat'; ! this file contains the refractive index per cell        
    else
        Write(*,'(a,a,a)')'Error : ', dielcomp_option, ' is an unknown dielectric decomposition option !!';
        stop 1
    endif    
        
    if (((trim(dielcomp_option) == 'fromdielcompositionfile') .OR. (trim(dielcomp_option) == 'fromshapefile')) &
        .AND. (SimScatterer%type_s .ne. 2)) then 
        Write(*,'(a,a)') 'Error : The requested dielectric decomposition option can only be used with type_scatterer = 2';
        stop 1        
    Endif     
    
    read(11,*); read(11,'(i1)'), Calc_EqSph;
    read(11,*);
        
    !! Parameters of the dicretization 
    read(11,*)
    read(11,*), ch_tmp;
    if (ch_tmp .eq. 'D') then
      read(11,*), Dlambda
    elseif (ch_tmp .eq. 'S') then
      read(11,*), Sc;
      if ((freq_unit == 'THz') .OR. (freq_unit == 'GHz')) then 
        SimScatterer%Sc = 1e-6*Sc;  
      else
        SimScatterer%Sc = Sc;
      endif
            
      Dlambda = 1;
    else
        If (rank ==0) Then
            Write(*,'(a)') 'Error when reading discretization parameter'
        endif
        go to 30;           
    endif
    read(11,*); read(11,*), Adapt_mesh
    read(11,*)
    
    !!!! Parameters of the Applied methods
    read(11,*)
    read(11,*);read(11,*), Nber_methods
    read(11,*);read(11,*), leng_meth
    Allocate(character(leng_meth) :: methods_names(Nber_methods))
    read(11,*)
    CBFM=0; MLCBFM=0; MoM=0; RGE=0;
    DO ii=1,Nber_methods
        read(11,*), methods_names(ii)
        !! The integers CBFM; MLCBFM and MoM represent the position of each method in the array methods_names 
        !! if this method is applied, and is equal to 0 otherwise        
        If (methods_names(ii) == 'CBFM-E') Then
            CBFM = ii
        Elseif (methods_names(ii) == 'MLCBFM-E') Then
            MLCBFM = ii
        Elseif (methods_names(ii) == 'RGE') Then
            RGE = ii   
        Elseif (methods_names(ii) == 'MoM') Then
            MoM = ii        
        Endif     
    Enddo
    read(11,*)
     
    !! Reading Transmitters ****************************************************************************************************************************************
    read(11,*)
    read(11,*);read(11,'(a2,a1,a2)'), NumIntType_t,ch_tmp,NumIntType_r
    if (NumIntType_r .eq. '') then 
        NumIntType_r = NumIntType_t;
    endif
    read(11,*);read(11,*), Ninc_sugg
    read(11,*);read(11,*), Nscat_sugg
    read(11,*);
    read(11,*);read(11,*),theta_init_trans_comp,theta_final_trans_comp,NTrTheta
    read(11,*);read(11,*),phi_init_trans_comp,phi_final_trans_comp,NTrPhi
    read(11,*)

    !! Reading Receivers ********************************************************************************************************
    read(11,*);read(11,*),theta_init_Recei,theta_final_Recei,NRxTheta
    read(11,*);read(11,*),phi_init_Recei,phi_final_Recei,NRxPhi
    
    ! Write Scattering matrix elements for each incident direction and Q per incident direction 
    read(11,*);read(11,*), wr_Sij
    read(11,*);read(11,*), wr_Qij    
    read(11,*)
    
    ! Get Transmitters/Scatterers depending on the type of the numerical integration used to average the scattering quantities
    ! over incident/scattering directions
    call get_trans_Receiv(Ninc_sugg,Nscat_sugg,sd_type,theta_init_trans_comp,theta_final_trans_comp,&
                          phi_init_trans_comp,phi_final_trans_comp, &
                          theta_init_Recei,theta_final_Recei, phi_init_Recei,phi_final_Recei,&
                          Transmitters_Comp,Receivers);
    
    ! if only 1 incident direction is used, we autmatically put wr_Sij and wr_Qij to 1
    if (NTr .eq. 1) then 
        wr_Sij=1;wr_Qij=1;
    endif
    
    ! Parameters of the numerical methods 
    
    !CBFM
    read(11,*)
    read(11,*);read(11,*), Navg_cells
    read(11,*);read(11,*), set_Nipws   !! if set_Nipws we will use setNipws in getParameters_CBFM.f90
    read(11,*);read(11,*), distr_ipws !! type of distribution for the N incident plane waves used 
                                      !! to generate the CBFs (see getTransmitters_CBFM for details)
    read(11,*);read(11,*), Nc_extended
    read(11,*);read(11,*), DR
    read(11,*);read(11,*), SR
    read(11,*), res_SR
    read(11,*);read(11,*), SR_Zc
    read(11,*) SR_Zc_type_ch
    read(11,*), Eps_SR_Zc
    read(11,*)        
         
    !ACA
    read(11,*);
    read(11,*);read(11,*), Use_ACA
    read(11,*);read(11,*), Nb_it_max                                                                                                                                                                                     
    read(11,*);read(11,*), Epsilon_ACA                                                                                                                                                                                           
    read(11,*);read(11,*), Vrb_ACA
    read(11,*) 
    
    ! decide SR_Zc_type from SR_Zc_type_ch
    If (trim(SR_Zc_type_ch) =='threshold') Then
        SR_Zc_type = 1; 
    ElseIf (trim(SR_Zc_type_ch) =='edistance') Then
        SR_Zc_type = 2;
    ElseIf (trim(SR_Zc_type_ch) =='spalgo_dz') Then
        SR_Zc_type = 3;
    EndIf
    
    If ((Use_ACA == 1) .and.(SR_Zc==1) .and. (SR_Zc_type .ne. 1)) Then  !! The use of ACA is available with only the first sparsity approach
        If (rank == 0) Then   
            Write(*,'(a)') 'ERROR : Wrong combination UseACA/Sparsity !!! Exit !!'            
        endif          
        go to 30;  !! stop 1 points any error in the chosed/entred simulation data           
    EndIf
    
    ! Save Sol Elements
    read(11,*);
    read(11,*);read(11,*), save_Zc
    read(11,*);read(11,*), save_Eint
    if (save_Eint == 1) then
        read(11,*), save_Eint_Nmax
    endif        
    
    !! close the dat file
    Close(11)
    
    !! Initialization ************************************************************************************************************
    !! ***************************************************************************************************************************
    

    If (shape_list .eq. 1) then
        ! First all the jobs will wait until job 0 check the existence and create if needed the SimShape.dat file
        if (rank == 0) then 
            inquire(file='inputs/SimShapes.dat', exist=fileExists)
            if (.not. fileExists) then
                if (Env_type .eq. 'WIND') Then 
                    CALL SYSTEM('dir /s /b inputs/shape.dat >> inputs/SimShapes.dat');              
                else
                    CALL SYSTEM('find '//trim(shape_folder_path)//' -type f -name inputs/shape.dat> inputs/SimShapes.dat');              
                endif         
            EndIf
        EndIf
        Call MPI_Barrier(MPI_COMM_WORLD,code);    
      
        ! Count Simulations ! To count the simulation files in the folder shape (Kuo Database) we search  Env_sep//'size'//Env_sep
        ! THESE CODE LINES WORK ONLY WITH KWO DATABASE SHAPE FOLDER STRUCTURE
        Open(12,File = 'inputs/SimShapes.dat') 
        eastat = 0; numlines =0;
        Do while (eastat .ge. 0)
            READ(12,'(a)',iostat=eastat) inputline
            ii = index(inputline,Env_sep//'size'//Env_sep);
            if ((eastat .ge. 0) .and. (ii .ne. 0)) then 
                numlines = numlines + 1
            endif                
        EndDo
        close(12);      
        NbSimulations = numlines;
        if (rank == 0) then 
            Write(*,'(a,i6)') 'The total number of Shape Simulations : ', NbSimulations
        endif
        Allocate(ShapesDirNames(NbSimulations));
        Open(12,File = 'SimShapes.dat') 
      
        ! Recover shape files pathes 
        ii = 1;
        Do while (ii .le. NbSimulations)
            READ(12,'(a)',iostat=eastat) inputline
            if (index(inputline,Env_sep//'size'//Env_sep) .ne. 0) then 
            ShapesDirNames(ii)= inputline;
            ii = ii + 1;
            endif
        EndDo
        close (12);     
    Else
            NbSimulations = 1;
    endif
    
    ! HERE START SCATTERER  
    Do Sim=1, NbSimulations
        If (shape_list .eq. 1) then 
            ShapeFilePath = ShapesDirNames(Sim);
            ii = index(ShapeFilePath,'shape'//Env_sep);
            SimScatterer%type_s = 2;
            if (ShapeFilePath(ii+6:ii+6) .eq. 'a') then 
                SimScatterer%info_s = ShapeFilePath(ii+6:ii+6)//ShapeFilePath(ii+8:ii+11);
                read(ShapeFilePath(ii+18:ii+30),'(a)') ap_str;
                read(ShapeFilePath(ii+18:ii+30),'(f13.6)') ap;
            elseif (ShapeFilePath(ii+6:ii+6) .eq. 'p') then 
                SimScatterer%info_s = ShapeFilePath(ii+6:ii+6)//ShapeFilePath(ii+8:ii+9);
                read(ShapeFilePath(ii+16:ii+28),'(a)') ap_str;
                read(ShapeFilePath(ii+16:ii+28),'(f13.6)') ap;
            else
                if (rank == 0) then
                    Write(*,'(a)') 'Something went wrong when reading ap from SimShapes.dat';
                endif                
                go to 30; 
            endif
            
            ap = ap/(10**lamb_mag); ! ap (m/mm/um)
            SimScatterer%a = ap/(10**lamb_mag)
            SimScatterer%dm = 2*ap/(10**lamb_mag)
            ap_str(1:1) = ap_str(3:3); ap_str(2:2)='.';
            ap_str(3:11) = ap_str(4:6)//ap_str(8:13);ap_str(12:13) = ''; !to prepare the name of the output folder
            
            if (rank == 0) Then 
                Write(*,'(a)')' '
                Write(*,'(a)')  '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++';
                Write(*,'(a,i5,a,i5,a)') '++ SIM ',Sim,' OUT OF ',NbSimulations,' ++++++++++++++++++++++++++++++++';
                Write(*,'(a)') '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++';
                Write(*,'(a)')' '
            endif            
        Else
            if (SimScatterer%ty_shape_in == 'cells') then 
                ShapeFilePath = 'inputs'//Env_sep//'Cells.dat'; 
            else
                ShapeFilePath = 'inputs'//Env_sep//'shape.dat';
            endif         
        Endif
        ! Scatterer Output Folder
        !Write(ap_str,'(f11.9)') ap;
        If (shape_list .eq. 1) then 
            if ((SimScatterer%type_s == 2) .OR. (SimScatterer%type_s == 6)) then 
                SimOutfld_name = trim(Outfld_name)//Env_sep//trim(SimScatterer%info_s)//'-ap='//trim(ap_str)//lamb_unit; 
            elseif (SimScatterer%type_s == 3) then
                SimOutfld_name = trim(Outfld_name)//Env_sep//trim(SimScatterer%info_s)//'-ac='//trim(ac_str)//lamb_unit//'-lc='//trim(lc_str)//lamb_unit;
            elseif (SimScatterer%type_s == 1) then
                SimOutfld_name = trim(Outfld_name)//Env_sep//'Sphere-ap='//trim(ap_str)//lamb_unit; 
            endif
        else
            if ((SimScatterer%type_s == 2) .OR. (SimScatterer%type_s == 6)) then 
                SimOutfld_name = trim(Outfld_name)//Env_sep//trim(SimScatterer%info_s)//'-ap='//trim(ap_str)//lamb_unit;
            elseif (SimScatterer%type_s == 3) then
                SimOutfld_name = trim(Outfld_name)//Env_sep//trim(SimScatterer%info_s)//'-ac='//trim(ac_str)//lamb_unit//'-lc='//trim(lc_str)//lamb_unit;
            elseif (SimScatterer%type_s == 1) then
                SimOutfld_name = trim(Outfld_name)//Env_sep//'Sphere-ap='//trim(ap_str)//lamb_unit;
            endif               
        endif        
        if (rank == 0) Then
            inquire(directory=trim(SimOutfld_name),exist=dirExists);        
            if (dirExists) Then
                call system('rm -r '//trim(SimOutfld_name))        
            EndIf
            call system('mkdir "'//trim(SimOutfld_name)//'"'); 
            if (Env_type == 'WIND') then
                call system('copy inputs\Simulation_data.dat "'//trim(SimOutfld_name)//'\"'); 
            else
                call system('cp inputs/Simulation_data.dat "'//trim(SimOutfld_name)//'/"'); 
            endif 
        endif
        Call MPI_Barrier(MPI_COMM_WORLD,code);  ! here all the jobs wait for the creation of the simulation folders 
        
        ! create analysis folder & store first moment (time reference)
        analysis_fold_name = trim(SimOutfld_name)//Env_sep//'Analysis';
        if (rank .eq. 0) then             
            CALL SYSTEM('mkdir "'//trim(analysis_fold_name)//'"');
        endif
        
        Do ii = 1,4
          Call MPI_Barrier(MPI_COMM_WORLD,code); 
          if (rank == ii) then
              ! create the S_files folder if needed
              inquire(directory=trim(analysis_fold_name),exist=dirExists);
              if (.not. dirExists)  Then
                  CALL SYSTEM('mkdir "'//trim(analysis_fold_name)//'"');
              EndIf
          endif
        EndDo   
        call print_allocate(14,'Reference(t=0)','NONE ',0);
    
        !! ICI GENERATION SCENE DE SIMULATION ET DISCRETISATION -------> Resultat : Nbc et tableau cellules
        !! PREPARING THE SIMULATION SCENE (the same for all the frequencies)
        ! To simplify, Let us discretize the simulation scene according to the higher considered frequency***************************
        ! thus the simulation scene is discretized and divided into blocks once !
        Lambda_w = maxval(Wavesle(:))/(10**lamb_mag)
        Freq_w = C0/Lambda_w;
        Omega_w = 2*Pi*Freq_w         !! angular frequency
        k_0 = (2*Pi)/Lambda_w 
    
        !! here we read a first time the m files to obtain, depending on the diel decomposition options 
        !! just to initialize SimScatterer%lambda. This code line is usefull if there is generation/discretization of sphere/cylinder/chebychev part ..., It
        !! is useless if the geometry is read from a shape file !
        !! Usefull for the sphere, cylinder or chebychev/GRD particle to be able to accurately discretize according to the higher frequency/refractive index, and this f
        call get_diel_values_lambdas(m_file_name,m_lambdas);
    
        ! to discretize for a multi-frequency simulation (applicable for type_part .ne. 2), for the moment we take into account the shortest lambda
        ! even if now the discretization is needed based on the geometrical modeling of the scatterer (from Kuo)
        rp_min = 1e2; rp_max = 0;
        SimScatterer%lambda_min = Lambda_w;
        SimScatterer%lambda_max = 0;    
        Do jj=1,Nfreq
            Do ii =1,Ndiel
                mrp = real(m_lambdas(ii,jj)) 
                mip = imag(m_lambdas(ii,jj))    
                rp = mrp**2-mip**2;  
                ip = 2*mrp*mip;
                
                if (rp .ge. rp_max) then 
                    rp_max = rp;
                    SimScatterer%Eps_max = rp+J*ip
                    SimScatterer%m_max = m_lambdas(ii,jj);
                endif 
                
                if (Wavesle(jj)/(10**lamb_mag*sqrt(rp)) .le. SimScatterer%lambda_min) then 
                    SimScatterer%lambda_min = Wavesle(jj)/((10**lamb_mag)*sqrt(rp));
                endif     
                
                if (rp .le. rp_min) then 
                    rp_min = rp;
                    SimScatterer%Eps_min = rp+J*ip
                    SimScatterer%m_min = m_lambdas(ii,jj);
                endif                 
                if (Wavesle(jj)/((10**lamb_mag)*sqrt(rp)) .ge. SimScatterer%lambda_max) then 
                    SimScatterer%lambda_max = Wavesle(jj)/((10**lamb_mag)*sqrt(rp));
                endif 
            EndDo       
        EndDo
        if ((trim(dielcomp_option) == 'fromonlymfile') .OR. &
            ((trim(dielcomp_option) == 'fromshapefile') .AND. (Ndiel .eq.1))) then
            homogs = 1;
        else
            homogs = 0;
        EndIf 
        
        !*****************************************************************************************************
        If (NbSimulations .eq. 1) then 
            if (rank == 0) then 
                Write (*,'(a)') '**************************************************************************************'
                Write (*,'(a)') '************* Computing of the Scattering by Complex-Shaped Scatterer *****************'
                Write (*,'(a)') '************************** CODE VIEM_MoM-CBFM_VoxelMesh ******************************'
                Write (*,'(a)') '*************************************************************************************'
            endif
        EndIf
        
        ! Simulation file   
        if (rank == 0) then 
    10      call date_and_time(date,time,zone,values);
            if (EqSph==0) then
            	if ((SimScatterer%type_s == 2) .OR. (SimScatterer%type_s == 6)) then 
                	file_name = trim(SimOutfld_name)//Env_sep//'Simulation_'//trim(SimScatterer%info_s)//'_'//date(5:6)//&
                    '-'//date(7:8)//'-'//date(1:4)//'_'//time(1:2)//'h'//time(3:4)//'.dat';   
            	elseif (SimScatterer%type_s == 1) then 
                	file_name = trim(SimOutfld_name)//Env_sep//'Simulation_Sphere_'//date(5:6)//&
                	'-'//date(7:8)//'-'//date(1:4)//'_'//time(1:2)//'h'//time(3:4)//'.dat';  
            	else
                	write(ch_tmp,'(i1)') SimScatterer%type_s;
                	file_name = trim(SimOutfld_name)//Env_sep//'Simulation_ty'//trim(ch_tmp)//'_'//date(5:6)//&
                	'-'//date(7:8)//'-'//date(1:4)//'_'//time(1:2)//'h'//time(3:4)//'.dat';  
            	endif        
            EndIf            
            Open(10,File = trim(file_name));

            Write (10,'(a)') '**************************************************************************************'
            Write (10,'(a)') '************* Computing of the Scattering by Complex-Shaped Scatterer *****************'
            Write (10,'(a)') '************************** CODE VIEM_MoM-CBFM_VoxelMesh ******************************'
            Write (10,'(a)') '*************************************************************************************'      
    
    
            !! Once generated, all these informations should be written in the output file and Simulation_data_out 
            Write(10,*) '' 
            write (10,'(a,i3,a)',advance='no') 'The number of frequencies = ', Nfreq,' : ['
            Do ii=1, Nfreq-1   
                Freq_w = C0/(Wavesle(ii)/(10**lamb_mag));
                a = nint(Freq_w/(10**freq_mag));
                if (a < 10) Then 
                    Allocate(character(4) ::stFreq)
                    ty = '(f4.2)';
                ElseIf (a < 100) Then
                    Allocate(character(5) ::stFreq)
                    ty = '(f5.2)';
                Else
                    Allocate(character(6) ::stFreq)
                    ty = '(f6.2)';
                EndIf    
                Write(stFreq,ty) Freq_w/(10**freq_mag)
                write (10,'(a,a)',advance='no') stFreq,'; '        
                Deallocate(stFreq);
            EndDo
            Freq_w = C0/(Wavesle(ii)/(10**lamb_mag)); 
            a = nint(Freq_w/(10**freq_mag));
            if (a < 10) Then 
                Allocate(character(4) ::stFreq)
                ty = '(f4.2)';
            ElseIf (a < 100) Then
                Allocate(character(5) ::stFreq)
                ty = '(f5.2)';
            Else
                Allocate(character(6) ::stFreq)
                ty = '(f6.2)';
            EndIf    
            Write(stFreq,ty) Freq_w/(10**freq_mag)  
            write (10,'(a,a,a)') stFreq,'] ',trim(freq_unit);
            Deallocate(stFreq);
    
            if (EqSph==1) then
                Write (*,'(a)') ' '
                Write (*,'(a)') '******************** VOLUME EQUIVALENT SPHERE SIMULATIONS **********************';
                Write (*,'(a)') '********************************************************************************';
                Write (*,'(a)') ' ';
                Deallocate(Cells,CBFM_Blocks,CBFM_Blocks_Ext,MLCBFM_BlDistr);
            endif
        endif
    
        if (EqSph==1) then
            Deallocate(Cells,CBFM_Blocks,CBFM_Blocks_Ext,MLCBFM_BlDistr);
        endif
        
        !! DISCRETIZATION & DIVISION INTO BLOCKS**************************************************************************
        ! Discretization
        Type_Par = SimScatterer%type_s
        info_p_fl = trim(SimScatterer%info_s);
        !Comp_time_disc = 0; values_init=0; values_final =0;
        call date_and_time(date_init,time_init,zone_init,values_init); 
        
        Call Discretization(SimScatterer,Cells,Ncells_SphDomains); 
        !! Remember that ap and dp refers to effective radius. ceci corrige quand necessaire ou garde la meme valeur si c bon
        if (Type_Par .ne. 1) then 
        	!SimScatterer%a = ((3*Nbc*SimScatterer%Sc**3.)/(4*pi))**(1./3.)  
            ! instead of N*d^3 we need to use sum(d^3) in case we are using different Sc (for adaptive mesh for example)!!
            SimScatterer%a = ((3*sum(Cells(1:Nbc)%Sc**3.))/(4*pi))**(1./3.)
        	SimScatterer%dm = 2.*SimScatterer%a; 
        endif 
        
        call date_and_time(date_final,time_final,zone_final,values_final)
        call Calcul_time_spent(values_init,values_final,Comp_time_disc)    
        
        ! Division into blocks     
        if ((CBFM .NE. 0) .OR. (MLCBFM .NE. 0)) Then 
            call date_and_time(date_init,time_init,zone_init,values_init); 
            ! Division into blocks depending on the type of scatterer  
            call Division_blocks(SimScatterer,Cells,Ncells_SphDomains,CBFM_Blocks,MLCBFM_BlDistr)    
            call date_and_time(date_final,time_final,zone_final,values_final)
            call Calcul_time_spent(values_init,values_final,Comp_time_div)     
                            
            ! EXTENTION OF THE CBFM BLOCKS TO AVOID THE IMPACT OF EDGE EFFECT ON THE ACCURACY OF THE SOLUTION
            ! This step is the same for all types of shapes
            !Comp_time_ext = 0; values_init=0; values_final =0;
            call date_and_time(date_init,time_init,zone_init,values_init); 
            Call Extend_blocks(SimScatterer,Cells,CBFM_Blocks,CBFM_Blocks_Ext)
            call date_and_time(date_final,time_final,zone_final,values_final)
            call Calcul_time_spent(values_init,values_final,Comp_time_ext)   
            
            !call Write_geometry_files(SimScatterer,Cells,CBFM_Blocks,CBFM_Blocks_Ext,'NEW');   !! TO CHECK AND MPI OPTIMIZE FROM MPI Code 1/29/2019
            
            ! I will discard this step for the moment for my MPI code !!
            ! Now get the Number of blocks for which the CBFs will be calculated and save the numbers 
            ! of the 'copy blocks' for each calculated block.
            !call Set_cp_CBFM_Blocks(Cells,CBFM_Blocks,CBFM_Blocks_Ext,cp_CBFM_Blocks);
            
            ! Check that all the jobs have the same (Division/Extension) Configuration
            allocate(all_NBlocks(nber_procs));
            call MPI_ALLGATHER(Nblocks,1,MPI_INTEGER,all_NBlocks,1,MPI_INTEGER,MPI_COMM_WORLD,code)
            call count_unique_vals(nber_procs,all_NBlocks,countU)
            if (countU .ne. 1) Then 
                if (rank == 0) Then 
                    Write(*,'(a)') 'Error : all the jobs have not the same division into blocks !';
                endif
                go to 30;            
            EndIf
            allocate(all_NbcBlocks(Nblocks,nber_procs));
            allocate(Nbc_blocks(Nblocks));
            Nbc_blocks = CBFM_Blocks(1:Nblocks)%Nbc_b;
            call MPI_ALLGATHER(Nbc_blocks,Nblocks,MPI_INTEGER,all_NbcBlocks,Nblocks,MPI_INTEGER,MPI_COMM_WORLD,code)
            deallocate(Nbc_blocks);
        
            Do ii=1, Nblocks
                allocate(Nbc_blocks(nber_procs));
                Nbc_blocks = all_NbcBlocks(ii,1:nber_procs);
                call count_unique_vals(nber_procs,Nbc_blocks,countU)
                if (countU .ne. 1) Then 
                    if (rank == 0) Then 
                        Write(*,'(a)') 'Error : all the jobs have not the same division into blocks !';
                    endif
                    go to 30;            
                EndIf 
                deallocate(Nbc_blocks);  
            EndDo
            deallocate(all_NBlocks,all_NbcBlocks);
            
            ! Blocks Distribution :  
            !! Divide up the CBFM blocks among the available MPI jobs
            Call MPI_distribution_blocks(CBFM_Blocks,MPI_CBFM_Blocks)
            if ((Adapt_mesh .eq. 0) .and. (nber_procs .gt. Nblocks)) then 
              if (rank == 0) then 
                Write(*,'(a,i4,a,i4,a)') 'Performance Error : Nprocs =',nber_procs,' > Nblocks =',Nblocks,' ! Please restart with fewer processors !';
              endif
              go to 30;
            endif   
      
            !if (Nbc .lt. 4e5) then 
            	call date_and_time(date_init,time_init,zone_init,values_init);             
            	call Write_geometry_files(SimScatterer,Cells,CBFM_Blocks,CBFM_Blocks_Ext,'NEW');   !! TO CHECK AND MPI OPTIMIZE  1/29/2019
            	call date_and_time(date_final,time_final,zone_final,values_final)
            	call Calcul_time_spent(values_init,values_final,Comp_time_write)     
	        !endif   
            
        Else 
            ! to avoid segmentation fault errors at the input of Compute_Electric_Fields.
            NBlocks = 1; Allocate(CBFM_Blocks(NBlocks));
            Nbc_ext =1; Allocate(CBFM_Blocks_Ext(NBlocks,Nbc_ext));
            Nblk_proc_max =1; Allocate(MPI_CBFM_Blocks(nber_procs,Nblk_proc_max+1));            
        EndIf     
        !! END DISCRETIZATION & DIVISION INTO BLOCKS***********************************************************************
        
        if (rank == 0) Then        
            Write(10,'(a)') ''
            Write(10,'(a)') 'Parameters of the scatterer : '
            Write(10,'(a,a)') 'P :    Type    a(mm)    dX(mm)   dY(mm)   dZ(mm)     X(m)      Y(m)      Z(m)         ',&
            '  m           lam(mm)   Sc(mm)' 
            p = SimScatterer%a*10**lamb_mag
            xp = SimScatterer%dx*10**lamb_mag
            yp = SimScatterer%dy*10**lamb_mag
            zp = SimScatterer%dz*10**lamb_mag
            ! type and effective radius
            Write(10,'(a,i1,a6,f10.4,f9.4,f9.4,f9.4)',advance='no') 'P : ',SimScatterer%type_s,trim(SimScatterer%info_s),p,xp,yp,zp        
            ! Position
            Write(10,'(e10.2)',advance='no') SimScatterer%xmin
            Write(10,'(e10.2)',advance='no') SimScatterer%ymin
            Write(10,'(e10.2)',advance='no') SimScatterer%zmin
            !Refractive index
            Write(10,'(a,f7.4,a,f7.4,a)',advance='no') '  (',real(SimScatterer%m_max),',',&
                imag(SimScatterer%m_max),')  '
            !Wavelength inside scatterer
            p = SimScatterer%lambda_min*10**lamb_mag
            Write(10,'(f8.3)',advance='no') p
            !Size of cell per scatterer (m/mm/um)
            Sc = SimScatterer%Sc*10**lamb_mag
            Write(10,'(f8.3)',advance='no') Sc        
    
            Write(10,*) ''
            Write(10,*) ''
            write (*,'(a)',advance='no') 'The dimensions of the scatterer = '
            write (*,'(f9.4,a)',advance='no') xp,'; ' 
            write (*,'(f9.4,a)',advance='no') yp,'; ' 
            write (*,'(f9.4,a,a)') zp, ' ',lamb_unit   
            write (*,'(a)',advance='no') 'The effective radius of the scatterer = '
            Write (*,'(f12.6,a,a)') SimScatterer%a*10**lamb_mag, ' ',lamb_unit 
            if (Adapt_mesh .eq. 0) then 
                Write(10,'(a,i7)') 'The Total Number of Cells =', Nbc
                Write(*,'(a,i7)') 'Total Number of Cells =', Nbc
            else
                Write(10,'(a,i7)') 'The Initial Total Number of Cells =', Nbc
                Write(*,'(a,i7)') 'Initial Total Number of Cells =', Nbc
            EndIf       
    
            Write(10,*) ''
            Write(10,*) ''
            Write(10,'(a,i2)') 'Nber_Methods = ',Nber_Methods
            Write(10,'(a)',advance='no') 'The applied methods are : '
            Do I=1,Nber_Methods-1
                Write(10,'(a,a)',advance='no') methods_names(I),'; '    
            Enddo
            Write(10,'(a)') methods_names(Nber_Methods)
    
            Write(*,*) '';Write(*,*) ''
            write (*,'(a)') '--------Numerical Methods-------- ' 
            Write(*,'(a,i2)') 'Nber_Methods = ',Nber_Methods
            Write(*,'(a)',advance='no') 'The applied methods are : '
            Do I=1,Nber_Methods-1
                Write(*,'(a,a)',advance='no') methods_names(I),'; '    
            Enddo
            Write(*,'(a)') methods_names(Nber_Methods)
    
    
            Write(*,*) ''; Write(*,*) ''
            Write(*,'(a)')  '------Transmitters/Receivers-----'
            Write(*,'(a,a,a,a)') 'Config of Tx/Rx = ', NumIntType_t,'/',NumIntType_r
            if ((NumIntType_t .eq. 'sd') .OR. (NumIntType_r .eq. 'sd')) then
                if (sd_type .eq. 1) Then
                    Write(*,'(a)') '-Type : FSU H&S St-d'                    
                ElseIf (sd_type .eq. 2) Then    
                    Write(*,'(a)') '-Type : W. St-d'             
                Else 
                    Write(*,'(a)') '-Type : W. Symm St-d'
                EndIf 
            endif
            Write(*,'(a,i5)') 'Nber of Tx = ', NTr
            if ((NumIntType_t .eq. 'aq') .OR. (NumIntType_t .eq. 'gl')) then 
                Write(*,'(a,i3,a,i3)') 'Ntheta = ',NTrTheta,'; Nphi = ',NTrPhi
            endif
            Write(*,'(a,i6)') 'Nber of Rx = ', NRx
            if ((NumIntType_r .eq. 'aq') .OR. (NumIntType_r .eq. 'gl')) then 
                Write(*,'(a,i3,a,i3)') 'Ntheta = ',NRxTheta,'; Nphi = ',NRxPhi
            endif      
        
            !*****************************************************************************************************     
            ! Write in the output file the transmitters/receivers
            Write(10,*) '';
            Write(10,'(a,a,a,a)') 'Config of Tx/Rx = ', NumIntType_t,'/',NumIntType_r
            if (NumIntType_t .eq. 'aq') then
                st_th_Tx=':  1.00';st_ph_Tx =':  1.00';st_th_Rx=':  1.00';st_ph_Rx =':  1.00';
                if (NTrTheta .gt. 1) then 
                    write(st_th_Tx,'(a,f6.2)') ':', (theta_final_trans_comp-theta_init_trans_comp)/(NTrTheta-1);
                endif            
                if (NTrPhi .gt. 1) then     
                    write(st_ph_Tx,'(a,f6.2)') ':', (phi_final_trans_comp-phi_init_trans_comp)/(NTrPhi-1);
                endif
                if (NRxTheta .gt. 1) then             
                    write(st_th_Rx,'(a,f6.2)') ':', (theta_final_Recei-theta_init_Recei)/(NRxTheta-1); 
                endif
                if (NRxPhi .gt. 1) then 
                    write(st_ph_Rx,'(a,f6.2)') ':', (phi_final_Recei-phi_init_Recei)/(NRxPhi-1);
                endif            
            else
                st_th_Tx = ''; st_th_Rx=''; st_ph_Tx = ''; st_ph_Rx='';
            
            endif
            Write(10,'(a,i5,a,f5.2,a,a,f6.2,a,f6.2,a,a,f6.2)') 'Nber_Transmitters =',NTr,'; Theta = ',theta_init_trans_comp,trim(st_th_Tx),':',&
                theta_final_trans_comp,'; Phi =',phi_init_trans_comp,trim(st_ph_Tx),':',phi_final_trans_comp
            Write(10,'(a,i5,a,f5.2,a,a,f6.2,a,f6.2,a,a,f6.2)') 'Nber_Receivers =',NRx,'; Theta = ',theta_init_Recei,trim(st_th_Rx),&
                ':',theta_final_Recei,'; Phi =',phi_init_Recei,trim(st_ph_Rx),':',phi_final_Recei
        
            ! Create and edit IncScattDirs file 
            Open(unit=41,File = trim(SimOutfld_name)//Env_sep//'IncScattDirs.dat');        
            Write(41,'(a)') 'INCIDENT DIRECTIONS : '
            Write(41,'(a,i5)') 'Ninc = ', NTr 
            if (NumIntType_t .ne. 'sd') then 
                Write(41,'(a,a)') 'Dist Type = ', NumIntType_t
            else
                if (sd_type .eq. 1) Then
                    Write(41,'(a)') 'Dist Type = sd (FSU H&S St-d)'                    
                ElseIf (sd_type .eq. 2) Then    
                    Write(41,'(a)') 'Dist Type = sd (W. St-d)'             
                Else 
                    Write(41,'(a)') 'Dist Type = sd (W. Symm St-d)'
                EndIf
            endif
            if ((NumIntType_t .eq. 'aq') .OR. (NumIntType_t .eq. 'gl')) Then 
                Write(41,'(a,i3,a,i3)') 'Ntheta = ',NTrTheta,'; Nphi = ',NTrPhi        
            EndIf   
            Write(41,'(a)') '    theta      phi ';    
            Do I=1,NTr
                Write(41,'(f9.2,f9.2)') Transmitters_Comp(I)%Theta, Transmitters_Comp(I)%Phi
            EndDo    
            Write(41,'(a)') ''
            Write(41,'(a)') 'SCATTERING DIRECTIONS : '
            Write(41,'(a,i5)') 'Nscat = ', NRx 
            if (NumIntType_r .ne. 'sd') then 
                Write(41,'(a,a)') 'Dist Type = ', NumIntType_r
            else
                if (sd_type .eq. 1) Then
                    Write(41,'(a)') 'Dist Type = sd (FSU H&S St-d)'                    
                ElseIf (sd_type .eq. 2) Then    
                    Write(41,'(a)') 'Dist Type = sd (W. St-d)'             
                Else 
                    Write(41,'(a)') 'Dist Type = sd (W. Symm St-d)'
                EndIf
            endif
            if ((NumIntType_r .eq. 'aq') .OR. (NumIntType_r .eq. 'gl')) Then 
                Write(41,'(a,i3,a,i3)') 'Ntheta = ',NRxTheta,'; Nphi = ',NRxPhi         
            EndIf   
            Write(41,'(a)') '    theta      phi ';    
            Do I=1,NRx
                Write(41,'(f9.2,f9.2)') Receivers(I)%Theta, Receivers(I)%Phi
            EndDo    
            Close(41);
    
            ! create the S_files folder if needed
            Sfold_name = trim(SimOutfld_name)//Env_sep//'S_files';
            inquire(directory=trim(Sfold_name),exist=dirExists);
            if (dirExists) Then
                call system('rm -r "'//trim(Sfold_name)//'"')        
            EndIf
            If (wr_Sij .eq. 1) then 
                call system('mkdir "'//trim(Sfold_name)//'"')
            EndIf  
            
            ! create the Q_files folder if needed
            if (EqSph==0) then
                Qfold_name = trim(SimOutfld_name)//Env_sep//'Q_files';
                inquire(directory=trim(Qfold_name),exist=dirExists);
                if (dirExists) Then
                    call system('rm -r "'//trim(Qfold_name)//'"')        
                EndIf
                If (wr_Qij .eq. 1) then 
                    call system('mkdir "'//trim(Qfold_name)//'"')
                EndIf  
            endif
        
            ! Create the Solution files (Eint, Zc, CBFs, Zcinv) if required
            If ((save_Eint .eq. 1) .OR. (save_Zc .eq. 1)) then 
                Solfold_name = trim(SimOutfld_name)//Env_sep//'Sol_files';
                inquire(directory=trim(Solfold_name),exist=dirExists);
                if (dirExists) Then
                    call system('rm -r "'//trim(Solfold_name)//'"')        
                EndIf
                call system('mkdir "'//trim(Solfold_name)//'"')            
            Endif
            if (Adapt_mesh .eq. 0) then
                if ((CBFM .NE. 0) .OR. (MLCBFM .NE. 0)) Then 
                    Write(10,*) ''; Write(10,*) ''
                    Write(10,'(a)') 'Division into blocks to apply the CBFM : '
                    Write(10,'(a,i6)') 'Total number of blocks = ',Nblocks        
                    Write(10,'(a)') 'Number of cells per Block  = '
                    Do I=1,Nblocks-1
	                    Write(10,'(i7,a)',advance='no') CBFM_Blocks(I)%Nbc_b,';'
                    Enddo
                    Write(10,'(i7)') CBFM_Blocks(Nblocks)%Nbc_b
                    Write(10,'(a)') 'Number of cells per Block after extension = '
                    Do I=1,Nblocks-1
	                    Write(10,'(i7,a)',advance='no') (CBFM_Blocks(I)%Nbc_b+CBFM_Blocks(I)%Nbc_ext),';'
                    Enddo
                    Write(10,'(i7)') (CBFM_Blocks(Nblocks)%Nbc_b+CBFM_Blocks(Nblocks)%Nbc_ext)
        
                    Write(*,*) ''; Write(*,*) ''
                    Write(*,'(a)')  '-------Division into blocks------'
                    Write(*,'(a,i6)') 'Total number of blocks = ',Nblocks
                    Write(*,'(a,f9.4,a,a)') 'Maximum block height = ', hBlock*10**lamb_mag, ' ',lamb_unit
                    Write(*,'(a,i6)') 'Maximum block length = ', maxval(CBFM_Blocks(1:NBlocks)%Nbc_b)
                    Write(*,'(a,i2)') 'Nbcells_ext =', Nc_extended
                    Write(*,'(a,i6)') 'Maximum block extension length = ', maxval(CBFM_Blocks(1:NBlocks)%Nbc_ext)    
                    if (define_use_Copies == 1) then 
                        Write(*,'(a,i4,a,i4,a)') 'Note that NcalBlks = ',NcalBlks,' out of total ',NBlocks,' blocks';
                    endif
                
                    Write(*,*) '';
                    Write (*,'(a,i2,a,i2,a,i2,a,i2,a)') '--> to discretize : ',Comp_time_disc(1),'j',Comp_time_disc(2)&
                    ,'h',Comp_time_disc(3),'min', Comp_time_disc(4),'sec'
                    Write (*,'(a,i2,a,i2,a,i2,a,i2,a)') '--> to divide into blocks: ',Comp_time_div(1),'j',Comp_time_div(2)&
                    ,'h',Comp_time_div(3),'min', Comp_time_div(4),'sec'
                    Write (*,'(a,i2,a,i2,a,i2,a,i2,a)') '--> to extend blocks: ',Comp_time_ext(1),'j',Comp_time_ext(2)&
                    ,'h',Comp_time_ext(3),'min', Comp_time_ext(4),'sec'
		            if (Nbc .lt. 100000) then
                     Write (*,'(a,i2,a,i2,a,i2,a,i2,a)') '--> to write geometry/blocks files : ',Comp_time_write(1),'j',Comp_time_write(2)&
                    ,'h',Comp_time_write(3),'min', Comp_time_write(4),'sec'
		            endif
                EndIf  
            EndIf
        EndIf
        ! here rank = 1 will quickly check if the folders are properly created  (I was having a weired problem of rank 0 not creating the folders !!)
        ! we will let all other tasks check if the folders ar created
        Do I = 1,4
          Call MPI_Barrier(MPI_COMM_WORLD,code); 
          if (rank == I) then
              ! create the S_files folder if needed
              Sfold_name = trim(SimOutfld_name)//Env_sep//'S_files';
              inquire(directory=trim(Sfold_name),exist=dirExists);
              if ((.not. dirExists) .and. (wr_Sij .eq. 1))  Then
                  call system('mkdir "'//trim(Sfold_name)//'"')    
              EndIf
              ! create the Q_files folder if needed
              if (EqSph==0) then
                  Qfold_name = trim(SimOutfld_name)//Env_sep//'Q_files';
                  inquire(directory=trim(Qfold_name),exist=dirExists);
                  if ((.not. dirExists) .and. (wr_Qij .eq. 1)) Then
                      call system('mkdir "'//trim(Qfold_name)//'"')
                  EndIf  
              endif
          
              ! Create the Solution files (Eint, Zc, CBFs, Zcinv) if required
              If ((save_Eint .eq. 1) .OR. (save_Zc .eq. 1)) then 
                  Solfold_name = trim(SimOutfld_name)//Env_sep//'Sol_files';
                  inquire(directory=trim(Solfold_name),exist=dirExists);
                  if (.not. dirExists) Then
                      call system('mkdir "'//trim(Solfold_name)//'"')      
                  EndIf                              
              Endif       
          endif  
        EndDo          
        Call MPI_Barrier(MPI_COMM_WORLD,code);          
        
        
        !! START THE COMPUTING OF THE ELECTRIC FIELDS DEPENDING ON THE FREQUENCY     
        Do ii=1, Nfreq   
            if (rank == 0) then 
                Write(10,*) ''; Write(10,*) ''
                Write(*,*) '';  Write(*,*) '';
            endif
            !! here we define the wavelength of the current experience  
            num_freq = ii;
            Lambda_w = Wavesle(ii)/10**lamb_mag
            Freq_w = C0/Lambda_w;
            Omega_w = 2*Pi*Freq_w         !! angular frequency
            !k_0 = Omega_w*sqrt(Eps0*Rmu0*Eps_air)
            k_0 = (2*Pi)/Lambda_w 
            ! refractive index
            SimScatterer%m_max = m_lambdas(1,ii)  
            mrp = real(m_lambdas(1,ii))
            mip = imag(m_lambdas(1,ii)) ; 
            rp = mrp**2-mip**2;
            ip = 2*mrp*mip;
            SimScatterer%Eps_max = rp+J*ip  
            SimScatterer%lambda_min = Lambda_w/sqrt(rp)
            ! Update scatterer Dlam
            SimScatterer%Dlamb = SimScatterer%lambda_min/SimScatterer%Sc;
            xmax = Pi*max(SimScatterer%dx,SimScatterer%dy,SimScatterer%dz)/Lambda_w; 
            xeq = 2*Pi*SimScatterer%a/Lambda_w; 
            
            Call DielComposition(m_lambdas,Cells); !!!! ATTENTION : TRAVAIL INACHEVE INPUT Mice AND Mwater + CALCUL EM En fction de Cells%m et pas Scatterer%m
                    
            ! here Set the cells parameter that depend on k_0 (so on lambda_p) and Eps_p
            old_Nbc = Nbc;
            Call SetCellsParams(SimScatterer,Cells,Upd_Cells);
            deallocate(Cells); Allocate(Cells(Nbc)); Cells=Upd_Cells; Deallocate(Upd_Cells);            
            
            if (rank == 0) then 
                Write(*,'(a,i3,a,i3,a)') '- SIMULATION ',ii,'/',Nfreq, ' : ****************************************'            
                a = nint(Freq_w/10.**freq_mag);
                if (a < 10) Then 
                    Allocate(character(5) ::stFreq); ty = '(f5.3)';
                ElseIf (a < 100) Then
                    Allocate(character(6) ::stFreq); ty = '(f6.3)';
                Else
                    Allocate(character(7) ::stFreq); ty = '(f7.3)';
                EndIf     
                Write(stFreq,ty) Freq_w/10.**freq_mag
                write (*,'(a,a,a,a)') 'The frequency of simulation = ',stFreq,' ',freq_unit 
                !write (*,'(a,F9.6,a,a)') ' -- > Wavelength = ',Lambda_w*10**lamb_mag,' ',lamb_unit
                write (*,*) ' -- > Wavelength = ',Lambda_w*10**lamb_mag,' ',lamb_unit    
                if (homogs .eq. 1) then  
                    !write (*,'(a,F9.6,a,a)') ' -- > Wavelength inside scatterer = ',SimScatterer%lambda_min*10**lamb_mag,' ',lamb_unit 
                    write (*,*) ' -- > Wavelength inside scatterer = ',SimScatterer%lambda_min*10**lamb_mag,' ',lamb_unit 
                    write (*,'(a,F7.4,a,ES10.3)') ' -- > m = ',mrp,' + j*',mip  
                    write (*,'(a,F7.4,a,ES10.3)') ' -- > Eps = ',rp,' + j*',ip  
                    deallocate(stFreq);            
                    Write(*,'(a,i4)') ' -- > Dlambda = ', SimScatterer%Dlamb
                    write (*,'(a,f6.4)') ' -- > d/aeff = ', SimScatterer%Sc/SimScatterer%a
                    write (*,'(a,f6.4)') ' -- > kd = ', k_0*SimScatterer%Sc
                    write (*,'(a,f6.4)') ' -- > |m|kd = ', abs(SimScatterer%m_min)*k_0*SimScatterer%Sc
                    
                    xeq_m = 2*Pi*SimScatterer%a/SimScatterer%lambda_min;
                    xmax_m = Pi*max(SimScatterer%dx,SimScatterer%dy,SimScatterer%dz)/SimScatterer%lambda_min;
                    write (*,'(a,f6.2)') ' -- > xeq =', xeq                
                    write (*,'(a,f6.2)') ' -- > xmax =', xmax 
                    write (*,'(a,f6.2)') ' -- > xeq_m =', xeq_m
                    write (*,'(a,f6.2)') ' -- > xmax_m =', xmax_m
                else
                    Allocate(vals(Nbc)); vals = Cells(1:Nbc)%lambda_n;
                    r_min = minval(vals); r_max = maxval(vals);
                    write (*,'(a,F9.6,a,F9.6,a,a)') ' -- > Wavelength inside scatterer = [',r_min*10**lamb_mag,' - ',r_max*10**lamb_mag,'] ',lamb_unit;
                    vals = real(Cells(1:Nbc)%m_n); r_min = minval(vals); r_max = maxval(vals);
                    vals = imag(Cells(1:Nbc)%m_n); i_min = minval(vals); i_max = maxval(vals);
                    write (*,'(a,F7.4,a,ES10.3,a,F7.4,a,ES10.3,a)') ' -- > m = [',r_min,' + j*',i_min,' - ',r_max,' + j*',i_max,']';
                    vals = real(Cells(1:Nbc)%Eps_n); r_min = minval(vals); r_max = maxval(vals); 
                    vals = imag(Cells(1:Nbc)%Eps_n); i_min = minval(vals); i_max = maxval(vals);
                    write (*,'(a,F7.4,a,ES10.3,a,F7.4,a,ES10.3,a)') ' -- > Eps = [',r_min,' + j*',i_min,' - ',r_max,' + j*',i_max,']'   
                    deallocate(stFreq); 
                    vals = Cells(1:Nbc)%Dlamb_n; r_min = minval(vals); r_max = maxval(vals); 
                    Write(*,'(a,f7.2,a,f7.2,a)') ' -- > Dlambda = [', r_min,' - ', r_max,']';
                    write (*,'(a,f6.4)') ' -- > d/aeff = ', SimScatterer%Sc/SimScatterer%a
                    write (*,'(a,f6.4)') ' -- > kd = ', k_0*SimScatterer%Sc
                    vals = abs(Cells(1:Nbc)%m_n); r_min = minval(vals)*k_0*SimScatterer%Sc; r_max = maxval(vals)*k_0*SimScatterer%Sc
                    write (*,'(a,f6.4,a,f6.4,a)') ' -- > |m|kd = [', r_min,' - ', r_max,']'
                    write (*,'(a,f6.2)') ' -- > xeq =', xeq                
                    write (*,'(a,f6.2)') ' -- > xmax =', xmax 
                    vals = Cells(1:Nbc)%lambda_n;
                    r_max = 2*Pi*SimScatterer%a/minval(vals); 
                    r_min = 2*Pi*SimScatterer%a/maxval(vals);
                    write (*,'(a,f6.2,a,f6.2,a)') ' -- > xeq_m = [',r_min,' - ',r_max,']';
                    r_max = Pi*max(SimScatterer%dx,SimScatterer%dy,SimScatterer%dz)/minval(vals); 
                    r_min = Pi*max(SimScatterer%dx,SimScatterer%dy,SimScatterer%dz)/maxval(vals);;
                    write (*,'(a,f6.2,a,f6.2,a)') ' -- > xmax_m = [',r_min,' - ',r_max,']';
                    deallocate(vals);
                endif
                                                   
                !*****************************************************************************************************
                !! Once generated, all these informations should be written in the output file  
                Write(10,*) ''
                Write(10,*) ''
                Write(10,'(a,i3,a,i3,a)') 'Simulation ',ii,'/',Nfreq, ' : ****************************************'             
                Write(10,*) ''
                Write(10,'(a,es12.2)') 'The frequency of simulation = ',Freq_w
                Write(10,'(a,F10.6,a,a)') 'The wavelength of simulation = ',Lambda_w*10**lamb_mag,' ',lamb_unit 
            endif
            
            if ((Adapt_mesh .eq. 1) .AND. (Nbc .ne. old_Nbc)) then 
            
                if (rank == 0) then 
                    Write(10,*) '';
                    Write(10,'(a,i7)') 'The updated Number of Cells After Adaptive Mesh =', Nbc
                    Write(*,*) ''; 
                    Write(*,'(a,i8)') 'After Adaptive Mesh, Total Number of Cells =',Nbc
                EndIf 
            
                if ((CBFM .NE. 0) .OR. (MLCBFM .NE. 0)) Then 
                    call Division_blocks(SimScatterer,Cells,Ncells_SphDomains,CBFM_Blocks,MLCBFM_BlDistr)    
                    Call Extend_blocks(SimScatterer,Cells,CBFM_Blocks,CBFM_Blocks_Ext)
                    call Write_geometry_files(SimScatterer,Cells,CBFM_Blocks,CBFM_Blocks_Ext,'UPD');   
                    
                    deallocate(MPI_CBFM_Blocks);                
                    Call MPI_distribution_blocks(CBFM_Blocks,MPI_CBFM_Blocks);
                    
                    if (rank == 0) then                    
                      Write(10,'(a)') 'Division into blocks to apply the CBFM : '
                      Write(10,'(a,i6)') 'Total number of blocks = ',Nblocks        
                      Write(10,'(a)') 'Number of cells per Block  = '
                      Do I=1,Nblocks-1
  	                    Write(10,'(i7,a)',advance='no') CBFM_Blocks(I)%Nbc_b,';'
                      Enddo
                      Write(10,'(i7)') CBFM_Blocks(Nblocks)%Nbc_b
                      Write(10,'(a)') 'Number of cells per Block after extension = '
                      Do I=1,Nblocks-1
  	                    Write(10,'(i7,a)',advance='no') (CBFM_Blocks(I)%Nbc_b+CBFM_Blocks(I)%Nbc_ext),';'
                      Enddo
                      Write(10,'(i7)') (CBFM_Blocks(Nblocks)%Nbc_b+CBFM_Blocks(Nblocks)%Nbc_ext)
          
                      Write(*,'(a,i6)') ' -- > Total number of blocks = ',Nblocks
                      Write(*,'(a,f6.3,a)') ' -- > Maximum block height = ', hBlock*1e3, ' mm'
                      Write(*,'(a,i6)') ' -- > Maximum block length = ', maxval(CBFM_Blocks(1:NBlocks)%Nbc_b)   
                      Write(*,'(a)') ''; Write(*,'(a)') ''
                    endif 
                EndIf                 
            EndIf
            
            ! check if Nprocs >= Nblocks
            if ((CBFM .NE. 0) .OR. (MLCBFM .NE. 0)) then 
                if (nber_procs .gt. Nblocks) then 
                    if (rank == 0) then 
                      Write(*,'(a,i4,a,i4,a)') 'Performance Error : Nprocs =',nber_procs,' > Nblocks =',Nblocks,' ! Please restart with fewer processors !';
                    endif
                    go to 30;   
                else
                    Call initializeNipws(SimScatterer);                                                                           
                Endif                                                                       
            Endif        
                       
            ! START COMPUTING OF THE ELECTRIC FIELDS
            Call Compute_Electric_Fields(SimScatterer,Cells,Transmitters_Comp,Receivers,methods_names,&
                CBFM_Blocks,CBFM_Blocks_Ext,MPI_CBFM_Blocks);      
        EndDo
        deallocate(m_lambdas)
    EndDo

    if (rank == 0) then 
        Write (*,'(a)') '*******************************************************************************'
        Write (*,'(a)') '************************* END OF SIMULATION, THANKS ***************************'
        Write (*,'(a)') '*******************************************************************************' 

        Write(10,'(a)') ''
        Write(10,'(a)') ''
        Write(10,'(a)') 'End of the simulation'; Write(10,'(a)') ''
        Write (10,'(a)') '*******************************************************************************'
        Write (10,'(a)') '*******************************************************************************'
        Close(10)
    endif
    Deallocate(Transmitters_Comp,Receivers,Cells)   
    Deallocate (methods_names)
    if ((CBFM .NE. 0) .OR. (MLCBFM .NE. 0)) Then 
        Deallocate(CBFM_Blocks); !! attention test 8/8/2018 comment/uncomment depending on test or no
    EndIf    
    
30  Call MPI_FINALIZE (code);
        
End PROGRAM Main_Scattering


