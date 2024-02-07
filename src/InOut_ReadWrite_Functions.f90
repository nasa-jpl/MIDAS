!! SUBROUTINES : 
!! - Get_InputData
!! - Read_ShapeFile
!! - Read_CellsFile
!! - Write_geometry_files
!! - Write_jobs_sim_info
!! - Write_memory_info : system_mem_usage & print_allocate

!! READ IMPUT FILE --------------------------------------------------------
SUBROUTINE Get_InputData(SimScatterer,Wavesle,methods_names,m_file_name,Transmitters_Comp,Receivers,error)
    USE Initialization
    USE common_variables
    USE iso_fortran_env
    USE DiverseUtil
    USE MPI
    implicit none
    
    ! IN/OUT
    type (Scatterer), INTENT(INOUT) :: SimScatterer
    Real(kind=8), Dimension(:), allocatable, INTENT(OUT) :: Wavesle
    character(:) ,allocatable, INTENT(OUT) :: methods_names(:)
    character(250),allocatable, INTENT(OUT) :: m_file_name(:)
    type(Dipole), Dimension(:),allocatable, INTENT(OUT) :: Transmitters_Comp
    type(Dipole), Dimension(:),allocatable, INTENT(OUT) :: Receivers
    Integer, INTENT(OUT) :: error
    
    ! Local 
    Integer :: Nwave, ind_line, ii,Ninc_sugg, Nscat_sugg
    Integer, Dimension(:), allocatable :: diel_comp_perc
    Real(kind=8) :: freq_min, freq_max, Wave_min, Wave_max, step_wave, step_freq
    Real(kind=8) :: p, r_cyl, l_cyl, ap, v_freq, Sc
    Character(1) :: ch_tmp
    character(250) :: m_file_name_0,inputline
    character(200) :: file_name
    character(9) :: wave_descr,SR_Zc_type_ch 
    character(6) :: tdata
    character(11) :: freq_unit_tmp
    character(10) :: lamb_unit_tmp
    
    INTERFACE 
        SUBROUTINE get_trans_Receiv(Ninc_in,Nscat_in,Transmitters_Comp,Receivers);
            USE Initialization
            USE common_variables
    
            IMPLICIT NONE
    
            !! IN/OUT ******************************************************************
    
            Integer, INTENT(IN) :: Ninc_in,Nscat_in
            type (Dipole), Dimension(:), allocatable, INTENT(OUT) :: Transmitters_Comp
            type (Dipole), Dimension(:), allocatable, INTENT(OUT) :: Receivers                   
        END SUBROUTINE get_trans_Receiv
    END INTERFACE
    
    !! Initialize error to 0
    error = 0
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
                error = 1
            endif
            error = 1
            go to 40
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
            error = 1
            go to 40
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
            error = 1
            go to 40
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
            error = 1
            go to 40
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
    read(11,'(a)'), Outfld_name
    
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
        if (rank == 0) then 
            Write(*,'(a,a,a)')'Error : ', dielcomp_option, ' is an unknown dielectric decomposition option !!'
        endif
        error = 1
        go to 40
    endif    
        
    if (((trim(dielcomp_option) == 'fromdielcompositionfile') .OR. (trim(dielcomp_option) == 'fromshapefile')) &
        .AND. (SimScatterer%type_s .ne. 2)) then
        if (rank == 0) then
            Write(*,'(a,a)') 'Error : The requested dielectric decomposition option can only be used with type_scatterer = 2';
        endif
        error = 1
        go to 40        
    Endif     
    read(11,*)
        
    !! Parameters of the dicretization 
    read(11,*)
    read(11,*), ch_tmp;
    if (ch_tmp .eq. 'S') then
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
        error = 1
        go to 40            
    endif
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
    read(11,*);read(11,*),NPolBeta
    beta_init_Pol= 0; beta_final_Pol=360.; ! we only consider 0-2Pi polar rotation 
    
    ! Write Scattering matrix elements for each incident direction and Q per incident direction 
    read(11,*);read(11,*), wr_Sij
    read(11,*);read(11,*), wr_Qij    
    read(11,*)
    
    ! Get Transmitters/Scatterers depending on the type of the numerical integration used to average the scattering quantities
    ! over incident/scattering directions
    call get_trans_Receiv(Ninc_sugg,Nscat_sugg,Transmitters_Comp,Receivers);
    ! if only 1 incident direction is used, we autmatically put wr_Sij and wr_Qij to 1
    if (NTr .eq. 1) then 
        wr_Sij=1;wr_Qij=1;
    endif
    
    
    ! Parameters of the numerical methods 
    
    !CBFM
    read(11,*)
    read(11,*);read(11,*), div_type
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
         
    !ACA !! 
    ! As we are not using the ACA for the MPI version yet, I deleted these lines 
    ! and simply initialized the ACA params to 0
    !read(11,*);
    !read(11,*);read(11,*), Use_ACA
    !read(11,*);read(11,*), Nb_it_max                                                                                                                                                                                     
    !read(11,*);read(11,*), Epsilon_ACA                                                                                                                                                                                           
    !read(11,*);read(11,*), Vrb_ACA
    !read(11,*)
    Use_ACA = 0; Nb_it_max= 50; Epsilon_ACA = 1E-4; Vrb_ACA = 0; 
    
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
        error = 1
        go to 40          
    EndIf
    
    ! get Far fiel approximation params
    ! in practice, FFA = 1 for precipitation particles & FFA = 0 for asteroid simulation 
    read(11,*);
    read(11,*);read(11,*), FFA
    read(11,*);read(11,*), Rso
    read(11,*);
    
    ! Save Sol Elements
    read(11,*);
    read(11,*);read(11,*), save_Zc
    read(11,*);read(11,*), save_Eint
    read(11,*), save_Eint_Nmax  ! used only if save_Eint=1
    read(11,*);read(11,*), save_Einc ! incident field (useful for inversion algorithms)        
    
    !! close the dat file
40  Close(11)
    
END SUBROUTINE Get_InputData

SUBROUTINE Read_ShapeFile(info_p,pr_NBcels,pr_lattice)

    USE Initialization
    USE common_variables
    
    Implicit NONE

    !IN/OUT 
    Character, INTENT(IN) :: info_p
    Integer, INTENT(OUT) :: pr_NBcels
    Integer, Dimension(:,:), allocatable, INTENT(OUT):: pr_lattice
    
    !! LOCAL
    Integer :: ii,jj
    Real(kind=8), Dimension(3) :: a1,a2
    Character ll
    Character(47) cc
        
    ! Read the shape file 
    Open(11,File = trim(ShapeFilePath))
    ! Inter-Dipole Distance
    !if ((info_p == 'a') .or. (info_p == 's')) Then
    !    read(11,'(a38,e12.8e3,a)') cc,Int_Dist,ll
    !ElseIf (info_p == 'p') Then
    !    read(11,'(a47,e12.8e3,a)') cc,Int_Dist,ll
    !EndIf    
    !Int_Dist = Int_Dist * 1e-6;
    read(11,*),cc
    
    ! NB_cells 
    read(11,*) pr_NBcels,ll
    
    ! a1 and a2
    read(11,'(f9.4,f9.4,f9.4,a)') a1(1),a1(2),a1(3),ll
    read(11,'(f9.4,f9.4,f9.4,a)') a2(1),a2(2),a2(3),ll
    
    ! For the moment we neglect the folowing 3 lines 
    read(11,*),cc;read(11,*),cc;read(11,*),cc;
    
    ! Read the positions in the lattice of the pr_NBcels cells (previously dipoles)
    Allocate(pr_lattice(pr_NBcels,6));
    Do ii= 1,pr_NBcels
        read(11,'(i7,i5,i5,i5,i5,i5,i5)') jj,pr_lattice(ii,1),pr_lattice(ii,2),pr_lattice(ii,3),&
            pr_lattice(ii,4),pr_lattice(ii,5),pr_lattice(ii,6) 
        ! To use only if you need to read a shape file copied from a 'target.out'file
        !read(11,'(i7,i5,i4,i4,i2,i2,i2)') jj,pr_lattice(ii,1),pr_lattice(ii,2),pr_lattice(ii,3),&
        !    pr_lattice(ii,4),pr_lattice(ii,5),pr_lattice(ii,6) 
    EndDo   

END SUBROUTINE Read_ShapeFile

SUBROUTINE Read_CellsFile(N,Cells_xyz_Sc,Cells_m_ind)

    USE Initialization
    USE common_variables
    
    Implicit NONE

    !IN/OUT 
    Integer, INTENT(OUT) :: N
    Real(kind=8), Dimension(:,:), allocatable, INTENT(OUT):: Cells_xyz_Sc
    Integer, Dimension(:), allocatable, INTENT(OUT):: Cells_m_ind
    
    !! LOCAL
    Integer :: ii,num_B
    character*1 :: cc
        
    ! Read the shape file 
    N = 0 
    OPEN(1,File = trim(ShapeFilePath))
    DO 
      READ (1,*, END=10) 
      N = N + 1 
    END DO 
    10 CLOSE (1) 
    IF (rank == 0) THEN
        write(*,*) 'N = ',N
    ENDIF 
    
    ! Read the positions in the lattice of the pr_NBcels cells (previously dipoles)
    Allocate(Cells_xyz_Sc(N,4));
    Allocate(Cells_m_ind(N));
    OPEN(1,File = trim(ShapeFilePath))
    Do ii= 1,N
        read(1,'(f12.6,a,f12.6,a,f12.6,a,f12.6,a,i8,a,i6)') Cells_xyz_Sc(ii,1),cc,Cells_xyz_Sc(ii,2),cc,Cells_xyz_Sc(ii,3),cc,&
            Cells_xyz_Sc(ii,4),cc,Cells_m_ind(ii),cc,num_B  
    EndDo   

END SUBROUTINE Read_CellsFile


SUBROUTINE Write_geometry_files(SimScatterer,Cells,CBFM_Blocks,CBFM_Blocks_Ext,option)

    !! 9-18-2019 : option added to choose the type of file writing depending on the value of adaptive mesh and the comparison between old_Nbc and Nbc
    USE Initialization
    USE common_variables
    USE MPI
    
    Implicit NONE
    
    !IN/OUT 
    type (Scatterer), INTENT(IN) :: SimScatterer
    type (Cell), Dimension(Nbc), INTENT(INOUT) :: Cells
    type(CBFM_Block), Dimension(Nblocks), INTENT(IN) :: CBFM_Blocks
    Integer, Dimension(Nblocks,Nbc_ext), INTENT(IN) :: CBFM_Blocks_Ext
    CHARACTER(3), INTENT(IN) :: option 

    Integer :: Type_sca,ii,jj,Nbcelsx,Nbcelsy,Nbcelsz
    Integer :: Nbcels_Dp,Nbcels_ext,Nbcels_int
    Real(kind=8) :: x0,y0,z0,Sc,ap
    
    Integer, Dimension(:,:), allocatable :: Part_in_lat
    Integer, Dimension(:), allocatable :: Nbcels_x,Nbcels_y,Nbcels_z
    
    CHARACTER(240) file_name
    CHARACTER(:), allocatable :: num_freq_str

    Type_sca = SimScatterer%type_s;
    Sc = SimScatterer%Sc; 
    ap = SimScatterer%dm/2.;
    
    ! Cells.dat file 
    if (rank == 0) then 
        ! Once Cells is reorganized, we can save it 
        ! The file Cellules.dat is used later to plot the 3D simulation scene 
        if (trim(option) .eq. 'NEW') then 
            If (EqSph ==0) then 
                file_name = trim(SimOutfld_name)//Env_sep//'Cells.dat';
            Else
                file_name = trim(SimOutfld_name)//Env_sep//'CellsES.dat';
            EndIf
        elseif (trim(option) .eq. 'UPD') then 
            if (Nfreq .eq. 1) then 
                file_name = trim(SimOutfld_name)//Env_sep//'Cells_adm.dat';
            else
                if (num_freq < 10) Then
                    Allocate(character(1)::num_freq_str)
                    Write(num_freq_str,'(i1)') num_freq
                ElseIf (num_freq < 100) Then
                    Allocate(character(2)::num_freq_str)
                    Write(num_freq_str,'(i2)') num_freq
                Else
                    Allocate(character(3)::num_freq_str)
                    Write(num_freq_str,'(i3)') num_freq
                EndIf 
                file_name = trim(SimOutfld_name)//Env_sep//'Cells_adm_f'//num_freq_str//'.dat';
            endif
        endif
    
        ! cell%num_diel to be written in Cells.dat file
        If (trim(dielcomp_option) == 'fromonlymfile') then 
            Cells(1:Nbc)%n_diel = 1;
        ElseIf ((trim(dielcomp_option) == 'fromdielcompositionfile') .OR. (trim(dielcomp_option) == 'random1')) then
            Cells(1:Nbc)%n_diel = (/1:Nbc/); ! for the other option ('fromshapefile') is read from shape file !         
        endif
        Open(14,File = trim(file_name))
        Do ii=1, Nbc
            !Write(14,'(f12.6,a,f12.6,a,f12.6,a,f12.6,a,i4,a,i4,a,f7.4,a,ES10.3,a,f7.4,a,ES10.3)') Cells(ii)%Xc,';',Cells(ii)%Yc, &
            !    ';',Cells(ii)%Zc,';',Cells(ii)%Sc,';',Cells(ii)%num_block,';',Cells(ii)%num_diel,';',&
            !    real(Cells(ii)%m_cell),' + j*',imag(Cells(ii)%m_cell),';', real(Cells(ii)%Eps_cell),' + j*',imag(Cells(ii)%Eps_cell);
        
            Write(14,'(f12.6,a,f12.6,a,f12.6,a,f12.6,a,i8,a,i6)') Cells(ii)%Xc,';',Cells(ii)%Yc, &
                ';',Cells(ii)%Zc,';',Cells(ii)%Sc,';',Cells(ii)%n_diel,';',Cells(ii)%n_block
        EndDo
        Close(14);
    endif
    
    if ((rank == 1) .AND. (trim(option) .eq. 'NEW')) then           
        If ((Type_sca .eq. 2)) Then    ! simply copy shape file 
            file_name=ShapeFilePath; 
            if (rank ==0 ) then 
                call system('cp '//trim(file_name)//' "'//trim(SimOutfld_name)//'/"');  
            endif

        else
            !! generate the file shape.dat for DDSCat simulations
            ! Of course here, for the moment,we suppose that we have only 1 scatterer
            ! we need the shape file for the validation of the MoM/CBFM results in comparison to DDScat and FEKO
            Nbcels_Dp = nint(SimScatterer%dm/Sc);
            Nbcels_ext = ceiling((ap - (0.8*ap/sqrt(3.)))/Sc) ; !(ap_cheb -c/2)/Sc + 0.8 is for "security" to be sure that 
                                                                ! this cube belongs entirely to the chebyshev particle
            Nbcels_int = Nbcels_Dp - 2*Nbcels_ext;          
            Nbcels_x = [Nbcels_int,Nbcels_Dp,Nbcels_ext,Nbcels_ext,Nbcels_int,Nbcels_int,Nbcels_Dp]; 
            Nbcels_y = [Nbcels_int,Nbcels_Dp,Nbcels_Dp,Nbcels_Dp,Nbcels_ext,Nbcels_ext,Nbcels_Dp];
            Nbcels_z = [Nbcels_int,Nbcels_ext,Nbcels_int,Nbcels_int,Nbcels_int,Nbcels_int,Nbcels_ext];
            Allocate(Part_in_lat(Nbc,6));
            If ((Type_sca == 1)) Then
                x0 = (1-0.5-Nbcels_Dp/2.)*Sc
                y0 = (1-0.5-Nbcels_Dp/2.)*Sc
                z0 = (1-0.5-Nbcels_Dp/2.)*Sc       
            ElseIf ((Type_sca == 6)) Then
                x0 = SimScatterer%xmin
                y0 = SimScatterer%ymin
                z0 = SimScatterer%zmin      
            ElseIf ((Type_sca == 3) .or. (Type_sca == 4))Then
                x0 = anint(((1-0.5-Nbcelsx/2.)*Sc)*10**Round_S)/10**Round_S
                y0 = anint(((1-0.5-Nbcelsy/2.)*Sc)*10**Round_S)/10**Round_S
                z0 = anint(((1-0.5)*Sc)*10**Round_S)/10**Round_S
            EndIf
            
            Do ii=1,Nbc
                Part_in_lat(ii,1) = nint((Cells(ii)%Xc - x0)/Sc) 
                Part_in_lat(ii,2) = nint((Cells(ii)%Yc - y0)/Sc)
                Part_in_lat(ii,3) = nint((Cells(ii)%Zc - z0)/Sc)             
            EndDo
            Part_in_lat(1:Nbc,4:6) = 1;
            !if (rank == 0) then    
            !    if (EqSph ==0) then 
            !        file_name = trim(SimOutfld_name)//Env_sep//'shape.dat';
            !    else
            !        file_name = trim(SimOutfld_name)//Env_sep//'shapeES.dat';
            !    endif        
            !    Open(14,File = trim(file_name))
            !    Open(12,File = 'inputs/Shape_head.dat')
            !    read(12,'(a)'), fline
            !    Write(14,'(a)') fline
            !    Write(14,'(i7,a)') Nbc,' # Number of Dipoles'
            !
            !    Do ii=1, 5
            !        read(12,'(a)'), fline;
            !        Write(14,'(a)') fline;            
            !    EndDo
            !    Close(12)        
            !    Do ii=1, Nbc
            !        Write(14,'(i7,i5,i5,i5,i5,i5,i5)') ii,Part_in_lat(ii,1), &
            !        Part_in_lat(ii,2),Part_in_lat(ii,3),Part_in_lat(ii,4), &
            !        Part_in_lat(ii,5),Part_in_lat(ii,6)
            !    EndDo
            !    Close(14);  
            !EndIf
        EndIf 
    endif
    
    ! Once Cells is reorganized, we can save it 
    ! The file Cellules.dat is used later to plot the 3D simulation scene 
    if (trim(option) .eq. 'NEW') then
        if (EqSph==0) Then 
            file_name = trim(SimOutfld_name)//Env_sep//'Blocks.dat';
        Else
            file_name = trim(SimOutfld_name)//Env_sep//'BlocksES.dat';
        EndIf   
    Elseif (trim(option) .eq. 'UPD') then 
        if (Nfreq .eq. 1) then 
            file_name = trim(SimOutfld_name)//Env_sep//'Blocks_adm.dat';
        else
            file_name = trim(SimOutfld_name)//Env_sep//'Blocks_adm_f'//num_freq_str//'.dat';
        endif
    Endif
    if (rank == 2) then 
        Open(14,File = trim(file_name))
        Do ii=1, Nblocks
            Write(14,'(a,i4)') '****** Block ',ii
            Write(14,'(a,i5)') 'Nbc = ',CBFM_Blocks(ii)%Nbc_b
            Write(14,'(a,i5)') 'Nbc_ext = ',CBFM_Blocks(ii)%Nbc_ext
            Write(14,'(a)') 'Num_cells_ext ='
            Do jj=1,CBFM_Blocks(ii)%Nbc_ext
                Write(14,'(i8)') CBFM_Blocks_Ext(ii,jj)
            EndDo
            Write(14,'(a)') ' ';
        EndDo
        Close(14)
    endif   
    
END SUBROUTINE Write_geometry_files


SUBROUTINE Write_Sfiles(nom_methode,Transmitters,Receivers,S_total)

    ! PHDF5 instructions copied from ph5example.f90 : https://github.com/mokus0/hdf5/blob/master/fortran/examples/ph5example.f90
    
    ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    !   Copyright by The HDF Group.                                               *
    !   Copyright by the Board of Trustees of the University of Illinois.         *
    !   All rights reserved.                                                      *
    !                                                                             *
    !   This file is part of HDF5.  The full HDF5 copyright notice, including     *
    !   terms governing use, modification, and redistribution, is contained in    *
    !   the files COPYING and Copyright.html.  COPYING can be found at the root   *
    !   of the source code distribution tree; Copyright.html can be found at the  *
    !   root level of an installed copy of the electronic HDF5 document set and   *
    !   is linked from the top-level documents page.  It can also be found at     *
    !   http://hdfgroup.org/HDF5/doc/Copyright.html.  If you do not have          *
    !   access to either file, you may request a copy from help@hdfgroup.org.     *
    ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    !
    ! Fortran parallel example.  Copied from Tutorial's example program of
    ! dataset.f90.

     !PROGRAM DATASET
     
    USE Initialization
    USE common_variables
    USE MPI
    !USE HDF5 ! This module contains all necessary modules
    
    IMPLICIT NONE

    !IN/OUT 
    character(8), INTENT(IN):: nom_methode
    type (Dipole), Dimension(NTr), INTENT(IN) :: Transmitters
    type (Dipole), Dimension(NRx_tot), INTENT(IN) :: Receivers
    COMPLEX(real64), Dimension(NRx_tot,4*NTr), INTENT(IN):: S_total


    ! Local 
    Integer :: a,kkt,kkr,NTr_wr_proc,Nths,Nphs
    Real(kind=8) :: th_i,ph_i
    CHARACTER(6) :: ty, kkt_st 
    CHARACTER(200) :: file_name_s, Sfold_name
    
    
    Real(kind=8), Dimension(:), allocatable :: Thetas,Phis,RecThetasVals,RecPhisVals
    Real(kind=8), Dimension(:,:), allocatable :: S_towrite
    CHARACTER(:), allocatable:: nom_meth_exact,stFreq,sim_name     
    COMPLEX(real64), Dimension(:,:), allocatable :: Svv_2d,Shh_2d,Svh_2d,Shv_2d
    COMPLEX(real64), Dimension(:), allocatable :: Svv_1d,Shh_1d,Svh_1d,Shv_1d
    
    ! PHDF5
    CHARACTER(LEN=10), PARAMETER :: default_fname = "sds.h5"  ! Default name
    CHARACTER(LEN=9), PARAMETER :: dsetname = "Smatrices" ! Dataset name

    CHARACTER(LEN=100) :: filename  ! File name
    INTEGER        :: fnamelen	     ! File name length
    !INTEGER(HID_T) :: file_id       ! File identifier
    !INTEGER(HID_T) :: dset_id       ! Dataset identifier
    !INTEGER(HID_T) :: filespace     ! Dataspace identifier in file
    !INTEGER(HID_T) :: plist_id      ! Property list identifier

!    INTEGER(HSIZE_T), DIMENSION(2) :: dimsf = (/5,8/) ! Dataset dimensions.
!     INTEGER, DIMENSION(7) :: dimsfi = (/5,8,0,0,0,0,0/)
!     INTEGER(HSIZE_T), DIMENSION(2) :: dimsfi = (/5,8/)
!    INTEGER(HSIZE_T), DIMENSION(2) :: dimsfi

    INTEGER, ALLOCATABLE :: data(:,:)   ! Data to write
    INTEGER :: data_rank = 2 ! Dataset rank

    INTEGER :: error, error_n  ! Error flags
    INTEGER :: ii, jj
    
    Sfold_name = trim(SimOutfld_name)//Env_sep//'S_files';
    
    NTr_wr_proc = (NTr/nber_procs)+1;
    
    If (nom_methode=='CBFM-E  ') Then
        Allocate(character(6) ::nom_meth_exact)
        nom_meth_exact = trim(nom_methode)
    Elseif ((nom_methode=='MoM     ') .OR. (nom_methode=='RGE     ')) Then 
        Allocate(character(3) ::nom_meth_exact)
        nom_meth_exact = trim(nom_methode)
    Endif
     
    a = nint(Freq_w/10**freq_mag);
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

    If (Nfreq == 1) Then
        Allocate(character(1)::sim_name)
        sim_name= ''
    Else
        if (num_freq < 10) Then
            Allocate(character(5)::sim_name)
            Write(sim_name,'(a,i1,a)') 'Sim', num_freq, '_'
        ElseIf (num_freq < 100) Then
            Allocate(character(6)::sim_name)
            Write(sim_name,'(a,i2,a)') 'Sim', num_freq, '_'
        Else
            Allocate(character(7)::sim_name)
            Write(sim_name,'(a,i3,a)') 'Sim', num_freq, '_'
        EndIf        
    EndIf 
    Write(stFreq,ty) Freq_w/10**freq_mag
    
    if (EqSph == 0) then
        file_name_s = trim(Sfold_name)//Env_sep//sim_name//'Smtable_'//stFreq//trim(freq_unit)//'_'//nom_meth_exact//'.h5';
    else
        file_name_s = trim(Sfold_name)//Env_sep//sim_name//'SmtableES_'//stFreq//trim(freq_unit)//'_'//nom_meth_exact//'.h5';
    endif

    !COMPLEX(real64), Dimension(NRx_tot,4*NTr), INTENT(IN):: S_total
    
    Allocate(S_towrite(8*NTr,NRx));
    Do ii=1,NTr
        Do jj=1,NRx
            S_towrite(8*(ii-1)+1,jj) = real(S_total(jj,4*(ii-1)+1));
            S_towrite(8*(ii-1)+2,jj) = imag(S_total(jj,4*(ii-1)+1));
            S_towrite(8*(ii-1)+3,jj) = real(S_total(jj,4*(ii-1)+2));
            S_towrite(8*(ii-1)+4,jj) = imag(S_total(jj,4*(ii-1)+2));
            S_towrite(8*(ii-1)+5,jj) = real(S_total(jj,4*(ii-1)+3));
            S_towrite(8*(ii-1)+6,jj) = imag(S_total(jj,4*(ii-1)+3));
            S_towrite(8*(ii-1)+7,jj) = real(S_total(jj,4*(ii-1)+4));
            S_towrite(8*(ii-1)+8,jj) = imag(S_total(jj,4*(ii-1)+4));
        EndDo
    Enddo
    !S_towrite(1:NRx,1:8*NTr) = 0.27       
     
     !dimsf(1) = 8*NTr
     !dimsf(2) = NRx
     
     !! **********************************************************
     !! **********************************************************
     !! COMMENTED WAITING for PHDF5
     
!     !
!     ! Initialize FORTRAN interface
!     !
!     if (rank == 0) then 
!     Write(*,*) 'Here 1'
!     endif
!     CALL h5open_f(error)
!     if (rank == 0) then 
!     Write(*,*) 'Here 2'
!     endif
!     !
!     ! Setup file access property list with parallel I/O access.
!     !
!     CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
!     if (rank == 0) then 
!     Write(*,*) 'Here 3'
!     endif
!     CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)
!     if (rank == 0) then 
!     Write(*,*) 'Here 4'
!     endif
!     !
!     ! Figure out the filename to use.  If your system does not support
!     ! getenv, comment that statement with this,
!     ! filename = ""
!!     CALL getenv("HDF5_PARAPREFIX", filename)
!!     fnamelen = LEN_TRIM(filename)
!!     if ( fnamelen == 0 ) then
!!	filename = default_fname
!!     else
!!	filename = filename(1:fnamelen) // "/" // default_fname
!!     endif
!!     print *, "Using filename = ", filename
!
!     !
!     ! Create the file collectively.
!     !
!     if (rank == 0) then 
!        Write(*,*) 'file_name_s = ',file_name_s
!     endif
!     CALL h5fcreate_f(file_name_s, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
!     if (rank == 0) then 
!     Write(*,*) 'Here 5'
!     endif
!     CALL h5pclose_f(plist_id, error)
!     !
!     ! Create the data space for the  dataset.
!     !
!     if (rank == 0) then 
!     Write(*,*) 'Here 6'
!     endif
!     CALL h5screate_simple_f(data_rank, dimsf, filespace, error)
!     if (rank == 0) then 
!     Write(*,*) 'Here 7'
!     endif
!     !
!     ! Create the dataset with default properties.
!     !
!     if (rank == 0) then 
!     Write(*,*) 'Here 8'
!     endif
!     CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, filespace, &
!                      dset_id, error)
!     !
!     ! Create property list for collective dataset write
!     !
!     if (rank == 0) then 
!     Write(*,*) 'Here 9'
!     endif
!     CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
!     if (rank == 0) then 
!     Write(*,*) 'Here 10'
!     endif
!     CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
!     if (rank == 0) then 
!     Write(*,*) 'Here 11'
!     endif
!     !
!     ! For independent write use
!     ! CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
!     !
!
!     !
!     ! Write the dataset collectively.
!     !
!     if (rank == 0) then 
!     Write(*,*) 'Here 12'
!     endif
!     CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, S_towrite, dimsfi, error, &
!                      xfer_prp = plist_id)
!     !
!     ! Deallocate data buffer.
!     !
!     if (rank == 0) then 
!     Write(*,*) 'Here 13'
!     endif
!     DEALLOCATE(S_towrite)
!
!     !
!     ! Close resources.
!     !
!     CALL h5sclose_f(filespace, error)
!     CALL h5dclose_f(dset_id, error)
!     CALL h5pclose_f(plist_id, error)
!     CALL h5fclose_f(file_id, error)
!     ! Attempt to remove the data file.  Remove the line if the compiler
!     ! does not support it.
!     !CALL unlink(filename)
!
!     !
!     ! Close FORTRAN interface
!     !
!     CALL h5close_f(error)
     !! COMMENTED BECAUSE WINDOWS
     !! **********************************************************
     !CALL MPI_FINALIZE(mpierror)   

END SUBROUTINE Write_Sfiles

SUBROUTINE Write_txt_Sfiles(nom_methode,Transmitters,Receivers,S_total)

    USE Initialization
    USE common_variables
    USE MPI
    
    IMPLICIT NONE

    !IN/OUT 
    character(8), INTENT(IN):: nom_methode
    type (Dipole), Dimension(NTr), INTENT(IN) :: Transmitters
    type (Dipole), Dimension(NRx_tot), INTENT(IN) :: Receivers
    COMPLEX(real64), Dimension(NRx_tot,4*NTr), INTENT(IN):: S_total
    
    ! local 
    Integer :: ii,jj,kkt,kkr,a,Nths,Nphs
    Integer :: id,nthreads,p,d,NTr_wr_proc 
    Real(kind=8) :: th_i,ph_i
    COMPLEX(real64) :: Vv, Vh, Hv, Hh
    
    CHARACTER(:), allocatable:: nom_meth_exact,stFreq,sim_name
    CHARACTER(200) :: file_name_s,Sfold_name
    CHARACTER(6) :: ty,kkt_st
    
    Sfold_name = trim(SimOutfld_name)//Env_sep//'S_files';
    
    NTr_wr_proc = (NTr/nber_procs)+1;
    If (nom_methode=='CBFM-E  ') Then
        Allocate(character(6) ::nom_meth_exact)
        nom_meth_exact = trim(nom_methode)
    Else ! MoM or RGE 
        Allocate(character(3) ::nom_meth_exact)
        nom_meth_exact = trim(nom_methode)
    Endif
     
    
    a = nint(Freq_w/10.**freq_mag);
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

    If (Nfreq == 1) Then
        Allocate(character(1)::sim_name)
        sim_name= ''
    Else
        if (num_freq < 10) Then
            Allocate(character(5)::sim_name)
            Write(sim_name,'(a,i1,a)') 'Sim', num_freq, '_'
        ElseIf (num_freq < 100) Then
            Allocate(character(6)::sim_name)
            Write(sim_name,'(a,i2,a)') 'Sim', num_freq, '_'
        Else
            Allocate(character(7)::sim_name)
            Write(sim_name,'(a,i3,a)') 'Sim', num_freq, '_'
        EndIf        
    EndIf 
    Write(stFreq,ty) Freq_w/10.**freq_mag
    

    DO kkt=1,NTr
          th_i =  Transmitters(kkt)%theta;
          ph_i = Transmitters(kkt)%phi;
          
              
              ! Comment here if PHDF5 S_files successful
              If (wr_Sij .eq. 1) Then
                  If ((kkt .gt. (rank*NTr_wr_proc)) .and. (kkt .le. (rank+1)*NTr_wr_proc)) then 
                    Write(kkt_st,'(a,i4.4)') 'kt',kkt;
                    if (EqSph == 0) then
                        file_name_s = trim(Sfold_name)//Env_sep//sim_name//'Smtable_'//stFreq//freq_unit//'_'//trim(kkt_st)//'_'//nom_meth_exact//'.dat';
                    else
                        file_name_s = trim(Sfold_name)//Env_sep//sim_name//'SmtableES_'//stFreq//freq_unit//'_'//trim(kkt_st)//'_'//nom_meth_exact//'.dat';
                    endif
                    Open(unit=21+rank,File = file_name_s)  
                    Write(21+rank, '(a,f10.4,a,f10.4)') 'THETA =',  Transmitters(kkt)%theta, '; PHI =',  Transmitters(kkt)%phi  
                    Write(21+rank,'(a,a)') '      theta       phi       Re(Svv)        Im(Svv)         Re(Svh)       Im(Svh) ',&
                                    '        Re(Shv)       Im(Shv)        Re(Shh)        Im(Shh) '
                    Do kkr=1, NRx  
                        Vv = S_total(kkr,4*(kkt-1)+1);
                        Vh = S_total(kkr,4*(kkt-1)+2); 
                        Hv = S_total(kkr,4*(kkt-1)+3); 
                        Hh = S_total(kkr,4*(kkt-1)+4);
                          
                        Write(21+rank,'(f9.2,a,f9.2,a,e12.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4)') &
                        Receivers(kkr)%theta,';  ',Receivers(kkr)%phi,';  ',Real(Vv),';  ',Imag(Vv),';  ',Real(Vh),&
                      ';  ',Imag(Vh),';  ', Real(Hv),';  ',Imag(Hv),';  ',Real(Hh),';  ',Imag(Hh)
                    EndDo
                    Close(21+rank);  
                  Endif  
              EndIf  
    EndDo
END SUBROUTINE Write_txt_Sfiles

SUBROUTINE Write_jobs_sim_info(CBFM_Blocks,MPI_CBFM_Blocks,K_patchs_all)

    USE Initialization
    USE common_variables
    USE MPI
    
    Implicit NONE
    
    !IN/OUT 
    type (CBFM_Block), Dimension(Nblocks), INTENT(IN) :: CBFM_Blocks
    Integer, Dimension(nber_procs,Nblk_proc_max+1), INTENT(IN) :: MPI_CBFM_Blocks
    Integer, Dimension(NBlocks), INTENT(IN) :: K_patchs_all
    
    ! local 
    Integer ii, Nblocks_proc, Ncells_proc, K_proc, a 
    CHARACTER(300) file_name,analysis_fold_name
    CHARACTER(:) ,allocatable :: stFreq
    CHARACTER(6) :: ty
    logical :: dirExists
    
    analysis_fold_name = trim(SimOutfld_name)//Env_sep//'Analysis';
    
    Call MPI_Barrier(MPI_COMM_WORLD,code);
    if (rank == 0) then         
        a = nint(Freq_w/(10**freq_mag));
        if (a < 10) Then 
            Allocate(character(5) ::stFreq); ty = '(f5.3)';
        ElseIf (a < 100) Then
            Allocate(character(6) ::stFreq); ty = '(f6.3)';
        Else
            Allocate(character(7) ::stFreq); ty = '(f7.3)';
        EndIf     
        Write(stFreq,ty) Freq_w/1E9
    
        file_name = trim(analysis_fold_name)//Env_sep//'MPIjobs_loadinfo_'//stFreq//trim(freq_unit)//'.dat'; 
    
        Open(14,File = trim(file_name))
        Write(14,'(a)') '   rank       Nblocks      Nbc          K';
        Do ii=1, nber_procs
            Nblocks_proc = MPI_CBFM_Blocks(ii,1);
            Ncells_proc = sum(CBFM_Blocks(MPI_CBFM_Blocks(ii,2:1+Nblocks_proc))%Nbc_b);
            K_proc = sum(K_patchs_all(MPI_CBFM_Blocks(ii,2:1+Nblocks_proc)));
            Write(14,'(i6,i12,i12,i12)') (ii-1), Nblocks_proc, Ncells_proc,K_proc
        EndDo
        Close(14)
    endif   
    Call MPI_Barrier(MPI_COMM_WORLD,code);
    
    !!!! Let's recap
    !! quelles sont les allocations majeures qui font la diff entre les MPI tasks/jobs
    !! on remarque que les deux valeurs qui menent essentiellement la dance sont Nbc_proc et Ktot_proc_max
    !Allocate(C_job_patchs(3*Nbc_proc,2*NTr_CBFM));
    !Allocate(C_job_patchs_tmp(3*Nbc_proc,Kmax_job));
    !Allocate(C_job_patchs(3*Nbc_proc,Kmax_job));    
    !
    !Allocate(Zreduite(Ktot_proc_max,K_total)); 
    !
    !Allocate(C_trans_patchs(3*Nbc_proc_max,K_max));
    !Allocate(curs_B_Cpatch_Next(NBlocks_EffNextjob));
    !
    !Allocate(Vreduit(Ktot_proc_max,2*NTr));    
    !
    !Allocate(ZredLoc(Mlocal,Nlocal),VredLoc(Mlocal,NRHSlocal));    
    !Allocate(AlphaProc(Ktot_proc_max,2*NTr));
    !
    !Allocate(E_total(3*Nbc_proc,2*NTr));
    !
    !
    !
    !!! mnt les allocations communes et/ou qui n'ont pas vraiment de poids majeurs sont :
    !Allocate(curs_B_cel(NBlocks),curs_B_Cpatch(MyNBlocks),K_patchs(MyNBlocks),K_patchs_all(NBlocks))
    !Allocate(blocks_to_sort(NBlocks),OrderBlocksProcs(Nblocks));
    !Ntests = 4; Allocate(testedBlks(Ntests))
    !Allocate(ExtSizes(NBlocks));
    !Allocate(fSR_blocks(MyNBlocks),nnz_blocks(MyNBlocks),spr_perc_blocks(MyNBlocks));
    !Allocate(Cells_Block(size));
    !
    !Allocate(EREFpatch_e(3*size,2*NTr_CBFM));
    !Allocate(Epatch_e(3*size,2*NTr_CBFM));
    !Allocate(Zpatch_e(2*klu+1,3*size));
    !Allocate(Zpatch_e_spr(nnz),row_sprZ(3*size+1),col_sprZ(nnz));
    !Allocate(Zpatch_e(3*size,3*size));
    !Allocate( AA(M,N), S(MIN(M,N)),U(M,M),VT(N,N), WW(MIN(M,N)-1));
    !Allocate(Cpatch_e(3*size, K));
    !
    !ALLOCATE(nb_elements_rec(nber_procs), deplts(nber_procs));
    !Allocate(vect_tmp(NBlocks));
    !Allocate(Ktot_procs(nber_procs));
    !Allocate(curs_row_Z_blocks(MyNBlocks));
    !Allocate(curs_col_Z_blocks(NBlocks)) ;
    !
    !Allocate(Cells_Block_ii(size1));
    !Allocate(Matrice1(NbreLig_mat1,NbreCol_mat1));
    !Allocate(Cells_Block_jj(size3));
    !Allocate(Matrice2(NbreLig_mat1,NbreLig_mat3));
    !Allocate(Matrice3(NbreLig_mat3,NbreCol_mat3));
    !Allocate(Mat_Inter(NbreCol_mat1,NbreLig_mat3));
    !Allocate(MatProduit(NbreCol_mat1,NbreCol_mat3));
    !
    !Allocate(Cells_Block_ii(size));
    !Allocate(E_ref_incident(NbreLig_mat1,2*NTr));
    !Allocate(Matrice1(NbreLig_mat1,NbreCol_mat1));
    !Allocate(Vect(NbreLig_mat1));
    !Allocate(VectProduit(NbreCol_mat1));
    !
    !Allocate(DESCA(9),DESCB(9));
    !Allocate(Curs_Kpatchs_all(NBlocks+1)); ! equivalent to curs_B_cel for Ni;
    !Allocate(K_patchs_eff(NBlocks_eff));
    !Allocate(IPIV(Mlocal+M_B));
    !
    !Allocate(K_patchs_eff(NBlocks_eff));
    !Allocate(VectProduit(BlockSize)); 
    
END SUBROUTINE Write_jobs_sim_info

subroutine system_mem_usage(valueRSS)

    use ifport !if on intel compiler

    ! You should know that : RSS is Resident Set Size (physically resident memory - 
    ! this is currently occupying space in the machine's physical memory),
    ! and VSZ is Virtual Memory Size (address space allocated - this has addresses 
    ! allocated in the process's memory map, but there isn't necessarily any actual 
    ! memory behind it all right now).

    ! This subroutine is from the following stackoverflow discussion :
    ! http://stackoverflow.com/questions/22028571/track-memory-usage-in-fortran-90?
    ! answertab=oldest#tab-top

    implicit none

    integer, intent(out) :: valueRSS

    character(len=200):: filename=' '
    character(len=80) :: line
    character(len=8)  :: pid_char=' '
    integer :: pid
    logical :: ifxst

    valueRSS=-1    ! return negative number if not found

    !--- get process ID

    pid=getpid()
    write(pid_char,'(I8)') pid
    filename='/proc/'//trim(adjustl(pid_char))//'/status'

    !--- read system file

    inquire (file=filename,exist=ifxst)
    if (.not.ifxst) then
      write (*,*) 'system file does not exist'
      return
    endif

    open(unit=100, file=filename, action='read')
    do
      read (100,'(a)',end=120) line
      if (line(1:6).eq.'VmRSS:') then
         read (line(7:),*) valueRSS
         exit
      endif
    enddo
    120 continue
    close(100)

    return
end subroutine system_mem_usage
    
subroutine print_allocate(Nchar,allocate_str,type_str,size)

    USE Initialization
    USE common_variables
    USE MPI
    
    IMPLICIT NONE

    !IN/OUT
    Integer, INTENT(IN) :: size,Nchar 
    character(Nchar), INTENT(IN) :: allocate_str
    character(5), INTENT(IN) ::type_str ! D for Double and S for Single REAL, COMP or INTG
        
    ! local
    Real(kind=8) :: size_MB
    character(300) :: analysis_fold_name,file_name
    character(19) :: time_allocate
    character(:), allocatable :: rank_str 
    character(8)  :: date
    character(10) :: time
    character(5)  :: zone
    integer,dimension(8) :: values
    
    if ((track_memory == 1) .and. (rank .lt. Njob_max)) then 
        call date_and_time(date,time,zone,values);
        time_allocate = date(5:6)//'-'//date(7:8)//'-'//date(1:4)//'_'//time(1:2)//':'//time(3:4)//':'//time(5:6);
    
        analysis_fold_name = trim(SimOutfld_name)//Env_sep//'Analysis';
    
        if (rank .lt. 10) then 
            allocate(character(1) ::rank_str);
            Write(rank_str,'(i1)') rank;
        elseif (rank .lt. 100) then 
            allocate(character(2) ::rank_str);
            Write(rank_str,'(i2)') rank;
        elseif (rank .lt. 1000) then 
            allocate(character(3) ::rank_str);
            Write(rank_str,'(i3)') rank;
        elseif (rank .lt. 10000) then 
            allocate(character(4) ::rank_str);
            Write(rank_str,'(i14)') rank;
        endif
    
    
        file_name = trim(analysis_fold_name)//Env_sep//'TrackAllocate_j'//rank_str//'.dat'; 
    
        if (trim(allocate_str) .eq. 'Reference(t=0)') then ! reference print allocate
            Open(30+rank,File = trim(file_name)); 
        else        
            Open(30+rank,File = trim(file_name), status = 'old', position = 'append'); 
        endif
    
        if (type_str .eq. 'DCOMP') then 
            size_MB = 64.*2.*size/1e6;
        elseif (type_str .eq. 'DREal') then 
            size_MB = 64.*size/1e6;
        else
            size_MB = 32.*size/1e6       
        endif
    
    
        Write(30+rank,'(a,i16,a10,f12.3,a,a)') time_allocate, size,type_str,size_MB,'    ',allocate_str;   
        Close(30+rank);
    endif
    
    

End Subroutine print_allocate
    
    
    
    
