SUBROUTINE Compute_EFields_MoM(Cells,Transmitters,Receivers,S_total,C_ext,C_abs)  
    
    ! Created on 8-15-2020 to track the high m value problem with spherical particle
    ! Now Compute_EFields includes calculating the scattering matrix and ExtAbsCsec fromIntField
    
    USE Initialization
    USE common_variables
    USE lapack95
    USE DiverseUtil
    USE MPI

    IMPLICIT NONE
    
    !! IN/OUT ******************************************************************
    type (Cell), Dimension(Nbc), INTENT(IN) :: Cells
    type (Dipole), Dimension(NTr), INTENT(IN) ::Transmitters
    type (Dipole), Dimension(NRx_tot), INTENT(IN) :: Receivers
    COMPLEX(real64), Dimension(NRx_tot,4*NTr), INTENT(OUT) :: S_total
    COMPLEX(real64), Dimension(NTr), INTENT(OUT) :: C_ext,C_abs
    
    !! Local 
    Integer :: ii, jj, kk, Nbc_proc_max,Nd,Nr,num_cel_min,num_cel_max
    Integer :: prev_job, next_job,IGlob_min,IGlob_max, Nelts
    Integer :: size2, pp, Eff_rank, iLoc_proc 
    DOUBLE COMPLEX :: Ztot_i_j, Vtot_i_j
    
    character(8)  :: date_init_N1, date_final_N1
    character(10) :: time_init_N1, time_final_N1
    character(5)  :: zone_init_N1, zone_final_N1
    Integer,dimension(8) :: values_init_N1, values_final_N1
    Integer, dimension(4):: time_calcul_N1
    character(200) :: file_name
    
    COMPLEX(real64), Dimension(:,:), allocatable :: E_total
    
    !! Scalapack :
    Integer :: icontxt,ICSRC,IRSRC,IA,JA,IB,JB,LWORK,LIWORK
    INTEGER :: INDXL2G,iloc,jloc,iGlob,jGlob
    INTEGER, EXTERNAL :: numroc
    REAL(kind=8), EXTERNAL :: pzlange
    REAL, EXTERNAL :: pslamch
    INTEGER, parameter :: ROW_SRC = 0, COL_SRC = 0
    INTEGER, parameter :: NDIMS = 2
    INTEGER, dimension(MPI_STATUS_SIZE) :: status 
    INTEGER,dimension(1:NDIMS) :: dims
    INTEGER, dimension(:), allocatable :: DESCA, DESCB, IPIV
    DOUBLE COMPLEX, Dimension(:,:), allocatable :: ZLoc,VLoc
    Real(kind=8) :: RCOND,EPSMCH,ANORM,ERRBD
    Real(kind=8), Dimension(:), allocatable :: WORK      
    Integer, Dimension(:), allocatable  :: IWORK 
    
    !! Local Scattered Fields
    Integer :: I,Ic,cel_beg,cel_end,num_emetteur, num_capteur
    Integer, Dimension(nber_procs) :: all_Nbc_procs    
    Real(kind=8) :: theta_capteur, phi_capteur
    Real(kind=8) :: Cext_e_V,Cext_e_H, Cabs_e_V,Cabs_e_H;
    COMPLEX(real64), Dimension(:,:), allocatable :: Green_dt_app,S_total_all
    COMPLEX(real64), Dimension(:), allocatable ::S_total_capteur,S_total_capteur_all
    COMPLEX(real64), Dimension(3) :: E_v, E_h, E_v_p, E_h_p
    COMPLEX(real64) :: Vv, Vh, Hv, Hh     
    COMPLEX(real64), Dimension(:), allocatable :: C_ext_all,C_abs_all
    COMPLEX(real64), Dimension(:,:),allocatable:: E_ref_incident
    type(Cell), Dimension(:), allocatable :: Cells_proc     
       
    
    !! Etot inside scatterer **********************************************************************************************
    !! Etot inside scatterer **********************************************************************************************
    !! Etot inside scatterer **********************************************************************************************
    !! Etot inside scatterer **********************************************************************************************
   
   call MPI_BARRIER(MPI_COMM_WORLD,code);
   If (rank == 0) Then 
        Write (*,*) ''
        Write (*,'(a)') '>>>>>>>>>>>>>>>>>>>>>>>>> Method Of Moments <<<<<<<<<<<<<<<<<<<<<<<<<'
        Write (*,*) ''; Write (*,*) ''
        Write (*,'(a)') '----------------- Total Fields Inside Scatterer ---------------------'
        Write (*,*) ''
    EndIf 
    
    Nbc_proc = Nbc/nber_procs;
    Nbc_proc_max = Nbc_proc+1;
        
    if (rank .lt. mod(Nbc,nber_procs)) then 
        Nbc_proc = Nbc_proc + 1;
    endif
    
    if (rank .eq. 0) then 
        Write(*,'(a,i4)') 'Nbc_proc ~= ', Nbc_proc
    endif         
    
    ! will be used more than once, so better create it here
    Nd = Nbc/nber_procs; 
    Nr = mod(Nbc,nber_procs);
    all_Nbc_procs(1:nber_procs) = Nd;
    Do ii=1, nber_procs        
        If (ii .le. Nr) then 
            all_Nbc_procs(ii) = all_Nbc_procs(ii) + 1;      
        Endif        
    EndDo 
    
    if (rank == 0) then 
        Write(*,'(a,i10)') 'Resolution of the MoM system of linear equations of 3Nbc =',3*Nbc
        time_calcul_N1 = 0
        call date_and_time(date_init_N1,time_init_N1,zone_init_N1,values_init_N1)  
    endif        
    
    !! Zloc 
    Allocate(DESCA(9),DESCB(9));
    DESCA = 0; DESCB = 0;
    IRSRC = 0; ICSRC =0;

    CALL MPI_DIMS_CREATE(nber_procs,NDIMS,dims,INFO)
    NPROW = dims(1); NPCOL=dims(2);  
    
    CALL BLACS_GET(-1,0,icontxt);
    CALL BLACS_GRIDINIT(icontxt,'Row-major', NPROW, NPCOL);
    CALL BLACS_GRIDINFO(icontxt,NPROW,NPCOL,MYROW,MYCOL);
    
    ! the following three lines are equivalent to : (coded to check/understand how numroc operates)
    ! call numroc_2d_opt(3*Nbc,3*Nbc,2*NTr,NPROW,NPCOL,Myrow,Mycol,Mlocal,Nlocal,NRHSlocal);        
    Mlocal = numroc(3*Nbc,M_B,Myrow,ROW_SRC,NPROW);
    Nlocal = numroc(3*Nbc,N_B,Mycol,COL_SRC,NPCOL);
    NRHSlocal = numroc(2*NTr,N_B,Mycol,COL_SRC,NPCOL); 
        
!    !if (rank .eq. 0) then 
!        Write(*,'(i3,a,i6,a,i6,a,i6)') rank, ' : Mloc ~= ', Mlocal, ', Nloc ~= ', Nlocal,' and NRHSlocal ~= ',NRHSlocal
!    !endif   
    
    
    Allocate(ZLoc(Mlocal,Nlocal),VLoc(Mlocal,NRHSlocal));
    call print_allocate(19,'ZLoc(Mlocal,Nlocal)','DCOMP',Mlocal*Nlocal);
    call print_allocate(23,'VLoc(Mlocal,NRHSlocal)','DCOMP',Mlocal*NRHSlocal);    
    
    !! here I need to fill out ZLoc and VredLoc
    ! from Block_Cyclic_Distribution_v1
    call MPI_BARRIER(MPI_COMM_WORLD,code);
    
    Do iLoc=1, Mlocal
        iGlob = INDXL2G(iLoc,M_B,Myrow,0,NPROW);     
        ! ZLoc
        Do jLoc=1, Nlocal
            jGlob = INDXL2G(jLoc,N_B,Mycol,0,NPCOL)
            call compute_Green_s_tr_single_element(Cells,iGlob,jGlob,Ztot_i_j);
            ZLoc(iLoc,jLoc) = Ztot_i_j;
        End do  
        Do jLoc=1, NRHSlocal
            jGlob = INDXL2G(jLoc,N_B,Mycol,0,NPCOL)
            call Incident_Field_single_element(Cells,Transmitters,iGlob,jGlob,Vtot_i_j)
            VLoc(iLoc,jLoc) = Vtot_i_j; 
        End do 
    End do      
    
    ! RESOLUTION WITH SCALAPACK ------------------------------------------------    
    CALL DESCINIT(DESCA,3*Nbc,3*Nbc,M_B,N_B,IRSRC,ICSRC,icontxt,Mlocal,INFO)
    CALL DESCINIT(DESCB,3*Nbc,2*NTr,M_B,N_B,IRSRC,ICSRC,icontxt,Mlocal,INFO)
            
    IA=1;JA=1;IB=1;JB=1;       
    ! Machine precision
    EPSMCH = PSLAMCH(icontxt,'E'); 
    
    ! LWORK >= 2*LOCr(N+MOD(IA-1,MB_A)) +  MAX( 2, MAX(NB_A*CEIL(NPROW-1,NPCOL),LOCc(N+MOD(JA-1,NB_A)) + NB_A*CEIL(NPCOL-1,NPROW)) ).
    ! LRWORK >= MAX( 1, 2*LOCc(N+MOD(JA-1,NB_A)) ).
    LWORK = 2*3*Nbc; ! Here 2*Matrix_size is enough for now (same as for CBFM-E), but needs to compute the exact LWORK and LIWORK for more accuracy/robustness
    LIWORK = 2*3*Nbc; !K_total; 
    Allocate(WORK(LWORK),IWORK(LIWORK));
    
    ! get Infinity NORM of ZLoc 
    ANORM = PZLANGE( 'I', 3*Nbc,3*Nbc, ZLoc, IA,JA,DESCA,WORK);
    Allocate(IPIV(Mlocal+M_B));
    if (rank == 0) then 
        Write(*,'(a)') 'PZGESV in progress ...' 
    endif
    CALL PZGESV(3*Nbc,2*NTr,ZLoc,IA,JA,DESCA,IPIV,VLoc,IB,JB,DESCB,INFO);
    if (rank .eq. 0) then
        Write(*,'(a,i4)') 'INFO = ',INFO
    endif
    
    if (INFO .GT. 0) then 
        if (rank .eq. 0) then 
            Write(*,'(a)') 'SINGULAR MATRIX';
        endif            
    elseif (3*Nbc .gt. 0) then 
        ! Get Reciprocal condition number RCOND of Zc
        !CALL PZGECON( 'I', 3*Nbc,ZLoc,IA,JA,DESCA, ANORM, RCOND,WORK,LWORK,IWORK,LIWORK,INFO);
        !RCOND = max(RCOND,EPSMCH);
        !ERRBD = EPSMCH/ RCOND
        if (rank == 0) then 
          Write(*,'(a,es10.3)') '--> ANORM(ZLoc) = ',ANORM  
          !Write(*,'(a,es12.5)') '--> RCOND = ',RCOND
          !Write(*,'(a,es12.5,a,es12.5)') '--> With EPSMCH =',EPSMCH,'; ERRBD =',ERRBD
        endif
    endif  
    
    CALL BLACS_GRIDEXIT(icontxt);
    Deallocate(ZLoc);         
    
    !! Now reorganize back the solution directly in E_total ********************************************
    Allocate(E_total(3*Nbc_proc_max,2*NTr));
    call print_allocate(26,'E_total(3*Nbc_proc,2*NTr)','DCOMP',3*Nbc_proc*2*NTr);
    E_total = 0D0
        
    prev_job = mod(nber_procs+rank-1,nber_procs)
    next_job = mod(rank+1,nber_procs)
    size2 = 3*Nbc_proc_max*2*NTr
   
    Do pp=1, nber_procs
        Call MPI_SENDRECV_REPLACE(E_total,2*size2,MPI_COMPLEX,prev_job,tag,next_job,tag,MPI_COMM_WORLD,status,code)
        
        Eff_rank = mod(rank+pp,nber_procs);         
        num_cel_min = sum(all_Nbc_procs(1:Eff_rank))+1; !Eff_rank*Nd+1 + min(Nr,Eff_rank);   
        num_cel_max = sum(all_Nbc_procs(1:Eff_rank+1)); !(Eff_rank+1)*Nd + min(Nr,Eff_rank);             
        IGlob_min = 3*(num_cel_min-1)+1;          
        IGlob_max = 3*num_cel_max ;
	    
        Do iLoc=1, Mlocal
    	    iGlob = INDXL2G(iLoc,M_B,Myrow,0,NPROW)
            iLoc_proc = iGlob - IGlob_min + 1; 
            if ((iGlob .ge. IGlob_min) .AND. (iGlob .le. IGlob_max)) then
                Do jLoc=1, NRHSlocal
            		jGlob = INDXL2G(jLoc,N_B,Mycol,0,NPCOL)
            		E_total(iLoc_proc,jGlob) = VLoc(iLoc,jLoc)
                End do
            endif
        enddo        
    End do    
                
    if (rank == 0) then 
        call date_and_time(date_final_N1,time_final_N1,zone_final_N1,values_final_N1)
        call Calcul_time_spent (values_init_N1,values_final_N1, time_calcul_N1)
        Write (*,'(a,i2,a,i2,a,i2,a,i2,a)') '--> to resolve Z * E = Einc ',time_calcul_N1(1),'j',time_calcul_N1(2)&
        ,'h',time_calcul_N1(3),'min', time_calcul_N1(4),'sec';
    endif     
    
    If ((save_Eint .eq. 1) .and. (Nbc .le. save_Eint_Nmax)) then         
        file_name = trim(SimOutfld_name)//Env_sep//'Ein_MPI.dat';
        if (rank == 0) then 
            Write(*,'(a)',advance='no') '--> to write Ein_tot ' 
            time_calcul_N1 = 0
            call date_and_time(date_init_N1,time_init_N1,zone_init_N1,values_init_N1);
        endif      
        
        if (rank == 0) then
            open(unit=14, file = trim(file_name));!, form = 'unformatted')
            Do ii=1,3*Nbc_proc
                Do kk=1,2*NTr          
                    Write(14,'(es16.8,a,es16.8)') real(E_total(ii,kk)),';',imag(E_total(ii,kk)); 
                EndDo
            EndDo
            close(14)
        endif
        call MPI_BARRIER(MPI_COMM_WORLD,code);
        do pp = 1, nber_procs - 1
            if(rank == pp) then
                open(unit = 14, file = trim(file_name), status = 'old', position = 'append')
                Do ii=1,3*Nbc_proc        
                    Do kk=1,2*NTr                           
                        Write(14,'(es16.8,a,es16.8)') real(E_total(ii,kk)),';',imag(E_total(ii,kk)); 
                    EndDo
                EndDo
                close(14)
            endif
            call MPI_BARRIER(MPI_COMM_WORLD,code);
        enddo      
        
        if (rank == 0) then 
            call date_and_time(date_final_N1,time_final_N1,zone_final_N1,values_final_N1)
            call Calcul_time_spent (values_init_N1,values_final_N1, time_calcul_N1)
            Write (*,'(i2,a,i2,a,i2,a,i2,a)') time_calcul_N1(2),'h', time_calcul_N1(3),'min', time_calcul_N1(4),'sec';
        endif
    endif       
    
    
    !! SCATTERING MATRIX **********************************************************************************************
    !! SCATTERING MATRIX **********************************************************************************************
    !! SCATTERING MATRIX **********************************************************************************************
    !! SCATTERING MATRIX **********************************************************************************************
    
    If (rank == 0) Then 
        Write (*,*) ''
        Write (*,'(a)') '----------------------- Scattred Fields --------------------------'
        Write (*,*) ''
        
        time_calcul_N1 = 0
        call date_and_time(date_init_N1,time_init_N1,zone_init_N1,values_init_N1)  
    EndIf
    
    
    ! now every proc focus on its cells/part of E_tot 
    Allocate(Cells_proc(Nbc_proc));        
    cel_beg = sum(all_Nbc_procs(1:rank))+1; 
    cel_end = sum(all_Nbc_procs(1:rank+1)); 
    Cells_proc(1:Nbc_proc) = Cells(cel_beg:cel_end);         
    
    DO num_capteur =1,NRx_tot  
	    Allocate(S_total_capteur(4*NTr),S_total_capteur_all(4*NTr));   
        Allocate(Green_dt_app(3,3*Nbc_proc))
        !! Dyade de Greene singuliere
        theta_capteur = Receivers(num_capteur)%theta
        phi_capteur = Receivers(num_capteur)%phi
        
        Call Green_ff_dt(Nbc_proc,Cells_proc,Receivers,num_capteur,theta_capteur,phi_capteur,Green_dt_app)            
        DO num_emetteur=1,NTr
            E_v = 0
            E_h = 0      
            DO I=1,3*Nbc_proc ! This the difference with the serial OpenMP code, Here each process computes S based on its part of E_total
                E_v(1)=E_v(1)+Green_dt_app(1,I)*E_total(I,num_emetteur)
                E_v(2)=E_v(2)+Green_dt_app(2,I)*E_total(I,num_emetteur)
                E_v(3)=E_v(3)+Green_dt_app(3,I)*E_total(I,num_emetteur)
            
                E_h(1)=E_h(1)+Green_dt_app(1,I)*E_total(I,num_emetteur+NTr)
                E_h(2)=E_h(2)+Green_dt_app(2,I)*E_total(I,num_emetteur+NTr)
                E_h(3)=E_h(3)+Green_dt_app(3,I)*E_total(I,num_emetteur+NTr)    
            ENDDO
            
            !! ---------------------------------------------------------------------------------!!
            !! ------------------------------Polarisation Vv -----------------------------------!!
            !! ---------------------------------------------------------------------------------!! 
            Vv= - E_v(1)*sin(theta_capteur*Pi/180.)+E_v(2)*cos(theta_capteur*Pi/180.)*cos(phi_capteur*Pi/180.) &
            +E_v(3)*cos(theta_capteur*Pi/180.)*sin(phi_capteur*Pi/180.);
                
    
            !! ---------------------------------------------------------------------------------!!
            !! ------------------------------Polarisation Vh -----------------------------------!!
            !! ---------------------------------------------------------------------------------!! 
            Vh= - E_v(2)*sin(phi_capteur*Pi/180.)+E_v(3)*cos(phi_capteur*Pi/180.);
                
    
            !! ---------------------------------------------------------------------------------!!
            !! ------------------------------Polarisation Hv -----------------------------------!!
            !! ---------------------------------------------------------------------------------!! 
            Hv= - E_h(1)*sin(theta_capteur*Pi/180.)+E_h(2)*cos(theta_capteur*Pi/180.)*cos(phi_capteur*Pi/180.) &
            + E_h(3)*cos(theta_capteur*Pi/180.)*sin(phi_capteur*Pi/180.);
                
            !! ---------------------------------------------------------------------------------!!
            !! ------------------------------Polarisation Hh -----------------------------------!!
            !! ---------------------------------------------------------------------------------!! 
            Hh= - E_h(2)*sin(phi_capteur*Pi/180.) + E_h(3)*cos(phi_capteur*Pi/180.);
             
            !! ---------------------------------------------------------------------------------!!
            !! -----------------------Remplissage du vecteur S_total ---------------------------!!
            !! ---------------------------------------------------------------------------------!!
            S_total_capteur(4*(num_emetteur-1)+1)  = Vv; 
            S_total_capteur(4*(num_emetteur-1)+2)  = Vh; 
            S_total_capteur(4*(num_emetteur-1)+3)  = Hv; 
            S_total_capteur(4*(num_emetteur-1)+4)  = Hh; 
        ENDDO      
        Deallocate(Green_dt_app);        
        Nelts = 4*NTr;
        Call MPI_ALLREDUCE(S_total_capteur,S_total_capteur_all,Nelts,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,code);
        S_total(num_capteur,1:4*NTr) = S_total_capteur_all(1:4*NTr);
        deallocate(S_total_capteur,S_total_capteur_all);        
    Enddo
    
   
    !! Cext & Cabs **********************************************************************************************
    !! Cext & Cabs ***********************************************************************************************
    !! Cext & Cabs *********************************************************************************************
    !! Cext & Cabs *********************************************************************************************
    
    if (rank == 0) then
        Write (*,*) '' 
        Write (*,'(a)') '-------- Extinction and Absorption from Internal Fields ----------'  
        Write (*,*) ''
    endif
    
    ! All procs recover again this important information 
    DO num_emetteur=1,NTr
        Cext_e_V = 0;Cext_e_H = 0;
        Cabs_e_V =0;Cabs_e_H =0;
        
        Allocate(E_ref_incident(3*Nbc_proc,2));
        Call Incident_Field(1,Nbc_proc,Cells_proc,NTr,Transmitters,num_emetteur,num_emetteur,E_ref_incident);
    
        DO I=1,Nbc_proc
            Cabs_e_V = Cabs_e_V + imag(Cells_proc(I)%parameter_Ce)*abs(sum(E_total(3*(I-1)+1:3*I,num_emetteur)))**2.*Cells_proc(I)%Sc**3. ;  
            Cabs_e_H = Cabs_e_H + imag(Cells_proc(I)%parameter_Ce)*abs(sum(E_total(3*(I-1)+1:3*I,num_emetteur+ &
                NTr)))**2.*Cells_proc(I)%Sc**3. ;  
            
            Cext_e_V = Cext_e_V + imag(Cells_proc(I)%parameter_Ce*sum(E_total(3*(I-1)+1:3*I,num_emetteur))&
                *conjg(sum(E_ref_incident(3*(I-1)+1:3*I,1))))*Cells_proc(I)%Sc**3. ;  
            Cext_e_H = Cext_e_H + imag(Cells_proc(I)%parameter_Ce*sum(E_total(3*(I-1)+1:3*I,num_emetteur+NTr))*&
                conjg(sum(E_ref_incident(3*(I-1)+1:3*I,2))))*Cells_proc(I)%Sc**3. ;            
        ENDDO 
        ! pas de 4pi ici car j'ai simplifie par le 4pi de Xi a l'interieur de la somme
        C_ext(num_emetteur) = K_air*(Cext_e_V+Cext_e_H)/2.  
        C_abs(num_emetteur) = K_air*(Cabs_e_V+Cabs_e_H)/2. 
        
        Deallocate(E_ref_incident)
    ENDDO      
    deallocate(Cells_proc);
    
    Allocate(C_ext_all(NTr),C_abs_all(NTr));
        
    Call MPI_BARRIER(MPI_COMM_WORLD ,code);
    Call MPI_ALLREDUCE(C_ext,C_ext_all,NTr,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD ,code);
    Call MPI_ALLREDUCE(C_abs,C_abs_all,NTr,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD ,code);    
    
    C_ext(1:NTr) = C_ext_all(1:NTr);
    C_abs(1:NTr) = C_abs_all(1:NTr);
    deallocate(C_ext_all,C_abs_all);
    
    if (rank == 0) then 
        call date_and_time(date_final_N1,time_final_N1,zone_final_N1,values_final_N1)
        call Calcul_time_spent (values_init_N1,values_final_N1, time_calcul_N1)
        Write (*,'(a,i2,a,i2,a,i2,a,i2,a)') '--> to compute scattered fields, Cext and Cabs from Etot ',time_calcul_N1(1),'j',time_calcul_N1(2)&
        ,'h',time_calcul_N1(3),'min', time_calcul_N1(4),'sec';
    endif  
      
    
END SUBROUTINE Compute_EFields_MoM





SUBROUTINE compute_Green_s_tr_single_element(Cells,ii,jj,Zmom_i_j) 
 
    USE Initialization
    USE common_variables
  
    IMPLICIT NONE
    !! IN/OUT ******************************************************************
    type (Cell), Dimension(Nbc), INTENT(IN) :: Cells
    INTEGER, INTENT(IN) :: ii,jj
    COMPLEX(KIND=8), INTENT(OUT) :: Zmom_i_j
    
    ! local
    INTEGER Is, Isc,Io, Ioc     
    COMPLEX(KIND=8) :: Term1, Term2, Term3
    COMPLEX(KIND=8) :: Fxx, Fyy, Fzz, Fxy, Fyz,Fxz,Fyx,Fzy,Fzx
    REAL(kind=8) :: Distx, Disty, Distz, Dist
    
    !! commencons par connaitre la cellule concernee par la ligne ii, colonne jj  de la matrice MoM 
    Ioc = mod(ii,3)
    If (Ioc .NE. 0) Then
        Io = ii/3 + 1
    Else 
        Io = ii/3
    Endif  
    
    Isc = mod(jj,3)
    If (Isc .NE. 0) Then
        Is = jj/3 + 1
    Else 
        Is = jj/3
    Endif  
    
    
    if (Io .eq. Is)then
        If (Ioc .eq. Isc) Then
            Zmom_i_j = 1 - Cells(Is)%parameter_Sing
        Else
            Zmom_i_j = 0;
        EndIf                 
	else
  
        Distx = Cells(Io)%Xc - Cells(Is)%Xc
        Disty = Cells(Io)%Yc - Cells(Is)%Yc
        Distz = Cells(Io)%Zc - Cells(Is)%Zc
        Dist = sqrt(Distx*Distx+Disty*Disty+Distz*Distz)
    
        Term1 = exp(J*K_air*Dist)/Dist**2.
        Term2 = J*K_air-1/Dist
        Term3 = Cells(Is)%parameter_Const*Cells(Is)%parameter_Ce
      
        If (Isc == 1) Then !! .X
            if (Ioc == 1) Then !! XX
                Fxx = (Term2 + K_air**2.*(Dist - Distx**2./Dist)-3.*Term2*Distx**2./Dist**2.)*Term1;
                Zmom_i_j = - Fxx*Term3; 
            Elseif (Ioc == 2) Then !! YX
                Fxy = (Distx*Disty*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1;
                Zmom_i_j = - Fxy*Term3 
            Else                !! ZX
                Fxz = (Distx*Distz*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1;
                Zmom_i_j = - Fxz*Term3; 
            EndIf
            
        ElseIf (Isc==2) Then !! .Y
            if (Ioc == 1) Then !! XY
                Fxy = (Distx*Disty*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
                Zmom_i_j = - Fxy*Term3;
            Elseif (Ioc == 2) Then !! YY
                Fyy = (Term2 + K_air**2.*(Dist - Disty**2./Dist)-3.*Term2*Disty**2./Dist**2.)*Term1
                Zmom_i_j = - Fyy*Term3;
            Else                !! ZY
                Fyz = (Disty*Distz*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
                Zmom_i_j = - Fyz*Term3;   
            EndIf  
                      
        Else  !! .Z
            if (Ioc == 1) Then !! XZ
                Fxz = (Distx*Distz*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
                Zmom_i_j = -Fxz*Term3;            
            Elseif (Ioc == 2) Then !! YZ
                Fyz = (Disty*Distz*(-K_air**2.-3.*Term2/Dist)/Dist)*Term1
                Zmom_i_j = -Fyz*Term3;             
            Else                !! ZZ
                Fzz = (Term2 + K_air**2.*(Dist - Distz**2./Dist)-3.*Term2*Distz**2./Dist**2.)*Term1
                Zmom_i_j = -Fzz*Term3;            
            EndIf
            
        EndIf
    EndIf 
          

 END SUBROUTINE compute_Green_s_tr_single_element
 
 
 
 
 
 SUBROUTINE Incident_Field_single_element(Cells,Transmitters,ii,kk,E_inc_i_k)

    USE Initialization
    USE common_variables
    USE iso_fortran_env
    Implicit none
    
    !IN/OUT
    type (Cell), Dimension(Nbc), INTENT(IN) :: Cells   
    type (Dipole), Dimension(NTr), INTENT(IN) :: Transmitters
    INTEGER, INTENT(IN) :: ii,kk
    COMPLEX(KIND=8), INTENT(OUT) :: E_inc_i_k
    
    ! local
    Integer :: Ir,num_trans,num_cel,V_or_H
    Complex(kind=8) :: K11x, K11y, K11z, Ex, Ey, Ez
    real(kind=8) :: theta_transmit, phi_transmit, Rx, Ry, Rz

    
    ! num_cel and X/Y/Z
    Ir = mod(ii,3)
    If (Ir .NE. 0) Then
        num_cel = ii/3 + 1
    Else 
        num_cel = ii/3
    Endif  
    
    ! num_trans and V/H polar
    If (kk .le. NTr) Then
        V_or_H = 1; ! Vertical Poloarization
        num_trans = kk;
    Else 
        V_or_H = 2; ! Horizontal Poloarization
        num_trans = kk - NTr; 
    Endif  
    
    theta_transmit = Transmitters(num_trans)%theta
    phi_transmit = Transmitters(num_trans)%phi

    K11x = K_air*cos(theta_transmit*Pi/180.); 
    K11y = K_air*sin(theta_transmit*Pi/180.)*cos(phi_transmit*Pi/180.)
    K11z = K_air*sin(theta_transmit*Pi/180.)*sin(phi_transmit*Pi/180.) 

    Rx = Cells(num_cel)%Xc
    Ry = Cells(num_cel)%Yc
    Rz = Cells(num_cel)%Zc
    
    if (V_or_H .eq. 1) then 

      !!--------------------------------Polarisation Verticale---------------------------------
        If (Ir .eq. 1) Then ! X
            E_inc_i_k = - exp(J*K11x*Rx)*exp(J*K11y*Ry)*exp(J*K11z*Rz)*sin(theta_transmit*Pi/180.); ! Ex  
        Elseif (Ir .eq. 2) Then ! Y
            E_inc_i_k = exp(J*K11x*Rx)*exp(J*K11y*Ry)*exp(J*K11z*Rz)*cos(theta_transmit*Pi/180.)*cos(phi_transmit*Pi/180.); ! Ey
        Else               ! Z
            E_inc_i_k = exp(J*K11x*Rx)*exp(J*K11y*Ry)*exp(J*K11z*Rz)*cos(theta_transmit*Pi/180.)*sin(phi_transmit*Pi/180.); ! Ez
        EndIf    
      
    Elseif (V_or_H .eq. 2) then              

    !!-------------------------------Polarisation Horizontale-------------------------------- 
        If (Ir .eq. 1) Then ! X 
            E_inc_i_k = 0.0;  
        ElseIf (Ir .eq. 2) Then ! Y                                                                      ! Ex
            E_inc_i_k = - exp(J*K11x*Rx)*exp(J*K11y*Ry)*exp(J*K11z*Rz)*sin(phi_transmit*Pi/180.);   ! Ey
        Else
            E_inc_i_k = exp(J*K11x*Rx)*exp(J*K11y*Ry)*exp(J*K11z*Rz)*cos(phi_transmit*Pi/180.);     ! Ez
        EndIf
    Endif   
End Subroutine Incident_Field_single_element



!! USE IF NEEDED TO WRITE ZLOC 

!! character(4) :: rank_str
        
    ! solve system with PSGESV ; the solution X overwrites the RHS B
!    if (rank == 0) then 
!        Write(*,*) '3*Nbc = ',3*Nbc
!        Write(*,*) '2*NTr = ',2*NTr
!        Write(*,*) 'IA = ',IA
!        Write(*,*) 'JA = ',JA
!        Write(*,*) 'DESCA = ',DESCA
!        Write(*,*) 'DESCB = ',DESCB
!        Write(*,*) 'IB = ',IB
!        Write(*,*) 'JB = ',JB
!    endif
!    
!    !! ***********************************************************************
!    !! ***********************************************************************
!    !! ***********************************************************************
!    !! Write ZLoc to compare with MoM matrix form ST MoM  
!    if (rank .lt. 10) Then 
!        Write(rank_str,'(i1)') rank
!    elseif (rank .lt. 100) Then 
!        Write(rank_str,'(i2)') rank
!    elseif (rank .lt. 1000) Then
!        Write(rank_str,'(i3)') rank
!    else
!        Write(rank_str,'(i4)') rank        
!    endif
!    file_name = trim(SimOutfld_name)//Env_sep//'ZLoc_'//trim(rank_str)//'.dat' 
!    Open(unit=60+rank,File = trim(file_name)); 
!    Write(60+rank,'(a,i4,a,i4)')  'Mlocal = ',Mlocal,'; Nlocal = ',Nlocal
!    Do ii=1, Mlocal
!      Do jj=1, Nlocal
!        Write(60+rank,'(e12.4,a,e12.4)')  Real(ZLoc(ii,jj)),'; ',Imag(ZLoc(ii,jj))
!      EndDo
!    EndDo
!    Close(60+rank);      
!    file_name = trim(SimOutfld_name)//Env_sep//'VLoc_'//trim(rank_str)//'.dat' 
!    Open(unit=60+rank,File = trim(file_name)); 
!    Write(60+rank,'(a,i4,a,i4)')  'Mlocal = ',Mlocal,'; NRHSlocal = ',NRHSlocal
!    Do ii=1, Mlocal
!      Do jj=1, NRHSlocal
!        Write(60+rank,'(e12.4,a,e12.4)')  Real(VLoc(ii,jj)),'; ',Imag(VLoc(ii,jj))
!      EndDo
!    EndDo
!    Close(60+rank);      
!    !! ***********************************************************************
!    !! ***********************************************************************
!    !! ***********************************************************************