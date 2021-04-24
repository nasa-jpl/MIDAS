SUBROUTINE Compute_EFields_ST_MoM(Cells,Transmitters,Receivers,S_total,C_ext,C_abs)

    ! Created on 9-2020 to track the error in MPI MoM and to run comet point target simulations 

    USE Initialization
    USE common_variables
    USE iso_fortran_env
    
    USE f95_precision
    USE lapack95
    

    Implicit none
    
    !IN/OUT
    type (Cell), Dimension(Nbc), INTENT(IN) :: Cells
    type (Dipole), Dimension(NTr), INTENT(IN) ::Transmitters
    type (Dipole), Dimension(NRx_tot), INTENT(IN) :: Receivers
    COMPLEX(real64), Dimension(NRx_tot,4*NTr), INTENT(OUT) :: S_total
    COMPLEX(real64), Dimension(NTr), INTENT(OUT) :: C_ext,C_abs
    
    !Local
    Integer :: I,jj,II,K,kk,ll,lig,Ic,ix,iy,iz
    Integer :: NMB,iii,jjj,cur_i_st,cur_i_end,cur_j_st,cur_j_end
    Integer :: taille_block_diffu, NBlocksEmetteurs, num_cel_fichier, num_fichier
    REAL(kind=8) :: RCOND
    CHARACTER(:), allocatable::exp_name
    character(200) :: file_name
    CHARACTER(3), allocatable::number_chars(:)
    character (len=8):: date
    character (len=10):: time
    character (len=5):: zone
    Integer, Dimension(8) :: values
    Integer, Dimension(:),allocatable :: ipiv
    Integer :: alloc_stat
    
    Integer, Dimension(:),allocatable :: curs_blocks_Emetteurs_Etotal
    Integer :: num_emetteur, num_capteur
    Real(kind=8) :: Cext_e_V,Cext_e_H, Cabs_e_V,Cabs_e_H
    Real(kind=8) :: theta_capteur, phi_capteur
    COMPLEX(real64), Dimension(:,:),allocatable::E_ref_incident
    COMPLEX(real64), Dimension(:,:),allocatable:: Mat_Green,Mat_Green_dr
    COMPLEX(real64) :: ff_coef, Vv, Vh, Hv, Hh
    COMPLEX(real64), Dimension(:), allocatable :: ff_coeffs
    COMPLEX(real64), Dimension(3) :: E_v, E_h
    COMPLEX(real64), Dimension(3*Nbc,2*NTr) :: E_total
    
    
    
    
    ! Time performances
    character(8)  :: date_init, date_final
    character(10) :: time_init, time_final
    character(5)  :: zone_init, zone_final    
    Integer,dimension(8) :: values_init, values_final
    Integer, dimension(4):: Comp_time  



    Write (*,*) ''
    Write (*,*) ''
    Write (*,'(a)') '-------------------------------------------------------------------------------'
    Write (*,'(a)') '-----Conventional MoM to Compute the Electric Fields Inside the Scatterer------'
    Write (*,'(a)') '-------------------------------------------------------------------------------'           
    Write (*,*) ''

    Write(10,*) '' 
    Write(10,'(a)') 'Conventional Method of Moments'

    Comp_time = 0
    call date_and_time(date_init,time_init,zone_init,values_init);
        
    !! Incident Field 
    Allocate(E_ref_incident(3*Nbc,2*NTr))
    Allocate(Mat_Green(3*Nbc,3*Nbc))

    Call Incident_Field(1,Nbc,Cells,NTr,Transmitters,1,NTr,E_ref_incident)
    
    !! Green Function
    Call Green_s_tr_total(Cells,Mat_Green) 
   
   
    !! ***********************************************************************
    ! Write Mat_Green (par block de M_B,N_B = 32)
    !file_name = trim(SimOutfld_name)//Env_sep//'MoM_M_B.dat' 
    !Open(unit=61,File = trim(file_name)); 
    !NMB = 3*Nbc/M_B;
    !Do ii=1, NMB
    !    cur_i_st = (ii-1)*M_B+1;
    !    cur_i_end = ii*M_B;
    !    Do jj=1, NMB
    !        cur_j_st = (jj-1)*M_B+1;
    !        cur_j_end = jj*M_B;
    !        Write(61,'(a,i2,a,i2,a)') ' ';
    !        Write(61,'(a,i2,a,i2,a,i2,a,i2,a)') 'Block : ',ii,',',jj,' (',M_B,';',M_B,')';
    !        Do iii=cur_i_st, cur_i_end
    !            Do jjj= cur_j_st,cur_j_end
    !                Write(61,'(e12.4,a,e12.4)')  Real(Mat_Green(iii,jjj)),'; ',Imag(Mat_Green(iii,jjj))
    !            EndDo
    !        EndDo                 
    !    EndDo
    !    Write(61,'(a,i2,a,i2,a)') ' ';
    !    Write(61,'(a,i2,a,i2,a,i2,a,i2,a)') 'Block : ',ii,',',jj,' (',M_B,';',3*Nbc-NMB*M_B,')';
    !    cur_j_st = cur_j_end + 1;
    !    Do iii=cur_i_st, cur_i_end
    !        Do jjj= cur_j_st,3*Nbc
    !            Write(61,'(e12.4,a,e12.4)')  Real(Mat_Green(iii,jjj)),'; ',Imag(Mat_Green(iii,jjj))
    !        EndDo
    !    EndDo 
    !EndDo
    !!! Last row
    !cur_i_st = cur_i_end+1;
    !cur_i_end = 3*Nbc;
    !Do jj=1, NMB
    !    cur_j_st = (jj-1)*M_B+1;
    !    cur_j_end = jj*M_B;
    !    Write(61,'(a,i2,a,i2,a)') ' ';
    !    Write(61,'(a,i2,a,i2,a,i2,a,i2,a)') 'Block : ',ii,',',jj,' (',3*Nbc-NMB*M_B,';',M_B,')';
    !    Do iii=cur_i_st, cur_i_end
    !        Do jjj= cur_j_st,cur_j_end
    !            Write(61,'(e12.4,a,e12.4)')  Real(Mat_Green(iii,jjj)),'; ',Imag(Mat_Green(iii,jjj))
    !        EndDo
    !    EndDo                 
    !EndDo
    !Write(61,'(a,i2,a,i2,a)') ' ';
    !Write(61,'(a,i2,a,i2,a,i2,a,i2,a)') 'Block : ',ii,',',jj,' (',3*Nbc-NMB*M_B,';',3*Nbc-NMB*M_B,')';
    !cur_j_st = cur_j_end + 1;
    !Do iii=cur_i_st, cur_i_end
    !    Do jjj= cur_j_st,3*Nbc
    !        Write(61,'(e12.4,a,e12.4)')  Real(Mat_Green(iii,jjj)),'; ',Imag(Mat_Green(iii,jjj))
    !    EndDo
    !EndDo 
    !Close(61);  
    
    !! Write MoM matrix (3Nbc,3Nbc) as is to check the previous file
    !file_name = trim(SimOutfld_name)//Env_sep//'MoM_tot.dat' 
    !Open(unit=61,File = trim(file_name)); 
    !Do iii=1, 3*Nbc
    !  Do jjj=1, 3*Nbc
    !    Write(61,'(e12.4,a,e12.4)')  Real(Mat_Green(iii,jjj)),'; ',Imag(Mat_Green(iii,jjj))
    !  EndDo
    !EndDo
    !Close(61);
    
    !! ***********************************************************************
    !! ***********************************************************************
    !! ***********************************************************************
    
    
        
    Write(*,*) ''
    Write(*,'(a,i12)') 'Resolution of the original EM problem of size 3*Nbc =',3*Nbc
    Write(*,*) ''

    !! MoM 
    Call gesvx(Mat_Green,E_ref_incident,E_total,RCOND=RCOND)

    Write(*,'(a,e12.3)') 'RCOND of the Green matrix = ', RCOND
    Write(*,*) ''
    Write(*,*) ''

    Write(10,*) ''
    Write(10,'(a,e12.3)') 'RCOND of the Green matrix = ', RCOND
    Write(10,*) ''

    Deallocate(Mat_Green,E_ref_incident)   
    call date_and_time(date_final,time_final,zone_final,values_final)
    call Calcul_time_spent(values_init,values_final,Comp_time)
    
    Write (*, '(a)') '';
    Write (*, '(a)') 'The total time to compute the internal electric field with Single-Task MoM';
    Write (*, '(a,i2,a,i2,a,i2,a,i2,a)')'is ', Comp_time(1),'j',Comp_time(2),'h',Comp_time(3),&
        'min', Comp_time(4),'sec'
    Write (10, '(a,i2,a,i2,a,i2,a,i2,a)') 'The total time to compute the internal electric field &
        &with Single-Task MoM is ',&
        Comp_time(1),'j',Comp_time(2),'h',Comp_time(3),'min', Comp_time(4),'sec'      
        
     
    ! Write Etot inside the scatterer 
    
    If ((save_Eint .eq. 1) .and. (Nbc .le. save_Eint_Nmax)) then        
        file_name = trim(SimOutfld_name)//Env_sep//'Ein.dat';  
        open(unit = 14, file = trim(file_name))
        Do ii=1,3*Nbc        
            Do kk=1,2*NTr                           
                Write(14,'(es16.8,a,es16.8)') real(E_total(ii,kk)),';',imag(E_total(ii,kk)); 
            EndDo
        EndDo
        close(14);            
    endif  
    
    ! Compute scattered fields ************************************************** 
    Comp_time = 0; call date_and_time(date_init,time_init,zone_init,values_init); 
    
    Write (*,*) ''
    Write (*,*) '----------------------- Scattred Fields --------------------------'
    
    
    DO num_capteur =1,NRx_tot
            
        theta_capteur = Receivers(num_capteur)%theta
        phi_capteur = Receivers(num_capteur)%phi
        
        Allocate(ff_coeffs(Nbc))
        Call GetFFieldCoeff(Nbc,Cells,theta_capteur,phi_capteur,ff_coeffs)    
                             
        DO num_emetteur=1,NTr
            E_v = 0
            E_h = 0      
            DO Ic=1,Nbc           
                ff_coef = ff_coeffs(Ic);
                ix = 3*(Ic-1)+1
                iy = 3*(Ic-1)+2
                iz = 3*Ic;
                
                E_v(1)=E_v(1)+ ff_coef*E_total(ix,num_emetteur)
                E_v(2)=E_v(2)+ ff_coef*E_total(iy,num_emetteur)
                E_v(3)=E_v(3)+ ff_coef*E_total(iz,num_emetteur)
            
                E_h(1)=E_h(1)+ ff_coef*E_total(ix,num_emetteur+NTr)
                E_h(2)=E_h(2)+ ff_coef*E_total(iy,num_emetteur+NTr)
                E_h(3)=E_h(3)+ ff_coef*E_total(iz,num_emetteur+NTr)                    
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
            S_total(num_capteur,4*(num_emetteur-1)+1) = Vv
            S_total(num_capteur,4*(num_emetteur-1)+2) = Vh
            S_total(num_capteur,4*(num_emetteur-1)+3) = Hv
            S_total(num_capteur,4*(num_emetteur-1)+4) = Hh               
        ENDDO            
        Deallocate(ff_coeffs);        
    Enddo
    
    
    !Now Compute_ExtAbsCsec_fromIntField       
    DO num_emetteur=1,NTr
        Cext_e_V = 0;Cext_e_H = 0;
        Cabs_e_V =0;Cabs_e_H =0;
        
        Allocate(E_ref_incident(3*Nbc,2));
        Call Incident_Field(1,Nbc,Cells,NTr,Transmitters,num_emetteur,num_emetteur,E_ref_incident);
    
        DO I=1,Nbc
            Cabs_e_V = Cabs_e_V + imag(Cells(I)%Che_n)*abs(sum(E_total(3*(I-1)+1:3*I,num_emetteur)))**2.*Cells(I)%Sc**3. ;  
            Cabs_e_H = Cabs_e_H + imag(Cells(I)%Che_n)*abs(sum(E_total(3*(I-1)+1:3*I,num_emetteur+ &
                NTr)))**2.*Cells(I)%Sc**3. ;  
            
            Cext_e_V = Cext_e_V + imag(Cells(I)%Che_n*sum(E_total(3*(I-1)+1:3*I,num_emetteur))&
                *conjg(sum(E_ref_incident(3*(I-1)+1:3*I,1))))*Cells(I)%Sc**3. ;  
            Cext_e_H = Cext_e_H + imag(Cells(I)%Che_n*sum(E_total(3*(I-1)+1:3*I,num_emetteur+NTr))*&
                conjg(sum(E_ref_incident(3*(I-1)+1:3*I,2))))*Cells(I)%Sc**3. ;            
        ENDDO 
        ! pas de 4pi ici car j'ai simplifie par le 4pi de Xi a l'interieur de la somme
        C_ext(num_emetteur) = k_0*(Cext_e_V+Cext_e_H)/2.  
        C_abs(num_emetteur) = k_0*(Cabs_e_V+Cabs_e_H)/2. 
        
        deallocate(E_ref_incident);
    ENDDO
    
    call date_and_time(date_final,time_final,zone_final,values_final)
    call Calcul_time_spent(values_init,values_final,Comp_time)
    Write (*,*) ''
    Write (*, '(a,i2,a,i2,a,i2,a,i2,a)') 'The calculation time for the scattered field is ',&
    Comp_time(1),'j',Comp_time(2),'h',Comp_time(3),'min', Comp_time(4),'sec'
    Write (*,*) ''; 
    Write (10, '(a,i2,a,i2,a,i2,a,i2,a)') 'The calculation time for the scattered field is ',&
    Comp_time(1),'j',Comp_time(2),'h',Comp_time(3),'min', Comp_time(4),'sec'
    Write (10,*) ''
    Write (10,*) ''  
    
    
    
    
    
End Subroutine Compute_EFields_ST_MoM