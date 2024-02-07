SUBROUTINE Compute_EFields_ST_MoM_InWork(nom_methode,Cells,Transmitters,Receivers,S_total,C_ext,C_abs)

    ! Created on 9-2020 to track the error in MPI MoM and to run comet point target simulations

    USE Initialization
    USE common_variables
    USE iso_fortran_env

    USE f95_precision
    USE lapack95


    Implicit none

    !IN/OUT
    character(8), INTENT(IN):: nom_methode
    type (Cell), Dimension(Nbc), INTENT(IN) :: Cells
    type (Dipole), Dimension(NTr), INTENT(IN) ::Transmitters
    type (Dipole), Dimension(NRx_tot), INTENT(IN) :: Receivers
    COMPLEX(real64), Dimension(NRx_tot,4*NTr), INTENT(OUT) :: S_total
    COMPLEX(real64), Dimension(NTr), INTENT(OUT) :: C_ext,C_abs

    !Local
    Integer :: I,jj,II,K,kk,ll,lig,Ic,ix,iy,iz,a,kkt,kkt_abs
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
    Real(kind=8) :: theta_capteur, phi_capteur,Beta,step_beta
    COMPLEX(real64), Dimension(:,:),allocatable::E_ref_incident
    COMPLEX(real64), Dimension(:,:),allocatable:: Mat_Green,Mat_Green_dr
    COMPLEX(real64) :: ff_coef, Vv, Vh, Hv, Hh,Vv_pol, Vh_pol, Hv_pol, Hh_pol
    COMPLEX(real64), Dimension(:), allocatable :: ff_coeffs
    COMPLEX(real64), Dimension(3) :: E_v, E_h, E_v_pol, E_h_pol
    COMPLEX(real64), Dimension(3*Nbc,2*NTr) :: E_total
    ! debug 
    Integer :: num_pol
    COMPLEX(real64), Dimension(:,:), allocatable :: S_total_pol
    CHARACTER(:), allocatable:: nom_meth_exact,stFreq,sim_name
    CHARACTER(200) :: file_name_s,Sfold_name
    CHARACTER(6) :: ty,kkt_st

    ! Time performances
    character(8)  :: date_init, date_final
    character(10) :: time_init, time_final
    character(5)  :: zone_init, zone_final
    Integer,dimension(8) :: values_init, values_final
    Integer, dimension(4):: Comp_time
    
    
    If (nom_methode=='CBFM-E  ') Then
        Allocate(character(6) ::nom_meth_exact)
        nom_meth_exact = trim(nom_methode)
    Elseif ((nom_methode=='MoM     ') .OR. (nom_methode=='RGE     ')) Then
        Allocate(character(3) ::nom_meth_exact)
        nom_meth_exact = trim(nom_methode)
    Endif
    Sfold_name = trim(SimOutfld_name)//Env_sep//'S_files';
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

    if (NPolBeta > 1) then
        step_beta = (beta_final_Pol-beta_init_Pol)/(NPolBeta-1);
    else
        step_beta = 0;
        
    endif

    if (rank == 0) then
      Write (*,*) ''
      Write (*,*) ''
      Write (*,'(a)') '-------------------------------------------------------------------------------'
      Write (*,'(a)') '-----Conventional MoM to Compute the Electric Fields Inside the Scatterer------'
      Write (*,'(a)') '-------------------------------------------------------------------------------'
      Write (*,*) ''

      Write(10,*) ''
      Write(10,'(a)') 'Conventional Method of Moments'
    endif
    Comp_time = 0
    call date_and_time(date_init,time_init,zone_init,values_init);

    !! Incident Field
    Allocate(E_ref_incident(3*Nbc,2*NTr))
    Allocate(Mat_Green(3*Nbc,3*Nbc))

    Call Incident_Field(1,Nbc,Cells,NTr,Transmitters,1,NTr,E_ref_incident)

    !! Green Function
    Call Green_s_tr_total(Cells,Mat_Green)

    if (rank == 0) then
    Write(*,*) ''
    Write(*,'(a,i12)') 'Resolution of the original EM problem of size 3*Nbc =',3*Nbc
    Write(*,*) ''
    endif

    !! MoM
    Call gesvx(Mat_Green,E_ref_incident,E_total,RCOND=RCOND)

    if (rank == 0) then
    Write(*,'(a,e12.3)') 'RCOND of the Green matrix = ', RCOND
    Write(*,*) ''
    Write(*,*) ''
    Write(10,*) ''
    Write(10,'(a,e12.3)') 'RCOND of the Green matrix = ', RCOND
    Write(10,*) ''
    endif

    Deallocate(Mat_Green,E_ref_incident)
    call date_and_time(date_final,time_final,zone_final,values_final)
    call Calcul_time_spent(values_init,values_final,Comp_time)

    if (rank == 0) then
    Write (*, '(a)') '';
    Write (*, '(a)') 'The total time to compute the internal electric field with Single-Task MoM';
    Write (*, '(a,i2,a,i2,a,i2,a,i2,a)')'is ', Comp_time(1),'j',Comp_time(2),'h',Comp_time(3),&
        'min', Comp_time(4),'sec'
    Write (10, '(a,i2,a,i2,a,i2,a,i2,a)') 'The total time to compute the internal electric field &
        &with Single-Task MoM is ',&
        Comp_time(1),'j',Comp_time(2),'h',Comp_time(3),'min', Comp_time(4),'sec'
    endif

    ! Write Etot inside the scatterer

    If ((save_Eint .eq. 1) .and. (Nbc .le. save_Eint_Nmax)) then
        file_name = trim(SimOutfld_name)//Env_sep//'Ein_MoM.dat';
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

    if (rank == 0) then
    Write (*,*) ''
    Write (*,*) '----------------------- Scattered Fields --------------------------'
    endif

    Allocate(S_total_pol(NRx_tot,4*NTr*NPolBeta)); ! Let's start by storing the entire S_total_pol(NRx_tot,4*NTr*NPolBeta), 
                                                                                    ! later we can do something similar to whats done in Compute_scattering_Matrices_InWork for debug
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
            
            Do num_pol = 1,NPolBeta
                Beta =  beta_init_Pol + (num_pol-1)*step_beta
                
                ! E_v_pol and E_h_pol (x, y and z) (the scattered field resulting from a polarized incident wave is the linear sum of the solution for Ev and Eh)
                E_v_pol = cos(Beta*Pi/180.)*E_v + sin(Beta*Pi/180.)*E_h
                E_h_pol = -sin(Beta*Pi/180.)*E_v + cos(Beta*Pi/180.)*E_h

                !! ---------------------------------------------------------------------------------!!
                !! ------------------------------Polarisation Vv -----------------------------------!!
                !! ---------------------------------------------------------------------------------!!
                Vv= - E_v(1)*sin(theta_capteur*Pi/180.)+E_v(2)*cos(theta_capteur*Pi/180.)*cos(phi_capteur*Pi/180.) &
                +E_v(3)*cos(theta_capteur*Pi/180.)*sin(phi_capteur*Pi/180.);
                
                ! Insert Beta
                ! The Vv above is Vvbeta (beta=0), so you can store Vv for Q calculations and write Vvbeta, same for the three other quantities
                Vv_pol = - E_v_pol(1)*cos(Beta*Pi/180.)*sin(theta_capteur*Pi/180.) &
                         + E_v_pol(2)*(cos(Beta*Pi/180.)*cos(theta_capteur*Pi/180.)*cos(phi_capteur*Pi/180.) - sin(Beta*Pi/180.)*sin(phi_capteur*Pi/180.)) &
                         + E_v_pol(3)*(cos(Beta*Pi/180.)*cos(theta_capteur*Pi/180.)*sin(phi_capteur*Pi/180.) + sin(Beta*Pi/180.)*cos(phi_capteur*Pi/180.));  
    
    
                !! ---------------------------------------------------------------------------------!!
                !! ------------------------------Polarisation Vh -----------------------------------!!
                !! ---------------------------------------------------------------------------------!!
                Vh= - E_v(2)*sin(phi_capteur*Pi/180.)+E_v(3)*cos(phi_capteur*Pi/180.);
    
                Vh_pol = + E_v_pol(1)* sin(Beta*Pi/180.)*sin(theta_capteur*Pi/180.) &
                         - E_v_pol(2)* (sin(Beta*Pi/180.)*cos(theta_capteur*Pi/180.)*cos(phi_capteur*Pi/180.)+cos(Beta*Pi/180.)*sin(phi_capteur*Pi/180.))&
                         + E_v_pol(3)* (cos(Beta*Pi/180.)*cos(phi_capteur*Pi/180.) - sin(Beta*Pi/180.)*cos(theta_capteur*Pi/180.)*sin(phi_capteur*Pi/180.))  ;      
        
                !! ---------------------------------------------------------------------------------!!
                !! ------------------------------Polarisation Hv -----------------------------------!!
                !! ---------------------------------------------------------------------------------!!
                Hv= - E_h(1)*sin(theta_capteur*Pi/180.)+E_h(2)*cos(theta_capteur*Pi/180.)*cos(phi_capteur*Pi/180.) &
                + E_h(3)*cos(theta_capteur*Pi/180.)*sin(phi_capteur*Pi/180.);
    
                Hv_pol = - E_h_pol(1)*cos(Beta*Pi/180.)*sin(theta_capteur*Pi/180.) &
                         + E_h_pol(2)*(cos(Beta*Pi/180.)*cos(theta_capteur*Pi/180.)*cos(phi_capteur*Pi/180.) - sin(Beta*Pi/180.)*sin(phi_capteur*Pi/180.)) &
                         + E_h_pol(3)*(cos(Beta*Pi/180.)*cos(theta_capteur*Pi/180.)*sin(phi_capteur*Pi/180.) + sin(Beta*Pi/180.)*cos(phi_capteur*Pi/180.));  
        
                !! ---------------------------------------------------------------------------------!!
                !! ------------------------------Polarisation Hh -----------------------------------!!
                !! ---------------------------------------------------------------------------------!!
                Hh= - E_h(2)*sin(phi_capteur*Pi/180.) + E_h(3)*cos(phi_capteur*Pi/180.);
    
                Hh_pol = + E_h_pol(1)* sin(Beta*Pi/180.)*sin(theta_capteur*Pi/180.) &
                         - E_h_pol(2)* (sin(Beta*Pi/180.)*cos(theta_capteur*Pi/180.)*cos(phi_capteur*Pi/180.)+cos(Beta*Pi/180.)*sin(phi_capteur*Pi/180.))&
                         + E_h_pol(3)* (cos(Beta*Pi/180.)*cos(phi_capteur*Pi/180.) - sin(Beta*Pi/180.)*cos(theta_capteur*Pi/180.)*sin(phi_capteur*Pi/180.))  ;      
         
    
                !! ---------------------------------------------------------------------------------!!
                !! -----------------------Remplissage du vecteur S_total ---------------------------!!
                !! ---------------------------------------------------------------------------------!!
                S_total(num_capteur,4*(num_emetteur-1)+1) = Vv
                S_total(num_capteur,4*(num_emetteur-1)+2) = Vh
                S_total(num_capteur,4*(num_emetteur-1)+3) = Hv
                S_total(num_capteur,4*(num_emetteur-1)+4) = Hh
                
                S_total_pol(num_capteur,4*NPolBeta*(num_emetteur-1)+4*(num_pol-1)+1)  = Vv_pol; 
                S_total_pol(num_capteur,4*NPolBeta*(num_emetteur-1)+4*(num_pol-1)+2)  = Vh_pol; 
                S_total_pol(num_capteur,4*NPolBeta*(num_emetteur-1)+4*(num_pol-1)+3)  = Hv_pol; 
                S_total_pol(num_capteur,4*NPolBeta*(num_emetteur-1)+4*(num_pol-1)+4)  = Hh_pol; 
            ENDDO
                
        ENDDO
        Deallocate(ff_coeffs);
    Enddo
    
    ! Write S_total_pol
    DO kkt=1,NTr 
      If ((wr_Sij .eq. 1)) Then
          Do num_pol = 1,NPolBeta 
              Beta =  beta_init_Pol + (num_pol-1)*step_beta
              kkt_abs = (kkt-1)*NPolBeta+num_pol;
              
              Write(kkt_st,'(a,i4.4)') 'kt',kkt_abs;
              if (EqSph == 0) then
                  file_name_s = trim(Sfold_name)//Env_sep//sim_name//'SmtableInWork_'//stFreq//trim(freq_unit)//'_'//trim(kkt_st)//'_'//nom_meth_exact//'.dat';
              else
                  file_name_s = trim(Sfold_name)//Env_sep//sim_name//'SmtableES_'//stFreq//trim(freq_unit)//'_'//trim(kkt_st)//'_'//nom_meth_exact//'.dat';
              endif
              
              Open(unit=21,File = file_name_s);
              Write(21+rank, '(a,f10.4,a,f10.4,a,f10.4)') 'THETA =',  Transmitters(kkt)%theta, '; PHI =',  Transmitters(kkt)%phi,'; BETA =', Beta
              Write(21+rank,'(a,a)') '      theta       phi      Re(Svv)        Im(Svv)         Re(Svh)       Im(Svh) ',&
                          '        Re(Shv)       Im(Shv)        Re(Shh)        Im(Shh) '
              
              Do num_capteur = 1, NRx                      
                Vv_pol = S_total_pol(num_capteur,4*NPolBeta*(kkt-1)+4*(num_pol-1)+1);
                Vh_pol = S_total_pol(num_capteur,4*NPolBeta*(kkt-1)+4*(num_pol-1)+2); 
                Hv_pol = S_total_pol(num_capteur,4*NPolBeta*(kkt-1)+4*(num_pol-1)+3); 
                Hh_pol = S_total_pol(num_capteur,4*NPolBeta*(kkt-1)+4*(num_pol-1)+4);
                
		theta_capteur = Receivers(num_capteur)%theta
        	phi_capteur = Receivers(num_capteur)%phi
         
                Write(21+rank,'(f10.4,a,f10.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4)') &
                theta_capteur,';  ',phi_capteur,';  ',Real(Vv_pol),';  ',Imag(Vv_pol),';  ',Real(Vh_pol),&
                  ';  ',Imag(Vh_pol),';  ', Real(Hv_pol),';  ',Imag(Hv_pol),';  ',Real(Hh_pol),';  ',Imag(Hh_pol)
              EndDo
              Close(21);
          EndDo
      EndIf
    EndDo
    Deallocate(S_total_pol);

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
    if (rank == 0) then
    Write (*,*) ''
    Write (*, '(a,i2,a,i2,a,i2,a,i2,a)') 'The calculation time for the scattered field is ',&
    Comp_time(1),'j',Comp_time(2),'h',Comp_time(3),'min', Comp_time(4),'sec'
    Write (*,*) '';
    Write (10, '(a,i2,a,i2,a,i2,a,i2,a)') 'The calculation time for the scattered field is ',&
    Comp_time(1),'j',Comp_time(2),'h',Comp_time(3),'min', Comp_time(4),'sec'
    Write (10,*) ''
    Write (10,*) ''
    endif




End Subroutine Compute_EFields_ST_MoM_InWork
