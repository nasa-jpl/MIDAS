 SUBROUTINE Compute_Scattering_Matrices_InWork(nom_methode,Cells,E_total,CBFM_Blocks,MPI_CBFM_Blocks,Transmitters,Receivers,S_total_out)
    
    USE Initialization
    USE common_variables
    USE iso_fortran_env
    USE MPI

    IMPLICIT NONE

    !IN/OUT
    character(8), INTENT(IN):: nom_methode
    type (Cell), Dimension(Nbc), INTENT(IN) :: Cells
    COMPLEX(real64), Dimension(3*Nbc_proc,2*NTr), INTENT(IN):: E_total
    ! we will need to know the blocks and MPI blocks repartition to have the good corresponding Cells_proc
    type (CBFM_Block), Dimension(Nblocks), INTENT(IN) :: CBFM_Blocks
    Integer, Dimension(nber_procs,Nblk_proc_max+1), INTENT(IN):: MPI_CBFM_Blocks
    type (Dipole), Dimension(NTr), INTENT(IN) :: Transmitters
    type (Dipole), Dimension(NRx_tot), INTENT(IN) :: Receivers
    !COMPLEX(real64), Dimension(NPolBeta*NRx_tot,4*NTr), INTENT(OUT) :: S_total     ! change 11/10/2021  ! originally Dimension(NRx_tot,4*NTr)
    COMPLEX(real64), Dimension(NRx_tot,4*NTr), INTENT(OUT) :: S_total_out
    ! Local 
    Integer :: ii,K,I,Ic,ix,iy,iz,num_pol,kkt,kkt_abs,a
    Integer :: Lig,id,nthreads,p,d,Nelts,cel_beg,cel_end
    Integer :: kk_job, kk, curs_cel,Nbc_b,NTr_wr_proc,NTr_WR_tot_
    Integer, Dimension(nber_procs) :: all_Nbc_procs        
    character(200) :: file_name
    COMPLEX(real64), Dimension(:,:), allocatable :: S_total_all
    COMPLEX(real64), Dimension(:), allocatable :: ff_coeffs, S_total_capteur, S_total_capteur_all
    COMPLEX(real64), Dimension(3) :: E_v, E_h, E_v_pol, E_h_pol
    COMPLEX(real64) :: ff_coef, Vv, Vh, Hv, Hh, Vv_pol, Vh_pol, Hv_pol, Hh_pol
    Integer :: num_emetteur, num_capteur
    Real(kind=8) :: theta_capteur, phi_capteur, Beta,step_beta
    type(Cell), Dimension(:), allocatable :: Cells_proc
    CHARACTER(:), allocatable:: nom_meth_exact,stFreq,sim_name
    CHARACTER(200) :: file_name_s,Sfold_name
    CHARACTER(6) :: ty,kkt_st
  
   
    If (rank == 0) Then 
        Write (*,*) ''
        Write (*,*) '----------------------- Scattering Matrices --------------------------'
    EndIf
    
   
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
    ! All procs recover again this important information 
    call MPI_ALLGATHER (Nbc_proc,1,MPI_INTEGER,all_Nbc_procs,1,MPI_INTEGER,MPI_COMM_WORLD,code);    
    
    ! now every proc focus on its cells/part of E_tot 
    Allocate(Cells_proc(Nbc_proc));
    MyNBlocks = MPI_CBFM_Blocks(rank+1,1);
    curs_cel = 1;
    Do kk_job = 1,MyNBlocks
        kk = MPI_CBFM_Blocks(rank+1,1+kk_job);
        Nbc_b = CBFM_Blocks(kk)%Nbc_b
        cel_beg = sum(CBFM_Blocks(1:kk-1)%Nbc_b)+1;
        cel_end = sum(CBFM_Blocks(1:kk)%Nbc_b);
        Cells_proc(curs_cel:curs_cel+Nbc_b-1) = Cells(cel_beg:cel_end);
        curs_cel = curs_cel + Nbc_b; 
    EndDo
        
    DO num_capteur =1,NRx_tot  
	    Allocate(S_total_capteur(4*NTr*NPolBeta),S_total_capteur_all(4*NTr*NPolBeta));   
        Allocate(ff_coeffs(Nbc_proc)); 
        !! Dyade de Greene singuliere
        theta_capteur = Receivers(num_capteur)%theta
        phi_capteur = Receivers(num_capteur)%phi
        
        Call GetFFieldCoeff(Nbc_proc,Cells_proc,theta_capteur,phi_capteur,ff_coeffs)          
        DO num_emetteur=1,NTr
            E_v = 0
            E_h = 0      
            DO Ic=1,Nbc_proc ! This the difference with the serial OpenMP code, Here each process computes S based on its part of E_total            
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
                
                ! remember in the scattering matrix you are linking Ei and Es ! whatever their polarization is !
                
                ! Look at page 3 of convention code particle : you will have to update the coordinates of Evs and Ehs, then calculating the 4 components of the S matrix will be quite easy
              
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
                
                S_total_capteur(4*NPolBeta*(num_emetteur-1)+4*(num_pol-1)+1)  = Vv_pol; 
                S_total_capteur(4*NPolBeta*(num_emetteur-1)+4*(num_pol-1)+2)  = Vh_pol; 
                S_total_capteur(4*NPolBeta*(num_emetteur-1)+4*(num_pol-1)+3)  = Hv_pol; 
                S_total_capteur(4*NPolBeta*(num_emetteur-1)+4*(num_pol-1)+4)  = Hh_pol; 

            ENDDO
        ENDDO      
        Deallocate(ff_coeffs); 
        Nelts = 4*NPolBeta*NTr;
        Call MPI_ALLREDUCE(S_total_capteur,S_total_capteur_all,Nelts,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,code);
        
        ! for each num_emetteur 1) save to send to  Compute_Scattering_Quantities and 2) write with Beta rotation before deleting 
        DO kkt=1,NTr 
            ! first get from S_total_capteur_all the Smatrix values that will be sent to Compute_Scattering_Quantities for the calculation of efficiency factors (4 values at Beta = 0 for each Tx = 1:NTr)  
            ! keep this size juste le remplissage sera different 
            S_total_out(num_capteur,4*(kkt-1)+1:4*kkt) = S_total_capteur_all(4*NPolBeta*(kkt-1)+1:4*NPolBeta*(kkt-1)+4);    !this corresponds to the first 4 elements (num_pol=1) calculated for each Transmitter 
        
        
            !! New strategy 7/9/2022 (needs enhancement when we will output to hf5 files )  : write here on the fly before deleting S_total_capteur_all (not in Compute_Scattering_Quantities)
            !! *******************************************************************************************************************************************
            ! The idea is to write on the fly 
            If ((wr_Sij .eq. 1) .and. (num_capteur .le. NRx)) Then
                  NTr_WR_tot_ = NTr*NPolBeta;
                  !NTr_wr_proc = (NTr/nber_procs)+1;
                  NTr_wr_proc = (NTr_WR_tot_/nber_procs)+1;
                  Do num_pol = 1,NPolBeta 
                      Beta =  beta_init_Pol + (num_pol-1)*step_beta
                      kkt_abs = (kkt-1)*NPolBeta+num_pol;
                      If ((kkt_abs .gt. (rank*NTr_wr_proc)) .and. (kkt_abs .le. (rank+1)*NTr_wr_proc)) then
                                            
                            Vv_pol = S_total_capteur_all(4*NPolBeta*(kkt-1)+4*(num_pol-1)+1);
                            Vh_pol = S_total_capteur_all(4*NPolBeta*(kkt-1)+4*(num_pol-1)+2); 
                            Hv_pol = S_total_capteur_all(4*NPolBeta*(kkt-1)+4*(num_pol-1)+3); 
                            Hh_pol = S_total_capteur_all(4*NPolBeta*(kkt-1)+4*(num_pol-1)+4);
                                                                                                          
                            Write(kkt_st,'(a,i4.4)') 'kt',kkt_abs;
                            if (EqSph == 0) then
                                file_name_s = trim(Sfold_name)//Env_sep//sim_name//'SmtableInWork_'//stFreq//trim(freq_unit)//'_'//trim(kkt_st)//'_'//nom_meth_exact//'.dat';
                            else
                                file_name_s = trim(Sfold_name)//Env_sep//sim_name//'SmtableES_'//stFreq//trim(freq_unit)//'_'//trim(kkt_st)//'_'//nom_meth_exact//'.dat';
                            endif
                            if (num_capteur .gt. 1) then
                                Open(unit=21+rank,File = file_name_s, Access='Append', Status='old');
                            else
                                Open(unit=21+rank,File = file_name_s);
                                Write(21+rank, '(a,f10.4,a,f10.4,a,f10.4)') 'THETA =',  Transmitters(kkt)%theta, '; PHI =',  Transmitters(kkt)%phi,'; BETA =', Beta
                                Write(21+rank,'(a,a)') '      theta       phi      Re(Svv)        Im(Svv)         Re(Svh)       Im(Svh) ',&
                                            '        Re(Shv)       Im(Shv)        Re(Shh)        Im(Shh) '
                            endif
                                
                            Write(21+rank,'(f10.4,a,f10.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4)') &
                            theta_capteur,';  ',phi_capteur,';  ',Real(Vv_pol),';  ',Imag(Vv_pol),';  ',Real(Vh_pol),&
                              ';  ',Imag(Vh_pol),';  ', Real(Hv_pol),';  ',Imag(Hv_pol),';  ',Real(Hh_pol),';  ',Imag(Hh_pol)
                            
                            Close(21+rank);
                      Endif
                  EndDo
              EndIf
        EndDo
        deallocate(S_total_capteur,S_total_capteur_all);        
    Enddo
    Deallocate(Cells_proc);
    
    ! I was having insufficient memory problem (buffer) with this approach (calculate all S_total and then in 1 communicate calculate the sum)
    !Allocate(S_total_all(NRx_tot,4*NTr));    
    !Call MPI_BARRIER(MPI_COMM_WORLD,code);
    !Nelts = 4*NRx_tot*NTr;
    !Call MPI_ALLREDUCE(S_total,S_total_all,Nelts,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,code);
        
    !S_total(1:NRx_tot,1:4*NTr) = S_total_all(1:NRx_tot,1:4*NTr);
    !deallocate(S_total_all);  
    
End Subroutine Compute_Scattering_Matrices_InWork