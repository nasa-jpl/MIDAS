SUBROUTINE Compute_Scattered_Fields(nom_methode,Cells,E_total,CBFM_Blocks,MPI_CBFM_Blocks,Transmitters,Receivers,Es_total)
    
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
    COMPLEX(real64), Dimension(NRx_tot,4*NTr), INTENT(OUT) :: Es_total
    
    
    ! Local 
    Integer :: ii,K,I,Ic,ix,iy,iz
    Integer :: Lig,id,nthreads,p,d,Nelts,cel_beg,cel_end
    Integer :: kk_job, kk, curs_cel,Nbc_b
    Integer, Dimension(nber_procs) :: all_Nbc_procs        
    character(200) :: file_name
    COMPLEX(real64), Dimension(:), allocatable :: ff_coeffs, Es_capteur, Es_capteur_all
    COMPLEX(real64), Dimension(3) :: E_v, E_h, E_v_p, E_h_p
    COMPLEX(real64) :: ff_coef, Vv, Vh, Hv, Hh
    Integer :: num_emetteur, num_capteur
    Real(kind=8) :: theta_capteur, phi_capteur, Beta,step_beta
    type(Cell), Dimension(:), allocatable :: Cells_proc
    COMPLEX(real64), Dimension(:,:), allocatable :: Green_s_dt_mat
    
    ! to write Es total 
    Integer ::  jj,dd,cc,kkt,kkr,NTr_wr_proc,a,Nths,Nphs
    character(200) :: file_name_s, Esfold_name
    CHARACTER(6) :: ty,kkt_st 
    CHARACTER(:), allocatable:: nom_meth_exact,stFreq,sim_name
    COMPLEX(real64), Dimension(:,:), allocatable :: Es_vv_2d,Es_hh_2d,Es_vh_2d,Es_hv_2d
    COMPLEX(real64), Dimension(:), allocatable :: Es_vv_1d,Es_hh_1d,Es_vh_1d,Es_hv_1d
    Real(kind=8), Dimension(:), allocatable :: Thetas,Phis,RecThetasVals,RecPhisVals
        
  
   
    If (rank == 0) Then 
        Write (*,*) ''
        Write (*,*) '----------------------- Scattered Fields --------------------------'
    EndIf
    
    
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
	    Allocate(Es_capteur(4*NTr),Es_capteur_all(4*NTr));   
        Allocate(Green_s_dt_mat(3,3*Nbc_proc)) 
        !! Dyade de Greene singuliere
        theta_capteur = Receivers(num_capteur)%theta
        phi_capteur = Receivers(num_capteur)%phi        
        Call Green_s_dt(Nbc_proc,Cells_proc,theta_capteur,phi_capteur,Green_s_dt_mat)  
          
        DO num_emetteur=1,NTr
            E_v = 0
            E_h = 0      
            DO I=1,3*Nbc_proc ! This the difference with the serial OpenMP code, Here each process computes S based on its part of E_total            
                           
                E_v(1)=E_v(1)+ Green_s_dt_mat(1,I)*E_total(I,num_emetteur)
                E_v(2)=E_v(2)+ Green_s_dt_mat(2,I)*E_total(I,num_emetteur)
                E_v(3)=E_v(3)+ Green_s_dt_mat(3,I)*E_total(I,num_emetteur)
            
                E_h(1)=E_h(1)+ Green_s_dt_mat(1,I)*E_total(I,num_emetteur+NTr)
                E_h(2)=E_h(2)+ Green_s_dt_mat(2,I)*E_total(I,num_emetteur+NTr)
                E_h(3)=E_h(3)+ Green_s_dt_mat(3,I)*E_total(I,num_emetteur+NTr)                    
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
              !! -----------------------Remplissage du vecteur Es_total ---------------------------!!
              !! ---------------------------------------------------------------------------------!!
              
              Es_capteur(4*(num_emetteur-1)+1)  = Vv; 
              Es_capteur(4*(num_emetteur-1)+2)  = Vh; 
              Es_capteur(4*(num_emetteur-1)+3)  = Hv; 
              Es_capteur(4*(num_emetteur-1)+4)  = Hh; 
        ENDDO      
        Deallocate(Green_s_dt_mat); 
        Nelts = 4*NTr;
        Call MPI_ALLREDUCE(Es_capteur,Es_capteur_all,Nelts,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,code);
        Es_total(num_capteur,1:4*NTr) = Es_capteur_all(1:4*NTr);
        deallocate(Es_capteur,Es_capteur_all);        
    Enddo
    Deallocate(Cells_proc);
    
    
    !!**************************************************************************************************
    ! Once we have the entire Es_total(NRx_tot,4*NTr), we write it (can be improved later using PHDF5) 
    !!**************************************************************************************************
    
    Esfold_name = trim(SimOutfld_name)//Env_sep//'Es_files';
    NTr_wr_proc = (NTr/nber_procs)+1
    
    ! preparation 
    If (nom_methode=='CBFM-E  ') Then
        Allocate(character(6) ::nom_meth_exact)
        nom_meth_exact = trim(nom_methode)
    Elseif ((nom_methode=='MoM     ') .OR. (nom_methode=='RGE     ')) Then 
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
    
    if ((NumIntType_r .eq. 'sd') .OR. (NumIntType_r .eq. 'lb')) then
        Allocate(Es_vv_1d(NRx),Es_hh_1d(NRx));
        Allocate(Es_vh_1d(NRx),Es_hv_1d(NRx));            
    else
        Allocate(Es_vv_2d(NRxTheta,NRxPhi),Es_hh_2d(NRxTheta,NRxPhi))
        Allocate(Es_vh_2d(NRxTheta,NRxPhi),Es_hv_2d(NRxTheta,NRxPhi))
        
        ! Receivers theta and phi vals 
        Allocate(Thetas(NRx),Phis(NRx))
        Allocate(RecThetasVals(NRxTheta),RecPhisVals(NRxPhi))         
            
        Thetas(:) = Receivers(:)%theta;
        Phis(:) = Receivers(:)%phi;
        call Unique1DArray_D(NRx,Nths,Thetas)
        call Unique1DArray_D(NRx,Nphs,Phis)
    
        RecThetasVals= Thetas(1:Nths); RecPhisVals= Phis(1:Nphs);
        deallocate(Thetas,Phis)
    endif
    
    DO dd=1,NTrPhi
        Do cc=1,NTrTheta
            kkt=(dd-1)* NTrTheta + cc;
            
            if ((NumIntType_r .eq. 'sd') .OR. (NumIntType_r .eq. 'lb')) then 
                Do kkr=1, NRx      
                    Es_vv_1d(kkr)= Es_total(kkr,4*(kkt-1)+1); 
                    Es_vh_1d(kkr) = Es_total(kkr,4*(kkt-1)+2); 
                    Es_hv_1d(kkr) = Es_total(kkr,4*(kkt-1)+3); 
                    Es_hh_1d(kkr)= Es_total(kkr,4*(kkt-1)+4);      
                EndDo
                               
                If ((kkt .gt. (rank*NTr_wr_proc)) .and. (kkt .le. (rank+1)*NTr_wr_proc)) then 
                  Write(kkt_st,'(a,i4.4)') 'kt',kkt;
                  file_name_s = trim(Esfold_name)//Env_sep//sim_name//'Esca_'//stFreq//freq_unit//'_'//trim(kkt_st)//'_'//nom_meth_exact//'.dat';
                  Open(unit=21+rank,File = file_name_s)    
                  Write(21+rank,'(a,a)') '      theta       phi       Re(Evv)        Im(Evv)         Re(Evh)       Im(Ehv) ',&
                                  '        Re(Ehv)       Im(Ehv)        Re(Ehh)        Im(Ehh) '
                  Do jj =1, NRxPhi  
                      Do ii =1, NRxTheta 
                          Write(21+rank,'(f9.2,a,f9.2,a,e12.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4)') &
                          Receivers(ii)%theta,';  ',Receivers(ii)%phi,';  ',Real(Es_vv_1d(ii)),';  ',Imag(Es_vv_1d(ii)),';  ',Real(Es_vh_1d(ii)),&
                        ';  ',Imag(Es_vh_1d(ii)),';  ', Real(Es_hv_1d(ii)),';  ',Imag(Es_hv_1d(ii)),';  ',Real(Es_hh_1d(ii)),';  ',Imag(Es_hh_1d(ii))
                      EndDo
                  EndDo
                  Close(21+rank);  
                endif        
      
            else                                    
                Do jj=1, NRxPhi                   
                    Do ii=1, NRxTheta                        
                        kkr = (jj-1)*NRxTheta + ii;                                
                        Es_vv_2d(ii,jj) = Es_total(kkr,4*(kkt-1)+1); 
                        Es_vh_2d(ii,jj) = Es_total(kkr,4*(kkt-1)+2); 
                        Es_hv_2d(ii,jj) = Es_total(kkr,4*(kkt-1)+3); 
                        Es_hh_2d(ii,jj) = Es_total(kkr,4*(kkt-1)+4);    
                    EndDo        
                EndDo      
                
                If ((kkt .gt. (rank*NTr_wr_proc)) .and. (kkt .le. (rank+1)*NTr_wr_proc)) then 
                    Write(kkt_st,'(a,i4.4)') 'kt',kkt;
                    file_name_s = trim(ESfold_name)//Env_sep//sim_name//'Esca_'//stFreq//freq_unit//'_'//trim(kkt_st)//'_'//nom_meth_exact//'.dat';
                    
                    Open(unit=21+rank,File = file_name_s)    
                    Write(21+rank,'(a,a)') '      theta       phi       Re(Evv)        Im(Evv)         Re(Evh)       Im(Ehv) ',&
                                    '        Re(Ehv)       Im(Ehv)        Re(Ehh)        Im(Ehh) '
                    Do jj =1, NRxPhi  
                        Do ii =1, NRxTheta 
                            Write(21+rank,'(f9.2,a,f9.2,a,e12.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4)') &
                            RecThetasVals(ii),';  ',RecPhisVals(jj),';  ',Real(Es_vv_2d(ii,jj)),';  ',Imag(Es_vv_2d(ii,jj)),';  ',Real(Es_vh_2d(ii,jj)),&
                                ';  ',Imag(Es_vh_2d(ii,jj)),';  ', Real(Es_hv_2d(ii,jj)),';  ',Imag(Es_hv_2d(ii,jj)),';  ',&
                                Real(Es_hh_2d(ii,jj)),';  ',Imag(Es_hh_2d(ii,jj))
                        EndDo
                    EndDo
                    Close(21+rank);   
                EndIf
                   
            EndIf
        Enddo             
    EndDo
       
End Subroutine Compute_Scattered_Fields