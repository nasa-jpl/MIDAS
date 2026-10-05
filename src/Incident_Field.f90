SUBROUTINE Incident_Field(cel_init,size_Cells,Cells_in,Nb_transmitters,Transmitters,num_tr_start, num_tr_end,E_ref_incident)

    USE Initialization
    USE common_variables
    USE iso_fortran_env
    Implicit none
    
    !IN/OUT
    Integer, INTENT(IN) :: cel_init,size_Cells
    type (Cell), Dimension(size_Cells), INTENT(IN) :: Cells_in
    Integer, INTENT(IN) :: Nb_transmitters ! attention here could be the commen variable NTr or Nipws when calculating the CBFs 
    type (Dipole), Dimension(Nb_transmitters), INTENT(IN) :: Transmitters
    Integer, INTENT(IN) :: num_tr_start, num_tr_end
    COMPLEX(real64), Dimension(3*size_Cells,2*(num_tr_end-num_tr_start+1)), INTENT(OUT)::E_ref_incident

    ! local
    Complex :: K11x, K11y, K11z, Ex, Ey, Ez
    real(kind=8) :: theta_transmit, phi_transmit, Rx, Ry, Rz
    Integer :: NcalcTr,num_trans,num_Eref_v, num_Eref_h,num_cel, sol,curs_cel

    E_ref_incident = 0
    NcalcTr = num_tr_end-num_tr_start+1;

    DO num_trans=num_tr_start, num_tr_end
        num_Eref_v = num_trans - num_tr_start + 1; 
        num_Eref_h = num_Eref_v + NcalcTr;
        
        theta_transmit = Transmitters(num_trans)%theta
        phi_transmit = Transmitters(num_trans)%phi
  
	    ! The incident wave is propagating in positive x (if theta_i=phi_i=0)
        K11x = k_0*cos(theta_transmit*Pi/180.); 
        K11y = k_0*sin(theta_transmit*Pi/180.)*cos(phi_transmit*Pi/180.)
        K11z = k_0*sin(theta_transmit*Pi/180.)*sin(phi_transmit*Pi/180.) 
  
        curs_cel = 1
        DO num_cel = cel_init, cel_init+size_Cells-1
            Rx = Cells_in(num_cel)%Xc
            Ry = Cells_in(num_cel)%Yc
            Rz = Cells_in(num_cel)%Zc

            !!--------------------------------Polarisation Verticale---------------------------------
            Ex = - exp(J*K11x*Rx)*exp(J*K11y*Ry)*exp(J*K11z*Rz)*sin(theta_transmit*Pi/180.)  
            Ey = exp(J*K11x*Rx)*exp(J*K11y*Ry)*exp(J*K11z*Rz)*cos(theta_transmit*Pi/180.)*cos(phi_transmit*Pi/180.)
            Ez = exp(J*K11x*Rx)*exp(J*K11y*Ry)*exp(J*K11z*Rz)*cos(theta_transmit*Pi/180.)*sin(phi_transmit*Pi/180.)
    
	        E_ref_incident(3*(curs_cel-1)+1,num_Eref_v)= Ex
            E_ref_incident(3*(curs_cel-1)+2,num_Eref_v)= Ey
            E_ref_incident(3*(curs_cel-1)+3,num_Eref_v)= Ez 

            !!-------------------------------Polarisation Horizontale--------------------------------  
            Ex = 0.0
            Ey = - exp(J*K11x*Rx)*exp(J*K11y*Ry)*exp(J*K11z*Rz)*sin(phi_transmit*Pi/180.)
            Ez = exp(J*K11x*Rx)*exp(J*K11y*Ry)*exp(J*K11z*Rz)*cos(phi_transmit*Pi/180.)

            E_ref_incident(3*(curs_cel-1)+1, num_Eref_h)= Ex
            E_ref_incident(3*(curs_cel-1)+2, num_Eref_h)= Ey
            E_ref_incident(3*(curs_cel-1)+3, num_Eref_h)= Ez 
    
            curs_cel = curs_cel + 1
        Enddo  
    Enddo 
End Subroutine Incident_Field
    
SUBROUTINE Incident_Field_Spherical(cel_init,size_Cells,Cells_in,Nb_transmitters,Transmitters,num_tr_start, num_tr_end,E_ref_incident)

    USE Initialization
    USE common_variables
    USE iso_fortran_env
    Implicit none
    
    !IN/OUT
    Integer, INTENT(IN) :: cel_init,size_Cells
    type (Cell), Dimension(size_Cells), INTENT(IN) :: Cells_in
    Integer, INTENT(IN) :: Nb_transmitters ! attention here could be the commen variable NTr or Nipws when calculating the CBFs 
    type (Dipole), Dimension(Nb_transmitters), INTENT(IN) :: Transmitters
    Integer, INTENT(IN) :: num_tr_start, num_tr_end
    COMPLEX(real64), Dimension(3*size_Cells,2*(num_tr_end-num_tr_start+1)), INTENT(OUT)::E_ref_incident

    ! local
    Complex :: K11x, K11y, K11z, Ex, Ey, Ez
    real(kind=8) :: theta_transmit, phi_transmit, Rx, Ry, Rz, Rbl
    real(kind=8) :: Dp, Rsource, Xs, Ys, Zs, Dx, Dy, Dz, R, ux, uy, uz, theta, phi
    real(kind=8) :: Cbl(3)
    Complex :: phase
    Integer :: NcalcTr,num_trans,num_Eref_v, num_Eref_h,num_cel, sol,curs_cel
    logical, save :: Rsource_printed = .false.

    ! bounding sphere of the block: centroid + max distance (within ~2x of the minimum
    ! enclosing sphere, which is all Rsource needs; no recursion, no N x 3 stack array)
    Cbl(1) = sum(Cells_in(:)%Xc)/size_Cells
    Cbl(2) = sum(Cells_in(:)%Yc)/size_Cells
    Cbl(3) = sum(Cells_in(:)%Zc)/size_Cells
    Rbl = sqrt(maxval((Cells_in(:)%Xc-Cbl(1))**2 + (Cells_in(:)%Yc-Cbl(2))**2 + (Cells_in(:)%Zc-Cbl(3))**2))

    ! source distance from the block centre: 10 block radii (at least 10 wavelengths)
    Rsource = 10.d0*max(Rbl, 2.d0*Pi/k_0)
    if ((rank .eq. 0) .and. (.not. Rsource_printed)) then
        write(*,'(a,f10.4,a)') ' -- > spherical IWs: Rsource = 10 x block radius (first block: ',Rsource*1e3,' mm)'
        Rsource_printed = .true.
    endif
    E_ref_incident = 0
    NcalcTr = num_tr_end-num_tr_start+1;

    DO num_trans=num_tr_start, num_tr_end
        num_Eref_v = num_trans - num_tr_start + 1; 
        num_Eref_h = num_Eref_v + NcalcTr;
        
        theta_transmit = Transmitters(num_trans)%theta
        phi_transmit = Transmitters(num_trans)%phi
  
        ! source on a sphere of radius Rsource around the block centre Cbl, in direction (theta, phi)
        Xs = Cbl(1) + Rsource*cos(theta_transmit*Pi/180.);
        Ys = Cbl(2) + Rsource*sin(theta_transmit*Pi/180.)*cos(phi_transmit*Pi/180.)
        Zs = Cbl(3) + Rsource*sin(theta_transmit*Pi/180.)*sin(phi_transmit*Pi/180.)
  
        curs_cel = 1
        DO num_cel = cel_init, cel_init+size_Cells-1
            Rx = Cells_in(num_cel)%Xc
            Ry = Cells_in(num_cel)%Yc
            Rz = Cells_in(num_cel)%Zc
            
            ! Vector from source to observation point
            Dx = Rx - Xs
            Dy = Ry - Ys
            Dz = Rz - Zs

            R = sqrt(Dx*Dx + Dy*Dy + Dz*Dz)

            ! I shouldn't have this issue for now, as I am choosing Rsource >> particle size (I need to reinforce it later)
            ! IF (R .EQ. 0.0) CYCLE
            
            ! Local radial unit vector
            ux = Dx / R
            uy = Dy / R
            uz = Dz / R

            ! Spherical wave phase + amplitude
            phase = exp(J * k_0 * R) / R

            ! Local angles (for polarization definition)
            theta = acos(ux)
            phi   = atan2(uz, uy)

            !!--------------------------------Polarisation Verticale---------------------------------
            Ex = - phase * sin(theta)  
            Ey = phase * cos(theta)*cos(phi)
            Ez = phase * cos(theta)*sin(phi)
    
	        E_ref_incident(3*(curs_cel-1)+1,num_Eref_v)= Ex
            E_ref_incident(3*(curs_cel-1)+2,num_Eref_v)= Ey
            E_ref_incident(3*(curs_cel-1)+3,num_Eref_v)= Ez 

            !!-------------------------------Polarisation Horizontale--------------------------------  
            Ex = 0.0
            Ey = - phase * sin(phi)
            Ez = phase * cos(phi)

            E_ref_incident(3*(curs_cel-1)+1, num_Eref_h)= Ex
            E_ref_incident(3*(curs_cel-1)+2, num_Eref_h)= Ey
            E_ref_incident(3*(curs_cel-1)+3, num_Eref_h)= Ez 
    
            curs_cel = curs_cel + 1
        Enddo  
    Enddo 
End Subroutine Incident_Field_Spherical  

SUBROUTINE Incident_Field_at_Rx(nom_methode,Transmitters,Receivers,E_incident_at_Rx)
    
    USE Initialization
    USE common_variables
    USE iso_fortran_env
    Implicit none
    
    !IN/OUT
    character(8), INTENT(IN):: nom_methode
    type (Dipole), Dimension(NTr), INTENT(IN) :: Transmitters
    type (Dipole), Dimension(NRx_tot), INTENT(IN) :: Receivers
    COMPLEX(real64), Dimension(NRx_tot,4*NTr), INTENT(OUT) :: E_incident_at_Rx
    

    ! local
    Integer :: NTr_wr_proc,num_trans,num_capteur
    real(kind=8) :: theta_transmit, phi_transmit
    real(kind=8) :: theta_capteur, phi_capteur, Rx, Ry, Rz
    Complex :: K11x, K11y, K11z, Ex_v, Ey_v, Ez_v, Ex_h, Ey_h, Ez_h
    Complex :: Vv, Vh, Hv, Hh
    
    
    ! to write Einc 
    Integer ::  ii,jj,dd,cc,kkt,kkr, Nths,Nphs,a
    character(200) :: file_name_s, Eifold_name
    CHARACTER(6) :: ty,kkt_st
    CHARACTER(:), allocatable:: nom_meth_exact,stFreq,sim_name
    COMPLEX(real64), Dimension(:,:), allocatable :: Ei_vv_2d,Ei_hh_2d,Ei_vh_2d,Ei_hv_2d
    COMPLEX(real64), Dimension(:), allocatable :: Ei_vv_1d,Ei_hh_1d,Ei_vh_1d,Ei_hv_1d
    Real(kind=8), Dimension(:), allocatable :: Thetas,Phis,RecThetasVals,RecPhisVals

    NTr_wr_proc = (NTr/nber_procs)+1
    E_incident_at_Rx = 0
    
    DO num_trans=1, NTr  
        If ((num_trans .gt. (rank*NTr_wr_proc)) .and. (num_trans .le. (rank+1)*NTr_wr_proc)) then    
          theta_transmit = Transmitters(num_trans)%theta
          phi_transmit = Transmitters(num_trans)%phi
    
  	   ! The incident wave is propagating in positive x (if theta_i=phi_i=0)
          K11x = k_0*cos(theta_transmit*Pi/180.); 
          K11y = k_0*sin(theta_transmit*Pi/180.)*cos(phi_transmit*Pi/180.)
          K11z = k_0*sin(theta_transmit*Pi/180.)*sin(phi_transmit*Pi/180.)   
          
          DO num_capteur= 1,NRx_tot
              theta_capteur = Receivers(num_capteur)%theta
              phi_capteur = Receivers(num_capteur)%phi
          
              Rx = Rso*cos(theta_capteur*Pi/180.) 
              Ry = Rso*sin(theta_capteur*Pi/180.)*cos(phi_capteur*Pi/180.)
              Rz = Rso*sin(theta_capteur*Pi/180.)*sin(phi_capteur*Pi/180.)
      
              !!--------------------------------Polarisation Verticale---------------------------------
              Ex_v = - exp(J*K11x*Rx)*exp(J*K11y*Ry)*exp(J*K11z*Rz)*sin(theta_transmit*Pi/180.)  
              Ey_v = exp(J*K11x*Rx)*exp(J*K11y*Ry)*exp(J*K11z*Rz)*cos(theta_transmit*Pi/180.)*cos(phi_transmit*Pi/180.)
              Ez_v = exp(J*K11x*Rx)*exp(J*K11y*Ry)*exp(J*K11z*Rz)*cos(theta_transmit*Pi/180.)*sin(phi_transmit*Pi/180.)
      
              !!-------------------------------Polarisation Horizontale--------------------------------  
              Ex_h = 0.0
              Ey_h = - exp(J*K11x*Rx)*exp(J*K11y*Ry)*exp(J*K11z*Rz)*sin(phi_transmit*Pi/180.)
              Ez_h = exp(J*K11x*Rx)*exp(J*K11y*Ry)*exp(J*K11z*Rz)*cos(phi_transmit*Pi/180.)            
              
              Vv= - Ex_v*sin(theta_capteur*Pi/180.)+ Ey_v*cos(theta_capteur*Pi/180.)*cos(phi_capteur*Pi/180.) &
                +Ez_v*cos(theta_capteur*Pi/180.)*sin(phi_capteur*Pi/180.)
              Vh= - Ey_v*sin(phi_capteur*Pi/180.)+Ez_v*cos(phi_capteur*Pi/180.);
              Hv= - Ex_h*sin(theta_capteur*Pi/180.)+Ey_h*cos(theta_capteur*Pi/180.)*cos(phi_capteur*Pi/180.) &
                + Ez_h*cos(theta_capteur*Pi/180.)*sin(phi_capteur*Pi/180.);
              Hh= - Ey_h*sin(phi_capteur*Pi/180.) + Ez_h*cos(phi_capteur*Pi/180.);

              E_incident_at_Rx(num_capteur,4*(num_trans-1)+1)= Vv
              E_incident_at_Rx(num_capteur,4*(num_trans-1)+2)= Vh
              E_incident_at_Rx(num_capteur,4*(num_trans-1)+3)= Hv
              E_incident_at_Rx(num_capteur,4*(num_trans-1)+4)= Hh
            Enddo 
        EndIf          
    Enddo
    
    !! ***********************************************************************************
    ! Write E_incident_at_Rx
    !! ***********************************************************************************
    
    Eifold_name = trim(SimOutfld_name)//Env_sep//'Ei_files';
    
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
        Allocate(Ei_vv_1d(NRx),Ei_hh_1d(NRx));
        Allocate(Ei_vh_1d(NRx),Ei_hv_1d(NRx));            
    else
        Allocate(Ei_vv_2d(NRxTheta,NRxPhi),Ei_hh_2d(NRxTheta,NRxPhi))
        Allocate(Ei_vh_2d(NRxTheta,NRxPhi),Ei_hv_2d(NRxTheta,NRxPhi))
        
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
                    Ei_vv_1d(kkr)= E_incident_at_Rx(kkr,4*(kkt-1)+1); 
                    Ei_vh_1d(kkr) = E_incident_at_Rx(kkr,4*(kkt-1)+2); 
                    Ei_hv_1d(kkr) = E_incident_at_Rx(kkr,4*(kkt-1)+3); 
                    Ei_hh_1d(kkr)= E_incident_at_Rx(kkr,4*(kkt-1)+4);      
                EndDo
                               
                If ((kkt .gt. (rank*NTr_wr_proc)) .and. (kkt .le. (rank+1)*NTr_wr_proc)) then 
                  Write(kkt_st,'(a,i4.4)') 'kt',kkt;
                  file_name_s = trim(Eifold_name)//Env_sep//sim_name//'Einc_'//stFreq//freq_unit//'_'//trim(kkt_st)//'_'//nom_meth_exact//'.dat';
                  Open(unit=21+rank,File = file_name_s)    
                  Write(21+rank,'(a,a)') '      theta       phi       Re(Evv)        Im(Evv)         Re(Evh)       Im(Ehv) ',&
                                  '        Re(Ehv)       Im(Ehv)        Re(Ehh)        Im(Ehh) '
                  Do jj =1, NRxPhi  
                      Do ii =1, NRxTheta 
                          Write(21+rank,'(f9.2,a,f9.2,a,e12.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4)') &
                          Receivers(ii)%theta,';  ',Receivers(ii)%phi,';  ',Real(Ei_vv_1d(ii)),';  ',Imag(Ei_vv_1d(ii)),';  ',Real(Ei_vh_1d(ii)),&
                        ';  ',Imag(Ei_vh_1d(ii)),';  ', Real(Ei_hv_1d(ii)),';  ',Imag(Ei_hv_1d(ii)),';  ',Real(Ei_hh_1d(ii)),';  ',Imag(Ei_hh_1d(ii))
                      EndDo
                  EndDo
                  Close(21+rank);  
                endif        
      
            else                                    
                Do jj=1, NRxPhi                   
                    Do ii=1, NRxTheta                        
                        kkr = (jj-1)*NRxTheta + ii;                                
                        Ei_vv_2d(ii,jj) = E_incident_at_Rx(kkr,4*(kkt-1)+1); 
                        Ei_vh_2d(ii,jj) = E_incident_at_Rx(kkr,4*(kkt-1)+2); 
                        Ei_hv_2d(ii,jj) = E_incident_at_Rx(kkr,4*(kkt-1)+3); 
                        Ei_hh_2d(ii,jj) = E_incident_at_Rx(kkr,4*(kkt-1)+4);    
                    EndDo        
                EndDo      
                
                If ((kkt .gt. (rank*NTr_wr_proc)) .and. (kkt .le. (rank+1)*NTr_wr_proc)) then 
                    Write(kkt_st,'(a,i4.4)') 'kt',kkt;
                    file_name_s = trim(Eifold_name)//Env_sep//sim_name//'Einc_'//stFreq//freq_unit//'_'//trim(kkt_st)//'_'//nom_meth_exact//'.dat';
                    
                    Open(unit=21+rank,File = file_name_s)    
                    Write(21+rank,'(a,a)') '      theta       phi       Re(Evv)        Im(Evv)         Re(Evh)       Im(Ehv) ',&
                                    '        Re(Ehv)       Im(Ehv)        Re(Ehh)        Im(Ehh) '
                    Do jj =1, NRxPhi  
                        Do ii =1, NRxTheta 
                            Write(21+rank,'(f9.2,a,f9.2,a,e12.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4)') &
                            RecThetasVals(ii),';  ',RecPhisVals(jj),';  ',Real(Ei_vv_2d(ii,jj)),';  ',Imag(Ei_vv_2d(ii,jj)),';  ',Real(Ei_vh_2d(ii,jj)),&
                                ';  ',Imag(Ei_vh_2d(ii,jj)),';  ', Real(Ei_hv_2d(ii,jj)),';  ',Imag(Ei_hv_2d(ii,jj)),';  ',&
                                Real(Ei_hh_2d(ii,jj)),';  ',Imag(Ei_hh_2d(ii,jj))
                        EndDo
                    EndDo
                    Close(21+rank);   
                EndIf
                   
            EndIf
        Enddo             
    EndDo 
     

END SUBROUTINE Incident_Field_at_Rx
