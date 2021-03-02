    SUBROUTINE Compute_Scattering_Quantities_1(nom_methode,SimScatterer,Transmitters,Receivers,S_total,C_ext,C_abs)
    
        USE Initialization
        USE common_variables
        USE iso_fortran_env
        USE MPI

        !! "Implicit Statement" 
        IMPLICIT NONE

        !IN/OUT
        character(8), INTENT(IN):: nom_methode
        type (Scatterer), INTENT(IN) :: SimScatterer
        type (Dipole), Dimension(NTr), INTENT(IN) :: Transmitters
        type (Dipole), Dimension(NRx_tot), INTENT(IN) :: Receivers
        COMPLEX(real64), Dimension(NRx_tot,4*NTr), INTENT(IN):: S_total
        COMPLEX(real64), Dimension(NTr), INTENT(IN):: C_ext,C_abs

        ! Local 
        !Integer :: param_name_length,taille_block_diffu, NBlocksEmetteurs
        !Integer, Dimension(:),allocatable :: curs_blocks_Emetteurs_Etotal
        Integer :: ii,jj,cc,dd,kkr,kkt,a,Nang_tot,Nths,Nphs,num_trans,Npts
        Integer :: ind_fwd,ind_bkw,curs,tdistr_sca,Step_in_cosTh
        Integer :: NTr_wr_proc
        
        CHARACTER(:), allocatable:: nom_meth_exact,stFreq,sim_name
        CHARACTER(200) :: file_name_s,file_name_q
        CHARACTER(200) :: Sfold_name,Qfold_name
        CHARACTER(6) :: ty
        CHARACTER(6) :: kkt_st        
        CHARACTER(11) :: Qextintfst,Qabsintfst
        
        Real(kind=8) :: TrTh_min, TrTh_max, TrPh_min, TrPh_max
        Real(kind=8) :: RecTh_min, RecTh_max, RecPh_min, RecPh_max
        
        Real(kind=8) :: ap,X,q_int,s,c,ang,cg,s_it,Sn,simpsonRule
        Real(kind=8) :: th_i,ph_i,Kix,Kiy,Kiz,Ki,th_s,ph_s,Ksx,Ksy,Ksz,Ks
        Real(kind=8) :: qe_int,qa_int,qs_int,qb_int,g_int,c_sca,q_sd
        Real(kind=8) :: Q_ext_av, Q_abs_av, Q_sca_av, Q_bks_av,g_av
        COMPLEX(real64) :: valVV, valHH, valVH, valHV
        COMPLEX(real64), Dimension(:,:), allocatable :: Svv_2d,Shh_2d,Svh_2d,Shv_2d
        COMPLEX(real64), Dimension(:), allocatable :: Svv_1d,Shh_1d,Svh_1d,Shv_1d
        Real(kind=8), Dimension(:), allocatable :: Thetas,Phis,RecThetasVals,RecPhisVals
        Real(kind=8), Dimension(:), allocatable :: TrThetasVals,TrPhisVals,s_t
        Real(kind=8), Dimension(:), allocatable :: xth, wth, xph, wph
        Real(kind=8), Dimension(:), allocatable :: x_leb,y_leb,z_leb,w,fxy_1d,fxy_g_1d
        Real(kind=8), Dimension(:,:), allocatable :: fxy_2d,fy,fxy_g_2d,fy_g,fy_ext,fy_abs,fy_sca,fy_bks
        Real(kind=8), Dimension(:,:), allocatable :: Q_ext,Q_sca,Q_abs,Q_bks,g,Q_ext_intf,Q_abs_intf
        
        if (rank .eq. 0) then 
            Write (*,*) ''
            Write (*,*) '-------------------- Scattered Quantities (',NumIntType_t,'/',NumIntType_r,') --------------------' 
            Write (*,*) ''
        EndIf
        ap = SimScatterer%dm/2.;
        X = K_air*SimScatterer%dm/2.
        Sfold_name = trim(SimOutfld_name)//Env_sep//'S_files';
        Qfold_name = trim(SimOutfld_name)//Env_sep//'Q_files';

        NTr_wr_proc = (NTr/nber_procs)+1;
        If (nom_methode=='CBFM-E  ') Then
            Allocate(character(6) ::nom_meth_exact)
            nom_meth_exact = trim(nom_methode)
        Elseif ((nom_methode=='MoM     ') .OR. (nom_methode=='RGE     ')) Then 
            Allocate(character(3) ::nom_meth_exact)
            nom_meth_exact = trim(nom_methode)
        Endif
        
        Allocate(Q_ext(NTrTheta,NTrPhi))
        Allocate(Q_sca(NTrTheta,NTrPhi))
        Allocate(Q_abs(NTrTheta,NTrPhi))
        Allocate(Q_bks(NTrTheta,NTrPhi))
        Allocate(g(NTrTheta,NTrPhi))
        Allocate(Q_ext_intf(NTrTheta,NTrPhi))
        Allocate(Q_abs_intf(NTrTheta,NTrPhi));
        
        ! Transmitters theta and phi values
        Allocate(Thetas(NTr),Phis(NTr))
        Allocate(TrThetasVals(NTrTheta),TrPhisVals(NTrPhi))
        kkt=1;
        Do jj=1, NTrPhi        
            Do ii=1, NTrTheta                        
                Thetas(kkt) = Transmitters(kkt)%theta;
                Phis(kkt) = Transmitters(kkt)%phi;
                kkt=kkt+1;
            EndDo
        EndDo
        call Unique1DArray_D(NTr,Nths,Thetas)
        call Unique1DArray_D(NTr,Nphs,Phis)
        
        TrThetasVals= Thetas(1:Nths); TrPhisVals= Phis(1:Nphs);
        Deallocate(Thetas,Phis)
        TrTh_min = TrThetasVals(1)*Pi/180.; TrTh_max = TrThetasVals(Nths)*Pi/180.;
        TrPh_min = TrPhisVals(1)*Pi/180.; TrPh_max =TrPhisVals(Nphs)*Pi/180.;     
        
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
        
        
        if ((NumIntType_r .eq. 'sd') .OR. (NumIntType_r .eq. 'lb')) then
            Allocate(Svv_1d(NRx),Shh_1d(NRx));
            Allocate(Svh_1d(NRx),Shv_1d(NRx));            
                 
            if (NumIntType_r == 'lb') Then !! Lebedev Qudrature
              ! prepare weight if Lebedev (otherwise w = 1 everywhere)
              Npts = NRx;
              Allocate(w(Npts))   
              Allocate(x_leb(Npts),y_leb(Npts),z_leb(Npts))
              call ld_by_order (Npts,x_leb,y_leb,z_leb,w) ;
              deallocate(x_leb,y_leb,z_leb);
            endif
        else
            Allocate(Svv_2d(NRxTheta,NRxPhi),Shh_2d(NRxTheta,NRxPhi))
            Allocate(Svh_2d(NRxTheta,NRxPhi),Shv_2d(NRxTheta,NRxPhi))
            
            ! Receivers theta and phi vals 
            Allocate(Thetas(NRx),Phis(NRx))
            Allocate(RecThetasVals(NRxTheta),RecPhisVals(NRxPhi))         
                
            Thetas(:) = Receivers(:)%theta;
            Phis(:) = Receivers(:)%phi;
            call Unique1DArray_D(NRx,Nths,Thetas)
            call Unique1DArray_D(NRx,Nphs,Phis)
        
            RecThetasVals= Thetas(1:Nths); RecPhisVals= Phis(1:Nphs);
            RecTh_min = RecThetasVals(1)*Pi/180.; RecTh_max = RecThetasVals(Nths)*Pi/180.;
            RecPh_min = RecPhisVals(1)*Pi/180.; RecPh_max =RecPhisVals(Nphs)*Pi/180.;
        
            ! check Step_in_cosTh
            if (nint(Thetas(2) - Thetas(1)) .eq. nint((Thetas(Nths)-Thetas(1))/(Nths-1))) then
                Step_in_cosTh = 0;
            else 
                Step_in_cosTh = 1;
            endif
            Deallocate(Thetas,Phis);
            
            ! needed to calculate the scattering extinction coefficients
            Sn = max((RecPh_max-RecPh_min),1.)*(cos(RecTh_min)-cos(RecTh_max));  ! if theta= 0:Pi and Phi=0:2Pi; Sn = 4Pi        
        endif
                               
        DO dd=1,NTrPhi
            Do cc=1,NTrTheta
                kkt=(dd-1)* NTrTheta + cc;
                
                th_i =  Transmitters(kkt)%theta;
                ph_i = Transmitters(kkt)%phi;
                Kix = - sin(th_i*Pi/180.)*cos(ph_i*Pi/180.)
                Kiy = - sin(th_i*Pi/180.)*sin(ph_i*Pi/180.)
                Kiz = - cos(th_i*Pi/180.) 
                !ki = sqrt(Kix**2.+Kiy**2.+Kiz**2.); !! logically equal to 1 but added just to check g 
                     
                ! Recovering S matrix elements taking into account
                ! ddscat and Mie code angle convention (forward direction <-->theta=0)
                if ((NumIntType_r .eq. 'sd') .OR. (NumIntType_r .eq. 'lb')) then 
                    Do kkr=1, NRx      !! LOOP ON RECEIVERS/SCATTERERS
                        Svv_1d(kkr)= S_total(kkr,4*(kkt-1)+1); 
                        Svh_1d(kkr) = S_total(kkr,4*(kkt-1)+2); 
                        Shv_1d(kkr) = S_total(kkr,4*(kkt-1)+3); 
                        Shh_1d(kkr)= S_total(kkr,4*(kkt-1)+4);      
                    EndDo 
                    
                    If (wr_Sij .eq. 1) Then
                        If ((kkt .gt. (rank*NTr_wr_proc)) .and. (kkt .le. (rank+1)*NTr_wr_proc)) then 
                          Write(kkt_st,'(a,i4.4)') 'kt',kkt;
                          if (EqSph == 0) then
                              file_name_s = trim(Sfold_name)//Env_sep//sim_name//'Smtable_'//stFreq//freq_unit//'_'//kkt_st//'_'//nom_meth_exact//'.dat';
                          else
                              file_name_s = trim(Sfold_name)//Env_sep//sim_name//'SmtableES_'//stFreq//freq_unit//'_'//kkt_st//'_'//nom_meth_exact//'.dat';
                          endif
                          Open(unit=21+rank,File = file_name_s)    
                          Write(21+rank,'(a,a)') '      theta       phi       Re(Svv)        Im(Svv)         Re(Svh)       Im(Shv) ',&
                                          '        Re(Shv)       Im(Shv)        Re(Shh)        Im(Shh) '
                          Do jj =1, NRxPhi  
                              Do ii =1, NRxTheta 
                                  Write(21+rank,'(f9.2,a,f9.2,a,e12.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4)') &
                                  Receivers(ii)%theta,';  ',Receivers(ii)%phi,';  ',Real(Svv_1d(ii)),';  ',Imag(Svv_1d(ii)),';  ',Real(Svh_1d(ii)),&
                                ';  ',Imag(Svh_1d(ii)),';  ', Real(Shv_1d(ii)),';  ',Imag(Shv_1d(ii)),';  ',Real(Shh_1d(ii)),';  ',Imag(Shh_1d(ii))
                              EndDo
                          EndDo
                          Close(21+rank);  
                        endif  
                    EndIf                  
                else                                    
                    Do jj=1, NRxPhi                   
                        Do ii=1, NRxTheta                        
                            kkr = (jj-1)*NRxTheta + ii;                                
                            Svv_2d(ii,jj) = S_total(kkr,4*(kkt-1)+1); 
                            Svh_2d(ii,jj) = S_total(kkr,4*(kkt-1)+2); 
                            Shv_2d(ii,jj) = S_total(kkr,4*(kkt-1)+3); 
                            Shh_2d(ii,jj) = S_total(kkr,4*(kkt-1)+4);    
                        EndDo        
                    EndDo 
                
                    If (wr_Sij .eq. 1) Then 
                        If ((kkt .gt. (rank*NTr_wr_proc)) .and. (kkt .le. (rank+1)*NTr_wr_proc)) then 
                        Write(kkt_st,'(a,i4.4)') 'kt',kkt;
                        if (EqSph == 0) then
                            file_name_s = trim(Sfold_name)//Env_sep//sim_name//'Smtable_'//stFreq//freq_unit//'_'//kkt_st//'_'//nom_meth_exact//'.dat';
                        else
                            file_name_s = trim(Sfold_name)//Env_sep//sim_name//'SmtableES_'//stFreq//freq_unit//'_'//kkt_st//'_'//nom_meth_exact//'.dat';
                        endif
                        Open(unit=21+rank,File = file_name_s)    
                        Write(21+rank,'(a,a)') '      theta       phi       Re(Svv)        Im(Svv)         Re(Svh)       Im(Shv) ',&
                                        '        Re(Shv)       Im(Shv)        Re(Shh)        Im(Shh) '
                        Do jj =1, NRxPhi  
                            Do ii =1, NRxTheta 
                                Write(21+rank,'(f9.2,a,f9.2,a,e12.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4)') &
                                RecThetasVals(ii),';  ',RecPhisVals(jj),';  ',Real(Svv_2d(ii,jj)),';  ',Imag(Svv_2d(ii,jj)),';  ',Real(Svh_2d(ii,jj)),&
                                    ';  ',Imag(Svh_2d(ii,jj)),';  ', Real(Shv_2d(ii,jj)),';  ',Imag(Shv_2d(ii,jj)),';  ',&
                                    Real(Shh_2d(ii,jj)),';  ',Imag(Shh_2d(ii,jj))
                            EndDo
                        EndDo
                        Close(21+rank);   
                        EndIf
                    EndIf     
                EndIf
                    
                ! Computing The scattering coefficients 
                ! Q_ext ************************************************************************************************************************
                ! FWD Direction 
                th_s = Transmitters(kkt)%Theta; ph_s = Transmitters(kkt)%phi;
                call getTxRxIndex(NRx_tot,Receivers,th_s,ph_s,ind_fwd);
                if (ind_fwd .eq. 0) then 
                    Write(*,'(a)') 'Error when calculating extenction cross section, fwd scattering direction not found !'
                    stop 0;
                endif                
                Q_ext(cc,dd) = (2/(X*X))*(abs(imag(S_total(ind_fwd,4*(kkt-1)+4)))+abs(imag(S_total(ind_fwd,4*(kkt-1)+1)))); ! SHH & SVV                  
                ! HERE Q_ext from internal Fiels
                Q_ext_intf(cc,dd) = C_ext(kkt)/(Pi*ap**2.);
                
                
                ! Q_scat & g ******************************************************************************************************************
                if ((NumIntType_r .eq. 'sd') .OR. (NumIntType_r .eq. 'lb')) then 
                    Allocate(fxy_1d(NRx));Allocate(fxy_g_1d(NRx));     
                    if (NumIntType_r == 'lb') then !LEBEDEV             
                        Do kkr =1, NRx
                            ! Q_scat
                            fxy_1d(kkr) = w(kkr)*(abs(Shh_1d(kkr))**2.+abs(Svv_1d(kkr))**2. &
                            +abs(Svh_1d(kkr))**2.+abs(Shv_1d(kkr))**2.);  
                                    
                            ! g
                            th_s = Receivers(kkr)%theta;  
                            ph_s = Receivers(kkr)%phi;
                            Ksx = -sin(th_s*Pi/180.)*cos(ph_s*Pi/180.)
                            Ksy = -sin(th_s*Pi/180.)*sin(ph_s*Pi/180.)
                            Ksz = -cos(th_s*Pi/180.)   
                            cg = kix*ksx+kiy*ksy+kiz*ksz 
                            fxy_g_1d(kkr) =fxy_1d(kkr)*cg;
                        EndDo  
                
                        ! Qscat 
                        q_sd = sum(fxy_1d); 
                        c_sca = 1./(2*K_air**2)*q_sd; 
                        Q_sca(cc,dd) = (4*Pi)/(2*Pi*X*X)*q_sd;
     
                        ! g                            
                        g(cc,dd) = 1./(2.*K_air**2.*c_sca)*sum(fxy_g_1d); 
                    elseif (NumIntType_r == 'sd') then ! SPHERICAL DESIGN
                        Do kkr =1, NRx
                            ! Q_scat
                            fxy_1d(kkr) = (abs(Shh_1d(kkr))**2.+abs(Svv_1d(kkr))**2. &
                            +abs(Svh_1d(kkr))**2.+abs(Shv_1d(kkr))**2.);  
                                    
                            ! g
                            th_s = Receivers(kkr)%theta;  
                            ph_s = Receivers(kkr)%phi;
                            Ksx = -sin(th_s*Pi/180.)*cos(ph_s*Pi/180.)
                            Ksy = -sin(th_s*Pi/180.)*sin(ph_s*Pi/180.)
                            Ksz = -cos(th_s*Pi/180.)   
                            cg = kix*ksx+kiy*ksy+kiz*ksz 
                            fxy_g_1d(kkr) =fxy_1d(kkr)*cg;
                        EndDo  
                
                        ! Qscat 
                        q_sd = (4*Pi/NRx)*sum(fxy_1d); 
                        c_sca = 1./(2*K_air**2)*q_sd; 
                        Q_sca(cc,dd) = 1/(2*Pi*X*X)*q_sd;
                        ! g           
                        g(cc,dd) = 1./(2.*K_air**2.*c_sca)*(4*Pi/NRx)*sum(fxy_g_1d); 
                    endif      
                    Deallocate(fxy_1d,fxy_g_1d);
                else
                    Q_sca(cc,dd)= 0;
                    g(cc,dd)= 0;                
                    Allocate(fxy_2d(NRxTheta,NRxPhi));Allocate(fxy_g_2d(NRxTheta,NRxPhi));
                    ! prepare the integrand fxy   
                    Do ii =1, NRxTheta                    
                        th_s =  RecThetasVals(ii);    
                        s = sin(th_s*Pi/180.);
                        Do jj =1, NRxPhi 
                            ! Q_scat
                            s_it = abs(Shh_2d(ii,jj))**2.+abs(Svv_2d(ii,jj))**2.+abs(Svh_2d(ii,jj))**2.+abs(Shv_2d(ii,jj))**2.;     
                            fxy_2d(ii,jj) = s_it*s;                          
                            ! g
                            ph_s = RecPhisVals(jj); ! ph_s used to calculate the scalar product ki.ks (for the asymetry parameter g)
                            Ksx = -sin(th_s*Pi/180.)*cos(ph_s*Pi/180.)
                            Ksy = -sin(th_s*Pi/180.)*sin(ph_s*Pi/180.)
                            Ksz = -cos(th_s*Pi/180.)                          
                            cg = (kix*ksx+kiy*ksy+kiz*ksz);!/(Ki*Ks); ! cg = cos(angle between Ki and Ks)
                            fxy_g_2d(ii,jj) =s_it*s*cg; 
                        EndDo
                    EndDo                 
                
                    ! now apply the Gausse Legendre or the Adaptive Quadrature
                    if (NumIntType_r == 'gl') Then ! Gausse Legendre
                        Allocate(xth(NRxTheta),wth(NRxTheta));
                        Allocate(xph(NRxPhi),wph(NRxPhi));
                        call getGaussLegendreQuadPts3(NRxTheta,RecTh_min,RecTh_max,xth,wth) ! LEGENDRE_RULE_FAST
                        if (NRxPhi > 1) Then
                            call getGaussLegendreQuadPts3(NRxPhi,RecPh_min,RecPh_max,xph,wph);
                        else
                            wph(1) = 1.;
                        endif                    
                    
                        wth = wth(NRxTheta:1:-1);
                        wph = wph(NRxPhi:1:-1);
                    
                        ! Qscat
                        q_int = 0;
                        Do ii =1, NRxTheta 
                            Do jj =1, NRxPhi     
                                q_int = q_int + wph(jj)*wth(ii)*fxy_2d(ii,jj);        
                            EndDo
                        EndDo 
                        c_sca = 1./(2.*K_air**2.)*q_int; 
                        Q_sca(cc,dd) = 2./(Sn*X*X)*q_int;  
                        
                        !g
                        q_int = 0;
                        Do ii =1, NRxTheta 
                            Do jj =1, NRxPhi     
                                q_int = q_int + wph(jj)*wth(ii)*fxy_g_2d(ii,jj);        
                            EndDo
                        EndDo 
                        g(cc,dd) = 1./(2*K_air**2.*c_sca)*q_int;   
                        Deallocate(xth,wth,xph,wph);  
                    
                    elseif (NumIntType_r == 'aq') Then ! Adaptive Quadrature
                        q_int = 0;
                    
                        !Q_scat
                        Allocate(fy(NRxPhi,1));
                        Do jj =1, NRxPhi
                            call getAdaptiveQuadPts2(NRxTheta,NRxPhi,NRxTheta,&
                                NRxPhi,fxy_2d,'th',jj,RecTh_min,RecTh_max,fy(jj,1))   
                        EndDo
                        if (NRxPhi > 1) Then
                            call getAdaptiveQuadPts2(NRxTheta,NRxPhi,NRxPhi,1,fy,&
                                    'ph',1,RecPh_min,RecPh_max,q_int);               
                        else
                            q_int = fy(1,1);                                                   
                        endif
                    
                        Q_sca(cc,dd) = 2./(Sn*X*X)*q_int;
                        c_sca = 1./(2*K_air**2)*q_int; 
                        deallocate(fy);    
                        
                        !g
                        q_int=0
                        allocate(fy_g(NRxPhi,1));
                        Do jj =1, NRxPhi
                            call getAdaptiveQuadPts2(NRxTheta,NRxPhi,NRxTheta,&
                                NRxPhi,fxy_g_2d,'th',jj,RecTh_min,RecTh_max,fy_g(jj,1))   
                        EndDo
                        if (NRxPhi > 1) Then 
                            call getAdaptiveQuadPts2(NRxTheta,NRxPhi,NRxPhi,1,fy_g,&
                                    'ph',1,RecPh_min,RecPh_max,q_int); 
                        else
                            q_int = fy_g(1,1);                        
                        endif
                        g(cc,dd) = 1./(2.*K_air**2.*c_sca)*q_int;   
                        deallocate(fy_g);                   
                    
                    elseif (NumIntType_r == 'tr') Then
                        ! Q_scat
                        call trap_2Dc(fxy_2d,RecTh_min,RecTh_max,RecPh_min,RecPh_max,q_int,NRxTheta,NRxPhi);
                        c_sca = 1./(2*K_air**2)*q_int; 
                        Q_sca(cc,dd) = 2./(Sn*X*X)*q_int; 
                        ! g 
                        call trap_2Dc(fxy_g_2d,RecTh_min,RecTh_max,RecPh_min,RecPh_max,q_int,NRxTheta,NRxPhi);
                        g(cc,dd) = 1./(2.*K_air**2.*c_sca)*q_int;                
                    endIf
                    Deallocate(fxy_2d,fxy_g_2d);
                EndIf
                                
                ! Q_abs = Qext - Qsca **************************************************************************************************
                Q_abs(cc,dd) = Q_ext(cc,dd) - Q_sca(cc,dd);
                ! and from internal field 
                Q_abs_intf(cc,dd) = C_abs(kkt)/(Pi*ap**2.);
                
                
                ! Q_backscattering *****************************************************************************************************
                ! BKW Direction 
                th_s = (180.-Transmitters(kkt)%Theta); ph_s = mod(Transmitters(kkt)%Phi+180.,360.);
                call getTxRxIndex(NRx_tot,Receivers,th_s,ph_s,ind_bkw);
                if (ind_bkw .eq. 0) then 
                    Write(*,'(a)') 'Error when calculating backscattering cross section, bkw scattering direction not found !'
                    stop 0;
                endif 
                Q_bks(cc,dd) = (1./(2.*Pi*X*X))*(abs(S_total(ind_bkw,4*(kkt-1)+1))**2.+abs(S_total(ind_bkw,4*(kkt-1)+4))**2.& ! SVV & SHH & SVH & SHV
                +abs(S_total(ind_bkw,4*(kkt-1)+2))**2.+abs(S_total(ind_bkw,4*(kkt-1)+3))**2.);   
                ! sphere : Qback = (1/(Pi*X*X))*abs(S1(2*Nang-1))^2.; ! TO CHECK and VALIDATE
                !Q_bks(cc,dd) = (1./(Pi*X*X))*abs(S_total(ind_bkw,4*(kkt-1)+4))**2.
            EndDo            
         EndDo
         if (NumIntType_r .eq. 'lb') then
             Deallocate(Svv_1d,Shh_1d,Svh_1d,Shv_1d,w); 
         elseif (NumIntType_r .eq. 'sd') then   
             Deallocate(Svv_1d,Shh_1d,Svh_1d,Shv_1d);               
         else
             Deallocate(RecThetasVals,RecPhisVals,Svv_2d,Shh_2d,Svh_2d,Shv_2d);
         endif
         
        
        ! IF DISPLAY Qext and Qabs from INTERNAL FIELD (Yurkin & Hoekstra)
        if (QextIFDisp == 1) then 
            Q_ext=0D0; Q_ext = Q_ext_intf;
            Qextintfst = ' (IntField)';
        else
            Qextintfst = ' (FarField)';            
        endif        
        Q_abs = Q_abs_intf;         
        
        ! WRT Q per id ----------------------------------------------------------------------------------------------------
        If ((rank .eq. 0) .and. (wr_Qij .eq. 1)) Then
            if (EqSph == 0) then
                file_name_q = trim(Qfold_name)//Env_sep//sim_name//'Qidtable_'//stFreq//freq_unit//'_'//nom_meth_exact//'.dat';
            else
                file_name_q = trim(Qfold_name)//Env_sep//sim_name//'QidtableES_'//stFreq//freq_unit//'_'//nom_meth_exact//'.dat';
            endif
            kkt=0;
            Open(unit=41,File = trim(file_name_q));    
            Write(41,'(a,a)') '     theta    phi     ',&
                'Q_ext        Q_abs      Q_scat      Q_bk      g(1)=<cos> '
            Do jj =1, NTrPhi  
                Do ii =1, NTrTheta 
                    kkt = kkt + 1 ;
                    Write(41,'(f9.2,f9.2,es12.4,es12.4,es12.4,es12.4,es12.4)') Transmitters(kkt)%theta,Transmitters(kkt)%phi, &
                        Q_ext(ii,jj), Q_abs(ii,jj), Q_sca(ii,jj),Q_bks(ii,jj),g(ii,jj)
                EndDo
            EndDo
            Close(41);    
        EndIf 
        !!-------------------------------------------------------------------------------------------------------------------
        
               
        !********************************************************************************************************
        ! HERE WE COMPUTE THE AVARAGED SCATTERING QUANTITIES USING Q_ext,Q_sca,Q_abs,Q_bks per direction
        ! of propagation (transmitter) 
        !********************************************************************************************************
        if (NTr > 1) Then
            
            Sn = max((TrPh_max-TrPh_min),1.)*(cos(TrTh_min)-cos(TrTh_max));
            
            Do ii=1, NTrTheta    
                s = sin(Transmitters(ii)%theta*Pi/180.);    
                Do jj=1, NTrPhi 
                    Q_ext(ii,jj) = Q_ext(ii,jj) * s;
                    Q_sca(ii,jj) = Q_sca(ii,jj) * s;
                    Q_bks(ii,jj) = Q_bks(ii,jj) * s;
                    Q_abs(ii,jj) = Q_abs(ii,jj) * s;
                    g(ii,jj) = g(ii,jj) * s;                                                                                     
                EndDo
            EndDo
                    
            Q_ext_av = 0; Q_abs_av = 0; 
            Q_sca_av = 0; Q_bks_av = 0; 
            g_av=0;
            
            !! GAUSS-LEGENDRE
            If (NumIntType_t == 'gl') Then            
                Allocate(xth(NTrTheta),wth(NTrTheta));
                Allocate(xph(NTrPhi),wph(NTrPhi));
                call getGaussLegendreQuadPts3(NTrTheta,TrTh_min,TrTh_max,xth,wth) ! LEGENDRE_RULE_FAST
                if (NTrPhi > 1) Then
                    call getGaussLegendreQuadPts3(NTrPhi,TrPh_min,TrPh_max,xph,wph);
                else
                    wph(1) = 1.
                endif                
            
                wth = wth(NTrTheta:1:-1);
                wph = wph(NTrPhi:1:-1);
                
                qe_int = 0;qa_int =0; qs_int = 0;qb_int = 0; g_int=0;
                Do ii =1, NTrTheta 
                    Do jj =1, NTrPhi     
                        qe_int = qe_int + wph(jj)*wth(ii)*Q_ext(ii,jj);     
                        qa_int = qa_int + wph(jj)*wth(ii)*Q_abs(ii,jj);     
                        qs_int = qs_int + wph(jj)*wth(ii)*Q_sca(ii,jj);
                        qb_int = qb_int + wph(jj)*wth(ii)*Q_bks(ii,jj);   
                        g_int = g_int + wph(jj)*wth(ii)*g(ii,jj);   
                    EndDo
                EndDo 
                
            elseif (NumIntType_t == 'aq') Then 
                ! Adaptive quadrature to compute the averaged scattering quantities                
                allocate(fy_ext(NTrPhi,1),fy_sca(NTrPhi,1));
                allocate(fy_bks(NTrPhi,1),fy_g(NTrPhi,1));
                allocate(fy_abs(NTrPhi,1))
                Do jj =1, NTrPhi
                    call getAdaptiveQuadPts2(NTrTheta,NTrPhi,NTrTheta,NTrPhi,&
                        Q_ext,'th',jj,TrTh_min,TrTh_max,fy_ext(jj,1)) 
                    call getAdaptiveQuadPts2(NTrTheta,NTrPhi,NTrTheta,NTrPhi,&
                        Q_abs,'th',jj,TrTh_min,TrTh_max,fy_abs(jj,1)) 
                    call getAdaptiveQuadPts2(NTrTheta,NTrPhi,NTrTheta,NTrPhi,&
                        Q_sca,'th',jj,TrTh_min,TrTh_max,fy_sca(jj,1)) 
                    call getAdaptiveQuadPts2(NTrTheta,NTrPhi,NTrTheta,NTrPhi,&
                        Q_bks,'th',jj,TrTh_min,TrTh_max,fy_bks(jj,1))
                    call getAdaptiveQuadPts2(NTrTheta,NTrPhi,NTrTheta,NTrPhi,&
                        g,'th',jj,TrTh_min,TrTh_max,fy_g(jj,1))
                EndDo
                if (NTrPhi > 1) Then
                    call getAdaptiveQuadPts2(NTrTheta,NTrPhi,NTrPhi,1,fy_ext,&
                            'ph',1,TrPh_min,TrPh_max,qe_int); 
                    call getAdaptiveQuadPts2(NTrTheta,NTrPhi,NTrPhi,1,fy_abs,&
                        'ph',1,TrPh_min,TrPh_max,qa_int);
                    call getAdaptiveQuadPts2(NTrTheta,NTrPhi,NTrPhi,1,fy_sca,&
                        'ph',1,TrPh_min,TrPh_max,qs_int);
                    call getAdaptiveQuadPts2(NTrTheta,NTrPhi,NTrPhi,1,fy_bks,&
                        'ph',1,TrPh_min,TrPh_max,qb_int);
                    call getAdaptiveQuadPts2(NTrTheta,NTrPhi,NTrPhi,1,fy_g,&
                        'ph',1,TrPh_min,TrPh_max,g_int);                    
                else
                    qe_int=fy_ext(1,1); qa_int=fy_abs(1,1);
                    qs_int=fy_sca(1,1); qb_int=fy_bks(1,1); g_int=fy_g(1,1);
                endif                
                deallocate(fy_ext,fy_sca,fy_bks,fy_g);
                
            elseif (NumIntType_t == 'tr') Then ! cubature trapezoid rule
                call trap_2Dc(Q_ext,TrTh_min,TrTh_max,TrPh_min,TrPh_min,qe_int,NTrTheta,NTrPhi);
                call trap_2Dc(Q_abs,TrTh_min,TrTh_max,TrPh_min,TrPh_min,qa_int,NTrTheta,NTrPhi);
                call trap_2Dc(Q_sca,TrTh_min,TrTh_max,TrPh_min,TrPh_min,qs_int,NTrTheta,NTrPhi);
                call trap_2Dc(Q_bks,TrTh_min,TrTh_max,TrPh_min,TrPh_min,qb_int,NTrTheta,NTrPhi);
                call trap_2Dc(g,TrTh_min,TrTh_max,TrPh_min,TrPh_min,g_int,NTrTheta,NTrPhi);                
            EndIf 
            
            ! Average over Sn
            Q_ext_av = 1./Sn*qe_int;
            Q_abs_av = 1./Sn*qa_int;
            Q_sca_av = 1./Sn*qs_int;
            Q_bks_av = 1./Sn*qb_int;
            g_av = 1./Sn*g_int; 
        Else
            Q_ext_av = Q_ext(1,1); Q_abs_av = Q_abs(1,1); 
            Q_sca_av = Q_sca(1,1); Q_bks_av = Q_bks(1,1); g_av = g(1,1);
        EndIf  
        
        Deallocate(Q_ext,Q_sca,Q_abs,Q_bks)    
        
        if (rank .eq. 0) then 
            if (EqSph == 0) then
                file_name_q = trim(SimOutfld_name)//Env_sep//'qtable_'//nom_meth_exact//'.dat'
            else
                file_name_q = trim(SimOutfld_name)//Env_sep//'qtableES_'//nom_meth_exact//'.dat'            
            endif
        
            If (num_freq > 1) Then
                Open(unit=41,File = trim(file_name_q), Access='Append', Status='old')    
            Else
                Open(unit=41,File = trim(file_name_q)) 
            EndIf
    
            If (num_freq == 1) Then
                Write(41, '(a,i8)') 'Number of Cells = ',Nbc 
                Write(41, '(a,i5,a,a,a)') 'Results averaged over ',NTr,' incident directions (',NumIntType_t,') :'
                Write(41, '(a,f6.2,a,f6.2,a,i4)') 'Theta_i_min = ',TrTh_min/pi*180.,'; Theta_i_max = ',TrTh_max/pi*180. ,'; NTheta_i = ',NTrTheta
                Write(41, '(a,f6.2,a,f6.2,a,i4)') 'Phi_i_min   = ',TrPh_min/pi*180.,'; Phi_i_max   = ', TrPh_max/pi*180.,'; NPhi_i   = ',NTrPhi
                Write(41, '(a,i5,a,a,a)') 'Results calculated with ',NRx,' scattering directions (',NumIntType_r,') :'
                if ((NumIntType_r == 'aq') .OR. (NumIntType_r == 'gl') .OR. (NumIntType_r == 'tr')) then 
                    Write(41, '(a,f6.2,a,f6.2,a,i4)') 'Theta_s_min = ',RecTh_min/pi*180.,'; Theta_s_max = ',RecTh_max/pi*180.,'; NTheta_s = ',NRxTheta
                    Write(41, '(a,f6.2,a,f6.2,a,i4)') 'Phi_s_min   = ',RecPh_min/pi*180.,'; Phi_s_max   = ',RecPh_max/pi*180.,'; NPhi_s   = ',NRxPhi 
                endif
                Write(41,'(a,a,a)') '      freq       aeff        wave',&
                    '       Q_ext        Q_abs      Q_scat      Q_bk      g(1)=<cos>   Ncels'        
            EndIf 
            if ((freq_unit == 'THz') .OR. (freq_unit == 'GHz')) then        
                Write(41,'(es12.4,es12.4,es12.4,es12.4,es12.4,es12.4,es12.4,es12.4,i8)') &
                    Freq_w,ap*1E6, lambda_w*1E6, Q_ext_av, Q_abs_av, Q_sca_av, &
                    Q_bks_av,g_av, Nbc
            else
                Write(41,'(es12.4,es12.4,es12.4,es12.4,es12.4,es12.4,es12.4,es12.4,i8)') &
                    Freq_w,ap, lambda_w, Q_ext_av, Q_abs_av, Q_sca_av, &
                    Q_bks_av,g_av, Nbc
            endif
    
            Close(41)
            Write(*,'(a,es12.4,a)') '- Q_ext = ', Q_ext_av, trim(Qextintfst)
            Write(*,'(a,es12.4)') '- Q_scat = ', Q_sca_av
            Write(*,'(a,es12.4)') '- Q_abs = ', Q_abs_av
            Write(*,'(a,es12.4)') '- Q_bks  = ', Q_bks_av
            Write(*,'(a,es12.4)') '- g(1)   = ', g_av
            Write(*,'(a,es12.4)') ' ' 
        EndIf        

    End Subroutine Compute_Scattering_Quantities_1
    
    
    
    SUBROUTINE Compute_Scattering_Quantities_2(nom_methode,SimScatterer,Transmitters,Receivers,S_total,C_ext,C_abs)
    
        USE Initialization
        USE common_variables
        USE MPI

        !! "Implicit Statement" 
        IMPLICIT NONE

        !IN/OUT
        character(8), INTENT(IN):: nom_methode
        type (Scatterer), INTENT(IN) :: SimScatterer
        type (Dipole), Dimension(NTr), INTENT(IN) :: Transmitters
        type (Dipole), Dimension(NRx_tot), INTENT(IN) :: Receivers
        COMPLEX(real64), Dimension(NRx_tot,4*NTr), INTENT(IN):: S_total
        COMPLEX(real64), Dimension(NTr), INTENT(IN):: C_ext,C_abs
        
        ! Local 
        Integer :: ii,jj,cc,dd,kkr,kkt,a,N,ind_fw_dir,ind_bkw_dir,Nths,Nphs
        Integer :: Nber_scattering_dirs,start_at,end_at,Npts
        Integer :: id,nthreads,p,d,NTr_wr_proc 
        CHARACTER(:), allocatable:: nom_meth_exact,stFreq,sim_name
        CHARACTER(200) :: file_name_s,file_name_q
        CHARACTER(200) :: Sfold_name,Qfold_name
        CHARACTER(6) :: ty
        CHARACTER(6) :: kkt_st
        CHARACTER(11) :: Qextintfst,Qabsintfst
        
        Real(kind=8) :: ap,X,ang,cg,s_it,Sn
        Real(kind=8) :: th_i,ph_i,Kix,Kiy,Kiz,th_s,ph_s,Ksx,Ksy,Ksz
        Real(kind=8) :: RecTh_min, RecTh_max, RecPh_min, RecPh_max
        Real(kind=8) :: q_sd,q_int,c_sca,Q_ext_av,Q_abs_av,Q_sca_av,Q_bks_av,g_av,s
        COMPLEX(real64) :: valVV, valHH, valVH, valHV
        COMPLEX(real64), Dimension(:), allocatable :: Svv_1d,Shh_1d,Svh_1d,Shv_1d
        COMPLEX(real64), Dimension(:,:), allocatable :: Svv_2d,Shh_2d,Svh_2d,Shv_2d
        Real(kind=8), Dimension(:), allocatable :: Thetas,Phis,RecThetasVals,RecPhisVals
        Real(kind=8), Dimension(:), allocatable :: x_leb,y_leb,z_leb,w
        Real(kind=8), Dimension(:), allocatable :: fxy_1d,fxy_g_1d
        Real(kind=8), Dimension(:,:), allocatable :: fxy_2d,fxy_g_2d,fy,fy_g
        Real(kind=8), Dimension(:), allocatable :: Q_ext,Q_sca,Q_abs,Q_bks,g,Q_ext_intf,Q_abs_intf  
        
        Write (*,*) ''
        Write (*,*) '-------------------- Scattred Quantities (',NumIntType_t,'/',NumIntType_r,') --------------------' 
        Write (*,*) ''
        
        ap = SimScatterer%dm/2.;
        X = K_air*SimScatterer%dm/2.
        
        Sfold_name = trim(SimOutfld_name)//Env_sep//'S_files';
        Qfold_name = trim(SimOutfld_name)//Env_sep//'Q_files';
        
        NTr_wr_proc = (NTr/nber_procs)+1;
        
        If (nom_methode=='CBFM-E  ') Then
            Allocate(character(6) ::nom_meth_exact)
            nom_meth_exact = trim(nom_methode)
        Elseif (nom_methode=='MoM     ') Then 
            Allocate(character(3) ::nom_meth_exact)
            nom_meth_exact = trim(nom_methode)
        Endif
        
        Allocate(Q_ext(NTr))
        Allocate(Q_ext_intf(NTr))
        Allocate(Q_sca(NTr))
        Allocate(Q_abs(NTr))
        Allocate(Q_abs_intf(NTr))
        Allocate(Q_bks(NTr))
        Allocate(g(NTr))
        
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
        
        ! prepare weight for Receivers if Lebedev (otherwise w = 1 everywhere)
        if (NumIntType_r == 'lb') Then !! Lebedev Qudrature
          Npts = NRx;
          Allocate(w(Npts)); 
          Allocate(x_leb(Npts),y_leb(Npts),z_leb(Npts))
          call ld_by_order (Npts,x_leb,y_leb,z_leb,w) ;
          deallocate(x_leb,y_leb,z_leb);
        endif
                
        if ((NumIntType_r == 'lb') .OR. (NumIntType_r == 'sd')) then 
            Allocate(Svv_1d(NRx),Shh_1d(NRx))
            Allocate(Svh_1d(NRx),Shv_1d(NRx))  
        else
            Allocate(Thetas(NRx),Phis(NRx))
            Allocate(RecThetasVals(NRxTheta),RecPhisVals(NRxPhi))         
                
            Thetas(1:NRx) = Receivers(1:NRx)%theta;
            Phis(1:NRx) = Receivers(1:NRx)%phi;
            call Unique1DArray_D(NRx,Nths,Thetas)
            call Unique1DArray_D(NRx,Nphs,Phis)
        
            RecThetasVals= Thetas(1:Nths); RecPhisVals= Phis(1:Nphs);
            RecTh_min = RecThetasVals(1)*Pi/180.; RecTh_max = RecThetasVals(Nths)*Pi/180.;
            RecPh_min = RecPhisVals(1)*Pi/180.; RecPh_max =RecPhisVals(Nphs)*Pi/180.;

	    ! needed to calculate the scattering extinction coefficients
            Sn = max((RecPh_max-RecPh_min),1.)*(cos(RecTh_min)-cos(RecTh_max));  ! if theta= 0:Pi and Phi=0:2Pi; Sn = 4Pi    
            
            Allocate(Svv_2d(NRxTheta,NRxPhi),Shh_2d(NRxTheta,NRxPhi))
            Allocate(Svh_2d(NRxTheta,NRxPhi),Shv_2d(NRxTheta,NRxPhi))
        endif
    
                       
        Do kkt=1,NTr  !! LOOP ON TRANSMITTERS
        
            th_i =  Transmitters(kkt)%theta;
            ph_i = Transmitters(kkt)%phi;
            Kix = - sin(th_i*Pi/180.)*cos(ph_i*Pi/180.)
            Kiy = - sin(th_i*Pi/180.)*sin(ph_i*Pi/180.)
            Kiz = - cos(th_i*Pi/180.)                         
                 
            ! Recovering S matrix elements taking into account
            ! ddscat and Mie code angle convention (forward direction <-->theta=0)
            if ((NumIntType_r .eq. 'sd') .OR. (NumIntType_r .eq. 'lb')) then 
                Do kkr=1, NRx      !! LOOP ON RECEIVERS/SCATTERERS
                    Svv_1d(kkr)= S_total(kkr,4*(kkt-1)+1); 
                    Svh_1d(kkr) = S_total(kkr,4*(kkt-1)+2); 
                    Shv_1d(kkr) = S_total(kkr,4*(kkt-1)+3); 
                    Shh_1d(kkr)= S_total(kkr,4*(kkt-1)+4);      
                EndDo 
            
                If (wr_Sij .eq. 1) Then
                    If ((kkt .gt. (rank*NTr_wr_proc)) .and. (kkt .le. (rank+1)*NTr_wr_proc)) then 
                      Write(kkt_st,'(a,i3.3)') 'kt',kkt;
                      if (EqSph == 0) then
                          file_name_s = trim(Sfold_name)//Env_sep//sim_name//'Smtable_'//stFreq//trim(freq_unit)//'_'//kkt_st//'_'//nom_meth_exact//'.dat';
                      else
                          file_name_s = trim(Sfold_name)//Env_sep//sim_name//'SmtableES_'//stFreq//trim(freq_unit)//'_'//kkt_st//'_'//nom_meth_exact//'.dat';
                      endif
                  
                      Open(unit=21+rank,File = file_name_s)    
                      Write(21+rank,'(a,a)') '      theta       phi       Re(Svv)        Im(Svv)         Re(Svh)       Im(Shv) ',&
                                      '        Re(Shv)       Im(Shv)        Re(Shh)        Im(Shh) '
                      Do ii =1, NRx 
                          Write(21+rank,'(f9.2,a,f9.2,a,e12.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4)') &
                          Receivers(ii)%theta,';  ',Receivers(ii)%phi,';  ',Real(Svv_1d(ii)),';  ',Imag(Svv_1d(ii)),';  ',Real(Svh_1d(ii)),&
                            ';  ',Imag(Svh_1d(ii)),';  ', Real(Shv_1d(ii)),';  ',Imag(Shv_1d(ii)),';  ',Real(Shh_1d(ii)),';  ',Imag(Shh_1d(ii))
                      EndDo
                      Close(21+rank);    
                    endif
                EndIf       
            Else
                Do jj=1, NRxPhi                   
                    Do ii=1, NRxTheta                        
                        kkr = (jj-1)*NRxTheta + ii;                                
                        Svv_2d(ii,jj) = S_total(kkr,4*(kkt-1)+1); 
                        Svh_2d(ii,jj) = S_total(kkr,4*(kkt-1)+2); 
                        Shv_2d(ii,jj) = S_total(kkr,4*(kkt-1)+3); 
                        Shh_2d(ii,jj) = S_total(kkr,4*(kkt-1)+4);    
                    EndDo        
                EndDo 
                
                If (wr_Sij .eq. 1) Then 
                    Write(kkt_st,'(a,i4.4)') 'kt',kkt;
                    if (EqSph == 0) then
                        file_name_s = trim(Sfold_name)//Env_sep//sim_name//'Smtable_'//stFreq//trim(freq_unit)//'_'//kkt_st//'_'//nom_meth_exact//'.dat';
                    else
                        file_name_s = trim(Sfold_name)//Env_sep//sim_name//'SmtableES_'//stFreq//trim(freq_unit)//'_'//kkt_st//'_'//nom_meth_exact//'.dat';
                    endif
                    Open(unit=20+id,File = file_name_s)    
                    Write(20+id,'(a,a)') '      theta       phi       Re(Svv)        Im(Svv)         Re(Svh)       Im(Shv) ',&
                                    '        Re(Shv)       Im(Shv)        Re(Shh)        Im(Shh) '
                    Do jj =1, NRxPhi  
                        Do ii =1, NRxTheta 
                            Write(20+id,'(f9.2,a,f9.2,a,e12.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4)') &
                            RecThetasVals(ii),';  ',RecPhisVals(jj),';  ',Real(Svv_2d(ii,jj)),';  ',Imag(Svv_2d(ii,jj)),';  ',Real(Svh_2d(ii,jj)),&
                                ';  ',Imag(Svh_2d(ii,jj)),';  ', Real(Shv_2d(ii,jj)),';  ',Imag(Shv_2d(ii,jj)),';  ',&
                                Real(Shh_2d(ii,jj)),';  ',Imag(Shh_2d(ii,jj))
                        EndDo
                    EndDo
                    Close(20+id);    
                EndIf    
            EndIf
                         
            ! Computing The scattering coefficients 
            ! Q_ext *******
            if ((NTr .eq. NRx) .and. ((NumIntType_r == 'lb') .OR. (NumIntType_r == 'sd'))) then 
                ind_fw_dir = kkt ;! the forward direction is the receiver = the current transmitter 
            else
                ind_fw_dir = NRx + kkt 
            endif            
            Q_ext(kkt) = (2/(X*X))*(abs(imag(S_total(ind_fw_dir,4*(kkt-1)+4)))+abs(imag(S_total(ind_fw_dir,4*(kkt-1)+1))));
            ! calculated from internal field
            Q_ext_intf(kkt) = C_ext(kkt)/(Pi*ap**2.);
            
            ! Q_scat & g ******
            if (NumIntType_r == 'lb') then !LEBEDEV  
                Allocate(fxy_1d(NRx));Allocate(fxy_g_1d(NRx));     
                Do kkr =1, NRx
                    ! Q_scat
                    fxy_1d(kkr) = w(kkr)*(abs(Shh_1d(kkr))**2.+abs(Svv_1d(kkr))**2. &
                    +abs(Svh_1d(kkr))**2.+abs(Shv_1d(kkr))**2.);  
                                    
                    ! g
                    th_s = Receivers(kkr)%theta;  
                    ph_s = Receivers(kkr)%phi;
                    Ksx = -sin(th_s*Pi/180.)*cos(ph_s*Pi/180.)
                    Ksy = -sin(th_s*Pi/180.)*sin(ph_s*Pi/180.)
                    Ksz = -cos(th_s*Pi/180.)   
                    cg = kix*ksx+kiy*ksy+kiz*ksz 
                    fxy_g_1d(kkr) =fxy_1d(kkr)*cg;
                EndDo                                 
                ! Qscat 
                q_sd = sum(fxy_1d); 
                c_sca = 1./(2*K_air**2)*q_sd; 
                Q_sca(kkt) = (4*Pi)/(2*Pi*X*X)*q_sd;
     
                ! g                            
                g(kkt) = 1./(2.*K_air**2.*c_sca)*sum(fxy_g_1d);
		Deallocate(fxy_1d,fxy_g_1d); 
            elseif (NumIntType_r == 'sd') then ! SPHERICAL DESIGN
                Allocate(fxy_1d(NRx));Allocate(fxy_g_1d(NRx));     
                Do kkr =1, NRx
                    ! Q_scat
                    fxy_1d(kkr) = (abs(Shh_1d(kkr))**2.+abs(Svv_1d(kkr))**2. &
                    +abs(Svh_1d(kkr))**2.+abs(Shv_1d(kkr))**2.);  
                                    
                    ! g
                    th_s = Receivers(kkr)%theta;  
                    ph_s = Receivers(kkr)%phi;
                    Ksx = -sin(th_s*Pi/180.)*cos(ph_s*Pi/180.)
                    Ksy = -sin(th_s*Pi/180.)*sin(ph_s*Pi/180.)
                    Ksz = -cos(th_s*Pi/180.)   
                    cg = kix*ksx+kiy*ksy+kiz*ksz 
                    fxy_g_1d(kkr) =fxy_1d(kkr)*cg;
                EndDo  
                
                ! Qscat 
                q_sd = (4*Pi/NRx)*sum(fxy_1d); 
                c_sca = 1./(2*K_air**2)*q_sd; 
                Q_sca(kkt) = 1/(2*Pi*X*X)*q_sd;
                ! g           
                g(kkt) = 1./(2.*K_air**2.*c_sca)*(4*Pi/NRx)*sum(fxy_g_1d); 
                
                Deallocate(fxy_1d,fxy_g_1d);
            else
                Allocate(fxy_2d(NRxTheta,NRxPhi));Allocate(fxy_g_2d(NRxTheta,NRxPhi));
                ! prepare the integrand fxy   
                Do ii =1, NRxTheta                    
                    th_s =  RecThetasVals(ii);    
                    s = sin(th_s*Pi/180.);
                    Do jj =1, NRxPhi 
                        ! Q_scat
                        s_it = abs(Shh_2d(ii,jj))**2.+abs(Svv_2d(ii,jj))**2.+abs(Svh_2d(ii,jj))**2.+abs(Shv_2d(ii,jj))**2.;     
                        fxy_2d(ii,jj) = s_it*s;                          
                        ! g
                        ph_s = RecPhisVals(jj); ! ph_s used to calculate the scalar product ki.ks (for the asymetry parameter g)
                        Ksx = -sin(th_s*Pi/180.)*cos(ph_s*Pi/180.)
                        Ksy = -sin(th_s*Pi/180.)*sin(ph_s*Pi/180.)
                        Ksz = -cos(th_s*Pi/180.)                          
                        cg = (kix*ksx+kiy*ksy+kiz*ksz);!/(Ki*Ks); ! cg = cos(angle between Ki and Ks)
                        fxy_g_2d(ii,jj) =s_it*s*cg; 
                    EndDo
                EndDo                 
                
                q_int = 0;                    
                !Q_scat
                allocate(fy(NRxPhi,1));
                Do jj =1, NRxPhi
                    call getAdaptiveQuadPts2(NRxTheta,NRxPhi,NRxTheta,&
                        NRxPhi,fxy_2d,'th',jj,RecTh_min,RecTh_max,fy(jj,1))   
                EndDo
                if (NRxPhi > 1) Then
                    call getAdaptiveQuadPts2(NRxTheta,NRxPhi,NRxPhi,1,fy,&
                            'ph',1,RecPh_min,RecPh_max,q_int);               
                else
                    q_int = fy(1,1);                                                   
                endif
                    
                Q_sca(kkt) = 2./(Sn*X*X)*q_int;
                c_sca = 1./(2*K_air**2)*q_int; 
                                        
                !g
                q_int=0
                allocate(fy_g(NRxPhi,1));
                Do jj =1, NRxPhi
                    call getAdaptiveQuadPts2(NRxTheta,NRxPhi,NRxTheta,&
                        NRxPhi,fxy_g_2d,'th',jj,RecTh_min,RecTh_max,fy_g(jj,1))   
                EndDo
                if (NRxPhi > 1) Then 
                    call getAdaptiveQuadPts2(NRxTheta,NRxPhi,NRxPhi,1,fy_g,&
                            'ph',1,RecPh_min,RecPh_max,q_int); 
                else
                    q_int = fy_g(1,1);                        
                endif
                g(kkt) = 1./(2.*K_air**2.*c_sca)*q_int;   
                deallocate(fy,fy_g); 
                Deallocate(fxy_2d,fxy_g_2d);
            endif  
            
            ! Q_abs
            Q_abs(kkt) = Q_ext(kkt) - Q_sca(kkt);
            ! from internal field 
            Q_abs_intf(kkt) = C_abs(kkt)/(Pi*ap**2.);
            
            ! bkwd direction
            if ((NTr .eq. NRx) .and. ((NumIntType_r == 'lb') .OR. (NumIntType_r == 'sd'))) then 
                ind_bkw_dir = NRx + kkt ; ! see get_trans_Receiv to understand how the receivers are cerated and organized
            else
                ind_bkw_dir = NRx + NTr + kkt; 
            endif               
            Q_bks(kkt) = (1./(2.*Pi*X*X))*(abs(S_total(ind_bkw_dir,4*(kkt-1)+1))**2.+abs(S_total(ind_bkw_dir,4*(kkt-1)+4))**2.&
            +abs(S_total(ind_bkw_dir,4*(kkt-1)+2))**2.+abs(S_total(ind_bkw_dir,4*(kkt-1)+3))**2.);
            ! sphere : Qback = (1/(Pi*X*X))*abs(S1(2*Nang-1))^2.; ! TO CHECK and VALIDATE
            !Q_bks(cc,dd) = (1./(Pi*X*X))*abs(Shh(ind_bkw_rec_th,ind_bkw_rec_ph))**2.
        EndDo  !! LOOP ON TRANSMITTERS
        

	if (NumIntType_r == 'lb') then 
            Deallocate(Svv_1d,Shh_1d,Svh_1d,Shv_1d,w);
        elseif (NumIntType_r == 'sd') then 
            Deallocate(Svv_1d,Shh_1d,Svh_1d,Shv_1d);
        else
            Deallocate(Svv_2d,Shh_2d,Svh_2d,Shv_2d);
	EndIf
        
        
        ! IF DISPLAY Qext and Qabs from INTERNAL FIELD (Yurkin & Hoekstra)
        if (QextIFDisp == 1) then 
            Q_ext=0D0; Q_ext = Q_ext_intf;
            Qextintfst = ' (IntField)';
        else
            Qextintfst = ' (FarField)';            
        endif
        Q_abs=0D0; Q_abs = Q_abs_intf;
                    
        !********************************************************************************************************
        ! HERE WE COMPUTE THE AVARAGED SCATTERING QUANTITIES USING Q_ext,Q_sca,Q_abs,Q_bks per direction
        ! of propagation (transmitter) 
        !********************************************************************************************************
        If (NTr > 1) Then    
            
            N = NTr; 
            
            if (NumIntType_t == 'sd') then        
                Q_ext_av = sum(Q_ext)/N; ! sd : (4*Pi/N)*sum(Q_ext) + divided by solid angle/(4*Pi);
                Q_abs_av = sum(Q_abs)/N; 
                Q_sca_av = sum(Q_sca)/N; 
                Q_bks_av = sum(Q_bks)/N;  
                g_av= sum(g)/N;      
            else
                Allocate(w(N));  
                Allocate(x_leb(N),y_leb(N),z_leb(N))
                call ld_by_order (N,x_leb,y_leb,z_leb,w) ;
                deallocate(x_leb,y_leb,z_leb);
                Q_ext_av=0;  Q_abs_av=0; Q_sca_av=0; Q_bks_av=0;g_av=0;
                do ii=1,N
                    Q_ext_av = Q_ext_av + w(ii)*Q_ext(ii);
                    Q_abs_av = Q_abs_av + w(ii)*Q_abs(ii);
                    Q_sca_av = Q_sca_av + w(ii)*Q_sca(ii);
                    Q_bks_av = Q_bks_av + w(ii)*Q_bks(ii);
                    g_av= g_av+w(ii)*g(ii);              
                enddo            
            endif
        Else
            Q_ext_av = Q_ext(1); Q_abs_av = Q_abs(1); Q_sca_av = Q_sca(1); Q_bks_av = Q_bks(1); g_av = g(1);
        EndIf 
                
        ! WRT Q per id 
        If ((rank .eq. 0) .and. (wr_Qij .eq. 1)) Then
            if (EqSph == 0) then
                file_name_q = trim(Qfold_name)//Env_sep//sim_name//'Qidtable_'//stFreq//trim(freq_unit)//'_'//nom_meth_exact//'.dat'
            else
                file_name_q = trim(Qfold_name)//Env_sep//sim_name//'QidtableES_'//stFreq//trim(freq_unit)//'_'//nom_meth_exact//'.dat'
            endif
            Open(unit=41,File = trim(file_name_q));    
             Write(41,'(a,a)') '     theta    phi     ',&
                'Q_ext        Q_abs      Q_scat      Q_bk      g(1)=<cos> '
            Do ii =1, NTr 
                Write(41,'(f9.2,f9.2,es12.4,es12.4,es12.4,es12.4,es12.4)') Transmitters(ii)%theta,Transmitters(ii)%phi, &
                    Q_ext(ii), Q_abs(ii), Q_sca(ii),Q_bks(ii),g(ii)
            EndDo
            Close(41);    
        EndIf            
        Deallocate(Q_ext,Q_sca,Q_abs,Q_bks,Q_ext_intf,Q_abs_intf);    
        
        if (rank .eq. 0) then 
            if (EqSph == 0) then 
                file_name_q = trim(SimOutfld_name)//Env_sep//'qtable_'//nom_meth_exact//'.dat'
            else
                file_name_q = trim(SimOutfld_name)//Env_sep//'qtableES_'//nom_meth_exact//'.dat'            
            endif
        
            If (num_freq > 1) Then
                Open(unit=41,File = trim(file_name_q), Access='Append', Status='old')    
            Else
                Open(unit=41,File = trim(file_name_q)) 
            EndIf
    
            If (num_freq == 1) Then
                Write(41, '(a,i8)') 'Number of Cells = ',Nbc 
                Write(41, '(a,i5,a,a,a)') 'Results averaged over ',NTr,' incident directions (',NumIntType_t,') :'
                Write(41, '(a,i5,a,a,a)') 'Results calculated with ',NRx,' scattering directions (',NumIntType_r,') :'
                if ((NumIntType_r == 'aq') .OR. (NumIntType_r == 'gl') .OR. (NumIntType_r == 'tr')) then 
                    Write(41, '(a,f6.2,a,f6.2,a,i4)') 'Theta_s_min = ',RecTh_min/pi*180.,'; Theta_s_max = ',RecTh_max/pi*180.,'; NTheta_s = ',NRxTheta
                    Write(41, '(a,f6.2,a,f6.2,a,i4)') 'Phi_s_min   = ',RecPh_min/pi*180.,'; Phi_s_max   = ',RecPh_max/pi*180.,'; NPhi_s   = ',NRxPhi 
                endif
                Write(41,'(a,a,a)') '      freq       aeff        wave',&
                    '       Q_ext        Q_abs      Q_scat      Q_bk      g(1)=<cos>   Ncels'        
            EndIf    
            
            if ( (freq_unit == 'THz') .OR. (freq_unit == 'GHz')) then   
                Write(41,'(es12.4,es12.4,es12.4,es12.4,es12.4,es12.4,es12.4,es12.4,i8)') &
                    Freq_w,ap*1E6, lambda_w*1E6, Q_ext_av, Q_abs_av, Q_sca_av, &
                    Q_bks_av,g_av, Nbc
            else
                Write(41,'(es12.4,es12.4,es12.4,es12.4,es12.4,es12.4,es12.4,es12.4,i8)') &
                Freq_w,ap, lambda_w, Q_ext_av, Q_abs_av, Q_sca_av, &
                 Q_bks_av,g_av, Nbc
            endif
    
            Close(41)
            Write(*,'(a,es12.4,a)') '- Q_ext = ', Q_ext_av, trim(Qextintfst)
            Write(*,'(a,es12.4)') '- Q_scat = ', Q_sca_av
            Write(*,'(a,es12.4)') '- Q_abs = ', Q_abs_av
            Write(*,'(a,es12.4)') '- Q_bks  = ', Q_bks_av
            Write(*,'(a,es12.4)') '- g(1)   = ', g_av
            Write(*,'(a,es12.4)') ' ';
        EndIf

    End Subroutine Compute_Scattering_Quantities_2
    
    subroutine Unique1DArray_D(N_in,N_out,Arr_a)
        ! Similar function as matlab unique to remove redundent elements.
        ! Author: Kong, kinaxj@gmail.com
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N_in
        INTEGER, INTENT(OUT) :: N_out
        real*8,DIMENSION(N_in),INTENT(INOUT):: Arr_a
        
        real*8,DIMENSION(:), allocatable:: Arr_b
        LOGICAL,DIMENSION(:), allocatable::mask
        INTEGER,DIMENSION(:),allocatable::index_vector,indexSos
        INTEGER::i,j,num
        num=size(Arr_a);  ALLOCATE(mask(num)); mask = .TRUE.
        DO i=num,2,-1
            mask(i)=.NOT.(ANY(Arr_a(:i-1)==Arr_a(i)))
        END DO
        ! Make an index vector
        allocate(indexSos(size(PACK([(i,i=1,num)],mask))))
        ALLOCATE(index_vector(size(indexSos))); index_vector=PACK([(i,i=1,num)],mask)
    
        ! Now copy the unique elements of a into b
        ALLOCATE(Arr_b(size(index_vector)))
        Arr_b=Arr_a(index_vector)
        Arr_a = 0D0;
        
        N_out= size(index_vector);
        Arr_a(1:N_out) = Arr_b;
        deallocate(Arr_b);        
    end subroutine Unique1DArray_D

    