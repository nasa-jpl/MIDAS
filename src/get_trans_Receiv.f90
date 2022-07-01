SUBROUTINE get_trans_Receiv(Ninc_in,Nscat_in,Transmitters_Comp,Receivers);

    ! Calculates Transmitters and Receivers directions (theta_i, Phi_i) and (Theta_s, Phi_s) respectively
    ! Inputs :
    ! Ninc_in  (Integer) : number of incident directions
    ! Nscat_in (Integer) : number of scattering directions
    ! thitrans,thftrans, phitrans,phftrans (Real kind=8) : limit values of theta_i, phi_i, theta_s and phi_s (theta_i/s = 0:180 and phi_i/s=0:360 when averaging over the entire sphere surface)

    ! Outputs
    ! Transmitters_Comp (Dipole type) : transmitters (theta_i,phi_i)
    ! Receivers         (Dipole type) : Receivers (theta_s,phi_s)

    ! modif 2-18-2020 : correction of the transformation (X,Y,Z) --> (theha, phi), It impacts the SD and LB configurations (line 118)
    ! We also replaced Phi+Pi by mod(Phi+2Pi,2Pi). thinking that +Pi doesn't change the outcoming result was wrong, mod(Phi+2Pi,2Pi) is the correct way to convert -Pi<Phi<Pi to 0<Phi<2Pi
    ! modif 2-24-2020 : correction of the call to ld_by_order. a subroutine lb_get_closer_Npts was created in sphere_lebedev_rule.f90 et the changes impacted get_trans_Receiv.f90,
    ! Compute_Scattering_Quantities.f90 and getParamaters_CBFM.f90

    USE Initialization
    USE common_variables

    IMPLICIT NONE

    !! IN/OUT ******************************************************************

    Integer, INTENT(IN) :: Ninc_in,Nscat_in
    type (Dipole), Dimension(:), allocatable, INTENT(OUT) :: Transmitters_Comp
    type (Dipole), Dimension(:), allocatable, INTENT(OUT) :: Receivers

    ! local
    Integer :: ii,Ind, I,K,order,Npts,NRx_extra,Step_in_cosTh
    Integer :: RxExt_exists, RxBks_exists

    !! Transmitters/Receivers
    Real(kind=8) :: thftrans,thitrans,phftrans,phitrans
    Real(kind=8) :: thfRecei,thiRecei,phfRecei,phiRecei
    Real(kind=8) :: Step_theta_trans_comp,step_phi_trans_comp
    Real(kind=8) :: step_theta_Recei,step_phi_Recei
    Real(kind=8) :: margin_trans_theta,margin_trans_phi,cosdth_init
    Real(kind=8) :: margin_rec_phi,margin_rec_theta
    Real(kind=8) :: step_theta_CBFM,step_phi_CBFM,theta, phi
    type (Dipole), Dimension(:), allocatable :: Receivers_tmp

    ! local Gauss Legendre
    Real(kind=8), Dimension(:), allocatable :: xth, wth, xph, wph

    ! spherical design & lebedev
    Real(kind=8), dimension(:),allocatable :: hypotxy
    Real(kind=8), dimension(:,:), allocatable :: XYZstd
    Real(kind=8), dimension(:),allocatable :: thetas_rd, phis_rd

    ! local Lebedev
    integer ( kind = 4 ), parameter :: nmax = 65
    integer ( kind = 4 ), parameter :: mmax = &
      ( ( nmax * 2 + 3 ) * ( nmax * 2 + 3 ) / 3 )

    real ( kind = 8 ) w(mmax)
    real ( kind = 8 ) x(mmax)
    real ( kind = 8 ) y(mmax)
    real ( kind = 8 ) z(mmax)

    ! initialize from Common_variables
    thitrans = theta_init_trans_comp; thftrans = theta_final_trans_comp
    phitrans = phi_init_trans_comp; phftrans = phi_final_trans_comp
    thiRecei = theta_init_Recei; thfRecei = theta_final_Recei
    phiRecei = phi_init_Recei; phfRecei = phi_final_Recei

    ! Step_in_cosTh = 1 if the angular step is considered in cosTh instead of Th
    ! (to avoid the problem of concentration around the poles)
    Step_in_cosTh = 0;

    sd_type = 2 ; ! we implemented different versions of spehrical t-design
                 ! sd_type = 1 Hardin and Sloane Spherical Designs : http://people.sc.fsu.edu/~jburkardt%20/f_src/sphere_design_rule/sphere_design_rule.html
                 ! sd_type = 2 Efficient Spherical T-Designs : https://web.maths.unsw.edu.au/~rsw/Sphere/EffSphDes/index.html
                 ! sd_type = 3 Symmetric Efficient Spherical T-Designs : https://web.maths.unsw.edu.au/~rsw/Sphere/EffSphDes/index.html

    if ((NumIntType_t .eq. 'sm') .and. &
        ((mod(NTrTheta,2) .ne. 1) .OR. (mod(NTrPhi,2) .ne. 1) .OR. &
        (mod(NRxTheta,2) .ne. 1) .OR. (mod(NRxPhi,2) .ne.1))) then
        NumIntType_t = 'aq'; ! switch automatically to adaptive quadrature
        EndIf


    !! For the moment (5/8/2019) we focus on the most appropriate and usefull NumIntType_t/NumIntType_r combined configuration
    !! The uniform config (NumIntType_t=NumIntType_r) are all available.
    if ((NumIntType_t .eq. 'gl') .OR. (NumIntType_t .eq. 'sm') .OR. (Step_in_cosTh .eq. 1)) then
        NumIntType_r = NumIntType_t;
    EndIf

    !!****************************************************************************
    !! TR ************************************************************************
    !!****************************************************************************
    ! read from IncScattDirs.dat file
    If (NumIntType_t .eq. 'rf') Then

        !! Open the dat file and read the simulation parameters.
        Open(41,File = 'inputs'//Env_sep//'InScattDirs.dat')
        Do I=1,5
            Read(41,*)
        EndDo

        NTr = NTrTheta*NTrPhi
        Allocate(Transmitters_Comp(NTr))
        Ind = 1
        DO I=1, NTrPhi
            DO K=1, NTrTheta
                Read(41,'(f9.4,f9.4)') theta, phi
                Transmitters_Comp(Ind) = Dipole(theta,phi)
                Ind = Ind + 1;
            Enddo
        Enddo

        Close(41);
        
    ! Spherical T_designs
    Elseif (NumIntType_t .eq. 'sd') Then
        ! Spherical Design
        if (sd_type .eq. 1) then
            call design_closer_order (Ninc_in,order,Npts);
        elseif (sd_type .eq. 2) then
            call ESD_design_closer_order (Ninc_in,order,Npts);
        elseif (sd_type .eq. 3) then
            call ESSD_design_closer_order (Ninc_in,order,Npts);
        else
            write(*,'(a)') 'invalid value of sd_type; set automatically to 2'
            sd_type = 2;
            call ESD_design_closer_order (Ninc_in,order,Npts);
        endif
        NTr = Npts;

        Allocate(XYZstd(3,Npts));
        if (sd_type .eq. 1) then
            call design_points(order,Npts,XYZstd);
        elseif (sd_type .eq. 2) then
            call ESD_design_points(order,Npts,XYZstd);
        elseif (sd_type .eq. 3) then
            call ESSD_design_points(order,Npts,XYZstd);
        endif

        !!! Test 2/11/2020
        !Open(21,File = 'ESSTD_XYZ.dat');
        !Do ii=1,Npts
        !    Write(21,'(es25.16,a,es25.16,a,es25.16,a)')  XYZstd(1,ii),',', XYZstd(2,ii),',',XYZstd(3,ii),', &'
        !EndDo
        !Close(21);
        !!!***********************************************************************************



        Allocate(hypotxy(Npts));
        ! convert cartesian to spherical
        hypotxy = hypot(XYZstd(2,:),XYZstd(3,:)); !hypotxy = sqrt(abs(x).^2 + abs(y).^2)

        Allocate(thetas_rd(Npts),phis_rd(Npts));
        thetas_rd = Pi/2.-atan2(XYZstd(1,:),hypotxy);
        phis_rd = atan2(XYZstd(3,:),XYZstd(2,:)); ! az = atan2(y,x);

        !!! Test 2/11/2020
        !Open(21,File = 'ESSTD_angles.dat');
        !Do ii=1,Npts
        !    Write(21,'(f10.3,f10.3)')  180.*thetas_rd(ii)/Pi,180.*phis_rd(ii)/Pi;
        !EndDo
        !Close(21);
        !!!***********************************************************************************

        Allocate(Transmitters_Comp(Npts));
        DO Ind=1, Npts
            theta = thetas_rd(Ind)/Pi*180.;
            phi =  (mod(phis_rd(Ind)+2.*Pi,2.*Pi))/Pi*180.
            Transmitters_Comp(Ind) = Dipole(theta,phi)
        Enddo
        deallocate(thetas_rd,phis_rd,XYZstd,hypotxy);

    !! Lebedev Quadrature
    ElseIf (NumIntType_t .eq. 'lb') Then

        ! Lebedev quadrature
        ! get the closer available order to the suggested Ninc_in
        call lb_get_closer_Npts(Ninc_in,Npts);
        call ld_by_order (Npts,x,y,z,w) ; ! for LB order = Npts
        NTr = Npts;

        Allocate(XYZstd(3,Npts));
        XYZstd(1,1:Npts) = x(1:Npts) ;
        XYZstd(2,1:Npts) = y(1:Npts) ;
        XYZstd(3,1:Npts) = z(1:Npts) ;

        Allocate(hypotxy(Npts));
        ! convert cartesian to spherical
        hypotxy = hypot(XYZstd(2,:),XYZstd(3,:)); !hypotxy = sqrt(abs(x).^2 + abs(y).^2)

        Allocate(thetas_rd(Npts),phis_rd(Npts));
        thetas_rd = Pi/2.-atan2(XYZstd(1,:),hypotxy);
        phis_rd = atan2(XYZstd(3,:),XYZstd(2,:)); ! az = atan2(y,x);

        Allocate(Transmitters_Comp(Npts));
        DO Ind=1, Npts
            theta = thetas_rd(Ind)/Pi*180.;
            phi =  (mod(phis_rd(Ind)+2.*Pi,2.*Pi))/Pi*180. ; ! mod(Phi+2Pi,2Pi) to have phi between 0 and 2Pi instead of -Pi->Pi
            Transmitters_Comp(Ind) = Dipole(theta,phi)
        Enddo
        deallocate(thetas_rd,phis_rd,XYZstd,hypotxy);

    ! gl : Gauss Legendre
    Elseif (NumIntType_t .eq. 'gl') then

        ! Transmitters / Incident Directions *********************************
        Allocate(xth(NTrTheta),wth(NTrTheta));
        Allocate(xph(NTrPhi),wph(NTrPhi));
        call getGaussLegendreQuadPts3(NTrTheta,thitrans*Pi/180.,thftrans*Pi/180.,xth,wth) ! LEGENDRE_RULE_FAST
        call getGaussLegendreQuadPts3(NTrPhi,phitrans*Pi/180.,phftrans*Pi/180.,xph,wph) ! LEGENDRE_RULE_FAST

        xth = xth(NTrTheta:1:-1)/Pi*180.
        xph = xph(NTrPhi:1:-1)/Pi*180.

        NTr = NTrTheta*NTrPhi
        Allocate(Transmitters_Comp(NTr))
        Ind = 1;
        DO I=1, NTrPhi
            phi = xph(I);
            DO K=1, NTrTheta
                theta = xth(K);
                Transmitters_Comp(Ind) = Dipole(theta,phi)
                Ind = Ind + 1;
            Enddo
        Enddo
        Deallocate(xth,wth,xph,wph);

    ! aq : Adaptive Quadrature OR tr : cubature trapezoid rule OR sm : Simpson Rule
    Elseif (Step_in_cosTh .eq. 0) then ! here aq .or. tr .or. sm with step in Th

        ! Transmitters / Incident Directions *********************************
        NTr = NTrTheta*NTrPhi
        margin_trans_theta =  (thftrans - thitrans);
        margin_trans_phi =  (phftrans - phitrans);
        ! step_theta_trans_comp
        If (NTrTheta == 1) Then
            step_theta_trans_comp = 0;
        Else
            step_theta_trans_comp = margin_trans_theta/(NTrTheta-1);
        EndIf

        !step_phi_trans_comp
        If (NTrPhi == 1) Then
            step_phi_trans_comp = 0;
        Else
            step_phi_trans_comp = margin_trans_phi/(NTrPhi-1);
        EndIf

        Allocate(Transmitters_Comp(NTr))
        Ind = 1;
        theta = 0
        phi = 0
        DO I=1, NTrPhi
            phi = phitrans + (I-1) * step_phi_trans_comp
            DO K=1, NTrTheta
                theta = thitrans + (K-1) * step_theta_trans_comp
                Transmitters_Comp(Ind) = Dipole(theta,phi)
                Ind = Ind + 1;
            Enddo
        Enddo

    Else ! here gl .or. aq .or. tr .or. sm with step in CosTh
        ! Transmitters
        NTr = NTrTheta*NTrPhi;
        cosdth_init = cos(thitrans*Pi/180.);
        margin_trans_theta =  cos(thftrans*Pi/180.)-cos(thitrans*Pi/180.);
        margin_trans_phi =  (phftrans - phitrans);
        ! step_theta_trans_comp
        If (NTrTheta == 1) Then
            step_theta_trans_comp = 0;
        Else
            step_theta_trans_comp = margin_trans_theta/(NTrTheta-1);
        EndIf

        !step_phi_trans_comp
        If (NTrPhi == 1) Then
            step_phi_trans_comp = 0;
        Else
            step_phi_trans_comp = margin_trans_phi/(NTrPhi-1);
        EndIf

        Allocate(Transmitters_Comp(NTr))
        Ind = 1;
        theta = 0
        phi = 0
        DO I=1, NTrPhi
            phi = phitrans + (I-1) * step_phi_trans_comp
            DO K=1, NTrTheta
                theta = acos(1.*(K-1)*step_theta_trans_comp+cosdth_init)*180./Pi;

                Transmitters_Comp(Ind) = Dipole(theta,phi)
                Ind = Ind + 1;
            Enddo
        Enddo
    EndIf


    !!****************************************************************************
    !! Rx ************************************************************************
    !!****************************************************************************
    If (NumIntType_r .eq. 'rf') Then

        !! Open the dat file and read the simulation parameters.
        Open(41,File = 'inputs'//Env_sep//'InScattDirs.dat')
        NTr = NTrTheta*NTrPhi
        Do I=1,NTr+11 ! 11 = 5+1+5
            Read(41,*)
        EndDo

        NRx = NTr ;
        !NRx_tot = NRx
        NRx_tot = NRx+NTr; ! for this option, given the values I have right now from Mark I will just add the bkw directions
        Allocate(Receivers(NRx_tot));
        Ind = 1
        DO I=1, NRxPhi
            DO K=1, NRxTheta
                Read(41,'(f9.4,f9.4)') theta, phi
                Receivers(Ind) = Dipole(theta,phi)
                Ind = Ind + 1;
            Enddo
        Enddo

        ! Now add the NTr bkw directions for the NTr transmitters
        DO I=1, NTr
            ! bkw direction
            theta = (180.-Transmitters_Comp(I)%Theta);
            phi = mod(Transmitters_Comp(I)%Phi+180.,360.);
            Receivers(Ind) = Dipole(theta,phi);
            Ind = Ind + 1;
        EndDo

        Close(41);
    Elseif (NumIntType_r .eq. 'sd') Then

        ! Spherical Design
        if (sd_type .eq. 1) then
            call design_closer_order (Nscat_in,order,Npts);
        elseif (sd_type .eq. 2) then
            call ESD_design_closer_order (Nscat_in,order,Npts);
        elseif (sd_type .eq. 3) then
            call ESSD_design_closer_order (Nscat_in,order,Npts);
        endif

        NRx = Npts;
        if ((NTr .eq. NRx) .and. (NumIntType_t .eq. 'sd')) then
            NRx_tot = NRx+NTr; ! = Npts + NTr = 2*NTr as NTr = NRx = Npts
            Allocate(Receivers(NRx_tot));
            Receivers(1:Npts)= Transmitters_Comp(1:NTr)
            ! Now bkw directions
            Do Ind=1, NTr
                theta = (180.-Transmitters_Comp(Ind)%Theta);
                phi = mod(Transmitters_Comp(Ind)%Phi+180.,360.);
                Receivers(NTr+Ind) = Dipole(theta,phi)
            EndDo
        else
            ! this applies to (NumIntType_t = aq) or (NumIntType_t = sd) with (NTr .ne. NRx) including (NTr .eq. 1)
            NRx_tot = NRx+2*NTr; !the incident direction (0,0) and its bkw direction (180,180)

            Allocate(Receivers(NRx_tot));
            Allocate(XYZstd(3,Npts));
            if (sd_type .eq. 1) then
                call design_points(order,Npts,XYZstd);
            elseif (sd_type .eq. 2) then
                call ESD_design_points(order,Npts,XYZstd);
            elseif (sd_type .eq. 3) then
                call ESSD_design_points(order,Npts,XYZstd);
            endif

            Allocate(hypotxy(Npts));
            ! convert cartesian to spherical
            hypotxy = hypot(XYZstd(2,:),XYZstd(3,:)); !hypotxy = sqrt(abs(x).^2 + abs(y).^2)

            Allocate(thetas_rd(Npts),phis_rd(Npts));
            thetas_rd = Pi/2.-atan2(XYZstd(1,:),hypotxy);
            phis_rd = atan2(XYZstd(3,:),XYZstd(2,:)); ! az = atan2(y,x);

            DO Ind=1, Npts
                theta = thetas_rd(Ind)/Pi*180.;
                phi =  (mod(phis_rd(Ind)+2.*Pi,2.*Pi))/Pi*180. ; ! mod(Phi+2Pi,2Pi) to have phi between 0 and 2Pi instead of -Pi->Pi
                Receivers(Ind) = Dipole(theta,phi)
            Enddo
            deallocate(thetas_rd,phis_rd,XYZstd,hypotxy);

            ! Now add the NTr fwd and NTr bkw directions for the NTr transmitters
            DO Ind=1, NTr
                theta = Transmitters_Comp(Ind)%Theta
                phi = Transmitters_Comp(Ind)%Phi;
                Receivers(Npts+Ind) = Dipole(theta,phi)
                ! bkw direction
                theta = (180.-Transmitters_Comp(Ind)%Theta);
                phi = mod(Transmitters_Comp(Ind)%Phi+180.,360.);
                Receivers(Npts+NTr+Ind) = Dipole(theta,phi)
            EndDo
        endif


    ElseIf (NumIntType_r .eq. 'lb') Then
        call lb_get_closer_Npts(Nscat_in,Npts); !for LB order = Npts
        NRx = Npts;

        if ((NTr .eq. NRx) .and. (NumIntType_t .eq. 'lb')) then
            NRx_tot = Npts+NTr;
            Allocate(Receivers(NRx_tot));
            Receivers(1:Npts)= Transmitters_Comp(1:NTr)
            ! Now bkw directions
            Do Ind=1, NTr
                theta = (180.-Transmitters_Comp(Ind)%Theta);
                phi = mod(Transmitters_Comp(Ind)%Phi+180.,360.);
                Receivers(NTr+Ind) = Dipole(theta,phi)
            EndDo
        else
            ! this applies to (NumIntType_t = aq) or (NumIntType_t = sd) with (NTr .ne. NRx) including (NTr .eq. 1)
            NRx_tot = Npts+2*NTr; !the incident direction (0,0) and its bkw direction (180,180)
            Allocate(Receivers(NRx_tot));

            call ld_by_order (Npts,x,y,z,w) ;
            Allocate(XYZstd(3,Npts));
            XYZstd(1,1:Npts) = x(1:Npts) ;
            XYZstd(2,1:Npts) = y(1:Npts) ;
            XYZstd(3,1:Npts) = z(1:Npts) ;

            Allocate(hypotxy(Npts));
            ! convert cartesian to spherical
            hypotxy = hypot(XYZstd(2,:),XYZstd(3,:)); !hypotxy = sqrt(abs(x).^2 + abs(y).^2)

            Allocate(thetas_rd(Npts),phis_rd(Npts));
            thetas_rd = Pi/2.-atan2(XYZstd(1,:),hypotxy);
            phis_rd = atan2(XYZstd(3,:),XYZstd(2,:)); ! az = atan2(y,x);

            DO Ind=1, Npts
                theta = thetas_rd(Ind)/Pi*180.;
                phi =  (mod(phis_rd(Ind)+2.*Pi,2.*Pi))/Pi*180. ; ! mod(Phi+2Pi,2Pi) to have phi between 0 and 2Pi instead of -Pi->Pi
                Receivers(Ind) = Dipole(theta,phi)
            Enddo
            deallocate(thetas_rd,phis_rd,XYZstd,hypotxy)

            ! Now add the NTr fwd and NTr bkw directions for the NTr transmitters
            DO Ind=1, NTr
                theta = Transmitters_Comp(Ind)%Theta
                phi = Transmitters_Comp(Ind)%Phi;
                Receivers(Npts+Ind) = Dipole(theta,phi)
                ! bkw direction
                theta = (180.-Transmitters_Comp(Ind)%Theta);
                phi = mod(Transmitters_Comp(Ind)%Phi+180.,360.);
                Receivers(Npts+NTr+Ind) = Dipole(theta,phi)
            EndDo
        endif


    ElseIf (NumIntType_r .eq. 'gl') Then
        ! RECEIVERS / Scattering Directions
        Allocate(xth(NRxTheta),wth(NRxTheta));
        Allocate(xph(NRxPhi),wph(NRxPhi));
        call getGaussLegendreQuadPts3(NRxTheta,thiRecei*Pi/180.,thfRecei*Pi/180.,xth,wth) ! LEGENDRE_RULE_FAST
        call getGaussLegendreQuadPts3(NRxPhi,phiRecei*Pi/180.,phfRecei*Pi/180.,xph,wph) ! LEGENDRE_RULE_FAST

        xth = xth(NRxTheta:1:-1)/Pi*180.
        xph = xph(NRxPhi:1:-1)/Pi*180.

        NRx = NRxTheta*NRxPhi
        NRx_tot = NRx + 2*NTr ! temporarily
        Allocate(Receivers(NRx_tot))
        Ind = 1;
        DO I=1, NRxPhi
            phi = xph(I);
            DO K=1, NRxTheta
                theta = xth(K);
                Receivers(Ind) = Dipole(theta,phi)
                Ind = Ind + 1;
            Enddo
        Enddo

        !! add extra-scatterers needed to calculate Cext and Cbks for all
        ! the available Transmitters/Incident directions
        NRx_extra = 0;
        Do I =1, NTr
            ! Forward direction
            theta = Transmitters_Comp(I)%Theta;
            phi = Transmitters_Comp(I)%Phi;
            call getTxRxIndex(NRx,Receivers(1:NRx),theta,phi,RxExt_exists);
            if (RxExt_exists .eq. 0) then ! add corresponding Rx
                NRx_extra = NRx_extra + 1;
                Receivers(NRx+NRx_extra) = Dipole(theta,phi);
            endif
            ! backward direction
            theta = (180.-Transmitters_Comp(I)%Theta);
            phi = mod(Transmitters_Comp(I)%Phi+180.,360.);
            call getTxRxIndex(NRx,Receivers(1:NRx),theta,phi,RxBks_exists);
            if (RxBks_exists .eq. 0) then ! add corresponding Rx
                NRx_extra = NRx_extra + 1;
                Receivers(NRx+NRx_extra) = Dipole(theta,phi);
            endif
        EndDo

        NRx_tot = NRx + NRx_extra;
        Allocate(Receivers_tmp(NRx_tot));
        Receivers_tmp(1:NRx_tot) = Receivers(1:NRx_tot);
        deallocate(Receivers); allocate(Receivers(NRx_tot));
        Receivers(1:NRx_tot) = Receivers_tmp(1:NRx_tot);
        deallocate(Receivers_tmp);

    ! aq : Adaptive Quadrature OR tr : cubature trapezoid rule OR sm : Simpson Rule
    Elseif (Step_in_cosTh .eq. 0) then ! here aq .or. tr .or. sm with step in Th
        NRx = NRxTheta*NRxPhi
        margin_rec_theta =  (thfRecei - thiRecei);
        margin_rec_phi =  (phfRecei - phiRecei);
        ! step_theta_Recei
        if (NRxTheta ==1) then
            step_theta_Recei = 0;
        else
            step_theta_Recei = margin_rec_theta/(NRxTheta-1);
        endif
        ! step_phi_Recei
        if (NRxPhi == 1) then
            step_phi_Recei = 0.;
        else
            step_phi_Recei = margin_rec_phi/(NRxPhi-1);
        endif

        if ((NumIntType_t .eq. 'sd') .OR. (NumIntType_t .eq. 'lb')) Then
            NRx_tot = NRx + 2*NTr
            Allocate(Receivers(NRx_tot));
            Ind = 1;
            theta = 0
            phi = 0
            DO I=1, NRxPhi
                phi = phiRecei + (I-1) * step_phi_Recei
                DO K=1, NRxTheta
                    theta = thiRecei + (K-1) * step_theta_Recei
                    Receivers(Ind) = Dipole(theta,phi)
                    Ind = Ind + 1;
                Enddo
            Enddo

            !! Now the receivers correponding to the fwd and bkw directions of the transmitters
            DO Ind=1, NTr
                theta = Transmitters_Comp(Ind)%Theta
                phi = Transmitters_Comp(Ind)%Phi;
                Receivers(NRx+Ind) = Dipole(theta,phi)
                ! bkw direction
                theta = (180.-Transmitters_Comp(Ind)%Theta);
                phi = mod(Transmitters_Comp(Ind)%Phi+180.,360.);
                Receivers(NRx+NTr+Ind) = Dipole(theta,phi)
            EndDo

        ElseIf (NumIntType_t .eq. 'gl') Then
            Write(*,'(a)') 'This combination of Tr/Rx is not developped yet';
            Stop 1;
        Else
            NRx_tot = NRx + 2*NTr ! temporarily
            Allocate(Receivers(NRx_tot))
            Ind = 1;
            theta = 0
            phi = 0
            DO I=1, NRxPhi
                phi = phiRecei + (I-1) * step_phi_Recei
                DO K=1, NRxTheta
                    theta = thiRecei + (K-1) * step_theta_Recei
                    Receivers(Ind) = Dipole(theta,phi)
                    Ind = Ind + 1;
                Enddo
            Enddo

            !! add extra-scatterers needed to calculate Cext and Cbks for all
            ! the available Transmitters/Incident directions
            NRx_extra = 0;
            Do I =1, NTr
                ! Forward direction
                theta = Transmitters_Comp(I)%Theta;
                phi = Transmitters_Comp(I)%Phi;
                call getTxRxIndex(NRx,Receivers(1:NRx),theta,phi,RxExt_exists);
                if (RxExt_exists .eq. 0) then ! add corresponding Rx
                    NRx_extra = NRx_extra + 1;
                    Receivers(NRx+NRx_extra) = Dipole(theta,phi);
                endif
                ! backward direction
                theta = (180.-Transmitters_Comp(I)%Theta);
                phi = mod(Transmitters_Comp(I)%Phi+180.,360.);
                call getTxRxIndex(NRx,Receivers(1:NRx),theta,phi,RxBks_exists);
                if (RxBks_exists .eq. 0) then ! add corresponding Rx
                    NRx_extra = NRx_extra + 1;
                    Receivers(NRx+NRx_extra) = Dipole(theta,phi);
                endif
            EndDo

            NRx_tot = NRx + NRx_extra;
            Allocate(Receivers_tmp(NRx_tot));
            Receivers_tmp(1:NRx_tot) = Receivers(1:NRx_tot);
            deallocate(Receivers); allocate(Receivers(NRx_tot));
            Receivers(1:NRx_tot) = Receivers_tmp(1:NRx_tot);
            deallocate(Receivers_tmp);
        EndIf

    Else
        NRx = NRxTheta*NRxPhi

        cosdth_init = cos(thiRecei*Pi/180.);
        margin_rec_theta =  cos(thfRecei*Pi/180.)-cos(thiRecei*Pi/180.);
        margin_rec_phi =  (phfRecei - phiRecei);
        ! step_theta_Recei
        if (NRxTheta ==1) then
            step_theta_Recei = 0;
        else
            step_theta_Recei = margin_rec_theta/(NRxTheta-1);
        endif
        ! step_phi_Recei
        if (NRxPhi == 1) then
            step_phi_Recei = 0.;
        else
            step_phi_Recei = margin_rec_phi/(NRxPhi-1);
        endif

        NRx_tot = NRx + 2*NTr ! temporarily
        Allocate(Receivers(NRx_tot))
        Ind = 1;
        theta = 0
        phi = 0
        DO I=1, NRxPhi
            phi = phiRecei + (I-1) * step_phi_Recei
            DO K=1, NRxTheta
                theta = acos(1.*(K-1)*step_theta_Recei+cosdth_init)*180./Pi;

                Receivers(Ind) = Dipole(theta,phi)
                Ind = Ind + 1;
            Enddo
        Enddo

        !! add extra-scatterers needed to calculate Cext and Cbks for all
        ! the available Transmitters/Incident directions
        NRx_extra = 0;
        Do I =1, NTr
            ! Forward direction
            theta = Transmitters_Comp(I)%Theta;
            phi = Transmitters_Comp(I)%Phi;
            call getTxRxIndex(NRx,Receivers(1:NRx),theta,phi,RxExt_exists);
            if (RxExt_exists .eq. 0) then ! add corresponding Rx
                NRx_extra = NRx_extra + 1;
                Receivers(NRx+NRx_extra) = Dipole(theta,phi);
            endif
            ! backward direction
            theta = (180.-Transmitters_Comp(I)%Theta);
            phi = mod(Transmitters_Comp(I)%Phi+180.,360.);
            call getTxRxIndex(NRx,Receivers(1:NRx),theta,phi,RxBks_exists);
            if (RxBks_exists .eq. 0) then ! add corresponding Rx
                NRx_extra = NRx_extra + 1;
                Receivers(NRx+NRx_extra) = Dipole(theta,phi);
            endif
        EndDo

        NRx_tot = NRx + NRx_extra;
        Allocate(Receivers_tmp(NRx_tot));
        Receivers_tmp(1:NRx_tot) = Receivers(1:NRx_tot);
        deallocate(Receivers); allocate(Receivers(NRx_tot));
        Receivers(1:NRx_tot) = Receivers_tmp(1:NRx_tot);
        deallocate(Receivers_tmp);
    EndIf

ENDSUBROUTINE get_trans_Receiv

subroutine getAngleIndex(ang,N,thetas,ind)

        USE Initialization

        Integer, INTENT(IN) :: N
        Real(kind=8),INTENT(IN) :: ang
        Real(kind=8),dimension(N), INTENT(IN) :: thetas
        Integer, INTENT(OUT) :: ind

        Integer ii;
        Real(kind=8),Parameter :: ErrAngEps = 1e-2

        Do ii = 1,N
            if (abs(thetas(ii)-ang) .le. ErrAngEps) then
                ind = ii;
                Exit;
            endif
        EndDo
        if (ii .gt. N) then
            Write(*,'(a)') 'Error when calculating backscattering cross section, bkw direction angle not found !'
            ind = 0;
        endif
    endsubroutine getAngleIndex

    subroutine getTxRxIndex(N,TxRx,th,ph,ind)

        USE Initialization

        Integer, INTENT(IN) :: N
        type (Dipole), Dimension(N), INTENT(IN) :: TxRx
        Real(kind=8),INTENT(IN) :: th, ph
        Integer, INTENT(OUT) :: ind

        Integer ii;
        Real(kind=8),Parameter :: ErrAngEps = 1e-2

        Do ii = 1,N
            if (abs(TxRx(ii)%theta-th) .le. ErrAngEps) then
                if (abs(TxRx(ii)%phi-ph) .le. ErrAngEps) then
                    ind = ii;
                    Exit;
                endif
            endif
        EndDo
        if (ii .gt. N) then
            ind = 0;
        endif
    endsubroutine getTxRxIndex
