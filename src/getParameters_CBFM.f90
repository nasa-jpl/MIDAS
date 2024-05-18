    SUBROUTINE getTransmitters_CBFM(Transmitters_CBFM);

    USE Initialization
    USE common_variables
    USE iso_fortran_env
    IMPLICIT NONE

    !! IN/OUT ******************************************************************
    type (Dipole), Dimension(:), allocatable, INTENT(OUT) :: Transmitters_CBFM

    ! Local
    Integer :: Ind,I,K
    Integer :: order,Npts
    Real(kind=8) :: step_theta_CBFM,step_phi_CBFM,theta_dipole,phi_dipole
    Real(kind=8) :: cosdth_init,margin_th,margin_ph,a,x_l,y_l,z_l
    Real(kind=8) :: th_init,th_end,ph_init,ph_end,delta
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

    th_init = 0.;th_end = 180.
    ph_init = 0.;ph_end = 360.

    if (distr_ipws == 1) then
        ! step_theta_CBFM, step_phi_CBFM
        !(360/step+1)*(180/step+1) = Nipws
        ! --> (Nipws-1)*step**2-(360+180)*step-360*180 = 0
        ! delta = b**2-4ac and solution = (-b+sqrt(delta))/2a
        ! ou aussi 2*k^2+3*k+1-N = 0 si k =180/step
        delta = (360+180)**2.+4*(Nipws-1)*(360*180);
        step_theta_CBFM = ((360+180)+sqrt(delta))/(2*(Nipws-1))
        step_phi_CBFM = step_theta_CBFM;
       
        ! finalement ca revient a :
        !step_theta_CBFM = 720./(sqrt(8.*Nipws+1.)-3.)
        !step_phi_CBFM = 720./(sqrt(8.*Nipws+1.)-3.)
        NTrTheta_CBFM = nint(180./step_theta_CBFM)+1
        NTrPhi_CBFM = nint(360./step_phi_CBFM)+1
        NTr_CBFM =  NTrTheta_CBFM*NTrPhi_CBFM
        Nipws = NTr_CBFM;
    else
        NTr_CBFM = Nipws
    endif

    If (distr_ipws == 1) Then !! uniform step in cos(theta) and Phi
        margin_ph = ph_end - ph_init;
        step_phi_CBFM = margin_ph/(NTrPhi_CBFM-1);

        ! here uniform distribution on with a constant step on cosTheta and phi
        margin_th = cos(th_end*Pi/180.)-cos(th_init*Pi/180.)
        step_theta_CBFM = margin_th/(NTrTheta_CBFM-1)  ! here we discretize cos Theta between 0 and pi
        cosdth_init = cos(th_init*Pi/180.);

        Allocate(Transmitters_CBFM(NTr_CBFM));
        Ind = 1;
        theta_dipole = 0.
        phi_dipole = 0.
        DO I=1, NTrPhi_CBFM
            phi_dipole = (I-1) * step_phi_CBFM
            DO K=1, NTrTheta_CBFM
                theta_dipole = acos(1.*(K-1)*step_theta_CBFM+cosdth_init)*180./Pi;
                Transmitters_CBFM(Ind) = Dipole(theta_dipole,phi_dipole)
                Ind = Ind + 1;
            Enddo
        Enddo
    elseif (distr_ipws == 2) Then ! random uniform in cos(theta) and phi
        call DATE_AND_TIME;
        call RANDOM_NUMBER(harvest=a);

        cosdth_init = cos(th_init*Pi/180.); ! cos(0)
        margin_th = cos(th_end*Pi/180.)-cos(th_init*Pi/180.); !(cos(180)-cos(0))
        margin_ph = ph_end-ph_init;

        Allocate(Transmitters_CBFM(NTr_CBFM));
        DO Ind=1, NTr_CBFM
          call RANDOM_NUMBER(harvest=a);
          theta_dipole = acos(margin_th*a+cosdth_init)*180./Pi;
          call RANDOM_NUMBER(harvest=a);
          phi_dipole = ph_init+a*margin_ph;

          Transmitters_CBFM(Ind) = Dipole(theta_dipole,phi_dipole)
        Enddo
    elseif ((distr_ipws == 3) .OR. (distr_ipws == 4)) Then ! spherical T-design (x) or Lebedev quad points
        if (distr_ipws == 3) then
            !call design_closer_order (NTr_CBFM,order,Npts); ! update 2/23/2024
            call ESD_design_closer_order (NTr_CBFM,order,Npts);
            NTr_CBFM = Npts;
            Nipws = Npts;
            Allocate(XYZstd(3,Npts));
            call ESD_design_points(order,Npts,XYZstd)
        else
            call lb_get_closer_Npts(NTr_CBFM,Npts); !
            call ld_by_order (Npts,x,y,z,w) ;
            Allocate(XYZstd(3,Npts));
            XYZstd(1,1:Npts) = x(1:Npts) ;
            XYZstd(2,1:Npts) = y(1:Npts) ;
            XYZstd(3,1:Npts) = z(1:Npts) ;
        endif

        Allocate(hypotxy(Npts));
        ! convert cartesian to spherical ! while taking into account the propagation along x of the incident direction (theta_r and phi_r are defined % x not z)
        hypotxy = hypot(XYZstd(2,:),XYZstd(3,:)); !hypotxy = sqrt(abs(x).^2 + abs(y).^2)

        Allocate(thetas_rd(Npts),phis_rd(Npts));
        thetas_rd = Pi/2.-atan2(XYZstd(1,:),hypotxy);
        phis_rd = atan2(XYZstd(3,:),XYZstd(2,:)); ! az = atan2(y,x);
        
        Allocate(Transmitters_CBFM(NTr_CBFM));
        DO Ind=1, NTr_CBFM
            theta_dipole = thetas_rd(Ind)/Pi*180.;
            phi_dipole = phis_rd(Ind)/Pi*180.
            Transmitters_CBFM(Ind) = Dipole(theta_dipole,phi_dipole)
        Enddo
        deallocate(thetas_rd,phis_rd,XYZstd,hypotxy)
    endif

    ENDSUBROUTINE getTransmitters_CBFM


    SUBROUTINE setfSR(nb,size,Cells_Block,vrb,fct_SR_blk,nnz);

    USE Initialization
    USE common_variables

    IMPLICIT NONE

    !! IN/OUT ******************************************************************
    Integer, INTENT(IN) :: nb,size,vrb
    type(Cell), Dimension(size), INTENT(IN) :: Cells_Block
    Real(kind=8), INTENT(OUT) :: fct_SR_blk
    Integer, INTENT(OUT) :: nnz


    !! local
    Integer :: N,cc,dd,jj,NiterSR,spr_size
    Real(kind=8) :: spr_perc,fSR_n,Beta_error,res,val_1,norm_Z,val_2,norm_Zspr,pr_res,curr_fSR,perc_spr

    ! arrays
    Integer, Dimension(:), allocatable :: row_sprZ,col_sprZ,col_sprZ_tmp
    Real(kind=8), Dimension(:), allocatable :: fSR_tests
    COMPLEX(real64), Dimension(:), allocatable :: Zpatch_e_spr,Zpatch_e_spr_tmp
    COMPLEX(real64), Dimension(:,:),allocatable :: Zpatch_e

    !! to start with, we save the current fct_SR
    !curr_fSR = fct_SR;

    NiterSR = 11 ;
    Allocate(fSR_tests(NiterSR));
    fSR_tests = [1e2,1e3,2e3,3e3,4e3,5e3,6e3,7e3,8e3,9e3,1e4];

    ! Frobinues Norm of the full matrix
    call Green_s_tr_partial_FN(size,Cells_Block,size,Cells_Block,norm_Z)
    jj =1; pr_res = -1; res = 1e3;
    do while ((jj .le. NiterSR) .and. (res .gt. res_SR)) ! exit if the differnce % FN(Z) lower than 0.01 %
        pr_res = res ;
        fct_SR_blk = fSR_tests(jj);

        !spr_perc = 100; ! let's say that we keep spr_perc % of the initial matrix Zii
                           ! Remember also that you're keeping only the upper part of the matrix
        !spr_size = nint((spr_perc*9.*size**2.)/100.)
        !nnz = 0;  ! to tell the code that I don't know nnz yet, and I need to know it here to avoid reallocationg later when calculating the CBFs
        !Allocate(Zpatch_e_spr(spr_size),row_sprZ(3*size+1),col_sprZ(spr_size));
        !Call SR_Green_s_tr_partial(size,Cells_Block,fct_SR_blk,spr_size,nnz,Zpatch_e_spr,row_sprZ,col_sprZ);

        Call SR_Green_s_tr_partial_FN(size,Cells_Block,fct_SR_blk,nnz,norm_Zspr);

        !!! FN Zspr
        !val_2 = 0.
        !Do cc =1,nnz
        !        val_2 = val_2 + (abs(Zpatch_e_spr(cc)))**2.
        !EndDo
        !!if (homogs == 1) then
        !!    val_2 = 2*val_2; ! don't worry the FN of the full Zii for homo = 1 was calculated the same way !
        !!endif
        !norm_Zspr = sqrt(val_2)
        !Deallocate(Zpatch_e_spr,col_sprZ,row_sprZ)

        res = abs(norm_Z-norm_Zspr)/norm_Z*100.
        perc_spr = (100.*nnz)/(9.*size**2.);
        if (vrb == 1 .and. rank ==0) then
            Write(*,*) 'nnz= ',nnz
            Write(*,'(a,ES7.1E1,a,f5.2,a,f7.2,a,ES7.1E1,a)') 'fSR = ',fct_SR_blk,' : with ',perc_spr,&
            ' % of Zii, NormZspr = ',norm_Zspr,'; Res = ',res,' %'
        endif

        jj = jj + 1;
    EndDo

    if (vrb .eq. 1 .and. rank ==0) then
        Write(*,'(a)') ' '
    endif

    ENDSUBROUTINE setfSR

    SUBROUTINE setNipws(nb,size,Cells_Block,vrb);

    USE Initialization
    USE common_variables

    IMPLICIT NONE

    !! IN/OUT ******************************************************************
    Integer, INTENT(IN) :: nb,size,vrb
    type(Cell), Dimension(size), INTENT(IN) :: Cells_Block

    !! local
    Integer :: jj,spr_size,nnz,mtype,iparm3,Nipws_init
    Integer :: pr_rank,pr_Nipws,NiterNpw,Nmax,nb_iter_out,curr_Nipws
    Real(kind=8) :: spr_perc,res

    ! arrays
    Integer, Dimension(:), allocatable :: row_sprZ,col_sprZ,col_sprZ_tmp
    Integer, Dimension(:), allocatable :: Nipws_tests

    type (Dipole), Dimension(:),allocatable :: Transmitters_CBFM
    COMPLEX(real64), Dimension(:), allocatable :: Zpatch_e_spr,Zpatch_e_spr_tmp
    COMPLEX(real64), Dimension(:,:),allocatable :: Zpatch_e, EREFpatch_e, Epatch_e
    COMPLEX(real64), Dimension(:,:), allocatable :: Matrix_U, Matrix_V

    INTERFACE
        SUBROUTINE getTransmitters_CBFM(Transmitters_CBFM)
            USE Initialization
            USE common_variables
            IMPLICIT NONE
            !! IN/OUT ******************************************************************
            type (Dipole), Dimension(:), allocatable, INTENT(OUT) :: Transmitters_CBFM

        END SUBROUTINE getTransmitters_CBFM
    END INTERFACE

    if (vrb .eq. 1 .and. rank ==0) then
        Write(*,'(a)') '-> Set Nipws :'
    endif

    spr_perc = 100; spr_size = nint((spr_perc*9.*size**2.)/100.)
    Allocate(Zpatch_e_spr_tmp(spr_size),col_sprZ_tmp(spr_size));
    Allocate(row_sprZ(3*size+1));
    ! here I will use the maximum fct_SR
    Call SR_Green_s_tr_partial(size,Cells_Block,fct_SR,spr_size,nnz,Zpatch_e_spr_tmp,row_sprZ,col_sprZ_tmp)

    Allocate(Zpatch_e_spr(nnz),col_sprZ(nnz));
    Zpatch_e_spr(1:nnz) = Zpatch_e_spr_tmp(1:nnz); col_sprZ(1:nnz) = col_sprZ_tmp(1:nnz)
    Deallocate(Zpatch_e_spr_tmp, col_sprZ_tmp)

    curr_Nipws = Nipws;

    If (homogs == 1) Then !
      mtype = 6 !complex symmetric matrix
    Else
      mtype = 13 !complex nonsymmetric matrix
    EndIf

    nb_iter_out = 0;
    Nipws=0;
    NiterNpw = 10; !8
    ! I prefer this approach (%2*Nipws from Chen & al URSI 2017) because Nipws increases more slowly
    if (distr_ipws == 3) then ! ici cas particulier de sphere_design_rule
        allocate(Nipws_tests(NiterNpw));
        Nipws_tests = [84,94,108,120,144,156,180,204,216,240];! commencons ainsi avant d'avoir les valeurs superieurs de N
    elseif (distr_ipws == 4) then
        allocate(Nipws_tests(59));
        Nipws_tests = [86,110,146,170,194,230,266,302,350,386,434,482,530,590, &
         650,  698,  770,  830,  890,  974, 1046, 1118, 1202, 1274, &
        1358, 1454, 1538, 1622, 1730, 1814, 1910, 2030, 2126, 2222, &
        2354, 2450, 2558, 2702, 2810, 2930, 3074, 3182, 3314, 3470, &
        3590, 3722, 3890, 4010, 4154, 4334, 4466, 4610, 4802, 4934, &
        5090, 5294, 5438, 5606, 5810];
    else
        allocate(Nipws_tests(NiterNpw));
        Nipws_tests = [91,190,231,325,496,703,861,1225,1891,2701]
    endif
    !Nipws_init = 32;
    jj= 1; res = 1e3;

    do while ((jj .le. NiterNpw) .and. (floor(res) .gt. 1)) ! if the difference is larger than 1 %
        ! previous
        pr_Nipws = Nipws;
        pr_rank = nb_iter_out;

        Nipws = Nipws_tests(jj)  !!Nipws_init*2**jj

        call getTransmitters_CBFM(Transmitters_CBFM);

        !! Incident field used to compute the CBFs (different from the scattering problem incident field )
        Allocate(EREFpatch_e(3*size,2*Nipws))
        Call Incident_Field(1,size,Cells_Block,Transmitters_CBFM,Nipws,EREFpatch_e);

        Allocate(Epatch_e(3*size,2*Nipws))
        iparm3 = 0; ! here iparm3 is not used
        Call pardiso_solver(3*size,2*Nipws,nnz,mtype,iparm3,row_sprZ,col_sprZ,Zpatch_e_spr,EREFpatch_e,Epatch_e)

        Nmax = 3*size ;
        Allocate(Matrix_U(3*size,Nmax), Matrix_V(Nmax,2*Nipws));
        CALL Calcul_GenMatrix_ACA(3*size,2*Nipws,Epatch_e,Nmax,nb_iter_out,Matrix_U,Matrix_V);
        deallocate(Epatch_e,Matrix_U,Matrix_V);

        res = real(nb_iter_out-pr_rank)/real(Nipws-pr_Nipws)*100.

        if (vrb .eq. 1) then
            Write(*,'(a,i4,a,i3,a,f7.2,a)') 'Nipws = ',Nipws,'; rank = ',nb_iter_out,'; Res = ',res,' %'
        endif

        jj = jj + 1;
        deallocate(Transmitters_CBFM, EREFpatch_e);
    enddo

    Nipws = Nipws_tests(jj-2);
    !Nipws = pr_Nipws ;

    if (curr_Nipws .gt. Nipws) then !! in this case, get back our current Nipws determined from the previous considered block
        Nipws = curr_Nipws;
    EndIf

    if (vrb .eq. 1) then
        Write(*,'(a)') ' '
    endif
    ENDSUBROUTINE setNipws

    SUBROUTINE initializeNipws(SimScatterer)

        USE Initialization
        USE common_variables

        IMPLICIT NONE

        !IN/OUT
        type (Scatterer), INTENT(INOUT) :: SimScatterer

        ! local
        integer :: rr
        real(kind=8) :: r_lambda

        ! here decide Nipws for the generation of CBFs depending on hmax/lambda_s
        if ((CBFM .NE. 0) .OR. (MLCBFM .NE. 0)) Then
            If (set_Nipws==0) then
            r_lambda = (hBlock/2.)/SimScatterer%lambda_min
            If (distr_ipws .eq. 3) then   ! spherical design
                Do rr=1,8
                If (r_lambda .le. rr) Then
                    Nipws = SphDes_Nipws_f_rlamb(rr);
                    Exit;
                EndIf
                EndDo
            ElseIf (distr_ipws .eq. 4) then   ! Lebedev quadrature
                Do rr=1,8
                If (r_lambda .le. rr) Then
                    Nipws = LebQuad_Nipws_f_rlamb(rr);
                    Exit;
                EndIf
                EndDo
            Else ! for the other distributions use Nipws_f_rlamb coming from my PhD
                Do rr=1,8
                If (r_lambda .le. rr) Then
                    Nipws = Nipws_f_rlamb(rr);
                    Exit;
                EndIf
                EndDo
            EndIf
            If (rr .eq. 9) Then
                if (rank .eq. 0) then
                Write(*,'(a,f8.4,a)') 'WARNING : Problem when selecting Nipws :  r0/lambda =',r_lambda,'  > 8'
                endif
                Nipws = Nipws_f_rlamb(8);
                Nipws = Nipws_f_rlamb(4);
                !Stop 1;
            EndIf
            Else
            Nipws = 2; ! Initialization
            EndIf
        EndIf
    END SUBROUTINE initializeNipws
