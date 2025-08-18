SUBROUTINE Compute_EFields_CBFME(Cells,CBFM_Blocks,CBFM_Blocks_Ext,MPI_CBFM_Blocks,Transmitters,K_patchs,C_job_patchs,E_total)
    ! Modifs 12/10/2023 : check RCOND and ERROR BOUND using SCALAPACK PZGECON; a commented extra-lines of code can be used to do the same with the expert
    ! Modif 8/16/2025 : adding CBFM threshold and Nipws selection as function of the dielectric properties 
    USE Initialization
    USE common_variables
    USE lapack95
    USE DiverseUtil
    USE MPI

    IMPLICIT NONE

    !! IN/OUT ******************************************************************
    type (Cell), Dimension(Nbc), INTENT(IN) :: Cells
    type (CBFM_Block), Dimension(Nblocks), INTENT(IN) :: CBFM_Blocks
    Integer, Dimension(Nblocks,Nbc_ext), INTENT(IN) :: CBFM_Blocks_Ext
    !Integer, Dimension(2+Nccp_max,NcalBlks), INTENT(IN) :: cp_CBFM_Blocks
    Integer, Dimension(nber_procs,Nblk_proc_max+1), INTENT(IN):: MPI_CBFM_Blocks
    type (Dipole), Dimension(NTr), INTENT(IN) ::Transmitters
    Integer, Dimension(:), allocatable, INTENT(INOUT) :: K_patchs
    COMPLEX(real64), Dimension(:,:), allocatable, INTENT(INOUT) :: C_job_patchs
    COMPLEX(real64), Dimension(:,:), allocatable, INTENT(OUT):: E_total

    !! Local MPI
    INTEGER :: fh,NBlocks_jobii, curs, ii_job, kk_job,prev_kk,Eff_rank,onlyZ
    INTEGER :: shift0,shift1,shift2,bytes_in_dp,bytes_in_integer, nb_values
    INTEGER(kind= MPI_OFFSET_KIND) :: offset_file
    INTEGER, dimension(MPI_STATUS_SIZE) :: status
    Integer, Dimension(:), allocatable :: blocks_to_sort,OrderBlocksProcs,Ktot_procs,Nbc_proc_tot
    integer, dimension(:), allocatable :: values
    Real(kind=8) :: e_re,e_im

    !! Scalapack :
    Integer :: icontxt,ICSRC,IRSRC,IA,JA,IB,JB
    INTEGER :: NbLinesPrev,Ktot_proc_max, size2, Nkk
    INTEGER :: INDXL2G,INDXG2L,iloc,jloc,iLoc_proc,iGlob,jLoc_proc,jGlob
    Integer :: block_iGlob, pos_loc_block,NBlocks_eff
    Integer :: LWORK, LIWORK
    INTEGER, EXTERNAL :: numroc
    REAL(kind=8), EXTERNAL :: pzlange
    REAL, EXTERNAL :: pslamch
    INTEGER, parameter :: ROW_SRC = 0, COL_SRC = 0
    INTEGER, parameter :: NDIMS = 2
    INTEGER,dimension(1:NDIMS) :: dims
    INTEGER, dimension(:), allocatable :: Blocks_eff, K_patchs_eff, Curs_Kpatchs_all
    INTEGER, dimension(:), allocatable :: DESCA, DESCB, IPIV
    DOUBLE COMPLEX, Dimension(:,:), allocatable :: ZredLoc,VredLoc!,ZFredLoc
    DOUBLE COMPLEX, Dimension(:,:), allocatable :: AlphaProc!,AlphaLoc
    Real(kind=8) :: EPSMCH,ANORM,ERRBD
    Real(kind=8), Dimension(:), allocatable :: WORK!,R,C,FERR,BERR
    Integer, Dimension(:), allocatable  :: IWORK
    !CHARACTER :: EQUED

    !! LOCAL
    Integer :: I,id,nthreads,p,d,kk,ii,jj,dd,pp,klu_cel,klu,CBFs_for_all_freqs
    Integer :: size,curs_lig_E,M,N,K,Diag_Representation,Nbc_b,Ncps
    Integer :: Nbc_b_ext, K_total,tdis,ii1,ii2,m_max_relative
    Integer :: size1,size3,NbreLig_mat1,NbreLig_mat3,NbreCol_mat1,NbreCol_mat3
    Integer :: curseurf,BlockSize,num_cel_fichier,num_fichier
    Integer :: nb,mtype,iparm3,pack_size,ii_beg,ii_end,jj_beg,jj_end,fu_ii,fu_jj,pack_ii,loc_ii,loc_jj
    Integer :: SMWA,lev_SMW,no_aca,sp_band_thresh,ii_s,jj_s
    Integer :: NiterNpw,NiterSR,Nipw_n,Ntests,vrb_cbfm_param,track_time_Zc
    Integer :: Kmax_job,K_proc, jj_job, K_max, Nbc_proc_max, sizeTrans
    Integer :: prev_job, next_job,rank_EffNextjob,NBlocks_EffNextjob,jj_Nextjob,prev_jj
    Integer :: spr_size,nnz,nnz_elts_Zc

    Real(kind=8) :: Threshold_CBFM,RCOND,norme,S_val,z_max,spr_perc,spr_perc_loc,CR,perc_ii
    Real(kind=8) :: f_Zc,locThresh_Zc,bandThresh_Zc,abs_elt,spr_ii,fct_SR_blk
    character(8)  :: date_init_N1, date_final_N1, date_init_ii, date_final_ii
    character(10) :: time_init_N1, time_final_N1, time_init_ii, time_final_ii
    character(5)  :: zone_init_N1, zone_final_N1, zone_init_ii, zone_final_ii

    ! track Zc time
    Integer :: num_file_time_Zc
    character(300) :: analysis_fold_name,timeZc_file_name
    character(:), allocatable :: rank_str
    character(8)  :: date_init_Zc, date_final_Zc
    character(10) :: time_init_Zc, time_final_Zc
    character(5)  :: zone_init_Zc, zone_final_Zc
    Integer,dimension(8) :: values_init_Zc, values_final_Zc
    Integer, dimension(4):: time_calcul_Zc

    character(200) :: file_name, Zcfilename
    CHARACTER(3), allocatable::number_chars(:)
    CHARACTER(15) :: file_name_Z
    CHARACTER(3) :: kk_st
    COMPLEX(real64) :: alpha, beta
    ! vectors and matrices
    Integer,dimension(8) :: values_init_N1, values_final_N1, values_init_ii, values_final_ii
    Integer, dimension(4):: time_calcul_N1, time_calcul_ii
    Integer, Dimension(:), allocatable :: nnz_blocks,curs_B_cel,curs_B_Cpatch,K_patchs_all
    Integer, Dimension(:), allocatable :: curs_row_Z_blocks, curs_col_Z_blocks, curs_B_Cpatch_Next
    Integer, Dimension(:), allocatable :: row_sprZ,col_sprZ,testedBlks,ExtSizes
    Integer, Dimension(:), allocatable :: nb_elements_rec,deplts,vect_tmp
    type (Dipole), Dimension(:),allocatable :: Transmitters_CBFM
    type(Cell), Dimension(:), allocatable :: Cells_Block,Cells_Block_ii,Cells_Block_jj
    Real(kind=8), Dimension(:), allocatable :: spr_perc_blocks, fSR_blocks, S, WW,vals
    Real(kind=8),Dimension(:,:),allocatable :: abs_Zpatch_e
    COMPLEX(real64), Dimension(:),allocatable :: Zpatch_e_spr,Vect,VectProduit,Zred_pack
    COMPLEX(real64), Dimension(:,:),allocatable :: Zpatch_e, EREFpatch_e, Epatch_e, Cpatch_e
    COMPLEX(real64), Dimension(:,:), allocatable :: AA, U, VT, C_job_patchs_tmp, C_trans_patchs
    COMPLEX(real64), Dimension(:,:), allocatable :: Matrice1, Matrice2,Matrice3,Mat_Inter,MatProduit,Vect_Inter
    COMPLEX(real64), Dimension(:,:), allocatable :: Zreduite, E_ref_incident,Vreduit, Alphas,Uin,Uout

    ! test 1/9/24
    CHARACTER(200) :: file_name_svd
    CHARACTER(10) :: k_st

    INTERFACE
        SUBROUTINE getTransmitters_CBFM(Transmitters_CBFM)
            USE Initialization
            USE common_variables
            IMPLICIT NONE
            !! IN/OUT ******************************************************************
            type (Dipole), Dimension(:), allocatable, INTENT(OUT) :: Transmitters_CBFM
        END SUBROUTINE getTransmitters_CBFM
    END INTERFACE

    SMWA = 0;

    If (rank == 0) Then
        Write (*,*) ''
        Write (*,'(a)') '>>>>>>>>>>>>>>>>>>>>>>>>>>> 1-Level CBFM-E <<<<<<<<<<<<<<<<<<<<<<<<<<<'
        Write (*,*) ''; Write (*,*) ''
        Write (*,'(a)') '----------------- Total Fields Inside Scatterer ---------------------'
        Write (*,*) ''
    EndIf

    ! just Initialization, the used values will be set through setfSR and setNipws
    fct_SR = 1e2;
    ! to show extra information
    vrb_cbfm_param = 0;
    track_time_Zc = 0; !Nt_display = 20; ! show time for each step in the generation of Zc for 'Nt_display' MPI tasks   !! IL VAUT MIEUX PEUT ETRE ECRIRE UN FICHIER QUE JE POURRAI LIRE AVEC MATLAB
    
    If (homogs == 1) Then ! mtype for PRADISO used later to solve the sparse Zii and Zc
      mtype = 6 !complex symmetric matrix
    Else
      mtype = 13 !complex nonsymmetric matrix
    EndIf
    if (debug_mode .eq. 0) then ! automatically
        vrb_cbfm_param = 0; track_time_Zc = 0;
    endif
    ! HERE GET MY NBlocks ! attention to the difference with NBlocks_job that can use for any other job
    ! MyNBlocks is the NBlocks_job of the current job
    MyNBlocks = MPI_CBFM_Blocks(rank+1,1);

    !! ALLOCATION AND PREPARATION HERE
    ! these two vectors will be used both for the genration of the CBFs and of the reduced matrix
    Allocate(curs_B_cel(NBlocks),curs_B_Cpatch(MyNBlocks),K_patchs_all(NBlocks))
    ! curs_B_Cpatch for the MyNBlocks calculated by the current process
    curs_B_Cpatch(1)=1
    DO kk_job=2, MyNBlocks
        prev_kk = MPI_CBFM_Blocks(rank+1,kk_job);  ! equivalent to kk-1 1+kk_job-1
        curs_B_Cpatch(kk_job) = curs_B_Cpatch(kk_job-1) + 3*CBFM_Blocks(prev_kk)%Nbc_b
    EndDo
    ! curs_B_cel for all the NBlocks blocks
    curs_B_cel(1)=1
    DO kk=2, NBlocks
        curs_B_cel(kk) = curs_B_cel(kk-1) + CBFM_Blocks(kk-1)%Nbc_b
    EndDo

    ! the table OrderBlocksProcs generated here is usefull to gather informations from procs
    ! it shows the global order of the blocks. Example : OrderBlocksProcs = [1,4,2,3]
    !  with 2 jobs per proc <-> job0 has b1 and b3 and job1 has 4 and 2
    Allocate(blocks_to_sort(NBlocks),OrderBlocksProcs(Nblocks));
    curs = 1;
    Do ii = 1, nber_procs
      NBlocks_jobii = MPI_CBFM_Blocks(ii,1);
      blocks_to_sort(curs:curs+NBlocks_jobii-1)= MPI_CBFM_Blocks(ii,2:1+NBlocks_jobii);
      curs = curs + NBlocks_jobii;
    EndDo
    OrderBlocksProcs = (/(ii, ii=1,NBlocks)/);
    call sort_asc(blocks_to_sort,NBlocks,OrderBlocksProcs);
    Deallocate(blocks_to_sort);

    ! Added 6/22 to force the code to calculate CBFs for all freqs : need to be improved later (input by user maybe, or determined % frequency range/step)
    CBFs_for_all_freqs = 1
    if ((CBFs_for_all_freqs .eq. 1) .OR. (num_freq .eq. 1)) then
    if (num_freq .gt. 1) then 
      deallocate(K_patchs,C_job_patchs);
    endif
    Allocate(K_patchs(MyNBlocks))
    !! Here we start by defining the adequate fSR
    !! here we define the tested blocks for which we will determine adequate fSR and Nipws.
    !! The calculated fSR and Nipws for these blocks will be then applied to all the other blocks
    Ntests = 4; Allocate(testedBlks(Ntests))
    testedBlks(1) = maxloc(CBFM_Blocks(:)%BSphCont(4),1);
    testedBlks(2) = minloc(CBFM_Blocks(:)%BSphCont(4),1);
    Allocate(ExtSizes(NBlocks));
    ExtSizes(:) = CBFM_Blocks(:)%Nbc_b + CBFM_Blocks(:)%Nbc_ext
    testedBlks(3) = maxloc(ExtSizes(:),1);
    testedBlks(4) = minloc(ExtSizes(:),1);
    deallocate(ExtSizes);
    if (vrb_cbfm_param == 1 .and. rank == 0) then
      Write(*,'(a)') ' '
      Write(*,'(a)') 'Characteristic blocks : [maxh, minh, maxNbc, minNbc]'
      Write(*,*) testedBlks(1:4)
      Write(*,'(a)') ' '
    endif

    call MPI_BARRIER(MPI_COMM_WORLD,code);
    ! SET fSR FOR EACH BLOCK
    Allocate(fSR_blocks(MyNBlocks),nnz_blocks(MyNBlocks),spr_perc_blocks(MyNBlocks));
    fSR_blocks = 0.;
    nnz_blocks = 0;
    time_calcul_N1 = 0
    call date_and_time(date_init_N1,time_init_N1,zone_init_N1,values_init_N1);
    If ((SR == 1) .and. (SMWA .ne. 1)) Then
        if (rank == 0 .and. debug_mode == 1) then
            Write(*,'(a)') ''; Write(*,'(a)') '-> Set fSR :'
        endif
        ! To ensure a higher accuracy, let us test all the blocks, I do not think that this has
        ! a huge impact on the cpu time but strongly helps to be totally sure about fSR and spr_perc
        do ii_job =1,MyNBlocks
            ii = MPI_CBFM_Blocks(rank+1,1+ii_job);
            Nbc_b = CBFM_Blocks(ii)%Nbc_b
            Nbc_b_ext = CBFM_Blocks(ii)%Nbc_ext
            size = Nbc_b + Nbc_b_ext
            Allocate(Cells_Block(size))
            if (debug_mode ==1 .and. vrb_cbfm_param ==0) Write(*,'(a,i4,a,i6)') 'Block ',ii,' of size',size;
            Cells_Block(1:Nbc_b) = Cells(curs_B_cel(ii):curs_B_cel(ii)+Nbc_b-1)
            Do kk=1,Nbc_b_ext
                Cells_Block(Nbc_b+kk) = Cells(CBFM_Blocks_Ext(ii,kk))
            EndDo
            call setfSR(ii,size,Cells_Block,vrb_cbfm_param,fct_SR_blk,nnz);
            fSR_blocks(ii_job) = fct_SR_blk;
            nnz_blocks(ii_job) = nnz
            spr_perc_blocks(ii_job) = (100.*nnz)/(9.*size**2.)
            deallocate(Cells_Block);
        enddo
    EndIf

    ! To use in case if the user would not test all the blocks for fSR and spr_ii (which is not recommended)
    spr_perc_loc = maxval(spr_perc_blocks);
    Call MPI_ALLREDUCE(spr_perc_loc,spr_perc,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,code)
    fct_SR = maxval(fSR_blocks);     ! this f_SR max is also used to calculate Nipws in the next step
    deallocate(spr_perc_blocks);

    ! SET Nipws FOR SELECTED BLOCKS
    if (set_Nipws .ne. 0) then
        if (rank == 0) then
            Write(*,'(a)') ''; Write(*,'(a)') '-> Set Nipws :'
        endif
        !! to avoid the expensive cpu time cost of the resolution + aca
        ! in the open mp versin we test the results obtained with the 4 characteristic blocks,
        !to reduce the cpu time we parallelize this paragraph
        ! let us here for the mpi version test all the blocks to start with !
        ! we can later for example re-define these 4 blocks per job (or globally maybe since the job is not seeing what the other jobs got as Nipws ??!!)
        do ii =1,MyNBlocks !4
            !! lb1 ***************************************
            !nb = testedBlks(ii)
            nb = MPI_CBFM_Blocks(rank+1,1+ii_job);
            Nbc_b = CBFM_Blocks(nb)%Nbc_b
            Nbc_b_ext = CBFM_Blocks(nb)%Nbc_ext
            size = Nbc_b + Nbc_b_ext
            if (vrb_cbfm_param ==0) Write(*,'(a,i4,a,i6)') 'Block ',nb,' of size',size;
            Allocate(Cells_Block(size))
            Cells_Block(1:Nbc_b) = Cells(curs_B_cel(nb):curs_B_cel(nb)+Nbc_b-1)
            Do kk=1,Nbc_b_ext
                Cells_Block(Nbc_b+kk) = Cells(CBFM_Blocks_Ext(nb,kk))
            EndDo
            call setNipws(nb,size,Cells_Block,nnz_blocks(ii),vrb_cbfm_param);
            deallocate(Cells_Block);
        enddo
    endif

    call MPI_BARRIER(MPI_COMM_WORLD,code);
    If ((rank == 0) .and. (SR == 1) .and. (SMWA .ne. 1)) Then
        Write(*,'(a,ES7.1E1,a,f5.2)') ' -- > fSR for CBFM = ',fct_SR, ' -> spr % = ',spr_perc
    EndIf
    
    ! Threshold for the generation of the CBFs, decided as function of max(|m|)
    vals = abs(Cells(1:Nbc)%m_n) 
    m_max_relative = ceiling(maxval(vals)) - 2. ! we will add this difference to the power of the threshold
    Threshold_CBFM = 10.**(-3-m_max_relative)
    if (rank .eq. 0) then 
        write(*,'(a,ES7.1E1)') ' -- > Threshold_CBFM = ',Threshold_CBFM
    endif
    
    ! for distr_ipws = 1, try [91,190,231,325,496,703,861]; for distr_ipws=3 try among 289,366,482,579,723,842,926; for distr_ipws=4 try among [110,194,230,350,434,590,770,974,1202]
    distr_ipws = 3; Nipws = 482 ! fix it here for now until debugging 4/9/2024 
    call getTransmitters_CBFM(Transmitters_CBFM);
    if (rank == 0) then 
        if (distr_ipws .eq. 1) then 
            Write(*,'(a,i6,a)') ' -- > Nipws for CBFM = ',Nipws, ' (uniform in cos(theta) and phi)'
        elseif (distr_ipws .eq. 3) then 
            Write(*,'(a,i6,a)') ' -- > Nipws for CBFM = ',Nipws, ' (sd)'
        elseif (distr_ipws .eq. 4) then
            Write(*,'(a,i6,a)') ' -- > Nipws for CBFM = ',Nipws, ' (lb)'
        endif
        call date_and_time(date_final_N1,time_final_N1,zone_final_N1,values_final_N1)
        call Calcul_time_spent(values_init_N1,values_final_N1, time_calcul_N1)
        Write(*,'(a,i2,a,i2,a,i2,a,i2,a)') ' --> to set CBFM parameters : ',time_calcul_N1(1),'j',time_calcul_N1(2)&
        ,'h',time_calcul_N1(3),'min', time_calcul_N1(4),'sec'
        Write(*,'(a)') ''
    endif

    Allocate(C_job_patchs(3*Nbc_proc,2*NTr_CBFM));
    call print_allocate(35,'C_job_patchs(3*Nbc_proc,2*NTr_CBFM)','DCOMP',3*Nbc_proc*2*NTr_CBFM); !(Nchar,allocate_str,type_str,size)
    C_job_patchs = 0D0;

    call MPI_BARRIER(MPI_COMM_WORLD,code);
    ! Here generation of the CBFS *************************************************************************************************
    if (rank == 0 ) Then
        if (debug_mode .eq. 1) then 
            Write(10,*) '';Write(10,'(a,e14.5)') 'Singular Values generated by the SVD with a threshold equal to ', Threshold_CBFM
        endif
        Write(*,*) ''
        If (DR == 1) Then
            Write(*,'(a)') 'Calculation of the CBFs (+DR):'
        ElseIf (SMWA == 1) Then
            Write(*,'(a)',advance='no') 'Calculation of the CBFs (+SMWF'
            If (SR ==1) Then
                Write(*,'(a)') '+SR):'
            else
                Write(*,'(a)') '):'
            endif
        ElseIf (SR ==1) Then
            if (homogs ==1) then
                Write(*,'(a,ES7.1E1,a)') 'Calculation of the CBFs (+SR-PardisoSym; fct_SR_max = ',fct_SR,'):'
            else
                Write(*,'(a,ES7.1E1,a)') 'Calculation of the CBFs (+SR-PardisoNonSym; fct_SR_max = ',fct_SR,'):'
            endif
        Else
            Write(*,'(a)') 'Calculation of the CBFs :'
        EndIf
        time_calcul_N1 = 0
        call date_and_time(date_init_N1,time_init_N1,zone_init_N1,values_init_N1)
    endif

    ! START CBFS HERE
    call MPI_BARRIER(MPI_COMM_WORLD,code);
    curs_lig_E = 1; K_patchs= 0;
    DO kk_job=1,MyNBlocks
        kk = MPI_CBFM_Blocks(rank+1,1+kk_job);
	      ! define the cells belonging to the extended block kk
        Nbc_b = CBFM_Blocks(kk)%Nbc_b
        Nbc_b_ext = CBFM_Blocks(kk)%Nbc_ext
        size = Nbc_b + Nbc_b_ext
        Allocate(Cells_Block(size))
        Cells_Block(1:Nbc_b) = Cells(curs_B_cel(kk):curs_B_cel(kk)+Nbc_b-1)
        Do ii=1,Nbc_b_ext
            Cells_Block(Nbc_b+ii) = Cells(CBFM_Blocks_Ext(kk,ii))
        EndDo

        !! Incident field used to compute the CBFs (different from the scattering problem incident field )
        Allocate(EREFpatch_e(3*size,2*NTr_CBFM));
        call print_allocate(30,'EREFpatch_e(3*size,2*NTr_CBFM)','DCOMP',3*size*2*NTr_CBFM);

        Call Incident_Field(1,size,Cells_Block,NTr_CBFM,Transmitters_CBFM,1,NTr_CBFM,EREFpatch_e);
        !! Resolve locally the scattering problem (compute the electric field inside the block kk
        !! resulting from the CBFM incident field : EM plane waves from the entire space)
        Allocate(Epatch_e(3*size,2*NTr_CBFM));
        call print_allocate(27,'Epatch_e(3*size,2*NTr_CBFM)','DCOMP',3*size*2*NTr_CBFM);

	      time_calcul_ii = 0
        call date_and_time(date_init_ii,time_init_ii,zone_init_ii,values_init_ii);
        If (DR == 1) Then
            klu_cel= 0  ! finally since we're considering only Diagonal elements Klu_cel will always == 0
            klu = 3*(klu_cel+1)-1;
            Allocate(Zpatch_e(2*klu+1,3*size))
            Call DR_Green_s_tr_partial(size,Cells_Block,klu_cel,klu,Zpatch_e)
            Write(*,'(a,i4,a,i6,a,i6)') 'j',rank,': Solving blk ',kk,' of size',size
            call gbsvx(Zpatch_e,EREFpatch_e,Epatch_e,RCOND=RCOND)
            Deallocate(Zpatch_e)
        ElseIf ((SR == 1) .and. (SMWA .ne. 1)) Then
            fct_SR_blk = fSR_blocks(kk_job); nnz  = nnz_blocks(kk_job); spr_perc = 100.*nnz/(9*size**2.);
            Allocate(Zpatch_e_spr(nnz),row_sprZ(3*size+1),col_sprZ(nnz));
            call print_allocate(17,'Zpatch_e_spr(nnz)','DCOMP',nnz);
            call print_allocate(13,'col_sprZ(nnz)','SINTG',nnz);

            Call SR_Green_s_tr_partial(size,Cells_Block,fct_SR_blk,nnz,Zpatch_e_spr,row_sprZ,col_sprZ);
	        Write(*,'(a,i4,a,i6,a,i6,a,f5.2,a)',advance='no') 'j',rank,': Solving blk ',kk,' of size',size,': nnz = ',spr_perc,' % of Zii'

            iparm3 = 0; ! here iparm3 is not used
            Call pardiso_solver(3*size,2*NTr_CBFM,nnz,mtype,iparm3,row_sprZ,col_sprZ,Zpatch_e_spr,EREFpatch_e,Epatch_e)
            Deallocate(Zpatch_e_spr,row_sprZ,col_sprZ);
        ElseIf (SMWA == 1) Then
            lev_SMW = 1;
            Write(*,'(a,i4,a,i6,a,i6)') 'j',rank,': Solving blk ',kk,' of size',size
            Call SMWFB_algorithm(lev_SMW,Cells,Cells_Block,size,2*NTr_CBFM,EREFpatch_e,Epatch_e);
        Else
            Allocate(Zpatch_e(3*size,3*size))
            Call Green_s_tr_partial(size,Cells_Block,size,Cells_Block,Zpatch_e)
            Write(*,'(a,i4,a,i6,a,i6)') 'j',rank,': Solving blk ',kk,' of size',size
            Call gesvx(Zpatch_e,EREFpatch_e,Epatch_e,RCOND=RCOND)
            Deallocate(Zpatch_e)
        EndIf
        Deallocate(Cells_Block);
	    call date_and_time(date_final_ii,time_final_ii,zone_final_ii,values_final_ii)
        call Calcul_time_spent(values_init_ii,values_final_ii, time_calcul_ii)
        Write(*,'(a,i6,a)') ' : ',time_calcul_ii(1)*86400+time_calcul_ii(2)*3600+time_calcul_ii(3)*60+time_calcul_ii(4),' sec'

        !! Decomposition en valeurs singulieres
        M = 3*size; N = 2*NTr_CBFM
        Allocate(S(MIN(M,N)),U(M,M),VT(N,N), WW(MIN(M,N)-1))
        call print_allocate(16,'U(3*size,3*size)','DCOMP',3*size*3*size);
        call print_allocate(25,'VT(2*NTr_CBFM,2*NTr_CBFM)','DCOMP',2*NTr_CBFM*2*NTr_CBFM);

        CALL GESVD(A=Epatch_e,S=S,U=U, VT=VT, JOB='U')              

        !! Normalisation et comparaison au seuil, K designera le nombre de valeurs singulieres retenues (non nulles)
        K = 0
        norme = S(1)
        DO dd=1, MIN(M,N)
            S_val = S(dd)/norme
            If (S_val >= Threshold_CBFM) Then
                S(dd) = S_val
                K = K + 1
            Else
                S(dd) = 0
            EndIF
        ENDDO


        !! the local C_patch (K first columns of the matrix U)
        ! Fill C_tot_patchs from Cpatch_e
        ! we only take the CBFs corresponding to the original size of the block kk (not extended)
        C_job_patchs(curs_B_Cpatch(kk_job):curs_B_Cpatch(kk_job)+3*Nbc_b-1,1:K) = U(1:3*Nbc_b,1:K);
        K_patchs(kk_job) = K

        !! Deallocaton de tous les vecteurs propores au bloc (Interieur de la boucle)
        Deallocate(EREFpatch_e,Epatch_e)
        Deallocate(S,U,VT,WW)
    ENDDO

    Kmax_job = maxval(K_patchs);
    Do kk_job = 1, nber_procs !! re-allocating per job is all what I was able to do here to reduce memory use for the total node
        if ((rank-1) .eq. kk_job) then
            Allocate(C_job_patchs_tmp(3*Nbc_proc,Kmax_job)); C_job_patchs_tmp = 0D0;
            call print_allocate(37,'C_job_patchs_tmp(3*Nbc_proc,Kmax_job)','DCOMP',3*Nbc_proc*Kmax_job);

            C_job_patchs_tmp(1:3*Nbc_proc,1:Kmax_job) = C_job_patchs(1:3*Nbc_proc,1:Kmax_job);
            deallocate(C_job_patchs); Allocate(C_job_patchs(3*Nbc_proc,Kmax_job));
            C_job_patchs = C_job_patchs_tmp; deallocate(C_job_patchs_tmp);
        endif
        call MPI_BARRIER(MPI_COMM_WORLD,code);
    EndDo
    else
    	if (rank == 0) then
        Write(*,'(a)') 'Re-Using CBFs of previous SIM : '
        endif
        Kmax_job = maxval(K_patchs);
    endif


    K_proc = sum(K_patchs)
    Call MPI_ALLREDUCE(K_proc,K_total,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,code)
    CR = 3*Nbc/K_total

    ALLOCATE(nb_elements_rec(nber_procs), deplts(nber_procs))
    nb_elements_rec(1:nber_procs) = MPI_CBFM_Blocks(1:nber_procs,1)
    deplts(1)=0
    Do ii=2,nber_procs
      deplts(ii)= sum(MPI_CBFM_Blocks(1:ii-1,1));
    EndDo

    call MPI_ALLGATHERV (K_patchs,MyNBlocks, MPI_INTEGER ,K_patchs_all,nb_elements_rec,&
    deplts,MPI_INTEGER,MPI_COMM_WORLD ,code);
    ! Reorder the collected K_patchs_all using OrderBlocksProcs
    Allocate(vect_tmp(NBlocks))
    Do ii=1,NBlocks
      jj= OrderBlocksProcs(ii);
      vect_tmp(ii) = K_patchs_all(jj)
    EndDo
    K_patchs_all = vect_tmp;
    Deallocate(vect_tmp);
    Deallocate(OrderBlocksProcs);

    ! determine Ktot per processor (sum Ki for the blocks of each process)
    Allocate(Ktot_procs(nber_procs))
    Ktot_procs = 0;
    Do pp=1,nber_procs
       Do ii=1,MPI_CBFM_Blocks(pp,1)
            K = MPI_CBFM_Blocks(pp,1+ii)
            Ktot_procs(pp) = Ktot_procs(pp) + K_patchs_all(K)
       EndDo
    EndDo
    Ktot_proc_max = maxval(Ktot_procs);

    if (rank == 0) then
        if (debug_mode .eq. 1) then
            Write(10,'(a)') '--> K_blocks ='
            Do I=1,Nblocks-1
  	        Write(10,'(i4,a)', advance='no') K_patchs_all(I),';'
            Enddo
            Write(10,'(i4)') K_patchs_all(Nblocks)
            Write(10,*) ''
        endif
        Write(*,*)
        Write(*,'(a,i6)') 'Total number of CBFs = ',K_total
        Write(*,'(a,F8.1)') 'Compression Rate = ', CR
        Write(*,'(a,F10.4,a)') '--> We kept ', 100*(1./CR),' % of the initial MoM matrix'

        Write(*,*)
        call date_and_time(date_final_N1,time_final_N1,zone_final_N1,values_final_N1)
        call Calcul_time_spent(values_init_N1,values_final_N1, time_calcul_N1)
        Write(*,'(a,i2,a,i2,a,i2,a,i2,a)') '--> to generate the CBFs : ',time_calcul_N1(1),'j',time_calcul_N1(2)&
        ,'h',time_calcul_N1(3),'min', time_calcul_N1(4),'sec'
    endif

    !!! HERE for performance analysis purposes, I am going to store the information (Nblocks,Ncells,K) per MPI job
    !! these data will be used later jointly with the remora output data to better understand the computational performance of my MPI code
    if (debug_mode .eq. 1) then
        call Write_jobs_sim_info(CBFM_Blocks,MPI_CBFM_Blocks,K_patchs_all);
    endif

    !!********************************************************************************************************************************
    !! Now generation of the reduced matrix Zc ***************************************************************************************
    !! For the MPI code I will only consider the simplest case (in terms of programming) which is the heterogeneous simulation scene
    !! It is useless to program the homogenous case since a realistic particle / EM simulation scene is usually heterogeneous
    !! So each job will have few row bands of the total reduced matrix
    If (rank ==0) Then
      Write (*,*)
      Write (*,*)
      Write(*,'(a)',advance='no') 'Generation of the reduced matrix'
      Write(*,'(a)') ' (full storage)' ! for the moment Later we can consider to sparsify the reduced matrix
      time_calcul_N1 = 0
      call date_and_time(date_init_N1,time_init_N1,zone_init_N1,values_init_N1)
      Write(*,'(a)') 'InJob Interactions :'
    EndIf

    ! Prepare the cursors to indicate the positions (row and column) of Zc_ii,jj
    ! in the total band of Zc stored by each MPI job
    Allocate(curs_row_Z_blocks(MyNBlocks))
    Allocate(curs_col_Z_blocks(NBlocks))
    !curs_row_Z_blocks
    curs_row_Z_blocks(1) =1;
    Do kk_job=2, MyNBlocks
        prev_kk = MPI_CBFM_Blocks(rank+1,kk_job)
        curs_row_Z_blocks(kk_job) = curs_row_Z_blocks(kk_job-1) + K_patchs_all(prev_kk)
    Enddo
    !curs_col_Z_blocks
    curs_col_Z_blocks(1) =1;
    Do kk=2, NBlocks
        curs_col_Z_blocks(kk) = curs_col_Z_blocks(kk-1) + K_patchs_all(kk-1)
    Enddo

    ! Allocate Zreduite per MPI job
    !! this Zreduite is the (rank+1) th horizantal band of the total Zreduite; Here we allocate with Ktot_proc_max because Zreduite will receive later the Zreduite
    !! of the other processors when distributing the total Zc (see line 726)
    Allocate(Zreduite(Ktot_proc_max,K_total),Vreduit(Ktot_proc_max,2*NTr));
    call print_allocate(31,'Zreduite(Ktot_proc_max,K_total)','DCOMP',Ktot_proc_max*K_total);
    call print_allocate(34,'Vreduit(Ktot_proc_max,2*NTr)','DCOMP',Ktot_proc_max*2*NTr);

    Zreduite = 0D0;
    Vreduit= 0D0;

    !! Zloc
    ! Let us also allocate Zloc here to be sure that all this will fit in memory before following to the the genration of the reduced matrix
    ! create the processors grid
    Allocate(DESCA(9),DESCB(9))
    DESCA = 0;DESCB = 0;
    IRSRC = 0; ICSRC =0;

    CALL MPI_DIMS_CREATE(nber_procs,NDIMS,dims,INFO)
    NPROW = dims(1); NPCOL=dims(2);

    CALL BLACS_GET(-1,0,icontxt);
    CALL BLACS_GRIDINIT(icontxt,'Row-major', NPROW, NPCOL)
    CALL BLACS_GRIDINFO(icontxt,NPROW,NPCOL,MYROW,MYCOL)

    ! the following three lines are equivalent to : (coded to check/understand how numroc operates)
    !call numroc_2d_opt(K_total,K_total,2*NTr,NPROW,NPCOL,Myrow,Mycol,Mlocal,Nlocal,NRHSlocal);
    Mlocal = numroc(K_total,M_B,Myrow,ROW_SRC,NPROW)
    Nlocal = numroc(K_total,N_B,Mycol,COL_SRC,NPCOL)
    NRHSlocal = numroc(2*NTr,N_B,Mycol,COL_SRC,NPCOL)
    if (rank .eq. 0) then
        Write(*,'(a,i6,a,i6,a,i6)') 'Mloc ~= ', Mlocal, ', Nloc ~= ', Nlocal,'and NRHSlocal ~= ',NRHSlocal
    endif

    Allocate(ZredLoc(Mlocal,Nlocal),VredLoc(Mlocal,NRHSlocal))
    ZredLoc = 0.0e0; VredLoc = 0.0e0
    call print_allocate(22,'ZredLoc(Mlocal,Nlocal)','DCOMP',Mlocal*Nlocal);
    call print_allocate(26,'VredLoc(Mlocal,NRHSlocal)','DCOMP',Mlocal*NRHSlocal);
    !-------------------------------------------------------------------------------------------------------------------------------------


    if ((debug_mode .eq. 1) .and. (track_time_Zc .eq. 1)) then
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
            Write(rank_str,'(i4)') rank;
        endif
        file_name = trim(analysis_fold_name)//Env_sep//'TrackZctime_j'//rank_str//'.dat';
        num_file_time_Zc = 50+nber_procs+rank;
        Open(num_file_time_Zc,File = trim(file_name));
        time_calcul_Zc = 0
        call date_and_time(date_init_Zc,time_init_Zc,zone_init_Zc,values_init_Zc);
    endif

    ! HERE START THE INJOB INTERACTIONS
    alpha = 1.;beta = 0. ! for the matrix multiplication
    ! line up the MPI jobs
    call MPI_BARRIER(MPI_COMM_WORLD,code);
    ! Start with computing and storing interactions between blocks proper to current MPI job ----------------------------------------------
    !--------------------------------------------------------------------------------------------------------------------------------------
    Do ii_job=1, MyNBlocks
      ii = MPI_CBFM_Blocks(rank+1,1+ii_job);
      !Write(*,'(a,i3,a,i5)') 'Job ',rank,' - Interactions with Block ',ii
      size1 = CBFM_Blocks(ii)%Nbc_b
      NbreLig_mat1 = 3*size1
      NbreCol_mat1 = K_patchs_all(ii)

      Allocate(Cells_Block_ii(size1))
      Cells_Block_ii(1:size1) = Cells(curs_B_cel(ii):curs_B_cel(ii)+size1-1)
      Allocate(Matrice1(NbreLig_mat1,NbreCol_mat1));
      call print_allocate(35,'Matrice1(NbreLig_mat1,NbreCol_mat1)','DCOMP',NbreLig_mat1*NbreCol_mat1);
      Matrice1 = C_job_patchs(curs_B_Cpatch(ii_job):curs_B_Cpatch(ii_job)+NbreLig_mat1-1,1:NbreCol_mat1)

      DO jj_job=1, MyNBlocks
          jj = MPI_CBFM_Blocks(rank+1,1+jj_job);
          size3 = CBFM_Blocks(jj)%Nbc_b
          NbreLig_mat3 = 3*size3
          NbreCol_mat3 = K_patchs_all(jj)

          Allocate(Cells_Block_jj(size3))
          Cells_Block_jj(1:size3) = Cells(curs_B_cel(jj):curs_B_cel(jj)+size3-1)

          Allocate(Matrice2(NbreLig_mat1,NbreLig_mat3));
          call print_allocate(35,'Matrice2(NbreLig_mat1,NbreLig_mat3)','DCOMP',NbreLig_mat1*NbreLig_mat3);
          Call Green_s_tr_partial(size1,Cells_Block_ii,size3,Cells_Block_jj,Matrice2)

          Allocate(Matrice3(NbreLig_mat3,NbreCol_mat3));
          call print_allocate(35,'Matrice3(NbreLig_mat3,NbreCol_mat3)','DCOMP',NbreLig_mat3*NbreCol_mat3);
          Matrice3 = C_job_patchs(curs_B_Cpatch(jj_job):curs_B_Cpatch(jj_job)+NbreLig_mat3-1,1:NbreCol_mat3)

          !! Zc(ii,jj)**********************************
          Allocate(Mat_Inter(NbreCol_mat1,NbreLig_mat3));
          Allocate(MatProduit(NbreCol_mat1,NbreCol_mat3));
          call print_allocate(37,'Mat_Inter(NbreCol_mat1,NbreLig_mat3)','DCOMP',NbreCol_mat1*NbreLig_mat3);
          call print_allocate(38,'MatProduit(NbreCol_mat1,NbreCol_mat3)','DCOMP',NbreCol_mat1*NbreCol_mat3);

          CALL ZGEMM('T','N',NbreCol_mat1,NbreLig_mat3,NbreLig_mat1,alpha,Matrice1,NbreLig_mat1,Matrice2,&
           NbreLig_mat1,beta,Mat_Inter,NbreCol_mat1)
          CALL ZGEMM('N','N',NbreCol_mat1,NbreCol_mat3,NbreLig_mat3,alpha,Mat_Inter,NbreCol_mat1,Matrice3,&
            NbreLig_mat3,beta,MatProduit,NbreCol_mat1)
          Zreduite(curs_row_Z_blocks(ii_job):curs_row_Z_blocks(ii_job)+NbreCol_mat1-1,curs_col_Z_blocks(jj):&
          curs_col_Z_blocks(jj)+NbreCol_mat3-1) = MatProduit

          Deallocate(Cells_Block_jj,Matrice2,Matrice3,Mat_inter,MatProduit)
      ENDDO
      Deallocate(Cells_Block_ii,Matrice1)
    EndDo

    if (track_time_Zc .eq. 1) then
        call date_and_time(date_final_Zc,time_final_Zc,zone_final_Zc,values_final_Zc)
        call Calcul_time_spent(values_init_Zc,values_final_Zc, time_calcul_Zc)
        Write (num_file_time_Zc,'(a,i6,a)') '* Injob interactions    : ',86400*time_calcul_Zc(1)+3600*time_calcul_Zc(2)+60*time_calcul_Zc(3)+time_calcul_Zc(4),' sec'
    endif


    ! Prepare a C_trans_patchs that will be used to transfer in a cyclic way the CBFs
    Call MPI_ALLREDUCE(Nbc_proc,Nbc_proc_max,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,code)
    K_max = maxval(K_patchs_all);
    Allocate(C_trans_patchs(3*Nbc_proc_max,K_max));
    call print_allocate(37,'C_trans_patchs(3*Nbc_proc_max,K_max)','DCOMP',3*Nbc_proc_max*K_max);
    C_trans_patchs = 0;
    C_trans_patchs(1:3*Nbc_proc,1:Kmax_job) = C_job_patchs(1:3*Nbc_proc,1:Kmax_job);
    sizeTrans = 3*Nbc_proc_max*K_max

    ! Now Interactions between blocks of different MPI jobs ------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------------------------------------------------------------
    ! identify previous and following MPI job
    prev_job = mod(nber_procs+rank-1,nber_procs)
    next_job = mod(rank+1,nber_procs)
    !call MPI_BARRIER(MPI_COMM_WORLD,code);
    if (rank == 0) then
      Write(*,'(a)') 'InterJob Interactions :'
    EndIf
    Do kk = 1, nber_procs-1
      if (rank == 0) Then
        Write(*,'(a,i4)') 'CBFs Transition kk = ',kk
      endif

      if (track_time_Zc .eq. 1) then
          time_calcul_Zc = 0
          call date_and_time(date_init_Zc,time_init_Zc,zone_init_Zc,values_init_Zc);
      endif

      ! SEND/RECEIVE the CBFs from/to the prev/next MPI job
      !! Every job/process sends its 'C_tot_patch_Zc_trans' to the previous process and receive the one of the next process
      Call MPI_SENDRECV_REPLACE(C_trans_patchs,2*sizeTrans,MPI_COMPLEX,prev_job,tag,&
                              	next_job,tag,MPI_COMM_WORLD,status,code)

      if (track_time_Zc .eq. 1) then
          call date_and_time(date_final_Zc,time_final_Zc,zone_final_Zc,values_final_Zc)
          call Calcul_time_spent(values_init_Zc,values_final_Zc, time_calcul_Zc)
          Write (num_file_time_Zc,'(a,i5,a,i6,a)') '-T',kk,' SendRecv        : ',86400*time_calcul_Zc(1)+3600*time_calcul_Zc(2)+60*time_calcul_Zc(3)+time_calcul_Zc(4),'sec'
      endif

      !! ATTENTION : We should consider the change of the content of 'C_tot_patch_Zc_trans' through the different kk
      !! Example : at kk = 1; C_tot_patch_Zc_trans of rank 0 includes the CBFs calculated by 0
      !!           at kk = 2; C_tot_patch_Zc_trans of rank 0 includes the CBFs calculated by Nprocs
      !!           at kk = n; C_tot_patch_Zc_trans of rank 0 includes the CBFs calculated by Nprocs- (n-1) + 1
      !! So we should identify the effective job which created the CBFs included in the received C_trans_patchs for the current kk
      rank_EffNextjob = mod(rank+kk,nber_procs)
      NBlocks_EffNextjob = MPI_CBFM_Blocks(rank_EffNextjob+1,1) ! To interact with the 'MyNBlocks' blocks of the current job

      ! For each k, We need these vectors to localise the cells and CBFs for the next effective job
      Allocate(curs_B_Cpatch_Next(NBlocks_EffNextjob))
      curs_B_Cpatch_Next(1)=1
      DO jj_Nextjob=2, NBlocks_EffNextjob
        prev_jj = MPI_CBFM_Blocks(rank_EffNextjob+1,jj_Nextjob);  ! equivalent to kk-1 1+kk_job-1
        curs_B_Cpatch_Next(jj_Nextjob) = curs_B_Cpatch_Next(jj_Nextjob-1) + 3*CBFM_Blocks(prev_jj)%Nbc_b
      EndDo

      if (track_time_Zc .eq. 1) then
          time_calcul_Zc = 0
          call date_and_time(date_init_Zc,time_init_Zc,zone_init_Zc,values_init_Zc);
      endif

      ! So now Interactions between the 'MyNBlocks' blocks of the current job and the 'NBlocks_EffNextjob' blocks of the effective next job
      DO ii_job=1,MyNBlocks
        ii = MPI_CBFM_Blocks(rank+1,1+ii_job);
        !Write(*,'(a,i3,a,i5)') 'Job ',rank,' - Interactions with Block ',ii
        size1 = CBFM_Blocks(ii)%Nbc_b
        NbreLig_mat1 = 3*size1
        NbreCol_mat1 = K_patchs_all(ii)

        Allocate(Cells_Block_ii(size1))
        Cells_Block_ii(1:size1) = Cells(curs_B_cel(ii):curs_B_cel(ii)+size1-1)

        Allocate(Matrice1(NbreLig_mat1,NbreCol_mat1));
        call print_allocate(35,'Matrice1(NbreLig_mat1,NbreCol_mat1)','DCOMP',NbreLig_mat1*NbreCol_mat1);
        Matrice1 = C_job_patchs(curs_B_Cpatch(ii_job):curs_B_Cpatch(ii_job)+NbreLig_mat1-1,1:NbreCol_mat1);

        !! ATTENTION TU PEUX ICI EVITER DE REMPLIR DEUX FOIS MATRICE 1 EN AJOUTANT UNE CONDITION SUR kk :
        !! SI kk=1 EFFECTUER LES INTERACTIONS EN INTERNE IL SE PEUT QUE CA NE SERVE A RIEN CA AJOUTE UNE CONDITION PAR BLOC LOCAL DONC DU TEMPS ...

        DO jj_Nextjob=1, NBlocks_EffNextjob
            jj = MPI_CBFM_Blocks(rank_EffNextjob+1,1+jj_Nextjob);
            size3 = CBFM_Blocks(jj)%Nbc_b
            NbreLig_mat3 = 3*size3
            NbreCol_mat3 = K_patchs_all(jj)

            Allocate(Cells_Block_jj(size3))
            Cells_Block_jj(1:size3) = Cells(curs_B_cel(jj):curs_B_cel(jj)+size3-1)

            Allocate(Matrice2(NbreLig_mat1,NbreLig_mat3));
            call print_allocate(35,'Matrice2(NbreLig_mat1,NbreLig_mat3)','DCOMP',NbreLig_mat1*NbreLig_mat3);
            Call Green_s_tr_partial(size1,Cells_Block_ii,size3,Cells_Block_jj,Matrice2)

            Allocate(Matrice3(NbreLig_mat3,NbreCol_mat3));
            call print_allocate(35,'Matrice3(NbreLig_mat3,NbreCol_mat3)','DCOMP',NbreLig_mat3*NbreCol_mat3);
            Matrice3 = C_trans_patchs(curs_B_Cpatch_Next(jj_Nextjob):curs_B_Cpatch_Next(jj_Nextjob)+NbreLig_mat3-1,1:NbreCol_mat3)

            !! Zc(ii,jj)**********************************
            Allocate(Mat_Inter(NbreCol_mat1,NbreLig_mat3));
            Allocate(MatProduit(NbreCol_mat1,NbreCol_mat3));
            call print_allocate(37,'Mat_Inter(NbreCol_mat1,NbreLig_mat3)','DCOMP',NbreCol_mat1*NbreLig_mat3);
            call print_allocate(38,'MatProduit(NbreCol_mat1,NbreCol_mat3)','DCOMP',NbreCol_mat1*NbreCol_mat3);

            CALL ZGEMM('T','N',NbreCol_mat1,NbreLig_mat3,NbreLig_mat1,alpha,Matrice1,NbreLig_mat1,Matrice2,&
             NbreLig_mat1,beta,Mat_Inter,NbreCol_mat1)
            CALL ZGEMM('N','N',NbreCol_mat1,NbreCol_mat3,NbreLig_mat3,alpha,Mat_Inter,NbreCol_mat1,Matrice3,&
              NbreLig_mat3,beta,MatProduit,NbreCol_mat1)
            Zreduite(curs_row_Z_blocks(ii_job):curs_row_Z_blocks(ii_job)+NbreCol_mat1-1,curs_col_Z_blocks(jj):&
            curs_col_Z_blocks(jj)+NbreCol_mat3-1) = MatProduit

            Deallocate(Cells_Block_jj,Matrice2,Matrice3,Mat_inter,MatProduit)
        EndDo
        Deallocate(Cells_Block_ii,Matrice1)
      EndDo
      Deallocate(curs_B_Cpatch_Next)

      if (track_time_Zc .eq. 1) then
          call date_and_time(date_final_Zc,time_final_Zc,zone_final_Zc,values_final_Zc)
          call Calcul_time_spent(values_init_Zc,values_final_Zc, time_calcul_Zc)
          Write (num_file_time_Zc,'(a,i5,a,i6,a)') '+T',kk,' Zijs Calculation: ',86400*time_calcul_Zc(1)+3600*time_calcul_Zc(2)+60*time_calcul_Zc(3)+time_calcul_Zc(4),'sec'
      endif
    EndDo
    Deallocate(C_trans_patchs);
    call MPI_BARRIER(MPI_COMM_WORLD,code);
    if (track_time_Zc .eq. 1) then
        close(num_file_time_Zc);
    endif
    ! End Calculation of Reduced Matrix Zc ******************************************************************************************

    ! GENERATION OF THE REDUCED FIELD Vc*********************************************************************************************
    if (rank == 0) then
      call date_and_time(date_final_N1,time_final_N1,zone_final_N1,values_final_N1)
      call Calcul_time_spent(values_init_N1,values_final_N1, time_calcul_N1)
      Write (*,'(a,i2,a,i2,a,i2,a,i2,a)') '--> to generate Zc ',time_calcul_N1(1),'j',time_calcul_N1(2)&
      ,'h',time_calcul_N1(3),'min', time_calcul_N1(4),'sec'
      Write (*,*)
      Write(*,'(a)') 'Generation of the reduced Incident Field'
      time_calcul_N1 = 0
      call date_and_time(date_init_N1,time_init_N1,zone_init_N1,values_init_N1)
    EndIf

    DO ii_job=1, MyNBlocks
        ii = MPI_CBFM_Blocks(rank+1,1+ii_job);
        size= CBFM_Blocks(ii)%Nbc_b
        NbreLig_mat1 = 3*size;
        NbreCol_mat1 = K_patchs_all(ii)

        Allocate(Cells_Block_ii(size))
        Cells_Block_ii(1:size) = Cells(curs_B_cel(ii):curs_B_cel(ii)+size-1)

        Allocate(E_ref_incident(NbreLig_mat1,2*NTr));
        call print_allocate(35,'E_ref_incident(NbreLig_mat1,2*NTr)','DCOMP',NbreLig_mat1*2*NTr);
        Call Incident_Field(1,size,Cells_Block_ii,NTr,Transmitters,1,NTr,E_ref_incident);

        Allocate(Matrice1(NbreLig_mat1,NbreCol_mat1));
        call print_allocate(36,'Matrice1(NbreLig_mat1,NbreCol_mat1)','DCOMP',NbreLig_mat1*NbreCol_mat1);
        Matrice1 = C_job_patchs(curs_B_Cpatch(ii_job):curs_B_Cpatch(ii_job)+NbreLig_mat1-1,1:NbreCol_mat1)

        DO jj=1, 2*NTr
          Allocate(Vect(NbreLig_mat1))
          Vect = E_ref_incident(:,jj)
          Allocate(VectProduit(NbreCol_mat1))
          CALL ZGEMV('T',NbreLig_mat1,NbreCol_mat1,alpha,Matrice1,NbreLig_mat1,Vect,1,beta,VectProduit,1)
          Vreduit(curs_row_Z_blocks(ii_job):curs_row_Z_blocks(ii_job)+NbreCol_mat1-1,jj) = VectProduit
          Deallocate(Vect,VectProduit)
        ENDDO
        Deallocate(Cells_Block_ii,Matrice1,E_ref_incident)
    ENDDO

    If (rank==0) then
      call date_and_time(date_final_N1,time_final_N1,zone_final_N1,values_final_N1)
      call Calcul_time_spent(values_init_N1,values_final_N1, time_calcul_N1)
      Write (*,'(a,i2,a,i2,a,i2,a,i2,a)') '--> to compute Vc ',time_calcul_N1(1),'j',time_calcul_N1(2)&
      ,'h',time_calcul_N1(3),'min', time_calcul_N1(4),'sec'
    EndIf
    call MPI_BARRIER(MPI_COMM_WORLD,code);
    ! End Reduced Incident Field Vc*********************************************************************************************

    !!! Resolution du systeme matriciel reduit**********************************************************************************
    if (rank == 0) then
        Write(*,*) ''
        Write(*,'(a,i8)') 'Resolution of the final problem of size K_total =',K_total
        if (debug_mode .eq. 1) then
            Write(10,*) ''
            Write(10,'(a,i8)') 'Size of the reduced matrix Zc = ', K_total
        endif

        time_calcul_N1 = 0
        call date_and_time(date_init_N1,time_init_N1,zone_init_N1,values_init_N1)
    endif

    !** Distribute Compressed matrix and incident field vector *****************************
    if (rank == 0) then
        Write(*,'(a)') '+ Block Cyclic Data Distribution'
    endif

    Allocate(Curs_Kpatchs_all(NBlocks+1)); ! equivalent to curs_B_cel for Ni
    Curs_Kpatchs_all(1) = 1;
    Do kk=2,NBlocks+1
        Curs_Kpatchs_all(kk) = Curs_Kpatchs_all(kk-1)+K_patchs_all(kk-1);
    EndDo

    Call Block_Cyclic_Distribution_v1(Ktot_proc_max,K_total,2*NTr,Zreduite,Vreduit,K_patchs,K_patchs_all,Curs_Kpatchs_all, Ktot_procs,MPI_CBFM_Blocks,ZredLoc,VredLoc);
    Deallocate(Zreduite,Vreduit);

    if (rank == 0) then
        Write(*,'(a)') '+ Resolution with SCALAPACK'
    endif

    ! RESOLUTION WITH SCALAPACK ------------------------------------------------
    CALL DESCINIT(DESCA,K_total,K_total,M_B,N_B,IRSRC,ICSRC,icontxt,Mlocal,INFO)
    CALL DESCINIT(DESCB,K_total,2*NTr,M_B,N_B,IRSRC,ICSRC,icontxt,Mlocal,INFO)

    IA=1;JA=1;IB=1;JB=1;
    ! Machine precision
    EPSMCH = PSLAMCH(icontxt,'E');
    ! get Infinity NORM of ZredLoc

    ! LWORK >= 2*LOCr(N+MOD(IA-1,MB_A)) +  MAX( 2, MAX(NB_A*CEIL(NPROW-1,NPCOL),LOCc(N+MOD(JA-1,NB_A)) + NB_A*CEIL(NPCOL-1,NPROW)) ).
    ! LRWORK >= MAX( 1, 2*LOCc(N+MOD(JA-1,NB_A)) ).
    LWORK = 2*K_total; ! Here 2*K_total is enough for now, but needs to compute the exact LWORK and LIWORK for more accuracy/robustness
    LIWORK = 2*K_total; !K_total;
    Allocate(WORK(LWORK),IWORK(LIWORK));
    ANORM = PZLANGE( 'I', K_total,K_total, ZredLoc, IA,JA,DESCA,WORK);
    Allocate(IPIV(Mlocal+M_B))
    ! Solve system with PSGESV ; the solution X overwrites the RHS B
    CALL PZGESV(K_total,2*NTr,ZredLoc,IA,JA,DESCA,IPIV,VredLoc,IB,JB,DESCB,INFO);
    if (INFO .GT. 0) then
        Write(*,'(a)') 'SINGULAR MATRIX';
    elseif (K_total .gt. 0) then
        ! Get Reciprocal condition number RCOND of Zc
        !LWORK = K_total; LIWORK = K_total;
        !CALL PZGECON( 'I', K_total,ZredLoc,IA,JA,DESCA, ANORM, RCOND,WORK,LWORK,IWORK,LIWORK,INFO);
        CALL PZGECON( 'I', K_total,ZredLoc,IA,JA,DESCA, ANORM, RCOND,WORK,LWORK,IWORK,LIWORK,INFO);
        !RCOND = max(RCOND,EPSMCH);
        ERRBD = EPSMCH/ RCOND
        if (rank == 0) then
          Write(*,'(a,es10.3)') '--> ANORM(ZredLoc) = ',ANORM
          Write(*,'(a,es12.5)') '--> RCOND = ',RCOND
          Write(*,'(a,es12.5,a,es12.5)') '--> With EPSMCH =',EPSMCH,'; ERRBD =',ERRBD
        endif
    endif
    ! solve system with PSGESVX (Expert Driver)
    !Allocate(R(K_total),C(K_total),FERR(N_B),BERR(N_B));
    !Allocate(ZFredLoc(Mlocal,Nlocal),AlphaLoc(Mlocal,NRHSlocal))
    !LWORK = K_total;LIWORK = K_total;
    !CALL PZGESVX('E','N',K_total,2*NTr,ZredLoc,IA,JA,DESCA,ZFredLoc,IA,JA,DESCA,IPIV,EQUED,R,C,VredLoc,IB,JB,DESCB,AlphaLoc,IB,JB,DESCB,RCOND,FERR,BERR,WORK,LWORK,IWORK,LIWORK,INFO);
    !if (INFO .GT. 0) then
    !    Write(*,'(a)') 'SINGULAR MATRIX';
    !elseif (K_total .gt. 0) then
    !    RCOND = max(RCOND,EPSMCH);
    !    if (rank == 0) then
    !      Write(*,'(a,es12.5)') '--> RCOND = ',RCOND
    !      Write(*,'(a,es12.5)') '--> FERR =',FERR(1)
    !    endif
    !endif
    CALL BLACS_GRIDEXIT(icontxt);
    Deallocate(ZredLoc);!,VredLoc);
    !********************************************************************************


    ! Now reorganize back the solution alpha
    Allocate(AlphaProc(Ktot_proc_max,2*NTr));
    call print_allocate(31,'AlphaProc(Ktot_proc_max,2*NTr)','DCOMP',Ktot_proc_max*2*NTr);
    AlphaProc = 0D0
    prev_job = mod(nber_procs+rank-1,nber_procs)
    next_job = mod(rank+1,nber_procs)
    size2 = Ktot_proc_max*2*NTr

    Do pp=1, nber_procs
        Call MPI_SENDRECV_REPLACE(AlphaProc,2*size2,MPI_COMPLEX,prev_job,tag,next_job,tag,MPI_COMM_WORLD,status,code)
        Eff_rank = mod(rank+pp,nber_procs);
	    K_proc = Ktot_procs(Eff_rank+1);
        NBlocks_eff = MPI_CBFM_Blocks(Eff_rank+1,1)
        Allocate(Blocks_eff(NBlocks_eff), K_patchs_eff(NBlocks_eff));
        Blocks_eff(1:NBlocks_eff) = MPI_CBFM_Blocks(Eff_rank+1,2:1+NBlocks_eff)
        K_patchs_eff(1:NBlocks_eff) = K_patchs_all(Blocks_eff);
        Do iLoc=1, Mlocal
    	    iGlob = INDXL2G(iLoc,M_B,Myrow,0,NPROW)
            block_iGlob = minloc(Curs_Kpatchs_all, 1, mask = Curs_Kpatchs_all .gt. iGlob) - 1
            pos_loc_block = minloc(MPI_CBFM_Blocks(Eff_rank+1,2:1+NBlocks_eff), 1, mask = MPI_CBFM_Blocks(Eff_rank+1,2:1+NBlocks_eff) .eq. block_iGlob);
            if ((pos_loc_block .gt. 1) .OR. ((pos_loc_block .eq. 1) .AND. (MPI_CBFM_Blocks(Eff_rank+1,2) .eq. block_iGlob))) then
                iLoc_proc = sum(K_patchs_eff(1:pos_loc_block-1)) + (iGlob - sum(K_patchs_all(1:block_iGlob-1)));
                Do jLoc=1, NRHSlocal
            		jGlob = INDXL2G(jLoc,N_B,Mycol,0,NPCOL)
            		AlphaProc(iLoc_proc,jGlob) = VredLoc(iLoc,jLoc)
                End do
            endif
        enddo
        deallocate(Blocks_eff,K_patchs_eff);
    End do
    deallocate(Curs_Kpatchs_all,K_patchs_all);!,AlphaLoc);

    if (rank == 0) then
        call date_and_time(date_final_N1,time_final_N1,zone_final_N1,values_final_N1)
        call Calcul_time_spent (values_init_N1,values_final_N1, time_calcul_N1)
        Write (*,'(a,i2,a,i2,a,i2,a,i2,a)') '--> to resolve Zc * Alpha = Vc ',time_calcul_N1(1),'j',time_calcul_N1(2)&
        ,'h',time_calcul_N1(3),'min', time_calcul_N1(4),'sec'

        time_calcul_N1 = 0
        call date_and_time(date_init_N1,time_init_N1,zone_init_N1,values_init_N1)
    endif


    Allocate(E_total(3*Nbc_proc,2*NTr));
    call print_allocate(26,'E_total(3*Nbc_proc,2*NTr)','DCOMP',3*Nbc_proc*2*NTr);
    E_total = 0D0
    DO kk=1,2*NTr
        curseurf = 1
        curs_lig_E = 1
        DO ii_job=1,MyNBlocks
            ii = MPI_CBFM_Blocks(rank+1,1+ii_job);
            BlockSize = 3*CBFM_Blocks(ii)%Nbc_b
            DO jj=1,K_patchs(ii_job)
                Allocate(VectProduit(BlockSize))
                VectProduit = AlphaProc(curseurf,kk)*C_job_patchs(curs_B_Cpatch(ii_job):curs_B_Cpatch(ii_job)+BlockSize-1,jj)
                E_total(curs_lig_E:curs_lig_E+BlockSize-1,kk)= E_total(curs_lig_E:curs_lig_E+BlockSize-1,kk)+ VectProduit
                curseurf = curseurf + 1
                Deallocate(VectProduit)
            ENDDO
            curs_lig_E = curs_lig_E + BlockSize
        ENDDO
    ENDDO

    if (rank == 0) then
        call date_and_time(date_final_N1,time_final_N1,zone_final_N1,values_final_N1)
        call Calcul_time_spent (values_init_N1,values_final_N1, time_calcul_N1)
        Write (*,'(a,i2,a,i2,a,i2,a,i2,a)') '--> to compute E as linear combination of Alphas and CBFs ',time_calcul_N1(1)&
            ,'j',time_calcul_N1(2),'h', time_calcul_N1(3),'min', time_calcul_N1(4),'sec'
    endif

    Deallocate(AlphaProc);
    if (Nfreq .eq. 1) then
        Deallocate(C_job_patchs);
    endif

    If ((save_Eint .eq. 1) .and. (Nbc .le. save_Eint_Nmax)) then

        file_name = trim(SimOutfld_name)//Env_sep//'Ein_CBFME_MPI.dat';
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

        !! TOUT CA NE MARCHE PAS POUR LE MOMENT
        !file_name = trim(SimOutfld_name)//Env_sep//'Ein_MPI.bin';
        !if (rank == 0) then
        !    open(unit=14, file = trim(file_name), form = 'UNFORMATTED')
        !    Do ii=1,3*Nbc_proc
        !        Do kk=1,2*NTr
        !            Write(14) real(E_total(ii,kk)),imag(E_total(ii,kk))
        !        EndDo
        !    EndDo
        !    close(14)
        !endif
        !call MPI_BARRIER(MPI_COMM_WORLD,code);
        !do pp = 1, nber_procs - 1
        !    if(rank == pp) then
        !        open(unit = 14, file = trim(file_name), form = 'unformatted', status = 'old', position = 'append')
        !        Do ii=1,3*Nbc_proc
        !            Do kk=1,2*NTr
        !                Write(14) real(E_total(ii,kk)),imag(E_total(ii,kk))
        !            EndDo
        !        EndDo
        !        close(14)
        !    endif
        !    call MPI_BARRIER(MPI_COMM_WORLD,code);
        !enddo

        !nb_values = 10;
        !Allocate(values(nb_values));
        !values(:)= (/(i+rank*100,i=1,nb_values)/)
        !print *, "Write process",rank, ":",values(:)
        !
        !call MPI_FILE_OPEN( MPI_COMM_WORLD,"data.dat", MPI_MODE_WRONLY + MPI_MODE_CREATE, &
        !                    MPI_INFO_NULL,fh,code)
        !IF (code /= MPI_SUCCESS) THEN
        !    PRINT *, 'Error in opening file'
        !    CALL MPI_ABORT( MPI_COMM_WORLD, 42, code)
        !END IF
        !call MPI_TYPE_SIZE( MPI_INTEGER,bytes_in_integer,code)
        !
        !offset_file=rank*nb_values*bytes_in_integer
        !call MPI_FILE_WRITE_AT(fh,offset_file,values,nb_values, MPI_INTEGER, &
        !                        status,code)
        !IF (code /= MPI_SUCCESS) THEN
        !    PRINT *, 'Error in writing file'
        !    CALL MPI_ABORT( MPI_COMM_WORLD, 6, code)
        !END IF
        !
        !call MPI_FILE_CLOSE(fh,code);

        ! Open the file
        !file_name = trim(SimOutfld_name)//Env_sep//'Ein_MPI.bin';
        !call MPI_FILE_OPEN(MPI_COMM_WORLD,file_name,MPI_MODE_RDWR + MPI_MODE_CREATE, MPI_INFO_NULL,fh,code);
        !IF (code /= MPI_SUCCESS) THEN
        !    PRINT *, 'Error in opening file'
        !    CALL MPI_ABORT( MPI_COMM_WORLD, 42, code)
        !END IF
        !
        !! Prepare offset_file and nb_values then Write
        !call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION,bytes_in_dp,code)
        !Allocate(Nbc_proc_tot(nber_procs));
        !Do pp=1,nber_procs
        !    MyNBlocks = MPI_CBFM_Blocks(pp,1);
        !    Nbc_proc_tot(pp) =  sum(CBFM_Blocks(MPI_CBFM_Blocks(pp,2:1+MyNBlocks))%Nbc_b);
        !EndDo
        !nb_values = 1;
        !shift1 = 3*sum(Nbc_proc_tot(1:rank));
        !Write(*,*) 'I am rank ',rank,'; My 3Nbc = ',3*Nbc_proc,'; My shift1 = ', shift1
        !Do kk = 1,2*NTr
        !    shift0 = (kk-1)*3*Nbc
        !    Do ii=1,3*Nbc_proc
        !        shift2 = (ii-1)
        !        offset_file = 2*(shift0+shift1+shift2)*bytes_in_dp;
        !        e_re = dreal(E_total(ii,kk));
        !        e_im = dimag(E_total(ii,kk));
        !        call MPI_FILE_WRITE_AT(fh,offset_file,e_re,nb_values,MPI_DOUBLE_PRECISION,status,code);
        !        call MPI_FILE_WRITE_AT(fh,offset_file+bytes_in_dp,e_im,nb_values,MPI_DOUBLE_PRECISION,status,code);
        !    EndDo
        !EndDo
        !IF (code /= MPI_SUCCESS) THEN
        !    PRINT *, 'Error in writing file'
        !    CALL MPI_ABORT(MPI_COMM_WORLD, 6, code)
        !END IF
        !
        !call MPI_FILE_CLOSE(fh,code);

        if (rank == 0) then
            call date_and_time(date_final_N1,time_final_N1,zone_final_N1,values_final_N1)
            call Calcul_time_spent (values_init_N1,values_final_N1, time_calcul_N1)
            Write (*,'(i2,a,i2,a,i2,a,i2,a)') time_calcul_N1(2),'h', time_calcul_N1(3),'min', time_calcul_N1(4),'sec';
        endif
    endif

END SUBROUTINE Compute_EFields_CBFME
