SUBROUTINE Write_jobs_sim_info(CBFM_Blocks,MPI_CBFM_Blocks,K_patchs_all)

    USE Initialization
    USE common_variables
    USE MPI
    
    Implicit NONE
    
    !IN/OUT 
    type (CBFM_Block), Dimension(Nblocks), INTENT(IN) :: CBFM_Blocks
    Integer, Dimension(nber_procs,Nblk_proc_max+1), INTENT(IN) :: MPI_CBFM_Blocks
    Integer, Dimension(NBlocks), INTENT(IN) :: K_patchs_all
    
    ! local 
    Integer ii, Nblocks_proc, Ncells_proc, K_proc, a 
    CHARACTER(300) file_name,analysis_fold_name
    CHARACTER(:) ,allocatable :: stFreq
    CHARACTER(6) :: ty
    logical :: dirExists
    
    analysis_fold_name = trim(SimOutfld_name)//Env_sep//'Analysis';
    
    Call MPI_Barrier(MPI_COMM_WORLD,code);
    if (rank == 0) then         
        a = nint(Freq_w/1E9);
        if (a < 10) Then 
            Allocate(character(5) ::stFreq); ty = '(f5.3)';
        ElseIf (a < 100) Then
            Allocate(character(6) ::stFreq); ty = '(f6.3)';
        Else
            Allocate(character(7) ::stFreq); ty = '(f7.3)';
        EndIf     
        Write(stFreq,ty) Freq_w/1E9
    
        file_name = trim(analysis_fold_name)//Env_sep//'MPIjobs_loadinfo_'//stFreq//'GHz.dat'; 
    
        Open(14,File = trim(file_name))
        Write(14,'(a)') '   rank       Nblocks      Nbc          K';
        Do ii=1, nber_procs
            Nblocks_proc = MPI_CBFM_Blocks(ii,1);
            Ncells_proc = sum(CBFM_Blocks(MPI_CBFM_Blocks(ii,2:1+Nblocks_proc))%Nbc_b);
            K_proc = sum(K_patchs_all(MPI_CBFM_Blocks(ii,2:1+Nblocks_proc)));
            Write(14,'(i6,i12,i12,i12)') (ii-1), Nblocks_proc, Ncells_proc,K_proc
        EndDo
        Close(14)
    endif   
    Call MPI_Barrier(MPI_COMM_WORLD,code);
    
    !!!! Let's recap
    !! quelles sont les allocations majeures qui font la diff entre les MPI tasks/jobs
    !! on remarque que les deux valeurs qui menent essentiellement la dance sont Nbc_proc et Ktot_proc_max
    !Allocate(C_job_patchs(3*Nbc_proc,2*NTr_CBFM));
    !Allocate(C_job_patchs_tmp(3*Nbc_proc,Kmax_job));
    !Allocate(C_job_patchs(3*Nbc_proc,Kmax_job));    
    !
    !Allocate(Zreduite(Ktot_proc_max,K_total)); 
    !
    !Allocate(C_trans_patchs(3*Nbc_proc_max,K_max));
    !Allocate(curs_B_Cpatch_Next(NBlocks_EffNextjob));
    !
    !Allocate(Vreduit(Ktot_proc_max,2*NTr));    
    !
    !Allocate(ZredLoc(Mlocal,Nlocal),VredLoc(Mlocal,NRHSlocal));    
    !Allocate(AlphaProc(Ktot_proc_max,2*NTr));
    !
    !Allocate(E_total(3*Nbc_proc,2*NTr));
    !
    !
    !
    !!! mnt les allocations communes et/ou qui n'ont pas vraiment de poids majeurs sont :
    !Allocate(curs_B_cel(NBlocks),curs_B_Cpatch(MyNBlocks),K_patchs(MyNBlocks),K_patchs_all(NBlocks))
    !Allocate(blocks_to_sort(NBlocks),OrderBlocksProcs(Nblocks));
    !Ntests = 4; Allocate(testedBlks(Ntests))
    !Allocate(ExtSizes(NBlocks));
    !Allocate(fSR_blocks(MyNBlocks),nnz_blocks(MyNBlocks),spr_perc_blocks(MyNBlocks));
    !Allocate(Cells_Block(size));
    !
    !Allocate(EREFpatch_e(3*size,2*NTr_CBFM));
    !Allocate(Epatch_e(3*size,2*NTr_CBFM));
    !Allocate(Zpatch_e(2*klu+1,3*size));
    !Allocate(Zpatch_e_spr(nnz),row_sprZ(3*size+1),col_sprZ(nnz));
    !Allocate(Zpatch_e(3*size,3*size));
    !Allocate( AA(M,N), S(MIN(M,N)),U(M,M),VT(N,N), WW(MIN(M,N)-1));
    !Allocate(Cpatch_e(3*size, K));
    !
    !ALLOCATE(nb_elements_rec(nber_procs), deplts(nber_procs));
    !Allocate(vect_tmp(NBlocks));
    !Allocate(Ktot_procs(nber_procs));
    !Allocate(curs_row_Z_blocks(MyNBlocks));
    !Allocate(curs_col_Z_blocks(NBlocks)) ;
    !
    !Allocate(Cells_Block_ii(size1));
    !Allocate(Matrice1(NbreLig_mat1,NbreCol_mat1));
    !Allocate(Cells_Block_jj(size3));
    !Allocate(Matrice2(NbreLig_mat1,NbreLig_mat3));
    !Allocate(Matrice3(NbreLig_mat3,NbreCol_mat3));
    !Allocate(Mat_Inter(NbreCol_mat1,NbreLig_mat3));
    !Allocate(MatProduit(NbreCol_mat1,NbreCol_mat3));
    !
    !Allocate(Cells_Block_ii(size));
    !Allocate(E_ref_incident(NbreLig_mat1,2*NTr));
    !Allocate(Matrice1(NbreLig_mat1,NbreCol_mat1));
    !Allocate(Vect(NbreLig_mat1));
    !Allocate(VectProduit(NbreCol_mat1));
    !
    !Allocate(DESCA(9),DESCB(9));
    !Allocate(Curs_Kpatchs_all(NBlocks+1)); ! equivalent to curs_B_cel for Ni;
    !Allocate(K_patchs_eff(NBlocks_eff));
    !Allocate(IPIV(Mlocal+M_B));
    !
    !Allocate(K_patchs_eff(NBlocks_eff));
    !Allocate(VectProduit(BlockSize)); 
    
END SUBROUTINE Write_jobs_sim_info