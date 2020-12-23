SUBROUTINE Block_Cyclic_Distribution_v1(M_Z,N_Z,NRHS_Z,Zred,Vred,K_patchs,K_patchs_all,Curs_Kpatchs_all,Ktot_procs,MPI_CBFM_Blocks,ZredLoc,VredLoc)  
    
    USE Initialization
    USE common_variables
    USE lapack95
    USE MPI

    IMPLICIT NONE
    
    !IN/OUT
    Integer, INTENT(IN) :: M_Z, N_Z, NRHS_Z
    DOUBLE COMPLEX, Dimension(M_Z,N_Z), INTENT(IN) :: Zred
    DOUBLE COMPLEX, Dimension(M_Z,NRHS_Z), INTENT(IN) :: Vred  
    Integer, Dimension(MyNBlocks), INTENT(IN) :: K_patchs
    Integer, Dimension(NBlocks), INTENT(IN) :: K_patchs_all
    Integer, Dimension(NBlocks+1), INTENT(IN) :: Curs_Kpatchs_all
    Integer, Dimension(nber_procs), INTENT(IN) :: Ktot_procs
    Integer, Dimension(nber_procs,Nblk_proc_max+1), INTENT(IN):: MPI_CBFM_Blocks
    DOUBLE COMPLEX, Dimension(Mlocal,Nlocal), INTENT(INOUT) :: ZredLoc
    DOUBLE COMPLEX, Dimension(Mlocal,NRHSlocal), INTENT(INOUT) ::VredLoc
    
    ! local 
    Integer :: kk
    Integer :: iLoc, jLoc, iLoc_proc, iGlob, jGlob, block_iGlob, pos_loc_block
    INTEGER :: INDXL2G,INDXG2L, size1, size2, prev_job, next_job
    Integer :: Eff_rank, K_proc, NBlocks_eff
    INTEGER, dimension(MPI_STATUS_SIZE) :: status     
    INTEGER, dimension(:), allocatable :: K_patchs_eff
    
    
    ! here distribution of both Zred and Vred when needed
    
    !! Rank : Start with distributing my proper data---------------------------- 
    ! iLoc refers to line in ZredLoc, iGlob in total Zc and iLoc_proc in Zreduite (my band Zc for current proc)
    Do iLoc=1, Mlocal
        iGlob = INDXL2G(iLoc,M_B,Myrow,0,NPROW)
        block_iGlob = minloc(Curs_Kpatchs_all, 1, mask = Curs_Kpatchs_all .gt. iGlob) - 1
        pos_loc_block = minloc(MPI_CBFM_Blocks(rank+1,2:1+MyNBlocks), 1, mask = MPI_CBFM_Blocks(rank+1,2:1+MyNBlocks) .eq. block_iGlob);
        if ((pos_loc_block .gt. 1) .OR. ((pos_loc_block .eq. 1) .AND. (MPI_CBFM_Blocks(rank+1,2) .eq. block_iGlob))) then 
            iLoc_proc = sum(K_patchs(1:pos_loc_block-1))+ (iGlob - sum(K_patchs_all(1:block_iGlob-1)));
            ! ZredLoc
            Do jLoc=1, Nlocal
                jGlob = INDXL2G(jLoc,N_B,Mycol,0,NPCOL)
                ZredLoc(iLoc,jLoc) = Zred(iLoc_proc,jGlob)
            End do  
            ! VredLoc
	        Do jLoc=1, NRHSlocal
		        jGlob = INDXL2G(jLoc,N_B,Mycol,0,NPCOL)
                VredLoc(iLoc,jLoc) = Vred(iLoc_proc,jGlob)
            End do
        end if    
    End do  
    
    ! Now get the necessary to fill out my ZredLoc and VredLoc from the other processes 
    ! EST CE NORMAL QU'ICI JE NE VOIS PAS L'INFO sur l'origine de mon actuelle Zreduite ? 
    size1 = M_Z*N_Z ; ! Ktot_proc_max*K_total previously and Ktot_proc_max**2 now
    size2 = M_Z*NRHS_Z; !Ktot_proc_max*2*NTr
    
    
    ! identify previous and following MPI job 
    prev_job = mod(nber_procs+rank-1,nber_procs)
    next_job = mod(rank+1,nber_procs)
    Do kk = 1, nber_procs-1
        ! SEND/RECEIVE Zc from/to the prev/next MPI job
        !! Every job/process sends its 'Zc' to the previous process and receive the one of the next process     
        Call MPI_SENDRECV_REPLACE(Zred,2*size1,MPI_COMPLEX,prev_job,tag,next_job,tag,MPI_COMM_WORLD,status,code)
        Call MPI_SENDRECV_REPLACE(Vred,2*size2,MPI_COMPLEX,prev_job,tag,next_job,tag,MPI_COMM_WORLD,status,code)
                                
        !! ATTENTION : We should consider the change of the content of 'C_tot_patch_Zc_trans' through the different kk 
        !! So we should identify the effective job which created the CBFs included in the received C_trans_patchs for the current kk 
        Eff_rank = mod(rank+kk,nber_procs);            
        K_proc = Ktot_procs(Eff_rank+1);
        NBlocks_eff = MPI_CBFM_Blocks(Eff_rank+1,1)
        Allocate(K_patchs_eff(NBlocks_eff));
        K_patchs_eff = K_patchs_all(MPI_CBFM_Blocks(Eff_rank+1,2:1+NBlocks_eff));
        Do iLoc=1, Mlocal
    	iGlob = INDXL2G(iLoc,M_B,Myrow,0,NPROW)
        block_iGlob = minloc(Curs_Kpatchs_all, 1, mask = Curs_Kpatchs_all .gt. iGlob) - 1
        pos_loc_block = minloc(MPI_CBFM_Blocks(Eff_rank+1,2:1+NBlocks_eff), 1, mask = MPI_CBFM_Blocks(Eff_rank+1,2:1+NBlocks_eff) .eq. block_iGlob);
        if ((pos_loc_block .gt. 1) .OR. ((pos_loc_block .eq. 1) .AND. (MPI_CBFM_Blocks(Eff_rank+1,2) .eq. block_iGlob))) then
            iLoc_proc = sum(K_patchs_eff(1:pos_loc_block-1)) + (iGlob - sum(K_patchs_all(1:block_iGlob-1)));
            ! fill ZredLoc
        	Do jLoc=1, Nlocal
            	jGlob = INDXL2G(jLoc,N_B,Mycol,0,NPCOL);
            	ZredLoc(iLoc,jLoc) = Zred(iLoc_proc,jGlob);
        	End do  
        	! fill VredLoc
			Do jLoc=1, NRHSlocal
				jGlob = INDXL2G(jLoc,N_B,Mycol,0,NPCOL);
                VredLoc(iLoc,jLoc) = Vred(iLoc_proc,jGlob);
       		End do
        end if
        End do
        deallocate(K_patchs_eff);
    EndDo        
    END SUBROUTINE Block_Cyclic_Distribution_v1
    
    
    
    SUBROUTINE Block_Cyclic_Distribution_v2(M_Z,N_Z,NRHS_Z,Zred,Vred,dist_interaction,K_patchs,K_patchs_all,Curs_Kpatchs_all,Ktot_procs,MPI_CBFM_Blocks,ZredLoc,VredLoc)  
    
    USE Initialization
    USE common_variables
    USE lapack95
    USE MPI

    IMPLICIT NONE
    
    !IN/OUT
    Integer, INTENT(IN) :: M_Z, N_Z, NRHS_Z,dist_interaction
    DOUBLE COMPLEX, Dimension(M_Z,N_Z), INTENT(IN) :: Zred
    DOUBLE COMPLEX, Dimension(M_Z,NRHS_Z), INTENT(IN) :: Vred  
    Integer, Dimension(MyNBlocks), INTENT(IN) :: K_patchs
    Integer, Dimension(NBlocks), INTENT(IN) :: K_patchs_all
    Integer, Dimension(NBlocks+1), INTENT(IN) :: Curs_Kpatchs_all
    Integer, Dimension(nber_procs), INTENT(IN) :: Ktot_procs
    Integer, Dimension(nber_procs,Nblk_proc_max+1), INTENT(IN):: MPI_CBFM_Blocks
    DOUBLE COMPLEX, Dimension(Mlocal,Nlocal), INTENT(INOUT) :: ZredLoc
    DOUBLE COMPLEX, Dimension(Mlocal,NRHSlocal), INTENT(INOUT) ::VredLoc
    
    ! local 
    Integer :: kk,iLoc, iLoc_proc, iGlob, block_iGlob, posi_loc_block
    Integer :: NBlocks_j, jLoc, jLoc_proc,jGlob, block_jGlob, posj_loc_block
    INTEGER :: INDXL2G,INDXG2L, size1, size2, prev_job, next_job
    Integer :: Eff_rank, rank_j, K_proc, NBlocks_eff
    INTEGER, dimension(MPI_STATUS_SIZE) :: status     
    INTEGER, dimension(:), allocatable :: K_patchs_j,K_patchs_eff
    
    ! dist_interaction = 0 
    ! So here distribution of both Zred and Vred after calculation of proper interaction inside 
    ! each MPI job (this is what I mean by distance interaction = 0)
    if (dist_interaction .eq. 0) then 
        
        !! Rank : Start with distributing my proper data----------------------------
        rank_j = rank ; 
        NBlocks_j = MPI_CBFM_Blocks(rank_j+1,1);
        Allocate(K_patchs_j(NBlocks_j))
        K_patchs_j = K_patchs_all(MPI_CBFM_Blocks(rank_j+1,2:1+NBlocks_j));
        
        ! iLoc refers to line in ZredLoc, iGlob in total Zc and iLoc_proc in Zreduite (my band Zc for current proc)
        Do iLoc=1, Mlocal
            iGlob = INDXL2G(iLoc,M_B,Myrow,0,NPROW)
            block_iGlob = minloc(Curs_Kpatchs_all, 1, mask = Curs_Kpatchs_all .gt. iGlob) - 1
            posi_loc_block = minloc(MPI_CBFM_Blocks(rank+1,2:1+MyNBlocks), 1, mask = MPI_CBFM_Blocks(rank+1,2:1+MyNBlocks) .eq. block_iGlob);
            if ((posi_loc_block .gt. 1) .OR. ((posi_loc_block .eq. 1) .AND. (MPI_CBFM_Blocks(rank+1,2) .eq. block_iGlob))) then 
                iLoc_proc = sum(K_patchs(1:posi_loc_block-1))+ (iGlob - sum(K_patchs_all(1:block_iGlob-1)));
                ! ZredLoc
                Do jLoc=1, Nlocal
                    jGlob = INDXL2G(jLoc,N_B,Mycol,0,NPCOL)
                    block_jGlob = minloc(Curs_Kpatchs_all, 1, mask = Curs_Kpatchs_all .gt. jGlob) - 1
                    posj_loc_block = minloc(MPI_CBFM_Blocks(rank_j+1,2:1+NBlocks_j), 1, mask = MPI_CBFM_Blocks(rank_j+1,2:1+NBlocks_j) .eq. block_jGlob);
                    if ((posj_loc_block .gt. 1) .OR. ((posj_loc_block .eq. 1) .AND. (MPI_CBFM_Blocks(rank_j+1,2) .eq. block_jGlob))) then 
                        jLoc_proc = sum(K_patchs_j(1:posj_loc_block-1))+ (jGlob - sum(K_patchs_all(1:block_jGlob-1)));
                        ZredLoc(iLoc,jLoc) = Zred(iLoc_proc,jLoc_proc)
                    endif
                End do                
                
                ! VredLoc
	            Do jLoc=1, NRHSlocal
		            jGlob = INDXL2G(jLoc,N_B,Mycol,0,NPCOL)
                    VredLoc(iLoc,jLoc) = Vred(iLoc_proc,jGlob)
                End do
            end if    
        End do
        deallocate(K_patchs_j);
    
        ! Now get the necessary to fill out my ZredLoc and VredLoc from the other processes 
        ! EST CE NORMAL QU'ICI JE NE VOIS PAS L'INFO sur l'origine de mon actuelle Zreduite ? 
        size1 = M_Z*N_Z ; ! Ktot_proc_max*K_total previously and Ktot_proc_max**2 now
        size2 = M_Z*NRHS_Z; !Ktot_proc_max*2*NTr        
        
        ! identify previous and following MPI job 
        prev_job = mod(nber_procs+rank-1,nber_procs)
        next_job = mod(rank+1,nber_procs)
        Do kk = 1, nber_procs-1
            ! SEND/RECEIVE Zc from/to the prev/next MPI job
            !! Every job/process sends its 'Zc' to the previous process and receive the one of the next process     
            Call MPI_SENDRECV_REPLACE(Zred,2*size1,MPI_COMPLEX,prev_job,tag,next_job,tag,MPI_COMM_WORLD,status,code)
            Call MPI_SENDRECV_REPLACE(Vred,2*size2,MPI_COMPLEX,prev_job,tag,next_job,tag,MPI_COMM_WORLD,status,code)
                                
            !! ATTENTION : We should consider the change of the content of 'C_tot_patch_Zc_trans' through the different kk 
            !! So we should identify the effective job which created the CBFs included in the received C_trans_patchs for the current kk 
            Eff_rank = mod(rank+kk,nber_procs);            
            K_proc = Ktot_procs(Eff_rank+1);
            NBlocks_eff = MPI_CBFM_Blocks(Eff_rank+1,1)
            Allocate(K_patchs_eff(NBlocks_eff));
            K_patchs_eff = K_patchs_all(MPI_CBFM_Blocks(Eff_rank+1,2:1+NBlocks_eff));
          
            ! the Zii which is spinning is the result of the interaction of the eff-rank with itself (meaning du distance interaction = 0)
            rank_j = Eff_rank ; ! rank_j = Eff_rank + dist_interaction ! rank_j here is also an effective rank_j 
            NBlocks_j = MPI_CBFM_Blocks(rank_j+1,1);
            Allocate(K_patchs_j(NBlocks_j))
            K_patchs_j = K_patchs_all(MPI_CBFM_Blocks(rank_j+1,2:1+NBlocks_j));
        
            Do iLoc=1, Mlocal
    	        iGlob = INDXL2G(iLoc,M_B,Myrow,0,NPROW)
                block_iGlob = minloc(Curs_Kpatchs_all, 1, mask = Curs_Kpatchs_all .gt. iGlob) - 1
                posi_loc_block = minloc(MPI_CBFM_Blocks(Eff_rank+1,2:1+NBlocks_eff), 1, mask = MPI_CBFM_Blocks(Eff_rank+1,2:1+NBlocks_eff) .eq. block_iGlob);
                if ((posi_loc_block .gt. 1) .OR. ((posi_loc_block .eq. 1) .AND. (MPI_CBFM_Blocks(Eff_rank+1,2) .eq. block_iGlob))) then
                    iLoc_proc = sum(K_patchs_eff(1:posi_loc_block-1)) + (iGlob - sum(K_patchs_all(1:block_iGlob-1)));
                    ! fill ZredLoc
        	        Do jLoc=1, Nlocal
            	        jGlob = INDXL2G(jLoc,N_B,Mycol,0,NPCOL);
                        block_jGlob = minloc(Curs_Kpatchs_all, 1, mask = Curs_Kpatchs_all .gt. jGlob) - 1
                        posj_loc_block = minloc(MPI_CBFM_Blocks(rank_j+1,2:1+NBlocks_j), 1, mask = MPI_CBFM_Blocks(rank_j+1,2:1+NBlocks_j) .eq. block_jGlob);
                        if ((posj_loc_block .gt. 1) .OR. ((posj_loc_block .eq. 1) .AND. (MPI_CBFM_Blocks(rank_j+1,2) .eq. block_jGlob))) then
                            jLoc_proc = sum(K_patchs_j(1:posj_loc_block-1)) + (jGlob - sum(K_patchs_all(1:block_jGlob-1)));
                            ZredLoc(iLoc,jLoc) = Zred(iLoc_proc,jLoc_proc);
                        endif
        	        End do  
        	        ! fill VredLoc
			        Do jLoc=1, NRHSlocal
				        jGlob = INDXL2G(jLoc,N_B,Mycol,0,NPCOL);
                        VredLoc(iLoc,jLoc) = Vred(iLoc_proc,jGlob);
       		        End do
                endif
            Enddo
            deallocate(K_patchs_eff,K_patchs_j);
        EndDo    
    Else ! here the distribution after calculation of Z (Rank i <-> j), the distance intercation gives us an idea about the effective rank j
         ! for this part NO Vred is distributed since it is not concerned by the interactions inter-jobs 
        
        rank_j = mod(rank+dist_interaction,nber_procs)
        NBlocks_j = MPI_CBFM_Blocks(rank_j+1,1);
        Allocate(K_patchs_j(NBlocks_j))
        K_patchs_j = K_patchs_all(MPI_CBFM_Blocks(rank_j+1,2:1+NBlocks_j));
        
        !! Rank : Start with distributing my proper data---------------------------- 
        ! iLoc refers to line in ZredLoc, iGlob in total Zc and iLoc_proc in Zreduite (my band Zc for current proc)
        Do iLoc=1, Mlocal
            iGlob = INDXL2G(iLoc,M_B,Myrow,0,NPROW)
            block_iGlob = minloc(Curs_Kpatchs_all, 1, mask = Curs_Kpatchs_all .gt. iGlob) - 1
            posi_loc_block = minloc(MPI_CBFM_Blocks(rank+1,2:1+MyNBlocks), 1, mask = MPI_CBFM_Blocks(rank+1,2:1+MyNBlocks) .eq. block_iGlob);
            if ((posi_loc_block .gt. 1) .OR. ((posi_loc_block .eq. 1) .AND. (MPI_CBFM_Blocks(rank+1,2) .eq. block_iGlob))) then 
                iLoc_proc = sum(K_patchs(1:posi_loc_block-1))+ (iGlob - sum(K_patchs_all(1:block_iGlob-1)));
                ! ZredLoc
                Do jLoc=1, Nlocal
                    jGlob = INDXL2G(jLoc,N_B,Mycol,0,NPCOL)
                    block_jGlob = minloc(Curs_Kpatchs_all, 1, mask = Curs_Kpatchs_all .gt. jGlob) - 1
                    posj_loc_block = minloc(MPI_CBFM_Blocks(rank_j+1,2:1+NBlocks_j), 1, mask = MPI_CBFM_Blocks(rank_j+1,2:1+NBlocks_j) .eq. block_jGlob);
                    if ((posj_loc_block .gt. 1) .OR. ((posj_loc_block .eq. 1) .AND. (MPI_CBFM_Blocks(rank_j+1,2) .eq. block_jGlob))) then 
                        jLoc_proc = sum(K_patchs_j(1:posj_loc_block-1))+ (jGlob - sum(K_patchs_all(1:block_jGlob-1)));
                        ZredLoc(iLoc,jLoc) = Zred(iLoc_proc,jLoc_proc)
                    endif
                End do                
            end if    
        End do 
        deallocate(K_patchs_j);
    
        ! Now get the necessary to fill out my ZredLoc and VredLoc from the other processes 
        ! EST CE NORMAL QU'ICI JE NE VOIS PAS L'INFO sur l'origine de mon actuelle Zreduite ? 
        size1 = M_Z*N_Z ; ! Ktot_proc_max*K_total previously and Ktot_proc_max**2 now
           
        ! identify previous and following MPI job 
        prev_job = mod(nber_procs+rank-1,nber_procs)
        next_job = mod(rank+1,nber_procs)
        Do kk = 1, nber_procs-1
            ! SEND/RECEIVE Zc from/to the prev/next MPI job
            !! Every job/process sends its 'Zc' to the previous process and receive the one of the next process     
            Call MPI_SENDRECV_REPLACE(Zred,2*size1,MPI_COMPLEX,prev_job,tag,next_job,tag,MPI_COMM_WORLD,status,code)
                                
            !! ATTENTION : We should consider the change of the content of 'C_tot_patch_Zc_trans' through the different kk 
            !! So we should identify the effective job which created the CBFs included in the received C_trans_patchs for the current kk 
            Eff_rank = mod(rank+kk,nber_procs);            
            K_proc = Ktot_procs(Eff_rank+1);
            NBlocks_eff = MPI_CBFM_Blocks(Eff_rank+1,1)
            Allocate(K_patchs_eff(NBlocks_eff));
            K_patchs_eff = K_patchs_all(MPI_CBFM_Blocks(Eff_rank+1,2:1+NBlocks_eff));
          
            rank_j = mod(Eff_rank+dist_interaction,nber_procs)
            NBlocks_j = MPI_CBFM_Blocks(rank_j+1,1);
            Allocate(K_patchs_j(NBlocks_j))
            K_patchs_j = K_patchs_all(MPI_CBFM_Blocks(rank_j+1,2:1+NBlocks_j));          
          
            Do iLoc=1, Mlocal
    	        iGlob = INDXL2G(iLoc,M_B,Myrow,0,NPROW)
                block_iGlob = minloc(Curs_Kpatchs_all, 1, mask = Curs_Kpatchs_all .gt. iGlob) - 1
                posi_loc_block = minloc(MPI_CBFM_Blocks(Eff_rank+1,2:1+NBlocks_eff), 1, mask = MPI_CBFM_Blocks(Eff_rank+1,2:1+NBlocks_eff) .eq. block_iGlob);
                if ((posi_loc_block .gt. 1) .OR. ((posi_loc_block .eq. 1) .AND. (MPI_CBFM_Blocks(Eff_rank+1,2) .eq. block_iGlob))) then
                    iLoc_proc = sum(K_patchs_eff(1:posi_loc_block-1)) + (iGlob - sum(K_patchs_all(1:block_iGlob-1)));
                    ! fill ZredLoc
        	        Do jLoc=1, Nlocal
            	        jGlob = INDXL2G(jLoc,N_B,Mycol,0,NPCOL);
                        block_jGlob = minloc(Curs_Kpatchs_all, 1, mask = Curs_Kpatchs_all .gt. jGlob) - 1
                        posj_loc_block = minloc(MPI_CBFM_Blocks(rank_j+1,2:1+NBlocks_j), 1, mask = MPI_CBFM_Blocks(rank_j+1,2:1+NBlocks_j) .eq. block_jGlob);
                        if ((posj_loc_block .gt. 1) .OR. ((posj_loc_block .eq. 1) .AND. (MPI_CBFM_Blocks(rank_j+1,2) .eq. block_jGlob))) then
                            jLoc_proc = sum(K_patchs_j(1:posj_loc_block-1)) + (jGlob - sum(K_patchs_all(1:block_jGlob-1)));
                            ZredLoc(iLoc,jLoc) = Zred(iLoc_proc,jLoc_proc);
                        endif
        	        End do  
                end if
            End do
            deallocate(K_patchs_eff,K_patchs_j);
        EndDo        
    EndIf 
        
    END SUBROUTINE Block_Cyclic_Distribution_v2
    
    
    SUBROUTINE numroc_2d_opt(M,N,NRHS,NROWS,NCOLS,Iproc,Jproc,Mloc,Nloc,NRHSloc) 
    
      USE Initialization
      USE common_variables
      USE lapack95
      USE MPI
  
      IMPLICIT NONE
      
      !IN/OUT
      Integer, INTENT(IN) :: M,N,NRHS,NROWS,NCOLS,Iproc,Jproc
      Integer, INTENT(OUT) :: Mloc,Nloc, NRHSloc
      
      ! Local 
      Integer :: Md, Nd, M0, N0
      Integer :: rest_M,rest_M_d,rest_M_m,rest_N,rest_N_d,rest_N_m
      
      ! Remember : M_B & N_B the scalapack block sizes are defined in Initialization
      
      ! Note that this "2d" version is equivalent to calling 3 times 
      ! a simpler numroc_opt (still different from the scalapack numroc) with only three parameters (M .OR. N .OR. NRHS, Mrows .OR. MCols, Iproc .OR. Jproc)
      ! I've chosen to do the three in the same subroutine even if this means repeating the same lines ! 
      
      Md =  M/(M_B*NROWS);      
      ! this is the equal part that all the processors will have
      M0 = Md*M_B;          
      ! then the rest of M/N is smaller than M_A*NROWS/N_A*NCOLS 
      ! so will not be distributed over all the processors
      rest_M = M - (Md*(M_B*NROWS));
      rest_M_d = rest_M/M_B;  ! procs with Iproc < rest_M_d (equivalent to Iproc +1 <= rest_M_div) will receive an extra M_B
      rest_M_m = mod(rest_M,M_B); ! the proc with excatly Iproc = rest_M_d 
      if (Iproc .lt. rest_M_d) then
        Mloc = M0 + M_B;      
      elseif (Iproc .eq. rest_M_d) then 
        Mloc = M0 + rest_M_m;
      else   
        Mloc = M0 ;        
      endif       
      
      ! same for Nloc = f(N, N_B, NROWS, Jproc) 
      Nd =  N/(N_B*NCOLS);
      N0 = Nd*N_B;
      
      rest_N = N - (Nd*(N_B*NCOLS));          
      rest_N_d = rest_N/N_B;  
      rest_N_m = mod(rest_N,N_B);  
      if (Jproc .lt. rest_N_d) then
        Nloc = N0 + N_B;      
      elseif (Jproc .eq. rest_N_d) then 
        Nloc = N0 + rest_N_m;
      else   
        Nloc = N0 ;        
      endif   
      
      ! Now NRHS same as Nloc just replace N with NRHS and Nloc with NRHSloc  
      Nd =  NRHS/(N_B*NCOLS);
      N0 = Nd*N_B;
      
      rest_N = NRHS - (Nd*(N_B*NCOLS));          
      rest_N_d = rest_N/N_B;  
      rest_N_m = mod(rest_N,N_B);  
      if (Jproc .lt. rest_N_d) then
        NRHSloc = N0 + N_B;      
      elseif (Jproc .eq. rest_N_d) then 
        NRHSloc = N0 + rest_N_m;
      else   
        NRHSloc = N0 ;        
      endif      
    END SUBROUTINE numroc_2d_opt