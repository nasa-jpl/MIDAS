SUBROUTINE MPI_distribution_blocks(CBFM_Blocks,MPI_CBFM_Blocks)

    ! HERE Ditribution of the CBFM blocks among the available MPI jobs       
        
    USE Initialization
    USE common_variables
    USE iso_fortran_env
    USE DiverseUtil
    USE MPI
    
    Implicit NONE
    
    !IN/OUT 
    type (CBFM_Block), Dimension(Nblocks), INTENT(IN) :: CBFM_Blocks
    Integer, Dimension(:,:), allocatable, INTENT(INOUT) :: MPI_CBFM_Blocks
    
    ! local 
    Integer :: ii,jj,Nblks_min_proc,n_rest,N,nb!,Doption
    Integer :: curs_block,blks_round,pp,bb,nbb,pp_s,pp_e
    Integer, dimension(nber_procs) :: procs_nbc_sort,procs_nbc, p_sort_out
    Integer, dimension(:), allocatable :: blocks_nbc,b_sort_out
    Integer, dimension(:,:), allocatable ::MPI_CBFM_Blocks_tmp
    
    !Doption = 2; ! Distribution option1 guaranty successive blocks for each proc not option 2
    
    Allocate(MPI_CBFM_Blocks(nber_procs,Nblocks)); ! TMP
    MPI_CBFM_Blocks = 0;
    
    If (nber_procs .eq. 1) Then ! inefficient but to consider
        Write(*,'(a)') 'Please Note that Nb_procs == 1 !! ';
        MPI_CBFM_Blocks(1,1) = Nblocks;
        MPI_CBFM_Blocks(1,2:2+Nblocks-1) = (/(ii, ii=1,Nblocks)/);
         
    ElseIf (Nblocks .gt. nber_procs) Then 
      Nblks_min_proc = Nblocks/nber_procs;            
      !if (Doption == 1) then
      !    ! we start by attributing Nblks_min_proc to each proc
      !    MPI_CBFM_Blocks(1:nber_procs,1) = Nblks_min_proc;   
      !    ! then the procs which have the 'n_rest' lower Nbc tot will have one extra block each. Very simple !!!
      !    n_rest = Nblocks-Nblks_min_proc*nber_procs
      !    if (n_rest .ne. 0) then 
      !        !do pp=1,nber_procs
      !        !    pp_s = Nblks_min_proc*(pp-1)+1;
      !        !    pp_e = Nblks_min_proc*pp
      !        !    procs_nbc(pp) = sum(CBFM_Blocks(pp_s:pp_e)%Nbc_b)+sum(CBFM_Blocks(pp_s:pp_e)%Nbc_ext); 
      !        !enddo 
      !        !Allocate(b_sort_out(nber_procs))
      !        !b_sort_out = (/(ii, ii=1,nber_procs)/); 
      !        !call Sort_asc(procs_nbc, nber_procs,b_sort_out);
      !        MPI_CBFM_Blocks(1:n_rest,1) = Nblks_min_proc + 1; ! !! pour le moment il n'y aucun grand effort pour bien balancer la charge !!!
      !                                                          ! la solution commentee est fausse, il faut chercher un efacon efficace pour distribuer
      !                                                          ! equitablement des blocs successifs !
      !        !! tu peux peut etre penser plus tard a re-organiser les cellules une fois les blocks equitablement ditribues ? cad distribuer sans prendre en compte
      !        !!le successif puis re-reorganiser les cellules !
      !        
      !        !deallocate(b_sort_out);                      
      !    endif 
      !    ! once I know the number of blocks per proc, I attribute successive blocks to each proc and re-calculate procs_nbc. Simple !!
      !    do pp=1,nber_procs
      !      nb = MPI_CBFM_Blocks(pp,1);
      !      pp_s = sum(MPI_CBFM_Blocks(1:pp-1,1))+1;
      !      pp_e = sum(MPI_CBFM_Blocks(1:pp,1));
      !      MPI_CBFM_Blocks(pp,2:2+nb-1) =  (/(ii, ii=pp_s,pp_e)/);
      !      procs_nbc(pp) = sum(CBFM_Blocks(pp_s:pp_e)%Nbc_b); 
      !    enddo     
      !    
      !elseif (Doption == 2) then 
          ! initialize MPI_CBFM_Blocks 
          MPI_CBFM_Blocks(1:nber_procs,1) = 1;  
          
          Allocate(blocks_nbc(Nblocks),b_sort_out(Nblocks))
          blocks_nbc(1:Nblocks) = CBFM_Blocks(1:Nblocks)%Nbc_b+CBFM_Blocks(1:Nblocks)%Nbc_ext; 
      
          b_sort_out = (/(ii, ii=1,Nblocks)/); 
          call Sort_desc(blocks_nbc, Nblocks,b_sort_out);
          
          MPI_CBFM_Blocks(1:nber_procs,2) = CBFM_Blocks(b_sort_out(1:nber_procs))%num_block; 
          procs_nbc(1:nber_procs) = blocks_nbc(1:nber_procs);
      
          n_rest = Nblocks-nber_procs
          
          ! Now for each round give the largest block (from the rest of b_sort_out) to the processor with the least work
          ! for each round each process takes only one block
          Do ii = 1,n_rest
            p_sort_out = (/(jj, jj=1,nber_procs)/);
            procs_nbc_sort = procs_nbc;
            call Sort_asc(procs_nbc_sort, nber_procs,p_sort_out)        
            
            pp = p_sort_out(1);
            MPI_CBFM_Blocks(pp,1) = MPI_CBFM_Blocks(pp,1) + 1;
            nbb = MPI_CBFM_Blocks(pp,1)
            bb = b_sort_out(nber_procs+ii)
            MPI_CBFM_Blocks(pp,1+nbb) = bb
            procs_nbc(pp) = procs_nbc(pp) + CBFM_Blocks(bb)%Nbc_b + CBFM_Blocks(bb)%Nbc_ext; 
          EndDo
      !endif
    Else
      MPI_CBFM_Blocks(1:Nblocks,1) = 1; 
      MPI_CBFM_Blocks(Nblocks+1:nber_procs,1) = 0;
      MPI_CBFM_Blocks(1:Nblocks,2) = CBFM_Blocks(1:Nblocks)%num_block;
    EndIf 
    
    ! Nblk_proc_max needed to know the size of MPI_CBFM_Blocks
    ! We cannot keep only Nblk_proc since every proc needs to know what all
    ! the other procs have also as Nblk_proc 
    Nblk_proc_max = maxval(MPI_CBFM_Blocks(1:nber_procs,1));
    Allocate(MPI_CBFM_Blocks_tmp(nber_procs,Nblk_proc_max+1))
    MPI_CBFM_Blocks_tmp(1:nber_procs,1:Nblk_proc_max+1) = MPI_CBFM_Blocks(1:nber_procs,1:Nblk_proc_max+1);
    Deallocate(MPI_CBFM_Blocks);
    Allocate(MPI_CBFM_Blocks(nber_procs,Nblk_proc_max+1))
    MPI_CBFM_Blocks = MPI_CBFM_Blocks_tmp;
    deallocate(MPI_CBFM_Blocks_tmp);    
    
    ! Compute Nbc_proc (COMMON VARIABLE) 
    N = MPI_CBFM_Blocks(rank+1,1);
    Nbc_proc = 0; 
    Do ii=1, N
      bb = MPI_CBFM_Blocks(rank+1,1+ii)
      Nbc_proc = Nbc_proc + CBFM_Blocks(bb)%Nbc_b;       
    EndDo 
    call MPI_BARRIER(MPI_COMM_WORLD,code); 
    !Write(*,*) 'job ', rank, '; Nbc = ',Nbc_proc,' : ',MPI_CBFM_Blocks(rank+1,2:Nblk_proc_max+1);
    
END SUBROUTINE MPI_distribution_blocks