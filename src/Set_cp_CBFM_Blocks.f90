SUBROUTINE Set_cp_CBFM_Blocks(Cells,CBFM_Blocks,CBFM_Blocks_Ext,cp_CBFM_Blocks) 
    
    USE Initialization
    USE common_variables
    USE iso_fortran_env
        
    IMPLICIT NONE
    
    !! IN/OUT ******************************************************************
    type (Cell), Dimension(Nbc), INTENT(IN) :: Cells
    type (CBFM_Block), Dimension(Nblocks), INTENT(IN) :: CBFM_Blocks
    Integer, Dimension(Nblocks,Nbc_ext), INTENT(IN) :: CBFM_Blocks_Ext
    Integer, Dimension(:,:), allocatable, INTENT(OUT) :: cp_CBFM_Blocks
    
    !! Local 
    Integer :: ii,jj,m,n,Check_Copies,newBlock_curs,Nbc_b,Nbc_b_ext    
    Integer, Dimension(:), allocatable :: BlkintN,BlkintNext
    Integer, Dimension(:,:), allocatable :: cp_CBFM_Blocks_tmp
    
    ! I will use this param for the moment to check the efficiency and accuracy of the definition and use of the 
    ! copies of the calculated blocks instead of calculating the CBFs for all the blocks 
    Check_Copies = 1;
    
    ! Here we recover two major informations : 'NcalBlks' and the array 'cp_CBFM_Blocks' 
    if ((Check_Copies .eq. 1) .and. (NbintBl .gt. 1)) then 
        Allocate(BlkintN(NbintBl),BlkintNext(NbintBl));
        Allocate(cp_CBFM_Blocks(2+NbintBl,NbintBl));
        cp_CBFM_Blocks = 0;
        
        BlkintN(1:NbintBl) = CBFM_Blocks(1:NbintBl)%Nbc_b;
        BlkintNext(1:NbintBl) = CBFM_Blocks(1:NbintBl)%Nbc_ext;
        newBlock_curs = 0;
        do ii=1,NbintBl
            if (BlkintN(ii) .ne. 0) then 
                ! initialization since this block is a new one 
                newBlock_curs = newBlock_curs + 1
                cp_CBFM_Blocks(1,newBlock_curs) = ii;
                cp_CBFM_Blocks(2,newBlock_curs) = 1; 
                cp_CBFM_Blocks(3,newBlock_curs) = ii; 
                ! now check if there are other blocks that are just copie of this one 
                Nbc_b = BlkintN(ii); BlkintN(ii) = 0; ! Recover BlkintN(ii) then put it to 0 
                Nbc_b_ext= BlkintNext(ii);     
                
		        jj=0;
                do while (jj .ne. ii)
                    jj = minloc(BlkintN(ii:NbintBl),1,mask=(BlkintN(ii:NbintBl) == Nbc_b));
                    jj = ii+jj-1;
                    !if (BlkintNext(jj) .eq. Nbc_b_ext) then
                     if (jj .ne. ii) then
                        BlkintN(jj) = 0;
                        BlkintNext(jj) = 0;  
                        cp_CBFM_Blocks(2,newBlock_curs) = cp_CBFM_Blocks(2,newBlock_curs) + 1;
                        m = cp_CBFM_Blocks(2,newBlock_curs);
                        cp_CBFM_Blocks(2+m,newBlock_curs) = jj; 
                     endif
                    !endif
                enddo   
            endif            
        enddo
        
        Deallocate(BlkintN,BlkintNext);
        Nccp_max = maxval(cp_CBFM_Blocks(2,:))
        m = 2 + Nccp_max; 
        n = newBlock_curs;
        Allocate(cp_CBFM_Blocks_tmp(m,n));
        cp_CBFM_Blocks_tmp(1:m,1:n) = cp_CBFM_Blocks(1:m,1:n);
        deallocate(cp_CBFM_Blocks);
        
        NcalBlks = n+NBlocks-NbintBl;
        Allocate(cp_CBFM_Blocks(m,NcalBlks));
        cp_CBFM_Blocks = 0;
        cp_CBFM_Blocks(1:m,1:n)=cp_CBFM_Blocks_tmp(1:m,1:n);
        cp_CBFM_Blocks(1,n+1:NcalBlks) = [NbintBl+1:NBlocks];
        cp_CBFM_Blocks(2,n+1:NcalBlks) = 1;
        cp_CBFM_Blocks(3,n+1:NcalBlks) = [NbintBl+1:NBlocks]; 
        deallocate(cp_CBFM_Blocks_tmp);
        Write(*,'(a,i4,a,i4,a)') 'Note that NcalBlks = ',NcalBlks,' out of total',NBlocks,' blocks';
    else
        Allocate(cp_CBFM_Blocks(3,NBlocks));
        NcalBlks = NBlocks;
        Nccp_max = 1;
        cp_CBFM_Blocks(1,1:NBlocks) = [1:NBlocks];
        cp_CBFM_Blocks(2,1:NBlocks) = 1;
        cp_CBFM_Blocks(3,1:NBlocks) = [1:NBlocks]; 
    endif
    
END SUBROUTINE Set_cp_CBFM_Blocks