SUBROUTINE Division_blocks(SimScatterer,Cells,Ncells_SphDomains,CBFM_Blocks,MLCBFM_BlDistr)
    
    USE Initialization
    USE common_variables
    USE MPI
    
    !IN/OUT 
    type (Scatterer), INTENT(INOUT) :: SimScatterer
    type (Cell), Dimension(Nbc), INTENT(INOUT) :: Cells
    Integer, Dimension(7), INTENT(IN) :: Ncells_SphDomains
    type (CBFM_Block), Dimension(:), allocatable, INTENT(OUT) :: CBFM_Blocks
    Integer, Dimension(:,:), allocatable, INTENT(OUT):: MLCBFM_BlDistr
    
    ! local 
    integer :: Type_Par, err,m, I, N
    real(kind=8) :: dmin, dmax
    Integer, Dimension(:), allocatable :: Diff_avg
    Real(kind=8), Dimension(:), allocatable :: hB_test
    
    INTERFACE        
        SUBROUTINE Division_blocks_csh(SimScatterer,Cells,CBFM_Blocks,MLCBFM_BlDistr,error_division)
        
            USE Initialization
            USE common_variables
            USE f95_precision
            Implicit NONE
            
            !IN/OUT 
            type (Scatterer), INTENT(INOUT) :: SimScatterer
            type (Cell), Dimension(Nbc), INTENT(INOUT) :: Cells
            type (CBFM_Block), Dimension(:), allocatable, INTENT(OUT) :: CBFM_Blocks
            Integer, Dimension(:,:), allocatable, INTENT(OUT):: MLCBFM_BlDistr
            Integer, INTENT(INOUT) :: error_division            
        END SUBROUTINE Division_blocks_csh
        
        SUBROUTINE Division_blocks_sph(SimScatterer,Cells,Ncells_SphDomains,CBFM_Blocks,MLCBFM_BlDistr,error_division)

            ! The division into blocks depends on the type of the considered scatterer
            ! It is much simpler for the conventional shapes : Sphere, Cylinder ...
            ! The division here is for the complex geometries (from file) 
            ! STILL can be improved ...        
                
            USE Initialization
            USE common_variables
            USE f95_precision
            Implicit NONE
            
            !IN/OUT 
            type (Scatterer), INTENT(INOUT) :: SimScatterer
            type (Cell), Dimension(Nbc), INTENT(INOUT) :: Cells
            Integer, Dimension(7), INTENT(IN) :: Ncells_SphDomains
            type (CBFM_Block), Dimension(:), allocatable, INTENT(OUT) :: CBFM_Blocks
            Integer, Dimension(:,:), allocatable, INTENT(OUT):: MLCBFM_BlDistr
            Integer, INTENT(INOUT) :: error_division               
        END SUBROUTINE Division_blocks_sph 
    
    END INTERFACE
    
    Type_Par = SimScatterer%type_s 
    
    
    if ((Type_Par .eq. 1) .or. (Type_Par .eq. 6)) then              
        Call Division_blocks_sph(SimScatterer,Cells,Ncells_SphDomains,CBFM_Blocks,MLCBFM_BlDistr,err);
                
        if (err .ne. 0) then 
            Write(*,'(a)') 'Something went wrong when dividing into blocks !'
            stop 1;
        endif
                
        hBlock = 2*maxval(CBFM_Blocks(:)%BSphCont(4));
    else
        !! ADAPTIVE DIVISION INTO BLOCKS TO GET AS CLOSE AS POSSIBLE 
        !!TO NBlock (avgNcells) indicated in 'Simulation_data.dat'   
        !! Division into M blocks of height hBlock 
        !! To get closer to the average NBlock indicated in simulation_data, 
        !! hBlock = hBlock*hB_test_step then Redivide 
        dmin = min(SimScatterer%dx,SimScatterer%dy,SimScatterer%dz);
        dmax = max(SimScatterer%dx,SimScatterer%dy,SimScatterer%dz);
        hBlock = dmax;
                
        !! modif to transfer to MPI code
        m = 8; ! number of trials
        Allocate(Diff_avg(m+1),hB_test(m+1)); 
        Diff_avg(:) = Nbc; ! Initialization to avoid min = 0
        I=1; N = Nbc; 
        err = 0;
        ! the second condition is to stop trying depending on Navg_cells
        Do while ((I .le. m) .and. (N .gt. Navg_cells)) 
            Call Division_blocks_csh(SimScatterer,Cells,CBFM_Blocks,MLCBFM_BlDistr,err);
            if ((err .ne. 0) .AND. (I .eq. 1)) then
                Write(*,'(a)') 'Something went wrong in the division into blocks !!';
                STOP 1;
            elseif (err .ne. 0) then 
                exit ;
            endif
            hB_test(I) = hBlock; 
            N = sum(CBFM_Blocks(:)%Nbc_b)/Nblocks;
            Diff_avg(I) = abs(N-Navg_cells);
            ! next test 
            hBlock = hBlock*hB_test_step; 
            I = I+1;
        EndDo
               
        if (Nbc .le. 2*Navg_cells) then  ! if Nbc < 2*Navg_cells the scatterer is simply divided to 2 blocks
            hBlock = dmax;
        else              
            ! to obtain the final division into blocks, we recover hBlock which resulted in the minimum Diff_avg
            hBlock = hB_test(minloc(Diff_avg,1));           
        endif
        err = 0;
        Call Division_blocks_csh(SimScatterer,Cells,CBFM_Blocks,MLCBFM_BlDistr,err)
        ! Finally, we save the maximum height resulting from 
        ! the division into blocks in /Param_MoMCBFM/ hBlock
        hBlock = 2*maxval(CBFM_Blocks(:)%BSphCont(4));     
        Deallocate(Diff_avg,hB_test);                  
    EndIf 

    
END SUBROUTINE Division_blocks
    
    

    
SUBROUTINE Division_blocks_csh(SimScatterer,Cells,CBFM_Blocks,MLCBFM_BlDistr,error_division)

    ! The division into blocks depends on the type of the considered scatterer
    ! It is much simpler for the conventional shapes : Sphere, Cylinder ...
    ! The division here is for the complex geometries (from file) 
    ! STILL can be improved ...        
        
    USE Initialization
    USE common_variables
        
    Implicit NONE
    
    !IN/OUT 
    type (Scatterer), INTENT(INOUT) :: SimScatterer
    type (Cell), Dimension(Nbc), INTENT(INOUT) :: Cells
    type (CBFM_Block), Dimension(:), allocatable, INTENT(OUT) :: CBFM_Blocks
    Integer, Dimension(:,:), allocatable, INTENT(OUT):: MLCBFM_BlDistr
    Integer, INTENT(INOUT) :: error_division
    
    ! local 
    Integer :: ii,jj,bb,Type_Par,Nbc_bef,Nbc_aft, cel_part
    Integer :: ix,iy,iz,iB, div, NbBdir0, NbBdir1or2
    Integer :: NbBx, NbBy, NbBz, NbBl, Nb_block
    Integer :: idir0,curs_glb_Blk,cel,cel_init,cel_final
    Integer :: NcB_max, cc,Nbcels,ncp,numB_cel
    Integer :: num_cel,curs_new_cel,curs_new_cel_part,Bl_init,numB_cel_gl,Nbc_block
    Integer :: Nb_LargBlock,Nbmax,Nbmin,iiB, new_iiB,dir_div,new_NbBl
    Integer :: condExit,new_Nb_block,addedB,iiBend,up_floor,curs_blk
    Integer :: md,Nbc_toDist_B,Nbc_perBlk,cursB,currentB,numB_gl
    Integer :: dir0, dir1, dir2, pos_max, step_dir,ind_max,nbc1,nbc2,ref_error,tr,NbTr
    Integer :: NbBlL2,GB, Nbc_GB
    Real(kind=8) :: hbox,hB,hx,hy,hz, xc,yc,zc,Sc
    Real(kind=8) :: d_infperc,d_supperc,Rd,Rd_fn,error1,error2,perc_err
    Real(kind=8) :: hBldir0,hBldir1or2,hBldir1or2_init,dmin,hbySc,dstep,step_tr
    Real(kind=8) :: prev_error,prev_hBldir1or2
    Real(kind=8), dimension(:,:), allocatable :: positions
    Real(kind=8), dimension(:), allocatable :: positions_divdir
    Real(kind=8), dimension(3) :: pr_dims
    character(200) :: file_name
    character :: info_p_fl
    logical :: Div_is_possible,CanDiv,last_division,switchToMean
    type (Cell), Dimension(:), allocatable:: New_cells,cells_Block,new_cells_Block
    type (CBFM_Block), Dimension(:), allocatable:: CBFM_Blocks_tmp,CBFM_Blocks_p_tmp,CBFM_Blocks_p
    Integer, dimension(:,:), allocatable :: cells_in_blocks
    Integer, dimension(:), allocatable :: order_div, Nbc_Level2
    Integer, dimension(3) :: dist_grid
    logical, dimension(3) :: mask
    Integer, Dimension(1) :: ind_shift
   
    error_division = 0;
    
    !********************************************************************************
    !**********************GEOMETRY BASED DIVISION INTO BLOCKS***********************
    !********************************************************************************
    ! let us allocate a primary CBFM_Blocks and a new Cells to put in the blocks and 
    ! the reorganized cells while scanning the scattterer
    Allocate(CBFM_Blocks_tmp(NbBlock_s_max))
    Allocate(New_Cells(Nbc)); ! the cells will be reorganized depending in the position of the block
                                ! to which they belong
    curs_new_cel_part = 1;
    curs_new_cel = 1;
    Nblocks = 0   ! here we initialize the total number of block
    curs_glb_Blk = 1; ! index block
    iB = 0;
    
    !! THIS CODE IS NOW INTENDED FOR SINGLE SCATTERING PROPERTIES 
    Type_Par = SimScatterer%type_s   
    info_p_fl = SimScatterer%info_s    
      
    ncp = Nbc;
    Sc = SimScatterer%Sc;
        
    ! here we consider only the cells composing the current scatterer
    Allocate(CBFM_Blocks_p_tmp(NbBlock_s_max)); 
    Allocate(positions(3,ncp)); 
    positions(1,1:ncp) = Cells(1:ncp)%Xc
    positions(2,1:ncp) = Cells(1:ncp)%Yc
    positions(3,1:ncp) = Cells(1:ncp)%Zc      
            
      
    !! The idea is to identify the direction along which the scatterer has its maximum size --> dir 0
    !! we divide the scatterer along dir 0 
    !! and then in an iterative way, divide along dir1 and dir 2 (the other two directions) 
    !! depending on the maxim size of the block to divide !
    pr_dims(1) = SimScatterer%dx
    pr_dims(2) = SimScatterer%dy
    pr_dims(3) = SimScatterer%dz
    Allocate(order_div(3));
    Do jj=1,3
    pos_max = minloc(pr_dims,1,pr_dims==maxval(pr_dims))  
    pr_dims(pos_max) = 0.
    order_div(jj) = pos_max       
    EndDo
    pr_dims(1) = SimScatterer%dx; pr_dims(2) = SimScatterer%dy
    pr_dims(3) = SimScatterer%dz
      
    !Write(*,*) 'order_div =',order_div
    !! *****************************************************************************
    ! *************************** Start dividing along dir0 ***************************  
    !! *****************************************************************************
      
    dir0 = order_div(1)  
    hbox = pr_dims(dir0); 
    div = hbox/hBlock +1; NbBdir0 = div;  
    hBldir0 = hbox/div
    NbBl = NbBdir0 
    
    ! Order/Arrange the current CBFM blocks 
    ! and put them in CBFM_Blocks_p
    Do idir0=1,NbBdir0
        CBFM_Blocks_p_tmp(idir0)%num_block = curs_glb_Blk+idir0-1;  
    EndDo           

    Allocate(CBFM_Blocks_p(NbBl))
    CBFM_Blocks_p(:) = CBFM_Blocks_p_tmp(1:NbBl) 
    if (NbBl .gt. 1) Then 
        ! Now that the blocks are defined and ordred, we scan the cells 
        ! to distribute them 
        Allocate(cells_in_blocks(NbBl,ncp+1))
        if (dir0==1) Then
        dmin = SimScatterer%xmin; 
        elseif (dir0==2) Then
        dmin = SimScatterer%ymin; 
        Else
        dmin = SimScatterer%zmin; 
        EndIf
          
        !dstep = hBldir0;
        hbySc = hBldir0/Sc; 
        dstep = nint(hbySc)*Sc;  
        Allocate(positions_divdir(ncp))
        positions_divdir = positions(dir0,1:ncp)
        Call distributeCells(ncp,NbBl,positions_divdir,dmin,dstep,Cells,cells_in_blocks)
        Deallocate(positions_divdir);
        ! Now reorganize the cells depending on the block to which they belong 
        ! that will make the application of the CBFM and the extension of the blocks easier
        curs_new_cel =1;        
        Do numB_cel=1,NbBl
            numB_gl = curs_glb_Blk + numB_cel - 1;
            Nbcels = cells_in_blocks(numB_cel,1)
            CBFM_Blocks_p(numB_cel)%Nbc_b= Nbcels
            Do num_cel=1,Nbcels
                cel = cells_in_blocks(numB_cel,1+num_cel)
                New_Cells(curs_new_cel) = Cells(cel)
                New_Cells(curs_new_cel)%num_block = numB_gl;       
                curs_new_cel = curs_new_cel + 1
            EndDo            
        EndDo
          
        Do num_cel=1,ncp
            New_Cells(num_cel)%num_cell = num_cel                    
        EndDo
        Cells(1:ncp) = New_Cells(1:ncp);
        Deallocate(cells_in_blocks) 
    Else
        CBFM_Blocks_p(1)%Nbc_b = ncp;              
    EndIf 
    Deallocate(positions)       
      
    !! **********************HERE PREPARE for simple 2Levels MLCBFM ****************
    !! *****************************************************************************
    If (MLCBFM ==1) Then
    NbBlL2 = NbBl;
    If (NbBl .ge. 4) then !4 blocks or more
        NbBlksL2 = NbBl-2;
        Allocate(Nbc_Level2(NbBlksL2)) 
        Nbc_Level2(1) = sum(CBFM_Blocks_p(1:2)%Nbc_b); 
        Nbc_Level2(2:NbBl-3) = CBFM_Blocks_p(3:NbBl-2)%Nbc_b 
        Nbc_Level2(NbBl-2) = sum(CBFM_Blocks_p(NbBl-1:NbBl)%Nbc_b);
    ElseIf (NbBl .eq. 3) then ! 3 blocks 
        NbBlksL2 = 2;
        Allocate(Nbc_Level2(2))
        If (CBFM_Blocks_p(1)%Nbc_b .lt. CBFM_Blocks_p(3)%Nbc_b) Then 
            Nbc_Level2(1) = sum(CBFM_Blocks_p(1:2)%Nbc_b); 
            Nbc_Level2(2) = CBFM_Blocks_p(3)%Nbc_b
        Else
            Nbc_Level2(1) = CBFM_Blocks_p(1)%Nbc_b 
            Nbc_Level2(2) = sum(CBFM_Blocks_p(2:3)%Nbc_b); 
        EndIf            
    Else     !! Il faut trouver une solution ici au cas ou la division suivant dir0 donne un seul bloc !
        NbBlksL2 = 2;
        Allocate(Nbc_Level2(2)) ! 2 blocks
        Nbc_Level2(1:2) = CBFM_Blocks_p(1:2)%Nbc_b             
    EndIf 
    Else
        NbBlksL2 = 1;
        Allocate(Nbc_Level2(1)) ! 2 blocks
        Nbc_Level2(1) = CBFM_Blocks_p(1)%Nbc_b        
    EndIf       
    !! *****************************************************************************        
    !! *****************************************************************************
      
    if (Type_Par .eq. 3) then 
	NBlocks = NbBl
	Allocate(CBFM_Blocks(Nblocks));
	NbBl = 0; 
	Do ii= 1,NBlocks 
           if (CBFM_Blocks_p(ii)%Nbc_b .ne. 0) then 
               NbBl = NbBl + 1;
	       CBFM_Blocks(NbBl) = CBFM_Blocks_p(ii);
           endif
	EndDo
        NBlocks = NbBl
        NbintBl = 0; 

        Deallocate(CBFM_Blocks_p); Allocate(CBFM_Blocks_p(NBlocks));
        CBFM_Blocks_p(1:Nblocks) = CBFM_Blocks(1:Nblocks);
        deallocate(CBFM_Blocks); Allocate(CBFM_Blocks(Nblocks));
        CBFM_Blocks(1:Nblocks) = CBFM_Blocks_p(1:Nblocks);
        Deallocate(CBFM_Blocks_p);
    else
        
        !! *****************************************************************************
        ! *************************** Next,dividing along X and Y*********************** 
        !! *****************************************************************************
        Nbmax = minval(CBFM_Blocks_p(:)%Nbc_b) ! the ref size will be the smallest block obtained
        Nbmax = min(Nbmax,Navg_cells); !min(Nbmax,Nbcel_Blk_max);  !after divison along z
        Div_is_possible = .TRUE.;
        !switchToMean = .FALSE.;
        last_division = .False.
        dir1 = order_div(2); dir2 = order_div(3);    ! just to identify the other directions (= to eliminate the first one)
        Deallocate(order_div);
     
        Do while (Div_is_possible)! .OR. last_division) ! not sure but still possible maybe 
            !If (switchToMean) then
            !    Nbmax = sum(CBFM_Blocks_p(:)%Nbc_b)/(max(1,NbBl));
            !EndIf             
            Div_is_possible = .False. ! and if inside the current loop a division was done, Div_is_possible is set to .TRUE.           
            ! each block will be divided along its largest dimension (X or Y), only if it's a good candidate for division
            new_NbBl = NbBl;cursB = 0
            Do iiB=1,NbBl
                Nbc_block = CBFM_Blocks_p(iiB)%Nbc_b   
                cel_init = sum(CBFM_Blocks_p(1:iiB-1)%Nbc_b)+1;
                cel_final = sum(CBFM_Blocks_p(1:iiB)%Nbc_b);
              
                Allocate(cells_Block(Nbc_block));
                cells_Block(1:Nbc_block) = cells(cel_init:cel_final);
                Allocate(positions(3,Nbc_block)); 
                positions(1,1:Nbc_block) = cells_Block(1:Nbc_block)%Xc
                positions(2,1:Nbc_block) = cells_Block(1:Nbc_block)%Yc
                positions(3,1:Nbc_block) = cells_Block(1:Nbc_block)%Zc      
                           
                ! check if the current block is 'divisable'
                CanDiv = (Nbc_block .gt. Nbmax) .AND. &
                    (abs(Nbmax-Nbc_block) .gt. abs(Nbmax-Nbc_block/2)) ! Candidate to division
            
                if (CanDiv) Then                        
                    Div_is_possible = .TRUE.;
                    ! we focus on the cells of the current block that is candidate to division
                    Allocate(new_cells_Block(Nbc_block));                        
                    Call SetBlockCtrs(Type_Par,iiB,NbBl,CBFM_Blocks_p,Nbc_block,positions)
                    if (CBFM_Blocks_p(iiB)%BCubCont(3,dir1) .gt. CBFM_Blocks_p(iiB)%BCubCont(3,dir2)) then
                        dir_div = dir1; 
                    else 
                        dir_div = dir2; 
                    endif
                          
                    hbox = CBFM_Blocks_p(iiB)%BCubCont(3,dir_div) 
                    div = 2 ; !div = hbox/hBlock +1; ! the second option is not the optimal for load balancing 
                    NbBdir1or2 = div;  
                    hBldir1or2 = hbox/div    
                  
                    Allocate(cells_in_blocks(NbBdir1or2,Nbc_block+1))
                    dmin = CBFM_Blocks_p(iiB)%BCubCont(1,dir_div)-Sc/2.; 
                
                    !dstep = hBldir1or2; 
                    hbySc = hBldir1or2/Sc; 
                    dstep = nint(hbySc)*Sc;
                
                    Allocate(positions_divdir(Nbc_block))
                    positions_divdir = positions(dir_div,1:Nbc_block);
                    Call distributeCells(Nbc_block,NbBdir1or2,positions_divdir,&
                        dmin,dstep,cells_Block,cells_in_blocks)
                    ! little correction if needed (juste once to avoid complication of the code)
                    nbc1 = cells_in_blocks(1,1);
                    nbc2 = cells_in_blocks(2,1);
                
                    NbTr = 20; tr=1;
                    Step_Tr = hBldir1or2/2.;
                    ref_error = Nbc_block/2 !or !Nbmax
                    error1 = (1.*abs(nbc1- ref_error))/ref_error  
                    error2 = (1.*abs(nbc2- ref_error))/ref_error
                    prev_error = 2*error1;  
                    !perc_err = 0.1;                 
                    !Do while ((tr .le. NbTr) .AND. (((nbc1 .gt. Nbmax) .and. (nbc2 .lt. Nbmax)) .OR. ((nbc2 .gt. Nbmax) .and. (nbc1 .lt. Nbmax))))
                    !Do while ((tr .le. NbTr) .AND. (error1 .gt. perc_err) .AND. (error2 .gt. perc_err))
                    Do while ((error1 .le. prev_error) .and. (tr .le. NbTr))                    
                        prev_error = error1;
                        prev_hBldir1or2 = hBldir1or2;
                    
                        if (nbc1 .gt. nbc2) then
                            hBldir1or2 = hBldir1or2 - Step_Tr
                        else
                            hBldir1or2 = hBldir1or2 + Step_Tr
                        endif                       
                    
                        !dstep = hBldir1or2; 
                        hbySc = hBldir1or2/Sc; 
                        dstep = nint(hbySc)*Sc;
                    
                        cells_in_blocks = 0;
                      
                        Call distributeCells(Nbc_block,NbBdir1or2,positions_divdir,&
                        dmin,dstep,cells_Block,cells_in_blocks)
                        nbc1 = cells_in_blocks(1,1);
                        nbc2 = cells_in_blocks(2,1);
                    
                        error1 = (1.*abs(nbc1- ref_error))/ref_error
                        error2 = (1.*abs(nbc2- ref_error))/ref_error 
                    
                        tr = tr + 1;
                        Step_Tr = Step_Tr/2.;
                    EndDo
                    ! final distribution :
                    hbySc = prev_hBldir1or2/Sc; 
                    dstep = nint(hbySc)*Sc;                    
                    cells_in_blocks = 0;
                    Call distributeCells(Nbc_block,NbBdir1or2,positions_divdir,&
                        dmin,dstep,cells_Block,cells_in_blocks)
                    Deallocate(positions_divdir);
                
                    ! Now reorganize the cells depending on the block to which they belong 
                    curs_new_cel =1;        
                    Do numB_cel=1,NbBdir1or2
                        Nbcels = cells_in_blocks(numB_cel,1)
                        currentB = cursB+numB_cel;
                        numB_gl = curs_glb_Blk + currentB - 1;
                        CBFM_Blocks_p_tmp(currentB)%num_block = numB_gl;
                        CBFM_Blocks_p_tmp(currentB)%Nbc_b= Nbcels  
                       
                        Do num_cel=1,Nbcels
                            cel = cells_in_blocks(numB_cel,1+num_cel)   
                            New_cells_Block(curs_new_cel) = cells_Block(cel)
                            New_Cells_Block(curs_new_cel)%num_block = numB_gl;   
                            curs_new_cel = curs_new_cel + 1
                        EndDo            
                    EndDo
                  
                    New_Cells(cel_init:cel_final) = new_cells_Block(1:Nbc_block)
                    Deallocate(cells_in_blocks)   
                  
                    ! Pay attention here when updating New_Cells_p staring from New_cells_Block 
                    cursB = cursB + NbBdir1or2
                    new_NbBl = new_NbBl - 1 + NbBdir1or2;
                  
                    Deallocate(cells_Block,positions,new_cells_Block);
              
                Else
                    cursB = cursB +1; numB_gl = curs_glb_Blk + cursB -1;
                    if (cursB .le. NbBlock_s_max) then 
                        CBFM_Blocks_p_tmp(cursB) = CBFM_Blocks_p(iiB);                    
                        CBFM_Blocks_p_tmp(cursB)%num_block = numB_gl;   
                        Do num_cel=1,Nbc_block
                            cells_Block(num_cel)%num_block = numB_gl; 
                        EndDo   
                        New_Cells(cel_init:cel_final) = cells_Block(1:Nbc_block)
                        Deallocate(cells_Block,positions);
                    else
                        error_division = 1;
                        go to 10;
                    endif                
                EndIf            
            EndDo
            NbBl = new_NbBl
            Deallocate(CBFM_Blocks_p); Allocate(CBFM_Blocks_p(NbBl))
            CBFM_Blocks_p(:) = CBFM_Blocks_p_tmp(1:NbBl) 
              
            ! update the num_cell field, then put New_Cells_p in Cells_p 
            Do num_cel=1,ncp
                New_Cells(num_cel)%num_cell = num_cel                    
            EndDo
            Cells(1:ncp) = New_Cells(1:ncp);    
          
            If (last_division) Then  ! since in the last division we reduced Nbmax, Div_is_possible can become true so we need to prohibit it
                Div_is_possible = .False.
            EndIf               
            If (.NOT. Div_is_possible) Then  ! here comes the last division
                last_division = .NOT. last_division    
                Nbmax = sum(CBFM_Blocks_p(:)%Nbc_b)/(max(1,NbBl))                   
            EndIf  
            !if ((maxval(CBFM_Blocks_p(:)%Nbc_b)) .lt. (1.6*Nbmax)) then
            !    switchToMean = .TRUE.;
            !EndIf  
            !Write(*,*) 'CBFM_Blocks_p(:)%Nbc_b = ',CBFM_Blocks_p(:)%Nbc_b     
            !Write(*,*) 'Nbmax =',Nbmax              
        EndDo       
        ! Update the new_cells and CBFM_Blocks_tmp 
        CBFM_Blocks_tmp(curs_glb_Blk:curs_glb_Blk+NbBl-1) = CBFM_Blocks_p(1:NbBl); 
        curs_glb_Blk = curs_glb_Blk + NbBl; Nblocks = Nblocks + NbBl; 
        NbintBl = 0; ! not applicable to this type of division
          
        ! And finally update Cells 
        Cells = New_Cells;     
        Allocate(CBFM_Blocks(Nblocks))
        CBFM_Blocks(:) = CBFM_Blocks_tmp(1:Nblocks)
        Deallocate(New_Cells,CBFM_Blocks_tmp);
    endif
    If (MLCBFM ==1) Then
      !! For the moment, I am using a simple 2 Levels CBFM and assume that I have one single scatterer for the multilevel CBFM 
      ! test level 1 
      !NbBlksL2 = 1
      !NberLevels = 1;
      !Allocate(MLCBFM_BlDistr(1,NbBlksL2+1))
      !MLCBFM_BlDistr(1,1) = 1; 
      !MLCBFM_BlDistr(1,2) = Nblocks;
      If (Type_Par == 2) Then !.AND. (info_p_fl .ne. 's')) Then
        NberLevels = 2;
        Allocate(MLCBFM_BlDistr(1,NbBlksL2+1))
        MLCBFM_BlDistr(1,1) = NbBlksL2;
        ii=1; GB = 1;  
        Do while (ii .le. Nblocks)
            jj=0; 
            Nbc_GB = CBFM_Blocks(ii)%Nbc_b;
            Do while (Nbc_GB .lt. Nbc_Level2(GB))
                ii = ii+1; jj = jj + 1;  
                Nbc_GB = Nbc_GB + CBFM_Blocks(ii)%Nbc_b;          
            EndDo
            MLCBFM_BlDistr(1,1+GB)= jj + 1;    
            ii = ii+1;
            GB = GB + 1;
        EndDo
      Else
        ! let us start with this simple configuration 
        ! pour le moment je mets a 0 parce que je n'utilise pas la LLCBFM pour la sphere
        ! je veux developper pour commencer une repartition au 2eme niveau en 4 blocs 
        NberLevels = 2;
        NbBlksL2 = 4;
        Allocate(MLCBFM_BlDistr(1,NbBlksL2+1))
        MLCBFM_BlDistr(1,1) = NbBlksL2;
        MLCBFM_BlDistr(1,5) = 0;
      
      EndIf
    Else
      NberLevels = 1;
      NbBlksL2 = 1;
      Allocate(MLCBFM_BlDistr(1,NbBlksL2+1))
      MLCBFM_BlDistr(1,1) = NbBlksL2;
      MLCBFM_BlDistr(1,2) = 1;
    EndIf
    
    
    ! define a narrower contour for each block
    Do ii=1,Nblocks        
        Nbc_block = CBFM_Blocks(ii)%Nbc_b
        cel_init = sum(CBFM_Blocks(1:ii-1)%Nbc_b)+1;
        cel_final = sum(CBFM_Blocks(1:ii)%Nbc_b);
        ! cancel Rotation 
        Allocate(positions(3,CBFM_Blocks(ii)%Nbc_b));
        positions(1,:) = Cells(cel_init:cel_final)%Xc
        positions(2,:) = Cells(cel_init:cel_final)%Yc
        positions(3,:) = Cells(cel_init:cel_final)%Zc
          
        Call SetBlockCtrs(Type_Par,ii,NBlocks,CBFM_Blocks,Nbc_block,positions) 
        deallocate(positions);
    EndDo
        
10  return;    
    
END SUBROUTINE Division_blocks_csh

!******************************************************************************************
!******************************************************************************************
!******************************************************************************************

SUBROUTINE Division_blocks_sph(SimScatterer,Cells,Ncells_SphDomains,CBFM_Blocks,MLCBFM_BlDistr,error_division)

    ! The division into blocks depends on the type of the considered scatterer
    ! It is much simpler for the conventional shapes : Sphere, Cylinder ...
    ! The division here is for the complex geometries (from file) 
    ! STILL can be improved ...        
        
    USE Initialization
    USE common_variables
    
    Implicit NONE
    
    !IN/OUT 
    type (Scatterer), INTENT(INOUT) :: SimScatterer
    type (Cell), Dimension(Nbc), INTENT(INOUT) :: Cells
    Integer, Dimension(7), INTENT(IN) :: Ncells_SphDomains
    type (CBFM_Block), Dimension(:), allocatable, INTENT(OUT) :: CBFM_Blocks
    Integer, Dimension(:,:), allocatable, INTENT(OUT):: MLCBFM_BlDistr
    Integer, INTENT(INOUT) :: error_division
    
    ! local 
    Integer :: ii,jj,cc,d,N,N1,mod,cel_init,cel_final,cheb_l
    Integer :: Type_Par,Nbc_block,Nbcels_int,Nbcels_ext,Nbcels_Dp
    Integer :: curs_new_cel, NbBl, cel, Ix,Iy,Iz,numB_cel,num_cel,Nbcels
    Integer :: new_NbBl,cursB,iiB,dir_div,NbB_dirdiv,nbc1,nbc2
    Integer :: curB_extern,NbTr,tr,ref_error,Nbmax
    Integer, dimension(:,:), allocatable :: cells_in_blocks,pos_in_lattice
    Integer, Dimension(:), allocatable :: Nbcels_x,Nbcels_y,Nbcels_z
        
    Real(kind=8) :: h1,h,ap,Sc,ap_cheb,Dp_cheb,ap_int,ap_ext,cheb_eps
    Real(kind=8) :: hBl_dirdiv,hbox,dmin,hbySc,dstep
    Real(kind=8) :: Step_Tr,error1,error2,perc_err
    Real(kind=8), dimension(:,:), allocatable :: positions
    Real(kind=8), dimension(:), allocatable :: positions_dirdiv
    Real(kind=8), dimension(:), allocatable :: start_x,start_y,start_z
        
    type (Cell), Dimension(:), allocatable :: new_Cells,cells_Block,new_cells_Block
    type (CBFM_Block), Dimension(:), allocatable:: CBFM_Blocks_in
    
    CHARACTER(200) :: file_name
    CHARACTER tmp_str,ch1,ch2
    CHARACTER(:), allocatable::info_p_fl
    
    logical :: Div_is_possible,CanDiv,last_division
    
    type (CBFM_Block), Dimension(:), allocatable:: Blocks_extern,Blocks_extern_tmp
    
    error_division = 0;
    Type_Par = SimScatterer%type_s
    info_p_fl = trim(SimScatterer%info_s);
    Sc = SimScatterer%Sc
    
    if (Type_Par == 6) Then ! these information was used to discretize the scatterer
        Read(info_p_fl,'(a,i2,a,f5.2)') ch1,cheb_l,ch2,cheb_eps;        
        ap_cheb = (1.+cheb_eps)*SimScatterer%dm/2.;
        Dp_cheb = anint(2*ap_cheb*10**Round_S)/10**Round_S;
        Nbcels_Dp = nint(Dp_cheb/Sc);        
        ! put back ap to r0 
        ap = SimScatterer%dm/2.; 
        Nbcels_ext = ceiling((ap_cheb - (0.8*ap/sqrt(3.)))/Sc); 
        Nbcels_int = Nbcels_Dp - 2*Nbcels_ext;
        ap_int = (Nbcels_int*Sc)/2.; 
        ap_ext = ap_cheb;
    elseif (Type_Par == 1) then 
        ap = SimScatterer%dm/2.; 
        Nbcels_Dp = nint(SimScatterer%dm/Sc);
        Nbcels_ext = ceiling((ap - (0.8*ap/sqrt(3.)))/Sc); 
        Nbcels_int = Nbcels_Dp - 2*Nbcels_ext;
        ap_int = (Nbcels_int*Sc)/2.; 
        ap_ext = ap;
    endif    
        
    Allocate(start_x(7),start_y(7),start_z(7));
    Allocate(Nbcels_x(7),Nbcels_y(7),Nbcels_z(7))
    start_x = [-ap_int,-ap_ext,-ap_ext,ap_int,-ap_int,-ap_int,-ap_ext]; 
    Nbcels_x = [Nbcels_int,Nbcels_Dp,Nbcels_ext,Nbcels_ext,Nbcels_int,Nbcels_int,Nbcels_Dp]; 
    start_y = [-ap_int,-ap_ext,-ap_ext,-ap_ext,-ap_ext,ap_int,-ap_ext]; 
    Nbcels_y = [Nbcels_int,Nbcels_Dp,Nbcels_Dp,Nbcels_Dp,Nbcels_ext,Nbcels_ext,Nbcels_Dp];
    start_z = [-ap_int,-ap_ext,-ap_int,-ap_int,-ap_int,-ap_int,ap_int]; 
    Nbcels_z = [Nbcels_int,Nbcels_ext,Nbcels_int,Nbcels_int,Nbcels_int,Nbcels_int,Nbcels_ext];  
    
    curs_new_cel = 1;
    Allocate(New_Cells(Nbc));
    Allocate(CBFM_Blocks_in(NbBlock_s_max))
    
    ! For this type of scatterers, we have 7 domains to divide into blocks 
    ! Start with the Cube inside the spherical scatterer 
    N1= Ncells_SphDomains(1); 
    
    if (N1 .le. Navg_cells) then 
        NbBl = 1;
        Allocate(cells_in_blocks(NbBl,1+N1)); 
        cells_in_blocks = 0; cells_in_blocks(1,1) = N1
        cells_in_blocks(1,2:1+N1) = Cells(1:N1)%num_cell; ! or cel  
        Cells(1:N1)%num_block = 1;     
    elseif ((N1/2) .le. Navg_cells) then 
        NbBl = 2; N = N1/NbBl; h = ap_int;
        Allocate(cells_in_blocks(NbBl,1+2*N)); ! the first column is to track cc per block
        cells_in_blocks = 0;
        Do cel = 1, N1
            ! position along x
            ix =ceiling((Cells(cel)%Xc+ap_int)/h);   !-(-ap_int)     
                        
            ! The block containing this cell is numB_cel
            numB_cel = ix
            ! add this information to the table cells_in_blocks
            cc = cells_in_blocks(numB_cel,1) + 2; ! first column number of cells for the block numB_cel
            cells_in_blocks(numB_cel,1) = cells_in_blocks(numB_cel,1) +1
            cells_in_blocks(numB_cel,cc) = Cells(cel)%num_cell; ! or cel  
            Cells(cel)%num_block = numB_cel; 
        EndDo
    elseif ((N1/4) .le. Navg_cells) then 
        NbBl = 4; N = N1/NbBl; h = ap_int;d=2;
        Allocate(cells_in_blocks(NbBl,1+2*N)); ! the first column is to track cc per block
        cells_in_blocks = 0;
        Do cel = 1, N1
            ! position along x
            ix =ceiling((Cells(cel)%Xc+ap_int)/h);   !-(-ap_int)     
            ! position along y
            iy=ceiling((Cells(cel)%Yc+ap_int)/h);
                        
            ! The block containing this cell is numB_cel
            numB_cel = (iy-1)*d + ix
            ! add this information to the table cells_in_blocks
            cc = cells_in_blocks(numB_cel,1) + 2; ! first column number of cells for the block numB_cel
            cells_in_blocks(numB_cel,1) = cells_in_blocks(numB_cel,1) +1
            cells_in_blocks(numB_cel,cc) = Cells(cel)%num_cell; ! or cel  
            Cells(cel)%num_block = numB_cel; 
        EndDo
    else   
        N = N1; h1 = 2*ap_int; mod = 0;
        d = ceiling((real(N1)/real(Navg_cells))**(1./3.));
        N = ceiling(N1/d**3.); h = h1/d;
        NbBl = d**3;
        if (NbBl + 6 .le. NbBlock_s_max) then 
            Allocate(cells_in_blocks(NbBl,1+2*N)); ! the first column is to track cc per block ! 2*N here in case the division doesn't give equal blocks
            cells_in_blocks = 0;
            N1= Ncells_SphDomains(1);
            Do cel = 1, N1
                ! position along x
                ix =ceiling((Cells(cel)%Xc+ap_int)/h);   !-(-ap_int)     
                ! position along y
                iy =ceiling((Cells(cel)%Yc+ap_int)/h);
                ! position along z
                iz =ceiling((Cells(cel)%Zc+ap_int)/h);
            
                ! The block containing this cell is numB_cel
                numB_cel = (iz-1)*(d**2) + (iy-1)*d + ix
                ! add this information to the table cells_in_blocks
                cc = cells_in_blocks(numB_cel,1) + 2; ! first column number of cells for the block numB_cel
                cells_in_blocks(numB_cel,1) = cells_in_blocks(numB_cel,1) +1
                cells_in_blocks(numB_cel,cc) = Cells(cel)%num_cell; ! or cel  
                Cells(cel)%num_block = numB_cel; 
            EndDo
        Else
            error_division = 1;
            go to 10;  
        endif
    endif  
    
    ! Now reorganize the cells depending on the block to which they belong 
    ! that will make the application of the CBFM and the extension of the blocks easier
    Do numB_cel=1,NbBl
        Nbcels = cells_in_blocks(numB_cel,1)
        CBFM_Blocks_in(numB_cel)%num_Block= numB_cel
        CBFM_Blocks_in(numB_cel)%Nbc_b= Nbcels
        Do num_cel=1,Nbcels
            cel = cells_in_blocks(numB_cel,1+num_cel)
            Cells(cel)%num_cell = curs_new_cel 
            New_Cells(curs_new_cel) = Cells(cel)
            curs_new_cel = curs_new_cel + 1
        EndDo            
    EndDo       
    Deallocate(cells_in_blocks);          
    
    Cells(1:N1) = New_Cells(1:N1);
    Nblocks = NbBl;
    ! Now itertively divide into blocks the 6 other surrounding domains
    NbBl = 6;
    Allocate(Blocks_extern(NbBl)); Allocate(Blocks_extern_tmp(NbBlock_s_max-NbBl));
    Blocks_extern(:)%Nbc_b = Ncells_SphDomains(2:7);
    Div_is_possible = .TRUE.;
    Do while (Div_is_possible)
        Div_is_possible = .False.            
        ! each block will be divided along its largest dimension (X or Y or Z), only if it's a good candidate for division
        new_NbBl = NbBl; 
        cursB = 0; 
        Do iiB=1,NbBl
            Nbc_block = Blocks_extern(iiB)%Nbc_b   
            cel_init = Ncells_SphDomains(1) + sum(Blocks_extern(1:iiB-1)%Nbc_b)+1;
            cel_final = Ncells_SphDomains(1) + sum(Blocks_extern(1:iiB)%Nbc_b);
              
            Allocate(cells_Block(Nbc_block));
            cells_Block(1:Nbc_block) = cells(cel_init:cel_final);
            Allocate(positions(3,Nbc_block)); 
            positions(1,1:Nbc_block) = cells_Block(1:Nbc_block)%Xc
            positions(2,1:Nbc_block) = cells_Block(1:Nbc_block)%Yc
            positions(3,1:Nbc_block) = cells_Block(1:Nbc_block)%Zc      
                           
            ! check if the current block is 'divisable'
            CanDiv = (Nbc_block .gt. Navg_cells) .AND. &
                (abs(Navg_cells-Nbc_block) .gt. abs(Navg_cells-Nbc_block/2)) .AND. & 
                ((Nblocks + NbBl + 1) .le. NbBlock_s_max) ! Candidate to division
          
            if (CanDiv) Then                        
                Div_is_possible = .TRUE.;
                ! we focus on the cells of the current block that is candidate to division
                Allocate(new_cells_Block(Nbc_block));                        
                Call SetBlockCtrs(Type_Par,iiB,NbBl,Blocks_extern,Nbc_block,positions)
                dir_div = maxloc(Blocks_extern(iiB)%BCubCont(3,:),1);                          
                hbox = Blocks_extern(iiB)%BCubCont(3,dir_div) 
                ! always divide by 2
                NbB_dirdiv = 2;  
                hBl_dirdiv = hbox/2.    
                  
                Allocate(cells_in_blocks(NbB_dirdiv,Nbc_block+1))
                dmin = Blocks_extern(iiB)%BCubCont(1,dir_div)-Sc/2.; 
                
                hbySc = hBl_dirdiv/Sc; 
                dstep = nint(hbySc)*Sc;
                  
                Allocate(positions_dirdiv(Nbc_block))
                positions_dirdiv = positions(dir_div,1:Nbc_block);
                Call distributeCells(Nbc_block,NbB_dirdiv,positions_dirdiv,&
                    dmin,dstep,cells_Block,cells_in_blocks)
                ! little correction if needed (juste once to avoid complication of the code)
                nbc1 = cells_in_blocks(1,1);
                nbc2 = cells_in_blocks(2,1);
                  
                 
                NbTr = 100; tr=1;
                Step_Tr = hBl_dirdiv/2.;
                ref_error = Nbc_block/2 
                error1 = real(abs(nbc1- ref_error))/real(ref_error)  
                error2 = real(abs(nbc2- ref_error))/real(ref_error)  
                perc_err = 0.1;                 
                
                Do while ((tr .le. NbTr) .AND. (error1 .gt. perc_err) .AND. (error2 .gt. perc_err))
                    if (nbc1 .gt. nbc2) then
                        hBl_dirdiv = hBl_dirdiv - Step_Tr
                    else
                        hBl_dirdiv = hBl_dirdiv + Step_Tr
                    endif                       
                    
                    !dstep = hBldir1or2; 
                    hbySc = hBl_dirdiv/Sc; 
                    dstep = nint(hbySc)*Sc;
                    
                    cells_in_blocks = 0;
                      
                    Call distributeCells(Nbc_block,NbB_dirdiv,positions_dirdiv,&
                    dmin,dstep,cells_Block,cells_in_blocks)
                    nbc1 = cells_in_blocks(1,1);
                    nbc2 = cells_in_blocks(2,1);
                      
                    error1 = real(abs(nbc1- ref_error))/real(ref_error)
                    error2 = real(abs(nbc2- ref_error))/real(ref_error) 
                      
                    tr = tr + 1;
                    Step_Tr = Step_Tr/2.;
                EndDo
                Deallocate(positions_dirdiv);
                ! Now reorganize the cells depending on the block to which they belong 
                curs_new_cel =1;        
                Do numB_cel=1,NbB_dirdiv
                    Nbcels = cells_in_blocks(numB_cel,1)
                    curB_extern = cursB+numB_cel;
                    Blocks_extern_tmp(curB_extern)%num_block = Nblocks+curB_extern;
                    Blocks_extern_tmp(curB_extern)%Nbc_b= Nbcels  
                       
                    Do num_cel=1,Nbcels
                        cel = cells_in_blocks(numB_cel,1+num_cel)   
                        New_cells_Block(curs_new_cel) = cells_Block(cel)
                        New_Cells_Block(curs_new_cel)%num_block = Nblocks+curB_extern;  
                        curs_new_cel = curs_new_cel + 1
                    EndDo            
                EndDo
                  
                New_Cells(cel_init:cel_final) = new_cells_Block(1:Nbc_block)
                Deallocate(cells_in_blocks)   
                cursB = cursB + NbB_dirdiv
                new_NbBl = new_NbBl - 1 + NbB_dirdiv;                  
                Deallocate(cells_Block,positions,new_cells_Block);    
            Else
                cursB = cursB +1; 
                if (cursB .le. NbBlock_s_max) then 
                    Blocks_extern_tmp(cursB) = Blocks_extern(iiB);                    
                    Blocks_extern_tmp(cursB)%num_block = Nblocks+cursB;   
                    Do num_cel=1,Nbc_block
                        cells_Block(num_cel)%num_block = Nblocks+cursB; ! Don't forget the internal blocks num_block
                    EndDo   
                    New_Cells(cel_init:cel_final) = cells_Block(1:Nbc_block)
                    Deallocate(cells_Block,positions);
                else
                    error_division = 1;
                    go to 10;
                endif                
            EndIf            
        EndDo
        NbBl = new_NbBl
        Deallocate(Blocks_extern); Allocate(Blocks_extern(NbBl))
        Blocks_extern(:) = Blocks_extern_tmp(1:NbBl) 
              
        ! update the num_cell field, then put New_Cells_p in Cells_p 
        Do num_cel=N1+1,Nbc
            New_Cells(num_cel)%num_cell = num_cel                    
        EndDo
        Cells(N1+1:Nbc) = New_Cells(N1+1:Nbc);                         
    EndDo       
    
    Allocate(CBFM_Blocks(Nblocks+NbBl));
    CBFM_Blocks(1:Nblocks) = CBFM_Blocks_in(1:Nblocks);
    CBFM_Blocks(Nblocks+1:Nblocks+NbBl) = Blocks_extern(1:NbBl);
    deallocate(CBFM_Blocks_in,Blocks_extern);
    NbintBl = Nblocks;
    Nblocks = Nblocks + NbBl;
        
    N = sum(CBFM_Blocks(1:Nblocks)%Nbc_b);
    if (N .NE. Nbc) then 
        error_division = 2;
        go to 10;  
    endif
     
     
    Type_Par = SimScatterer%type_s 
    ! define a narrower contour for each block
    Do ii=1,Nblocks 
        CBFM_Blocks(ii)%num_block = ii;
        Nbc_block = CBFM_Blocks(ii)%Nbc_b
        cel_init = sum(CBFM_Blocks(1:ii-1)%Nbc_b)+1;
        cel_final = sum(CBFM_Blocks(1:ii)%Nbc_b);
        ! cancel Rotation 
        Allocate(positions(3,CBFM_Blocks(ii)%Nbc_b));
        positions(1,:) = Cells(cel_init:cel_final)%Xc
        positions(2,:) = Cells(cel_init:cel_final)%Yc
        positions(3,:) = Cells(cel_init:cel_final)%Zc
          
        Call SetBlockCtrs(Type_Par,ii,NBlocks,CBFM_Blocks,Nbc_block,positions) 
        deallocate(positions);
    EndDo     
    
    
    ! for the moment, the MLCBFM version is not implemented here 
    NberLevels = 1;
    NbBlksL2 = 1;
    Allocate(MLCBFM_BlDistr(1,NbBlksL2+1))
    MLCBFM_BlDistr(1,1) = NbBlksL2;
    MLCBFM_BlDistr(1,2) = 1;   
        
10  return; 
    
END SUBROUTINE Division_blocks_sph

    
SUBROUTINE BlockCenter(NbcBlk,Blk,pos,BlkCent)
  
      USE Initialization
      USE f95_precision
      Implicit NONE
  
      !IN/OUT
      Integer, INTENT(IN) :: NbcBlk
      type (CBFM_Block), INTENT(IN) :: Blk
      Real(kind=8), Dimension(3,NbcBlk), INTENT(IN) :: pos
      Real(kind=8), Dimension(3), INTENT(OUT) :: BlkCent
  
      ! Parameter for the histogram 
      Integer, Parameter:: Nr =100; ! Number of ranges 
  
      ! local 
      Integer :: ii,jj,iX,iY,iZ
      Real(kind=8) :: Xmin,Xmax,Ymin,Ymax,Zmin,Zmax
      Real(kind=8), Dimension(3,Nr) :: Range
      Integer, Dimension(3,Nr+1) :: Bucket
  
  
      ! Absolute min and max
      Xmin = Blk%BCubCont(1,1); Xmax = Blk%BCubCont(2,1)
      Ymin = Blk%BCubCont(1,2); Ymax = Blk%BCubCont(2,2)
      Zmin = Blk%BCubCont(1,3); Zmax = Blk%BCubCont(2,3)
  
      ! construct the range array 
      DO jj=1,Nr
          Range(1,jj) = Xmin + jj*(Xmax-Xmin)/Nr;
          Range(2,jj) = Ymin + jj*(Ymax-Ymin)/Nr
          Range(3,jj) = Zmin + jj*(Zmax-Zmin)/Nr        
      EndDo   
  
      Bucket(1:3,1:Nr+1) = 0;
      DO ii = 1, NbcBlk                       ! for each input score
          ! X
           DO jj = 1, Nr                    ! determine the bucket
              IF (pos(1,ii) < Range(1,jj)) THEN
                 Bucket(1,jj) = Bucket(1,jj) + 1
                 EXIT
              END IF               
           END DO 
           !Y
           DO jj = 1, Nr                    ! determine the bucket
              IF (pos(2,ii) < Range(2,jj)) THEN
                 Bucket(2,jj) = Bucket(2,jj) + 1
                 EXIT
              END IF               
           END DO 
           !Z    
           DO jj = 1, Nr                    ! determine the bucket
              IF (pos(3,ii) < Range(3,jj)) THEN
                 Bucket(3,jj) = Bucket(3,jj) + 1
                 EXIT
              END IF               
           END DO         
           ! don't forget the last bucket
           IF (pos(1,ii) >= Range(1,Nr))  Bucket(1,Nr+1) = Bucket(1,Nr+1)+1
           IF (pos(2,ii) >= Range(2,Nr))  Bucket(2,Nr+1) = Bucket(2,Nr+1)+1
           IF (pos(3,ii) >= Range(3,Nr))  Bucket(3,Nr+1) = Bucket(3,Nr+1)+1
      END DO
  
      ! So the center of the block is simply the center of the range with the maximum Bucket 
      iX = minloc(Bucket(1,:),1,mask=(Bucket(1,:) .eq. maxval(Bucket(1,:))));
      iY = minloc(Bucket(2,:),1,mask=(Bucket(2,:) .eq. maxval(Bucket(2,:))));
      iZ = minloc(Bucket(3,:),1,mask=(Bucket(3,:) .eq. maxval(Bucket(3,:))));
  
      BlkCent(1) = ((Xmin + (iX-1)*(Xmax-Xmin)/Nr) + (Xmin + iX*(Xmax-Xmin)/Nr))/2.
      BlkCent(2) = ((Ymin + (iY-1)*(Ymax-Ymin)/Nr) + (Ymin + iY*(Ymax-Ymin)/Nr))/2.
      BlkCent(3) = ((Zmin + (iZ-1)*(Zmax-Zmin)/Nr) + (Zmin + iZ*(Zmax-Zmin)/Nr))/2.  
End Subroutine BlockCenter

!******************************************************************************************
!******************************************************************************************
!******************************************************************************************
  
SUBROUTINE SetBlockCtrs(typep,numB,NbBl,CBFM_Blocks,Nbc_b,positions)
  
      USE Initialization
      USE f95_precision
      Implicit NONE
      
      !IN/OUT
      Integer, INTENT(IN) :: numB,Nbc_b,NbBl,typep
      type (CBFM_Block), Dimension(NbBl), INTENT(INOUT) :: CBFM_Blocks
      Real(kind=8), Dimension(3,Nbc_b), INTENT(IN) :: positions
  
      ! local
      Integer :: jj,cel_init,cel_final 
      Real(kind=8) :: Rd_fn,Rd
      
      CBFM_Blocks(numB)%BCubCont(1,1) = minval(positions(1,:));
      CBFM_Blocks(numB)%BCubCont(2,1) = maxval(positions(1,:));
      CBFM_Blocks(numB)%BCubCont(1,2) = minval(positions(2,:));
      CBFM_Blocks(numB)%BCubCont(2,2) = maxval(positions(2,:));
      CBFM_Blocks(numB)%BCubCont(1,3) = minval(positions(3,:));
      CBFM_Blocks(numB)%BCubCont(2,3) = maxval(positions(3,:));
        
      CBFM_Blocks(numB)%BCubCont(3,1) = CBFM_Blocks(numB)%BCubCont(2,1) - CBFM_Blocks(numB)%BCubCont(1,1)
      CBFM_Blocks(numB)%BCubCont(3,2) = CBFM_Blocks(numB)%BCubCont(2,2) - CBFM_Blocks(numB)%BCubCont(1,2)
      CBFM_Blocks(numB)%BCubCont(3,3) = CBFM_Blocks(numB)%BCubCont(2,3) - CBFM_Blocks(numB)%BCubCont(1,3)
      
      
      ! starting from the cubical contour of the block of cells, we define a spherical contour as a center (Xc,Yc and Zc)
      ! and Rb (radius). For a better definition of adjascent blocks later in the code, we will check the distribution of 
      ! the cordinates of the cells to better define the center of the block (instead of a simple (c_min + c_max)/2)
      ! A REVOIR
      CBFM_Blocks(numB)%BSphCont(1) = (CBFM_Blocks(numB)%BCubCont(1,1)+CBFM_Blocks(numB)%BCubCont(2,1))/2.;
      CBFM_Blocks(numB)%BSphCont(2) = (CBFM_Blocks(numB)%BCubCont(1,2)+CBFM_Blocks(numB)%BCubCont(2,2))/2.;
      CBFM_Blocks(numB)%BSphCont(3) = (CBFM_Blocks(numB)%BCubCont(1,3)+CBFM_Blocks(numB)%BCubCont(2,3))/2.;
      
      !Call BlockCenter(CBFM_Blocks(numB)%Nbc_b,CBFM_Blocks(numB),positions,CBFM_Blocks(numB)%BSphCont(1:3));
      
      Rd_fn = 0;
      Do jj=1,Nbc_b
          Rd = sqrt((positions(1,jj)-CBFM_Blocks(numB)%BSphCont(1))**2.+(positions(2,jj)-CBFM_Blocks(numB)%BSphCont(2))**2 &
              + (positions(3,jj)-CBFM_Blocks(numB)%BSphCont(3))**2)
          Rd_fn = max(Rd_fn,Rd);
      EndDo  
      CBFM_Blocks(numB)%BSphCont(4)=Rd_fn;
  
      
END SUBROUTINE SetBlockCtrs

!******************************************************************************************
!******************************************************************************************
!******************************************************************************************
  
SUBROUTINE distributeCells(Nbc_td,NbBl,positions_td,dim_min,dim_step,Cells_td,cells_in_blocks)
  
      USE Initialization
      USE f95_precision
      Implicit NONE
  
      ! IN/OUT
      INTEGER, INTENT(IN) :: Nbc_td,NbBl
      Real(kind=8), INTENT(IN) :: dim_min,dim_step
      Real(kind=8), Dimension(Nbc_td), INTENT(IN) :: positions_td
      type (Cell), Dimension(Nbc_td), INTENT(IN):: Cells_td
      Integer, dimension(NbBl,Nbc_td+1), INTENT(OUT) :: cells_in_blocks
      ! local
      
      Integer :: cel,i_dir,numB_cel,numB_cel_gl,cc
      Real(kind=8) :: cord     
  
      ! Now that the blocks are defined and ordred, we scan the cells 
      ! belonging to the current scatterer to distribute them 
      ! throughout the CBFM blocks (along Z for now)
      cells_in_blocks = 0;
      ! each row of cells_in_blocks contain in the first column the number of cells
      ! comprising this block, than the indexes of these cells in the following columns
      Do cel = 1, Nbc_td
          !if (cel== 956) then
          !    Write(*,*) 'cel 956';
          !endif
          ! position along z
          cord = positions_td(cel);
          !i_dir=0;
          !Do while (cord .gt. i_dir*dim_step+dim_min)
          !    i_dir = i_dir +1 ;                
          !EndDo  
          i_dir = ceiling((cord-dim_min)/dim_step);
          ! The block containing this cell is numB_cel
          numB_cel = min(max(i_dir,1),NbBl) ! to avoid problems with the cells of the boundaries
          
          ! add this information to the table cells_in_blocks
          cc = cells_in_blocks(numB_cel,1) + 2; ! first column number of cells for the block numB_cel
          cells_in_blocks(numB_cel,1) = cells_in_blocks(numB_cel,1) +1
          cells_in_blocks(numB_cel,cc) = cel;       
      EndDo
END SUBROUTINE distributeCells