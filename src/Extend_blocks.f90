SUBROUTINE Extend_blocks(SimScatterer,Cells,CBFM_Blocks,CBFM_Blocks_Ext)

    USE Initialization
    USE common_variables
    USE MPI

    Implicit NONE

    !IN/OUT 
    type(Scatterer), INTENT(IN) :: SimScatterer
    type(Cell), Dimension(Nbc), INTENT(IN) :: Cells
    type(CBFM_Block), Dimension(Nblocks), INTENT(INOUT) :: CBFM_Blocks
    Integer, Dimension(:,:), allocatable, INTENT(OUT) :: CBFM_Blocks_Ext
    
    !! LOCAL
    Integer :: ii,jj,Type_Par,BB,Nbext_max,Block_init,Block_final
    Integer :: BBp,NbBadj,Upd_NbBadj,numBlock, cel_init,cel_final
    Integer :: NbCelExt,NbCelExt_Badj,pp,ncb,cel_in_Badj
    Integer, Dimension(:,:), allocatable :: Extensions
    Integer, Dimension(:), allocatable :: Badj,Upd_Badj
    character(200) :: file_name
    Real(kind=8) :: ExtLength, Sc
    Real(kind=8) :: x_min1,x_max1,y_min1,y_max1,z_min1,z_max1
    Real(kind=8) :: x_min2,x_max2,y_min2,y_max2,z_min2,z_max2
    Real(kind=8) :: X_extZo_min,X_extZo_max,Y_extZo_min,Y_extZo_max
    Real(kind=8) :: Z_extZo_min,Z_extZo_max
    Real(kind=8) :: xc_b1,yc_b1,zc_b1,xc_b2,yc_b2,zc_b2,dis_12,thr_dis,lam_p
    Real(kind=8) :: X_ext_min,X_ext_max,Y_ext_min,Y_ext_max
    Real(kind=8) :: Z_ext_min,Z_ext_max,Xc,Yc,Zc
    Real(kind=8), dimension(:,:), allocatable :: positions
    
    
    
    ! Initialize the table Extensions
    Nbext_max = Fact_Nbext_max*maxval(CBFM_Blocks(:)%Nbc_b);
    Allocate(Extensions(Nblocks,Nbext_max));
    Extensions = 0;
            
    ! Scan the blocks and determine the cells present in the peripheral of the block
    pp = 1
    Type_Par = SimScatterer%type_s  
    
    Block_init = 1;
    Block_final = Nblocks;
    Sc = SimScatterer%Sc;
    
    if ((Block_init .ne. Block_final) .and. (Nc_extended .ne. 0)) Then         
        Do BB=Block_init,Block_final
            ExtLength = Nc_extended*Sc;
                  
            xc_b1 = CBFM_Blocks(BB)%BSphCont(1);
            yc_b1 = CBFM_Blocks(BB)%BSphCont(2);
            zc_b1 = CBFM_Blocks(BB)%BSphCont(3);
            ! determine the adjascent blocks based on the limit coordinate 
            ! of each block belonging to this scatterer
            Allocate(Badj(Nblocks));
            NbBadj = 0; Badj =0
          
            ! find the adjascent blocks 
            ! scan before BB
            Do BBp=Block_init,BB-1
                xc_b2 = CBFM_Blocks(BBp)%BSphCont(1);
                yc_b2 = CBFM_Blocks(BBp)%BSphCont(2);
                zc_b2 = CBFM_Blocks(BBp)%BSphCont(3);
                      
                dis_12 = sqrt((xc_b2-xc_b1)**2.+(yc_b2-yc_b1)**2.+(zc_b2-zc_b1)**2.)
                thr_dis = CBFM_Blocks(BB)%BSphCont(4)+CBFM_Blocks(BBp)%BSphCont(4);
                If (dis_12 .le. thr_dis) Then
                    NbBadj = NbBadj + 1;
                    Badj(NbBadj)=BBp;                    
                EndIf                              
            EndDo 
            ! scan after BB  
            Do BBp=BB+1,Block_final
                xc_b2 = CBFM_Blocks(BBp)%BSphCont(1);
                yc_b2 = CBFM_Blocks(BBp)%BSphCont(2);
                zc_b2 = CBFM_Blocks(BBp)%BSphCont(3);
                      
                dis_12 = sqrt((xc_b2-xc_b1)**2.+(yc_b2-yc_b1)**2.+(zc_b2-zc_b1)**2.)
                !Waiting the enhancment ... 
                thr_dis = (CBFM_Blocks(BB)%BSphCont(4)+CBFM_Blocks(BBp)%BSphCont(4));
                If (dis_12 .le. thr_dis) Then
                    NbBadj = NbBadj + 1;
                    Badj(NbBadj)=BBp;                    
                EndIf                        
            EndDo
          
            ! Now that we know the adjascent blocksto BB, let's determine the cells of these blocks
            ! belonging to the extension zone of the block BB 
            ! first Based on the contour of the block, let us define the extended zone 
            X_ext_min= CBFM_Blocks(BB)%BCubCont(1,1)-ExtLength-Sc/2.; 
            X_ext_max= CBFM_Blocks(BB)%BCubCont(2,1)+ExtLength+Sc/2.;
            Y_ext_min= CBFM_Blocks(BB)%BCubCont(1,2)-ExtLength-Sc/2.; 
            Y_ext_max= CBFM_Blocks(BB)%BCubCont(2,2)+ExtLength+Sc/2.; 
            Z_ext_min= CBFM_Blocks(BB)%BCubCont(1,3)-ExtLength-Sc/2.; 
            Z_ext_max= CBFM_Blocks(BB)%BCubCont(2,3)+ExtLength+Sc/2.; 
                  
          
            NbCelExt = 0;
            Upd_NbBadj = 0;
            Allocate(Upd_Badj(NbBadj)); ! update NbBadj and Badj depending on the number of cells 
            ! extending the block belonging to the current blocks considered as adjascent
            Do numBlock=1,NbBadj
                BBp = Badj(numBlock)
                ncb = CBFM_Blocks(BBp)%Nbc_b
                cel_init = sum(CBFM_Blocks(1:BBp-1)%Nbc_b)+1;
                cel_final = sum(CBFM_Blocks(1:BBp)%Nbc_b);
                      
                ! Need to cancel the rotation of the cells first (get them back to the vertical position)
                Allocate(Positions(3,ncb))
                positions(1,1:ncb) = Cells(cel_init:cel_final)%Xc
                positions(2,1:ncb) = Cells(cel_init:cel_final)%Yc
                positions(3,1:ncb) = Cells(cel_init:cel_final)%Zc
                                    
                NbCelExt_Badj = 0;    ! used to check if this block is actually adjascent
                cel_in_Badj= 1;
                Do ii=cel_init,cel_final                            
                    Xc = positions(1,cel_in_Badj); 
                    Yc = positions(2,cel_in_Badj);
                    Zc = positions(3,cel_in_Badj);
                    if ((X_ext_min .le. Xc) .AND. (Xc .le. X_ext_max) .AND. &
                        (Y_ext_min .le. Yc) .AND. (Yc .le. Y_ext_max) .AND. & 
                        (Z_ext_min .le. Zc) .AND. (Zc .le. Z_ext_max)) then
                      
                        NbCelExt = NbCelExt + 1;
                        NbCelExt_Badj = NbCelExt_Badj + 1; 
                        Extensions(BB,NbCelExt) = ii;  !or simply ii
                    EndIf 
                    cel_in_Badj = cel_in_Badj + 1;                
                EndDo 
                ! here update NbBadj and Badj depending on NbCelExt_Badj
                if (NbCelExt_Badj .ne. 0) Then
                    Upd_NbBadj = Upd_NbBadj + 1;
                    Upd_Badj(Upd_NbBadj) = BBp                        
                EndIf                    
                Deallocate(Positions)
            EndDo                     
                  
            CBFM_Blocks(BB)%Nbc_ext = NbCelExt
            ! Attention : we don't forget to add the obtained information about NbBadj and Badj to each CBFM block
            CBFM_Blocks(BB)%NbBadj = Upd_NbBadj;
            Allocate(CBFM_Blocks(BB)%Badj(Upd_NbBadj)); CBFM_Blocks(BB)%Badj = Upd_Badj; 
            Deallocate(Badj,Upd_Badj)        
        EndDo
    Else
        CBFM_Blocks(Block_init:Block_final)%Nbc_ext = 0;
        Extensions(Block_init:Block_final,:) = 0;
    EndIf 
    

    Nbc_ext = max(maxval(CBFM_Blocks(:)%Nbc_ext),1);    
    Allocate(CBFM_Blocks_Ext(Nblocks,Nbc_ext));
    CBFM_Blocks_Ext(1:Nblocks,1:Nbc_ext)=Extensions(1:Nblocks,1:Nbc_ext)         

END SUBROUTINE Extend_blocks 
