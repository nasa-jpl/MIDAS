RECURSIVE SUBROUTINE SMWFB_algorithm(lev_SMW,Cells,CellsBlock,n,nrhs,Ein,Eout)
     
    USE Initialization
    USE common_variables
    USE iso_fortran_env
    USE lapack95
    
    ! "Implicit Statement" 
    IMPLICIT NONE

    !! IN/OUT ******************************************************************
    type (Cell), Dimension(Nbc), INTENT(IN) :: Cells
    type (Cell), Dimension(n), INTENT(IN) :: CellsBlock
    INTEGER, INTENT(INOUT) :: lev_SMW
    INTEGER, INTENT(IN) ::n,nrhs
    !COMPLEX(real64), Dimension(n,n), INTENT(INOUT):: Z
    COMPLEX(real64), Dimension(3*n,nrhs), INTENT(INOUT):: Ein
    COMPLEX(real64), Dimension(3*n,nrhs), INTENT(OUT):: Eout
       
    ! Local
    Integer, Parameter :: Nber_lev_SMW = 2;
    Integer :: n1,n2,curs_lig,curs_col,nout1,nout2,curr_lev_SMW
    Integer :: ii,jj,curs_ligI,curs_colJ
    Integer :: n1_in, n2_in,spr_perc,spr_size,nnz
    Integer :: nrhs_tot,Nb_it_max_smwf,mtype,iparm3
    Integer :: div_dir,m,itr,n_iter
    Real(kind=8) :: RCOND,Epsilon,dmin,dstep,Sc,norm_Zspr
    Real(kind=8), dimension(3) :: blk_dims
    Real(kind=8), dimension(:), allocatable :: Cells_pos
    COMPLEX(real64) :: alpha, beta
    type (Cell), Dimension(:),allocatable :: CellsBlock1
    type (Cell), Dimension(:),allocatable :: CellsBlock2
    Integer, Dimension(:), allocatable :: row_sprZ,col_sprZ
    Integer, dimension(:,:), allocatable :: cells_in_blocks
    COMPLEX(real64), Dimension(:,:), allocatable :: U_Nmax,V_Nmax,U12,V12,U21,V21
    COMPLEX(real64), Dimension(:,:), allocatable :: Zf,E1,E2,Etot,Etoto,E1o,E2o,U12o,U21o
    COMPLEX(real64), Dimension(:,:), allocatable :: U1,U2,U1o,U2o
    COMPLEX(real64), Dimension(:,:), allocatable :: P12,P21,X1,X2,Y1,Y2,A,B,I,Mmult
    COMPLEX(real64), Dimension(:),allocatable :: Zpatch_e_spr_tmp
    COMPLEX(real64), Dimension(:),POINTER :: Zpatch_e_spr
    
    alpha = 1.;beta = 0.
    ! mtype for PRADISO
    If (homogs == 1) Then 
      mtype = 6 
    Else
      mtype = 13 
    EndIf   
    ! here properties of the ACA specific to the SMWF algorithm 
    Nb_it_max_smwf = 300;
    Epsilon = 1e-3;
    !**********************************************************
    
    !1)    
    If (lev_SMW == (Nber_lev_SMW+1)) Then
        if (SR==1) Then
            spr_perc = 50; ! let's say that we keep spr_perc % of the initial matrix Zii 
                           ! Remember also that you're keeping only the upper part of the matrix
            spr_size = nint((spr_perc*9.*n**2.)/100.)
            fct_SR = 1e3;
            Call SR_Green_s_tr_partial_FN(n,CellsBlock,fct_SR,nnz,norm_Zspr);
            
            Allocate(Zpatch_e_spr(nnz),row_sprZ(3*n+1),col_sprZ(nnz));
            Call SR_Green_s_tr_partial(n,CellsBlock,fct_SR,nnz,Zpatch_e_spr,row_sprZ,col_sprZ)
            
            Write(*,'(a,i5)',advance='no') 'solve Zf of block of size n =',n
            Write(*,'(a,f5.2,a)') ': nnz = ',100.*nnz/(9*n**2.),' % of Zii' 
          
            iparm3 = 0;
            Call pardiso_solver(3*n,nrhs,nnz,mtype,iparm3,row_sprZ,col_sprZ,Zpatch_e_spr,Ein,Eout)
             
            Deallocate(Zpatch_e_spr,row_sprZ,col_sprZ);
        else
          Allocate(Zf(3*n,3*n))
          Call Green_s_tr_partial(n,CellsBlock,n,CellsBlock,Zf) 
          !Write(*,'(a,i5)') 'solve Zf of block of size n =',n;
          call gesvx(Zf,Ein,Eout,RCOND=RCOND);     
        endif
                  
    Else
        lev_SMW = lev_SMW + 1;
        
        ! Efficient/convenient division into 2 blocks -----------
        blk_dims(1) = maxval(CellsBlock(1:n)%Xc)-minval(CellsBlock(1:n)%Xc);
        blk_dims(2) = maxval(CellsBlock(1:n)%Yc)-minval(CellsBlock(1:n)%Yc);
        blk_dims(3) = maxval(CellsBlock(1:n)%Zc)-minval(CellsBlock(1:n)%Zc);
        
        div_dir = maxloc(blk_dims,1) ;
        
        Allocate(Cells_pos(n))
        if (div_dir==1) then 
            Cells_pos(1:n) = CellsBlock(1:n)%Xc;
        elseif (div_dir==2) then 
            Cells_pos(1:n) = CellsBlock(1:n)%Yc;
        else
            Cells_pos(1:n) = CellsBlock(1:n)%Zc;
        endif
        m = 2;Sc = CellsBlock(1)%Sc
        dmin = minval(Cells_pos);
        ii = (blk_dims(div_dir)/m)/Sc;
        dstep = ii*Sc;
        
        Allocate(cells_in_blocks(m,n))
        
        Call distributeCells(n,m,Cells_pos,dmin,dstep,CellsBlock,cells_in_blocks)        
        n1 = cells_in_blocks(1,1);
        n2 = cells_in_blocks(2,1);
        n_iter = 10;itr=1;
        Do while ((itr .le. n_iter) .and. (abs(n1-n/2)/n .gt. 0.2))
            if (n1 .gt. n2) then
                dstep = dstep - dstep/n_iter; ! such as dstep/ntr itsel decreases from an iteration to another                
            else
                dstep = dstep + dstep/n_iter;                
            endif
            cells_in_blocks =0;
            Call distributeCells(n,m,Cells_pos,dmin,dstep,CellsBlock,cells_in_blocks)  
            n1 = cells_in_blocks(1,1);
            n2 = cells_in_blocks(2,1);
            itr = itr + 1;
        EndDo     
        
        !n1 = n/2
        !n2 = n-n1;
        Allocate(CellsBlock1(n1),CellsBlock2(n2));
        CellsBlock1(1:n1) = CellsBlock(cells_in_blocks(1,2:1+n1));
        CellsBlock2(1:n2) = CellsBlock(cells_in_blocks(2,2:1+n2));
        
        !CellsBlock1(1:n1) = CellsBlock(1:n1);
        !CellsBlock2(1:n2) = CellsBlock(n1+1:n);
        
        n1_in = n1; n1 = 3*n1
        n2_in = n2; n2 = 3*n2
        !2)                
        ! ACA applied to Z12 of size n1 x n2 -- > U12 and V12 
        Allocate(U_Nmax(n1,Nb_it_max_smwf), V_Nmax(Nb_it_max_smwf,n2))
        Call Calcul_Matrix_ACA_SMW(n1_in,n2_in,CellsBlock1,CellsBlock2,Nb_it_max_smwf,Epsilon,nout1,U_Nmax,V_Nmax) 
        nout1 = min(nout1,Nb_it_max_smwf);                
        Allocate(U12(n1,nout1), V12(nout1,n2))
        U12 = U_Nmax(1:n1,1:nout1); deallocate(U_Nmax)
        V12 = V_Nmax(1:nout1,1:n2); deallocate(V_Nmax)       
         
        
        !3)
        ! ACA applied to Z21 of size n2 x n1 
        Allocate(U_Nmax(n2,Nb_it_max_smwf), V_Nmax(Nb_it_max_smwf,n1))
        Call Calcul_Matrix_ACA_SMW(n2_in,n1_in,CellsBlock2,CellsBlock1,Nb_it_max_smwf,Epsilon,nout2,U_Nmax,V_Nmax)              
        Allocate(U21(n2,nout2), V21(nout2,n1))
        U21 = U_Nmax(1:n2,1:nout2); deallocate(U_Nmax)
        V21 = V_Nmax(1:nout2,1:n1); deallocate(V_Nmax) 
        
               
        curr_lev_SMW = lev_SMW; ! to avoid any impact of 4) on lev_SMW of 5)
        !4) handle Z11
        ! get Z11 and E1 from Z,Einc
        !Allocate(E1(n1,nrhs),E1o(n1,nrhs),U12o(n1,nout1))
        !E1(1:n1,1:nrhs) = Ein(1:n1,1:nrhs)
        nrhs_tot = nrhs + nout1;
        Allocate(Etot(n1,nrhs_tot),Etoto(n1,nrhs_tot))
        Etot(1:n1,1:nrhs) = Ein(1:n1,1:nrhs)
        Etot(1:n1,nrhs+1:nrhs_tot) = U12(1:n1,1:nout1);
        
        !! E1=Z11-1*E1
        !call SMWFB_algorithm_1(lev_SMW,Cells,CellsBlock1,n1_in,nrhs,E1,E1o);        
        !! U12=Z11-1*U12
        !call SMWFB_algorithm_1(curr_lev_SMW,Cells,CellsBlock1,n1_in,nout1,U12,U12o)
        ! update E1 and U12
        
        call SMWFB_algorithm(lev_SMW,Cells,CellsBlock1,n1_in,nrhs_tot,Etot,Etoto); 
        
        Allocate(E1(n1,nrhs));
        !E1 = E1o;
        !U12 = U12o       
        E1 = Etoto(1:n1,1:nrhs);
        U12 = Etoto(1:n1,nrhs+1:nrhs_tot);
        deallocate(Etot,Etoto);
             
        !5) handle Z22
        !Allocate(E2(n2,nrhs),E2o(n2,nrhs),U21o(n2,nout2))
        !E2(1:n2,1:nrhs) = Ein(n1+1:3*n,1:nrhs)
        nrhs_tot = nrhs + nout2;
        Allocate(Etot(n2,nrhs_tot),Etoto(n2,nrhs_tot))
        Etot(1:n2,1:nrhs) = Ein(n1+1:3*n,1:nrhs)
        Etot(1:n2,nrhs+1:nrhs_tot) = U21(1:n2,1:nout2)
        
        !!E2=Z22-1*E2
        !call SMWFB_algorithm_1(curr_lev_SMW,Cells,CellsBlock2,n2_in,nrhs,E2,E2o)
        !!U21=Z22-1*U21
        !call SMWFB_algorithm_1(curr_lev_SMW,Cells,CellsBlock2,n2_in,nout2,U21,U21o)
        
        call SMWFB_algorithm(curr_lev_SMW,Cells,CellsBlock2,n2_in,nrhs_tot,Etot,Etoto); 
        Allocate(E2(n2,nrhs));
        E2 = Etoto(1:n2,1:nrhs); 
        U21=Etot(1:n2,nrhs+1:nrhs_tot);
        deallocate(Etot,Etoto);
        
                
        !6) & 7)
        Allocate(P12(nout2,nout1),P21(nout1,nout2));
        Allocate(Y1(nout2,nrhs),Y2(nout1,nrhs))
        
        !P12= V21 x U12
        CALL ZGEMM('N','N',nout2,nout1,n1,alpha,V21,nout2,U12,n1,beta,P12,nout2)            
        !Y1 = V21 x E1
        CALL ZGEMM('N','N',nout2,nrhs,n1,alpha,V21,nout2,E1,n1,beta,Y1,nout2) 
        !P21= V12 x U21
        CALL ZGEMM('N','N',nout1,nout2,n2,alpha,V12,nout1,U21,n2,beta,P21,nout1)    
        !Y2 = V12 x E2
        CALL ZGEMM('N','N',nout1,nrhs,n2,alpha,V12,nout1,E2,n2,beta,Y2,nout1) 
        
        ! 8)
        Allocate(A(nout1,nout1),I(nout1,nout1))
        Allocate(B(nout1,nrhs));
        Allocate(X1(nout2,nrhs),X2(nout1,nrhs))
        
        I = 0;
        Do ii=1,nout1
            I(ii,ii)=1; 
        EndDo
        Allocate(Mmult(nout1,nout1))
        CALL ZGEMM('N','N',nout1,nout1,nout2,alpha,P21,nout1,P12,nout2,beta,Mmult,nout1) 
        A = I - Mmult; 
        deallocate(Mmult);
        
        Allocate(Mmult(nout1,nrhs))
        CALL ZGEMM('N','N',nout1,nrhs,nout2,alpha,P21,nout1,Y1,nout2,beta,Mmult,nout1) 
        B = Y2 - Mmult; 
        deallocate(Mmult)
        call gesvx(A,B,X2,RCOND=RCOND) ! X2
        
                
        Allocate(Mmult(nout2,nrhs))
        CALL ZGEMM('N','N',nout2,nrhs,nout1,alpha,P12,nout2,X2,nout1,beta,Mmult,nout2) 
        X1 = Y1 - Mmult; deallocate(Mmult) ! X1
                
        !9)
        Allocate(Mmult(n1,nrhs));
        !Mmult = matmul(U12,X2);
        CALL ZGEMM('N','N',n1,nrhs,nout1,alpha,U12,n1,X2,nout1,beta,Mmult,n1)         
        Eout(1:n1,1:nrhs) = E1-Mmult;
        Deallocate(Mmult); 
        Allocate(Mmult(n2,nrhs));
        CALL ZGEMM('N','N',n2,nrhs,nout2,alpha,U21,n2,X1,nout2,beta,Mmult,n2)   
        Eout(n1+1:3*n,1:nrhs) = E2-Mmult
        Deallocate(Mmult);       
        
    EndIf
    
END SUBROUTINE SMWFB_algorithm