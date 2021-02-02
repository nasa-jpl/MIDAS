!! SUBROUTINES : 
!! - Calcul_MatrixZij_ACA
!! - Calcul_GenMatrix_ACA
!! - Calcul_Matrix_ACA_SMW

SUBROUTINE Calcul_MatrixZij_ACA(Cells,I_Z,J_Z,curs_lig_Iz,curs_col_Jz,NbreLigMatZ,NbreColMatZ,nb_iter_out,Matrix_U,Matrix_V)  
     
    USE Initialization
    USE common_variables
    USE iso_fortran_env

    ! "Implicit Statement" 
    IMPLICIT NONE

    !! IN/OUT ******************************************************************
    type (Cell), Dimension(Nbc), INTENT(IN) :: Cells
    INTEGER, INTENT(IN) :: I_Z,J_Z,curs_lig_Iz,curs_col_Jz,NbreLigMatZ,NbreColMatZ
    INTEGER, INTENT(OUT) :: nb_iter_out
    COMPLEX(real64), Dimension(NbreLigMatZ,Nb_it_max), INTENT(OUT):: Matrix_U
    COMPLEX(real64), Dimension(Nb_it_max,NbreColMatZ), INTENT(OUT):: Matrix_V
        
    ! Local
    Integer :: curs_lig_green,curs_col_green,col,ii,jj,k,Ik,Jk,nj,l,ccol,llig
    Integer, Dimension(:), allocatable ::Indexes_Row, Indexes_Column
    Integer, Dimension(2) :: loc_max_i_j    
    Real(kind=8) :: Epsilon,norm_prec_Z,norm_Uk,norm_Vk,val_2,Res
    real (kind = 4), external :: ZLANGE
    Real(kind=8), Dimension(:), allocatable :: NormeF_Zk
    Real(kind=8), Dimension(1,1) ::sum_norm_UV
    Logical :: Continuons
    DOUBLE PRECISION, Dimension(max(NbreLigMatZ,NbreColMatZ)) :: WORK
    COMPLEX(real64) :: pivot
    COMPLEX(real64), Dimension(:,:), allocatable :: ProdUlVl_row, ProdUlVl_col,vect_R_tmp,Unj,Uk,Vnj,Vk
    COMPLEX(real64), Dimension(:,:), allocatable :: Appro_R, Z_rowI, Z_columnJ
    
    Epsilon = Epsilon_ACA
                
    ! Indexes_Row et Indexes_Column are the arrays containing orderly selected row and column indexes of the matrix Zmn
    Allocate(Indexes_Row(Nb_it_max+1)); Allocate(Indexes_Column(Nb_it_max+1))
    Indexes_Row = 0;Indexes_Column = 0;
    !! line or column from Zij depending on the prgress in the algorithm
    Allocate(Z_rowI(1,NbreColMatZ)); Allocate(Z_columnJ(NbreLigMatZ,1))
    !! Approximate Z matrix and Appproximate Error matrix
    Allocate(Appro_R(NbreLigMatZ,NbreColMatZ)); Appro_R = 0;
    Allocate(NormeF_Zk(Nb_it_max+1))
    Allocate(ProdUlVl_row(1,NbreColMatZ));Allocate(ProdUlVl_col(NbreLigMatZ,1))
    
    !! Initialisation 
    Matrix_U=0; Matrix_V=0
    Indexes_Row(1) = 1  !!I1
    !! Pour initialiser la premiere ligne de R_appro, j'appelle Green_s_tr_partielle et Green_ns_pi_tr_partielle pour calculer Zij(1,:); 1 etant la premiere ligne de la matrice Zij qui concerne le bloc i,j de la matrice totale
    !Z_rowI = Matrice2(Indexes_Row(1):Indexes_Row(1),1:NbreColMatZ)
    curs_lig_Green = curs_lig_Iz+Indexes_Row(1)-1
    Call computeBlockRow(Cells,curs_lig_Green,curs_col_Jz,NbreColMatZ,Z_rowI)  
    Appro_R(1:1,1:NbreColMatZ) = Z_rowI
    Deallocate(Z_rowI)
    loc_max_i_j = maxloc(abs(Appro_R(Indexes_Row(1):Indexes_Row(1),1:NbreColMatZ))) 
    Indexes_Column(1) = loc_max_i_j(2)   !!J1
    !! Initialisation de V
    pivot = Appro_R(Indexes_Row(1),Indexes_Column(1))
    Do col=1,NbreColMatZ
        Matrix_V(1,col) = Appro_R(Indexes_Row(1),col)/pivot
    EndDo 
                
    !Z_columnJ = Matrice2(1:NbreLigMatZ,Indexes_Column(1):Indexes_Column(1))
    curs_col_Green = curs_col_Jz+Indexes_Column(1)-1
    Call computeBlockCol(Cells,curs_lig_Iz,NbreLigMatZ,curs_col_Green,Z_columnJ)                                    
    Appro_R(1:NbreLigMatZ,Indexes_Column(1):Indexes_Column(1)) = Z_columnJ
    Deallocate(Z_columnJ)
    !! Initialisation de U
    Matrix_U(1:NbreLigMatZ,1:1) = Appro_R(1:NbreLigMatZ,Indexes_Column(1):Indexes_Column(1))
                
        
    !! Frobenius Norm of Z1
    !norm_prec_Z = 0.          
    !! Frobenius_Norm(Matrix_U(1:NbreLigMatZ,1:1))
    !norm_Uk = ZLANGE('f',NbreLigMatZ,1,Matrix_U(1:NbreLigMatZ,1:1),NbreLigMatZ,work)
    !! Frobenius_Norm(Matrix_V(1:1,1:NbreColMatZ))
    !norm_Vk = ZLANGE('f',1,NbreColMatZ,Matrix_V(1:1,1:NbreColMatZ),1,work)   
    
    norm_prec_Z = 0.
    ! Frobenius_Norm(Matrix_U(1:NbreLigMatZ,1:1))
    val_2 = 0.
    Do ii =1,NbreLigMatZ
        val_2 = val_2 + (abs(Matrix_U(ii,1)))**2
    EndDo
    norm_Uk = sqrt(val_2)
    
    ! Frobenius_Norm(Matrix_V(1:1,1:NbreColMatZ))
    val_2 = 0.
    Do jj=1,NbreColMatZ
        val_2 = val_2 + (abs(Matrix_V(1,jj)))**2
    EndDo
    norm_Vk = sqrt(val_2)
    
    NormeF_Zk(1) = sqrt(norm_prec_Z**2 + norm_Uk**2*norm_Vk**2)
    
    !! attention max(R(i,Jk)), i!=I1(=1)
    Allocate(vect_R_tmp(NbreLigMatZ,1))
    vect_R_tmp(1:NbreLigMatZ,1:1) = Appro_R(1:NbreLigMatZ,Indexes_Column(1):Indexes_Column(1))
    vect_R_tmp(Indexes_Row(1),1) = 0
    loc_max_i_j = maxloc(abs(vect_R_tmp(1:NbreLigMatZ,1:1))) 
    Indexes_Row(2) = loc_max_i_j(1)
    Deallocate(vect_R_tmp)
                
    k=2
    continuons = .TRUE.
    Do while (continuons .AND. k <= Nb_it_max)
        !!(1)
        Ik = Indexes_Row(k)
        ProdUlVl_row = 0.
        Do l=1,k-1
            ProdUlVl_row(1,1:NbreColMatZ) = ProdUlVl_row(1,1:NbreColMatZ) + Matrix_U(Ik,l)*Matrix_V(l,1:NbreColMatZ)
        EndDo
        !Appro_R(Ik,1:NbreColMatZ) = Matrice2(Ik,1:NbreColMatZ) - ProdUlVl_row(1,1:NbreColMatZ)
        Allocate(Z_rowI(1,NbreColMatZ))
        curs_lig_Green = curs_lig_Iz+Ik-1
        Call computeBlockRow(Cells,curs_lig_Green,curs_col_Jz,NbreColMatZ,Z_rowI)
        Appro_R(Ik,1:NbreColMatZ) = Z_rowI(1,1:NbreColMatZ) - ProdUlVl_row(1,1:NbreColMatZ) 
        Deallocate(Z_rowI)
        !!(2)
        Allocate(vect_R_tmp(1,NbreColMatZ))
        Do ccol=1,NbreColMatZ
            If (ANY(Indexes_Column(1:k-1) == ccol)) Then
                vect_R_tmp(1:1,ccol) =0
            Else
                vect_R_tmp(1:1,ccol) = Appro_R(Indexes_Row(k):Indexes_Row(k),ccol)
            EndIf
        EndDo                    
        loc_max_i_j = maxloc(abs(vect_R_tmp(1:1,1:NbreColMatZ))) 
        Indexes_Column(k) = loc_max_i_j(2)   !!J1
        Jk = Indexes_Column(k)     
        Deallocate(vect_R_tmp)
        !!(3)
        pivot = Appro_R(Ik,Jk)
        Do col=1,NbreColMatZ
            Matrix_V(k,col) = Appro_R(Ik,col)/pivot
        EndDo
        !!(4)
        ProdUlVl_col = 0
        Do l=1,k-1
            ProdUlVl_col(1:NbreLigMatZ,1) = ProdUlVl_col(1:NbreLigMatZ,1) + Matrix_V(l,Jk)*Matrix_U(1:NbreLigMatZ,l)
        EndDo
        !Appro_R(1:NbreLigMatZ,Jk) = Matrice2(1:NbreLigMatZ,Jk) - ProdUlVl_col(1:NbreLigMatZ,1)
        Allocate(Z_columnJ(NbreLigMatZ,1))
        curs_col_Green = curs_col_Jz+Jk-1
        Call computeBlockCol(Cells,curs_lig_Iz,NbreLigMatZ,curs_col_Green,Z_columnJ)
        Appro_R(1:NbreLigMatZ,Jk) = Z_columnJ(1:NbreLigMatZ,1) - ProdUlVl_col(1:NbreLigMatZ,1)
        Deallocate(Z_columnJ)
                    
        !!(5)
        Matrix_U(1:NbreLigMatZ,k) = Appro_R(1:NbreLigMatZ,Jk)
        !!(6)
        norm_prec_Z = NormeF_Zk(k-1)
        
        !! Frobenius_Norm with ZLANGE
        !!Frobenius_Norm(Matrix_U(1:NbreLigMatZ,k:k))
        !norm_Uk = ZLANGE('f',NbreLigMatZ,1,Matrix_U(1:NbreLigMatZ,k:k),NbreLigMatZ,work)
        !!Frobenius_Norm(Matrix_V(k:k,1:NbreColMatZ))
        !norm_Vk = ZLANGE('f',1,NbreColMatZ,Matrix_V(k:k,1:NbreColMatZ),1,work)
        
        !Frobenius_Norm(Matrix_U(1:NbreLigMatZ,k:k))
        val_2 = 0.
        Do ii =1,NbreLigMatZ
           val_2 = val_2 + (abs(Matrix_U(ii,k)))**2
        EndDo
        norm_Uk = sqrt(val_2)
        
        !Frobenius_Norm(Matrix_V(k:k,1:NbreColMatZ))
        val_2 = 0.
        Do jj=1,NbreColMatZ
            val_2 = val_2 + (abs(Matrix_V(k,jj)))**2
        EndDo
        norm_Vk = sqrt(val_2)
    
        sum_norm_UV = 0.
        Do nj=1,k-1
            Allocate(Unj(1,NbreLigMatZ));Allocate(Uk(NbreLigMatZ,1))
            Allocate(Vnj(NbreColMatZ,1));Allocate(Vk(1,NbreColMatZ))
            Unj = transpose(Matrix_U(1:NbreLigMatZ,nj:nj))        
            Uk = Matrix_U(1:NbreLigMatZ,k:k)
            Vnj = transpose(Matrix_V(nj:nj,1:NbreColMatZ))        
            Vk = Matrix_V(k:k,1:NbreColMatZ)
            sum_norm_UV(1:1,1:1) = sum_norm_UV(1:1,1:1) + abs(Matmul(Unj,Uk))*abs(Matmul(Vk,Vnj))
            Deallocate(Unj,Uk,Vnj,Vk)
        EndDo
        NormeF_Zk(k) = sqrt(norm_prec_Z**2 + 2*sum_norm_UV(1,1) + norm_Uk**2*norm_Vk**2)
        !!(7)                   
        Res = (norm_Uk*norm_Vk)/NormeF_Zk(k)
        If (Res <= Epsilon) Then
            Continuons = .FALSE.
        Else
            !!(8)
            Allocate(vect_R_tmp(NbreLigMatZ,1))
            Do llig=1,NbreLigMatZ
                If (ANY(Indexes_Row(1:k) == llig)) Then
                    vect_R_tmp(llig,1:1) = 0.
                Else
                    vect_R_tmp(llig,1:1) = Appro_R(llig,Indexes_Column(k):Indexes_Column(k))
                EndIf
            EndDo
            loc_max_i_j = maxloc(abs(vect_R_tmp(1:NbreLigMatZ,1:1)))
            Indexes_Row(k+1) = loc_max_i_j(1)
            k = k+1
            Deallocate(vect_R_tmp)
        EndIf
    EndDo 
    
    nb_iter_out = k
               
    Deallocate(Indexes_Row,Indexes_Column)
    Deallocate (Appro_R); Deallocate(NormeF_Zk)
    Deallocate(ProdUlVl_row, ProdUlVl_col)
    
    If (vrb_ACA == 1) Then
        Write(*,'(a,i4,a,i4,a,f12.6,a,i3,a)') 'I = ',I_Z,'; J = ',J_Z,' : Res = ',Res, ' after ', k, ' iterations.'
    EndIf

END SUBROUTINE Calcul_MatrixZij_ACA


SUBROUTINE Calcul_GenMatrix_ACA(m,n,Zin,nmax,nout,Matrix_U,Matrix_V)  
     
    USE Initialization
    USE common_variables
    USE iso_fortran_env

    ! "Implicit Statement" 
    IMPLICIT NONE

    !! IN/OUT ******************************************************************
    INTEGER, INTENT(IN) :: m,n,nmax
    COMPLEX(real64), Dimension(m,n), INTENT(IN):: Zin
    INTEGER, INTENT(OUT) :: nout
    COMPLEX(real64), Dimension(m,nmax), INTENT(OUT):: Matrix_U
    COMPLEX(real64), Dimension(nmax,n), INTENT(OUT):: Matrix_V
        
    ! Local
    Integer :: curs_lig_green,curs_col_green,col,ii,jj,k,Ik,Jk,nj,l,ccol,llig
    Integer, Dimension(:), allocatable ::Indexes_Row, Indexes_Column
    Integer, Dimension(2) :: loc_max_i_j    
    Real(kind=8) :: Epsilon,norm_prec_Z,norm_Uk,norm_Vk,val_2,Res
    real (kind = 4), external :: ZLANGE
    Real(kind=8), Dimension(:), allocatable :: NormeF_Zk
    Real(kind=8), Dimension(1,1) ::sum_norm_UV
    Logical :: Continuons
    DOUBLE PRECISION, Dimension(max(m,n)) :: WORK
    COMPLEX(real64) :: pivot
    COMPLEX(real64), Dimension(:,:), allocatable :: ProdUlVl_row, ProdUlVl_col,vect_R_tmp,Unj,Uk,Vnj,Vk
    COMPLEX(real64), Dimension(:,:), allocatable :: Appro_Z, Appro_R, Z_rowI, Z_columnJ
    
    Epsilon = Epsilon_ACA
                
    ! Indexes_Row et Indexes_Column are the arrays containing orderly selected row and column indexes of the matrix Zmn
    Allocate(Indexes_Row(nmax+1)); Allocate(Indexes_Column(nmax+1))
    Indexes_Row = 0;Indexes_Column = 0;
    !! line or column from Zij depending on the prgress in the algorithm
    Allocate(Z_rowI(1,n)); Allocate(Z_columnJ(m,1))
    !! Approximate Z matrix and Appproximate Error matrix
    Allocate (Appro_Z(m,n)); Allocate(Appro_R(m,n)); Appro_R = 0;
    Allocate(NormeF_Zk(nmax+1))
    Allocate(ProdUlVl_row(1,n));Allocate(ProdUlVl_col(m,1))
    
    !! Initialisation 
    Matrix_U=0; Matrix_V=0
    Indexes_Row(1) = 1  !!I1
    Appro_Z = 0
    !!
    Z_rowI = Zin(Indexes_Row(1):Indexes_Row(1),1:n)
    Appro_R(1:1,1:n) = Z_rowI
    Deallocate(Z_rowI)
    
    loc_max_i_j = maxloc(abs(Appro_R(Indexes_Row(1):Indexes_Row(1),1:n))) 
    Indexes_Column(1) = loc_max_i_j(2)   !!J1
    !! Initialisation de V
    pivot = Appro_R(Indexes_Row(1),Indexes_Column(1))
    Do col=1,n
        Matrix_V(1,col) = Appro_R(Indexes_Row(1),col)/pivot
    EndDo 
                
    Z_columnJ = Zin(1:m,Indexes_Column(1):Indexes_Column(1))
    Appro_R(1:m,Indexes_Column(1):Indexes_Column(1)) = Z_columnJ
    Deallocate(Z_columnJ)
    !! Initialisation de U
    Matrix_U(1:m,1:1) = Appro_R(1:m,Indexes_Column(1):Indexes_Column(1))
                
        
    !! Frobenius Norm of Z1
    !norm_prec_Z = 0.          
    !! Frobenius_Norm(Matrix_U(1:NbreLigMatZ,1:1))
    !norm_Uk = ZLANGE('f',NbreLigMatZ,1,Matrix_U(1:NbreLigMatZ,1:1),NbreLigMatZ,work)
    !! Frobenius_Norm(Matrix_V(1:1,1:NbreColMatZ))
    !norm_Vk = ZLANGE('f',1,NbreColMatZ,Matrix_V(1:1,1:NbreColMatZ),1,work)   
    
    norm_prec_Z = 0.
    ! Frobenius_Norm(Matrix_U(1:NbreLigMatZ,1:1))
    val_2 = 0.
    Do ii =1,m
        val_2 = val_2 + (abs(Matrix_U(ii,1)))**2
    EndDo
    norm_Uk = sqrt(val_2)
    
    ! Frobenius_Norm(Matrix_V(1:1,1:NbreColMatZ))
    val_2 = 0.
    Do jj=1,n
        val_2 = val_2 + (abs(Matrix_V(1,jj)))**2
    EndDo
    norm_Vk = sqrt(val_2)
    
    NormeF_Zk(1) = sqrt(norm_prec_Z**2 + norm_Uk**2*norm_Vk**2)
    
    !! attention max(R(i,Jk)), i!=I1(=1)
    Allocate(vect_R_tmp(m,1))
    vect_R_tmp(1:m,1:1) = Appro_R(1:m,Indexes_Column(1):Indexes_Column(1))
    vect_R_tmp(Indexes_Row(1),1) = 0
    loc_max_i_j = maxloc(abs(vect_R_tmp(1:m,1:1))) 
    Indexes_Row(2) = loc_max_i_j(1)
    Deallocate(vect_R_tmp)
                
    k=2
    continuons = .TRUE.
    Do while (continuons .AND. k <= Nb_it_max)
        !!(1)
        Ik = Indexes_Row(k)
        ProdUlVl_row = 0.
        Do l=1,k-1
            ProdUlVl_row(1,1:n) = ProdUlVl_row(1,1:n) + Matrix_U(Ik,l)*Matrix_V(l,1:n)
        EndDo
        Appro_R(Ik,1:n) = Zin(Ik,1:n) - ProdUlVl_row(1,1:n)
        
        !!(2)
        Allocate(vect_R_tmp(1,n))
        Do ccol=1,n
            If (ANY(Indexes_Column(1:k-1) == ccol)) Then
                vect_R_tmp(1:1,ccol) =0
            Else
                vect_R_tmp(1:1,ccol) = Appro_R(Indexes_Row(k):Indexes_Row(k),ccol)
            EndIf
        EndDo                    
        loc_max_i_j = maxloc(abs(vect_R_tmp(1:1,1:n))) 
        Indexes_Column(k) = loc_max_i_j(2)   !!J1
        Jk = Indexes_Column(k)     
        Deallocate(vect_R_tmp)
        !!(3)
        pivot = Appro_R(Ik,Jk)
        Do col=1,n
            Matrix_V(k,col) = Appro_R(Ik,col)/pivot
        EndDo
        !!(4)
        ProdUlVl_col = 0
        Do l=1,k-1
            ProdUlVl_col(1:m,1) = ProdUlVl_col(1:m,1) + Matrix_V(l,Jk)*Matrix_U(1:m,l)
        EndDo
        Appro_R(1:m,Jk) = Zin(1:m,Jk) - ProdUlVl_col(1:m,1)
                            
        !!(5)
        Matrix_U(1:m,k) = Appro_R(1:m,Jk)
        !!(6)
        norm_prec_Z = NormeF_Zk(k-1)
        
        !! Frobenius_Norm with ZLANGE
        !!Frobenius_Norm(Matrix_U(1:NbreLigMatZ,k:k))
        !norm_Uk = ZLANGE('f',NbreLigMatZ,1,Matrix_U(1:NbreLigMatZ,k:k),NbreLigMatZ,work)
        !!Frobenius_Norm(Matrix_V(k:k,1:NbreColMatZ))
        !norm_Vk = ZLANGE('f',1,NbreColMatZ,Matrix_V(k:k,1:NbreColMatZ),1,work)
        
        !Frobenius_Norm(Matrix_U(1:NbreLigMatZ,k:k))
        val_2 = 0.
        Do ii =1,m
           val_2 = val_2 + (abs(Matrix_U(ii,k)))**2
        EndDo
        norm_Uk = sqrt(val_2)
        
        !Frobenius_Norm(Matrix_V(k:k,1:NbreColMatZ))
        val_2 = 0.
        Do jj=1,n
            val_2 = val_2 + (abs(Matrix_V(k,jj)))**2
        EndDo
        norm_Vk = sqrt(val_2)
    
        sum_norm_UV = 0.
        Do nj=1,k-1
            Allocate(Unj(1,m));Allocate(Uk(m,1))
            Allocate(Vnj(n,1));Allocate(Vk(1,n))
            Unj = transpose(Matrix_U(1:m,nj:nj))        
            Uk = Matrix_U(1:m,k:k)
            Vnj = transpose(Matrix_V(nj:nj,1:n))        
            Vk = Matrix_V(k:k,1:n)
            sum_norm_UV(1:1,1:1) = sum_norm_UV(1:1,1:1) + abs(Matmul(Unj,Uk))*abs(Matmul(Vk,Vnj))
            Deallocate(Unj,Uk,Vnj,Vk)
        EndDo
        NormeF_Zk(k) = sqrt(norm_prec_Z**2 + 2*sum_norm_UV(1,1) + norm_Uk**2*norm_Vk**2)
        !!(7)                   
        Res = (norm_Uk*norm_Vk)/NormeF_Zk(k)
        If (Res <= Epsilon) Then
            Continuons = .FALSE.
        Else
            !!(8)
            Allocate(vect_R_tmp(m,1))
            Do llig=1,m
                If (ANY(Indexes_Row(1:k) == llig)) Then
                    vect_R_tmp(llig,1:1) = 0.
                Else
                    vect_R_tmp(llig,1:1) = Appro_R(llig,Indexes_Column(k):Indexes_Column(k))
                EndIf
            EndDo
            loc_max_i_j = maxloc(abs(vect_R_tmp(1:m,1:1)))
            Indexes_Row(k+1) = loc_max_i_j(1)
            k = k+1
            Deallocate(vect_R_tmp)
        EndIf
    EndDo 
    
    nout = k
               
    Deallocate(Indexes_Row,Indexes_Column)
    Deallocate (Appro_Z, Appro_R); Deallocate(NormeF_Zk)
    Deallocate(ProdUlVl_row, ProdUlVl_col)
    
    If (vrb_ACA == 1) Then
        Write(*,'(a,f12.6,a,i3,a)') 'Res = ',Res, ' after ', k, ' iterations.'
    EndIf

END SUBROUTINE Calcul_GenMatrix_ACA


SUBROUTINE Calcul_Matrix_ACA_SMW(nCells1,nCells2,CellsB1,CellsB2,Nb_it_max_smwf,Epsilon,nb_iter_out,Matrix_U,Matrix_V)  
     
    USE Initialization
    USE common_variables
    USE iso_fortran_env

    ! "Implicit Statement" 
    IMPLICIT NONE

    !! IN/OUT ******************************************************************
    !type (Cell), Dimension(Nbc), INTENT(IN) :: Cells
    INTEGER, INTENT(IN) :: nCells1,nCells2
    type (Cell), Dimension(nCells1), INTENT(IN) :: CellsB1
    type (Cell), Dimension(nCells2), INTENT(IN) :: CellsB2
    
    INTEGER, INTENT(IN) :: Nb_it_max_smwf
    Real(kind=8), INTENT(IN) :: Epsilon
    INTEGER, INTENT(OUT) :: nb_iter_out
    COMPLEX(real64), Dimension(3*nCells1,Nb_it_max_smwf), INTENT(OUT):: Matrix_U
    COMPLEX(real64), Dimension(Nb_it_max_smwf,3*nCells2), INTENT(OUT):: Matrix_V
        
    ! Local
    Integer :: curs_lig_green,curs_col_green,col,ii,jj,k,Ik,Jk,nj,l,ccol,llig
    Integer :: NbreLigMatZ, NbreColMatZ
    Integer, Dimension(:), allocatable ::Indexes_Row, Indexes_Column
    Integer, Dimension(2) :: loc_max_i_j    
    Real(kind=8) :: norm_prec_Z,norm_Uk,norm_Vk,val_2,Res
    real (kind = 4), external :: ZLANGE
    Real(kind=8), Dimension(:), allocatable :: NormeF_Zk
    Real(kind=8), Dimension(1,1) ::sum_norm_UV
    Logical :: Continuons
    DOUBLE PRECISION, Dimension(max(3*nCells1,3*nCells2)) :: WORK
    COMPLEX(real64) :: pivot
    COMPLEX(real64), Dimension(:,:), allocatable :: ProdUlVl_row, ProdUlVl_col,vect_R_tmp,Unj,Uk,Vnj,Vk
    COMPLEX(real64), Dimension(:,:), allocatable :: Appro_Z, Appro_R, Z_rowI, Z_columnJ
    
    vrb_ACA = 1;
        
    NbreLigMatZ = 3*nCells1;
    NbreColMatZ = 3*nCells2;
                
    ! Indexes_Row et Indexes_Column are the arrays containing orderly selected row and column indexes of the matrix Zmn
    Allocate(Indexes_Row(Nb_it_max_smwf+1)); Allocate(Indexes_Column(Nb_it_max_smwf+1))
    Indexes_Row = 0;Indexes_Column = 0;
    !! line or column from Zij depending on the prgress in the algorithm
    Allocate(Z_rowI(1,NbreColMatZ)); Allocate(Z_columnJ(NbreLigMatZ,1))
    !! Approximate Z matrix and Appproximate Error matrix
    Allocate (Appro_Z(NbreLigMatZ,NbreColMatZ)); Allocate(Appro_R(NbreLigMatZ,NbreColMatZ)); Appro_R = 0;
    Allocate(NormeF_Zk(Nb_it_max_smwf+1))
    Allocate(ProdUlVl_row(1,NbreColMatZ));Allocate(ProdUlVl_col(NbreLigMatZ,1))
    
    !! Initialisation 
    Matrix_U=0; Matrix_V=0
    Indexes_Row(1) = 1  !!I1
    Appro_Z = 0
    !! Pour initialiser la premiere ligne de R_appro, j'appelle Green_s_tr_partielle et Green_ns_pi_tr_partielle pour calculer Zij(1,:); 1 etant la premiere ligne de la matrice Zij qui concerne le bloc i,j de la matrice totale
    Call computeBlockRow_SMW(nCells1,nCells2,CellsB1,CellsB2,Indexes_Row(1),Z_rowI)  
    Appro_R(1:1,1:NbreColMatZ) = Z_rowI
    Deallocate(Z_rowI)
    loc_max_i_j = maxloc(abs(Appro_R(Indexes_Row(1):Indexes_Row(1),1:NbreColMatZ))) 
    Indexes_Column(1) = loc_max_i_j(2)   !!J1
    !! Initialisation de V
    pivot = Appro_R(Indexes_Row(1),Indexes_Column(1))
    Do col=1,NbreColMatZ
        Matrix_V(1,col) = Appro_R(Indexes_Row(1),col)/pivot
    EndDo 
                
    Call computeBlockCol_SMW(nCells1,nCells2,CellsB1,CellsB2,Indexes_Column(1),Z_columnJ)                                    
    Appro_R(1:NbreLigMatZ,Indexes_Column(1):Indexes_Column(1)) = Z_columnJ
    Deallocate(Z_columnJ)
    !! Initialisation de U
    Matrix_U(1:NbreLigMatZ,1:1) = Appro_R(1:NbreLigMatZ,Indexes_Column(1):Indexes_Column(1))
                
        
    !! Frobenius Norm of Z1
    !norm_prec_Z = 0.          
    !! Frobenius_Norm(Matrix_U(1:NbreLigMatZ,1:1))
    !norm_Uk = ZLANGE('f',NbreLigMatZ,1,Matrix_U(1:NbreLigMatZ,1:1),NbreLigMatZ,work)
    !! Frobenius_Norm(Matrix_V(1:1,1:NbreColMatZ))
    !norm_Vk = ZLANGE('f',1,NbreColMatZ,Matrix_V(1:1,1:NbreColMatZ),1,work)   
    
    norm_prec_Z = 0.
    ! Frobenius_Norm(Matrix_U(1:NbreLigMatZ,1:1))
    val_2 = 0.
    Do ii =1,NbreLigMatZ
        val_2 = val_2 + (abs(Matrix_U(ii,1)))**2
    EndDo
    norm_Uk = sqrt(val_2)
    
    ! Frobenius_Norm(Matrix_V(1:1,1:NbreColMatZ))
    val_2 = 0.
    Do jj=1,NbreColMatZ
        val_2 = val_2 + (abs(Matrix_V(1,jj)))**2
    EndDo
    norm_Vk = sqrt(val_2)
    
    NormeF_Zk(1) = sqrt(norm_prec_Z**2 + norm_Uk**2*norm_Vk**2)
    
    !! attention max(R(i,Jk)), i!=I1(=1)
    Allocate(vect_R_tmp(NbreLigMatZ,1))
    vect_R_tmp(1:NbreLigMatZ,1:1) = Appro_R(1:NbreLigMatZ,Indexes_Column(1):Indexes_Column(1))
    vect_R_tmp(Indexes_Row(1),1) = 0
    loc_max_i_j = maxloc(abs(vect_R_tmp(1:NbreLigMatZ,1:1))) 
    Indexes_Row(2) = loc_max_i_j(1)
    Deallocate(vect_R_tmp)
                
    k=2
    continuons = .TRUE.
    Do while (continuons .AND. k <= Nb_it_max_smwf)
        !!(1)
        Ik = Indexes_Row(k)
        ProdUlVl_row = 0.
        Do l=1,k-1
            ProdUlVl_row(1,1:NbreColMatZ) = ProdUlVl_row(1,1:NbreColMatZ) + Matrix_U(Ik,l)*Matrix_V(l,1:NbreColMatZ)
        EndDo
        
        Allocate(Z_rowI(1,NbreColMatZ))
        Call computeBlockRow_SMW(nCells1,nCells2,CellsB1,CellsB2,Ik,Z_rowI)  
        Appro_R(Ik,1:NbreColMatZ) = Z_rowI(1,1:NbreColMatZ) - ProdUlVl_row(1,1:NbreColMatZ) 
        Deallocate(Z_rowI)
        !!(2)
        Allocate(vect_R_tmp(1,NbreColMatZ))
        Do ccol=1,NbreColMatZ
            If (ANY(Indexes_Column(1:k-1) == ccol)) Then
                vect_R_tmp(1:1,ccol) =0
            Else
                vect_R_tmp(1:1,ccol) = Appro_R(Indexes_Row(k):Indexes_Row(k),ccol)
            EndIf
        EndDo                    
        loc_max_i_j = maxloc(abs(vect_R_tmp(1:1,1:NbreColMatZ))) 
        Indexes_Column(k) = loc_max_i_j(2)   !!J1
        Jk = Indexes_Column(k)     
        Deallocate(vect_R_tmp)
        !!(3)
        pivot = Appro_R(Ik,Jk)
        Do col=1,NbreColMatZ
            Matrix_V(k,col) = Appro_R(Ik,col)/pivot
        EndDo
        !!(4)
        ProdUlVl_col = 0
        Do l=1,k-1
            ProdUlVl_col(1:NbreLigMatZ,1) = ProdUlVl_col(1:NbreLigMatZ,1) + Matrix_V(l,Jk)*Matrix_U(1:NbreLigMatZ,l)
        EndDo
        
        Allocate(Z_columnJ(NbreLigMatZ,1))
        Call computeBlockCol_SMW(nCells1,nCells2,CellsB1,CellsB2,Jk,Z_columnJ)
        Appro_R(1:NbreLigMatZ,Jk) = Z_columnJ(1:NbreLigMatZ,1) - ProdUlVl_col(1:NbreLigMatZ,1)
        Deallocate(Z_columnJ)
                    
        !!(5)
        Matrix_U(1:NbreLigMatZ,k) = Appro_R(1:NbreLigMatZ,Jk)
        !!(6)
        norm_prec_Z = NormeF_Zk(k-1)
        
        !! Frobenius_Norm with ZLANGE
        !!Frobenius_Norm(Matrix_U(1:NbreLigMatZ,k:k))
        !norm_Uk = ZLANGE('f',NbreLigMatZ,1,Matrix_U(1:NbreLigMatZ,k:k),NbreLigMatZ,work)
        !!Frobenius_Norm(Matrix_V(k:k,1:NbreColMatZ))
        !norm_Vk = ZLANGE('f',1,NbreColMatZ,Matrix_V(k:k,1:NbreColMatZ),1,work)
        
        !Frobenius_Norm(Matrix_U(1:NbreLigMatZ,k:k))
        val_2 = 0.
        Do ii =1,NbreLigMatZ
           val_2 = val_2 + (abs(Matrix_U(ii,k)))**2
        EndDo
        norm_Uk = sqrt(val_2)
        
        !Frobenius_Norm(Matrix_V(k:k,1:NbreColMatZ))
        val_2 = 0.
        Do jj=1,NbreColMatZ
            val_2 = val_2 + (abs(Matrix_V(k,jj)))**2
        EndDo
        norm_Vk = sqrt(val_2)
    
        sum_norm_UV = 0.
        Do nj=1,k-1
            Allocate(Unj(1,NbreLigMatZ));Allocate(Uk(NbreLigMatZ,1))
            Allocate(Vnj(NbreColMatZ,1));Allocate(Vk(1,NbreColMatZ))
            Unj = transpose(Matrix_U(1:NbreLigMatZ,nj:nj))        
            Uk = Matrix_U(1:NbreLigMatZ,k:k)
            Vnj = transpose(Matrix_V(nj:nj,1:NbreColMatZ))        
            Vk = Matrix_V(k:k,1:NbreColMatZ)
            sum_norm_UV(1:1,1:1) = sum_norm_UV(1:1,1:1) + abs(Matmul(Unj,Uk))*abs(Matmul(Vk,Vnj))
            Deallocate(Unj,Uk,Vnj,Vk)
        EndDo
        NormeF_Zk(k) = sqrt(norm_prec_Z**2 + 2*sum_norm_UV(1,1) + norm_Uk**2*norm_Vk**2)
        !!(7)                   
        Res = (norm_Uk*norm_Vk)/NormeF_Zk(k)
        If (Res <= Epsilon) Then
            Continuons = .FALSE.
        Else
            !!(8)
            Allocate(vect_R_tmp(NbreLigMatZ,1))
            Do llig=1,NbreLigMatZ
                If (ANY(Indexes_Row(1:k) == llig)) Then
                    vect_R_tmp(llig,1:1) = 0.
                Else
                    vect_R_tmp(llig,1:1) = Appro_R(llig,Indexes_Column(k):Indexes_Column(k))
                EndIf
            EndDo
            loc_max_i_j = maxloc(abs(vect_R_tmp(1:NbreLigMatZ,1:1)))
            Indexes_Row(k+1) = loc_max_i_j(1)
            k = k+1
            Deallocate(vect_R_tmp)
        EndIf
    EndDo 
    
    nb_iter_out = min(Nb_it_max_smwf,k)
               
    Deallocate(Indexes_Row,Indexes_Column)
    Deallocate (Appro_Z, Appro_R); Deallocate(NormeF_Zk)
    Deallocate(ProdUlVl_row, ProdUlVl_col)
    
    If (vrb_ACA == 1) Then
        Write(*,'(a,f12.6,a,i3,a)') 'Res = ',Res, ' after ', k, ' iterations.'
    EndIf

END SUBROUTINE Calcul_Matrix_ACA_SMW