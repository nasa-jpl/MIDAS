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