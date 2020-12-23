INCLUDE 'mkl_pardiso.f90'
SUBROUTINE pardiso_solver(n,nrhs,nnz,mtype,iparm3,ia,ja,a,b,x)

  !*******************************************************************************   
  !   Copyright(C) 2004-2013 Intel Corporation. All Rights Reserved.
  !   
  !   The source code, information  and  material ("Material") contained herein is
  !   owned  by Intel Corporation or its suppliers or licensors, and title to such
  !   Material remains  with Intel Corporation  or its suppliers or licensors. The
  !   Material  contains proprietary information  of  Intel or  its  suppliers and
  !   licensors. The  Material is protected by worldwide copyright laws and treaty
  !   provisions. No  part  of  the  Material  may  be  used,  copied, reproduced,
  !   modified, published, uploaded, posted, transmitted, distributed or disclosed
  !   in any way  without Intel's  prior  express written  permission. No  license
  !   under  any patent, copyright  or  other intellectual property rights  in the
  !   Material  is  granted  to  or  conferred  upon  you,  either  expressly,  by
  !   implication, inducement,  estoppel or  otherwise.  Any  license  under  such
  !   intellectual  property  rights must  be express  and  approved  by  Intel in
  !   writing.
  !   
  !   *Third Party trademarks are the property of their respective owners.
  !   
  !   Unless otherwise  agreed  by Intel  in writing, you may not remove  or alter
  !   this  notice or  any other notice embedded  in Materials by Intel or Intel's
  !   suppliers or licensors in any way.
  !
  !*******************************************************************************
  !   Content : MKL PARDISO Fortran-90 example
  !
  !*******************************************************************************
  !----------------------------------------------------------------------
  ! This subroutine is a modified version of the Example program to show the 
  ! use of the "PARDISO" routine for symmetric linear systems
  ! /opt/intel/composer_xe_2013.5.192/mkl/examples/solverf/source
  !---------------------------------------------------------------------
  
  USE mkl_pardiso
  IMPLICIT NONE
  
  !IN/OUT
  INTEGER, INTENT(IN) :: n,nnz,nrhs,mtype,iparm3
  INTEGER, Dimension(n+1), INTENT(IN) :: ia
  INTEGER, Dimension(nnz), INTENT(IN) :: ja  
  COMPLEX, Dimension(nnz), INTENT(IN) :: a
  COMPLEX, Dimension(n,nrhs), INTENT(INOUT) :: b
  COMPLEX(kind=4), Dimension(n,nrhs), INTENT(INOUT) :: x    
  
  !.. Internal solver memory pointer 
  TYPE(MKL_PARDISO_HANDLE), ALLOCATABLE  :: pt(:)
  !.. All other variables
  INTEGER maxfct, mnum, phase, error, msglvl
  INTEGER error1
  INTEGER, ALLOCATABLE :: iparm( : )
  !INTEGER, ALLOCATABLE :: ia( : )
  !INTEGER, ALLOCATABLE :: ja( : )
  !COMPLEX, ALLOCATABLE :: a( : )
  !COMPLEX, ALLOCATABLE :: b( :,: )
  !COMPLEX, ALLOCATABLE :: x( :,: )
  INTEGER i,j,c, idum(1)
  COMPLEX ddum(1)
  !.. Fill all arrays containing matrix data.
  !n = 8 
  !nnz = 18
  !nrhs = 1 
  
  !Write(*,*) 'n =',n
  !Write(*,*) 'nnz =',nnz
  !Write(*,*) 'nrhs =',nrhs
  !Write(*,*) 'mtype = ',mtype
  !Write(*,*) 'iparm3 =',iparm3
  !Write(*,*) 'ia(1:10) =', ia(1:10)
  !Write(*,*) 'ja(1:10) =', ja(1:10)
  !Write(*,*) 'a(1:10) =', a(1:10)
  
  
  !Write(*,*) 'ia = '
  !Write(*,*) ia(1:n+1);
  !Write(*,*) 'ja = '
  !Write(*,*) ja(1:nnz);
  !Write(*,*) 'a = '
  !Do i=1,nnz
  !  Write(*,*) a(i);
  !EndDo
  !Write(*,*) 'b ='
  !DO i = 1, n
  !  Do c=1,nrhs
  !   WRITE(*,*) ' b(',i,',',c,') = ', b(i,c)
  !  End DO
  !END DO
  
  
  !ALLOCATE( ia ( n + 1 ) )
  !ia = (/ 1, 5, 8, 10, 12, 15, 17, 18, 19 /)
  !ALLOCATE( ja ( nnz ) )
  !ja = (/ 1,    3,       6, 7,    &
  !           2, 3,    5,          &
  !              3,             8, &
  !                 4,       7,    &
  !                    5, 6, 7,    &
  !                       6,    8, &
  !                          7,    &
  !                             8 /)
  !ALLOCATE( a ( nnz ) )
  !a = (/ 7.d0,        1.d0,             2.d0, 7.d0,        &
  !             -4.d0, 8.d0,       2.d0,                    &
  !                    1.d0,                         5.d0,  &
  !                          7.d0,             9.d0,        &
  !                                5.d0, 1.d0, 5.d0,        &
  !                                     -1.d0,       5.d0,  &
  !                                           11.d0,        &
  !                                                  5.d0 /)
  !ALLOCATE( b ( n,nrhs ) )
  !ALLOCATE( x ( n,nrhs ) )
  !..
  !.. Set up PARDISO control parameter
  !..
  ALLOCATE( iparm ( 64 ) )
  
  do i = 1, 64
     iparm(i) = 0
  end do 
  
  iparm(1) = 1 ! no solver default
  iparm(2) = 3 ! fill-in reordering from METIS
  iparm(3) = iparm3 ! numbers of processors, value of OMP_NUM_THREADS
  iparm(4) = 0 ! no iterative-direct algorithm
  iparm(5) = 0 ! no user fill-in reducing permutation
  iparm(6) = 0 ! =0 solution on the first n compoments of x
  iparm(7) = 0! Not in use
  iparm(8) = 20 ! numbers of iterative refinement steps
  iparm(10) = 13 ! perturbe the pivot elements with 1E-13
  iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
  iparm(13) = 1 ! maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm(13) = 1 in case of inappropriate accuracy
  iparm(14) = 0 ! Output: number of perturbed pivots
  iparm(15) = 0 ! Not in use
  iparm(16) = 0 ! Not in use
  iparm(17) = 0 ! Not in use
  iparm(18) = -1 ! Output: number of nonzeros in the factor LU
  iparm(19) = -1 ! Output: Mflops for LU factorization
  iparm(20) = 0 ! Output: Numbers of CG Iterations
  iparm(21) = 1;  !Apply 1x1 and 2x2 Bunch and Kaufman pivoting during the factorization process
  iparm(24) = 1;  !PARDISO uses new two - level factorization algorithm
  
  maxfct = 1 !Maximum number of numerical factorizations
  mnum = 1  !Which factorization to use
  
  error  = 0 ! initialize error flag
  msglvl = 0 ! print statistical information
  !mtype  = 6 !complex symmetric matrix
  !mtype  = 13 ! complex nonsymmetric matrix
  
  !.. Initiliaze the internal solver memory pointer. This is only
  ! necessary for the FIRST call of the PARDISO solver.
  
  ALLOCATE ( pt ( 64 ) )
  do i = 1, 64
     pt( i )%DUMMY =  0 
  end do
  
  !.. Reordering and Symbolic Factorization, This step also allocates
  ! all memory that is necessary for the factorization
  
  phase = 11 ! only reordering and symbolic factorization
  CALL pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja, &
  idum, nrhs, iparm, msglvl, ddum, ddum, error)
     
  !WRITE(*,*) 'Reordering completed ... '
  IF (error /= 0) THEN
     WRITE(*,*) 'The following ERROR was detected: ', error
     GOTO 1000
  END IF
  !WRITE(*,*) 'Number of nonzeros in factors = ',iparm(18)
  !WRITE(*,*) 'Number of factorization MFLOPS = ',iparm(19)
  
  !.. Factorization.
  phase = 22 ! only factorization
  CALL pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja, &
  idum, nrhs, iparm, msglvl, ddum, ddum, error)
  !WRITE(*,*) 'Factorization completed ... '
  IF (error /= 0) THEN
     WRITE(*,*) 'The following ERROR was detected: ', error
     GOTO 1000
  ENDIF
  
  !.. Back substitution and iterative refinement
  phase = 33 ! only factorization
  
  CALL pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja, &
  idum, nrhs, iparm, msglvl, b, x, error)
  !WRITE(*,*) 'Solve completed ... '
  IF (error /= 0) THEN
     WRITE(*,*) 'The following ERROR was detected: ', error
     GOTO 1000
  ENDIF
  !WRITE(*,*) 'The solution of the system is '
  !DO i = 1, n
  !  Do c=1,nrhs
  !   WRITE(*,*) ' x(',i,',',c,') = ', x(i,c)
  !  End DO
  !END DO
        
  1000 CONTINUE
  !.. Termination and release of memory
  phase = -1 ! release internal memory
  CALL pardiso (pt, maxfct, mnum, mtype, phase, n, ddum, idum, idum, &
  idum, nrhs, iparm, msglvl, ddum, ddum, error1)
  
  !IF ( ALLOCATED( ia ) )      DEALLOCATE( ia )
  !IF ( ALLOCATED( ja ) )      DEALLOCATE( ja )
  !IF ( ALLOCATED( a ) )       DEALLOCATE( a )
  !IF ( ALLOCATED( b ) )       DEALLOCATE( b )
  !IF ( ALLOCATED( x ) )       DEALLOCATE( x )
  IF ( ALLOCATED( iparm ) )   DEALLOCATE( iparm )
  
  IF (error1 /= 0) THEN
     WRITE(*,*) 'The following ERROR on release stage was detected: ', error1
  ENDIF
  
  !IF ( error /= 0 ) STOP 1
  !STOP 0

END SUBROUTINE 
