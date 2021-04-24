SUBROUTINE Compute_Scattered_Fields(nom_methode,Cells,E_total,CBFM_Blocks,MPI_CBFM_Blocks,Transmitters,Receivers,Es_total)
    
    USE Initialization
    USE common_variables
    USE iso_fortran_env
    USE MPI

    IMPLICIT NONE

    !IN/OUT
    character(8), INTENT(IN):: nom_methode
    type (Cell), Dimension(Nbc), INTENT(IN) :: Cells
    COMPLEX(real64), Dimension(3*Nbc_proc,2*NTr), INTENT(IN):: E_total
    ! we will need to know the blocks and MPI blocks repartition to have the good corresponding Cells_proc
    type (CBFM_Block), Dimension(Nblocks), INTENT(IN) :: CBFM_Blocks
    Integer, Dimension(nber_procs,Nblk_proc_max+1), INTENT(IN):: MPI_CBFM_Blocks
    type (Dipole), Dimension(NTr), INTENT(IN) :: Transmitters
    type (Dipole), Dimension(NRx_tot), INTENT(IN) :: Receivers
    COMPLEX(real64), Dimension(NRx_tot,4*NTr), INTENT(OUT) :: Es_total

    ! Local 
    Integer :: R, ii,K, num_capteur_fichier,num_file,I,Ic
    Integer :: Lig,id,nthreads,p,d,Nelts,cel_beg,cel_end
    Integer :: kk_job, kk, curs_cel,Nbc_b
    Integer, Dimension(nber_procs) :: all_Nbc_procs        
    character(200) :: file_name
    CHARACTER(:), allocatable::nom_meth_exact
    COMPLEX(real64), Dimension(:,:), allocatable :: Green_dt, Es_total_all
    COMPLEX(real64), Dimension(:), allocatable ::Es_total_capteur,Es_total_capteur_all
    COMPLEX(real64), Dimension(3) :: E_v, E_h, E_v_p, E_h_p
    COMPLEX(real64) :: Vv, Vh, Hv, Hh
    Integer :: num_emetteur, num_capteur
    Real(kind=8) :: theta_capteur, phi_capteur
    type(Cell), Dimension(:), allocatable :: Cells_proc
  
   
    ! R (m) is the distance between the scatterer and the receivers 
    R = 10000; ! 10 Km 
    
    If (rank == 0) Then 
        Write (*,*) ''
        Write (*,*) '----------------------- Scattred Fields  --------------------------'
    EndIf
    
    
    ! All procs recover again this important information 
    call MPI_ALLGATHER (Nbc_proc,1,MPI_INTEGER,all_Nbc_procs,1,MPI_INTEGER,MPI_COMM_WORLD,code);    
    !Write(*,*) 'Nbc_proc =',all_Nbc_procs
    
    ! now every proc focus on its cells/part of E_tot 
    Allocate(Cells_proc(Nbc_proc));
    MyNBlocks = MPI_CBFM_Blocks(rank+1,1);
    curs_cel = 1;
    Do kk_job = 1,MyNBlocks
        kk = MPI_CBFM_Blocks(rank+1,1+kk_job);
        Nbc_b = CBFM_Blocks(kk)%Nbc_b
        cel_beg = sum(CBFM_Blocks(1:kk-1)%Nbc_b)+1;
        cel_end = sum(CBFM_Blocks(1:kk)%Nbc_b);
        Cells_proc(curs_cel:curs_cel+Nbc_b-1) = Cells(cel_beg:cel_end);
        curs_cel = curs_cel + Nbc_b; 
    EndDo
        
    
    DO num_capteur =1,NRx_tot  
	Allocate(Es_total_capteur(4*NTr),Es_total_capteur_all(4*NTr));   
        
        theta_capteur = Receivers(num_capteur)%theta
        phi_capteur = Receivers(num_capteur)%phi
        
        !! Dyade de Green singuliere
        ! we store the 3*3*N here just to avoid recomputing each 3X3 matrix for each transmitter position 
        Allocate(Green_dt(3,3*Nbc_proc))  
        Call Green_s_dt(Nbc_proc,Cells_proc, R,theta_capteur,phi_capteur,Green_dt)
                   
        DO num_emetteur=1,NTr
            E_v = 0
            E_h = 0      
            DO I=1,3*Nbc_proc ! This the difference with the serial OpenMP code, Here each process computes S based on its part of E_total
                E_v(1)=E_v(1)+Green_dt(1,I)*E_total(I,num_emetteur)
                E_v(2)=E_v(2)+Green_dt(2,I)*E_total(I,num_emetteur)
                E_v(3)=E_v(3)+Green_dt(3,I)*E_total(I,num_emetteur)
            
                E_h(1)=E_h(1)+Green_dt(1,I)*E_total(I,num_emetteur+NTr)
                E_h(2)=E_h(2)+Green_dt(2,I)*E_total(I,num_emetteur+NTr)
                E_h(3)=E_h(3)+Green_dt(3,I)*E_total(I,num_emetteur+NTr)    
            ENDDO
            
            !! ---------------------------------------------------------------------------------!!
            !! ------------------------------Polarisation Vv -----------------------------------!!
            !! ---------------------------------------------------------------------------------!! 
            Vv= - E_v(1)*sin(theta_capteur*Pi/180.)+E_v(2)*cos(theta_capteur*Pi/180.)*cos(phi_capteur*Pi/180.) &
            +E_v(3)*cos(theta_capteur*Pi/180.)*sin(phi_capteur*Pi/180.);
                
    
            !! ---------------------------------------------------------------------------------!!
            !! ------------------------------Polarisation Vh -----------------------------------!!
            !! ---------------------------------------------------------------------------------!! 
            Vh= - E_v(2)*sin(phi_capteur*Pi/180.)+E_v(3)*cos(phi_capteur*Pi/180.);
                
    
            !! ---------------------------------------------------------------------------------!!
            !! ------------------------------Polarisation Hv -----------------------------------!!
            !! ---------------------------------------------------------------------------------!! 
            Hv= - E_h(1)*sin(theta_capteur*Pi/180.)+E_h(2)*cos(theta_capteur*Pi/180.)*cos(phi_capteur*Pi/180.) &
            + E_h(3)*cos(theta_capteur*Pi/180.)*sin(phi_capteur*Pi/180.);
                
            !! ---------------------------------------------------------------------------------!!
            !! ------------------------------Polarisation Hh -----------------------------------!!
            !! ---------------------------------------------------------------------------------!! 
            Hh= - E_h(2)*sin(phi_capteur*Pi/180.) + E_h(3)*cos(phi_capteur*Pi/180.);
             
            !! ---------------------------------------------------------------------------------!!
            !! -----------------------Remplissage du vecteur S_total ---------------------------!!
            !! ---------------------------------------------------------------------------------!!
            
            Es_total_capteur(4*(num_emetteur-1)+1)  = Vv; 
            Es_total_capteur(4*(num_emetteur-1)+2)  = Vh; 
            Es_total_capteur(4*(num_emetteur-1)+3)  = Hv; 
            Es_total_capteur(4*(num_emetteur-1)+4)  = Hh; 
        ENDDO      
        Deallocate(Green_dt);        
        Nelts = 4*NTr;
        Call MPI_ALLREDUCE(Es_total_capteur,Es_total_capteur_all,Nelts,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,code);
        Es_total(num_capteur,1:4*NTr) = Es_total_capteur_all(1:4*NTr);
        deallocate(Es_total_capteur,Es_total_capteur_all);        
    Enddo
    Deallocate(Cells_proc);
    
    ! I was having insufficient memory problem (buffer) with this approach (calculate all S_total and then in 1 communicate calculate the sum)
    !Allocate(S_total_all(NRx_tot,4*NTr));    
    !Call MPI_BARRIER(MPI_COMM_WORLD,code);
    !Nelts = 4*NRx_tot*NTr;
    !Call MPI_ALLREDUCE(S_total,S_total_all,Nelts,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,code);
        
    !S_total(1:NRx_tot,1:4*NTr) = S_total_all(1:NRx_tot,1:4*NTr);
    !deallocate(S_total_all);  
    
End Subroutine Compute_Scattered_Fields
