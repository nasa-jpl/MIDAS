SUBROUTINE Compute_Scattering_Matrices_init(nom_methode,Cells,E_total,CBFM_Blocks,MPI_CBFM_Blocks,Transmitters,Receivers,S_total)
    
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
    !COMPLEX(real64), Dimension(NPolBeta*NRx_tot,4*NTr), INTENT(OUT) :: S_total     ! change 11/10/2021  ! originally Dimension(NRx_tot,4*NTr)
    COMPLEX(real64), Dimension(NRx_tot,4*NTr), INTENT(OUT) :: S_total
    ! Local 
    Integer :: ii,K,I,Ic,ix,iy,iz
    Integer :: Lig,id,nthreads,p,d,Nelts,cel_beg,cel_end
    Integer :: kk_job, kk, curs_cel,Nbc_b
    Integer, Dimension(nber_procs) :: all_Nbc_procs        
    character(200) :: file_name
    CHARACTER(:), allocatable::nom_meth_exact
    COMPLEX(real64), Dimension(:,:), allocatable :: S_total_all
    COMPLEX(real64), Dimension(:), allocatable :: ff_coeffs, S_total_capteur, S_total_capteur_all
    COMPLEX(real64), Dimension(3) :: E_v, E_h
    COMPLEX(real64) :: ff_coef, Vv, Vh, Hv, Hh
    Integer :: num_emetteur, num_capteur
    Real(kind=8) :: theta_capteur, phi_capteur
    type(Cell), Dimension(:), allocatable :: Cells_proc
  
   
    If (rank == 0) Then 
        Write (*,*) ''
        Write (*,*) '----------------------- Scattering Matrices --------------------------'
    EndIf
    
    ! All procs recover again this important information 
    call MPI_ALLGATHER (Nbc_proc,1,MPI_INTEGER,all_Nbc_procs,1,MPI_INTEGER,MPI_COMM_WORLD,code);    
    
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
	    Allocate(S_total_capteur(4*NTr),S_total_capteur_all(4*NTr));   
        Allocate(ff_coeffs(Nbc_proc)); 
        !! Dyade de Greene singuliere
        theta_capteur = Receivers(num_capteur)%theta
        phi_capteur = Receivers(num_capteur)%phi
        
        Call GetFFieldCoeff(Nbc_proc,Cells_proc,theta_capteur,phi_capteur,ff_coeffs)          
        DO num_emetteur=1,NTr
            E_v = 0
            E_h = 0      
            DO Ic=1,Nbc_proc ! This the difference with the serial OpenMP code, Here each process computes S based on its part of E_total            
                ff_coef = ff_coeffs(Ic);
                ix = 3*(Ic-1)+1
                iy = 3*(Ic-1)+2
                iz = 3*Ic;
                
                E_v(1)=E_v(1)+ ff_coef*E_total(ix,num_emetteur)
                E_v(2)=E_v(2)+ ff_coef*E_total(iy,num_emetteur)
                E_v(3)=E_v(3)+ ff_coef*E_total(iz,num_emetteur)
            
                E_h(1)=E_h(1)+ ff_coef*E_total(ix,num_emetteur+NTr)
                E_h(2)=E_h(2)+ ff_coef*E_total(iy,num_emetteur+NTr)
                E_h(3)=E_h(3)+ ff_coef*E_total(iz,num_emetteur+NTr)                    
            ENDDO
                
            ! remember in the scattering matrix you are linking Ei and Es ! whatever their polarization is !
            
            ! Look at page 3 of convention code particle : you will have to update the coordinates of Evs and Ehs, then calculating the 4 components of the S matrix will be quite easy
          
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
            
            S_total_capteur(4*(num_emetteur-1)+1)  = Vv; 
            S_total_capteur(4*(num_emetteur-1)+2)  = Vh; 
            S_total_capteur(4*(num_emetteur-1)+3)  = Hv; 
            S_total_capteur(4*(num_emetteur-1)+4)  = Hh; 

        ENDDO      
        Deallocate(ff_coeffs); 
        Nelts = 4*NTr;
        Call MPI_ALLREDUCE(S_total_capteur,S_total_capteur_all,Nelts,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,code);
        S_total(num_capteur,1:4*NTr) = S_total_capteur_all(1:4*NTr);
        deallocate(S_total_capteur,S_total_capteur_all);        
    Enddo
    Deallocate(Cells_proc);
    
    ! I was having insufficient memory problem (buffer) with this approach (calculate all S_total and then in 1 communicate calculate the sum)
    !Allocate(S_total_all(NRx_tot,4*NTr));    
    !Call MPI_BARRIER(MPI_COMM_WORLD,code);
    !Nelts = 4*NRx_tot*NTr;
    !Call MPI_ALLREDUCE(S_total,S_total_all,Nelts,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,code);
        
    !S_total(1:NRx_tot,1:4*NTr) = S_total_all(1:NRx_tot,1:4*NTr);
    !deallocate(S_total_all);  
    
    
    ! write S files (txt for now for 3rd angle rotation debug)
    
    
End Subroutine Compute_Scattering_Matrices_init