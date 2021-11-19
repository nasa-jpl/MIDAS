SUBROUTINE Compute_Scattering_Matrices(nom_methode,Cells,E_total,CBFM_Blocks,MPI_CBFM_Blocks,Transmitters,Receivers,S_total)
    
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
    COMPLEX(real64), Dimension(3) :: E_v, E_h, E_v_p, E_h_p
    COMPLEX(real64) :: ff_coef, Vv, Vh, Hv, Hh
    Integer :: num_emetteur, num_capteur
    Real(kind=8) :: theta_capteur, phi_capteur, Beta,step_beta
    type(Cell), Dimension(:), allocatable :: Cells_proc
  
   
    If (rank == 0) Then 
        Write (*,*) ''
        Write (*,*) '----------------------- Scattring Matrices --------------------------'
    EndIf
    
    
    if (NPolBeta > 1) then
        step_beta = (beta_final_Pol-beta_init_Pol)/(NPolBeta-1);
    else
        step_beta = 0;
        
    endif
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
            
            !Do num_pol = 1:NPolBeta
            !    Beta =  beta_init_Pol + (num_pol-1)*step_beta
                
                ! need to figure out E_v_pol and E_h_pol (x, y and z)
                !E_v_pol = cos(Beta*Pi/180.)*E_v + sin(Beta*Pi/180.)*E_h
                !E_h_pol =  
                
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
    
End Subroutine Compute_Scattering_Matrices


SUBROUTINE GetFFieldCoeff(Nc,Cells_in,theta_capteur,phi_capteur,ff_coeffs)

    USE Initialization
    USE common_variables
    IMPLICIT NONE
    
    Integer, INTENT(IN) :: Nc
    type (Cell), Dimension(Nc), INTENT(IN) :: Cells_in
    Real(kind=8), INTENT(IN) :: theta_capteur,phi_capteur
    COMPLEX(real64), Dimension(Nc), INTENT(OUT) :: ff_coeffs
    
    !Local
    Integer Is
    Real(kind=8) :: xc,yc,zc,x_cap,y_cap,z_cap
    COMPLEX(real64) :: ffc, f_kapChe

    DO Is=1, Nc
             
        ! Cell in scatterer x, y & z         
        xc = Cells_in(Is)%Xc;
        yc = Cells_in(Is)%Yc;
        zc = Cells_in(Is)%Zc;
        
        !! Receiver x, y & z : since 3/13/2019 !
        x_cap = cos(theta_capteur*Pi/180.) 
        y_cap = sin(theta_capteur*Pi/180.)*cos(phi_capteur*Pi/180.) 
        z_cap = sin(theta_capteur*Pi/180.)*sin(phi_capteur*Pi/180.)
        
        ! translation theorem E2= E1*exp(-ik*delta.u) 
        ffc= exp(-J*k_0*(x_cap*xc+y_cap*yc+z_cap*zc))    ! ici ce n'est pas Gr_mn le terme qui depends de rmn est sorti a l'exterieur de S on l'a plus ici
                                                                ! ceci est le terme de dephasage du theoreme de translation E2(cell_i) = E1(0)*exp(-jk delta u)
                                                                ! le terme de green est mnt en fait a l'exterieur vu la definition de la matrice S !
                                                                ! en fait c'est une methode totalement differente de ce qu'on a ustilise pour le champ diffracte a r (inside or outside the scatterer)
        f_kapChe= Cells_in(Is)%Kappa_n*Cells_in(Is)%Che_n
          
        ff_coeffs(Is) = ffc*k_0**2.*f_kapChe /(4*Pi) * k_0    ! the last k_0 comes from the definition of the S matrix with DDSCAT  
        
    ENDDO

END SUBROUTINE GetFFieldCoeff