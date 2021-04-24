SUBROUTINE Compute_ExtAbsCsec_fromIntField(nom_methode,Cells,E_total,Transmitters,C_ext,C_abs)
    
    USE Initialization
    USE common_variables
    USE iso_fortran_env
    USE MPI

    IMPLICIT NONE

    !IN/OUT
    character(8), INTENT(IN):: nom_methode
    type (Cell), Dimension(Nbc), INTENT(IN) :: Cells
    COMPLEX(real64), Dimension(3*Nbc_proc,2*NTr), INTENT(IN):: E_total
    type (Dipole), Dimension(NTr), INTENT(IN) :: Transmitters
    COMPLEX(real64), Dimension(NTr), INTENT(OUT) :: C_ext,C_abs

    ! Local 
    Integer :: ii,K,I,Ic,Lig,id,nthreads,p,d 
    Integer :: cel_beg, cel_end
    character(200) :: file_name
    CHARACTER(:), allocatable::nom_meth_exact
    Integer :: num_emetteur, num_capteur
    Integer, Dimension(nber_procs) :: all_Nbc_procs
    Real(kind=8) :: Cext_e_V,Cext_e_H, Cabs_e_V,Cabs_e_H; 
    Real(kind=8) :: theta_capteur, phi_capteur
    COMPLEX(real64), Dimension(:), allocatable :: C_ext_all,C_abs_all
    COMPLEX(real64), Dimension(:,:),allocatable::E_ref_incident
    type(Cell), Dimension(:), allocatable :: Cells_proc
    
    if (rank == 0) then 
        Write (*,*) '-------- Extinction and Absorption from Internal Fields ----------'  
        Write (*,*) ''
    endif
    
    ! All procs recover again this important information 
    call MPI_ALLGATHER (Nbc_proc,1,MPI_INTEGER,all_Nbc_procs,1,MPI_INTEGER,MPI_COMM_WORLD,code);    
        
    ! now every proc focus on its cells/part of E_tot 
    Allocate(Cells_proc(Nbc_proc));
    cel_beg = sum(all_Nbc_procs(1:rank))+1;
    cel_end = sum(all_Nbc_procs(1:rank+1));
    Cells_proc(1:Nbc_proc) = Cells(cel_beg:cel_end)
    
               
    DO num_emetteur=1,NTr
        Cext_e_V = 0;Cext_e_H = 0;
        Cabs_e_V =0;Cabs_e_H =0;
        
        Allocate(E_ref_incident(3*Nbc_proc,2));
        Call Incident_Field(1,Nbc_proc,Cells_proc,NTr,Transmitters,num_emetteur,num_emetteur,E_ref_incident);
    
        DO I=1,Nbc_proc
            Cabs_e_V = Cabs_e_V + imag(Cells_proc(I)%Che_n)*abs(sum(E_total(3*(I-1)+1:3*I,num_emetteur)))**2.*Cells_proc(I)%Sc**3. ;  
            Cabs_e_H = Cabs_e_H + imag(Cells_proc(I)%Che_n)*abs(sum(E_total(3*(I-1)+1:3*I,num_emetteur+ &
                NTr)))**2.*Cells_proc(I)%Sc**3. ;  
            
            Cext_e_V = Cext_e_V + imag(Cells_proc(I)%Che_n*sum(E_total(3*(I-1)+1:3*I,num_emetteur))&
                *conjg(sum(E_ref_incident(3*(I-1)+1:3*I,1))))*Cells_proc(I)%Sc**3. ;  
            Cext_e_H = Cext_e_H + imag(Cells_proc(I)%Che_n*sum(E_total(3*(I-1)+1:3*I,num_emetteur+NTr))*&
                conjg(sum(E_ref_incident(3*(I-1)+1:3*I,2))))*Cells_proc(I)%Sc**3. ;            
        ENDDO 
        ! pas de 4pi ici car j'ai simplifie par le 4pi de Xi a l'interieur de la somme
        C_ext(num_emetteur) = k_0*(Cext_e_V+Cext_e_H)/2.  
        C_abs(num_emetteur) = k_0*(Cabs_e_V+Cabs_e_H)/2. 
        
        Deallocate(E_ref_incident)
    ENDDO      
    deallocate(Cells_proc);
    
    Allocate(C_ext_all(NTr),C_abs_all(NTr));
        
    Call MPI_BARRIER(MPI_COMM_WORLD ,code);
    Call MPI_ALLREDUCE(C_ext,C_ext_all,NTr,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD ,code);
    Call MPI_ALLREDUCE(C_abs,C_abs_all,NTr,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD ,code);    
    
    C_ext(1:NTr) = C_ext_all(1:NTr);
    C_abs(1:NTr) = C_abs_all(1:NTr);
    deallocate(C_ext_all,C_abs_all);

End Subroutine Compute_ExtAbsCsec_fromIntField
