SUBROUTINE Compute_Electric_Fields(SimScatterer,Cells,Transmitters,Receivers,methods_names,CBFM_Blocks,CBFM_Blocks_Ext,MPI_CBFM_Blocks)
 
    USE Initialization
    USE common_variables
    USE iso_fortran_env
    USE MPI
    Implicit none    
    
    !! IN/OUT
    type (Scatterer), INTENT(IN) :: SimScatterer
    type (Cell), Dimension(Nbc), INTENT(IN) :: Cells
    type (Dipole), Dimension(NTr), INTENT(IN) :: Transmitters
    type (Dipole), Dimension(NRx_tot), INTENT(IN) :: Receivers
    CHARACTER *leng_meth,dimension(Nber_methods), INTENT(IN):: methods_names
    type (CBFM_Block), Dimension(NBlocks), INTENT(IN):: CBFM_Blocks
    Integer, Dimension(NBlocks,Nbc_ext), INTENT(IN):: CBFM_Blocks_Ext
    Integer, Dimension(nber_procs,Nblk_proc_max+1), INTENT(IN):: MPI_CBFM_Blocks
        
    !! LOCAL
    Integer :: ii,I,K,Ind,RatioRLambda,rdm_ipws
    Real(kind=8) :: step_theta_CBFM,step_phi_CBFM,theta_dipole,phi_dipole
    Real(kind=8) :: cosdth_init,margin_th,margin_ph,a,x_l,y_l,z_l
    Real(kind=8) :: th_init,th_end,ph_init,ph_end
    COMPLEX(real64), Dimension(:,:),allocatable :: E_total,S_total
    COMPLEX(real64), Dimension(:),allocatable :: C_ext,C_abs
        
    ! Time performances
    character(8)  :: date_init, date_final
    character(10) :: time_init, time_final
    character(5)  :: zone_init, zone_final    
    Integer,dimension(8) :: values_init, values_final
    Integer, dimension(4):: Comp_time    
    
    INTERFACE 
        SUBROUTINE Compute_EFields_CBFME(Cells,CBFM_Blocks,CBFM_Blocks_Ext,MPI_CBFM_Blocks,Transmitters,E_total)    
 
            USE Initialization
            USE common_variables
            USE lapack95
            USE DiverseUtil
            USE MPI

            IMPLICIT NONE
    
            !! IN/OUT ******************************************************************
            type (Cell), Dimension(Nbc), INTENT(IN) :: Cells
            type (CBFM_Block), Dimension(Nblocks), INTENT(IN) :: CBFM_Blocks
            Integer, Dimension(Nblocks,Nbc_ext), INTENT(IN) :: CBFM_Blocks_Ext
            Integer, Dimension(nber_procs,Nblk_proc_max+1), INTENT(IN):: MPI_CBFM_Blocks
            type (Dipole), Dimension(NTr), INTENT(IN) ::Transmitters
            COMPLEX(real64), Dimension(:,:), allocatable, INTENT(OUT):: E_total
        END SUBROUTINE Compute_EFields_CBFME 
        
        SUBROUTINE Compute_EFields_CBFME_ACA(Cells,CBFM_Blocks,CBFM_Blocks_Ext,MPI_CBFM_Blocks,Transmitters,E_total)  
            USE Initialization
            USE common_variables
            USE lapack95
            USE DiverseUtil
            USE MPI

            IMPLICIT NONE
    
            !! IN/OUT ******************************************************************
            type (Cell), Dimension(Nbc), INTENT(IN) :: Cells
            type (CBFM_Block), Dimension(Nblocks), INTENT(IN) :: CBFM_Blocks
            Integer, Dimension(Nblocks,Nbc_ext), INTENT(IN) :: CBFM_Blocks_Ext
            !Integer, Dimension(2+Nccp_max,NcalBlks), INTENT(IN) :: cp_CBFM_Blocks
            Integer, Dimension(nber_procs,Nblk_proc_max+1), INTENT(IN):: MPI_CBFM_Blocks
            type (Dipole), Dimension(NTr), INTENT(IN) ::Transmitters
            COMPLEX(real64), Dimension(:,:), allocatable, INTENT(OUT):: E_total
        END SUBROUTINE Compute_EFields_CBFME_ACA        
    
            
    END INTERFACE    
    
        
    !**********************************************
    !***************  CBFM-E  *********************
    !**********************************************    
    if (CBFM .NE. 0) Then
        !! Compute Fields inside the scatterer *********************************
        Comp_time = 0
        call date_and_time(date_init,time_init,zone_init,values_init)
        !if (SR_Zc == 0) Then  
          if (Use_ACA==0) then
              Call Compute_EFields_CBFME(Cells,CBFM_Blocks,CBFM_Blocks_Ext,MPI_CBFM_Blocks,Transmitters,E_total); !,cp_CBFM_Blocks)
          else
              Call Compute_EFields_CBFME_ACA(Cells,CBFM_Blocks,CBFM_Blocks_Ext,MPI_CBFM_Blocks,Transmitters,E_total);  
          EndIf
        !Else   SR_Zc = 1 not implemented for now !                   
        !Endif
        
        call date_and_time(date_final,time_final,zone_final,values_final)
        call Calcul_time_spent(values_init,values_final,Comp_time)
        if (rank == 0) then 
            Write (*, '(a)') '';
            Write (*, '(a)') 'The total time to compute the internal electric field with 1L CBFM-E';
            Write (*, '(a,i2,a,i2,a,i2,a,i2,a)')'is ', Comp_time(1),'j',Comp_time(2),'h',Comp_time(3),&
                'min', Comp_time(4),'sec'
            Write (10, '(a,i2,a,i2,a,i2,a,i2,a)') 'The total time to compute the internal electric field &
                &with 1L CBFM-E is ',&
                Comp_time(1),'j',Comp_time(2),'h',Comp_time(3),'min', Comp_time(4),'sec'
        endif
                       
        ! Compute scattered fields ************************************************** 
        Comp_time = 0; call date_and_time(date_init,time_init,zone_init,values_init); 
        
        Allocate(S_total(NRx_tot,4*NTr),C_ext(NTr),C_abs((NTr)));
        Call Compute_Scattering_Matrices('CBFM-E  ',Cells,E_total,CBFM_Blocks,MPI_CBFM_Blocks,Transmitters,Receivers,S_total);
        Call Compute_ExtAbsCsec_fromIntField('CBFM-E  ',Cells,E_total,Transmitters,C_ext,C_abs);
        
        if (rank == 0) then 
            call date_and_time(date_final,time_final,zone_final,values_final)
            call Calcul_time_spent(values_init,values_final,Comp_time)
            Write (*,*) ''
            Write (*, '(a,i2,a,i2,a,i2,a,i2,a)') 'The calculation time for the scattered field is ',&
            Comp_time(1),'j',Comp_time(2),'h',Comp_time(3),'min', Comp_time(4),'sec'
            Write (*,*) ''; 
            Write (10, '(a,i2,a,i2,a,i2,a,i2,a)') 'The calculation time for the scattered field is ',&
            Comp_time(1),'j',Comp_time(2),'h',Comp_time(3),'min', Comp_time(4),'sec'
            Write (10,*) ''
            Write (10,*) ''  
        endif
        if ((NumIntType_t .eq. 'aq') .OR. (NumIntType_t .eq. 'gl') .OR. (NumIntType_t .eq. 'tr') .OR. (NumIntType_t .eq. 'sm')) Then
          call Compute_Scattering_Quantities_1('CBFM-E  ',SimScatterer,Transmitters,Receivers,S_total,C_ext,C_abs)
        else
          call Compute_Scattering_Quantities_2('CBFM-E  ',SimScatterer,Transmitters,Receivers,S_total,C_ext,C_abs)
        endif
        Deallocate(E_total,S_total,C_ext,C_abs)  
    Endif
    
    
    !**********************************************
    !*****************  MoM  **********************
    !**********************************************    
    if (MoM .NE. 0) Then 
    !if ((MoM .NE. 0) .and. (rank .EQ. 0)) Then ! if we want to use/test with the single Task MoM   
        !! Compute Fields inside the scatterer *********************************
        Comp_time = 0
        call date_and_time(date_init,time_init,zone_init,values_init)
        
        Allocate(S_total(NRx_tot,4*NTr),C_ext(NTr),C_abs((NTr)));
        Call Compute_EFields_MoM(Cells,Transmitters,Receivers,S_total,C_ext,C_abs);  ! MPI MoM
        !Call Compute_EFields_ST_MoM(Cells,Transmitters,Receivers,S_total,C_ext,C_abs) ! MPI single-task MoM (equivalent to OpenMP MoM)
        
        call date_and_time(date_final,time_final,zone_final,values_final)
        call Calcul_time_spent(values_init,values_final,Comp_time)
        
        if (rank == 0) then 
            Write (*, '(a)') '';
            Write (*, '(a)') 'The total time to compute the internal & Scattered fields with MoM';
            Write (*, '(a,i2,a,i2,a,i2,a,i2,a)')'is ', Comp_time(1),'j',Comp_time(2),'h',Comp_time(3),&
                'min', Comp_time(4),'sec'
            Write (10, '(a,i2,a,i2,a,i2,a,i2,a)') 'The total time to compute the internal & Scattered fields &
                &with MoM is ',&
                Comp_time(1),'j',Comp_time(2),'h',Comp_time(3),'min', Comp_time(4),'sec'
        endif
                       
        if ((NumIntType_t .eq. 'aq') .OR. (NumIntType_t .eq. 'gl') .OR. (NumIntType_t .eq. 'tr') .OR. (NumIntType_t .eq. 'sm')) Then
          call Compute_Scattering_Quantities_1('CBFM-E  ',SimScatterer,Transmitters,Receivers,S_total,C_ext,C_abs)
        else
          call Compute_Scattering_Quantities_2('CBFM-E  ',SimScatterer,Transmitters,Receivers,S_total,C_ext,C_abs)
        endif
        Deallocate(S_total,C_ext,C_abs)  
    Endif
    
END SUBROUTINE Compute_Electric_Fields
    
    