SUBROUTINE Write_geometry_files(SimParticle,Cells,CBFM_Blocks,CBFM_Blocks_Ext,option)

    !! 9-18-2019 : option added to choose the type of file writing depending on the value of adaptive mesh and the comparison between old_Nbc and Nbc
    USE Initialization
    USE common_variables
    USE MPI
    
    Implicit NONE
    
    !IN/OUT 
    type (Particle), INTENT(IN) :: SimParticle
    type (Cell), Dimension(Nbc), INTENT(INOUT) :: Cells
    type(CBFM_Block), Dimension(Nblocks), INTENT(IN) :: CBFM_Blocks
    Integer, Dimension(Nblocks,Nbc_ext), INTENT(IN) :: CBFM_Blocks_Ext
    CHARACTER(3), INTENT(IN) :: option 

    Integer :: Type_Par,ii,jj,Nbcelsx,Nbcelsy,Nbcelsz
    Integer :: Nbcels_Dp,Nbcels_ext,Nbcels_int
    Real(kind=8) :: x0,y0,z0,Sc,ap
    
    Integer, Dimension(:,:), allocatable :: Part_in_lat
    Integer, Dimension(:), allocatable :: Nbcels_x,Nbcels_y,Nbcels_z
    
    CHARACTER(240) file_name
    CHARACTER(:), allocatable :: num_freq_str

    Type_Par = SimParticle%type_p;
    Sc = SimParticle%Sc_p; 
    ap = SimParticle%Dp/2.;
    
    ! Cells.dat file 
    if (rank == 0) then 
        ! Once Cells is reorganized, we can save it 
        ! The file Cellules.dat is used later to plot the 3D simulation scene 
        if (trim(option) .eq. 'NEW') then 
            If (EqSph ==0) then 
                file_name = trim(SimOutfld_name)//Env_sep//'Cells.dat';
            Else
                file_name = trim(SimOutfld_name)//Env_sep//'CellsES.dat';
            EndIf
        elseif (trim(option) .eq. 'UPD') then 
            if (Nfreq .eq. 1) then 
                file_name = trim(SimOutfld_name)//Env_sep//'Cells_adm.dat';
            else
                if (num_freq < 10) Then
                    Allocate(character(1)::num_freq_str)
                    Write(num_freq_str,'(i1)') num_freq
                ElseIf (num_freq < 100) Then
                    Allocate(character(2)::num_freq_str)
                    Write(num_freq_str,'(i2)') num_freq
                Else
                    Allocate(character(3)::num_freq_str)
                    Write(num_freq_str,'(i3)') num_freq
                EndIf 
                file_name = trim(SimOutfld_name)//Env_sep//'Cells_adm_f'//num_freq_str//'.dat';
            endif
        endif
    
        ! cell%num_diel to be written in Cells.dat file
        If (trim(dielcomp_option) == 'fromonlymfile') then 
            Cells(1:Nbc)%num_diel = 1;
        ElseIf ((trim(dielcomp_option) == 'fromdielcompositionfile') .OR. (trim(dielcomp_option) == 'random1')) then
            Cells(1:Nbc)%num_diel = (/1:Nbc/); ! for the other option ('fromshapefile') is read from shape file !         
        endif
        Open(14,File = trim(file_name))
        Do ii=1, Nbc
            !Write(14,'(f12.6,a,f12.6,a,f12.6,a,f12.6,a,i4,a,i4,a,f7.4,a,ES10.3,a,f7.4,a,ES10.3)') Cells(ii)%Xc,';',Cells(ii)%Yc, &
            !    ';',Cells(ii)%Zc,';',Cells(ii)%Sc,';',Cells(ii)%num_block,';',Cells(ii)%num_diel,';',&
            !    real(Cells(ii)%m_cell),' + j*',imag(Cells(ii)%m_cell),';', real(Cells(ii)%Eps_cell),' + j*',imag(Cells(ii)%Eps_cell);
        
            Write(14,'(f12.6,a,f12.6,a,f12.6,a,f12.6,a,i6,a,i8)') Cells(ii)%Xc,';',Cells(ii)%Yc, &
                ';',Cells(ii)%Zc,';',Cells(ii)%Sc,';',Cells(ii)%num_block,';',Cells(ii)%num_diel;
        EndDo
        Close(14);
    endif
    
    if ((rank == 1) .AND. (trim(option) .eq. 'NEW')) then           
        If ((Type_Par .eq. 2)) Then    ! simply copy shape file 
            file_name=ShapeFilePath; 
            if (rank ==0 ) then 
                call system('cp '//trim(file_name)//' "'//trim(SimOutfld_name)//'/"');  
            endif

        else
            !! generate the file shape.dat for DDSCat simulations
            ! Of course here, for the moment,we suppose that we have only 1 particle
            ! we need the shape file for the validation of the MoM/CBFM results in comparison to DDScat and FEKO
            Nbcels_Dp = nint(SimParticle%Dp/Sc);
            Nbcels_ext = ceiling((ap - (0.8*ap/sqrt(3.)))/Sc) ; !(ap_cheb -c/2)/Sc + 0.8 is for "security" to be sure that 
                                                                ! this cube belongs entirely to the chebyshev particle
            Nbcels_int = Nbcels_Dp - 2*Nbcels_ext;          
            Nbcels_x = [Nbcels_int,Nbcels_Dp,Nbcels_ext,Nbcels_ext,Nbcels_int,Nbcels_int,Nbcels_Dp]; 
            Nbcels_y = [Nbcels_int,Nbcels_Dp,Nbcels_Dp,Nbcels_Dp,Nbcels_ext,Nbcels_ext,Nbcels_Dp];
            Nbcels_z = [Nbcels_int,Nbcels_ext,Nbcels_int,Nbcels_int,Nbcels_int,Nbcels_int,Nbcels_ext];
            Allocate(Part_in_lat(Nbc,6));
            If ((Type_Par == 1)) Then
                x0 = (1-0.5-Nbcels_Dp/2.)*Sc
                y0 = (1-0.5-Nbcels_Dp/2.)*Sc
                z0 = (1-0.5-Nbcels_Dp/2.)*Sc       
            ElseIf ((Type_Par == 6)) Then
                x0 = SimParticle%pr_xmin
                y0 = SimParticle%pr_ymin
                z0 = SimParticle%pr_zmin      
            ElseIf ((Type_Par == 3) .or. (Type_Par == 4))Then
                x0 = anint(((1-0.5-Nbcelsx/2.)*Sc)*10**Round_Sp)/10**Round_Sp
                y0 = anint(((1-0.5-Nbcelsy/2.)*Sc)*10**Round_Sp)/10**Round_Sp
                z0 = anint(((1-0.5)*Sc)*10**Round_Sp)/10**Round_Sp
            EndIf
            
            Do ii=1,Nbc
                Part_in_lat(ii,1) = nint((Cells(ii)%Xc - x0)/Sc) 
                Part_in_lat(ii,2) = nint((Cells(ii)%Yc - y0)/Sc)
                Part_in_lat(ii,3) = nint((Cells(ii)%Zc - z0)/Sc)             
            EndDo
            Part_in_lat(1:Nbc,4:6) = 1;
            !if (rank == 0) then    
            !    if (EqSph ==0) then 
            !        file_name = trim(SimOutfld_name)//Env_sep//'shape.dat';
            !    else
            !        file_name = trim(SimOutfld_name)//Env_sep//'shapeES.dat';
            !    endif        
            !    Open(14,File = trim(file_name))
            !    Open(12,File = 'inputs/Shape_head.dat')
            !    read(12,'(a)'), fline
            !    Write(14,'(a)') fline
            !    Write(14,'(i7,a)') Nbc,' # Number of Dipoles'
            !
            !    Do ii=1, 5
            !        read(12,'(a)'), fline;
            !        Write(14,'(a)') fline;            
            !    EndDo
            !    Close(12)        
            !    Do ii=1, Nbc
            !        Write(14,'(i7,i5,i5,i5,i5,i5,i5)') ii,Part_in_lat(ii,1), &
            !        Part_in_lat(ii,2),Part_in_lat(ii,3),Part_in_lat(ii,4), &
            !        Part_in_lat(ii,5),Part_in_lat(ii,6)
            !    EndDo
            !    Close(14);  
            !EndIf
        EndIf 
    endif
    
    ! Once Cells is reorganized, we can save it 
    ! The file Cellules.dat is used later to plot the 3D simulation scene 
    if (trim(option) .eq. 'NEW') then
        if (EqSph==0) Then 
            file_name = trim(SimOutfld_name)//Env_sep//'Blocks.dat';
        Else
            file_name = trim(SimOutfld_name)//Env_sep//'BlocksES.dat';
        EndIf   
    Elseif (trim(option) .eq. 'UPD') then 
        if (Nfreq .eq. 1) then 
            file_name = trim(SimOutfld_name)//Env_sep//'Blocks_adm.dat';
        else
            file_name = trim(SimOutfld_name)//Env_sep//'Blocks_adm_f'//num_freq_str//'.dat';
        endif
    Endif
    if (rank == 2) then 
        Open(14,File = trim(file_name))
        Do ii=1, Nblocks
            Write(14,'(a,i4)') '****** Block ',ii
            Write(14,'(a,i5)') 'Nbc = ',CBFM_Blocks(ii)%Nbc_b
            Write(14,'(a,i5)') 'Nbc_ext = ',CBFM_Blocks(ii)%Nbc_ext
            Write(14,'(a)') 'Num_cells_ext ='
            Do jj=1,CBFM_Blocks(ii)%Nbc_ext
                Write(14,'(i8)') CBFM_Blocks_Ext(ii,jj)
            EndDo
            Write(14,'(a)') ' ';
        EndDo
        Close(14)
    endif    
    
    
END SUBROUTINE Write_geometry_files