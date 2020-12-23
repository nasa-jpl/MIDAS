SUBROUTINE Particle_DielComposition(m_lambdas,SimParticle,Cells)
    
    USE Initialization
    USE common_variables
    USE iso_fortran_env
        
    Implicit NONE
    
    ! IN/OUT 
    Complex, Dimension(Ndiel,Nfreq), INTENT(IN) :: m_lambdas
    type (Particle), INTENT(INOUT) :: SimParticle
    type (Cell), Dimension(Nbc), INTENT(INOUT):: Cells
    
    ! local 
    Integer :: ii,ii_diel
    real(kind=8) :: mrp,mip,rp,ip
       
    
    !! It is in this subroutine the dielectric constant is selected for each cell ! Cell.m_cell and Cell.Eps_cell
    !! later in UpdateCellsParameters, depending on the dielectric properties decided here and on the wavelength 
    !! in consideration, lambda_s per cell will be determined ! Then, if required by the user, Sc for each cell can be adapted to its lambda_s
    
    !! ATTENTION : Cells(ii)%num_diel is set here as it is used just after to retrieve the good m value for each cell from m_lambdas
    ! cell%num_diel (in case if was not already initialized in Write_geometry_files.f90)
    ! if 'fromshapefile' : num_diel is read from shapefile and if 'random2' num_diel is initialized in get_diel_values_lambdas (below)
    If (trim(dielcomp_option) == 'fromonlymfile') then 
        Cells(1:Nbc)%num_diel = 1;
    ElseIf ((trim(dielcomp_option) == 'fromdielcompositionfile') .OR. (trim(dielcomp_option) == 'random1')) then
        Cells(1:Nbc)%num_diel = (/1:Nbc/); ! for the other option ('fromshapefile') is read from shape file ! 
    endif
    
    !! Since Cell.m_cell and Cell.Eps_cell can vary with the frequency this subroutine should be called inside the loop on lambda!
    if ((trim(dielcomp_option) == 'fromshapefile') .OR. (trim(dielcomp_option) == 'fromdielcompositionfile') .OR. (trim(dielcomp_option) == 'fromonlymfile')) then
        Do ii= 1,Nbc
            ii_diel = Cells(ii)%num_diel;
            Cells(ii)%m_cell = m_lambdas(ii_diel,num_freq);
            
            mrp = real(m_lambdas(ii_diel,num_freq))
            mip = imag(m_lambdas(ii_diel,num_freq)) ; 
            rp = mrp**2-mip**2;
            ip = 2*mrp*mip;
            Cells(ii)%Eps_cell = rp+J*ip  
            Cells(ii)%lambda_cell = Lambda_w/sqrt(rp)
            ! Update particle Dlam
            Cells(ii)%Dlamb_cell = Cells(ii)%lambda_cell/Cells(ii)%Sc;
        EndDo        
    else ! I think it will be the same as fromshapefile and fromdielcompositionfile as we will use also Cells(ii)%num_diel and m_lambdas properly generated for
        ! the two random options in get_diel_values_lambdas
        Write(*,*) 'This dielcomp_option is still under development ! Thank you for your patience!';
        stop 1;
    EndIf     
END SUBROUTINE Particle_DielComposition
    
SUBROUTINE get_diel_values_lambdas(m_file_name,m_lambdas)
    
    USE Initialization
    USE common_variables
    USE iso_fortran_env
    USE strings
    
    Implicit NONE
    
    ! IN/OUT 
    character(250), dimension(Ndiel), INTENT(IN) :: m_file_name
    Complex, Dimension(:,:), allocatable, INTENT(OUT) :: m_lambdas
    
    ! out : m_lambdas - Ndiel x Nfreq  ==> Ndiel is the number of mfiles for 'fromshapefile' OR 'fromonlymfile' OR 'random2'
    !                                      Ndiel is 2 (min & max values) for 'random1'
    !                                      Ndiel = Nbfilelines for 'fromdielcompositionfile'. Since we don't have any information yet about Nbc
                                                                  ! we will read all the available m values and compare later with Nbc.  
        
    ! local 
    Integer :: ii,jj,dd,nvals, nitems,n,nargs,Nfreq_dielfile
    Real(kind=8) :: lmbd, mrp, mip
    character(250) :: fname;
    Complex, Dimension(:), allocatable :: ms_mfile
    Real(kind=8), Dimension(:), allocatable :: lamb_ii
    Character(200) cc,frmt
    character(100),dimension(3) :: args1
    character(100),dimension(2) :: args2
        
    if (trim(dielcomp_option) == 'fromdielcompositionfile') then
        !Allocate(m_lambdas(Ndiel,Nfreq));
        !m_lambdas = 1.; ! Actually m_lambdas is useless for this option at this level as we will not need any information about m or lambda inside the scatterer
                        ! before having all the discretization info from the shape file. so we simply allocate 1x Nfreq and initialize it to 1 just to be 
                        ! consistent with regard to the other options
        
        Open(11,File = trim('inputs/dielcomposition.dat'));
        read(11,'(a)'),cc; read(11,'(a)'),cc;
        call parse(cc,'=',args1,nargs);
        read(args1(3),'(i)'), Nfreq_dielfile
        call parse(args1(2),';',args2,nargs);
        read(args2(1),'(i)'), Ndiel
        
        Allocate(m_lambdas(Ndiel,Nfreq));
        Allocate(lamb_ii(2*Nfreq_dielfile))
        read(11,*),cc
        if (Nfreq .le. Nfreq_dielfile) then 
            frmt = '(i8,';
            Do ii=1, Nfreq_dielfile-1
                frmt = trim(frmt)//'f9.4,f9.4,'
            EndDo 
            frmt = trim(frmt)//'f9.4,f9.4)';
                        
            Do ii= 1,Ndiel
                read(11,trim(frmt)) dd,lamb_ii(:)
                do jj = 1,Nfreq
                    m_lambdas(ii,jj) = lamb_ii(2*jj-1)+J*lamb_ii(2*jj);                    
                EndDo                
            EndDo   
        elseif (Nfreq_dielfile .eq. 1) then ! we use this single value for all calculated frequencies 
            Do ii= 1,Ndiel
                read(11,'(i8,f9.4,f9.4)') jj,mrp,mip
                m_lambdas(ii,:) = mrp+J*mip;                                   
            EndDo
        else ! too complicated to decide here -> error 
            stop 1;
            Write(*,'(a,a)') 'Nfreq_dielfile < Nfreq and .ne. to 1 ! Please use another dielcompositionfile !! '            
        endif
        
    else 
        !! Ndiel >=1 if (trim(dielcomp_option) == 'fromshapefile') .OR. (trim(dielcomp_option) == 'random2') : we read from Ndiel different m files 
        !! Ndiel = 2 if (trim(dielcomp_option) == 'random1') : we read m min and max values from m files, they can or not depend on frequency 
        !! Ndiel = 1 if (trim(dielcomp_option) == 'fromonlymfile') : ! we read m per frequency from the only available m file
        Allocate(m_lambdas(Ndiel,Nfreq));
        Allocate(ms_mfile(Nfreq));
        ! we read m values from m files, they can or not depend on frequency
        Do dd=1,Ndiel
            ms_mfile = 0;
            nvals = 0;
            ! fill out m_Wavesle(ii,:)
            fname = m_file_name(dd);
            Open(11,File = trim(fname)) 
                
            Do ii=1,3; read (11,*); enddo !3 first info lines in an m file  
            Do
                read (11,*, end=10),lmbd,mrp,mip
                ms_mfile(nvals+1) = mrp+J*mip
                nvals = nvals + 1;
            Enddo
10          close(11);
            if (nvals .eq. Nfreq) then
                m_lambdas(dd,1:Nfreq) = ms_mfile(1:Nfreq);
            elseif (nvals .eq. 1) then 
                m_lambdas(dd,1:Nfreq) = ms_mfile(1);
            else
                stop 1;
                Write(*,'(a,a)') 'Error while reading m values from ',fname
            endif
        endDo
        deallocate(ms_mfile);               
    endif       
    ENDSUBROUTINE get_diel_values_lambdas
    
integer function nitems(line)
    
    ! number of space-separated items in a line
    ! function nitems from https://www.tek-tips.com/viewthread.cfm?qid=1688013

    character line*(*)    
    logical back
    integer length
    
    back = .true.        
    length = len_trim(line)    
    k = index(line(1:length), ' ', back)
    if (k == 0) then
        nitems = 0
        return
    end if    
    
    nitems = 1
    do 
        ! starting with the right most blank space, 
        ! look for the next non-space character down
        ! indicating there is another item in the line
        do
            if (k <= 0) exit
            
            if (line(k:k) == ' ') then
                k = k - 1
                cycle
            else
                nitems = nitems + 1
                exit
            end if
            
        end do
        
        ! once a non-space character is found,
        ! skip all adjacent non-space character
        do
            if ( k<=0 ) exit
            
            if (line(k:k) /= ' ') then
                k = k - 1
                cycle
            end if
            
            exit
            
        end do
        
        if (k <= 0) exit
            
    end do
end function nitems     
 