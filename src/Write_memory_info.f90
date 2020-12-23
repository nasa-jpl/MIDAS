subroutine system_mem_usage(valueRSS)

    use ifport !if on intel compiler

    ! You should know that : RSS is Resident Set Size (physically resident memory - 
    ! this is currently occupying space in the machine's physical memory),
    ! and VSZ is Virtual Memory Size (address space allocated - this has addresses 
    ! allocated in the process's memory map, but there isn't necessarily any actual 
    ! memory behind it all right now).

    ! This subroutine is from the following stackoverflow discussion :
    ! http://stackoverflow.com/questions/22028571/track-memory-usage-in-fortran-90?
    ! answertab=oldest#tab-top

    implicit none

    integer, intent(out) :: valueRSS

    character(len=200):: filename=' '
    character(len=80) :: line
    character(len=8)  :: pid_char=' '
    integer :: pid
    logical :: ifxst

    valueRSS=-1    ! return negative number if not found

    !--- get process ID

    pid=getpid()
    write(pid_char,'(I8)') pid
    filename='/proc/'//trim(adjustl(pid_char))//'/status'

    !--- read system file

    inquire (file=filename,exist=ifxst)
    if (.not.ifxst) then
      write (*,*) 'system file does not exist'
      return
    endif

    open(unit=100, file=filename, action='read')
    do
      read (100,'(a)',end=120) line
      if (line(1:6).eq.'VmRSS:') then
         read (line(7:),*) valueRSS
         exit
      endif
    enddo
    120 continue
    close(100)

    return
end subroutine system_mem_usage
    
subroutine print_allocate(Nchar,allocate_str,type_str,size)

    USE Initialization
    USE common_variables
    USE MPI
    
    IMPLICIT NONE

    !IN/OUT
    Integer, INTENT(IN) :: size,Nchar 
    character(Nchar), INTENT(IN) :: allocate_str
    character(5), INTENT(IN) ::type_str ! D for Double and S for Single REAL, COMP or INTG
        
    ! local
    Real(kind=8) :: size_MB
    character(300) :: analysis_fold_name,file_name
    character(19) :: time_allocate
    character(:), allocatable :: rank_str 
    character(8)  :: date
    character(10) :: time
    character(5)  :: zone
    integer,dimension(8) :: values
    
    if ((track_memory == 1) .and. (rank .lt. Njob_max)) then 
        call date_and_time(date,time,zone,values);
        time_allocate = date(5:6)//'-'//date(7:8)//'-'//date(1:4)//'_'//time(1:2)//':'//time(3:4)//':'//time(5:6);
    
        analysis_fold_name = trim(SimOutfld_name)//Env_sep//'Analysis';
    
        if (rank .lt. 10) then 
            allocate(character(1) ::rank_str);
            Write(rank_str,'(i1)') rank;
        elseif (rank .lt. 100) then 
            allocate(character(2) ::rank_str);
            Write(rank_str,'(i2)') rank;
        elseif (rank .lt. 1000) then 
            allocate(character(3) ::rank_str);
            Write(rank_str,'(i3)') rank;
        elseif (rank .lt. 10000) then 
            allocate(character(4) ::rank_str);
            Write(rank_str,'(i14)') rank;
        endif
    
    
        file_name = trim(analysis_fold_name)//Env_sep//'TrackAllocate_j'//rank_str//'.dat'; 
    
        if (trim(allocate_str) .eq. 'Reference(t=0)') then ! reference print allocate
            Open(30+rank,File = trim(file_name)); 
        else        
            Open(30+rank,File = trim(file_name), status = 'old', position = 'append'); 
        endif
    
        if (type_str .eq. 'DCOMP') then 
            size_MB = 64.*2.*size/1e6;
        elseif (type_str .eq. 'DREal') then 
            size_MB = 64.*size/1e6;
        else
            size_MB = 32.*size/1e6       
        endif
    
    
        Write(30+rank,'(a,i16,a10,f12.3,a,a)') time_allocate, size,type_str,size_MB,'    ',allocate_str;   
        Close(30+rank);
    endif
    
    

endsubroutine print_allocate
    
    
    
    