SUBROUTINE Discretization(SimScatterer,Cells,Ncells_SphDomains)
    
    USE Initialization
    USE common_variables
    USE iso_fortran_env
    USE MPI
    
    Implicit NONE
    
    ! IN/OUT 
    type (Scatterer), INTENT(INOUT) :: SimScatterer
    type (Cell), Dimension(:), allocatable, INTENT(OUT):: Cells
    Integer, Dimension(7), INTENT(OUT), OPTIONAL:: Ncells_SphDomains

    ! Local
    Integer :: ii,idip,Nbcels_Dp,Qd,Ql,Dp_tmp,Dp_int,Lp_tmp,Lp_int
    Integer :: I,Nbc_init,Comp_cel,Ix,Iy,Iz,Nbc_p
    Integer :: dd,Ncelp_m,nb,Type_Par,ap_count,info_count
    Integer :: tot_count,pr_NBcels,Nbc_x,Nbc_y,Nbc_z
    Integer :: Nbcelsx, Nbcelsy, Nbcelsz, hcy_rcy, Nbcelsz_sphe
    Integer :: ncp,Nbc_bef,Nbc_aft,chgmt,Nbc_p_dom
    Integer :: cheb_l,Nbcels_int,Nbcels_ext
    Integer :: Nbcels_d, Nbcels_h
    Real(kind=8) :: Sc_max,Sc_tmp,Sc,Dp,ap,pas,lambda_p
    Real(kind=8) :: X,Y,Z,Rc,A,B,C,alphaz,alphay,hcy,rcy,z_ref
    Real(kind=8) :: x0,y0,z0
    Real(kind=8) :: cheb_eps,hypotxy,theta_rd,phi_rd,Rlim,ap_cheb,ap_int,Dp_cheb 
    Real(kind=8) :: h_cyl, d_cyl, Lp
    
    Integer, Dimension(:,:), allocatable :: Part_in_lat
    Integer, Dimension(:), allocatable :: Nbcels_x,Nbcels_y,Nbcels_z
    type (Cell), Dimension(:), allocatable :: TmpCells,TmpCellsPerm
    Real(kind=8), Dimension(3) :: a1,a2
    Real(kind=8), dimension(:,:), allocatable :: positions
    Real(kind=8), dimension(:), allocatable :: start_x, start_y, start_z
    COMPLEX(real64) :: Ac,Ai, Bi, Ci
    
    CHARACTER(:), allocatable::info_p_fl
    CHARACTER tmp_str,ch1,ch2
    CHARACTER(3) st_ap
    CHARACTER(9) wp_info
    CHARACTER(14) wp_ap
    CHARACTER(23) wp_tot
    CHARACTER(100) fline
    CHARACTER(240) file_name
    
      
    INTERFACE
        SUBROUTINE Read_ShapeFile(info_p,pr_NBcels,pr_lattice)
            USE Initialization
            USE common_variables
            !USE f95_precision
            Implicit NONE
            !IN/OUT            
            Character, INTENT(IN) :: info_p
            Integer, INTENT(OUT) :: pr_NBcels
            Integer, Dimension(:,:), allocatable, INTENT(OUT):: pr_lattice
            
        END SUBROUTINE Read_ShapeFile    
    END INTERFACE
    
    ! REMEMBER Type_Par indidates the type of geometry of the simulated scatterer :
    ! 1 : Simple Sphere : just 1 
    ! 2 : Read pristine or aggregate complex geometry from Shape file. Example :" 2a-0006 or 2p-08
    ! 3 : Cylinder
    ! 4 : Cylinder associated with two semi-spheres 
    ! 5
    ! 6 : Chebyshev particle. Example : 6c09-0.20 

    ! general initialization to cover both Spherical and Arbitray shapes
    Ncells_SphDomains = 0;
    
    ! STEP 1 : Determine Sc (Size of cell) for each scatterer depending on its size and lambda
    ! This will enable us to define an initial Nbc to allocate the array Cells 
    !***********************************************************************************************
    Nbc_init = 0
    
    Type_Par = SimScatterer%type_s
    if (Type_Par .eq. 2) then 
        Allocate(character(1) :: info_p_fl);       
    else
        Allocate(character(9) :: info_p_fl);
    endif
    info_p_fl = trim(SimScatterer%info_s);
    if (Type_Par == 6) Then !6c02-00.20
        Read(info_p_fl,'(a,i2,a,f5.2)') ch1,cheb_l,ch2,cheb_eps;
    endif    
    
    If (Type_par == 2) then 
        Nbc_init = Nbc_init + NBc_max_alloc;
    Else
        If ((Type_Par == 1) .or. (Type_Par == 5)) Then ! sphere 
            Dp = anint(SimScatterer%dm*10**Round_S)/10**Round_S;
        Elseif (Type_Par == 6) Then !6c02-00.20 ! chebychev particle
            Read(info_p_fl,'(a,i2,a,f4.2)') ch1,cheb_l,ch2,cheb_eps;
            ap_cheb = (1.+cheb_eps)*(SimScatterer%a);  
            Dp = anint(2*ap_cheb*10**Round_S)/10**Round_S;
        Elseif (Type_Par == 3) Then !cylinder
            Dp = SimScatterer%dm;
            Lp = SimScatterer%dx; 
        EndIf        
            
        if (Dlambda .gt. 1) then
          Sc_max = anint(10**Round_S*SimScatterer%lambda_min/Dlambda)&
            /(10**Round_S);            
          Nbcels_Dp = Dp/Sc_max;
          Sc_tmp = anint(Sc_max*10**(Round_S-1))/10**(Round_S-1); 
        else
          Sc_tmp = SimScatterer%Sc;
        endif               
         
        if (Type_Par == 3) Then !cylinder : here Sc should divide both r and l 
            Qd = nint(Dp*10.**Round_S)/nint(Sc_tmp*10.**Round_S)
            Dp_tmp = Qd*nint(Sc_tmp*10.**Round_S)
            Dp_int = (nint(Dp*10.**(Round_D+1.))*10.**Round_S)/10.**(Round_D+1.)
            
            Ql = nint(Lp*10.**Round_S)/nint(Sc_tmp*10.**Round_S)
            Lp_tmp = Ql*nint(Sc_tmp*10.**Round_S)
            Lp_int = (nint(Lp*10.**(Round_D+1.))*10.**Round_S)/10.**(Round_D+1.)
        
            chgmt = 0;
            Do while ((Dp_tmp .NE. Dp_int) .OR. (Lp_tmp .NE. Lp_int))
                pas = 5*(1./10**Round_S)
                Sc_tmp = Sc_tmp - pas 
                Qd = nint(Dp*10**Round_S)/nint(Sc_tmp*10**Round_S)
                Dp_tmp = Qd*nint(Sc_tmp*10**Round_S);
                
                Ql = nint(Lp*10**Round_S)/nint(Sc_tmp*10**Round_S)
                Lp_tmp = Ql*nint(Sc_tmp*10**Round_S);
                chgmt = 1;            
            EndDo 
                
            if ((Dlambda .eq. 1)  .and.  (chgmt .eq. 1)) then
                if (rank .eq. 0) then
                Write(*,'(a,f6.2,a)') 'Attention : Sc was automatically changed to Sc = ',Sc_tmp*1e6,' um' 
                endif
            endif
        
            SimScatterer%Sc = Sc_tmp
            SimScatterer%Dlamb = anint((SimScatterer%lambda_min)/SimScatterer%Sc)
            Qd = nint(Dp*10**Round_S)/nint(SimScatterer%Sc*10**Round_S)
            Ql = nint(Lp*10**Round_S)/nint(SimScatterer%Sc*10**Round_S)
        else
            Qd = nint(Dp*10**Round_S)/nint(Sc_tmp*10**Round_S)
            Dp_tmp = Qd*nint(Sc_tmp*10**Round_S)
            Dp_int = (nint(Dp*10**(Round_D+1))*10**Round_S)/10**(Round_D+1)
            
            chgmt = 0;
            Do while (Dp_tmp .NE. Dp_int)
                pas = 5*(1./10**Round_S)
                Sc_tmp = Sc_tmp - pas 
                Qd = nint(Dp*10**Round_S)/nint(Sc_tmp*10**Round_S)
                Dp_tmp = Qd*nint(Sc_tmp*10**Round_S);
                chgmt = 1;            
            EndDo 
                
            if ((Dlambda .eq. 1)  .and.  (chgmt .eq. 1)) then
                Write(*,'(a,f6.2,a)') 'Attention : Sc was automatically changed to Sc = ',Sc_tmp*1e6,' um' 
            endif
        
            SimScatterer%Sc = Sc_tmp
            SimScatterer%Dlamb = anint((SimScatterer%lambda_min)/SimScatterer%Sc)
            Qd = nint(Dp*10**Round_S)/nint(SimScatterer%Sc*10**Round_S)            
        EndIf
               
        if ((Type_Par == 1) .OR. (Type_Par == 6)) Then 
            Nbc_init = Nbc_init + Qd**3    
        elseIf (Type_Par == 3) Then 
            Nbc_init = Qd**2*Ql
        elseIf (Type_Par == 4) Then
            ! here we add the cells contained in the two half-sphere
            Nbc_init = Nbc_init + 4*Qd**2*(hcy_rcy*Qd) + 8*Qd**3
        EndIf
    EndIf          
    !***********************************************************************************************
    
    Allocate (TmpCells(Nbc_init))     
    ! Discretization
    !***********************************************************************************************
    Comp_cel = 0
    Nbc_p = 0;
    Type_Par = SimScatterer%type_s ! for Now, Type =1 if sphere; =2 if the scatterer is to be read 
    ! from a shape file (data set from DDSCAT)
    
        
    ! The descritization of the scatterer depends on its type 
    !! NOTE THAT THE WAY I GENERATE THE SPHERE HERE IS NOT SIMILAR TO WHAT DDSCAT DOES ! 
    ! IF I CONVERT TARGET.OUT TO SHAPE FILE WE ARE CLOSER TO MIE! SOOOOOOO A EXPLORER !!!!
    If ((Type_Par == 1) .or. (Type_Par == 5))Then ! SIMPLE SPHERE 
        Sc = SimScatterer%Sc
        ap = SimScatterer%dm/2.; 
        Nbcels_Dp = nint(SimScatterer%dm/Sc); 
        
        Nbcels_ext = ceiling((ap - (0.8*ap/sqrt(3.)))/Sc) ; !(ap_cheb -c/2)/Sc + 0.8 is for "security" to be sure that 
                                                                ! this cube belongs entirely to the chebyshev particle  
        Nbcels_int = Nbcels_Dp - 2*Nbcels_ext;
        ap_int = (Nbcels_int*Sc)/2.;
                
        Allocate(start_x(7),start_y(7),start_z(7));
        Allocate(Nbcels_x(7),Nbcels_y(7),Nbcels_z(7))
        start_x = [-ap_int,-ap,-ap,ap_int,-ap_int,-ap_int,-ap]; 
        Nbcels_x = [Nbcels_int,Nbcels_Dp,Nbcels_ext,Nbcels_ext,Nbcels_int,Nbcels_int,Nbcels_Dp]; 
        start_y = [-ap_int,-ap,-ap,-ap,-ap,ap_int,-ap]; 
        Nbcels_y = [Nbcels_int,Nbcels_Dp,Nbcels_Dp,Nbcels_Dp,Nbcels_ext,Nbcels_ext,Nbcels_Dp];
        start_z = [-ap_int,-ap,-ap_int,-ap_int,-ap_int,-ap_int,ap_int]; 
        Nbcels_z = [Nbcels_int,Nbcels_ext,Nbcels_int,Nbcels_int,Nbcels_int,Nbcels_int,Nbcels_ext];
        
        ! CUBE INSCRIT DANS LA SPHERE
        DO Iz=1,Nbcels_z(1)
            DO Iy=1,Nbcels_y(1)
                DO Ix=1,Nbcels_x(1)                               
                    
                    Comp_cel = Comp_cel + 1
                    Nbc_p = Nbc_p + 1
                    
                    TmpCells(Comp_cel)%num_cell = Comp_cel
                    TmpCells(Comp_cel)%Xc = start_x(1) + (Ix-0.5)*Sc
                    TmpCells(Comp_cel)%Yc = start_y(1) + (Iy-0.5)*Sc
                    TmpCells(Comp_cel)%Zc = start_z(1) + (Iz-0.5)*Sc 
                    TmpCells(Comp_cel)%Sc = Sc                   
                Enddo
            Enddo
        Enddo
        Ncells_SphDomains(1) = Nbc_p;
        ! Les blocs entourant le bloc inscrit dans la sphere
        Do ii=2,7 
            Nbc_p_dom = 0; 
            DO Iz=1,Nbcels_z(ii)
                DO Iy=1,Nbcels_y(ii)
                    DO Ix=1,Nbcels_x(ii)              
                    
                        x = start_x(ii) + (Ix-0.5)*Sc
                        y = start_y(ii) + (Iy-0.5)*Sc
                        z = start_z(ii) + (Iz-0.5)*Sc                    
                    
                        Rc = sqrt(x**2.+y**2.+z**2.);
                        If (Rc < ap) Then                                                                                    
                            Comp_cel = Comp_cel + 1
                            Nbc_p = Nbc_p + 1;
                            Nbc_p_dom = Nbc_p_dom + 1;
                    
                            TmpCells(Comp_cel)%num_cell = Comp_cel
                            TmpCells(Comp_cel)%Xc = x
                            TmpCells(Comp_cel)%Yc = y
                            TmpCells(Comp_cel)%Zc = z
                            TmpCells(Comp_cel)%Sc = Sc                                                     
                        EndIf
                    Enddo
                Enddo
            Enddo 
            Ncells_SphDomains(ii) = Nbc_p_dom;
        EndDo      
        Nbc = Nbc_p               
    ElseIf (Type_Par == 2 ) Then ! COMPLEX SHAPE PARTICLE (DDSCAT : read the shape file)
        
        !! A - According to the simulation data, find out the name of the shape file to read 
        if ((info_p_fl == 'a') .or. (info_p_fl == 's')) Then
            wp_info = '(a6,a6,a1'  
            info_count = 14;
        Else if (info_p_fl == 'p') Then
            wp_info = '(a6,a4,a1'
            info_count = 12;
        EndIf
        
        ap = SimScatterer%a*10**3;
        st_ap = '000'; ap_count = 25;
        If (ap .lt. 1) Then                
            wp_ap = ',a3,F10.6,a12)'
        ElseIf (ap .lt. 10) Then
            wp_ap = ',a2,F11.6,a12)'
        Else
            wp_ap = ',a1,F12.6,a12)'
        EndIf
        
        !wp_tot = wp_info//wp_ap; tot_count = info_count + ap_count;
        !Allocate(character(tot_count) :: file_name)
        !write(file_name,wp_tot) 'shape\',trim(SimScatterer%info_s),'\',st_ap,(SimScatterer%a*10**6),'um\shape.dat' 
        !file_name=ShapeFilePath; 
        !if (rank ==0 ) then 
        !    call system('cp '//trim(file_name)//' "'//trim(SimOutfld_name)//'/"');  
        !endif
                
        !! B - Read from the the shape file the number and positions of the cells/dipoles in the lattice
        Call Read_ShapeFile(info_p_fl,pr_NBcels,Part_in_lat);
        
        !! C - Now convert the previous information to a set of cells (with x, y and z positions)
        ! First of all, check if the loaded pr_NBcels is greater than the current Nbc_init, if yes 
        ! Copy TmpCells in greater array and deallocate it 
        If (pr_NBcels .ge. Nbc_init) Then 
            Nbc_init = pr_NBcels + Nbc_init ;
            Allocate(TmpCellsPerm(Nbc_init));TmpCellsPerm(1:Comp_cel) = TmpCells(1:Comp_cel);
            Deallocate(TmpCells); Allocate(TmpCells(Nbc_init));
            TmpCells(1:Comp_cel)=TmpCellsPerm(1:Comp_cel); Deallocate(TmpCellsPerm);                
        EndIf
        
        ! Start with updating the properties of the current scatterer according to the 
        ! information read from shape.dat (new Sc, Nbc_p and angle of rotation)
        Sc = ((4*Pi)/(3*pr_NBcels))**(1./3.)*SimScatterer%a;        
        SimScatterer%Sc = Sc;
        
        SimScatterer%Dlamb = anint((SimScatterer%lambda_min)/Sc)
        Nbc = pr_NBcels   
        
        ! Now, calculate the X, Y and Z coordiantes and other EM properties for each cell 
        ! we consider that the origin of the new coordinate system a1,a2 is the origin of the latice (0,0,0)
        
        Do idip = 1,pr_NBcels
            Comp_cel = Comp_cel + 1
            
            TmpCells(Comp_cel)%num_cell = Comp_cel
            TmpCells(Comp_cel)%Xc = Sc*Part_in_lat(idip,1) + Sc/2.
            TmpCells(Comp_cel)%Yc = Sc*Part_in_lat(idip,2) + Sc/2.
            TmpCells(Comp_cel)%Zc = Sc*Part_in_lat(idip,3) + Sc/2.
            
            TmpCells(Comp_cel)%Sc = Sc     
            TmpCells(Comp_cel)%num_diel = Part_in_lat(idip,4); ! we only consider isotropic scatterers for the moment         
        EndDo            
        
        ! we recall that the scatterer attributes pb_xmin, pb_xmax,pb_ymin, pb_ymax 
        ! pb_zmin and pb_zmax of the current scatterer will be assigned in the subroutine Division_blocks  
        
    ElseIf (Type_Par == 3) Then 
        ! the scatterer is a simple Cylinder of height h (mm) and radius r (mm) read from the simulation data input file
        Sc = SimScatterer%Sc;
        d_cyl = SimScatterer%dy; !or SimScatterer%dz
        
        DO Ix=1,Ql
            DO Iy=1,Qd
                DO Iz=1,Qd              
                    x = (Ix-0.5)*Sc
                    y = -d_cyl/2.+(Iy-0.5)*Sc
                    z = -d_cyl/2.+(Iz-0.5)*Sc
                
                    Rc = sqrt(y**2.+z**2.);
                    If (Rc < d_cyl/2.) Then                     
                                                                    
                        Comp_cel = Comp_cel + 1
                        Nbc_p = Nbc_p + 1
                    
                        TmpCells(Comp_cel)%num_cell = Comp_cel
                        TmpCells(Comp_cel)%Xc = x
                        TmpCells(Comp_cel)%Yc = y
                        TmpCells(Comp_cel)%Zc = z
                        TmpCells(Comp_cel)%Sc = Sc
                    EndIf
                Enddo
            Enddo
        Enddo        
    ElseIf (Type_Par == 5) Then ! The scatterer is a sphere attached to N cylinders 
    ElseIf (Type_Par == 6)Then ! Chebyshev particle 
        Sc = SimScatterer%Sc
        ap_cheb = (1.+cheb_eps)*SimScatterer%dm/2.;
        Dp_cheb = anint(2*ap_cheb*10**Round_S)/10**Round_S;
        Nbcels_Dp = nint(Dp_cheb/Sc); 
        
        ! put back ap to r0 
        ap = SimScatterer%dm/2.; 
        Nbcels_ext = ceiling((ap_cheb - (0.8*ap/sqrt(3.)))/Sc) ; !(ap_cheb -c/2)/Sc + 0.8 is for "security" to be sure that 
                                                                ! this cube belongs entirely to the chebyshev particle  
        Nbcels_int = Nbcels_Dp - 2*Nbcels_ext;
        ap_int = (Nbcels_int*Sc)/2.; 
        
        Allocate(start_x(7),start_y(7),start_z(7));
        Allocate(Nbcels_x(7),Nbcels_y(7),Nbcels_z(7))
        start_x = [-ap_int,-ap_cheb,-ap_cheb,ap_int,-ap_int,-ap_int,-ap_cheb]; 
        Nbcels_x = [Nbcels_int,Nbcels_Dp,Nbcels_ext,Nbcels_ext,Nbcels_int,Nbcels_int,Nbcels_Dp]; 
        start_y = [-ap_int,-ap_cheb,-ap_cheb,-ap_cheb,-ap_cheb,ap_int,-ap_cheb]; 
        Nbcels_y = [Nbcels_int,Nbcels_Dp,Nbcels_Dp,Nbcels_Dp,Nbcels_ext,Nbcels_ext,Nbcels_Dp];
        start_z = [-ap_int,-ap_cheb,-ap_int,-ap_int,-ap_int,-ap_int,ap_int]; 
        Nbcels_z = [Nbcels_int,Nbcels_ext,Nbcels_int,Nbcels_int,Nbcels_int,Nbcels_int,Nbcels_ext];
        
        ! CUBE INSCRIT DANS LA SPHERE
        DO Iz=1,Nbcels_z(1)
            DO Iy=1,Nbcels_y(1)
                DO Ix=1,Nbcels_x(1)                               
                    
                    Comp_cel = Comp_cel + 1
                    Nbc_p = Nbc_p + 1
                    
                    TmpCells(Comp_cel)%num_cell = Comp_cel
                    TmpCells(Comp_cel)%Xc = start_x(1) + (Ix-0.5)*Sc
                    TmpCells(Comp_cel)%Yc = start_y(1) + (Iy-0.5)*Sc
                    TmpCells(Comp_cel)%Zc = start_z(1) + (Iz-0.5)*Sc 
                    TmpCells(Comp_cel)%Sc = Sc                    
                Enddo
            Enddo
        Enddo
        Ncells_SphDomains(1) = Nbc_p;
        ! Les blocs entourant le bloc inscrit dans la sphere
        Do ii=2,7 
            Nbc_p_dom = 0; 
            DO Iz=1,Nbcels_z(ii)
                DO Iy=1,Nbcels_y(ii)
                    DO Ix=1,Nbcels_x(ii)              
                    
                        x = start_x(ii) + (Ix-0.5)*Sc
                        y = start_y(ii) + (Iy-0.5)*Sc
                        z = start_z(ii) + (Iz-0.5)*Sc                    
                    
                        hypotxy = hypot(x,y);
                        !theta_rd = Pi*nint(Pi/2.-atan2(z,hypotxy)/Pi*180.)/180.;
                        !phi_rd = Pi*nint(atan2(y,x)/Pi*180.)/180.; 
                        theta_rd = Pi/2.-atan2(z,hypotxy);
                        phi_rd = atan2(y,x);
                    
                        Rc = sqrt(x**2.+y**2.+z**2.);
                        Rlim = ap*(1.+cheb_eps*cos(cheb_l*theta_rd)*cos(cheb_l*phi_rd));
                                                            
                        If (Rc .le. Rlim) Then                                                                    
                            Comp_cel = Comp_cel + 1
                            Nbc_p = Nbc_p + 1;
                            Nbc_p_dom = Nbc_p_dom + 1;
                    
                            TmpCells(Comp_cel)%num_cell = Comp_cel
                            TmpCells(Comp_cel)%Xc = x
                            TmpCells(Comp_cel)%Yc = y
                            TmpCells(Comp_cel)%Zc = z
                            TmpCells(Comp_cel)%Sc = Sc                                                      
                        EndIf
                    Enddo
                Enddo
            Enddo 
            Ncells_SphDomains(ii) = Nbc_p_dom;
        EndDo
        If (Type_Par == 5) Then
            ! ici j'ajoute les cylindres 
        EndIf        
        Nbc = Nbc_p     
    EndIF
    
    Nbc = Comp_cel
    Allocate(Cells(Nbc));Cells(1:Nbc) = TmpCells(1:Nbc); Deallocate(TmpCells)
    
    ! We compute the maximum dimension of the scatterer in the 3 directions 
    ! (Scatterer xmin, xmax, ymin, ymax, zmin, zmax and Scatterer dx, dy and dz)
    ! these information will be used later inside Division_blocks.f90     
    Sc = SimScatterer%Sc;    
    ! before defining the box containing the scatterer, we should bring it back 
    ! to the vertical position (theta =0; Phi =0)
    SimScatterer%xmin= minval(Cells(1:Nbc)%Xc) - Sc/2.
    SimScatterer%xmax= maxval(Cells(1:Nbc)%Xc) + Sc/2.
    SimScatterer%ymin= minval(Cells(1:Nbc)%Yc) - Sc/2.
    SimScatterer%ymax= maxval(Cells(1:Nbc)%Yc) + Sc/2.
    SimScatterer%zmin= minval(Cells(1:Nbc)%Zc) - Sc/2.
    SimScatterer%zmax= maxval(Cells(1:Nbc)%Zc) + Sc/2.
    
    SimScatterer%dx= SimScatterer%xmax-SimScatterer%xmin
    SimScatterer%dy= SimScatterer%ymax-SimScatterer%ymin
    SimScatterer%dz= SimScatterer%zmax-SimScatterer%zmin    

END SUBROUTINE Discretization