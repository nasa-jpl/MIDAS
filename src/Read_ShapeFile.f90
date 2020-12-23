SUBROUTINE Read_ShapeFile(info_p,pr_NBcels,Int_Dist,a1,a2,pr_lattice)

    USE Initialization
    USE common_variables
    
    Implicit NONE

    !IN/OUT 
    Character, INTENT(IN) :: info_p
    Integer, INTENT(OUT) :: pr_NBcels
    Real(kind=8), INTENT(OUT) :: Int_Dist
    Real(kind=8), Dimension(3), INTENT(OUT):: a1,a2
    Integer, Dimension(:,:), allocatable, INTENT(OUT):: pr_lattice
    
    !! LOCAL
    Integer :: ii,jj
    Character ll
    Character(47) cc
        
    ! Read the shape file 
    Open(11,File = trim(ShapeFilePath))
    ! Inter-Dipole Distance
    if ((info_p == 'a') .or. (info_p == 's')) Then
        read(11,'(a38,e12.8e3,a)') cc,Int_Dist,ll
    ElseIf (info_p == 'p') Then
        read(11,'(a47,e12.8e3,a)') cc,Int_Dist,ll
    EndIf    
    Int_Dist = Int_Dist * 1e-6;
    
    ! NB_cells 
    read(11,*) pr_NBcels,ll
    
    ! a1 and a2
    read(11,'(f9.4,f9.4,f9.4,a)') a1(1),a1(2),a1(3),ll
    read(11,'(f9.4,f9.4,f9.4,a)') a2(1),a2(2),a2(3),ll
    
    ! For the moment we neglect the folowing 3 lines 
    read(11,*),cc;read(11,*),cc;read(11,*),cc;
    
    ! Read the positions in the lattice of the pr_NBcels cells (previously dipoles)
    Allocate(pr_lattice(pr_NBcels,6));
    Do ii= 1,pr_NBcels
        read(11,'(i7,i5,i5,i5,i5,i5,i5)') jj,pr_lattice(ii,1),pr_lattice(ii,2),pr_lattice(ii,3),&
            pr_lattice(ii,4),pr_lattice(ii,5),pr_lattice(ii,6) 
        ! To use only if you need to read a shape file copied from a 'target.out'file
        !read(11,'(i7,i5,i4,i4,i2,i2,i2)') jj,pr_lattice(ii,1),pr_lattice(ii,2),pr_lattice(ii,3),&
        !    pr_lattice(ii,4),pr_lattice(ii,5),pr_lattice(ii,6) 
    EndDo   

END SUBROUTINE Read_ShapeFile
