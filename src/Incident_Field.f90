SUBROUTINE Incident_Field(cel_init,size_Cells,Cells_in,Nb_transmitters,Transmitters,num_tr_start, num_tr_end,E_ref_incident)

    USE Initialization
    USE common_variables
    USE iso_fortran_env
    Implicit none
    
    !IN/OUT
    Integer, INTENT(IN) :: cel_init,size_Cells
    type (Cell), Dimension(size_Cells), INTENT(IN) :: Cells_in
    Integer, INTENT(IN) :: Nb_transmitters ! attention here could be the commen variable NTr or Nipws when calculating the CBFs 
    type (Dipole), Dimension(Nb_transmitters), INTENT(IN) :: Transmitters
    Integer, INTENT(IN) :: num_tr_start, num_tr_end
    COMPLEX(real64), Dimension(3*size_Cells,2*(num_tr_end-num_tr_start+1)), INTENT(OUT)::E_ref_incident

    ! local
    Complex :: K11x, K11y, K11z, Ex, Ey, Ez
    real(kind=8) :: theta_transmit, phi_transmit, Rx, Ry, Rz
    Integer :: NcalcTr,num_trans,num_Eref_v, num_Eref_h,num_cel, sol,curs_cel

    E_ref_incident = 0
    NcalcTr = num_tr_end-num_tr_start+1;

    DO num_trans=num_tr_start, num_tr_end
        num_Eref_v = num_trans - num_tr_start + 1; 
        num_Eref_h = num_Eref_v + NcalcTr;
        
        theta_transmit = Transmitters(num_trans)%theta
        phi_transmit = Transmitters(num_trans)%phi
  
	! The incident wave is propagating in positive x (if theta_i=phi_i=0)
        K11x = k_0*cos(theta_transmit*Pi/180.); 
        K11y = k_0*sin(theta_transmit*Pi/180.)*cos(phi_transmit*Pi/180.)
        K11z = k_0*sin(theta_transmit*Pi/180.)*sin(phi_transmit*Pi/180.) 
  
        curs_cel = 1
        DO num_cel = cel_init, cel_init+size_Cells-1
        Rx = Cells_in(num_cel)%Xc
        Ry = Cells_in(num_cel)%Yc
        Rz = Cells_in(num_cel)%Zc

        !!--------------------------------Polarisation Verticale---------------------------------
        Ex = - exp(J*K11x*Rx)*exp(J*K11y*Ry)*exp(J*K11z*Rz)*sin(theta_transmit*Pi/180.)  
        Ey = exp(J*K11x*Rx)*exp(J*K11y*Ry)*exp(J*K11z*Rz)*cos(theta_transmit*Pi/180.)*cos(phi_transmit*Pi/180.)
        Ez = exp(J*K11x*Rx)*exp(J*K11y*Ry)*exp(J*K11z*Rz)*cos(theta_transmit*Pi/180.)*sin(phi_transmit*Pi/180.)
    
	    E_ref_incident(3*(curs_cel-1)+1,num_Eref_v)= Ex
        E_ref_incident(3*(curs_cel-1)+2,num_Eref_v)= Ey
        E_ref_incident(3*(curs_cel-1)+3,num_Eref_v)= Ez 

        !!-------------------------------Polarisation Horizontale--------------------------------  
        Ex = 0.0
        Ey = - exp(J*K11x*Rx)*exp(J*K11y*Ry)*exp(J*K11z*Rz)*sin(phi_transmit*Pi/180.)
        Ez = exp(J*K11x*Rx)*exp(J*K11y*Ry)*exp(J*K11z*Rz)*cos(phi_transmit*Pi/180.)

        E_ref_incident(3*(curs_cel-1)+1, num_Eref_h)= Ex
        E_ref_incident(3*(curs_cel-1)+2, num_Eref_h)= Ey
        E_ref_incident(3*(curs_cel-1)+3, num_Eref_h)= Ez 
    
        curs_cel = curs_cel + 1
      Enddo  
    Enddo  


End Subroutine Incident_Field