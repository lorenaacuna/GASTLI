!      /'''''''''''''''''''''''''''\
!     /          FUNCTIONS          \
!     \............................./

! This module contains important functions that can be used everywhere in the code.

MODULE functions

  USE constants
  USE dimensions
  USE parameters
  USE mazevet
  !USE AQUA_EOS,only:DIM_OUTPUT,LOAD_TABLE_PT,INTERPOLATE_AQUA_PT
  !USE AQUA_EOS,only:LOAD_TABLE_RhoT,INTERPOLATE_AQUA_RHOT
  !USE AQUA_EOS,only:LOAD_TABLE_RhoU,INTERPOLATE_AQUA_RHOU

  IMPLICIT NONE

CONTAINS

!
! List of functions:
!
! interp2d_opt: 2d interpolation for tables
!
! EOS_chabrier: density and adiabatic gradient for H/He by Chabrier+ 21. Used by new_rho (subroutine)
! 
! U_chabrier: all derivatives for H/He by Chabrier+21. Used by function grun_HHe
!
! grun_mazevet: Gruneisen parameter for water by Mazevet+ 19. Used by new_T
!
! grun_AQUA: Gruneisen parameter for water by AQUA EOS. Only used in tests
!
! grun_HHe: Gruneisen parameter for H/He by Chabrier+ 21. Used by new_T
!
! grun_sesame: Gruneisen parameter for rock from SESAME's dry sand EOS. Used by new_T
!




pure function integrate(x, y) result(r)
    !! Calculates the integral of an array y with respect to x using the trapezoid
    !! approximation. Note that the mesh spacing of x does not have to be uniform.
    Double precision, intent(in)  :: x(:)         !! Variable x
    Double precision, intent(in)  :: y(size(x))   !! Function y(x)
    Double precision              :: r            !! Integral int(y(x) times dx)

    ! Integrate using the trapezoidal rule
    associate(n => size(x))
      r = sum((y(1+1:n-0) + y(1+0:n-1))*(x(1+1:n-0) - x(1+0:n-1)))/2
    end associate
  end function




  !-------------------------------------!
  !                                     !
  !         FUNCTION OtoH               !
  !                                     !
  !-------------------------------------!
  !

  FUNCTION OtoH(Zin,Pin,Tin)

    !------



    DOUBLE PRECISION, INTENT(IN)    :: Zin, &    ! Metallicity
                                       Pin, &    ! Pressure in Pa
                                       Tin       ! Temperature in K



    DOUBLE PRECISION :: mu_H,      &
                        mu_H2O,    &
                        mu_Z,      &
                        mu_X,      &
                        poly0,     &
                        poly1,     &
                        X,         &
                        PcritH, TcritH, &
                        OtoH


   !----------

    ! Molecular weights
    mu_H = 1.
    mu_H2O = 18.
    mu_Z = mu_H2O           ! Metals represented by water

    ! Critical point of H
    TcritH = 10**4.185
    PcritH = (10**11.788)/10

    ! Coefficients for transition
    poly0 = -14336053.846627245
    poly1 = 273871954709.12747


    X = (1-Zin)/1.379

    OtoH = 0.d0
    

    IF (Tin > TcritH) THEN 
        ! atomic
        mu_X = mu_H
        OtoH = (Zin/X)*(mu_X/mu_Z)/( 1 + 2*( (Zin/X)*(mu_X/mu_Z) ) )
    ELSE IF ((Tin < TcritH).AND.(Pin > poly0*Tin + poly1)) THEN
        ! atomic
        mu_X = mu_H
        OtoH = (Zin/X)*(mu_X/mu_Z)/( 1 + 2*( (Zin/X)*(mu_X/mu_Z) ) )
    ELSE 
        ! molecular
        mu_X = 2 * mu_H
        OtoH = (Zin / X) * (mu_X / mu_Z) / (2 + 2 * ((Zin / X) * (mu_X / mu_Z)))
    END IF

    RETURN

  END FUNCTION OtoH




  !-------------------------------------!
  !                                     !
  !         FUNCTION interp2d_opt       !
  !                                     !
  !-------------------------------------!
  !
  !  2d linear interpolation 
  !
  !  x_in: Value to evaluate in x-axis
  !  y_in: Value to evaluate in y-axis
  !  x_array: Array of x-axis grid. This one is the one that always changes i.e [a,b,c,...,a,b,c...]
  !                  Dimension is Nx.
  !  y_array: Array of y-axis grid. This one is the constant one i.e [d,d,d,...,e,e,e...]
  !                  Dimension is Ny
  !  z_array: Array with the parameter to interpolate. Dimension is Nx*Ny
  !  return: final interpolated value, z(x_in,y_in)

  FUNCTION interp2d_opt(x_in,y_in,x_array,y_array,z_array,Nx,Ny)

    !------


    INTEGER :: i_x1, &
               i_y1, &
               i_x2, &
               i_y2, &
               Nx, Ny, Ntot


    DOUBLE PRECISION, INTENT(IN)    :: x_in, &
                                       y_in, &
                                       x_array(Nx), &
                                       y_array(Ny), &
                                       z_array(Nx*Ny)


    DOUBLE PRECISION :: x1,  &
                        x2,  &
                        y1,  &
                        y2,  &
                        f11, & 
                        f22, &
                        f12, &
                        f21, &
                        cte, &
                        interp2d_opt 



   !----------


    i_x1 = MAXLOC(x_array, DIM = 1, MASK=(x_array <= x_in)) 
    i_x2 = i_x1 + 1
    i_y1 = MAXLOC(y_array, DIM = 1, MASK = (y_array <= y_in))
    i_y2 = i_y1 + 1

    IF (i_x2==1) THEN
       i_x2 = 2
       i_x1 = 1
    END IF

    IF (i_y2==1) THEN
       i_y2 = 2
       i_y1 = 1
    END IF


    IF (( (i_x1>Nx).OR.(i_x2>Nx) ).AND.( (i_y1>Ny).OR.(i_y2>Ny) )) THEN
       i_x1 = Nx
       i_x2 = Nx
       i_y1 = Ny
       i_y2 = Ny

       interp2d_opt = z_array((i_y2-1)*Nx+i_x1) 


    ELSE IF ((i_x1>Nx).OR.(i_x2>Nx)) THEN
       i_x1 = Nx
       i_x2 = Nx
      
       x1 = x_array(i_x1)
       x2 = x_array(i_x2)
       y1 = y_array(i_y1)
       y2 = y_array(i_y2)
 
       f11 = z_array((i_y1-1)*Nx+i_x1) ! -> this is your y0
       f12 = z_array((i_y2-1)*Nx+i_x1) ! -> this is your y1
       f21 = z_array((i_y1-1)*Nx+i_x2)
       f22 = z_array((i_y2-1)*Nx+i_x2)

!    WRITE(6,*) 'x-axis indexes > Nx'
!   
!    WRITE(6,*) 'i_x1 =',i_x1
!    WRITE(6,*) 'i_x2 =',i_x2
!    WRITE(6,*) 'i_y1 =',i_y1
!    WRITE(6,*) 'i_y2 =',i_y2
!
!    WRITE(6,*) 'x1 =',x1
!    WRITE(6,*) 'x2 =',x2
!    WRITE(6,*) 'y1 =',y1
!    WRITE(6,*) 'y2 =',y2
!
!    WRITE(6,*) 'f11 =',f11
!    WRITE(6,*) 'f12 =',f12
!    WRITE(6,*) 'f21 =',f21
!    WRITE(6,*) 'f22 =',f22


       interp2d_opt = f11 + (y_in - y1) * (f12-f11)/(y2-y1)

    ELSE IF ((i_y1>Ny).OR.(i_y2>Ny)) THEN
       i_y1 = Ny
       i_y2 = Ny

       x1 = x_array(i_x1)
       x2 = x_array(i_x2)
       y1 = y_array(i_y1)
       y2 = y_array(i_y2)
 
       f11 = z_array((i_y1-1)*Nx+i_x1) ! -> this is your y0
       f12 = z_array((i_y2-1)*Nx+i_x1) 
       f21 = z_array((i_y1-1)*Nx+i_x2) ! -> this is your y1
       f22 = z_array((i_y2-1)*Nx+i_x2)

!    WRITE(6,*) 'y-axis indexes > Ny'
!   
!    WRITE(6,*) 'i_x1 =',i_x1
!    WRITE(6,*) 'i_x2 =',i_x2
!    WRITE(6,*) 'i_y1 =',i_y1
!    WRITE(6,*) 'i_y2 =',i_y2
!
!    WRITE(6,*) 'x1 =',x1
!    WRITE(6,*) 'x2 =',x2
!    WRITE(6,*) 'y1 =',y1
!    WRITE(6,*) 'y2 =',y2
!
!    WRITE(6,*) 'f11 =',f11
!    WRITE(6,*) 'f12 =',f12
!    WRITE(6,*) 'f21 =',f21
!    WRITE(6,*) 'f22 =',f22


       interp2d_opt = f11 + (x_in - x1) * (f21-f11)/(x2-x1)



    ELSE


    x1 = x_array(i_x1)
    x2 = x_array(i_x2)
    y1 = y_array(i_y1)
    y2 = y_array(i_y2)

    ! Evaluate function
    f11 = z_array((i_y1-1)*Nx+i_x1)
    f12 = z_array((i_y2-1)*Nx+i_x1)
    f21 = z_array((i_y1-1)*Nx+i_x2)
    f22 = z_array((i_y2-1)*Nx+i_x2)

!    WRITE(6,*) 'Normal case'
!   
!    WRITE(6,*) 'i_x1 =',i_x1
!    WRITE(6,*) 'i_x2 =',i_x2
!    WRITE(6,*) 'i_y1 =',i_y1
!    WRITE(6,*) 'i_y2 =',i_y2
!
!    WRITE(6,*) 'x1 =',x1
!    WRITE(6,*) 'x2 =',x2
!    WRITE(6,*) 'y1 =',y1
!    WRITE(6,*) 'y2 =',y2
!
!    WRITE(6,*) 'f11 =',f11
!    WRITE(6,*) 'f12 =',f12
!    WRITE(6,*) 'f21 =',f21
!    WRITE(6,*) 'f22 =',f22

    cte = 1/((x2-x1)*(y2-y1))

!    WRITE(6,*) 'cte =',cte

    interp2d_opt = cte*( f11*(x2-x_in)*(y2-y_in) + f21*(x_in-x1)*(y2-y_in) + &
                        f12*(x2-x_in)*(y_in-y1) + f22*(x_in-x1)*(y_in-y1) )
    END IF

    RETURN

  END FUNCTION interp2d_opt




  !-------------------------------------!
  !                                     !
  !         FUNCTION EOS_chabrier       !
  !                                     !
  !-------------------------------------!
  !
  ! Bilinear interpolation of Chabrier EOS for H/He
  ! EOS_chabrier(1): density in kg/m3
  ! EOS_chabrier(2): Adiabatic gradient

  FUNCTION EOS_chabrier(logPchin,logTchin)

    !------

    DOUBLE PRECISION, INTENT(IN)          :: logPchin, logTchin

    DOUBLE PRECISION :: interp,        &
                        Y_astk,         &
                        Y_til,          &
                        rho_cd21,      &
                        logpress_conv, &
                        Vmix_value,    &
                        inv_rho_ideal, &
                        rho_ideal,     &
                        inv_rho_corr,  &
                        rho_corr

                        

    DOUBLE PRECISION, DIMENSION(2) :: EOS_chabrier 


   !----------

    Y_astk = 0.275d0
    Y_til = 0.245d0

    ! Density
    interp = interp2d_opt(logPchin,logTchin,logP_input,logT_input,logrho_ch,size(logP_input),size(logT_input))

    rho_cd21 = 10**interp         ! Reverse log. Density in g/cm3
    logpress_conv = logPchin + 10d0    ! log(P) in dyn/cm2
    Vmix_value = interp2d_opt(logTchin,logpress_conv,logT_HG,logP_HG,Vmix_HG,size(logT_HG),size(logP_HG))
    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                    !
    ! Original (correction at Y = 0.245) !
    !                                    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!    EOS_chabrier(1) = rho_cd21*1d3 ! Change of units from g/cm3 to kg/m3



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                    !
    !          No correction             !
    !                                    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !inv_rho_ideal = (1/rho_cd21) - Vmix_value*Y_til*(1-Y_til)
    !rho_ideal = 1d0 / inv_rho_ideal

    !EOS_chabrier(1) = rho_ideal*1d3 ! Change of units from g/cm3 to kg/m3



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                    !
    !   With correction at Y = 0.275     !
    !                                    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    inv_rho_corr = (1/rho_cd21) + Vmix_value * (Y_astk*(1-Y_astk) - Y_til*(1-Y_til))
    rho_corr = 1d0 / inv_rho_corr

    EOS_chabrier(1) = rho_corr*1d3 ! Change of units from g/cm3 to kg/m3


    EOS_chabrier(2) = 0.d0


    RETURN

  END FUNCTION EOS_chabrier


  !-------------------------------------!
  !                                     !
  !         FUNCTION U_chabrier         !
  !                                     !
  !-------------------------------------!
  !
  ! Bilinear interpolation of Chabrier internal energy for H/He
  ! 
  ! U_chabrier(1): pressure in GPa
  ! U_chabrier(2): internal energy in MJ/kg
  ! U_chabrier(3): entropy in MJ/kg/K
  ! U_chabrier(4): ( dlog rho/dlog T )_P 
  ! U_chabrier(5): ( dlog rho/dlog P )_T
  ! U_chabrier(6): ( dlog S/dlog T )_P
  ! U_chabrier(7): adiabatic gradient
  ! 

  FUNCTION U_chabrier(logRHOchin,logTchin)

    !------

    DOUBLE PRECISION, INTENT(IN)          :: logRHOchin, logTchin

    DOUBLE PRECISION :: logRHO1, logRHO2, logT1, logT2, f11, f22, f12, f21, cte, interp

    DOUBLE PRECISION, DIMENSION(7) :: U_chabrier 

    INTEGER :: i_rho1,i_t1,i_rho2,i_t2, Nrho


   !----------

    ! Pressure
    interp = interp2d_opt(logRHOchin,logTchin,logrho_input,logT_input,logP_ch,size(logrho_input),size(logT_input))

    U_chabrier(1) = 10**interp ! Reverse log. Units: GPa 



    ! Internal energy
    interp = interp2d_opt(logRHOchin,logTchin,logrho_input,logT_input,logU_ch,size(logrho_input),size(logT_input))

    U_chabrier(2) = 10**interp ! Reverse log. Units: MJ/kg



    ! Entropy
    interp = interp2d_opt(logRHOchin,logTchin,logrho_input,logT_input,logS_ch,size(logrho_input),size(logT_input))

    U_chabrier(3) = 10**interp ! Reverse log. Units: MJ/kg/K



    ! ( dlog rho/dlog T )_P
    interp = interp2d_opt(logRHOchin,logTchin,logrho_input,logT_input,dlrho_dlT_P,size(logrho_input),size(logT_input))

    U_chabrier(4) = interp ! Derivative



    ! ( dlog rho/dlog P )_T
    interp = interp2d_opt(logRHOchin,logTchin,logrho_input,logT_input,dlrho_dlP_T,size(logrho_input),size(logT_input))

    U_chabrier(5) = interp ! Derivative



    ! ( dlog S/dlog T )_P
    interp = interp2d_opt(logRHOchin,logTchin,logrho_input,logT_input,dlS_dlT_P,size(logrho_input),size(logT_input))

    U_chabrier(6) = interp ! Derivative



    ! Adiabatic gradient
    interp = interp2d_opt(logRHOchin,logTchin,logrho_input,logT_input,grad_ad,size(logrho_input),size(logT_input))

    U_chabrier(7) = interp ! Derivative


    RETURN

  END FUNCTION U_chabrier


  !------------------------------------!
  !                                    !
  !        FUNCTION grun_mazevet       !
  !                                    !
  !------------------------------------!
  !
  ! Function that calculates the Gruneisen parameter of supercritical water as
  ! in Mazevet et al. 2019
  ! 
  FUNCTION grun_mazevet(T_in,rho_in)

    !------

    DOUBLE PRECISION, INTENT(IN)         :: rho_in,   & ! Input density in kg.m3
                                            T_in        ! Input temperature in K

    DOUBLE PRECISION   :: grun_mazevet,    & ! Result: Gruneisen parameter 
                          toli,          & ! tolerance
                          T1,T2,         & ! Increased T
			  rho1,rho2,     & ! Increased rho
                          P1,P2,         & ! Increased P
                          u1,u2,         & ! Increased internal energy 
                          PnkT,FNkT,UNkT,CV,CHIT,CHIR,PMbar,USPEC    ! Output of Mazevet's function h2ofit    
   
  DOUBLE PRECISION :: ratio_w, ratio_oliv, gam_oliv


    !------
   
      toli = 0.2

      T1 = T_in - toli
      T2 = T_in + toli

      !rho1 = rho_in/1000d0
      !rho2 = rho_in*(1 + toli)/1000d0

      rho1 = rho_in/1000d0
      rho2 = rho_in/1000d0

      CALL h2ofit(rho1,T1,PnkT,FNkT,UNkT,CV,CHIT,CHIR,PMbar,USPEC)
      P1 = PMbar*1d11  ! Pressure in Pa
      u1 = (USPEC*rho1)*0.1d0
      CALL h2ofit(rho2,T2,PnkT,FNkT,UNkT,CV,CHIT,CHIR,PMbar,USPEC)
      P2 = PMbar*1d11  ! Pressure in Pa
      u2 = (USPEC*rho2)*0.1d0



      grun_mazevet = (P2-P1)/(u2-u1)
      RETURN

    !------

     END FUNCTION grun_mazevet


  !------------------------------------!
  !                                    !
  !        FUNCTION heat_cap_water     !
  !                                    !
  !------------------------------------!
  !
  ! Function that calculates the heat capacity Cv of supercritical water as
  ! in Mazevet et al. 2019
  ! 
  FUNCTION heat_cap_water(T_in,rho_in)

    !------

    DOUBLE PRECISION, INTENT(IN)         :: rho_in,   & ! Input density in kg.m3
                                            T_in        ! Input temperature in K

    DOUBLE PRECISION   :: toli, T1, T2,    & ! Temp
			  rho1, rho2,      & ! rho
                          Cv_kb, k_b, m_h, & 
                          F1, F2,          &
                          PnkT,FNkT,UNkT,CV,CHIT,CHIR,PMbar,USPEC    ! Output of Mazevet's function h2ofit    
   
    DOUBLE PRECISION, DIMENSION(2) :: heat_cap_water ! heat_cap_water(1): Cv, heat capacity in J/kg/K
                                                    ! heat_cap_water(2): S, entropy 

    !------

      k_b = 1.380649d-23   ! Boltzmann constant in J/K
      m_h = 1.6735d-24     ! Mass of hydrogen atom in grams

   
      ! Heat capacity calculation
      
      CALL h2ofit(rho_in/1000d0,T_in,PnkT,FNkT,UNkT,CV,CHIT,CHIR,PMbar,USPEC)
      Cv_kb = CV
      heat_cap_water(1) = Cv_kb * (k_b/m_h) * 1d3   ! In J/kg/K


     ! Entropy calculation

     toli = 0.2

     T1 = T_in - toli
     T2 = T_in + toli

     rho1 = rho_in/1000d0
     rho2 = rho_in/1000d0

     CALL h2ofit(rho1,T1,PnkT,FNkT,UNkT,CV,CHIT,CHIR,PMbar,USPEC)
     F1 = FNkT * k_b*T1/m_h                                       ! In J/kg
     CALL h2ofit(rho2,T2,PnkT,FNkT,UNkT,CV,CHIT,CHIR,PMbar,USPEC)
     F2 = FNkT * k_b*T2/m_h                                       ! In J/kg


     heat_cap_water(2) = (-1d0)*(F2-F1)/(T2-T1)                   ! In J/kg/K
     

      RETURN

    !------

     END FUNCTION heat_cap_water




  !------------------------------------!
  !                                    !
  !        FUNCTION grun_AQUA          !
  !                                    !
  !------------------------------------!
  !
  ! Function that calculates the Gruneisen parameter of supercritical water as
  ! in Haldemann et al. 2020
  ! 
!  FUNCTION grun_AQUA(T_in,rho_in)
!
!    !------
!
!    DOUBLE PRECISION, INTENT(IN)         :: rho_in,   & ! Input density in kg.m3
!                                             T_in        ! Input temperature in K
!
!    DOUBLE PRECISION   :: grun_AQUA,    & ! Result: Gruneisen parameter 
!                          toli,          & ! tolerance
!                          T1,T2,         & ! Increased T
!  			  rho1,rho2,     & ! Increased rho
!                          P1,P2,         & ! Increased P
!                          u1,u2            ! Increased internal energy 
! 
!    REAL(8)  :: AQUA_RES(DIM_OUTPUT)
!    CHARACTER(LEN=200) :: path
!    integer, save :: counter_eos = 0
!
!    !------
!   
!      toli = 0.2
!
!      T1 = T_in - toli
!      T2 = T_in + toli
!
!
!
!
!
!    !IF ((T_in>150d0).AND.(T_in<1d5).AND.(rho_in>0d0).AND.(rho_in<1d4)) THEN
!
! 
!        AQUA_RES = INTERPOLATE_AQUA_RHOT(max(1d-1, rho_in),max(151d0, T1))
!        P1 = AQUA_RES(1) 
!        u1 = rho_in*AQUA_RES(4)
!        AQUA_RES = INTERPOLATE_AQUA_RHOT(max(1d-1, rho_in),max(151d0,T2))
!        P2 = AQUA_RES(1) 
!        u2 = rho_in*AQUA_RES(4)
!
!    !ELSE
!    !    WRITE(6,*) 'Out of table'
!    !    WRITE(6,*) T_in, rho_in
!    !    rho1 = rho_in/1000d0
!    !    rho2 = rho_in/1000d0
!    !
!    !    CALL h2ofit(rho1,T1,PnkT,FNkT,UNkT,CV,CHIT,CHIR,PMbar,USPEC)
!    !    P1 = PMbar*1d11  ! Pressure in Pa
!    !    u1 = (USPEC*rho1)*0.1d0
!    !    CALL h2ofit(rho2,T2,PnkT,FNkT,UNkT,CV,CHIT,CHIR,PMbar,USPEC)
!    !    P2 = PMbar*1d11  ! Pressure in Pa
!    !    u2 = (USPEC*rho2)*0.1d0
!    !END IF
!
!
!      grun_AQUA = (P2-P1)/(u2-u1)
!      RETURN
!
!    !------
!
!     END FUNCTION grun_AQUA



  !------------------------------------!
  !                                    !
  !        FUNCTION grun_HHe           !
  !                                    !
  !------------------------------------!
  !
  ! Function that calculates the Gruneisen parameter of H/He
  ! with Chabrier EOS
  !
  ! 
  ! grun_HHe(1): Gruneisen parameter
  ! grun_HHe(2): Cv - Heat capacity at constant volume (J/kg/K)
  ! grun_HHe(3): Cv - Heat capacity at constant pressure (J/kg/K)
  ! 


  FUNCTION grun_HHe(T_in,rho_in)

    !------


    DOUBLE PRECISION, INTENT(IN)         ::  rho_in,   & ! Input density in kg.m3
                                             T_in        ! Input temperature in K


    DOUBLE PRECISION   :: press, rho_un, &
                          entr,          &
                          dlrho_dlT,     &
                          dlrho_dlP,     &
                          dlS_dlT,       &
                          chi_T,         &
                          chi_rho,       &
                          Cp, Cv,        &
                          delta_ad,      &
                          toli,          & ! tolerance
                          T1,T2,         & ! Increased T
  			  rho1,rho2,     & ! Increased rho
                          P1,P2,         & ! Increased P
                          u1,u2,         & ! Increased internal energy 
                          grun_scalar

                          
                   
  
    DOUBLE PRECISION, DIMENSION(7) :: output, output1, output2

    DOUBLE PRECISION, DIMENSION(3) :: grun_HHe       ! Result: Gruneisen parameter 
 

  
    !------

 
      rho_un = rho_in*1e-3                                ! In g/cm3
      output = U_chabrier(log10(rho_un),log10(T_in))      
      press = 1d9*output(1)                               ! Convert pressure from GPa to Pa
      entr = 1d6*output(3)                                ! Convert entropy from MJ/kg/K to J/kg/K
      dlrho_dlT = output(4)
      dlrho_dlP = output(5)
      dlS_dlT = output(6)
      delta_ad = output(7)


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                    !
      !         Method 2/Method 3          !
      !                                    !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      chi_T = dlrho_dlT / dlrho_dlP
      chi_rho = 1d0 / dlrho_dlP

      Cp = entr * dlS_dlT
      Cv = Cp - (press/(rho_in*T_in)) * (chi_T**2/chi_rho)
  
      ! old
      ! grun_HHe = -(press * chi_T)/(rho_in*T_in*Cp) 

       
      grun_scalar = -(press * chi_T)/(rho_in*T_in*Cv) 


      grun_HHe(1) = grun_scalar
      grun_HHe(2) = Cv
      grun_HHe(3) = Cp             ! In J/kg/K
     


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                    !
      !       As adiabatic gradient        !
      !                                    !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
!      grun_HHe = delta_ad




      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                    !
      !         Method 1   (dP/dU)         !
      !                                    !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!      toli = 0.2

!      T1 = T_in - toli
!      T2 = T_in + toli


!      output1 = U_chabrier(log10(rho_in*1e-3),log10(T1))      
!      P1 = 1d9*output1(1)                                    ! Convert pressure from GPa to Pa
!      u1 = (1d6*output1(2))*rho_in                           ! Convert internal energy from MJ/kg to J/kg, and then to J/m3
      

!      output2 = U_chabrier(log10(rho_in*1e-3),log10(T2))      
!      P2 = 1d9*output2(1)                                    ! Convert pressure from GPa to Pa
!      u2 = (1d6*output2(2))*rho_in                           ! Convert internal energy from MJ/kg to J/kg, and then to J/m3

!      grun_HHe = (P2-P1)/(u2-u1)


      RETURN



    !------

     END FUNCTION grun_HHe


  !------------------------------------!
  !                                    !
  !       FUNCTION grun_sesame         !
  !                                    !
  !------------------------------------!
  !
  ! Function that calculates the Gruneisen parameter of rock
  ! with SESAME's dry sand EOS
  !
  ! grun_sesame(1): Gruneisen parameter
  ! grun_sesame(2): Cv - Heat capacity at constant volume (J/kg/K)
  ! grun_sesame(3): Cp - Heat capacity at constant pressure (J/kg/K)
  ! 


  FUNCTION grun_sesame(T_in,P_in)

    !------


    DOUBLE PRECISION, INTENT(IN)         ::  P_in,   &   ! Input pressure in Pa
                                             T_in        ! Input temperature in K


    DOUBLE PRECISION   :: logPin,      &
                          logTin,      &
                          dlrho_dlT,   &
                          dlrho_dlP,   &
                          dlS_dlT,     &
                          logS_cgs,    &
                          entr,        &
                          logrho_cgs,  &
                          chi_T, chi_rho, &
                          Cv, Cp, P_cgs, rho_cgs, grun_scalar

                   
    DOUBLE PRECISION, DIMENSION(4) :: grun_sesame ! Result: Gruneisen parameter 


  
    !------

    ! Use same method as the Gruneisen parameter of H/He (method 2)

      logPin = log10(P_in*10)
      logTin = log10(T_in)

      dlrho_dlT = interp2d_opt(logTin,logPin,logT_sesame,logP_sesame, &
               dlrho_dlT_p_sesame,size(logT_sesame),size(logP_sesame))

      dlrho_dlP = interp2d_opt(logTin,logPin,logT_sesame,logP_sesame, &
               dlrho_dlP_t_sesame,size(logT_sesame),size(logP_sesame))

      dlS_dlT = interp2d_opt(logTin,logPin,logT_sesame,logP_sesame, &
               dlS_dlT_p_sesame,size(logT_sesame),size(logP_sesame))

      logS_cgs = interp2d_opt(logTin,logPin,logT_sesame,logP_sesame, &
               logS_sesame,size(logT_sesame),size(logP_sesame))

      entr = 10**logS_cgs

      logrho_cgs = interp2d_opt(logTin,logPin,logT_sesame,logP_sesame, &
               logrho_sesame,size(logT_sesame),size(logP_sesame))





      chi_T = dlrho_dlT / dlrho_dlP
      chi_rho = 1d0 / dlrho_dlP

      Cp = entr * dlS_dlT

      P_cgs = 10**logPin
      rho_cgs = 10**logrho_cgs

      Cv = Cp - (P_cgs/(rho_cgs*T_in)) * (chi_T**2/chi_rho)
  

!      WRITE(6,*) 'P_in = ', P_in  
!      WRITE(6,*) 'logPin = ', logPin
!      WRITE(6,*) 'logTin = ', logTin
!      WRITE(6,*) 'dlrho_dlT = ', dlrho_dlT
!      WRITE(6,*) 'dlrho_dlP = ', dlrho_dlP
!      WRITE(6,*) 'dlS_dlT = ', dlS_dlT
!      WRITE(6,*) 'logS_cgs = ', logS_cgs
!      WRITE(6,*) 'entr = ', entr
!      WRITE(6,*) 'logrho_cgs = ', logrho_cgs
!      WRITE(6,*) 'X_T = ', X_T
!      WRITE(6,*) 'Cv = ', Cv


     ! old
     ! grun_sesame = abs(   ( P_cgs * chi_T )/( rho_cgs*T_in*Cp )    )     


     grun_scalar = abs(   ( P_cgs * chi_T )/( rho_cgs*T_in*Cv )    ) 

     grun_sesame(1) = grun_scalar
     grun_sesame(2) = Cv*1d-4      ! From erg/g/K to J/kg/K
     grun_sesame(3) = Cp*1d-4      ! From erg/g/K to J/kg/K
     grun_sesame(4) = entr*1d-4    ! From erg/g/K to J/kg/K

     RETURN



    !------

     END FUNCTION grun_sesame





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                        !
!           Original functions           !
!                                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !-------------------------------------!
  !                                     !
  !          FUNCTION integral          !
  !                                     !
  !-------------------------------------!
  !
  ! Calculates the intagral of an array using the trapezoidal rule.
  !
  FUNCTION integral(lim_dn,lim_up,f,x)

    !------

    DOUBLE PRECISION                           :: integral  ! Result

    INTEGER,                        INTENT(IN) :: lim_dn, & ! Lower limit of the integral
                                                  lim_up    ! Upper limit of the integral

    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: f,      & ! Function to integrate
                                                  x         ! Variable of integration

    INTEGER                                    :: i         ! Index

    !------

    integral=0.d0

    !------

    IF (lim_up-lim_dn < 0) THEN
       DO i=lim_up, lim_dn-1
          integral = integral - (x(i+1)-x(i))*(f(i)+f(i+1))*0.5d0
       END DO
       
    ELSE IF (lim_up == lim_dn) THEN
       integral = 0.d0
       
    ELSE
       DO i=lim_dn, lim_up-1
          integral = integral + (x(i+1)-x(i))*(f(i)+f(i+1))*0.5d0
       END DO
    END IF

    RETURN

    !------
    
  END FUNCTION integral



  !----------------------------------!
  !                                  !
  !          FUNCTION ifNtZ          !
  !                                  !
  !----------------------------------!
  !
  ! Puts the input variable to zero if it is NaN (does nothing otherwise).
  !
  FUNCTION ifNtZ(var)

    !------

    DOUBLE PRECISION             :: ifNtZ ! Result

    DOUBLE PRECISION, INTENT(IN) :: var   ! Variable

    !------

    IF (isnan(var)) THEN
       ifNtZ = 0.d0
    ELSE
       ifNtZ = var
    END IF

    RETURN
       
    !------
    
  END FUNCTION ifNtZ



  !-------------------------------------!
  !                                     !
  !          FUNCTION mass_btw          !
  !                                     !
  !-------------------------------------!
  !
  ! Calculates the mass that is included between two given radii in the planet, or for a specific layer.
  !
  FUNCTION mass_btw(min,max,lay)

    !------

    DOUBLE PRECISION                   :: mass_btw ! Result

    INTEGER,          INTENT(IN)       :: min,   & ! Index of the lower radius
                                          max,   & ! Index of the upper radius
                                          lay      ! Layer to consider for the density

    INTEGER                            :: i,     & ! Spatial index
                                          k        ! Layer index

    DOUBLE PRECISION, DIMENSION(n_pts) :: arr      ! Temporary array

    !------

    mass_btw = 0.d0

    !------

    IF (min > max) THEN
       WRITE(6,*) '/!\ ERROR in function "mass_btw": lower radius must be smaller than upper radius.'
       WRITE(6,*) min
       WRITE(6,*) max
       WRITE(6,*) lay
       STOP
    END IF

    ! Use of the planet's density profile
    IF (lay==0) THEN
       DO k=1, n_lay
          arr = 0.d0

          IF ((min>=intrf(k).AND.min<intrf(k+1)) .AND. max>=intrf(k+1)) THEN
             FORALL(i=min:intrf(k+1)) arr(i) = rho(k,i) * r(i)**2.d0
            !mass_btw = mass_btw + 4.d0*pi*integral(min,intrf(k+1),arr,r)
            mass_btw = mass_btw + 4.d0*pi*integrate(r(min:intrf(k+1)),arr(min:intrf(k+1)))
            !mass_btw = mass_btw + 4.d0*pi*1.5d25
            !WRITE(6,*) 'case 1',integrate(r(min:intrf(k+1)),arr(min:intrf(k+1)))


          ELSE IF ((min>=intrf(k).AND.min<intrf(k+1)) .AND. (max>=intrf(k).AND.max<intrf(k+1))) THEN
             FORALL(i=min:max) arr(i) = rho(k,i) * r(i)**2.d0
             !mass_btw = mass_btw + 4.d0*pi*integral(min,max,arr,r)
             mass_btw = mass_btw + 4.d0*pi*integrate(r(min:max),arr(min:max))
             !mass_btw = mass_btw + 4.d0*pi*1.5d25
            !WRITE(6,*) 'case 2',integrate(r(min:max),arr(min:max))


          ELSE IF (min<intrf(k) .AND. (max>=intrf(k).AND.max<intrf(k+1))) THEN
             FORALL(i=intrf(k):max) arr(i) = rho(k,i) * r(i)**2.d0
             !mass_btw = mass_btw + 4.d0*pi*integral(intrf(k),max,arr,r)
             mass_btw = mass_btw + 4.d0*pi*integrate(arr(intrf(k):max),r(intrf(k):max))
             !mass_btw = mass_btw + 4.d0*pi*0.0
            !WRITE(6,*) 'case 3',integrate(arr(intrf(k):max),r(intrf(k):max))


          ELSE IF (min<intrf(k) .AND. max>=intrf(k+1)) THEN
             FORALL(i=intrf(k):intrf(k+1)) arr(i) = rho(k,i) * r(i)**2.d0
             !mass_btw = mass_btw + 4.d0*pi*integral(intrf(k),intrf(k+1),arr,r)
             mass_btw = mass_btw + 4.d0*pi*integrate(r(intrf(k):intrf(k+1)),arr(intrf(k):intrf(k+1)))
             !mass_btw = mass_btw + 4.d0*pi*1.5d25
            !WRITE(6,*) 'case 4',integrate(r(intrf(k):intrf(k+1)),arr(intrf(k):intrf(k+1)))

          
          ELSE
             CYCLE
          END IF
       END DO

       RETURN

    ! Use of a specific layer's density profile
    ELSE
       arr = rho(lay,:) * r**2.d0

       !mass_btw = 4.d0*pi*integral(min,max,arr,r)
       mass_btw = 4.d0*pi*integrate(r(min:max),arr(min:max))
       !WRITE(6,*) 'case 5',integrate(r(min:max),arr(min:max))
       !mass_btw = 4.d0*pi*1.5d25

       RETURN
    END IF
       
    !------
    
  END FUNCTION mass_btw



  !------------------------------------!
  !                                    !
  !          FUNCTION P_water          !
  !                                    !
  !------------------------------------!
  !
  ! P(T) equations of the different phase transition of water.
  !
  FUNCTION P_water(h,T1)

    !------

    DOUBLE PRECISION             :: P_water ! Result: pressure of the phase transition at T1

    INTEGER,          INTENT(IN) :: h       ! Curve index

    DOUBLE PRECISION, INTENT(IN) :: T1      ! Temperature

    DOUBLE PRECISION             :: Theta   ! Reduced temperature

    !------

    P_water = 0.d0

    !------

    ! Ih - Vapor (IAPWS 2011 Release)
    IF (h==1) THEN
       IF (T1<=Ttp(1)) THEN
          Theta = T1 / Ttp(1)
          P_water = Ptp(1) * exp(theta**(-1.d0) * (-0.212144006d2*theta**0.00333333333d0 &
               +0.273203819d2*theta**1.20666667d0 -0.610598130d1*theta**1.70333333d0))
          RETURN
       ELSE
          P_water = -1.d0
       END IF
       
    ! Ih - Liquid (IAPWS 2011 Release)
    ELSE IF (h==2) THEN
       IF (T1>=Ttp(2).AND.T1<=Ttp(1)) THEN
          Theta = T1 / Ttp(1)
          P_water = Ptp(1) * (1.d0 +0.119539337d7*(1.d0-theta**3.d0) +0.808183159d5*(1.d0-theta**25.75d0) &
               +0.333826860d4*(1.d0-theta**103.750d0))
          RETURN
       ELSE
          P_water = -1.d0
       END IF
       
    ! III - Liquid (IAPWS 2011 Release)
    ELSE IF (h==3) THEN
       IF (T1>=Ttp(2).AND.T1<=Ttp(3)) THEN
          Theta = T1 / Ttp(2)
          P_water = Ptp(2) * (1.d0 -0.299948d0*(1.d0-theta**60.d0))
          RETURN
       ELSE
          P_water = -1.d0
       END IF
       
    ! V - Liquid (Wagner & Pruss 2002)
    ELSE IF (h==4) THEN
       IF (T1>=Ttp(3).AND.T1<=Ttp(4)) THEN
          Theta = T1 / Ttp(3)
          P_water = Ptp(3) * (1.d0 -1.18721d0*(1.d0-theta**8.d0))
          RETURN
       ELSE
          P_water = -1.d0
       END IF
       
    ! VI - Liquid (Wagner & Pruss 2002)
    ELSE IF (h==5) THEN
       IF (T1>=Ttp(4).AND.T1<=Ttp(5)) THEN
          Theta = T1 / Ttp(4)
          P_water = Ptp(4) * (1.d0 -1.07476*(1.d0-theta**4.6d0))
          RETURN
       ELSE
          P_water = -1.d0
       END IF
       
    ! VII - Liquid (Frank et al. 2003)
    ELSE IF (h==6) THEN
       IF (T1>=Ttp(5)) THEN
          Theta = T1 / Ttp(5)
          P_water = Ptp(5) -0.764d9*(1.d0-theta**4.32d0)
          RETURN
       ELSE
          P_water = -1.d0
       END IF
       
    ! Liquid - Vapor (Wagner & Pruss 2002)
    ELSE IF (h==7) THEN
       IF (T1>=Ttp(1).AND.T1<=Ttp(6)) THEN
          Theta = T1 / Ttp(6)
          P_water = Ptp(6) * exp(theta**(-1.d0) * (-7.85951783d0*(1.d0-theta) +1.84408259d0*(1.d0-theta)**1.5d0 &
             -11.7866497d0*(1.d0-theta)**3.d0 +22.6807411*(1.d0-theta)**3.5d0 -15.9618719d0*(1.d0-theta)**4.d0 &
             +1.80122502d0*(1.d0-theta)**7.5d0))
          RETURN
       ELSE
          P_water = -1.d0
       END IF

    ! Unknown transition
    ELSE
       WRITE(6,*) '/!\ ERROR in function "P_water": the asked transition is unknown.'
       STOP
       
    END IF
       
    !------
    
  END FUNCTION P_water


  
  !------------------------------------!
  !                                    !
  !          FUNCTION dec2hex          !
  !                                    !
  !------------------------------------!
  !
  ! Converts a decimal integer "a" into its expression in hexadecimal.
  !
  FUNCTION dec2hex(a)

    !------

    CHARACTER(LEN=2)    :: dec2hex   ! Result (hexadecimal value)

    INTEGER, INTENT(IN) :: a         ! Number to convert from decimal to hexadecimal

    INTEGER             :: digit1, & ! Modulo of "a" by 16
                           digit2    ! Rest of "a" by 16

    CHARACTER           :: char1,  & ! First symbol of the hexadecimal value
                           char2     ! Second symbol of the hexadecimal value

    !------

    digit1 = a/16
    digit2 = mod(a,16)
    
    SELECT CASE (digit1)
    CASE (0)
       char1 = '0'
    CASE (1)
       char1 = '1'
    CASE (2)
       char1 = '2'
    CASE (3)
       char1 = '3'
    CASE (4)
       char1 = '4'
    CASE (5)
       char1 = '5'
    CASE (6)
       char1 = '6'
    CASE (7)
       char1 = '7'
    CASE (8)
       char1 = '8'
    CASE (9)
       char1 = '9'
    CASE (10)
       char1 = 'A'
    CASE (11)
       char1 = 'B'
    CASE (12)
       char1 = 'C'
    CASE (13)
       char1 = 'D'
    CASE (14)
       char1 = 'E'
    CASE (15)
       char1 = 'F'
    END SELECT

    SELECT CASE (digit2)
    CASE (0)
       char2 = '0'
    CASE (1)
       char2 = '1'
    CASE (2)
       char2 = '2'
    CASE (3)
       char2 = '3'
    CASE (4)
       char2 = '4'
    CASE (5)
       char2 = '5'
    CASE (6)
       char2 = '6'
    CASE (7)
       char2 = '7'
    CASE (8)
       char2 = '8'
    CASE (9)
       char2 = '9'
    CASE (10)
       char2 = 'A'
    CASE (11)
       char2 = 'B'
    CASE (12)
       char2 = 'C'
    CASE (13)
       char2 = 'D'
    CASE (14)
       char2 = 'E'
    CASE (15)
       char2 = 'F'
    END SELECT

    WRITE(dec2hex,'(2(A1))') char1, char2

    RETURN
    
    !------
    
  END FUNCTION dec2hex



END MODULE functions
