!      /''''''''''''''''''''''''''''\
!     /          PARAMETERS          \
!     \............................../

! This module lists:
!      (1) the variables of the program
!      (2) the input parameters, as well as their values (in the subroutines)
! in order to be able to use them everywhere in the code.

MODULE parameters

  USE constants
  USE dimensions

  IMPLICIT NONE
  

  ! VARIABLES
  ! ---------

  INTEGER                                   :: cnt_conv     = 0,    & ! Counter for unreached convergence
                                               water_phase  = 0       ! Water phase at the planet's surface

  INTEGER,           DIMENSION(n_lay+1)     :: intrf        = 0,    & ! Position of the lower interface of layers
                                               intrf_old    = 0,    & ! Position of the interfaces at step -1
                                               intrf_oold   = 0,    & ! Position of the interfaces at step -2
                                               intrf_ooold  = 0,    & ! Position of the interfaces at step -2
                                               intrf_oooold = 0,    & ! Position of the interfaces at step -4
                                               intrf_old5   = 0,    & ! Position of the interfaces at step -5
					                                intrf_old6   = 0,    & ! Position of the interfaces at step -6
					                                intrf_old7   = 0       ! Position of the interfaces at step -7

  INTEGER,           DIMENSION(n_pts)       :: layer        = 0       ! Layer

  DOUBLE PRECISION                          :: x_corec      = 0.d0, & ! Core mass fraction, calculated by the model
                                               x_H2Oc       = 0.d0, & ! Water mass fraction, calculated by the model
                                               M_Pc         = 0.d0, & ! Planet mass, calculated by the model [kg]
                                               R_Pc         = 0.d0, & ! Planet radius, calculated by the model [m]
                                               FeSic        = 0.d0, & ! Fe/Si fraction, calculated by the model
                                               C_iner_tot   = 0.d0, & ! Polar moment of inertia of the planet,
                                                                      ! calculated by the model [kg.m2]
                                               C_MoI1       = 0.d0, & ! C/MR2
                                               C_MoI2       = 0.d0, & ! Cm/C
                                               g_mmax       = 0.d0, & ! Overall maximal gravity acceleration [m.s-2]
                                               g_max        = 0.d0, & ! Maximal gravity acceleration [m.s-2]
                                               g_max_old    = 0.d0, & ! Previous maximal gravity acceleration
                                                                      ! [m.s-2]
                                               P_mmax       = 0.d0, & ! Overall maximal pressure [Pa]
                                               P_max        = 0.d0, & ! Maximal pressure [Pa]
                                               P_max_old    = 0.d0, & ! Previous maximal pressure [Pa]
                                               T_mmax       = 0.d0, & ! Overall maximal temperature [K]
                                               T_max        = 0.d0, & ! Maximal temperature [K]
                                               T_max_old    = 0.d0, & ! Previous maximal temperature [K]
                                               rho_mmax     = 0.d0, & ! Overall maximal density [kg.m-3]
                                               rho_max      = 0.d0, & ! Maximal density [kg.m-3]
                                               rho_max_old  = 0.d0, & ! Previous maximal density [kg.m-3]
 					                                rho_surf     = 0.d0, &
                                               rho_surf_HHe = 0.d0, &
                                               rho_surf_w   = 0.d0, &
                                               rho_surf_oliv = 0.d0, &
                                               inv_rho_surf = 0.d0, &
                                               inverse_rho         = 0.d0, &
                                               rho_min_mat  = 0.d0, &
                                               gam_surf_HHe = 0.d0
      
                                               
  DOUBLE PRECISION,  DIMENSION(n_lay)       :: M            = 0.d0, & ! Mass of the layers [kg]
                                               C_iner       = 0.d0    ! Moment of inertia of the layers [kg.m2]

  DOUBLE PRECISION,  DIMENSION(n_pts)       :: r            = 0.d0    ! Radius [m]

  DOUBLE PRECISION,  DIMENSION(n_lay,n_pts) :: g            = 0.d0, & ! Gravity acceleration [m.s-2]
                                               P            = 0.d0, & ! Pressure [Pa]
                                               T            = 0.d0, & ! Temperature [K]
                                               rho          = 0.d0, & ! Density [kg.m-3]
					                                rho_old      = 0.d0, & ! Density of previous iteration [kg.m-3]
					                                rho_oold     = 0.d0, &
                                               gam          = 0.d0    ! Gruneisen parameter
  
  DOUBLE PRECISION, DIMENSION(n_pts) :: rho_oliv_core       = 0.d0, & 
                                        rho_oliv_env        = 0.d0, & 
                                        rho_oliv            = 0.d0, & 
                                        rho_w_core          = 0.d0, &
                                        rho_w_env           = 0.d0, &
                                        rho_w               = 0.d0, &
                                        rho_HHe             = 0.d0, &
                                        gam_oliv_core       = 0.d0, &
                                        gam_oliv_env        = 0.d0, &
                                        gam_oliv            = 0.d0, &
                                        gam_w_core          = 0.d0, &
                                        gam_w_env           = 0.d0, &
                                        gam_w               = 0.d0, &
                                        gam_HHe             = 0.d0, &
                                        cv_HHe_arr          = 0.d0, &
                                        cv_H2O_arr          = 0.d0, &
                                        cv_rock_arr         = 0.d0, &
                                        cv_core             = 0.d0, &
                                        cv_env              = 0.d0, &
                                        cv_all              = 0.d0, &
                                        entropy_H2O         = 0.d0, &
                                        entropy_rock        = 0.d0, &
                                        entropy_HHe         = 0.d0, &
                                        S_all               = 0.d0, &
                                        S_core              = 0.d0, &
                                        S_env               = 0.d0


  
  ! PARAMETERS
  ! ----------

  INTEGER                                       :: j_max,               & ! Maximum number of iterations
                                                   cnt_conv_max,        & ! Maximum number of detected oscillations
                                                   chk_EOS,             & ! Index to check the validity range of EOS
                                                   corEOS                 ! Type of correction for specific EOS

  DOUBLE PRECISION                              :: pow_law,             & ! Power exponent for spatial grid
                                                   conv_prec,           & ! Precision of the convergence condition
                                                   M_P,                 & ! Planet mass [kg]
                                                   R_P,                 & ! Planet radius [m]
                                                   x_core,              & ! Core mass fraction
                                                   x_H2O,               & ! Water mass fraction
                                                   f_alloy,             & ! Fraction of alloy in the core
                                                   MgD,                 & ! Mg number (Mg#)
                                                   MgSi,                & ! Mg/Si fraction
                                                   FeSi,                & ! Fe/Si fraction
                                                   FeSit,               & ! Fe/Si fraction (theoretical value)
                                                   P_surf,              & ! Surface pressure [Pa]
                                                   T_surf                 ! Surface temperature [K]

  DOUBLE PRECISION,  DIMENSION(n_EOS)           :: EOS_lim_P              ! Limit superior of the EOS's validity
                                                                          ! range in pressure [Pa]
       

  INTEGER,           DIMENSION(n_lay)           :: ilayer,              & ! List of layers
                                                   use_lay,             & ! Presence of each layer in the planet
                                                   iEOS,                & ! EOS to use in each layer
                                                   EOS_th,              & ! Consideration of the thermal component
                                                   n_mat                  ! Number of different materials per layer

  DOUBLE PRECISION,  DIMENSION(n_lay)           :: rho_lim_EOS,         & ! Computed lower limit in density of the
                                                   del_T,               & ! Upper temperature drop of the layers [K]
                                                   Mmol,                & ! Molar mass [kg.mol-1]
                                                   Z_eff,               & ! Effective atomic charge
                                                   rho_0,               & ! Density at ambiant conditions [kg.m-3]
                                                   T_0,                 & ! Reference temperature [K]
                                                   K_0,                 & ! Reference bulk modulus [Pa]
                                                   Kp_0,                & ! Pressure derivative of bulk modulus
                                                   gam_0,               & ! Reference Gruneisen parameter
                                                   q,                   & ! Adiabatic power exponant
                                                   a_T,                 & ! 1st alpha parameter [K-1]
                                                   b_T,                 & ! 2nd alpha parameter [K-2]
                                                   c_T,                 & ! 3rd alpha parameter [K]
                                                   a_P,                 & ! Temperature derivative of bulk modulus
                                                   T_D0                   ! Reference Debye temperature [K]

  DOUBLE PRECISION,  DIMENSION(n_lay,n_mat_max) :: x                      ! Mole fraction of layer's materials



  DOUBLE PRECISION ::   Zenv,       & ! Metallicity of envelope
                        OtoH_env,   & !
                        O_to_H_sun, & !
                        metallicity






  ! COEFFICIENTS
  ! ---------

  DOUBLE PRECISION, DIMENSION(3,8)   :: table_0       ! Table for ideal gas coefficients
  DOUBLE PRECISION, DIMENSION(15,56) :: table_r       ! Table for residual part coefficients

  DOUBLE PRECISION,  DIMENSION(8)    :: n0i,gamma0i   ! Ideal gas coefficients
  DOUBLE PRECISION,  DIMENSION(56)   :: cr, &         ! Residual part coefficients
					dr, &
					tr, &
					nr, &
					alpha, &
					beta, &
					gamma, &
					epsilon, &
				 	ar, &
					br, &
					B, &
					C, &
					D, &
					A



  ! Chabrier EOS tables
    DOUBLE PRECISION, DIMENSION(121*441) ::  logrho_ch, grad_ad_PT

    
    DOUBLE PRECISION, DIMENSION(121*241) ::  logU_ch, logP_ch, &
                                                    logS_ch, dlrho_dlT_P,        &
                                                    dlrho_dlP_T, dlS_dlT_P, grad_ad
                                                   

                                                  
  DOUBLE PRECISION, DIMENSION(121) ::  logT_input 
                                           
  DOUBLE PRECISION, DIMENSION(441) ::  logP_input 

  DOUBLE PRECISION, DIMENSION(241) ::  logrho_input 


 ! HG23 EOS correction for H/He mixture

  DOUBLE PRECISION, DIMENSION(441) :: logP_HG

  DOUBLE PRECISION, DIMENSION(121) :: logT_HG

  DOUBLE PRECISION, DIMENSION(441*121) :: Vmix_HG, Smix_HG


 ! Mazevet EOS for water

  DOUBLE PRECISION, DIMENSION(441) :: P_maz_water

  DOUBLE PRECISION, DIMENSION(121) :: T_maz_water

  DOUBLE PRECISION, DIMENSION(441*121) :: rho_maz_water


 ! SESAME EOS for rock (dry sand)

   DOUBLE PRECISION, DIMENSION(64) :: logT_sesame

   DOUBLE PRECISION, DIMENSION(41) :: logP_sesame

   DOUBLE PRECISION, DIMENSION(41*64) ::   logrho_sesame,      &
                                           logS_sesame,        &
                                           dlrho_dlT_p_sesame, &
                                           dlS_dlT_p_sesame,   &
                                           dlrho_dlP_t_sesame 
 
 



  ! Mazevet+ 22 EOS for H and He: tables for data

  ! DOUBLE PRECISION, DIMENSION(1001) :: rho_maz_HHe

  ! DOUBLE PRECISION, DIMENSION(29) :: temp_maz_HHe


  ! DOUBLE PRECISION, DIMENSION(1001*29) :: press_maz_HHe, U_maz_HHe

  
CONTAINS


  !----------------------------------------------!
  !                                              !
  !          SUBROUTINE read_EOS_maz22           !
  !                                              !
  !----------------------------------------------!
  !
  ! This reads Chabrier's EOS
  !
!  SUBROUTINE read_EOS_maz22
!
    !------

!    DOUBLE PRECISION, DIMENSION(5,1001*29) :: table_maz_HHe

!    !------

!  OPEN(UNIT = 158,file="Input/Mazevet_HHe/rho_maz.dat", status='OLD', action = 'READ') 
!  READ(158,*) rho_maz_HHe 
!  CLOSE(UNIT = 158)


!  OPEN(UNIT = 159,file="Input/Mazevet_HHe/temp_maz.dat", status='OLD', action = 'READ') 
!  READ(159,*) temp_maz_HHe 
!  CLOSE(UNIT = 159)

   
!  OPEN(UNIT = 160,file="Input/Mazevet_HHe/table1.dat", status='OLD', action = 'READ') 
!  READ(160,*) table_maz_HHe 
!  CLOSE(UNIT = 160)

!  press_maz_HHe = table_maz_HHe(3,1:)
!  U_maz_HHe = table_maz_HHe(4,1:)




    !------



!  END SUBROUTINE read_EOS_maz22





  !----------------------------------------------!
  !                                              !
  !          SUBROUTINE calc_FeSi_theo           !
  !                                              !
  !----------------------------------------------!
  !
  ! Calculates theoretical FeSi
  !
  SUBROUTINE calc_FeSi_theo

    ! THEORETICAL Fe/Si RATIO OF THE PLANET
    ! -------------------------------------

    FeSit = MgSi/MgD * (1.d0 - MgD + (Mmol(2)/Mmol(1))*(x_core/(1.d0-x_core-x_H2O)))

  END SUBROUTINE calc_FeSi_theo



  !-----------------------------------------------------!
  !                                                     !
  !          SUBROUTINE reinitialize_variables          !
  !                                                     !
  !-----------------------------------------------------!
  !
  ! Re-initializes all variables (not the parameters, that are read in a file).
  !
  SUBROUTINE reinitialize_variables

    INTEGER :: k, & ! Layer index
               l    ! Coefficients index

    !------

    ! Integer
    cnt_conv     = 0
    water_phase  = 0

    intrf        = 0
    intrf_old    = 0
    intrf_oold   = 0
    intrf_ooold  = 0
    intrf_oooold = 0

    layer        = 0

    ! Double precision
    x_corec      = 0.d0
    x_H2Oc       = 0.d0
    M_Pc         = 0.d0
    R_Pc         = 0.d0
    FeSic        = 0.d0
    g_mmax       = 0.d0
    g_max        = 0.d0
    g_max_old    = 0.d0
    P_mmax       = 0.d0
    P_max        = 0.d0
    P_max_old    = 0.d0
    T_mmax       = 0.d0
    T_max        = 0.d0
    T_max_old    = 0.d0
    rho_mmax     = 0.d0
    rho_max      = 0.d0
    rho_max_old  = 0.d0

    M            = 0.d0

    r            = 0.d0

    g            = 0.d0
    P            = 0.d0
    T            = 0.d0
    rho          = 0.d0
    gam          = 0.d0
    
    rho_old      = 0.d0
    rho_surf      = 0.d0
    rho_oold     = 0.d0

    !------

  END SUBROUTINE reinitialize_variables

  



  !----------------------------------------------!
  !                                              !
  !          SUBROUTINE read_coefficients        !
  !                                              !
  !----------------------------------------------!
  !
  ! This reads the coefficients for the Gruneisen parameter
  !
  SUBROUTINE read_coefficients

  INTEGER :: i  ! Index
 
  !OPEN(UNIT = 170,file='Input/IAPWS95-0.dat', status='OLD', action = 'READ')
  !READ(170,*) table_0
  !CLOSE(UNIT = 170)

  !n0i = table_0(2,:)
  !gamma0i = table_0(3,:)

  n0i(1) = -8.320446484d0
  n0i(2) = 6.683210528d0
  n0i(3) = 3.00632d0
  n0i(4) = 0.012436d0
  n0i(5) = 0.97315d0
  n0i(6) = 1.2795
  n0i(7) = 0.96956
  n0i(8) = 0.24873

  gamma0i(1) = 0.d0
  gamma0i(2) = 0.d0
  gamma0i(3) = 0.d0
  gamma0i(4) = 1.28728967d0
  gamma0i(5) = 3.53734222d0
  gamma0i(6) = 7.74073708d0
  gamma0i(7) = 9.24437796d0
  gamma0i(8) = 27.5075105d0

  !OPEN(UNIT = 171,file='Input/IAPWS95-r.dat', status='OLD', action = 'READ')
  !READ(171,*) table_r
  !CLOSE(UNIT = 171)

  !cr = table_r(2,:)
  !dr = table_r(3,:)
  !tr = table_r(4,:)
  !nr = table_r(5,:)
  !alpha = table_r(6,:)
  !beta = table_r(7,:)
  !gamma = table_r(8,:)
  !epsilon = table_r(9,:)
  !ar = table_r(10,:)
  !br = table_r(11,:)
  !B = table_r(12,:)
  !C = table_r(13,:)
  !D = table_r(14,:)
  !A = table_r(15,:)

  DO i = 1,7
     cr(i) = 0.d0
  END DO
  
  DO i = 8,22
     cr(i) = 1.d0
  END DO
  
  DO i = 23,42
     cr(i) = 2.d0
  END DO
 
  DO i = 43,46
    cr(i) = 3.d0
  END DO
  
  cr(47) = 4.d0
  
  DO i = 48,51
     cr(i) = 6.d0
  END DO
  
  DO i = 52,56
     cr(i) = 0.d0
  END DO

  DO i = 1,3
     dr(i) = 1.d0
  END DO
 
  dr(4) = 2.d0
  dr(5) = 2.d0
  dr(6) = 3.d0
  dr(7) = 4.d0
     
  DO i = 8,10
     dr(i) = 1.d0
  END DO

  dr(11) = 2.d0
  dr(12) = 2.d0
  dr(13) = 2.d0
  dr(14) = 4.d0
  dr(15) = 4.d0
  dr(16) = 5.d0
  dr(17) = 7.d0
  dr(18) = 9.d0
  dr(19) = 10.d0
  dr(20) = 11.d0
  dr(21) = 13.d0
  dr(22) = 15.d0
  dr(23) = 1.d0
  dr(24) = 2.d0
  dr(25) = 2.d0
  dr(26) = 2.d0
  dr(27) = 3.d0
 
  DO i = 28,30
     dr(i) = 4.d0
  END DO

  dr(31) = 5.d0
  dr(32) = 6.d0
  dr(33) = 6.d0
  dr(34) = 7.d0
  
  DO i = 35,39
     dr(i) = 9.d0
  END DO
  
  dr(40) = 10.d0
  dr(41) = 10.d0
  dr(42) = 12.d0
  dr(43) = 3.d0
  dr(44) = 4.d0
  dr(45) = 4.d0
  dr(46) = 5.d0
  dr(47) = 14.d0
  dr(48) = 3.d0
  
  DO i = 49,51
     dr(i) = 6.d0
  END DO
  
  DO i = 52,54
     dr(i) = 3.d0
  END DO

  dr(55) = 0.d0
  dr(56) = 0.d0

  tr(1) = -0.5d0
  tr(2) = 0.875d0
  tr(3) = 1.d0
  tr(4) = 0.5d0
  tr(5) = 0.75d0
  tr(6) = 0.375d0
  tr(7) = 1.d0
  tr(8) = 4.d0
  tr(9) = 6.d0
  tr(10) = 12.d0
  tr(11) = 1.d0
  tr(12) = 5.d0
  tr(13) = 4.d0
  tr(14) = 2.d0
  tr(15) = 13.d0
  tr(16) = 9.d0
  tr(17) = 3.d0
  tr(18) = 4.d0
  tr(19) = 11.d0
  tr(20) = 4.d0
  tr(21) = 13.d0
  tr(22) = 1.d0
  tr(23) = 7.d0
  tr(24) = 1.d0
  tr(25) = 9.d0
  tr(26) = 10.d0
  tr(27) = 10.d0
  tr(28) = 3.d0
  tr(29) = 7.d0
  tr(30) = 10.d0
  tr(31) = 10.d0
  tr(32) = 6.d0
  tr(33) = 10.d0
  tr(34) = 10.d0
  tr(35) = 1.d0
  tr(36) = 2.d0
  tr(37) = 3.d0
  tr(38) = 4.d0
  tr(39) = 8.d0
  tr(40) = 6.d0
  tr(41) = 9.d0
  tr(42) = 8.d0
  tr(43) = 16.d0
  tr(44) = 22.d0
  tr(45) = 23.d0
  tr(46) = 23.d0
  tr(47) = 10.d0
  tr(48) = 50.d0
  tr(49) = 44.d0
  tr(50) = 46.d0
  tr(51) = 50.d0
  tr(52) = 0.d0
  tr(53) = 1.d0
  tr(54) = 4.d0
  tr(55) = 0.d0
  tr(56) = 0.d0

  nr(1) = 0.12533547935523d-1
  nr(2) = 0.78957634722828d1
  nr(3) = -0.87803203303561d1
  nr(4) = 0.31802509345418d0
  nr(5) = -0.26145533859358d0
  nr(6) = -0.78199751687981d-2
  nr(7) = 0.88089493102134d-2
  nr(8) = -0.66856572307965d0
  nr(9) = 0.20433810950965d0
  nr(10) = -0.66212605039687d-4
  nr(11) = -0.19232721156002d0
  nr(12) = -0.25709043003438d0
  nr(13) = 0.16074868486251d0
  nr(14) = -0.40092828925807d-1
  nr(15) = 0.39343422603254d-6 
  nr(16) = -0.75941377088144d-5
  nr(17) = 0.56250979351888d-3
  nr(18) = -0.15608652257135d-4
  nr(19) = 0.11537996422951d-8
  nr(20) = 0.36582165144204d-6
  nr(21) = -0.13251180074668d-11
  nr(22) = -0.62639586912454d-9
  nr(23) = -0.10793600908932d0
  nr(24) = 0.17611491008752d-1
  nr(25) = 0.22132295167546d0
  nr(26) = -0.40247669763528d0
  nr(27) = 0.58083399985759d0
  nr(28) = 0.49969146990806d-2
  nr(29) = -0.31358700712549d-1
  nr(30) = -0.74315929710341d0
  nr(31) = 0.47807329915480d0
  nr(32) = 0.20527940895948d-1
  nr(33) = -0.13636435110343d0
  nr(34) = 0.14180634400617d-1
  nr(35) = 0.83326504880713d-2
  nr(36) = -0.29052336009585d-1
  nr(37) = 0.38615085574206d-1
  nr(38) = -0.20393486513704d-1
  nr(39) = -0.16554050063734d-2
  nr(40) = 0.19955571979541d-2
  nr(41) = 0.15870308324157d-3
  nr(42) = -0.16388568342530d-4
  nr(43) = 0.43613615723811d-1
  nr(44) = 0.34994005463765d-1
  nr(45) = -0.76788197844621d-1
  nr(46) = 0.22446277332006d-1
  nr(47) = -0.62689710414685d-4
  nr(48) = -0.55711118565645d-9
  nr(49) = -0.19905718354408d0
  nr(50) = 0.31777497330738d0
  nr(51) = -0.11841182425981d0
  nr(52) = -0.31306260323435d2
  nr(53) = 0.31546140237781d2
  nr(54) = -0.25213154341695d4
  nr(55) = -0.14874640856724d0
  nr(56) = 0.31806110878444d0

  DO i = 1,51
    alpha(i) = 0.d0
    beta(i) = 0.d0
    gamma(i) = 0.d0
    epsilon(i) = 0.d0
  END DO
  
  DO i = 52,54
     alpha(i) = 20.d0
  END DO
 
  alpha(55) = 0.d0
  alpha(56) = 0.d0
 
  beta(52) = 150.d0
  beta(53) = 150.d0
  beta(54) = 250.d0
  beta(55) = 0.3d0
  beta(56) = 0.3d0

  gamma(52) = 1.21d0
  gamma(53) = 1.21d0
  gamma(54) = 1.25d0
  gamma(55) = 0.d0
  gamma(56) = 0.d0
  
  epsilon(52) = 1.d0
  epsilon(53) = 1.d0
  epsilon(54) = 1.d0
  epsilon(55) = 0.d0
  epsilon(56) = 0.d0

  DO i=1,54
     ar(i) = 0.d0
     br(i) = 0.d0
     B(i) = 0.d0
     C(i) = 0.d0
     D(i) = 0.d0
     A(i) = 0.d0
  END DO
  
  ar(55) = 3.5d0
  ar(56) = 3.5d0
  
  br(55) = 0.85d0
  br(56) = 0.95d0
  
  B(55) = 0.2d0
  B(56) = 0.2d0
 
  C(55) = 28.d0
  C(56) = 32.d0
  
  D(55) = 700.d0
  D(56) = 800.d0

  A(55) = 0.32d0
  A(56) = 0.32d0

    !------

  END SUBROUTINE read_coefficients

    
END MODULE parameters
