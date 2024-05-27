
MODULE GASTLI_interior

  USE funcmazevet
  USE mazevet
  USE constants
  USE dimensions
  USE parameters
  USE functions
  USE modEOS
  USE subroutines


  IMPLICIT NONE

CONTAINS


SUBROUTINE GASTLI_interior_subroutine(j_maxin, cnt_conv_maxin, conv_precin,pow_lawin,chk_EOSin, EOS_lim_Pin, &
           corEOSin, Ttpin, Ptpin, &
           f_alloyin, MgDin, MgSiin, ilayerin, use_layin, iEOSin, EOS_thin, n_matin, del_Tin, xin, Mmolin, Z_effin, &
           rho_0in, T_0in, K_0in, Kp_0in, gam_0in, qin, a_Tin, b_Tin, c_Tin, a_Pin,T_D0in, &
           logT_sesamein,        & ! ok from here ...
           logP_sesamein,        &
           logrho_sesamein,      & 
           logS_sesamein,        &
           dlrho_dlT_p_sesamein, &
           dlS_dlT_p_sesamein,   &
           dlrho_dlP_t_sesamein, & ! ... to here
           P_maz_waterin,        & ! ok from here ...
           T_maz_waterin,        &
           rho_maz_waterin,      & ! ... to here
           logT_inputin,         & ! ok from here ...
           logP_inputin,         &
           logrho_inputin,       &
           logrho_chin,          &
           logU_chin,            &
           logP_chin,            &
           logS_chin,            &
           dlrho_dlT_Pin,        &
           dlrho_dlP_Tin,        &
           dlS_dlT_Pin,          &
           grad_adin,            &
           grad_ad_PTin,         & ! ... to here
           logP_HGin,            & ! ok from here...
           logT_HGin,            &
           Vmix_HGin,            & 
           Smix_HGin,            & ! ... to here
           M_Pin,x_corein,x_H2Oin, T_surfin,P_surfin, &
           Zenvin,                                    & ! NEW MODEL PLANET INPUT: ENVELOPE METALLICITY
           M_Pout,x_coreout,x_H2Oout,R_Pout,FeSiout,rhoout,rout,interfaceout,intrf_hist,iter_num, &
           gout,Tout,Pout,rhoarrout, cvout, Sout, &
           OtoHout,metalout)                            ! NEW MODEL PLANET OUTPUT: [Fe/H] with respect to Sun value






    ! INPUT

    DOUBLE PRECISION, INTENT(IN) ::   M_Pin,        &  ! Mass of planet [M_E]
                                      x_corein,     &  ! Core mass fraction 
                                      x_H2Oin,      &  ! Water mass fraction
                                      f_alloyin,    &  ! Alloy percentage in core
                                      MgDin,        &  ! # Mg number
                                      MgSiin,       &  ! Mg/Si mole ratio
                                      T_surfin,     &  ! Surface temperature [K]
                                      P_surfin,     &  ! Surface pressure [Pa]
                                      Zenvin,       &  ! Envelope metallicity
                                      pow_lawin,    &  ! Power exponent for spatial grid
                                      conv_precin      ! Precision of the convergence condition


    INTEGER, INTENT(IN)                         :: j_maxin,               & ! Maximum number of iterations
                                                   cnt_conv_maxin,        & ! Maximum number of detected oscillations
                                                   chk_EOSin,             & ! Index to check the validity range of EOS
                                                   corEOSin                 ! Type of correction for specific EOS


    DOUBLE PRECISION, DIMENSION(n_EOS), INTENT(IN)     :: EOS_lim_Pin   ! Limit superior of the EOS's validity
    
    DOUBLE PRECISION, DIMENSION(6), INTENT(IN) ::  Ttpin,                          & ! Temperature of water's triple points [K]
                                                   Ptpin                             ! Pressure of water's triple points [Pa]

    INTEGER, DIMENSION(n_lay), INTENT(IN)         :: ilayerin,              & ! List of layers
                                                     use_layin,             & ! Presence of each layer in the planet
                                                     iEOSin,                & ! EOS to use in each layer
                                                     EOS_thin,              & ! Consideration of the thermal component
                                                     n_matin                  ! Number of different materials per layer

    DOUBLE PRECISION, DIMENSION(n_lay), INTENT(IN):: del_Tin,               & ! Upper temperature drop of the layers [K]
                                                     Mmolin,                & ! Molar mass [kg.mol-1]
                                                     Z_effin,               & ! Effective atomic charge
                                                     rho_0in,               & ! Density at ambiant conditions [kg.m-3]
                                                     T_0in,                 & ! Reference temperature [K]
                                                     K_0in,                 & ! Reference bulk modulus [Pa]
                                                     Kp_0in,                & ! Pressure derivative of bulk modulus
                                                     gam_0in,               & ! Reference Gruneisen parameter
                                                     qin,                   & ! Adiabatic power exponant
                                                     a_Tin,                 & ! 1st alpha parameter [K-1]
                                                     b_Tin,                 & ! 2nd alpha parameter [K-2]
                                                     c_Tin,                 & ! 3rd alpha parameter [K]
                                                     a_Pin,                 & ! Temperature derivative of bulk modulus
                                                     T_D0in                   ! Reference Debye temperature [K]

    DOUBLE PRECISION,  DIMENSION(n_lay,n_mat_max), INTENT(IN) :: xin                      ! Mole fraction of layer's materials



    DOUBLE PRECISION, DIMENSION(64), INTENT(IN) :: logT_sesamein

    DOUBLE PRECISION, DIMENSION(41), INTENT(IN) :: logP_sesamein

    DOUBLE PRECISION, DIMENSION(41*64), INTENT(IN) ::    logrho_sesamein,      &
                                                         logS_sesamein,        &
                                                         dlrho_dlT_p_sesamein, &
                                                         dlS_dlT_p_sesamein,   &
                                                         dlrho_dlP_t_sesamein 


    DOUBLE PRECISION, DIMENSION(441), INTENT(IN) :: P_maz_waterin

    DOUBLE PRECISION, DIMENSION(121), INTENT(IN) :: T_maz_waterin

    DOUBLE PRECISION, DIMENSION(441*121), INTENT(IN) :: rho_maz_waterin


    DOUBLE PRECISION, DIMENSION(121*441), INTENT(IN) ::  logrho_chin, grad_ad_PTin

    
    DOUBLE PRECISION, DIMENSION(121*241), INTENT(IN) ::  logU_chin, logP_chin, &
                                                    logS_chin, dlrho_dlT_Pin,        &
                                                    dlrho_dlP_Tin, dlS_dlT_Pin, grad_adin
                                               

                                                  
    DOUBLE PRECISION, DIMENSION(121), INTENT(IN) ::  logT_inputin
                                           
    DOUBLE PRECISION, DIMENSION(441), INTENT(IN) ::  logP_inputin

    DOUBLE PRECISION, DIMENSION(241), INTENT(IN) ::  logrho_inputin 


    DOUBLE PRECISION, DIMENSION(441), INTENT(IN) :: logP_HGin

    DOUBLE PRECISION, DIMENSION(121), INTENT(IN) :: logT_HGin

    DOUBLE PRECISION, DIMENSION(441*121), INTENT(IN) :: Vmix_HGin, Smix_HGin



                                                                         
    ! OUTPUT

    DOUBLE PRECISION, INTENT(OUT) ::   M_Pout,     &   ! Mass of planet (re-check)
                                       x_coreout,  &   ! Core mass fraction (re-check)
                                       x_H2Oout,   &   ! Water mass fraction (re-check)
                                       R_Pout,     &   ! Radius of planet [R_E]
                                       FeSiout,    &   ! Fe/Si mole ratio
                                       rhoout,     &   ! Density of planet [g.cm-3]
                                       OtoHout,    &   !
                                       metalout


    DOUBLE PRECISION,  DIMENSION(n_pts), INTENT(OUT)       :: rout,        &   ! Radius [m]
                                                              gout,        &   !
                                                              Tout,        &   ! 
                                                              Pout,        &   !
                                                              rhoarrout,   &   !
                                                              cvout,       &   !
                                                              Sout
 


    INTEGER, DIMENSION(n_lay+1), INTENT(OUT)     :: interfaceout               ! Position of the lower interface of layers
    INTEGER, DIMENSION(n_lay+1,100), INTENT(OUT) :: intrf_hist                 ! interfaceout but for all iterations
    INTEGER, DIMENSION(100), INTENT(OUT) :: iter_num





  
  ! DECLARATION
  ! -----------

  INTEGER :: i, & ! Spatial index
             j, & ! Iteration index
             k, & ! Layer index
             my_unit, &
             i_center, i_core, i_totalpl

  DOUBLE PRECISION :: start, finish, logrho_rock_surf, rho_rock_surf, &
                      Ptest, Ttest, grun_test, rock_test, water_test, &
                      cv_sesame, cv_HHe, cv_water, S_cd21, Smix_value, &
                      Y_astk, Y_til


  LOGICAL :: test ! Test for loop exit


  CHARACTER(LEN=21)   :: name ! Name of the file


  DOUBLE PRECISION, DIMENSION(2) :: output_ch, HHe_test, output_heat_w
  DOUBLE PRECISION, DIMENSION(3) :: HHe_out2

  DOUBLE PRECISION, DIMENSION(4) :: sesame_out,  &
                                    sesame_out2


  DOUBLE PRECISION, DIMENSION(7) :: output_ch_S

  

  
  ! INITIALIZATION
  ! --------------

  CALL reinitialize_variables



!  call cpu_time(start)
       

  ! Model and EOS control parameters
  j_max = j_maxin 
  cnt_conv_max = cnt_conv_maxin 
  conv_prec = conv_precin
  pow_law = pow_lawin
  chk_EOS = chk_EOSin
  EOS_lim_P = EOS_lim_Pin 
  corEOS = corEOSin

  ! 'constants.in' file parameters
  Ttp = Ttpin
  Ptp = Ptpin

  ! Secondary compositional parameters
  f_alloy = f_alloyin
  MgD = MgDin
  MgSi = MgSiin


  ! 'layer.in' files parameters
  ilayer = ilayerin
  use_lay = use_layin
  iEOS = iEOSin 
  EOS_th = EOS_thin
  n_mat = n_matin
  del_T = del_Tin
  x = xin

  ! Calculated from 'layer.in' parameters
  Mmol = Mmolin                                                  
  Z_eff = Z_effin
  rho_0 = rho_0in
  !WRITE(6,*) 'RHO_0 IN =', rho_0in
  T_0 = T_0in
  K_0 = K_0in
  Kp_0 = Kp_0in
  gam_0 = gam_0in
  q = qin
  a_T = a_Tin
  b_T = b_Tin
  c_T = c_Tin
  a_P = a_Pin
  T_D0 = T_D0in  

  ! Primary compositional parameters, mass and surface conditions      
  M_P   = M_Pin   * M_E
  x_core = x_corein
  x_H2O = x_H2Oin
  T_surf = T_surfin
  P_surf = P_surfin
  Zenv = Zenvin


  ! EOS data parameters

  ! SESAME
  logT_sesame        = logT_sesamein
  logP_sesame        = logP_sesamein
  logrho_sesame      = logrho_sesamein
  logS_sesame        = logS_sesamein
  dlrho_dlT_p_sesame = dlrho_dlT_p_sesamein
  dlS_dlT_p_sesame   = dlS_dlT_p_sesamein
  dlrho_dlP_t_sesame = dlrho_dlP_t_sesamein

  ! Mazevet for water 
  P_maz_water   = P_maz_waterin
  T_maz_water   = T_maz_waterin
  rho_maz_water = rho_maz_waterin

  ! Chabrier+23 for H/He 
  logT_input     = logT_inputin
  logP_input     = logP_inputin
  logrho_input   = logrho_inputin
  logrho_ch      = logrho_chin
  logU_ch        = logU_chin
  logP_ch        = logP_chin
  logS_ch        = logS_chin
  dlrho_dlT_P    = dlrho_dlT_Pin
  dlrho_dlP_T    = dlrho_dlP_Tin
  dlS_dlT_P      = dlS_dlT_Pin
  grad_ad        = grad_adin
  grad_ad_PT     = grad_ad_PTin

  ! HG23 correction
  logP_HG = logP_HGin
  logT_HG = logT_HGin
  Vmix_HG = Vmix_HGin
  Smix_HG = Smix_HGin

 
  ! First subroutines
  CALL density_limit_EOS

!  call cpu_time(finish)
!  WRITE(6,*) 'DATA LOADING TIME = ',finish-start



  ! Envelope metallicity
  ! Zenv = 0.029d0                         ! Jupiter's value  -> Now this is a free input param.
  OtoH_env = OtoH(Zenv,P_surf,T_surf)      ! O:H ratio

  O_to_H_sun = 4.898e-4                    ! From Lodders+ 03
  metallicity = OtoH_env/O_to_H_sun

 ! WRITE(6,*) '                        '
 ! WRITE(6,*) '------------------------'
 ! WRITE(6,*) '      METALLICITY       '
 ! WRITE(6,*) 'Zenv = ', Zenv
 ! WRITE(6,*) 'O:H ratio =', OtoH_env
 ! WRITE(6,*) 'Metallicity = ', metallicity
 ! WRITE(6,*) '                        '
 ! WRITE(6,*) '------------------------'

  
  ! WRITE(6,*) 'Test: chabrier at logP=1, logT=4', EOS_chabrier(1.d0,4.d0)

  ! Opening log file if asked
  !IF (log_type==1) OPEN(6, FILE='log.txt', STATUS='new')

  ! Surface conditions
  !IF (x_H2O>0.d0) THEN
  !   CALL chk_water_phase
  !END IF

  ! Boundary conditions
  g(1,1)     = 0.d0   ! Central gravity acceleration
  P(n_lay,:) = P_surf ! Surface pressure
  T(n_lay,:) = T_surf ! Surface temperature

  !WRITE(6,*) 'FLAG1 = ', g(layer(1199),1199)
  ! Initialization of iteration counter
  j = 0

  ! Radius
  FORALL(i = 1:n_pts) r(i) = M_P**pow_law * dble(i-1)/dble(n_pts-1)

  
  !write(6,*) rho_0(1),rho_0(2),rho_0(3)


   ! CALCULATE SURFACE DENSITIES
   ! Olivine
   !WRITE(6,*) 'RHO_0 ARRAY = ', rho_0
   !WRITE(6,*) 'FLAGRHO FIRST = ', rho(k,1199)
   FORALL(k = 1:n_lay) rho(k,:) = rho_0(k)
   !WRITE(6,*) 'FLAGRHO AAA = ', rho(k,1199)
   rho_surf_oliv = rho_0(2)

   ! sesame dry sand
   !WRITE(6,*) 'logP = ', log10(P_surf*10)
   !WRITE(6,*) 'logT = ', log10(T_surf)
   !WRITE(6,*) 'logrho_sesame = ', logrho_sesame(1)

   logrho_rock_surf = interp2d_opt(log10(T_surf),log10(P_surf*10),logT_sesame,logP_sesame, &
               logrho_sesame,size(logT_sesame),size(logP_sesame))
   rho_rock_surf = (10**logrho_rock_surf)*1d3
  
   !WRITE(6,*) logrho_rock_surf

   ! Test for Mazevet EOS from file
   !WRITE(6,*) ''
   !WRITE(6,*) 'Test for Mazevet EOS for water (read from file) '
   !WRITE(6,*) 'P in GPa = ', P_surf/1e9
   !WRITE(6,*) 'T in K = ', T_surf
   !WRITE(6,*) 'Density in kg/m3 = ', interp2d_opt(P_surf/1e9,T_surf,P_maz_water,T_maz_water, &
   !            rho_maz_water,size(P_maz_water),size(T_maz_water))
   !WRITE(6,*) ''


   ! Test for all EOS - Jupiter's central P and T
   !WRITE(6,*) ''
   !WRITE(6,*) 'Test for all EOS - Jupiters central P and T '
   Ptest = 14.316e13
   Ttest = 18335.
   !WRITE(6,*) 'P in GPa = ', Ptest/1e9
   !WRITE(6,*) 'T in K = ', Ttest
   !WRITE(6,*) ''
   !WRITE(6,*) '-- Water --'
   !water_test = interp2d_opt(Ptest/1e9,Ttest,P_maz_water,T_maz_water, &
   !            rho_maz_water,size(P_maz_water),size(T_maz_water))
   !WRITE(6,*) 'Density in kg/m3 = ', water_test
   !WRITE(6,*) ''
   !WRITE(6,*) '-- Rock --'
   rock_test = interp2d_opt(log10(Ttest),log10(Ptest*10),logT_sesame,logP_sesame, &
                             logrho_sesame,size(logT_sesame),size(logP_sesame))
   !WRITE(6,*) 'Density in kg/m3 = ', (10**rock_test)*1d3
   !WRITE(6,*) ''
   !WRITE(6,*) '-- H/He --'
   !HHe_test = EOS_chabrier(log10(Ptest/1e9),log10(Ttest))
   !WRITE(6,*) 'Density in kg/m3 = ', HHe_test(1)
   !WRITE(6,*) ''

  ! WRITE(6,*) 'logT_sesame =', logT_sesame
  ! WRITE(6,*) 'logP_sesame =', logP_sesame
  ! WRITE(6,*) 'logrho_sesame =', logrho_sesame



   ! Water
   !CALL rho_surface_w    ! This calculates the parameter rho_surf_w 
                          ! automatically
   rho_surf_w  = interp2d_opt(P_surf/1e9,T_surf,P_maz_water,T_maz_water, &
               rho_maz_water,size(P_maz_water),size(T_maz_water))
   
   ! H/He
   output_ch = EOS_chabrier(log10(P_surf/1e9),log10(T_surf))
   rho_surf_HHe = output_ch(1)
   gam_surf_HHe = output_ch(2)
 
   
   ! INITIALISE DENSITY ARRAYS
   ! Olivine
   !rho_oliv(:) = rho_surf_oliv

   ! Water
   !rho_w(:) = rho_surf_w

   ! H/He
   !rho_HHe(:) = rho_surf_HHe


   ! Layer 1: core
   inverse_rho = 0.5d0/rho_surf_w + 0.5d0/rho_surf_oliv
   rho_0(1) =  1.d0/inverse_rho
   !rho(1,:) = rho_0(1)

   ! Layer 2: envelope   
   inverse_rho = (1-Zenv)/rho_surf_HHe + Zenv/rho_surf_w
   rho_0(2) = 1.d0/inverse_rho

   !rho_0(2) = rho_surf_w
   !rho(2,:) = rho_0(2)


   ! SURFACE DENSITY OF THE MIXTURE
   rho_surf = rho_0(2)

   !rho_0(3) = 0.d0

  IF (isnan(rho_surf)) THEN
             WRITE(6,*) 'Error in interior structure model (Fortran): surface density is NaN'
             WRITE(6,*) 'This is likely due to the surface pressure being so low that it is '
             WRITE(6,*) 'out of the EOSs validity range.'
             WRITE(6,*) 'We recommend to increase the surface pressure to 1000 bar ideally,'
             WRITE(6,*) 'or to 10 bar (minimum).'             
             STOP 
             !CYCLE
   END IF


   ! CALCULATE MINIMUM DENSITY TO HELP CONVERGENCE
   rho_min_mat = min(rho_surf_HHe,rho_surf_oliv,rho_surf_w)
   !rho_min_mat = 240d0

   ! Print this for debugging
   !WRITE(6,*) 'Surface densities are in kg/m3'
   !WRITE(6,*) 'Surface density of rock: ', rho_surf_oliv
   !WRITE(6,*) 'Surface density of rock (SESAME): ', rho_rock_surf
   !WRITE(6,*) 'Surface density of water: ', rho_surf_w
   !WRITE(6,*) 'Surface density of H/He: ', rho_surf_HHe
   !WRITE(6,*) 'Surface density of envelope mixture: ', rho_surf
   !WRITE(6,*) 'Overall density: ', rho_min_mat

   !WRITE(6,*) ""
   !WRITE(6,*) 'Test for SESAMEs Gruneisen parameter: '
   Ptest = 2.182728436018975d9
   Ttest = 5d3
   sesame_out = grun_sesame(Ttest,Ptest)
   grun_test = sesame_out(1)
   !WRITE(6,*) 'Gruneisen param = ', grun_test
   !WRITE(6,*) ""


   !WRITE(6,*) ""
   !WRITE(6,*) 'Test for heat capacities: '
   ! Rock
   sesame_out2 = grun_sesame(T_surf,P_surf)
   cv_sesame = sesame_out2(2)
   
   ! H/He
   HHe_out2 = grun_HHe(T_surf,rho_surf_HHe)
   cv_HHe = HHe_out2(2)

   ! Water
   output_heat_w = heat_cap_water(T_surf,rho_surf_w)
   cv_water = output_heat_w(1)


   !WRITE(6,*) 'Rock - Cv [J/kg/K] =  ', cv_sesame
   !WRITE(6,*) 'H/He - Cv [J/kg/K] =  ', cv_HHe
   !WRITE(6,*) 'Water - Cv [J/kg/K] = ', cv_water

   !WRITE(6,*) ""



   ! Surface rock density with sesame
   ! Comment this line for Vinet EOS for olivine
   rho_surf_oliv = rho_rock_surf

 !  IF (x_H2O>0d0) THEN 
 !      rho_surf = 0.5d0/rho_surf_HHe + 0.5d0/rho_0(1)
 !  ELSE 
 !      rho_surf = 0.5d0/rho_surf_w + 0.5d0/rho_0(1)
 !  END IF
  

 !  WRITE(name,'(A18)') 'Output/rho_w00.dat'
 !  my_unit = 200
 !
 !
 !  OPEN(my_unit, FILE=name, STATUS='unknown')
 !  WRITE(my_unit,'(A5)') 'rho_w'
 !  DO i=1, n_pts
 !    WRITE(my_unit,'(ES16.4E2)') rho_w(i)
 !  END DO
 !  CLOSE(my_unit)


!   WRITE(name,'(A20)') 'Output/rho_HHe00.dat'
!   my_unit = 200

!   OPEN(my_unit, FILE=name, STATUS='unknown')
!   WRITE(my_unit,'(A100)') 'i  r  rho_oliv  rho_w  rho_HHe  rho(1,i)  rho(2,i)  rho(3,i)  rho_planet'
!   DO i=1, n_pts
!     WRITE(my_unit,'(I4,8ES16.4E2)') i, r(i), rho_oliv(i), rho_w(i), rho_HHe(i), rho(1,i), rho(2,i), rho(3,i), rho(layer(i),i)
!   END DO
!   CLOSE(my_unit)




  ! Interfaces
  FORALL(k = 1:n_lay) intrf(k) = (k-1)*n_pts/20 +1
  intrf(n_lay+1) = n_pts + 1
  CALL layering

  ! Display input
!  CALL display_input

  
  ! Gravity
  !WRITE(6,*) 'FLAG2 = ', g(layer(1199),1199)
  !WRITE(6,*) 'FLAGRHO = ', rho(k,1199)
  CALL new_g
  !WRITE(6,*) 'g(layer(1199),1199) = ', g(layer(1199),1199)


  !WRITE(6,*) 'GRAVITY OK'

  ! Pressure
  CALL new_P

  !WRITE(6,*) 'PRESSURE OK'


  ! Temperature

  !  call cpu_time(start)
    CALL new_T         
  !  call cpu_time(finish)
  !  WRITE(6,*)'("Time = ",f6.3," seconds.")',finish-start


  !WRITE(6,*) 'TEMPERATURE OK'


  ! Density
  !rho_old = rho


  CALL new_rho

  !WRITE(6,*) 'DENSITY OK'



  !rho_new = 0.7d0 * rho + 0.3d0 * rho_old
  !rho_new 
  !rho = rho_new


  !WRITE(6,*) 'flag1'

  !WRITE(6,*) intrf
  ! Current planet mass and radius
  M_Pc = mass_btw(1,intrf(n_lay),0)
  R_Pc = r(intrf(n_lay))

  ! Results
  !CALL display_log(j)
  !CALL write_profile(j)

!  WRITE(11,'(I4)',advance='no') j
!  DO k=1, n_lay
!     WRITE(11,'(ES27.15E2)',advance='no') r(intrf(k))/R_E
!  END DO
!  WRITE(11,*) ''
  

!   WRITE(name,'(A20)') 'Output/rho_HHe01.dat'
!   my_unit = 200

!   OPEN(my_unit, FILE=name, STATUS='unknown')
!   WRITE(my_unit,'(A34)') 'rho_oliv  rho_w  rho_HHe  rho(2,i)'
!   DO i=1, n_pts
!     WRITE(my_unit,'(4ES16.4E2)') rho_oliv(i), rho_w(i), rho_HHe(i), rho(2,i)
!   END DO
!   CLOSE(my_unit)



  ! ITERATION PROCESS
  ! -----------------

  DO
     j = j+1
     
     ! Interfaces
     !WRITE(6,*) 'main'
     !WRITE(6,*) intrf
     !WRITE(6,*) j
     CALL interfaces          
     CALL layering

     ! Gravity
     !WRITE(6,*) 'g(layer(1199),1199) = ', g(layer(1199),1199)
     CALL new_g
     !WRITE(6,*) 'g(layer(1199),1199) = ', g(layer(1199),1199)

     
     ! Pressure
     CALL new_P

     ! Temperature
  !  call cpu_time(start)
    CALL new_T         
  !  call cpu_time(finish)
  !  WRITE(6,*) 'TIME = ',finish-start


     rho_oold = rho_old
     rho_old = rho
     CALL new_rho
     IF (j.ge.3) rho = sqrt(0.6d0*rho**2.d0 + 0.3d0*rho_old**2.d0 + 0.1d0*rho_oold**2.d0)



  !   WRITE(name,'(A12,I3.3,A4)') 'Output/rho_w', j, '.dat'
  !   my_unit = j + 200
  !
  !
  !   OPEN(my_unit, FILE=name, STATUS='unknown')
  !   WRITE(my_unit,'(A5)') 'rho_w'
  !   DO i=1, n_pts
  !     WRITE(my_unit,'(ES16.4E2)') rho_w(i)
  !   END DO
  !   CLOSE(my_unit)

!     WRITE(name,'(A14,I3.3,A4)') 'Output/rho_HHe', j+1, '.dat'
!     my_unit = j + 200


!     OPEN(my_unit, FILE=name, STATUS='unknown')
!     WRITE(my_unit,'(A34)') 'rho_oliv  rho_w  rho_HHe  rho(2,i)'
!     DO i=1, n_pts
!        WRITE(my_unit,'(4ES16.4E2)') rho_oliv(i), rho_w(i), rho_HHe(i), rho(2,i)
!     END DO
!     CLOSE(my_unit)

 ! WRITE(6,*) 'flag2'

     !WRITE(6,*) intrf
     ! Current planet mass and radius
     M_Pc = mass_btw(1,intrf(n_lay),0)
     R_Pc = r(intrf(n_lay))

     ! Results
!     CALL display_log(j)
!     CALL write_profile(j)
     
!     IF (mod(j-1,rec)==0 .AND. plot_type1/=0) THEN
!        CALL plot_profile(j,plot_type1)
!        WRITE(6,*) '[i] Plotting profile'
!     END IF
     
!     WRITE(11,'(I4)',advance='no') j
!     DO k=1, n_lay
!        WRITE(11,'(ES27.15E2)',advance='no') r(intrf(k))/R_E
!     END DO
!     WRITE(11,*) ''
        
intrf_hist(:,j) = intrf
iter_num(j) = j

     ! Test for loop exit
     test = (abs(g_max-g_max_old) < conv_prec*g_max) .AND. (abs(P_max-P_max_old) < conv_prec*P_max) .AND. &
          (abs(T_max-T_max_old) < conv_prec*T_max) .AND. (abs(rho_max-rho_max_old) < conv_prec*rho_max)

     DO k=1, n_lay
        test = test .AND. (intrf(k)==intrf_old(k))
     END DO

     CALL check_convergence


     IF (test) THEN ! Case 1: convergence reached
        WRITE(6,*) '[i] Convergence reached.'
        EXIT
     ELSE IF (cnt_conv>=cnt_conv_max) THEN ! Case 2: detected oscillations around convergence
        WRITE(6,*) '[i] Too many oscillations detected: exiting iterations loop.'
        WRITE(6,*) ''
        !WRITE(6,*) '/!\ WARNING: convergence not reached.'
        EXIT
     ELSE IF (j>=j_max) THEN ! Case 3: reached maximum number of iterations
        WRITE(6,*) '[i] Allowed maximum number of iterations reached.'
        WRITE(6,*) ''
        !WRITE(6,*) '/!\ WARNING: convergence not reached.'
        EXIT
     END IF

     

     
  END DO




  ! FINAL COMPUTATION
  ! -----------------

  CALL mass_layers
!  CALL FeSi_ratio
!  CALL mom_inertia
  x_corec = M(1)/M_Pc
  x_H2Oc  = M(2)/M_Pc


  M_Pout = M_Pc/M_E
  x_coreout = x_corec
  x_H2Oout = x_H2Oc
  R_Pout = R_Pc/R_E
  FeSiout = FeSic
  rhoout = 3.d0*M_Pc/(4.d0*pi*R_Pc**3.d0)/1.d3

  OtoHout = OtoH_env
  metalout = metallicity

  rout = r
  interfaceout = intrf

  DO i = 1,n_pts

     gout(i) = g(layer(i),i)
     Tout(i) = T(layer(i),i)
     Pout(i) = P(layer(i),i)
     rhoarrout(i) = rho(layer(i),i)
 
  END DO



  !WRITE(6,*) 'g(layer(1199),1199) = ', g(layer(1199),1199)


  DO i = 1,n_pts

     gout(i) = g(layer(i),i)
     Tout(i) = T(layer(i),i)
     Pout(i) = P(layer(i),i)
     rhoarrout(i) = rho(layer(i),i)
 
  END DO


  ! Calculation of Cv
  ! Core = water + rock

  i_center = 1
  i_core = intrf(2)


  DO i = i_center,i_core
     sesame_out2 = grun_sesame(Tout(i),Pout(i))
     cv_rock_arr(i) = sesame_out2(2)
     entropy_rock(i) = sesame_out2(4)

     output_heat_w = heat_cap_water(Tout(i),rho_w_core(i))
     cv_h2o_arr(i) = output_heat_w(1)
     entropy_H2O(i) = output_heat_w(2)
 
  END DO

  cv_core =  0.5d0*cv_rock_arr + 0.5d0*cv_h2o_arr
  S_core  =  0.5d0*entropy_rock + 0.5d0*entropy_h2o


  ! Envelope = water + H/He
  i_totalpl = intrf(3)

  cv_h2o_arr = 0.d0

  Y_astk = 0.275d0
  Y_til = 0.245d0

  DO i = i_core,i_totalpl

     output_heat_w = heat_cap_water(Tout(i),rho_w_env(i))
     cv_h2o_arr(i) = output_heat_w(1)
     entropy_H2O(i) = output_heat_w(2)

     HHe_out2 = grun_HHe(Tout(i),rho_HHe(i))
     cv_HHe_arr(i) = HHe_out2(2)

     output_ch_S = U_chabrier(log10(rho_HHe(i)/1d3),log10(Tout(i)))
     S_cd21 = 1d6 * output_ch_S(3)                                                                                      ! Correction at Y = 0.24, J/kg/K
     Smix_value = interp2d_opt(log10(Tout(i)),log10(Pout(i))+10d0,logT_HG,logP_HG,Smix_HG,size(logT_HG),size(logP_HG))  ! Smix in erg/g/K
     entropy_HHe(i) = S_cd21 + Smix_value*1d-4 * ((1-Y_astk)*Y_astk - (1-Y_til)*Y_til)                                  ! Correction at Y = 0.275, J/kg/K
         
 
  END DO

  cv_env = Zenv*cv_h2o_arr + (1d0-Zenv)*cv_HHe_arr
  S_env  = Zenv*entropy_H2O + (1d0-Zenv)*entropy_HHe


  ! Final cv
  cv_all(i_center:i_core) = cv_core(i_center:i_core)
  cv_all(i_core:i_totalpl) = cv_env(i_core:i_totalpl)

  S_all(i_center:i_core) = S_core(i_center:i_core)
  S_all(i_core:i_totalpl) = S_env(i_core:i_totalpl)


  cvout = cv_all
  Sout = S_all




 !  OPEN(87, FILE='Output/gruneisen.dat', STATUS='unknown')
 !  WRITE(87,'(A100)',advance='no') 'r(i)  gam(1,i)  gam_oliv_c  gam_w_c  rho_oliv_c  rho_w_c'
 !  WRITE(87,'(A100)',advance='no') 'gam(2,i)  gam_oliv_e  gam_w_e  rho_oliv_e  rho_w_e  gam_hhe  rho_hhe  gam(layer(i),i)'
 !  WRITE(87,'(A100)')               'S_HHe  S_water  S_rock  S_all '
 !  DO i=1, n_pts
 !    WRITE(87,'(18ES16.4E2)') r(i), gam(1,i), gam_oliv_core(i), gam_w_core(i), rho_oliv_core(i), rho_w_core(i), &
 !                                   gam(2,i), gam_oliv_env(i),  gam_w_env(i), rho_oliv_env(i), rho_w_env(i), gam_HHe(i), &
 !                                   rho_HHe(i), gam(layer(i),i), entropy_HHe(i), entropy_H2O(i), entropy_rock(i), S_all(i)
 !  END DO
 !  CLOSE(87)



  ! DISPLAY RESULTS
  ! ---------------
  
!  IF (plot_type1/=0) CALL plot_profile(j,plot_type1)
!  IF (plot_type2/=0) CALL plot_interfaces(plot_type2)
!  IF (plot_type3/=0) CALL plot_ternary(plot_type3)
!  IF (x_H2O>0.d0.AND.plot_type4/=0) CALL plot_water_phase(plot_type4)
!  CALL display_results(j)
!  CALL plot_interior(1)
  
  
  ! CLOSING FILES
  ! -------------

!  IF (log_type==1) CLOSE(6)
!  CLOSE(11)


  ! END PROGRAM
  ! -----------

END SUBROUTINE GASTLI_interior_subroutine

END MODULE GASTLI_interior
