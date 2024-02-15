!      /''''''''''''''''''''''''\
!     /          modEOS          \
!     \........................../

! This module contains all equations of state (EOS) that can be used by the program. One can add more of
! them if needed. The EOS are all contained in one function, that chooses which EOS to use depending on
! the considered layer.

MODULE modEOS

  USE constants
  USE dimensions
  USE parameters
  USE functions
  USE mazevet
!  USE AQUA_EOS,only:DIM_OUTPUT,LOAD_TABLE_PT,INTERPOLATE_AQUA_PT
!  USE AQUA_EOS,only:LOAD_TABLE_RhoT,INTERPOLATE_AQUA_RHOT
!  USE AQUA_EOS,only:LOAD_TABLE_RhoU,INTERPOLATE_AQUA_RHOU

  IMPLICIT NONE

CONTAINS

  !--------------------------------!
  !                                !
  !          FUNCTION EOS          !
  !                                !
  !--------------------------------!
  !
  ! Function that chooses which EOS to use in the program, among the following ones:
  !      EOS #1: BM3 (Third order Birch-Murnaghan)
  !      EOS #2: MGD (Mie-Gruneisen-Debye)
  !      EOS #3: Vinet + Debye thermal correction
  !      EOS #4: HSvdw (Hard Sphere Vaan der Waals)
  !	 EOS #5: GE (Gibbs Energy)  
  !      EOS #6: LJ potential (Duan & Zhang 2006) for supercrit.
  !      EOS #7: H02 (Holzapfel simple) + Debye thermal correction
  !      EOS #8: H12 (Holzapfel accurate) + Debye thermal correction  
  !      EOS #9: IAPWS95 for supercrit.
  !      EOS #10: Mazevet et al. 2019 for supercrit.
  !      EOS #11: AQUA EOS (Haldemann et al. 2020)
  !      EOS #12: Chabrier EOS for H/He
  !      EOS #13: Mazevet+22 EOS for H/He
  ! /!\ Always verify that the number of EOS is equal to 'n_EOS' in 'dimensions.f95'.
  !
  FUNCTION EOS(rho1,T1_in,k,wEOS)

    !------

    ! Variables of the function
    DOUBLE PRECISION                     :: EOS       ! Result

    INTEGER,          INTENT(IN)         :: k,      & ! Layer index
                                            wEOS      ! Index of EOS to use

    DOUBLE PRECISION, INTENT(IN)         :: rho1,   & ! Density
                                            T1_in        ! Temperature

    ! Variables common to all EOS
    DOUBLE PRECISION                     :: x         ! Density to uncompressed density ratio

    ! Variables of Debye thermal correction
    INTEGER                              :: l         ! Integration index

    DOUBLE PRECISION                     :: gam1,   & ! Gruneisen parameter
                                            T_D       ! Debye temperature
    
    DOUBLE PRECISION, DIMENSION(n_int)   :: y,      & ! Variable of integration
                                            y_0,    & ! Variable of integration
                                            arr,    & ! Function to integrate
                                            arr_0     ! Function to integrate

    ! Parameter for BM3 EOS
    DOUBLE PRECISION                     :: K_T0      ! Temperature-corrected bulk modulus
                                            
    ! Correction of several EOS
    DOUBLE PRECISION, DIMENSION(n_lay,3) :: alphaC, & ! Correction on reference bulk modulus for several EOS [Pa]
                                            betaC     ! Correction on pressure derivative of bulk modulus for EOS
                                            
    ! Parameters for the H12 EOS
    DOUBLE PRECISION                     :: P_TF0,  & ! Pressure of a Fermi gas at ambiant conditions [Pa]
                                            C10,    & ! First coefficient of H12 equation
                                            C12       ! Second coefficient of H12 equation
           
    ! Coefficients for LJ potential EOS
    DOUBLE PRECISION                     :: a1_1, &
					    a2_1, &
                                            a3_1, &
                                            a4_1, &
					    a5_1, &
					    a6_1, &
					    a7_1, &
					    a8_1, &
					    a9_1, &
					    a10_1, &
					    a11_1, &
					    a12_1, &
					    alpha_lj_1, &
					    beta_lj_1, &
					    gamma_lj_1, &
 	                                    a1_2, &
					    a2_2, &
                                            a3_2, &
                                            a4_2, &
					    a5_2, &
					    a6_2, &
					    a7_2, &
					    a8_2, &
					    a9_2, &
					    a10_2, &
					    a11_2, &
					    a12_2, &
					    alpha_lj_2, &
					    beta_lj_2, &
					    gamma_lj_2 	

    ! Coefficients for EOS 4 HSVdw
    DOUBLE PRECISION 			:: aa1, &
					   aa2, &
					   aa3, &
					   aa4, &
					   aa5, &
					   aa6, &
					   aa7, &
					   aa8, &
					   aa9, &
					   aa10, &
					   aa11, &
					   aa12, &
					   alpha_eos, &
					   beta_eos, &
					   gam_eos
    
    ! Critical temperature and pressure for LJ potential EOS and Gibbs
    DOUBLE PRECISION                     :: Tc, &     ! Critical temperature of water [K]
					    Pc        ! Critical pressure of water [cm3.mol-1 or Pa]
    DOUBLE PRECISION :: Rp, &
                        Vc, &
                        Tre, &
                        F_pot, &
                        E_pot, &
                        D_pot, &
                        B_pot, &
                        C_pot, &
                        cf_Mmol, &
                        V, &
                        Z_pot, &
                        EOS_bar, &
                        vap_a, &
                        vap_b, &
                        aZ    

    ! Variables for IAPWS95 
    DOUBLE PRECISION   :: rhocrit, &     ! Critical rho [kg.m-3]    
                          tau, &
                          delta, &       ! Reduced temperature and density
                          Rcrit, &       ! R constant [kJ.kg-1.K-1]
                          phird, &       ! phi^r_delta for pressure
                          crochet3, &    ! intermediate
                          ddd_delta, &
                          dps_delta, &   ! intermediate
                          der_dd, &      ! intermediate
                          dd, &
                          ps, &          ! intermediate
                          th

    ! Variables for Mazevet et al. 2019
   DOUBLE PRECISION   :: PnkT,FNkT,UNkT,CV,CHIT,CHIR,PMbar,USPEC   !, &
                          !rho_gcc



    ! Variables for AQUA, Haldemann et al. 2020
    ! REAL(8)  :: AQUA_RES(DIM_OUTPUT)
    ! CHARACTER(LEN=200) :: path
    ! integer, save :: counter_eos = 0


    ! Variables for Chabrier
    DOUBLE PRECISION   :: rho_un
                          
    DOUBLE PRECISION, DIMENSION(7) :: output


    ! Variables for Mazevet+22
    DOUBLE PRECISION   :: press_GPa


    !------


    ! Coefficients for the correction of thermodynamical parameters on several EOS in the following density ranges:

    !    1. rho/rho_0(k) in [1:1.5] (Valencia et al. 2007)
 
    alphaC(1,1) = 0.d0
    alphaC(2,1) = 0.d0

    betaC(1,1) = 1.01d0
    betaC(2,1) = 1.d0

    !    2. rho/rho_0(k) in [1:5]

    alphaC(1,2) = -52.d9
    alphaC(2,2) =  0.d0

    betaC(1,2) = 1.25d0
    betaC(2,2) = 1.d0

    !    3. rho/rho_0(k) in [1:10]

    alphaC(1,3) = -94.d9
    alphaC(2,3) = 0.d0

    betaC(1,3) = 1.5d0
    betaC(2,3) = 1.d0

    !------

    ! Coefficients for LJ potential EOS (Duan & Zhang 2006, table 4)
    ! 0.0 - 0.2 GPa
    a1_1 = 4.38269941d-2
    a2_1 = -1.68244362d-1
    a3_1 = -2.36923373d-1
    a4_1 = 1.13027462d-2
    a5_1 = -7.67764181d-2
    a6_1 = 9.71820593d-2
    a7_1 = 6.62674916d-5
    a8_1 = 1.06637349d-3
    a9_1 = -1.23265258d-3 
    a10_1 = -8.93953948d-6
    a11_1 = -3.88124606d-5
    a12_1 = 5.61510206d-5
    alpha_lj_1 = 7.51274488d-3
    beta_lj_1 = 2.51598931d0
    gamma_lj_1 = 3.94000000d-2

    ! 0.2 - 10.0 GPa
    a1_2 = 4.68071541d-2
    a2_2 = -2.81275941d-1
    a3_2 = -2.43926365d-1
    a4_2 = 1.10016958d-2
    a5_2 = -3.86603525d-2
    a6_2 = 9.30095461d-2
    a7_2 = -1.15747171d-5
    a8_2 = 4.19873848d-4
    a9_2 = -5.82739501d-4 
    a10_2 = 1.00936000d-6
    a11_2 = -1.01713593d-5
    a12_2 = 1.63934213d-5
    alpha_lj_2 = -4.49505919d-2
    beta_lj_2 = -3.15028174d-1
    gamma_lj_2 = 1.25000000d-2

    !------


    ! Case #0: no EOS (void)
    IF (wEOS==0) THEN
       EOS = 0.d0
       
       RETURN

       
       
    ! EOS #1: BM3
    ELSE IF (wEOS==1) THEN
       IF (EOS_th(k)==0) THEN ! Without consideration of the thermal component
          x = rho1 / rho_0(k)
          K_T0 = K_0(k)
       ELSE IF (EOS_th(k)==1) THEN ! With consideration of the thermal component
          x = rho1 / rho_0(k) * exp(a_T(k)*(T_0(k)-T1_in) + 0.5d0*b_T(k)*(T_0(k)**2-T1_in**2) &
               + c_T(k)*(1.d0/T_0(k)-1.d0/T1_in))
          K_T0 = K_0(k) + a_P(k)*(T1_in-T_0(k))
       END IF

       EOS = 3.d0*K_T0/2.d0 * (x**(7.d0/3.d0)-x**(5.d0/3.d0)) &
            * (1.d0-3.d0/4.d0*(4.d0-Kp_0(k))*(x**(2.d0/3.d0)-1.d0))
       
     RETURN

       
       
    ! EOS #2: MGD   
    ELSE IF (wEOS==2) THEN
       x = rho1 / rho_0(k)
       
       ! Part 1: Mie-Gruneisen formula
       EOS = 3.d0*K_0(k)/2.d0 * (x**(7.d0/3.d0)-x**(5.d0/3.d0)) &
            * (1.d0-3.d0/4.d0*(4.d0-Kp_0(k))*(x**(2.d0/3.d0)-1.d0))
       
       ! Part 2: Debye thermal correction (if consideration of the thermal component is asked)
       IF (EOS_th(k)==1) THEN
          gam1 = gam_0(k) * x**(-q(k))
          T_D = T_D0(k) * x**gam1
          
          FORALL(l=1:n_int)
             y(l) = dble(l-1)*T_D/T1_in/dble(n_int-1)
             y_0(l) = dble(l-1)*T_D/T_0(k)/dble(n_int-1)
          END FORALL
          arr(1) = 0.d0
          arr_0(1) = 0.d0
          FORALL(l=2:n_int)
             arr(l) = y(l)**3.d0 / (exp(y(l))-1.d0)
             arr_0(l) = y_0(l)**3.d0 / (exp(y_0(l))-1.d0)
          END FORALL
          
          EOS = EOS + 9.d0*gam1*rho1*R_g/Mmol(k)/(T_D**3.d0) &
          !     * (T1_in**4.d0*integral(1,n_int,arr,y)-T_0(k)**4.d0*integral(1,n_int,arr_0,y_0))
           * (T1_in**4.d0*integrate(y(1:n_int),arr(1:n_int))-T_0(k)**4.d0*integrate(y_0(1:n_int),arr_0(1:n_int)))

       END IF

     RETURN

       

    ! EOS #3: Vinet
    ELSE IF (wEOS==3) THEN
       x = rho1 / rho_0(k)
       
       ! Part 1: Vinet formula
       IF (corEOS==0) THEN ! without correction
          EOS = 3.d0*K_0(k) * (x**(2.d0/3.d0)-x**(1.d0/3.d0)) &
               * exp(3.d0/2.d0*(Kp_0(k)-1.d0)*(1.d0-x**(-1.d0/3.d0)))
       ELSE ! with correction
          EOS = 3.d0*(K_0(k)+alphaC(k,corEOS)) * (x**(2.d0/3.d0)-x**(1.d0/3.d0)) &
               * exp(3.d0/2.d0*(Kp_0(k)*betaC(k,corEOS)-1.d0)*(1.d0-x**(-1.d0/3.d0)))
       END IF
       
       ! Part 2: Debye thermal correction (if consideration of the thermal component is asked)
       IF (EOS_th(k)==1) THEN
          gam1 = gam_0(k) * x**(-q(k))
          T_D = T_D0(k) * x**gam1
          
          FORALL(l=1:n_int)
             y(l) = dble(l-1)*T_D/T1_in/dble(n_int-1)
             y_0(l) = dble(l-1)*T_D/T_0(k)/dble(n_int-1)
          END FORALL
          arr(1) = 0.d0
          arr_0(1) = 0.d0
          FORALL(l=2:n_int)
             arr(l) = y(l)**3.d0 / (exp(y(l))-1.d0)
             arr_0(l) = y_0(l)**3.d0 / (exp(y_0(l))-1.d0)
          END FORALL
          
          EOS = EOS + 9.d0*gam1*rho1*R_g/Mmol(k)/(T_D**3.d0) &
          !     * (T1_in**4.d0*integral(1,n_int,arr,y)-T_0(k)**4.d0*integral(1,n_int,arr_0,y_0))
          !* (T1_in**4.d0*integrate(y(1:n_int),arr(1:n_int))-T_0(k)**4.d0*integrate(y_0(1:n_int),arr_0(1:n_int)))
          * (T1_in**4.d0*(1+1)-T_0(k)**4.d0*(1+1))

       END IF

    RETURN



    ! EOS #4: HSvdw
    ELSE IF (wEOS==4) THEN 

       aa1 = 4.38269941e-02
       aa2 = -1.68244362e-01
       aa3 = -2.36923373e-01
       aa4 = 1.13027462e-02
       aa5 = -7.67764181e-02
       aa6 = 9.71820593e-02
       aa7 = 6.62674916e-05
       aa8 = 1.06637349e-03
       aa9 = -1.23265258e-03
       aa10 = -8.93953948e-06
       aa11 = -3.88124606e-05
       aa12 = 5.61510206e-05
       alpha_eos = 7.51274488e-03
       beta_eos = 2.51598931e+00
       gam_eos = 3.94e-02
    
       Tre = T1_in / 647.15
	        
       V = (0.01801528 / rho1) * 1e+6
       Vc = 83.14467 * 647.25 / 221.19
	        
       aZ = 1 + (aa1 + aa2/(Tre**2) + aa3/(Tre**3))*(Vc/V) + (aa4 + aa5/(Tre**2) + aa6/(Tre**3))*(Vc/V)**2 &
          + (aa7 + aa8/(Tre**2) + aa9/(Tre**3))*(Vc/V)**4 + (aa10 + aa11/(Tre**2) + aa12/(Tre**3))*(Vc/V)**5 &
          + ((alpha_eos/(Tre**3))*(Vc/V)**2)*(beta_eos + gam_eos *(Vc/V)**2)*exp(- gam_eos * (Vc/V)**2)   
    
       EOS = aZ * 83.14467 * T1_in / V
       EOS = EOS * 1e+05
    
    RETURN


    ! EOS #5: GE - Simplified water vapor
    ELSE IF (wEOS==5) THEN

       cf_Mmol = 18.01528d-3
       V = (cf_Mmol/rho1)   

       ! Critical temperature and pressure (Suarez-Arriaga 2014). K and Pa
       Tc = 647.27d0
       Pc = 221.15d5

       ! Gas constant in J.mol-1.K-1
       Rp = 8.314472d0

       vap_a = (27.d0/64.d0) * Rp**2 * Tc**2/Pc
       vap_b = Rp * Tc/(8*Pc)

       Z_pot = V/(V-vap_b) - vap_a/(Rp*T1_in*V)

       EOS = Z_pot * Rp * T1_in/V

       RETURN



       ! EOS #6: LJ potential
       ELSE IF (wEOS==6) THEN

       ! Critical temperature and pressure (Duan & Zhang 2006). K and cm3.mol-1 (or bar)
       Tc = 647.25d0
       Pc = 221.19d0
	  
       ! Gas constant in cm3.bar.mol-1.K-1
       Rp = 83.14467d0


       ! Parameters for EOS
       Vc = Rp*Tc/Pc
       Tre = T1_in/Tc
       F_pot = alpha_lj_2/Tre**3
       E_pot = a10_2 + a11_2/Tre**2 + a12_2/Tre**3
       D_pot = a7_2 + a8_2/Tre**2 + a9_2/Tre**3
       C_pot = a4_2 + a5_2/Tre**2 + a6_2/Tre**3
       B_pot = a1_2 + a2_2/Tre**2 + a3_2/Tre**3

       ! Volume in cm3
       cf_Mmol = 18.01528d-3
       V = (cf_Mmol/rho1)*1d6

      ! Compressibility factor
      Z_pot = 1.d0 + B_pot*Vc/V + C_pot*(Vc/V)**2.d0 + D_pot*(Vc/V)**4.d0 + E_pot*(Vc/V)**5.d0 + F_pot*(Vc/V)**2.d0 &
            * (beta_lj_2 + gamma_lj_2*(Vc/V)**2.d0) * exp(-gamma_lj_2*(Vc/V)**2.d0)

      EOS_bar = Z_pot*Rp*T1_in/V


      ! Pressure in Pa
      !EOS = EOS_bar*1.d5
      !RETURN

      IF (EOS_bar*1.d5 > 2.d8) THEN
          EOS = EOS_bar*1.d5
          RETURN
      ELSE

       ! Parameters for EOS
       Vc = Rp*Tc/Pc
       Tre = T1_in/Tc
       F_pot = alpha_lj_1/Tre**3
       E_pot = a10_1 + a11_1/Tre**2 + a12_1/Tre**3
       D_pot = a7_1 + a8_1/Tre**2 + a9_1/Tre**3
       C_pot = a4_1 + a5_1/Tre**2 + a6_1/Tre**3
       B_pot = a1_1 + a2_1/Tre**2 + a3_1/Tre**3

       ! Volume in cm3
       cf_Mmol = 18.01528d-3
       V = (cf_Mmol/rho1)*1d6

      ! Compressibility factor
      Z_pot = 1.d0 + B_pot*Vc/V + C_pot*(Vc/V)**2.d0 + D_pot*(Vc/V)**4.d0 + E_pot*(Vc/V)**5.d0 + F_pot*(Vc/V)**2.d0 &
           * (beta_lj_1 + gamma_lj_1*(Vc/V)**2.d0) * exp(-gamma_lj_1*(Vc/V)**2.d0)

      !EOS_bar = Z_pot*Rp*T1_in/V

      ! Pressure in Pa
      EOS = EOS_bar*1d5

      RETURN
    END IF              
	        

    ! EOS #7: H02
    ELSE IF (wEOS==7) THEN
       x = rho1 / rho_0(k)

       ! Part 1: Holzapfel 02 formula
       IF (corEOS==0) THEN ! without correction
          EOS = 3.d0*K_0(k) * (x**(5.d0/3.d0)-x**(4.d0/3.d0)) &
               * exp(3.d0/2.d0*(Kp_0(k)-3.d0)*(1.d0-x**(-1.d0/3.d0)))
       ELSE ! with correction
          EOS = 3.d0*(K_0(k)+alphaC(k,corEOS)) * (x**(5.d0/3.d0)-x**(4.d0/3.d0)) &
               * exp(3.d0/2.d0*((Kp_0(k)*betaC(k,corEOS))-3.d0)*(1.d0-x**(-1.d0/3.d0)))
       END IF
       
       ! Part 2: Debye thermal correction (if consideration of the thermal component is asked)
       IF (EOS_th(k)==1) THEN
          gam1 = gam_0(k) * x**(-q(k))
          T_D = T_D0(k) * x**gam1
          
          FORALL(l=1:n_int)
             y(l) = dble(l-1)*T_D/T1_in/dble(n_int-1)
             y_0(l) = dble(l-1)*T_D/T_0(k)/dble(n_int-1)
          END FORALL
          arr(1) = 0.d0
          arr_0(1) = 0.d0
          FORALL(l=2:n_int)
             arr(l) = y(l)**3.d0 / (exp(y(l))-1.d0)
             arr_0(l) = y_0(l)**3.d0 / (exp(y_0(l))-1.d0)
          END FORALL
          
          EOS = EOS + 9.d0*gam1*rho1*R_g/Mmol(k)/(T_D**3.d0) &
          !     * (T1_in**4.d0*integral(1,n_int,arr,y)-T_0(k)**4.d0*integral(1,n_int,arr_0,y_0))
          * (T1_in**4.d0*integrate(y(1:n_int),arr(1:n_int))-T_0(k)**4.d0*integrate(y_0(1:n_int),arr_0(1:n_int)))

       END IF

       RETURN

       

    ! EOS #8: H12
    ELSE IF (wEOS==8) THEN
       x = rho1 / rho_0(k)

       ! Coefficients for H12 formula
       P_TF0 = 2.337d-38 * (N_A*rho_0(k)*Z_eff(k)/Mmol(k))**(5.d0/3.d0)
       IF (corEOS==0) THEN ! without correction
          C10 = log(P_TF0/(3.d0*K_0(k)))
          C12 = 1.d0/2.d0 * (3.d0*Kp_0(k) - 2.d0*C10 - 9.d0)
       ELSE ! with correction
          C10 = log(P_TF0/(3.d0*(K_0(k)+alphaC(k,corEOS))))
          C12 = 1.d0/2.d0 * (3.d0*(Kp_0(k)*betaC(k,corEOS)) - 2.d0*C10 - 9.d0)
       END IF

       ! Part 1: Holzapfel 12 formula
       IF (corEOS==0) THEN ! without correction
          EOS = 3.d0*K_0(k) * (x**(5.d0/3.d0)-x**(4.d0/3.d0)) &
               * exp((C10+C12)*(1.d0-x**(-1.d0/3.d0)) - C12*(1.d0-x**(-1.d0/3.d0))**2.d0)
       ELSE ! with correction
          EOS = 3.d0*(K_0(k)+alphaC(k,corEOS)) * (x**(5.d0/3.d0)-x**(4.d0/3.d0)) &
               * exp((C10+C12)*(1.d0-x**(-1.d0/3.d0)) - C12*(1.d0-x**(-1.d0/3.d0))**2.d0)
       END IF
       
       ! Part 2: Debye thermal correction (if consideration of the thermal component is asked)
       IF (EOS_th(k)==1) THEN
          gam1 = gam_0(k) * x**(-q(k))
          T_D = T_D0(k) * x**gam1
          
          FORALL(l=1:n_int)
             y(l) = dble(l-1)*T_D/T1_in/dble(n_int-1)
             y_0(l) = dble(l-1)*T_D/T_0(k)/dble(n_int-1)
          END FORALL
          arr(1) = 0.d0
          arr_0(1) = 0.d0
          FORALL(l=2:n_int)
             arr(l) = y(l)**3.d0 / (exp(y(l))-1.d0)
             arr_0(l) = y_0(l)**3.d0 / (exp(y_0(l))-1.d0)
          END FORALL
          
          EOS = EOS + 9.d0*gam1*rho1*R_g/Mmol(k)/(T_D**3.d0) &
          !     * (T1_in**4.d0*integral(1,n_int,arr,y)-T_0(k)**4.d0*integral(1,n_int,arr_0,y_0))
          * (T1_in**4.d0*integrate(y(1:n_int),arr(1:n_int))-T_0(k)**4.d0*integrate(y_0(1:n_int),arr_0(1:n_int)))

       END IF

       RETURN

    ! EOS IAPWS95
    ELSE IF (wEOS==9) THEN	
    
    rhocrit = 322.0d0 
    Rcrit = 0.46151805d0 

    IF (T1_in < 1.d0) THEN
        tau = 647.096d0/1.d0
    ELSE 
        tau = 647.096d0/T1_in
    END IF
       
    
    delta = rho1/rhocrit

    CALL read_coefficients
 
     ! phird

     phird = 0.d0


     ! First term
     DO l=1,7
         phird = phird + nr(l) * dr(l) * delta**(dr(l)-1.d0) * tau**tr(l)
     END DO
     
     ! Second term
     DO l=8,51
         phird = phird + nr(l) * exp(-delta**cr(l)) & 
               * (delta**(dr(l)-1.d0)*tau**tr(l)*(dr(l)-cr(l)*delta**cr(l)))
     END DO

     ! Third term
     DO l=52,54
         crochet3 = dr(l)/delta - 2*alpha(l)*(delta-epsilon(l))
         phird = phird + nr(l) * delta**dr(l) * tau**tr(l) &
               * exp(-alpha(l)*(delta-epsilon(l))**2.d0-beta(l)*(tau-gamma(l))**2.d0) * crochet3
     END DO

     ! Fourth term
     DO l=55,56
         ps = exp(-C(l)*(delta-1.d0)**2.d0 - D(l)*(tau-1.d0)**2.d0)
         th = (1.d0-tau) + A(l)*((delta-1.d0)**2.d0)**(1/(2*beta(l)))
         dd = th**2.d0 + B(l)*((delta-1.d0)**2.d0)**ar(l)
         dps_delta = - 2.d0 * C(l) * (delta-1.d0) * ps
         der_dd = (delta-1.d0) * (A(l)*th*(2/beta(l))*((delta-1.d0)**2.d0)**((1/(2*beta(l)))-1) &
              + 2*B(l)*ar(l)*((delta-1.d0)**2.d0)**(ar(l)-1.d0))
         ddd_delta = br(l) * dd**(br(l)-1.d0) * der_dd
         phird = phird + nr(l) * (dd**br(l)*(ps+delta*dps_delta)+ddd_delta*delta*ps)
     END DO

     
     EOS = 1d3*(1 + delta*phird)*rho1*Rcrit*T1_in
     RETURN

    ! EOS Mazevet et al. 2019
    ELSE IF (wEOS==10) THEN

    !rho_gcc = rho1/1000d0 ! Density in g/cc


    CALL h2ofit(rho1/1000d0,T1_in,PnkT,FNkT,UNkT,CV,CHIT,CHIR,PMbar,USPEC)

    EOS = PMbar*1d11  ! Pressure in Pa
    RETURN


    ! AQUA EOS (Haldemann et al. 2020)
!    ELSE IF (wEOS==11) THEN
!    
!    IF (counter_eos==0) THEN
!
!        path = 'Tables'
!        CALL LOAD_TABLE_PT(path)
!        CALL LOAD_TABLE_RhoT(path)
!        CALL LOAD_TABLE_RhoU(path)
!        counter_eos = counter_eos + 1
!
!     END IF
!
!     !IF ((T1_in>150d0).AND.(T1_in<1d5).AND.(rho1>0d0).AND.(rho1<1d4)) THEN
!
!        AQUA_RES = INTERPOLATE_AQUA_RHOT(max(1d-1, rho1),max(151d0,T1_in))
!        !WRITE(6,*) 'AQUA_RES(1)'
!        !WRITE(6,*) AQUA_RES(1)
!        EOS = AQUA_RES(1) 
!        RETURN
!
!     !ELSE
!        !WRITE(6,*) 'Out of table'
!        !CALL h2ofit(rho1/1000d0,T1_in,PnkT,FNkT,UNkT,CV,CHIT,CHIR,PMbar,USPEC)
!        !WRITE(6,*) 'EOS'
!        !WRITE(6,*) PMbar*1d11
!        !EOS = PMbar*1d11  ! Pressure in Pa
!        !RETURN
!     !END IF
           
    ! Chabrier EOS for H/He
    ELSE IF (wEOS==12) THEN

 
    rho_un = rho1*1e-3           ! Density in gcc
    output = U_chabrier(log10(rho_un),log10(T1_in))      

    EOS = 1d9*output(1)          ! Pressure in Pa
    RETURN
           
    ! Mazevet+22 EOS for H/He
 !   ELSE IF (wEOS==13) THEN
 !
 !
 !   rho_un = rho1*1e-3           ! Density in gcc
 !   press_GPa = interp2d_opt(rho_un,T1_in,rho_maz_HHe,temp_maz_HHe,press_maz_HHe,size(rho_maz_HHe),size(temp_maz_HHe))    
 !
 !
 !   EOS = 1d9*press_GPa          ! Pressure in Pa
 !   RETURN


    ! Wrong choice of EOS
    ELSE
       WRITE(6,*) '/!\ ERROR in function "EOS": the asked EOS does not exist.'
       STOP

       
    END IF


    !------

  END FUNCTION EOS

END MODULE modEOS
