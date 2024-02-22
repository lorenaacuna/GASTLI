!      /'''''''''''''''''''''''''''''\
!     /          SUBROUTINES          \
!     \.............................../

! This module contains subroutines designed for calculation.

MODULE subroutines

  USE constants
  USE dimensions
  USE parameters
  USE functions
  USE modEOS
  USE OMP_LIB

  IMPLICIT NONE

CONTAINS

  !---------------------------------------!
  !                                       !
  !          SUBROUTINE layering          !
  !                                       !
  !---------------------------------------!
  !
  ! Assigns a layer to each spatial index using the position of layer interfaces.
  !
  SUBROUTINE layering

    !------

    INTEGER :: i, & ! Spatial index
               k    ! Layer index

    !------
 
    !$OMP PARALLEL DO
    DO k=1, n_lay
       FORALL(i = intrf(k):intrf(k+1)-1) layer(i)=k
    END DO
    !$OMP END PARALLEL DO

    !------
    
  END SUBROUTINE layering



  !------------------------------------------------!
  !                                                !
  !          SUBROUTINE density_limit_EOS          !
  !                                                !
  !------------------------------------------------!
  !
  ! Computes the lower limit in density of the chosen EOS for all layers.
  !
  SUBROUTINE density_limit_EOS

    !------

    INTEGER :: k ! Layer index

    !------

    DO k=1, n_lay
       rho_lim_EOS(k) = 0.d0

       IF (EOS_th(k)==1) THEN
          DO WHILE (isnan(EOS(rho_lim_EOS(k),T_0(k),k,iEOS(k))))
             rho_lim_EOS(k) = rho_lim_EOS(k) + 10.d0
             
             IF (rho_lim_EOS(k)>rho_0(k)) THEN
                WRITE(6,*) '/!\ ERROR in subroutine "density_limit_EOS": the chosen EOS cannot be used for this layer.'
                WRITE(6,*) '     k=', k, 'EOS=', iEOS(k)
                STOP
             END IF
          END DO
       END IF
    END DO
    
    !------
    
  END SUBROUTINE density_limit_EOS



  !----------------------------------------------!
  !                                              !
  !          SUBROUTINE chk_water_phase          !
  !                                              !
  !----------------------------------------------!
  !
  ! Checks to which water phase the planet's surface conditions (T,P) correspond.
  !
  SUBROUTINE chk_water_phase

    !------

    IF (T_surf<=Ttp(2)) THEN
       IF (P_surf<=P_water(1,T_surf))                              water_phase = 2 ! Vapor
       IF (P_surf>P_water(1,T_surf).AND.P_surf<=Ptp(2))            water_phase = 4 ! Ih
       IF (P_surf>Ptp(2))                                          water_phase = 6 ! Other
    ELSE IF (T_surf>Ttp(2).AND.T_surf<=Ttp(3)) THEN
       IF (P_surf<=P_water(1,T_surf))                              water_phase = 2 ! Vapor
       IF (P_surf>P_water(1,T_surf).AND.P_surf<=P_water(2,T_surf)) water_phase = 4 ! Ih
       IF (P_surf>P_water(2,T_surf).AND.P_surf<=P_water(3,T_surf)) water_phase = 3 ! Liquid
       IF (P_surf>P_water(3,T_surf))                               water_phase = 6 ! Other
    ELSE IF (T_surf>Ttp(3).AND.T_surf<=Ttp(1)) THEN
       IF (P_surf<=P_water(1,T_surf))                              water_phase = 2 ! Vapor
       IF (P_surf>P_water(1,T_surf).AND.P_surf<=P_water(2,T_surf)) water_phase = 4 ! Ih
       IF (P_surf>P_water(2,T_surf).AND.P_surf<=P_water(4,T_surf)) water_phase = 3 ! Liquid
       IF (P_surf>P_water(4,T_surf))                               water_phase = 6 ! Other
    ELSE IF (T_surf>Ttp(1).AND.T_surf<=Ttp(4)) THEN
       IF (P_surf<=P_water(7,T_surf))                              water_phase = 2 ! Vapor
       IF (P_surf>P_water(7,T_surf).AND.P_surf<=P_water(4,T_surf)) water_phase = 3 ! Liquid
       IF (P_surf>P_water(4,T_surf))                               water_phase = 6 ! Other
    ELSE IF (T_surf>Ttp(4).AND.T_surf<=Ttp(5)) THEN
       IF (P_surf<=P_water(7,T_surf))                              water_phase = 2 ! Vapor
       IF (P_surf>P_water(7,T_surf).AND.P_surf<=P_water(5,T_surf)) water_phase = 3 ! Liquid
       IF (P_surf>P_water(5,T_surf))                               water_phase = 6 ! Other
    ELSE IF (T_surf>Ttp(5).AND.T_surf<=Ttp(6)) THEN
       IF (P_surf<=P_water(7,T_surf))                              water_phase = 2 ! Vapor
       IF (P_surf>P_water(7,T_surf).AND.P_surf<=P_water(6,T_surf)) water_phase = 3 ! Liquid
       IF (P_surf>P_water(6,T_surf))                               water_phase = 5 ! VII
    ELSE IF (T_surf>Ttp(6)) THEN
       IF (P_surf<=Ptp(6))                                         water_phase = 2 ! Vapor
       IF (P_surf>Ptp(6).AND.P_surf<=P_water(6,T_surf))            water_phase = 1 ! Supercritical
       IF (P_surf>P_water(6,T_surf))                               water_phase = 5 ! VII
    END IF

    !IF (water_phase==1) THEN
       !WRITE(6,*) '/!\ ERROR in subroutine "chk_water_phase": supercritical water is not included.'
       !STOP
    IF (water_phase==2) THEN
       WRITE(6,*) '/!\ ERROR in subroutine "chk_water_phase": water vapor is not included.'
       !STOP
    ELSE IF (water_phase==4) THEN
       WRITE(6,*) '/!\ ERROR in subroutine "chk_water_phase": water ice Ih is not included in this version.'
       STOP
    ELSE IF (water_phase==6) THEN
       WRITE(6,*) '/!\ ERROR in subroutine "chk_water_phase": water ices other than Ih and VII are not included.'
       STOP
    ELSE IF (water_phase==3) THEN
       WRITE(6,*) '/!\ ERROR in subroutine "chk_water_phase": liquid water is not included in this version.'
       STOP
    ELSE IF (water_phase==5) THEN
       WRITE(6,*) '/!\ ERROR in subroutine "chk_water_phase": VII ice is not included in this version.'
       STOP
    END IF
    
    !------
    
  END SUBROUTINE chk_water_phase



  !------------------------------------!
  !                                    !
  !          SUBROUTINE new_g          !
  !                                    !
  !------------------------------------!
  !
  ! Computes the new gravity acceleration in each layer using the Gauss theorem.
  !
  SUBROUTINE new_g

    !------

    INTEGER                            :: i, & ! Spatial index
                                          k    ! Layer index

    DOUBLE PRECISION, DIMENSION(n_pts) :: arr  ! Temporary array

    !------

    ! Calculation
    DO k=1, n_lay
       arr = rho(k,:)*r**2
       !WRITE(6,*) 'rho(k,1199) = ', rho(k,1199)
       
       IF(k==1) THEN
          !$OMP PARALLEL DO
          DO i=2, n_pts
             !g(k,i) = g(k,intrf(k))*(r(intrf(k))/r(i))**2 + 4.d0*pi*G_grav*integral(intrf(k),i,arr,r)/r(i)**2
             g(k,i) = g(k,intrf(k))*(r(intrf(k))/r(i))**2 + 4.d0*pi*G_grav*integrate(r(intrf(k):i),arr(intrf(k):i))/r(i)**2
             !g(k,i) = g(k,intrf(k))*(r(intrf(k))/r(i))**2 + 4.d0*pi*G_grav*1.5d23/r(i)**2
             !WRITE(6,*) 'g(k,i)', integrate(r(intrf(k):i),arr(intrf(k):i))

          END DO
          !$OMP END PARALLEL DO
       ELSE
          !$OMP PARALLEL DO
          DO i=2, n_pts
             !g(k,i) = g(k-1,intrf(k))*(r(intrf(k))/r(i))**2 + 4.d0*pi*G_grav*integral(intrf(k),i,arr,r)/r(i)**2
             g(k,i) = g(k-1,intrf(k))*(r(intrf(k))/r(i))**2 + 4.d0*pi*G_grav*integrate(r(intrf(k):i),arr(intrf(k):i))/r(i)**2
             !g(k,i) = g(k-1,intrf(k))*(r(intrf(k))/r(i))**2 + 4.d0*pi*G_grav*1.5d23/r(i)**2
             !WRITE(6,*) 'g(k,i)', integrate(r(intrf(k):i),arr(intrf(k):i))
          END DO
          !$OMP END PARALLEL DO
       END IF
    END DO

    ! Save and update max values
    g_max_old = g_max
    g_max = 0.d0
    
    !$OMP PARALLEL DO
    DO i=1, n_pts
       g_max  = max(g_max, g(layer(i),i))
       g_mmax = max(g_mmax,g(layer(i),i))
    END DO
    !$OMP END PARALLEL DO

    !------
    
  END SUBROUTINE new_g



  !------------------------------------!
  !                                    !
  !          SUBROUTINE new_P          !
  !                                    !
  !------------------------------------!
  !
  ! Computes the new pressure in each layer using the hydrostatic equilibrium.
  !
  SUBROUTINE new_P

    !------

    INTEGER                            :: i,   & ! Spatial index
                                          k      ! Layer index


    INTEGER,          DIMENSION(n_lay) :: isPOOR ! Marker of pressure out of EOS validity range

    DOUBLE PRECISION, DIMENSION(n_pts) :: arr    ! Temporary array

    !------

    ! Initialization
    isPOOR = 0

    ! Calculation
    DO k=n_lay-1, 1, -1
       arr = rho(k,:)*g(k,:)
       
       !$OMP PARALLEL 
       !WRITE(6,*) 'Max number of threads = ', OMP_GET_MAX_THREADS()
       !WRITE(6,*) 'Thread number =', OMP_GET_THREAD_NUM()
       DO i=n_pts, 1, -1
          !P(k,i) = P(k+1,intrf(k+1)) + integral(i,intrf(k+1),arr,r)
          P(k,i) = P(k+1,intrf(k+1)) + integrate(r(i:intrf(k+1)),arr(i:intrf(k+1)))
          !P(k,i) = P(k+1,intrf(k+1)) + 1.5d12
          !WRITE(6,*) 'P(k,i) ', integrate(r(i:intrf(k+1)),arr(i:intrf(k+1)))
       END DO
       !$OMP END PARALLEL 

       ! Test for EOS validity range
       IF (P(k,intrf(k)) > EOS_lim_P(iEOS(k))) isPOOR(k)=k
    END DO

    IF (chk_EOS==1 .AND. ANY(isPOOR/=0)) THEN
       WRITE(6,*) '/!\ WARNING: the pressure is out of the EOS validity range'
       WRITE(6,*) '     k=', pack(isPOOR, isPOOR/=0)
    END IF

    ! Save and update max values
    P_max_old = P_max
    P_max = 0.d0
    
    !$OMP PARALLEL DO
    DO i=1, n_pts
       P_max  = max(P_max, P(layer(i),i))
       P_mmax = max(P_mmax,P(layer(i),i))
    END DO
    !$OMP END PARALLEL DO

    !------
        
  END SUBROUTINE new_P


  
  !------------------------------------!
  !                                    !
  !          SUBROUTINE new_T          !
  !                                    !
  !------------------------------------!
  !
  ! Computes the new temperature in each layer using the adiabatic assumption.
  !
  SUBROUTINE new_T

    !------

    INTEGER                            :: i,  & ! Spatial index
                                          k     ! Layer index

    !DOUBLE PRECISION, DIMENSION(n_pts) :: arr, gam_oliv, gam_w, gam_HHe  ! Array 

    DOUBLE PRECISION, DIMENSION(n_pts) :: arr


    DOUBLE PRECISION :: inv_gam

    DOUBLE PRECISION, DIMENSION(3) :: HHe_output

    DOUBLE PRECISION, DIMENSION(4) :: sesame_output


 
    !------

    ! Initialization
    arr   = 0.d0
    
    ! Calculation
    
    FORALL(k=1:n_lay, i=1:n_pts) gam(k,i) = gam_0(k) * (rho_0(k)/rho(k,i))**q(k)

 !   ! Consistent calculation of the Gruneisen parameter for pure olivine
 !   DO i=1,n_pts 
 !        gam(1,i) = gam_0(1) * (rho_0(1)/rho_oliv(i))**q(1)
 !   END DO

   !WRITE(6,*) 'HELLO T'


    ! Layer 1: core

    !$OMP PARALLEL DO PRIVATE(inv_gam)
    DO i=1,n_pts 
         ! Olivine
         ! Vinet EOS
         !!!!gam_oliv_core(i) = gam_0(1) * (rho_0(1)/rho_oliv_core(i))**q(1)
         ! SESAME EOS
         sesame_output = grun_sesame(T(1,i),P(1,i))
         gam_oliv_core(i) = sesame_output(1)


         ! Water
         gam_w_core(i) = grun_mazevet(T(1,i),rho_w_core(i))
         !gam_w(i) = grun_AQUA(T(1,i),rho_w(i))

         inv_gam = 0.5d0/gam_w_core(i) + 0.5d0/gam_oliv_core(i)
         gam(1,i) = 1.d0/inv_gam
    END DO
    !$OMP END PARALLEL DO


   !WRITE(6,*) 'CORE T OK'

    ! Layer 2: envelope

    !$OMP PARALLEL DO PRIVATE(inv_gam)
    DO i=1,n_pts 
         ! Olivine  
         !gam_oliv_env(i) = gam_0(2) * (rho_0(2)/rho_oliv_env(i))**q(2)
     

         ! H/He
         HHe_output = grun_HHe(T(2,i),rho_HHe(i))
         gam_HHe(i) = HHe_output(1)
         

         ! Water
         gam_w_env(i) = grun_mazevet(T(2,i),rho_w_env(i))
         !gam_w(i) = grun_AQUA(T(2,i),rho_w(i))


         ! Mixture
         ! Metals represented by water
         inv_gam = (1-Zenv)/gam_HHe(i) + Zenv/gam_w_env(i)

         ! Metals represented by water+olivine
         ! inv_gam = (1-Zenv)/gam_HHe(i) + Zenv*0.5d0/gam_w_env(i) + Zenv*0.5d0/gam_oliv_env(i)

         gam(2,i) = 1.d0/inv_gam

    END DO
    !$OMP END PARALLEL DO



   !WRITE(6,*) 'ENVELOPE T OK'



    DO k=n_lay-1, 1, -1
       arr(1) = gam(k,1)*g(k,1) * (rho(k,2)-rho(k,1))/(P(k,2)-P(k,1))
       arr(n_pts) = gam(k,n_pts)*g(k,n_pts) * (rho(k,n_pts)-rho(k,n_pts-1))/(P(k,n_pts)-P(k,n_pts-1))
       FORALL(i=2:n_pts-1) arr(i) = gam(k,i)*g(k,i) * (rho(k,i+1)-rho(k,i-1))/(P(k,i+1)-P(k,i-1))
       WHERE(isnan(arr)) arr = 0.d0
       
       !$OMP PARALLEL DO
       DO i=n_pts, 1, -1
          !T(k,i) = (del_T(k)+T(k+1,intrf(k+1))) * exp(integral(i,intrf(k+1),arr,r))   
          T(k,i) = (del_T(k)+T(k+1,intrf(k+1))) * exp(integrate(r(i:intrf(k+1)),arr(i:intrf(k+1))))  
          !T(k,i) = (del_T(k)+T(k+1,intrf(k+1))) * exp(0d0) 
          !WRITE(6,*) 'T(k,i) ', integrate(r(i:intrf(k+1)),arr(i:intrf(k+1)))
 
          T(k,i) = max(T(k,i),T_surf)

       END DO
       !$OMP END PARALLEL DO
    END DO

    !write(6,*) T(2,:)

    ! Save and update max values
    T_max_old = T_max
    T_max = 0.d0

    !$OMP PARALLEL DO
    DO i=1, n_pts
       T_max  = max(T_max, T(layer(i),i))
       T_mmax = max(T_mmax,T(layer(i),i))
    END DO
    !$OMP END PARALLEL DO


!   OPEN(87, FILE='Output/gruneisen.dat', STATUS='unknown')
!   WRITE(87,'(A100)') 'P(1,i)  T(1,i)  rho(1,i)  r(i)  gamma(i)  gam_oliv  gam_w  rho_oliv  rho_w'
!   DO i=1, n_pts
!     WRITE(87,'(9ES16.4E2)') P(1,i), T(1,i), rho(1,i), r(i), gam(1,i), gam_oliv(i), gam_w(i), rho_oliv(i), rho_w(i)
!   END DO
!   CLOSE(87)
!
!   OPEN(87, FILE='Output/gruneisen_layer2.dat', STATUS='unknown')
!   WRITE(87,'(A117)') 'P(2,i)  T(2,i)  rho(2,i)  r(i)  gamma(i)  gam_oliv  gam_HHe  rho_oliv  rho_HHe  gam(layer(i),i)'
!   DO i=1, n_pts
!     WRITE(87,'(10ES16.4E2)') P(2,i), T(2,i), rho(2,i), r(i), gam(2,i), gam_oliv(i), &
!           gam_HHe(i), rho_oliv(i), rho_HHe(i), gam(layer(i),i) 
!   END DO
!   CLOSE(87)



    !------

  END SUBROUTINE new_T



  !--------------------------------------!
  !                                      !
  !          SUBROUTINE Newton           !
  !                                      !
  !--------------------------------------!
  !
  ! Newton method to compute new density
  !
  SUBROUTINE Newton(k, Pin, Tin, EOSid, last, rhoout)

    !------

    INTEGER, INTENT(IN) :: k, EOSid

    DOUBLE PRECISION, DIMENSION(n_pts), INTENT(IN) :: Pin, Tin

    DOUBLE PRECISION, DIMENSION(n_pts), INTENT(OUT) :: rhoout

    INTEGER, INTENT(OUT) :: last

    
    INTEGER          :: i,         & ! Spatial index
                        it_Ntn,    & ! Iteration counter of the Newton method
                        max_Ntn      ! Allowed maximum of iterations
                        
    DOUBLE PRECISION :: eps,       & ! Precision on the density
                        rhoN,      & ! Density used by the Newton method
                        drhoN,     & ! Small variation of the density
                        A_q,         & ! Value of the EOS at rhoN
                        B_q,         & ! Value of the EOS at rhoN+drhoN
                        Der,       & ! Derivative of the EOS at rhoN
                        old_rhoN,  & ! Previous value of rhoN
                        oold_rhoN, & ! Previous (twice) value of rhoN
                        
                        ! Parameters needed for the initial guess of EOS H12
                        P_TF0,     & ! Pressure of a Fermi gas at ambiant conditions [Pa]
                        C10,       & ! First coefficient of H12 equation
                        C12,       & ! Second coefficient of H12 equation
                        cf_a,      & ! First coefficient of polynomial Taylor expansion
                        cf_b,      & ! Second coefficient of polynomial Taylor expansion
                        cf_d,      & ! Third coefficient of polynomial Taylor expansion
                        cf_p,      & ! First coefficient of reduced polynomial Taylor expansion
                        cf_q,      & ! Second coefficient of reduced polynomial Taylor expansion
                        delta        ! Discriminant of the reduced polynomial Taylor expansion

    LOGICAL          :: t_rho,     & ! Test on the convergence
                        t_cnt        ! Test on the number of iterations

    DOUBLE PRECISION, DIMENSION(2) :: rho_nn 

    DOUBLE PRECISION :: logrho_rock
                        
    !------


    rhoout = 0.d0
    last = 0

    !WRITE(6,*) 'In Newton_subrout'
    !WRITE(6,*) 'k =', k
    !WRITE(6,*) 'iEOS(k) =', iEOS(k)

    ! Initialization
    eps     = 1.d-4
    max_Ntn = 100

    ! You need to define the layer before calling Newton  
    drhoN = rho_0(k) * 1.d-4
       
   !!!!!!!!!!!!$OMP PARALLEL PRIVATE(it_Ntn, old_rhoN, oold_rhoN, t_rho, t_cnt, A_q, B_q, Der, rhoN)
    points: DO i=1, n_pts

    ! Cut


!             rhoN = 433.22d0
!
!             ! Checking the range of validity of the EOS
!             IF (rhoN < rho_lim_EOS(k)) THEN
!
!                 IF (i >= intrf(k+1)) THEN
!                    EXIT points
!                 ELSE IF (i < intrf(k)) THEN
!                    CYCLE points
!                 END IF
!             END IF
!
!
!          rhoout(i) = rhoN     
!          last = i

     ! Cut


          ! Initial guess (from first-order Taylor expansion)
          IF (EOSid==1 .OR. EOSid==2) THEN ! For BM3 and MGD EOSs
             rhoN = rho_0(k) * (Kp_0(k) - 5.d0 + sqrt(1.d0 + 2.d0*P(k,i)/K_0(k)*(Kp_0(k)-4.d0))) / (Kp_0(k)-4.d0)
             IF (isnan(rhoN)) rhoN = rho_0(k)
          ELSE IF (EOSid==3 .OR. EOSid==14) THEN ! For Vinet EOS
             rhoN = rho_0(k) * (Kp_0(k) - 2.d0 + sqrt(1.d0 + 2.d0*P(k,i)/K_0(k)*(Kp_0(k)-1.d0))) / (Kp_0(k)-1.d0)
             IF (isnan(rhoN)) rhoN = rho_0(k)
          ELSE IF (EOSid==7) THEN ! For H02 EOS
             rhoN = rho_0(k)
             rhoN = rho_0(k) * (Kp_0(k) - 4.d0 + sqrt(1.d0 + 2.d0*P(k,i)/K_0(k)*(Kp_0(k)-3.d0))) / (Kp_0(k)-3.d0)
             IF (isnan(rhoN)) rhoN = rho_0(k)
          ELSE IF (EOSid==8) THEN ! For H12 EOS
             P_TF0 = 2.337d-38 * (N_A*rho_0(k)*Z_eff(k)/Mmol(k))**(5.d0/3.d0)
             C10 = log(P_TF0/(3.d0*K_0(k)))
             C12 = 1.d0/2.d0 * (3.d0*Kp_0(k) - 2.d0*C10 - 9.d0)

             cf_a = - C12 / 9.d0
             cf_b = (C10+C12) / 3.d0
             cf_d = - P(k,i) / K_0(k)
             cf_p = - cf_b**2.d0/(3.d0*cf_a**2.d0) + 1.d0/cf_a
             cf_q = 2.d0*cf_b**3.d0/(27.d0*cf_a**3.d0) - cf_b/(3.d0*cf_a**2.d0) + cf_d/cf_a
             delta = cf_q**2.d0 + 4.d0*cf_p**3.d0/27.d0
             
             rhoN = rho_0(k) * (1.d0 - cf_b/(3.d0*cf_a) + (0.5d0*(-cf_q-sqrt(delta)))**(1.d0/3.d0) &
                  + (0.5d0*(-cf_q+sqrt(delta)))**(1.d0/3.d0))
             IF (isnan(rhoN)) rhoN = rho_0(k)
          ELSE IF(EOSid==6) THEN ! LJ potential  
	     rhoN=736.d0
          ELSE IF(EOSid==9) THEN ! IAPWS95  
             rhoN = 1000d0
             IF (isnan(rhoN)) rhoN = rho_0(k) 
          ELSE IF(EOSid==10) THEN ! Mazevet+19  
             rhoN = 1000d0
          ELSE IF(EOSid==11) THEN ! Haldemann+20 and Mazevet+19  
             rhoN = 1000d0
          ELSE IF(EOSid==12) THEN ! Chabrier EOS for H/He  
             rhoN = 200d0
          ELSE IF(EOSid==13) THEN ! Mazevet+22 EOS for H/He  
             rhoN = 200d0
          ELSE IF(EOSid==5) THEN ! Gibbs energy  
	     rhoN=rho_0(k)
             IF (isnan(rhoN)) rhoN = rho_0(k)
          ELSE IF (i==1) THEN ! No initial guess for other EOSs, using current value of rho
             rhoN = rho(k,i)
             WRITE(6,*) '/!\ WARNING in subroutine "new_rho": no initial guess for the used EOS.'
          END IF



          ! Initialization
          it_Ntn = 0
          
          old_rhoN  = 0.d0
          oold_rhoN = 0.d0

          t_rho = .TRUE.
          t_cnt = .TRUE.

          newtonloop: DO WHILE (t_rho .AND. t_cnt)




             it_Ntn = it_Ntn + 1

             ! Value and derivative of the EOS at current value of rho (or initial guess if possible)
             !A_q = EOS(rhoN,T(k,i),k,iEOS(k))
             !B_q = EOS(rhoN+drhoN,T(k,i),k,iEOS(k))

             A_q = EOS(rhoN,Tin(i),k,EOSid)
             B_q = EOS(rhoN+drhoN,Tin(i),k,EOSid)


             !IF (isnan(A_q) .OR. isnan(B_q)) THEN
                !WRITE(6,*) '/!\ ERROR in subroutine "new_rho": error on the value of new density.'
                !WRITE(6,*) '     k=', k, 'i=', i, 'rhoN=', rhoN
                !WRITE(6,*) A_q
                !WRITE(6,*) B_q
                !WRITE(6,*) T(k,i)
                !WRITE(6,*) P(k,i)
                !WRITE(6,*) iEOS(k)
                !STOP
             !END IF
             
             Der = (B_q-A_q) / drhoN
             
             ! Save previous values of rho
             oold_rhoN = old_rhoN
             old_rhoN  = rhoN

             ! New value of the density 
             !rhoN = rhoN - (A_q-P(k,i))/Der

             rhoN = rhoN - (A_q-Pin(i))/Der

             ! Checking the new value of rho
             IF (rhoN > 1.d100) THEN
                rhoN = old_rhoN
                !WRITE(6,*) '/!\ WARNING in subroutine "new_rho": new density too high, using less precise result.'
                EXIT newtonloop
             END IF

             !!!!!!!$OMP CRITICAL
             ! Checking the range of validity of the EOS
             IF (rhoN < rho_lim_EOS(k)) THEN

                 IF (i >= intrf(k+1)) THEN
                    EXIT points
                 ELSE IF (i < intrf(k)) THEN
                    CYCLE points
                 END IF
             END IF
             !!!!!!!$OMP END CRITICAL

 
             ! Test on possible oscillations of the Newton method
             IF (oold_rhoN == rhoN) THEN
                rhoN = (rhoN + old_rhoN) / 2.d0
             END IF

             ! Test of the convergence
             t_rho = abs(1.d0-old_rhoN/rhoN) .GT. eps
             t_cnt = it_Ntn .LT. max_Ntn



          END DO newtonloop

          rhoout(i) = rhoN     
          last = i



     END DO points
     !!!!!!!!!!!!$OMP END PARALLEL





    !------

  END SUBROUTINE Newton


  !--------------------------------------!
  !                                      !
  !        SUBROUTINE mod_Newton         !
  !                                      !
  !--------------------------------------!
  !
  ! Newton method to compute new density
  !
  SUBROUTINE mod_Newton(k, Pin, Tin, EOSid, last, rhoout)

    !------

    INTEGER, INTENT(IN) :: k, EOSid

    DOUBLE PRECISION, DIMENSION(441*121), INTENT(IN) :: Pin, Tin

    DOUBLE PRECISION, DIMENSION(441*121), INTENT(OUT) :: rhoout

    INTEGER, INTENT(OUT) :: last

    
    INTEGER          :: i,         & ! Spatial index
                        it_Ntn,    & ! Iteration counter of the Newton method
                        max_Ntn      ! Allowed maximum of iterations
                        
    DOUBLE PRECISION :: eps,       & ! Precision on the density
                        rhoN,      & ! Density used by the Newton method
                        drhoN,     & ! Small variation of the density
                        A_q,         & ! Value of the EOS at rhoN
                        B_q,         & ! Value of the EOS at rhoN+drhoN
                        Der,       & ! Derivative of the EOS at rhoN
                        old_rhoN,  & ! Previous value of rhoN
                        oold_rhoN, & ! Previous (twice) value of rhoN
                        
                        ! Parameters needed for the initial guess of EOS H12
                        P_TF0,     & ! Pressure of a Fermi gas at ambiant conditions [Pa]
                        C10,       & ! First coefficient of H12 equation
                        C12,       & ! Second coefficient of H12 equation
                        cf_a,      & ! First coefficient of polynomial Taylor expansion
                        cf_b,      & ! Second coefficient of polynomial Taylor expansion
                        cf_d,      & ! Third coefficient of polynomial Taylor expansion
                        cf_p,      & ! First coefficient of reduced polynomial Taylor expansion
                        cf_q,      & ! Second coefficient of reduced polynomial Taylor expansion
                        delta        ! Discriminant of the reduced polynomial Taylor expansion

    LOGICAL          :: t_rho,     & ! Test on the convergence
                        t_cnt        ! Test on the number of iterations

    DOUBLE PRECISION, DIMENSION(2) :: rho_nn 

    DOUBLE PRECISION :: logrho_rock
                        
    !------


    rhoout = 0.d0
    last = 0

    !WRITE(6,*) 'In Newton_subrout'
    !WRITE(6,*) 'k =', k
    !WRITE(6,*) 'iEOS(k) =', iEOS(k)

    ! Initialization
    eps     = 1.d-4
    max_Ntn = 100

    ! You need to define the layer before calling Newton  
    drhoN = rho_0(k) * 1.d-4
       
   !!!!!!!!!!!!$OMP PARALLEL PRIVATE(it_Ntn, old_rhoN, oold_rhoN, t_rho, t_cnt, A_q, B_q, Der, rhoN)
    points: DO i=1, 441*121



          ! Initial guess (from first-order Taylor expansion)
          IF (EOSid==1 .OR. EOSid==2) THEN ! For BM3 and MGD EOSs
             rhoN = rho_0(k) * (Kp_0(k) - 5.d0 + sqrt(1.d0 + 2.d0*P(k,i)/K_0(k)*(Kp_0(k)-4.d0))) / (Kp_0(k)-4.d0)
             IF (isnan(rhoN)) rhoN = rho_0(k)
          ELSE IF (EOSid==3 .OR. EOSid==14) THEN ! For Vinet EOS
             rhoN = rho_0(k) * (Kp_0(k) - 2.d0 + sqrt(1.d0 + 2.d0*P(k,i)/K_0(k)*(Kp_0(k)-1.d0))) / (Kp_0(k)-1.d0)
             IF (isnan(rhoN)) rhoN = rho_0(k)
          ELSE IF (EOSid==7) THEN ! For H02 EOS
             rhoN = rho_0(k)
             rhoN = rho_0(k) * (Kp_0(k) - 4.d0 + sqrt(1.d0 + 2.d0*P(k,i)/K_0(k)*(Kp_0(k)-3.d0))) / (Kp_0(k)-3.d0)
             IF (isnan(rhoN)) rhoN = rho_0(k)
          ELSE IF (EOSid==8) THEN ! For H12 EOS
             P_TF0 = 2.337d-38 * (N_A*rho_0(k)*Z_eff(k)/Mmol(k))**(5.d0/3.d0)
             C10 = log(P_TF0/(3.d0*K_0(k)))
             C12 = 1.d0/2.d0 * (3.d0*Kp_0(k) - 2.d0*C10 - 9.d0)

             cf_a = - C12 / 9.d0
             cf_b = (C10+C12) / 3.d0
             cf_d = - P(k,i) / K_0(k)
             cf_p = - cf_b**2.d0/(3.d0*cf_a**2.d0) + 1.d0/cf_a
             cf_q = 2.d0*cf_b**3.d0/(27.d0*cf_a**3.d0) - cf_b/(3.d0*cf_a**2.d0) + cf_d/cf_a
             delta = cf_q**2.d0 + 4.d0*cf_p**3.d0/27.d0
             
             rhoN = rho_0(k) * (1.d0 - cf_b/(3.d0*cf_a) + (0.5d0*(-cf_q-sqrt(delta)))**(1.d0/3.d0) &
                  + (0.5d0*(-cf_q+sqrt(delta)))**(1.d0/3.d0))
             IF (isnan(rhoN)) rhoN = rho_0(k)
          ELSE IF(EOSid==6) THEN ! LJ potential  
	     rhoN=736.d0
          ELSE IF(EOSid==9) THEN ! IAPWS95  
             rhoN = 1000d0
             IF (isnan(rhoN)) rhoN = rho_0(k) 
          ELSE IF(EOSid==10) THEN ! Mazevet+19  
             rhoN = 1000d0
          ELSE IF(EOSid==11) THEN ! Haldemann+20 and Mazevet+19  
             rhoN = 1000d0
          ELSE IF(EOSid==12) THEN ! Chabrier EOS for H/He  
             rhoN = 200d0
          ELSE IF(EOSid==13) THEN ! Mazevet+22 EOS for H/He  
             rhoN = 200d0
          ELSE IF(EOSid==5) THEN ! Gibbs energy  
	     rhoN=rho_0(k)
             IF (isnan(rhoN)) rhoN = rho_0(k)
          ELSE IF (i==1) THEN ! No initial guess for other EOSs, using current value of rho
             rhoN = rho(k,i)
             WRITE(6,*) '/!\ WARNING in subroutine "new_rho": no initial guess for the used EOS.'
          END IF



          ! Initialization
          it_Ntn = 0
          
          old_rhoN  = 0.d0
          oold_rhoN = 0.d0

          t_rho = .TRUE.
          t_cnt = .TRUE.

          newtonloop: DO WHILE (t_rho .AND. t_cnt)




             it_Ntn = it_Ntn + 1

             ! Value and derivative of the EOS at current value of rho (or initial guess if possible)
             !A_q = EOS(rhoN,T(k,i),k,iEOS(k))
             !B_q = EOS(rhoN+drhoN,T(k,i),k,iEOS(k))

             A_q = EOS(rhoN,Tin(i),k,EOSid)
             B_q = EOS(rhoN+drhoN,Tin(i),k,EOSid)


             !IF (isnan(A_q) .OR. isnan(B_q)) THEN
                !WRITE(6,*) '/!\ ERROR in subroutine "new_rho": error on the value of new density.'
                !WRITE(6,*) '     k=', k, 'i=', i, 'rhoN=', rhoN
                !WRITE(6,*) A_q
                !WRITE(6,*) B_q
                !WRITE(6,*) T(k,i)
                !WRITE(6,*) P(k,i)
                !WRITE(6,*) iEOS(k)
                !STOP
             !END IF
             
             Der = (B_q-A_q) / drhoN
             
             ! Save previous values of rho
             oold_rhoN = old_rhoN
             old_rhoN  = rhoN

             ! New value of the density 
             !rhoN = rhoN - (A_q-P(k,i))/Der

             rhoN = rhoN - (A_q-Pin(i))/Der

             ! Checking the new value of rho
             IF (rhoN > 1.d100) THEN
                rhoN = old_rhoN
                !WRITE(6,*) '/!\ WARNING in subroutine "new_rho": new density too high, using less precise result.'
                EXIT newtonloop
             END IF


 
             ! Test on possible oscillations of the Newton method
             IF (oold_rhoN == rhoN) THEN
                rhoN = (rhoN + old_rhoN) / 2.d0
             END IF

             ! Test of the convergence
             t_rho = abs(1.d0-old_rhoN/rhoN) .GT. eps
             t_cnt = it_Ntn .LT. max_Ntn



          END DO newtonloop

          rhoout(i) = rhoN     
          last = i



     END DO points
     !!!!!!!!!!!!$OMP END PARALLEL





    !------

  END SUBROUTINE mod_Newton



  
  !--------------------------------------!
  !                                      !
  !          SUBROUTINE new_rho          !
  !                                      !
  !--------------------------------------!
  !
  ! Computes the new density in each layer using a specific equation of state (EOS) via the Newton method.
  !
  SUBROUTINE new_rho

    !------
    
    INTEGER          :: i,         & ! Spatial index
                        k,         & ! Layer index
                        EOS_id,layer_id,last
                       
    DOUBLE PRECISION, DIMENSION(2) :: rho_nn 

    DOUBLE PRECISION, DIMENSION(n_pts) :: inverse_rho  

    ! Variables for AQUA, Haldemann et al. 2020
    ! REAL(8)  :: AQUA_RES(DIM_OUTPUT)


    INTEGER, DIMENSION(4) :: layer_array, & 
                             EOS_array,   &
                             last_array

    DOUBLE PRECISION,  DIMENSION(n_pts) :: P_mat, T_mat, rho_out        
                   
    DOUBLE PRECISION,  DIMENSION(4,n_pts) :: rho_mat   

    DOUBLE PRECISION :: rhoN,logrho_rock


     
    !------


    layer_array = (/1, 1, 2, 2/)
    !EOS_array = (/3, 10, 10, 12/)
    ! EOS_array = (/rock, water, water, HHe/)
    EOS_array = (/14, 10, 10, 12/)




    ! Calculation
    !$OMP PARALLEL DO PRIVATE(layer_id, EOS_id, P_mat, T_mat, rho_out, last)
    layers: DO k=1, 4
 
          layer_id = layer_array(k)
          EOS_id = EOS_array(k)
          
          P_mat = P(layer_id,:)
          T_mat = T(layer_id,:)


          !!!!! CALCULATION OF DENSITY !!!!!!!!!

         
          ! H/He
          IF (EOS_id==12) THEN

             points1: DO i=1, n_pts

                 rho_nn = EOS_chabrier(log10(P_mat(i)/1e9),log10(T_mat(i)))
                 rhoN = rho_nn(1)   


                 IF (rhoN < rho_lim_EOS(layer_id)) THEN
                 
                     IF (i >= intrf(layer_id+1)) THEN
                        EXIT points1
                     ELSE IF (i < intrf(layer_id)) THEN
                        CYCLE points1
                     END IF

                 END IF

                 rho_out(i) = rhoN     
                 last = i

             END DO points1

          ! rock
          ELSE IF (EOS_id==14) THEN

             points2: DO i=1, n_pts

                 logrho_rock = interp2d_opt(log10(T_mat(i)),log10(P_mat(i)*10),logT_sesame,logP_sesame, &
                             logrho_sesame,size(logT_sesame),size(logP_sesame))

                 IF (isnan(logrho_rock)) THEN
                     logrho_rock = -0.77286
                 END IF


                 rhoN = (10**logrho_rock)*1d3

                 IF (rhoN < rho_lim_EOS(layer_id)) THEN
                 
                     IF (i >= intrf(layer_id+1)) THEN
                        EXIT points2
                     ELSE IF (i < intrf(layer_id)) THEN
                        CYCLE points2
                     END IF

                 END IF

                 rho_out(i) = rhoN     
                 last = i

             END DO points2

          ! Mazevet for water
          ELSE IF (EOS_id==10) THEN

             points3: DO i=1, n_pts


                rhoN = interp2d_opt(P_mat(i)/1e9,T_mat(i),P_maz_water,T_maz_water, &
                   rho_maz_water,size(P_maz_water),size(T_maz_water))

                 IF (rhoN < rho_lim_EOS(layer_id)) THEN
                 
                     IF (i >= intrf(layer_id+1)) THEN
                        EXIT points3
                     ELSE IF (i < intrf(layer_id)) THEN
                        CYCLE points3
                     END IF

                 END IF

                 rho_out(i) = rhoN     
                 last = i

             END DO points3


          ! Other EOSs
          ELSE

             CALL Newton(layer_id, P_mat, T_mat, EOS_id, last, rho_out)
             

          END IF 

         
          !!!!! END OF CALCULATION OF DENSITY !!!!!!!!!

          
   
          rho_mat(k,:) = rho_out
          last_array(k) = last    

    END DO layers
    !$OMP END PARALLEL DO

    
    FORALL(i=last_array(1)+1:n_pts) rho_mat(1,i) = rho_surf_oliv
    FORALL(i=last_array(2)+1:n_pts) rho_mat(2,i) = rho_surf_w
    FORALL(i=last_array(3)+1:n_pts) rho_mat(3,i) = rho_surf_w
    FORALL(i=last_array(4)+1:n_pts) rho_mat(4,i) = rho_surf_HHe

    rho_oliv_core = rho_mat(1,:)
    rho_w_core = rho_mat(2,:)
    rho_w_env = rho_mat(3,:)
    rho_HHe = rho_mat(4,:)

    ! Mixture for core
    inverse_rho = 0.5d0/rho_oliv_core + 0.5d0/rho_w_core
    rho(1,:) = 1.d0/inverse_rho              


    ! Mixture for envelope
    ! Metals represented by water only
    inverse_rho = Zenv/rho_w_env + (1d0-Zenv)/rho_HHe

    ! Metals represented by water+olivine
    !inverse_rho = Zenv*0.5d0/rho_oliv_env + Zenv*0.5d0/rho_w_env + (1d0-Zenv)/rho_HHe

    rho(2,:) = max(1.d0/inverse_rho,rho_surf)  


    ! Save and update max values
    rho_max_old = rho_max
    rho_max = 0.d0

    !$OMP PARALLEL DO
    DO i=1, n_pts
       rho_max  = max(rho_max, rho(layer(i),i))
       rho_mmax = max(rho_mmax,rho(layer(i),i))
    END DO
    !$OMP END PARALLEL DO

    !------

  END SUBROUTINE new_rho

  !--------------------------------------!
  !                                      !
  !          SUBROUTINE rho_surface_w    !
  !                                      !
  !--------------------------------------!
  !
  ! Computes the surface density
  !
  SUBROUTINE rho_surface_w

    !------



    
    INTEGER          :: i,         & ! Spatial index
                        k,         & ! Layer index
                        it_Ntn,    & ! Iteration counter of the Newton method
                        max_Ntn,   & ! Allowed maximum of iterations
                        last2         ! Last position where the EOS is used
                        
    DOUBLE PRECISION :: eps,         & ! Precision on the density
                        rhoNew,      & ! Density used by the Newton method
                        drhoNew,     & ! Small variation of the density
                        A,           & ! Value of the EOS at rhoN
                        B,           & ! Value of the EOS at rhoN+drhoN
                        Deriva,      & ! Derivative of the EOS at rhoN
                        old_rhoNew,  & ! Previous value of rhoN
                        oold_rhoNew    ! Previous (twice) value of rhoN

                        

    LOGICAL          :: t_rho,     & ! Test on the convergence
                        t_cnt        ! Test on the number of iterations

     ! Variables for AQUA, Haldemann et al. 2020
     ! REAL(8)  :: AQUA_RES(DIM_OUTPUT)

                        
    !------

    ! Initialization
    eps     = 1.d-4
    max_Ntn = 100

    k = 2
    iEOS(k) = 10 ! Supercritical - Mazevet EOS
    !iEOS(k) = 11  ! For AQUA EOS




 !  IF (iEOS(k) == 11) THEN
 !
 !       CALL LOAD_TABLE_PT('Tables')
 !       CALL LOAD_TABLE_RhoT('Tables')
 !       CALL LOAD_TABLE_RhoU('Tables')
 !
 ! 
 !       AQUA_RES = INTERPOLATE_AQUA_PT(P_surf,T_surf)
 !       WRITE(6,*) 'AQUA_RES(1)'
 !       WRITE(6,*) AQUA_RES(1)
 !       rhoNew = AQUA_RES(1) 
 !       rho_surf = rhoNew
 ! 
 !   ELSE 
     




    ! Calculation       
       drhoNew = rho_0(k) * 1.d-4
       

          IF(iEOS(k)==6) THEN ! LJ potential  
	     rhoNew=736.d0
             IF (isnan(rhoNew)) rhoNew = 736.d0 
          ELSE IF(iEOS(k)==9) THEN ! IAPWS95  
             rhoNew = 1000d0
             IF (isnan(rhoNew)) rhoNew = rho_0(k) 
          ELSE IF(iEOS(k)==5) THEN ! Gibbs energy  
	     rhoNew=rho_0(k)
             IF (isnan(rhoNew)) rhoNew = rho_0(k)
          ELSE IF(iEOS(k)==10) THEN ! Mazevet+19  
             rhoNew = 1000d0
          ELSE IF(iEOS(k)==11) THEN ! Haldemann+20 and Mazevet+19  
             rhoNew = 1000d0
          !ELSE IF (i==1) THEN ! No initial guess for other EOSs, using current value of rho
             !rhoNew = rho(k,i)
             !WRITE(6,*) '/!\ WARNING in subroutine "new_rho": no initial guess for the used EOS.'
          END IF

          ! Initialization
          it_Ntn = 0
          
          old_rhoNew  = 0.d0
          oold_rhoNew = 0.d0

          t_rho = .TRUE.
          t_cnt = .TRUE.

          newton: DO WHILE (t_rho .AND. t_cnt)
             it_Ntn = it_Ntn + 1

             ! Value and derivative of the EOS at current value of rho (or initial guess if possible)
             A = EOS(rhoNew,T_surf,k,iEOS(k))
             B = EOS(rhoNew+drhoNew,T_surf,k,iEOS(k))

             IF (isnan(A) .OR. isnan(B)) THEN
                WRITE(6,*) '/!\ ERROR in subroutine "rho_surface": error on the value of new density.'
                WRITE(6,*) '     k=', k, 'i=', i, 'rhoN=', rhoNew
                WRITE(6,*) A
                WRITE(6,*) B
                WRITE(6,*) T_surf
                WRITE(6,*) P_surf
                !WRITE(6,*) iEOS(k)
                STOP
             END IF
             
             Deriva = (B-A) / drhoNew
             
             ! Save previous values of rho
             oold_rhoNew = old_rhoNew
             old_rhoNew  = rhoNew

             ! New value of the density 
             rhoNew = rhoNew - (A-P_surf)/Deriva
             

             ! Checking the new value of rho
             IF (rhoNew > 1.d100) THEN
                rhoNew = old_rhoNew
                WRITE(6,*) '/!\ WARNING in subroutine "new_rho": new density too high, using less precise result.'
                !EXIT newton
             END IF

 
             ! Test on possible oscillations of the Newton method
             IF (oold_rhoNew == rhoNew) THEN
                rhoNew = (rhoNew + old_rhoNew) / 2.d0
             END IF

             ! Test of the convergence
             t_rho = abs(1.d0-old_rhoNew/rhoNew) .GT. eps
             t_cnt = it_Ntn .LT. max_Ntn

          END DO newton

          !rho_surf_w = rhoNew/10.d0
          rho_surf_w = rhoNew

  !    END IF

    !------

  END SUBROUTINE rho_surface_w


  !--------------------------------------!
  !                                      !
  !          SUBROUTINE rho_surface_HHe  !
  !                                      !
  !--------------------------------------!
  !
  ! Computes the surface density
  !
  SUBROUTINE rho_surface_HHe

    !------



    
    INTEGER          :: i,         & ! Spatial index
                        k,         & ! Layer index
                        it_Ntn,    & ! Iteration counter of the Newton method
                        max_Ntn,   & ! Allowed maximum of iterations
                        last2         ! Last position where the EOS is used
                        
    DOUBLE PRECISION :: eps,         & ! Precision on the density
                        rhoNew,      & ! Density used by the Newton method
                        drhoNew,     & ! Small variation of the density
                        A,           & ! Value of the EOS at rhoN
                        B,           & ! Value of the EOS at rhoN+drhoN
                        Deriva,      & ! Derivative of the EOS at rhoN
                        old_rhoNew,  & ! Previous value of rhoN
                        oold_rhoNew    ! Previous (twice) value of rhoN

                        

    LOGICAL          :: t_rho,     & ! Test on the convergence
                        t_cnt        ! Test on the number of iterations


                        
    !------

    ! Initialization
    eps     = 1.d-4
    max_Ntn = 100

    k = 2
    iEOS(k) = 12 ! Chabrier EOS for H/He




    ! Calculation       
       drhoNew = rho_0(k) * 1.d-4
       

          IF(iEOS(k)==6) THEN ! LJ potential  
	     rhoNew=736.d0
             IF (isnan(rhoNew)) rhoNew = 736.d0 
          ELSE IF(iEOS(k)==9) THEN ! IAPWS95  
             rhoNew = 1000d0
             IF (isnan(rhoNew)) rhoNew = rho_0(k) 
          ELSE IF(iEOS(k)==5) THEN ! Gibbs energy  
	     rhoNew=rho_0(k)
             IF (isnan(rhoNew)) rhoNew = rho_0(k)
          ELSE IF(iEOS(k)==10) THEN ! Mazevet+19  
             rhoNew = 1000d0
          ELSE IF(iEOS(k)==11) THEN ! Haldemann+20 and Mazevet+19  
             rhoNew = 1000d0
          ELSE IF(iEOS(k)==12) THEN ! Chabrier EOS 
             rhoNew = 200d0
          !ELSE IF (i==1) THEN ! No initial guess for other EOSs, using current value of rho
             !rhoNew = rho(k,i)
             !WRITE(6,*) '/!\ WARNING in subroutine "new_rho": no initial guess for the used EOS.'
          END IF

          ! Initialization
          it_Ntn = 0
          
          old_rhoNew  = 0.d0
          oold_rhoNew = 0.d0

          t_rho = .TRUE.
          t_cnt = .TRUE.

          newton: DO WHILE (t_rho .AND. t_cnt)
             it_Ntn = it_Ntn + 1

             ! Value and derivative of the EOS at current value of rho (or initial guess if possible)
             A = EOS(rhoNew,T_surf,k,iEOS(k))
             B = EOS(rhoNew+drhoNew,T_surf,k,iEOS(k))

             IF (isnan(A) .OR. isnan(B)) THEN
                WRITE(6,*) '/!\ ERROR in subroutine "rho_surface": error on the value of new density.'
                WRITE(6,*) '     k=', k, 'i=', i, 'rhoN=', rhoNew
                WRITE(6,*) A
                WRITE(6,*) B
                WRITE(6,*) T_surf
                WRITE(6,*) P_surf
                !WRITE(6,*) iEOS(k)
                STOP
             END IF
             
             Deriva = (B-A) / drhoNew
             
             ! Save previous values of rho
             oold_rhoNew = old_rhoNew
             old_rhoNew  = rhoNew

             ! New value of the density 
             rhoNew = rhoNew - (A-P_surf)/Deriva
             

             ! Checking the new value of rho
             IF (rhoNew > 1.d100) THEN
                rhoNew = old_rhoNew
                WRITE(6,*) '/!\ WARNING in subroutine "new_rho": new density too high, using less precise result.'
                !EXIT newton
             END IF

 
             ! Test on possible oscillations of the Newton method
             IF (oold_rhoNew == rhoNew) THEN
                rhoNew = (rhoNew + old_rhoNew) / 2.d0
             END IF

             ! Test of the convergence
             t_rho = abs(1.d0-old_rhoNew/rhoNew) .GT. eps
             t_cnt = it_Ntn .LT. max_Ntn

          END DO newton

          rho_surf_HHe = rhoNew


    !------

  END SUBROUTINE rho_surface_HHe
  

  


  !-----------------------------------------!
  !                                         !
  !          SUBROUTINE interfaces          !
  !                                         !
  !-----------------------------------------!
  !
  ! Computes the new interface positions between the layers.
  !
  SUBROUTINE interfaces

    !------

    INTEGER                            :: i,  & ! Spatial index
                                          k,lll,i_o,i_oo,i_int     ! Layer index

    DOUBLE PRECISION                   :: M_lay ! Mass of a layer

    DOUBLE PRECISION, DIMENSION(n_pts)                   :: arr    
    !------

    ! Saving previous positions
    intrf_old7   = intrf_old6
    intrf_old6   = intrf_old5
    intrf_old5   = intrf_oooold
    intrf_oooold = intrf_ooold
    intrf_ooold  = intrf_oold
    intrf_oold   = intrf_old
    intrf_old    = intrf

    ! New positions
    DO k=1, n_lay
       ! If layer k is the core
       IF (ilayer(k)==1) THEN
          intrf(k+1) = intrf(k)
           
          DO WHILE (mass_btw(intrf(k),intrf(k+1),k) < M_P*x_core)
             intrf(k+1) = intrf(k+1) + 1

             IF (intrf(k+1) > n_pts) THEN
                WRITE(6,*) '/!\ ERROR in subroutine "interfaces": the new interface radius is out of range. Core'
                STOP
             END IF
          END DO
          


       ! If layer k is the supercritical water layer
       ELSE IF (ilayer(k)==2) THEN
          intrf(k+1) = intrf(k)
  
          !M_lay = M_P*x_H2O 

          !WRITE(6,*) 'M_lay before calculation'
          !WRITE(6,*) M_lay
          !WRITE(6,*) 'M_lay 1st calculation'


          !M_lay = M_P - mass_btw(1,intrf(k),0)
          !WRITE(6,*) M_lay


          IF (x_H2O.eq.0.d0.or.M_lay.le.0.d0) CYCLE
  
          
          IF (.true.) THEN  
          lll = n_pts
  
          DO WHILE (rho(k,lll).le.rho_surf)  ! This loop finds the point at which rho=rho_surf
             lll = lll-1
          END DO
  
          intrf(k+1) = lll+1  ! Planet radius goes until rho=rho_surf
  
          DO WHILE(intrf(k+1).lt.intrf(k)) ! This moves the profile outwards if Rp < interface mantle-supercrit
             arr(:) = rho(k,:)
             !$OMP PARALLEL DO
             DO i=2,n_pts
                rho(k,i) = arr(i-1)
             END DO
             !$OMP END PARALLEL DO
             rho(k,1) = arr(1)
             intrf(k+1) = min(intrf(k+1)+1,n_pts)
          END DO


  
          M_lay = M_P - mass_btw(1,intrf(k),0) ! Update of the mass of the supercrit layer with new intrf 
          
          
          !WRITE(6,*) 'M_lay 2nd calculation'
          !WRITE(6,*) mass_btw(1,intrf(k),0)
          !WRITE(6,*) intrf(k)
          !WRITE(6,*) intrf
          !WRITE(6,*) M_lay

          !WRITE(6,*) mass_btw(1,intrf(k),1)
          !WRITE(6,*) r(1:2)
          

         IF (isnan(M_lay)) THEN
             WRITE(6,*) 'Error in interior structure model (Fortran): Mass of core layer is NaN'
             WRITE(6,*) 'This is likely due to a spatial resolution of the radius grid that is'
             WRITE(6,*) 'too small to resolve the size of the core'
             WRITE(6,*) 'Increase the resolution of the grid by setting a lower value for the'
             WRITE(6,*) 'input parameter pow_law_formass'             
             STOP 
         END IF



         lll = 0
         i_o = -1
         i_oo = -2
         i_int = 1  ! These are counters
  
         DO WHILE(lll.lt.n_pts) ! This is to avoid intrf(k+1) > maximum number points
            i_oo = i_o
            i_o = intrf(k+1)
            lll = lll+1
           
            !WRITE(6,*) 'interfaces'
            !WRITE(6,*) intrf
            !WRITE(6,*) mass_btw(intrf(k),intrf(k+1),k)
            !WRITE(6,*) M_lay

            

            IF (mass_btw(intrf(k),intrf(k+1),k) < M_lay) THEN ! Move profile outwards until mass of supercritical
               arr(:) = rho(k,:)                              ! layer is filled
               !$OMP PARALLEL DO
               DO i=2,n_pts
                  rho(k,i) = arr(i-1)
               END DO
               !$OMP END PARALLEL DO
               rho(k,1) = arr(1)
               intrf(k+1) = min(intrf(k+1)+1,n_pts)
            ELSE           ! This moves profile inwards if the mass of the supercrit layer is too much 
               arr(:) = rho(k,:) 
               !$OMP PARALLEL DO
               DO i=1,n_pts-1
                  rho(k,i) = arr(i+1)
               END DO
               !$OMP END PARALLEL DO
               rho(k,n_pts) = arr(n_pts)
               intrf(k+1) = max(intrf(k+1)-1,1)
               !intrf(k+1) = max(intrf(k+1)-1,intrf(k))
            END IF
            
            IF (intrf(k+1).eq.i_oo.or.intrf(k+1).eq.i_o) EXIT  ! This detects the bouncing between two i positions 
          END DO                                               ! of the supercritical - void interface
             
          CYCLE    ! This jumps to next layer (void)
          END IF

!          !##########################################

!              DO WHILE ((mass_btw(intrf(k),intrf(k+1),k) < M_lay) .AND. (intrf(k+1) < n_pts))
!                  intrf(k+1) = intrf(k+1) + 1
!
!                  IF (intrf(k+1) > n_pts) THEN
!                     WRITE(6,*) '/!\ ERROR in subroutine "interfaces": the new interface radius is out of range. Supercritical'
!                     STOP
!                  END IF
!              END DO
          


       ! If layer k is the void
       ELSE IF (ilayer(k)==3) THEN          
               CYCLE ! Do nothing

       ! If layer k is not defined
       ELSE
          WRITE(6,*) '/!\ ERROR in subroutine "interfaces": this type of layer is not defined.'

          
       END IF
    END DO
    
    
    !------
        
  END SUBROUTINE interfaces



  !------------------------------------------------!
  !                                                !
  !          SUBROUTINE check_convergence          !
  !                                                !
  !------------------------------------------------!
  !
  ! Computes the mass of each layer in the planet.
  !
  SUBROUTINE check_convergence
  
    INTEGER :: prec_intrf

    !------

    
    !prec_intrf = 0

    !IF ((ANY(intrf/=intrf_old) .AND. ALL(abs(intrf-intrf_oold)<=prec_intrf)) .OR. &
         !(ANY(intrf/=intrf_old) .AND. ANY(intrf/=intrf_oold) .AND. ALL(abs(intrf-intrf_ooold)<=prec_intrf)) .OR. &
         !(ANY(intrf/=intrf_old) .AND. ANY(intrf/=intrf_oold) .AND. ANY(intrf/=intrf_ooold) .AND. &
         !ALL(abs(intrf-intrf_oooold)<=prec_intrf)) .OR. (ANY(intrf/=intrf_old) .AND. ANY(intrf/=intrf_oold) .AND. &
         !ANY(intrf/=intrf_ooold) .AND. ANY(intrf/=intrf_oooold) .AND. ALL(abs(intrf-intrf_old5)<=prec_intrf)) .OR. &
         !(ANY(intrf/=intrf_old) .AND. ANY(intrf/=intrf_oold) .AND. ANY(intrf/=intrf_ooold) .AND. &
         !ANY(intrf/=intrf_oooold) .AND. ANY(intrf/=intrf_old5) .AND. ALL(abs(intrf-intrf_old6)<=prec_intrf)) .OR. &
         !(ANY(intrf/=intrf_old) .AND. ANY(intrf/=intrf_oold) .AND. ANY(intrf/=intrf_ooold) .AND. &
         !ANY(intrf/=intrf_oooold) .AND. ANY(intrf/=intrf_old5) .AND. ANY(intrf/=intrf_old6) .AND. &
         !ALL(abs(intrf-intrf_old7)<=prec_intrf))) THEN


    IF ((ANY(intrf/=intrf_old) .AND. ALL(intrf==intrf_oold)) .OR. &
         (ANY(intrf/=intrf_old) .AND. ANY(intrf/=intrf_oold) .AND. ALL(intrf==intrf_ooold)) .OR. &
         (ANY(intrf/=intrf_old) .AND. ANY(intrf/=intrf_oold) .AND. ANY(intrf/=intrf_ooold) .AND. &
         ALL(intrf==intrf_oooold))) THEN
       cnt_conv = cnt_conv + 1
    END IF


       !cnt_conv = cnt_conv + 1

       !WRITE(6,*) cnt_conv       

    !END IF
    
    !------

  END SUBROUTINE check_convergence


  !------------------------------------------!
  !                                          !
  !          SUBROUTINE mass_layers          !
  !                                          !
  !------------------------------------------!
  !
  ! Computes the mass of each layer in the planet.
  !
  SUBROUTINE mass_layers

    !------

    INTEGER :: k ! Layer index
    
    !------

    !$OMP PARALLEL DO
    DO k=1, n_lay
       !WRITE(6,*) 'mass_layers'
       !WRITE(6,*) intrf
       M(k) = mass_btw(intrf(k),intrf(k+1),k)
    END DO
    !$OMP END PARALLEL DO

    IF (M(n_lay)<0.d0) M(n_lay) = 0.d0

    !------
    
  END SUBROUTINE mass_layers



  !-----------------------------------------!
  !                                         !
  !          SUBROUTINE FeSi_ratio          !
  !                                         !
  !-----------------------------------------!
  !
  ! Computes the Fe/Si mole ratio in the planet obtained with the program.
  !
  SUBROUTINE FeSi_ratio

    !------

    INTEGER                            :: k       ! Layer index

    DOUBLE PRECISION                   :: Fe,   & ! Total amount (mole) of Fe in the planet
                                          Si      ! Total amount (mole) of Si in the planet

    DOUBLE PRECISION, DIMENSION(n_lay) :: x_Fe, & ! Amount (mole) of Fe in the layers
                                          x_Si    ! Amount (mole) of Si in the layers
    
    !------

    Fe = 0.d0
    Si = 0.d0
    
    ! For adaptations: x(number layer, number material in that layer)

    x_Fe(1) = 0.d0*x(1,1) + 2.d0*x(1,2) 
    x_Fe(2) = 0.d0*x(2,1)


    x_Si(1) = 1.d0*x(1,1) + 1.d0*x(1,2) 
    x_Si(2) = 0.d0*x(2,1)


    DO k=1, n_lay-1
       Fe = Fe + x_Fe(k)*M(k)/Mmol(k)
       Si = Si + x_Si(k)*M(k)/Mmol(k)
    END DO

    FeSic = Fe/Si

    !------
    
  END SUBROUTINE FeSi_ratio



  !------------------------------------------!
  !                                          !
  !          SUBROUTINE mom_inertia          !
  !                                          !
  !------------------------------------------!
  !
  ! Computes the polar moment of inertia of the planet (and its portion due to the solid outer shell)
  ! obtained with the program.
  !
  SUBROUTINE mom_inertia

    !------

    INTEGER                            :: k, & ! Layer index
                                          i    ! Radius index

    DOUBLE PRECISION, DIMENSION(n_pts) :: arr  ! Temporary array
    
    !------

    C_iner_tot = 0.d0

    DO k=1, n_lay
       DO i=1, n_pts
          arr(i) = rho(k,i)*r(i)**4.d0
       END DO
       
       C_iner(k) = 8.d0/3.d0*pi * integral(intrf(k),intrf(k+1),arr,r)

       C_iner_tot = C_iner_tot + C_iner(k)
    END DO

    ! Computing C/MR2 and Cm/C

    !C_MoI1 = C_iner_tot / (M_Pc * R_Pc**2.d0)
    !C_MoI2 = (C_iner(2)+C_iner(3)) / C_iner_tot

    !------
    
  END SUBROUTINE mom_inertia


END MODULE subroutines
