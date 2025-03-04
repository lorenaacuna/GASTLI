!      /''''''''''''''''''''''''''''\
!     /          PARAMETERS          \
!     \............................../

! This module lists:
!      (1) the variables of the program
!      (2) the input parameters, as well as their values (in the subroutines)
! in order to be able to use them everywhere in the code.

  
MODULE fort_input

  IMPLICIT NONE

  
CONTAINS


  !----------------------------------------------!
  !                                              !
  !          SUBROUTINE read_layer_parameters    !
  !                                              !
  !----------------------------------------------!
  !
  ! Initializes all input parameters by reading the values in specific files.
  !
  SUBROUTINE read_layer_parameters(path_to_file,n_lay,n_mat_max, ilayer,use_lay,iEOS,EOS_th,del_T,n_mat,x,&
             cf_Mmol,cf_Z_eff, cf_rho_0,cf_T_0,cf_K_0,cf_Kp_0,cf_gam_0,cf_q,cf_a_T,cf_b_T,cf_c_T,cf_a_P,cf_T_D0)

    !------

    INTEGER, INTENT(IN)  :: n_lay, n_mat_max
    
    CHARACTER(len=*), INTENT(IN) :: path_to_file

    CHARACTER(LEN=40)                      :: name        ! Name of the "layer*.in" files

    INTEGER                                :: k,        & ! Layer index
                                              unit,     & ! Unit of the "layer*.in" files
                                              h,        & ! Lines read index
                                              l           ! Coefficients index

    INTEGER,          INTENT(OUT)           ::      ilayer(n_lay),              & ! List of layers
                                                    use_lay(n_lay),             & ! Presence of each layer in the planet
                                                    iEOS(n_lay),                & ! EOS to use in each layer
                                                    EOS_th(n_lay),              & ! Consideration of the thermal component
                                                    n_mat(n_lay)                  ! Number of different materials per layer


    DOUBLE PRECISION,   INTENT(OUT)           :: del_T(n_lay)                  ! Upper temperature drop of the layers [K]

    DOUBLE PRECISION,   INTENT(OUT) :: x(n_lay,n_mat_max)                      ! Mole fraction of layer's materials

    DOUBLE PRECISION,  INTENT(OUT)  :: cf_Mmol(n_lay,n_mat_max),  & ! Coefficients for molar mass
                                                   cf_Z_eff(n_lay,n_mat_max), & ! Coefficients for effective atomic charge
                                                   cf_rho_0(n_lay,n_mat_max), & ! Coefficients for reference density
                                                   cf_T_0(n_lay,n_mat_max),   & ! Coefficients for reference temperature
                                                   cf_K_0(n_lay,n_mat_max),   & ! Coefficients for reference bulk modulus
                                                   cf_Kp_0(n_lay,n_mat_max),  & ! Coefficients for pressure derivative of bulk mod.
                                                   cf_gam_0(n_lay,n_mat_max), & ! Coefficients for reference Gruneisen parameter
                                                   cf_q(n_lay,n_mat_max),     & ! Coefficients for adiabatic power exponent
                                                   cf_a_T(n_lay,n_mat_max),   & ! Coefficients for 1st alpha parameter
                                                   cf_b_T(n_lay,n_mat_max),   & ! Coefficients for 2nd alpha parameter
                                                   cf_c_T(n_lay,n_mat_max),   & ! Coefficients for 3rd alpha parameter
                                                   cf_a_P(n_lay,n_mat_max),   & ! Coefficients for temp. derivative of bulk modulus
                                                   cf_T_D0(n_lay,n_mat_max)     ! Coefficients for reference Debye temperature

    CHARACTER(LEN=12)     :: lay_name(n_lay)                ! Name of the layers
    CHARACTER(LEN=7)      :: lay_col(n_lay)                ! Color of the layers, to use in plot
    CHARACTER(LEN=12) :: mat_name(n_lay,n_mat_max)              ! Name of the materials in a layer
    CHARACTER(LEN=12) :: mol_name(n_lay,n_mat_max)               ! Name of the molecules in a layer



    !------


       cf_Mmol  = 0.d0
       cf_Z_eff = 0.d0
       cf_rho_0 = 0.d0
       cf_T_0   = 0.d0
       cf_K_0   = 0.d0
       cf_Kp_0  = 0.d0
       cf_gam_0 = 0.d0
       cf_q     = 0.d0
       cf_a_T   = 0.d0
       cf_b_T   = 0.d0
       cf_c_T   = 0.d0
       cf_a_P   = 0.d0
       cf_T_D0  = 0.d0

    DO k=1, n_lay

       x(k,:)   = 0.d0

       ! Reinitalize characters and parameters
       DO l=1, n_mat_max
          mat_name(k,l) = '------------'
          mol_name(k,l) = '------------'
       END DO


       ! Opening file
       WRITE(name,'(A11,I1.1,A3)') 'Input/layer', k, '.in'
       unit = 30 + k

       OPEN(unit, FILE=TRIM(ADJUSTL(path_to_file))//name, STATUS='old', ACTION='read')

       ! WRITE(6,*) name

       ! Skip header
       DO h=1, 7
          READ(unit,*)
       END DO

       READ(unit,*) ilayer(k)

       IF(ilayer(k)/=k) THEN
          WRITE(6,'(A72,I1,A1)') ' /!\ ERROR in subroutine "read_parameters": wrong index value for layer ', k, '.'
          STOP
       END IF

       READ(unit,'(A12)') lay_name(k)
       READ(unit,'(A7)')  lay_col(k)
       READ(unit,*)       use_lay(k)

       IF(use_lay(k)/=1) THEN
          WRITE(6,*) '/!\ ERROR in subroutine "read_parameters": the not-use of a layer is not implemented yet.'
          STOP
       END IF

       READ(unit,*) iEOS(k)
       READ(unit,*) EOS_th(k)
       READ(unit,*) del_T(k)
       READ(unit,*) n_mat(k)

       ! Skip to section "Composition"
       DO h=1, 4
          READ(unit,*)
       END DO

       READ(unit,*) (mat_name(k,l), l=1,n_mat(k))
       READ(unit,*) (mol_name(k,l), l=1,n_mat(k))
       READ(unit,*) (x(k,l),        l=1,n_mat(k))
       READ(unit,*)
       READ(unit,*) (cf_Mmol(k,l),    l=1,n_mat(k))
       READ(unit,*) (cf_Z_eff(k,l),   l=1,n_mat(k))
       READ(unit,*) (cf_rho_0(k,l),   l=1,n_mat(k))
       READ(unit,*)
       READ(unit,*) (cf_T_0(k,l),     l=1,n_mat(k))
       READ(unit,*) (cf_K_0(k,l),     l=1,n_mat(k))
       READ(unit,*) (cf_Kp_0(k,l),    l=1,n_mat(k))
       READ(unit,*) (cf_gam_0(k,l),   l=1,n_mat(k))
       READ(unit,*) (cf_q(k,l),       l=1,n_mat(k))
       READ(unit,*)
       READ(unit,*) (cf_a_T(k,l),     l=1,n_mat(k))
       READ(unit,*) (cf_b_T(k,l),     l=1,n_mat(k))
       READ(unit,*) (cf_c_T(k,l),     l=1,n_mat(k))
       READ(unit,*) (cf_a_P(k,l),     l=1,n_mat(k))
       READ(unit,*)
       READ(unit,*) (cf_T_D0(k,l),    l=1,n_mat(k))


       ! Closing file
       CLOSE(unit)


    END DO

    !------
        
  END SUBROUTINE read_layer_parameters


  !---------------------------------------------!
  !                                             !
  !          SUBROUTINE read_constants          !
  !                                             !
  !---------------------------------------------!
  !
  ! Reads the value of some constants in a specific file.
  !
  SUBROUTINE read_constants(path_to_file,Ttp,Ptp)

    !------

    CHARACTER(len=*), INTENT(IN) :: path_to_file


    ! Thermodynamical constants (read)

    DOUBLE PRECISION, DIMENSION(6), INTENT(OUT) :: Ttp,                          & ! Temperature of water's triple points [K]
                                                   Ptp                             ! Pressure of water's triple points [Pa]


    INTEGER :: h, & ! Lines read index
               i    ! Triple points index
    
    !------

    ! FILE "constants.in"
    ! -------------------
    
    OPEN(29, FILE=TRIM(ADJUSTL(path_to_file))//'Input/constants.in', STATUS='old', ACTION='read')

    ! Skip header to section "Water triple points"
    DO h=1, 10
       READ(29,*)
    END DO

    DO i=1, 6
       READ(29,*) Ttp(i), Ptp(i)
    END DO

    CLOSE(29)

    !------
        
  END SUBROUTINE read_constants







  !----------------------------------------------!
  !                                              !
  !         SUBROUTINE read_EOS_sesame           !
  !                                              !
  !----------------------------------------------!
  !
  ! This reads SESAME's EOS for rock
  !
  SUBROUTINE read_EOS_sesame(path_to_file, logT_sesame,logP_sesame,logrho_sesame,logS_sesame, &
                                                         dlrho_dlT_p_sesame,    &
                                                         dlS_dlT_p_sesame,      &
                                                         dlrho_dlP_t_sesame )

    !------

    CHARACTER(len=*), INTENT(IN) :: path_to_file


    ! SESAME EOS for rock (dry sand)

    DOUBLE PRECISION, DIMENSION(64), INTENT(OUT) :: logT_sesame

    DOUBLE PRECISION, DIMENSION(41), INTENT(OUT) :: logP_sesame

    DOUBLE PRECISION, DIMENSION(41*64), INTENT(OUT) ::   logrho_sesame,      &
                                                         logS_sesame,        &
                                                         dlrho_dlT_p_sesame, &
                                                         dlS_dlT_p_sesame,   &
                                                         dlrho_dlP_t_sesame 

    DOUBLE PRECISION, DIMENSION(7,64*41) :: table_sesame
    

    OPEN(UNIT = 158,file=TRIM(ADJUSTL(path_to_file))//"Input/SESAME/logPcgs_sesame.dat", status='OLD', action = 'READ') 
    READ(158,*) logP_sesame
    CLOSE(UNIT = 158)


    OPEN(UNIT = 160,file=TRIM(ADJUSTL(path_to_file))//"Input/SESAME/sesame_sand_eos.dat", status='OLD', action = 'READ') 
    READ(160,*) table_sesame
    CLOSE(UNIT = 160)

    logrho_sesame = table_sesame(3,1:)
    logS_sesame = table_sesame(4,1:)
    dlrho_dlT_p_sesame = table_sesame(5,1:)
    dlS_dlT_p_sesame = table_sesame(6,1:)
    dlrho_dlP_t_sesame = table_sesame(7,1:)

    logT_sesame = (/1.89,1.94,1.99,2.04,2.09,2.14,2.19,2.24,2.29, &
    2.34,2.39,2.44,2.49,2.54,2.59,2.64,2.69,2.74,2.79,2.84,2.89,2.94,2.99,3.04, &
    3.09,3.14,3.19,3.24,3.29,3.34,3.39,3.44,3.49,3.54,3.59,3.64,3.69,3.74,3.79, &
    3.84,3.89,3.94,3.99,4.04,4.09,4.14,4.19,4.24,4.29,4.34,4.39,4.44,4.49,4.54, &
    4.59,4.64,4.69,4.74,4.79,4.84,4.89,4.94,4.99,5.04/)



    !------


  END SUBROUTINE read_EOS_sesame






  !----------------------------------------------!
  !                                              !
  !       SUBROUTINE read_mazevet_water          !
  !                                              !
  !----------------------------------------------!
  !
  ! This reads Mazevet EOS for water (file only for density)
  !
  SUBROUTINE read_mazevet_water(path_to_file,P_maz_water,T_maz_water,rho_maz_water)  

    !------

    CHARACTER(len=*), INTENT(IN) :: path_to_file


    DOUBLE PRECISION, DIMENSION(441), INTENT(OUT) :: P_maz_water

    DOUBLE PRECISION, DIMENSION(121), INTENT(OUT) :: T_maz_water

    DOUBLE PRECISION, DIMENSION(441*121), INTENT(OUT) :: rho_maz_water


    DOUBLE PRECISION, DIMENSION(3,441*121) :: table_mazevet_water

    DOUBLE PRECISION :: logT_maz, logP_maz

    INTEGER :: Nt, Np, l


    !------

  Nt = 121
  Np = 441

 
  OPEN(UNIT = 158,file=TRIM(ADJUSTL(path_to_file))//"Input/Mazevet_water/EOS_water_Mazevet.dat", status='OLD', action = 'READ') 

  READ(158,*) table_mazevet_water

  CLOSE(UNIT = 158)

  rho_maz_water = table_mazevet_water(3,1:)


  DO l = 1, Nt
     logT_maz = 2.0+(l-1)*0.05 
     T_maz_water(l) = 10**logT_maz    
  END DO 


  DO l = 1, Np
     logP_maz = -9.0+(l-1)*0.05    
     P_maz_water(l) = 10**logP_maz         
  END DO 



    !------

  END SUBROUTINE read_mazevet_water  





  !----------------------------------------------!
  !                                              !
  !          SUBROUTINE read_EOS_ch              !
  !                                              !
  !----------------------------------------------!
  !
  ! This reads Chabrier's EOS
  !
  SUBROUTINE read_EOS_ch(path_to_file, logT_input,logP_input,logrho_input, &
                         logrho_ch, logU_ch, logP_ch,        &
                         logS_ch, dlrho_dlT_P,               &
                         dlrho_dlP_T, dlS_dlT_P, grad_ad, grad_ad_PT)

    !------

    CHARACTER(len=*), INTENT(IN) :: path_to_file


    DOUBLE PRECISION, DIMENSION(121*441), INTENT(OUT) ::  logrho_ch, grad_ad_PT

    
    DOUBLE PRECISION, DIMENSION(121*241), INTENT(OUT) ::  logU_ch, logP_ch, &
                                                    logS_ch, dlrho_dlT_P,        &
                                                    dlrho_dlP_T, dlS_dlT_P, grad_ad
                                               

                                                  
    DOUBLE PRECISION, DIMENSION(121), INTENT(OUT) ::  logT_input 
                                           
    DOUBLE PRECISION, DIMENSION(441), INTENT(OUT) ::  logP_input 

    DOUBLE PRECISION, DIMENSION(241), INTENT(OUT) ::  logrho_input 


    INTEGER                                :: i,j,      & ! 
                                              l,        & ! 
                                              unit,     & ! 
                                              Nt, Np, Ntot, c, Nrho, NtotU


    DOUBLE PRECISION, DIMENSION(10)   :: line   

    !------



  Nt = 121
  Np = 441
  Nrho = 241

  Ntot = Nt*Np
  NtotU = Nt*Nrho

  ! Table 1

  !ALLOCATE(logrho_ch(Ntot))
  !ALLOCATE(grad_ad_PT(Ntot))



  ! Table 2

  !ALLOCATE(logP_ch(NtotU))             ! Column 2
  !ALLOCATE(logU_ch(NtotU))             ! Column 4
  !ALLOCATE(logS_ch(NtotU))             ! Column 5
  !ALLOCATE(dlrho_dlT_P(NtotU))         ! Column 6
  !ALLOCATE(dlrho_dlP_T(NtotU))         ! Column 7
  !ALLOCATE(dlS_dlT_P(NtotU))           ! Column 8
  !ALLOCATE(grad_ad(NtotU))             ! Column 10


  ! For rho(P,T) table - DENSITY 

  unit = 65 


  OPEN(unit, FILE=TRIM(ADJUSTL(path_to_file))//'Input/Chabrier/TABLEEOS_2021_TP_Y0275_v1', STATUS='old', ACTION='read')

  c = 1

  ! Skip header
  READ(unit,*)

  DO j = 1, Nt

     ! Skip # iT... line
     READ(unit,*)

     ! Read isotherm
     DO i = 1, Np

          READ(unit,*) (line(l),     l=1,10)

          !logT_ch(c) = line(1)
          !logP_ch(c) = line(2)
          logrho_ch(c) = line(3)
          grad_ad_PT(c) = line(10)

          c = c+1
                   
     END DO

  END DO

  CLOSE(unit)



  ! For U(rho,T) table - INTERNAL ENERGY

  unit = 65 

  OPEN(unit, FILE=TRIM(ADJUSTL(path_to_file))//'Input/Chabrier/TABLEEOS_2021_Trho_Y0275_v1', STATUS='old', ACTION='read')

  c = 1

  ! Skip header
  READ(unit,*)

  DO j = 1, Nt

     ! Skip # iT... line
     READ(unit,*)

     ! Read isotherm
     DO i = 1, Nrho

          READ(unit,*) (line(l),     l=1,10)

          logP_ch(c) = line(2)
          logU_ch(c) = line(4)
          logS_ch(c) = line(5)
          dlrho_dlT_P(c) = line(6)
          dlrho_dlP_T(c) = line(7)
          dlS_dlT_P(c) = line(8)
          grad_ad(c) = line(10)

          c = c+1
                   
     END DO

  END DO

  CLOSE(unit)



  DO l = 1, Nt
     logT_input(l) = 2.0+(l-1)*0.05      
  END DO 


  DO l = 1, Np
     logP_input(l) = -9.0+(l-1)*0.05      
  END DO 

  DO l = 1, Nrho
     logrho_input(l) = -6.0+(l-1)*0.05      
  END DO 


    !------



  END SUBROUTINE read_EOS_ch



  !----------------------------------------------!
  !                                              !
  !          SUBROUTINE read_HG23_corr           !
  !                                              !
  !----------------------------------------------!
  !
  ! This reads HG23 corrections for EOS of H and He
  !
  SUBROUTINE read_HG23_corr(path_to_file,logP_HG,logT_HG,Vmix_HG,Smix_HG)

    !------

    CHARACTER(len=*), INTENT(IN) :: path_to_file


    ! HG23 EOS correction for H/He mixture

    DOUBLE PRECISION, DIMENSION(441), INTENT(OUT) :: logP_HG

    DOUBLE PRECISION, DIMENSION(121), INTENT(OUT) :: logT_HG

    DOUBLE PRECISION, DIMENSION(441*121), INTENT(OUT) :: Vmix_HG, Smix_HG


    DOUBLE PRECISION, DIMENSION(4,441*121) :: table_HG

    !------

 
  OPEN(UNIT = 158,file=TRIM(ADJUSTL(path_to_file))//"Input/Chabrier/logP_HG23_corr.dat", status='OLD', action = 'READ') 
  READ(158,*) logP_HG 
  CLOSE(UNIT = 158)

  OPEN(UNIT = 158,file=TRIM(ADJUSTL(path_to_file))//"Input/Chabrier/logT_HG23_corr.dat", status='OLD', action = 'READ') 
  READ(158,*) logT_HG 
  CLOSE(UNIT = 158)

  OPEN(UNIT = 158,file=TRIM(ADJUSTL(path_to_file))//"Input/Chabrier/HG23_Vmix_Smix.csv", status='OLD', action = 'READ') 

  ! Skip header
  READ(158,*)

  READ(158,*) table_HG

  CLOSE(UNIT = 158)

  Vmix_HG = table_HG(3,1:)
  Smix_HG = table_HG(4,1:)


    !------

  END SUBROUTINE read_HG23_corr







  !----------------------------------------------!
  !                                              !
  !          SUBROUTINE calc_parameters          !
  !                                              !
  !----------------------------------------------!
  !
  ! Initializes all input parameters by reading the values in specific files.
  !
  SUBROUTINE calc_parameters(n_lay,n_mat_max,n_mat,xin,cf_Mmol,cf_Z_eff, cf_rho_0, &
             cf_T_0,cf_K_0,cf_Kp_0,cf_gam_0,cf_q,cf_a_T,cf_b_T,cf_c_T,cf_a_P,cf_T_D0, &
             f_alloy,MgD,MgSi,Mmol,Z_eff,rho_0,T_0,K_0,Kp_0,gam_0,q,a_T,b_T,c_T,a_P,T_D0,xout)

    !------

    INTEGER, INTENT(IN)  :: n_lay, n_mat_max

    INTEGER,          INTENT(IN)           ::  n_mat(n_lay)                  ! Number of different materials per layer

    DOUBLE PRECISION,   INTENT(IN) :: xin(n_lay,n_mat_max)                      ! Mole fraction of layer's materials (input)

    DOUBLE PRECISION,  INTENT(IN)  ::              cf_Mmol(n_lay,n_mat_max),  & ! Coefficients for molar mass
                                                   cf_Z_eff(n_lay,n_mat_max), & ! Coefficients for effective atomic charge
                                                   cf_rho_0(n_lay,n_mat_max), & ! Coefficients for reference density
                                                   cf_T_0(n_lay,n_mat_max),   & ! Coefficients for reference temperature
                                                   cf_K_0(n_lay,n_mat_max),   & ! Coefficients for reference bulk modulus
                                                   cf_Kp_0(n_lay,n_mat_max),  & ! Coefficients for pressure derivative of bulk mod.
                                                   cf_gam_0(n_lay,n_mat_max), & ! Coefficients for reference Gruneisen parameter
                                                   cf_q(n_lay,n_mat_max),     & ! Coefficients for adiabatic power exponent
                                                   cf_a_T(n_lay,n_mat_max),   & ! Coefficients for 1st alpha parameter
                                                   cf_b_T(n_lay,n_mat_max),   & ! Coefficients for 2nd alpha parameter
                                                   cf_c_T(n_lay,n_mat_max),   & ! Coefficients for 3rd alpha parameter
                                                   cf_a_P(n_lay,n_mat_max),   & ! Coefficients for temp. derivative of bulk modulus
                                                   cf_T_D0(n_lay,n_mat_max)     ! Coefficients for reference Debye temperature

    DOUBLE PRECISION, INTENT(IN)                  :: f_alloy,             & ! Fraction of alloy in the core
                                                     MgD,                 & ! Mg number (Mg#)
                                                     MgSi                   ! Mg/Si fraction


    INTEGER                                :: l,        & ! Coefficients index
                                              k,        & ! Layer index 
                                              counter

    DOUBLE PRECISION                       :: sum_x,    & ! Sum of mole fractions in a layer
                                              sum_w,    & ! Sum of weight fractions in a layer
                                              inv_rho,  & ! Inverse of density
                                              sum_Vr,   & ! Sum of reduced volumes
                                              sum_VrtK, & ! Sum of reduced volumes times bulk moduli
                                              sum_VroK, & ! Sum of reduced volumes over bulk moduli
                                              K_V0,     & ! Voigt bulk modulus
                                              K_R0        ! Reuss bulk modulus

    DOUBLE PRECISION  :: x(n_lay,n_mat_max)                      ! Mole fraction of layer's materials (internal)

                                              
    DOUBLE PRECISION, DIMENSION(n_mat_max) :: w,        & ! Weight fraction of layer's materials
                                              Vr_0        ! "Reduced" reference volume of the materials

    DOUBLE PRECISION,  DIMENSION(n_lay), INTENT(OUT) :: Mmol,                & ! Molar mass [kg.mol-1]
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
  
    DOUBLE PRECISION,   INTENT(OUT) :: xout(n_lay,n_mat_max)                      ! Mole fraction of layer's materials (output)
 
    
    !------



    Mmol  = 0.d0
    Z_eff = 0.d0
    rho_0 = 0.d0
    T_0   = 0.d0
    K_0   = 0.d0
    Kp_0  = 0.d0
    gam_0 = 0.d0
    q     = 0.d0                                       
    a_T   = 0.d0
    b_T   = 0.d0
    c_T   = 0.d0
    a_P   = 0.d0
    T_D0  = 0.d0

    ! WRITE(6,*) 'cf_Mmol(1,:) = ', cf_Mmol(1,:)


    x = 0.d0
 
    DO k=1, n_lay

      ! WRITE(6,*) x(k,:) 
      ! WRITE(6,*) 'x (before) = ', x

      x(k,:) = xin(k,:)

      w     = 0.d0
      sum_x = 0.d0
      sum_w = 0.d0
      inv_rho = 0.d0
      Vr_0  = 0.d0
      sum_Vr = 0.d0
      sum_VrtK = 0.d0
      sum_VroK = 0.d0
      K_V0 = 0.d0
      K_R0 = 0.d0


     ! WRITE(6,*) 'x = ', x

     ! WRITE(6,*) 'sum_w =', sum_w
     ! WRITE(6,*) 'w = ', w
     ! WRITE(6,*) 'sum_x = ', sum_x


       ! Attribution of the mole fractions if asked
       IF (ANY(x(k,:)<-0.5d0) .AND. ANY(x(k,:)>0.d0)) THEN
          WRITE(6,*) '/!\ ERROR in subroutine "read_parameters": all mole fractions must be set to "-1".'
          STOP
       ELSE IF (ANY(x(k,:) < -0.5d0)) THEN

    !      IF (k==1) THEN ! If layer k is the core
    !        x(k,1) = 1.d0 - f_alloy ! Fe
    !        x(k,2) = f_alloy        ! FeS
    !        x(k,3) = 0.d0           ! FeO (absent from original model)
    !        x(k,4) = 0.d0           ! FeSi (absent from original model)
    !        x(k,5) = 0.d0           ! Fe7C3 (absent from original model)
    !      ELSE IF (k==2) THEN ! If layer k is the lower mantle
    !        x(k,1) = MgD/MgSi * MgD               ! MgSiO3
    !        x(k,2) = MgD/MgSi * (1.d0-MgD)        ! FeSiO3
    !        x(k,3) = (1.d0-MgD/MgSi) * MgD        ! MgO
    !        x(k,4) = (1.d0-MgD/MgSi) * (1.d0-MgD) ! FeO
    !      ELSE IF (k==3) THEN ! If layer k is the upper mantle

           IF (k==1) THEN ! If layer k is the core
             !x(k,1) = 2.d0*(1.d0-MgD/MgSi) * MgD         ! Mg2SiO4
             !x(k,2) = 2.d0*(1.d0-MgD/MgSi) * (1.d0-MgD) ! Fe2SiO4
             x(k,1) = 0.9d0     ! Mg2SiO4
             x(k,2) = 0.1d0     ! Fe2SiO4


    !      ELSE IF (k==4) THEN ! If layer k is the supercritical water layer
    !        x(k,1) = 1.d0 ! H2O supercritical
    !      ELSE IF (k==5) THEN ! If layer k is the void 




           ELSE IF (k==2) THEN ! If layer k is the envelope
             !x(k,1) = 2.d0*(1.d0-MgD/MgSi) * MgD         ! Mg2SiO4
             !x(k,2) = 2.d0*(1.d0-MgD/MgSi) * (1.d0-MgD) ! Fe2SiO4
             x(k,1) = 0.9d0     ! Mg2SiO4
             x(k,2) = 0.1d0     ! Fe2SiO4


           ELSE IF (k==3) THEN ! If layer k is the void 
             x(k,1) = 0.d0 ! Void
          END IF

       ELSE IF (ALL(x(k,:)>-0.5d0)) THEN
          DO l=1, n_mat(k)
             sum_x = sum_x + x(k,l)
          END DO
          IF (sum_x>1.d0-1.d-10 .AND. sum_x<1.d0+1.d-10) THEN
             DO l=1, n_mat(k)
                x(k,l) = x(k,l) / sum_x
             END DO
          ELSE IF (sum_x<1.d0 .OR. sum_x>1.d0) THEN
             WRITE(6,*) '/!\ ERROR in subroutine "read_parameters": the sum of all mole fractions is not equal to unity.'
             STOP
          END IF
       END IF

       IF (ANY(x(k,:)<0.d0) .OR. ANY(x(k,:)>1.d0)) THEN
          WRITE(6,*) '/!\ ERROR in subroutine "read_parameters": please check mole fractions or/and composition parameters.'
          STOP
       END IF


       ! Computation of the layer's parameters using the material coefficients
       ! Molar mass and effective atomic number
       DO l=1, n_mat(k)
          Mmol(k)  = Mmol(k)  + x(k,l)*cf_Mmol(k,l)
          Z_eff(k) = Z_eff(k) + x(k,l)*cf_Z_eff(k,l)
       END DO
       ! Reference density
       DO l=1, n_mat(k)
          sum_w = sum_w + cf_Mmol(k,l) * x(k,l)
          !WRITE(6,*) 'k = ',k
          !WRITE(6,*) 'l = ',l
          !WRITE(6,*) 'sum_w = ', sum_w
       END DO
       
       DO l=1, n_mat(k)
          w(l) = cf_Mmol(k,l) * x(k,l) / sum_w
          ! WRITE(6,*) 'k = ',k
          ! WRITE(6,*) 'l = ',l
          ! WRITE(6,*) 'w = ', w
       END DO

       DO l=1, n_mat(k)
          inv_rho = inv_rho + w(l)/cf_rho_0(k,l)
       END DO

       !WRITE(6,*) inv_rho
       !WRITE(6,*) w
       !WRITE(6,*) sum_w


       rho_0(k) = 1.d0 / inv_rho

       !IF (k==8) rho_0(k) = 0.d0 ! Case of the void
       !IF (k==5) rho_0(k) = 0.d0 ! Case of the void
       IF (k==3) rho_0(k) = 0.d0 ! Case of the void

       ! Reference bulk modulus
       DO l=1, n_mat(k)
          Vr_0(l) = cf_Mmol(k,l) * x(k,l) / cf_rho_0(k,l)

          sum_Vr = sum_Vr + Vr_0(l)

          sum_VrtK = sum_VrtK + Vr_0(l)*cf_K_0(k,l)
          sum_VroK = sum_VroK + Vr_0(l)/cf_K_0(k,l)
       END DO

       K_V0 = sum_VrtK / sum_Vr
       K_R0 = sum_Vr / sum_VroK

       K_0(k) = 0.5d0 * (K_V0 + K_R0)
       !WRITE(6,*) 'K_0 =', K_0
       ! Other parameters (simple approximation of their value for a mixture of materials)
       DO l=1, n_mat(k)
          T_0(k)   = T_0(k)   + x(k,l)*cf_T_0(k,l)
          Kp_0(k)  = Kp_0(k)  + x(k,l)*cf_Kp_0(k,l)
          gam_0(k) = gam_0(k) + x(k,l)*cf_gam_0(k,l)
          q(k)     = q(k)     + x(k,l)*cf_q(k,l)

          a_T(k)   = a_T(k)   + x(k,l)*cf_a_T(k,l)
          b_T(k)   = b_T(k)   + x(k,l)*cf_b_T(k,l)
          c_T(k)   = c_T(k)   + x(k,l)*cf_c_T(k,l)
          a_P(k)   = a_P(k)   + x(k,l)*cf_a_P(k,l)

          T_D0(k)  = T_D0(k)  + x(k,l)*cf_T_D0(k,l)
       END DO

     xout = x

  !WRITE(6,*) 'x (after) = ', x
       
    END DO

    !------
        
  END SUBROUTINE calc_parameters

END MODULE fort_input

