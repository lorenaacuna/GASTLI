!      /'''''''''''''''''''''''''''\
!     /          CONSTANTS          \
!     \............................./

! This module lists some useful physical and mathematical constants for the code.

MODULE constants

  IMPLICIT NONE

                                    ! Universal constants
  DOUBLE PRECISION, PARAMETER    :: pi        = 3.141592d0,       & ! Value of Pi
                                    R_g       = 8.3144621d0,      & ! Universal gas constant [J.mol-1.K-1]
                                    N_A       = 6.022d23,         & ! Avogadro number [mol-1]
                                    G_grav    = 6.67408d-11,      & ! Gravitational constant [N.m2.kg-2]
                                    unit_AU   = 1.49597870700d11, & ! Astronomical Unit [m]
                                    unit_day  = 8.64d4,           & ! Day [s]
                                    M_Su      = 1.9884d30,        & ! Solar mass [kg]
                                    R_Su      = 6.96342d8,        & ! Solar radius [m]
                                    T_eff_Su  = 5.772d3,          & ! Solar effective temperature [K]
                                    M_E       = 5.97237d24,       & ! Earth mass [kg]
                                    R_E       = 6.371d6,          & ! Earth radius [m] (mean value)
                                    rho_E     = 5.513d3,          & ! Earth mean density [kg.m-3]
                                    x_core_E  = 0.325d0,          & ! Core mass fraction of Earth
                                    x_H2O_E   = 0.0005d0,         & ! Water mass fraction of Earth
                                    f_alloy_E = 0.13d0,           & ! Fraction of alloy in the Earth
                                    MgD_E     = 0.9d0,            & ! Mg number (Mg#) of the Earth
                                    MgSi_E    = 1.131d0,          & ! Mg/Si ratio of the Earth, corrected from solar
                                    FeSi_E    = 0.986d0,          & ! Fe/Si ratio of the Earth, corrected from solar
                                    C_iner_tot_E = 8.034d37,      & ! Polar moment of inertia [kg.m2] (Lambeck 1980)
                                    Alb_E     = 0.306d0,          & ! Earth bond albedo (De Pater & Lissauer 2015)
                                    D_orb_E   = 1.d0,             & ! Earth orbital distance [AU]
                                    P_orb_E   = 365.25696d0,      & ! Earth orbital period [day]
                                    P_surf_E  = 1.d5,             & ! Earth surface pressure [Pa] (1 bar = 1e5 Pa)
                                    T_surf_E  = 288d0,            & ! Earth surface temperature [K] (288 K)
                                    T_eq_E    = 255d0               ! Earth equilibrium temperature [K] (255 K)

                                    ! Thermodynamical constants (input)
  DOUBLE PRECISION, DIMENSION(6) :: Ttp,                          & ! Temperature of water's triple points [K]
                                    Ptp                             ! Pressure of water's triple points [Pa]


END MODULE constants
