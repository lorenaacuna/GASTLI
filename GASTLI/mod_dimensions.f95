!      /''''''''''''''''''''''''''''\
!     /          DIMENSIONS          \
!     \............................../

! This module lists the different lengths used for arrays in the code (i.e. their dimensions),
! which may be changed if needed.

MODULE dimensions

  IMPLICIT NONE

  INTEGER, PARAMETER :: n_lay     = 3,     & ! Number of layers in the planet (plus void)
                        n_EOS     = 10,     & ! Number of available EOS
                        n_pts     = 10000, & ! Number of points (steps for the spatial integration)
                        n_int     = 100,   & ! Number of points used for integrations not on space
                        n_mat_max = 5        ! Maximum number of different materials in a layer

END MODULE dimensions
