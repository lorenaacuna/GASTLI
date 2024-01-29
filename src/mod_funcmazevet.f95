!      /'''''''''''''''''''''''''''\
!     /          FUNCmazevet        \
!     \............................./

! This module contains the function used for Mazevet+2019 EOS.

MODULE funcmazevet


  IMPLICIT NONE

CONTAINS


  !---------------------------------------!
  !                                       !
  !          FUNCTION FERMIF              !
  !                                       !
  !---------------------------------------!
  !
  ! Fermi function
  !

  FUNCTION FERMIF(X) ! Fermi function

     !------

    DOUBLE PRECISION  :: FERMIF, & ! Result
			 F

    DOUBLE PRECISION, INTENT(IN)  :: X

     !------

      IF (X > 40.d0) THEN
         F = 0.d0
      ELSE IF (X < -40.d0) THEN
         F = 1.d0
      ELSE
         F = 1.d0/(dexp(X)+1.d0)
      END IF

      FERMIF = F
      RETURN

     !------
      
      END FUNCTION FERMIF












END MODULE funcmazevet
