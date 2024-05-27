
MODULE fort_spec_func

  IMPLICIT NONE

  
CONTAINS


function integ_parab(x,y,z,fx,fy,fz,a,b)

  implicit none
  ! I/O
  double precision :: x,y,z,fx,fy,fz,a,b
  double precision :: integ_parab
  ! Internal
  double precision :: c1,c2,c3

  c3 = ((fz-fy)/(z-y)-(fz-fx)/(z-x))/(y-x)
  c2 = (fz-fx)/(z-x)-c3*(z+x)
  c1 = fx-c2*x-c3*x**2d0

  integ_parab = c1*(b-a)+c2*(b**2d0-a**2d0)/2d0+c3*(b**3d0-a**3d0)/3d0

end function integ_parab




END MODULE fort_spec_func
