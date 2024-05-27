
MODULE fort_spec

  USE fort_spec_func


  IMPLICIT NONE

  
CONTAINS






subroutine calc_radius(struc_len,press,gravity,rho,P0_cgs, &
     R_pl,var_grav,radius)

  implicit none
  ! I/O
  INTEGER, intent(in)                         :: struc_len
  DOUBLE PRECISION, intent(in)                :: P0_cgs
  DOUBLE PRECISION, intent(in)                :: press(struc_len), &
       rho(struc_len)
  DOUBLE PRECISION, intent(in)                :: gravity, R_pl
  LOGICAL, intent(in)                         :: var_grav
  DOUBLE PRECISION, intent(out)               :: radius(struc_len)

  ! Internal
  INTEGER                                     :: i_str
  DOUBLE PRECISION                            :: R0, inv_rho(struc_len) !, integ_parab

  inv_rho = 1d0/rho

  radius = 0d0
  R0=0d0
  if (var_grav) then

     !write(*,*) '####################### VARIABLE GRAVITY'
     !write(*,*) '####################### VARIABLE GRAVITY'
     !write(*,*) '####################### VARIABLE GRAVITY'
     !write(*,*) '####################### VARIABLE GRAVITY'

     ! Calculate radius with vertically varying gravity, set up such that at P=P0, i.e. R=R_pl
     ! the planet has the predefined scalar gravity value
     do i_str = struc_len-1, 1, -1
        if ((press(i_str+1) > P0_cgs) .AND. (press(i_str) <= P0_cgs)) then
           if (i_str <= struc_len-2) then
              R0 = radius(i_str+1) - integ_parab(press(i_str),press(i_str+1),press(i_str+2), &
                   inv_rho(i_str),inv_rho(i_str+1),inv_rho(i_str+2),P0_cgs,press(i_str+1))/gravity &
                   /R_pl**2d0
           else
              R0 = radius(i_str+1)-(1d0/rho(i_str)+1d0/rho(i_str+1))/(2d0*gravity)* &
                   (press(i_str+1)-P0_cgs)/R_pl**2d0
           end if
        end if
        if (i_str <= struc_len-2) then
           radius(i_str) = radius(i_str+1) - integ_parab(press(i_str),press(i_str+1),press(i_str+2), &
                inv_rho(i_str),inv_rho(i_str+1),inv_rho(i_str+2),press(i_str),press(i_str+1))/gravity &
                /R_pl**2d0
        else
           radius(i_str) = radius(i_str+1)-(1d0/rho(i_str)+1d0/rho(i_str+1))/(2d0*gravity)* &
                (press(i_str+1)-press(i_str))/R_pl**2d0
        end if
     end do
     R0 = 1d0/R_pl -R0
     radius = radius + R0
     radius = 1d0/radius

  else

     !write(*,*) '####################### CONSTANT GRAVITY'
     !write(*,*) '####################### CONSTANT GRAVITY'
     !write(*,*) '####################### CONSTANT GRAVITY'
     !write(*,*) '####################### CONSTANT GRAVITY'


     ! Calculate radius with vertically constant gravity
     do i_str = struc_len-1, 1, -1
        if ((press(i_str+1) > P0_cgs) .AND. (press(i_str) <= P0_cgs)) then
           if (i_str <= struc_len-2) then
              R0 = radius(i_str+1) + integ_parab(press(i_str),press(i_str+1),press(i_str+2), &
                   inv_rho(i_str),inv_rho(i_str+1),inv_rho(i_str+2),P0_cgs,press(i_str+1))/gravity
           else
              R0 = radius(i_str+1)+(1d0/rho(i_str)+1d0/rho(i_str+1))/(2d0*gravity)* &
                   (press(i_str+1)-P0_cgs)
           end if
        end if
        if (i_str <= struc_len-2) then
           radius(i_str) = radius(i_str+1) + integ_parab(press(i_str),press(i_str+1),press(i_str+2), &
                inv_rho(i_str),inv_rho(i_str+1),inv_rho(i_str+2),press(i_str),press(i_str+1))/gravity
        else
           radius(i_str) = radius(i_str+1)+(1d0/rho(i_str)+1d0/rho(i_str+1))/(2d0*gravity)* &
                (press(i_str+1)-press(i_str))
        end if
     end do

     R0 = R_pl-R0
     radius = radius + R0
!!$     write(*,*) R0, P0_cgs, gravity, R_pl, press(20), rho(20)

  end if

end subroutine calc_radius




END MODULE fort_spec
