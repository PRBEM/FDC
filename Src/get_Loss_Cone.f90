! get Loss Cone by xudong,201112
! Assumping the magnet is dipolar

subroutine get_Loss_cone(L, Alpha_LC)
use constants,only:RAD
! L, L-shells
!Alpha_LC, the pitch angle at the edge of Loss cone, unit in Deg
implicit none
real*8,intent(in)::L
real*8,intent(out)::Alpha_LC
real*8::temp, sinAlpha_lc_2, sinAlpha_lc


temp=(L**3.0)*sqrt(4.0-3.0/L)
	if(temp==0.0) then
		write(*,*)"the Alpha_LC is error!"
		stop
	end if
sinAlpha_lc_2=1.0/temp

sinAlpha_lc=sqrt(sinAlpha_lc_2)
Alpha_LC=asin(sinAlpha_lc)
! Rad to Deg
Alpha_LC=Alpha_LC/RAD

return

end subroutine get_Loss_Cone
