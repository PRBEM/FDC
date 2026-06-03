! local diffusion rate at one latitude

subroutine D_local_n(alpha,E,L,nres,lambda,Bw,Theta_function,DDN)
use constants
implicit none

  real*8,intent(in)::alpha,E,L,lambda,Bw,Theta_function(7)
  integer,intent(in)::nres
  real*8,intent(out)::DDN(3)


  real*8::WeightValue(3)

  real*8::theta_min
  real*8::theta_max
 ! real*8::theta_m
  !real*8::theta_w
  real*8::theta

  !common /group2/ theta_min,theta_max,theta_m,theta_w

  integer::nint_tmp

  real*8::BB1
  real*8::BB2
  real*8::BB3

  integer::i

  !theta_m=Theta_function(1)
  !theta_w=Theta_function(2)
  theta_min=Theta_function(5)
  theta_max=Theta_function(6)


!write(*,*)Theta_max

  BB1=(theta_max-theta_min)*RAD
  BB2=1.0/(real(nint_X)-1.0)
  BB3=BB1*BB2

  DDN(1)=0.0
  DDN(2)=0.0
  DDN(3)=0.0

!    if(lambda*180.0/3.1415926>44.6 .and. theta*180.0/PI>43.0) then
!                write(*,*)theta_max
!                write(*,*)theta_min
!      end if

  !write(*,*)Theta_function(2)

  do i=1,nint_X
     theta=theta_min/180.0*PI+(real(i)-1.0)*BB3

	 WeightValue(1)=0.0
	 WeightValue(2)=0.0
	 WeightValue(3)=0.0

     call WeightFunction(alpha,E,L,lambda,Bw,Theta_function,nres,theta,WeightValue)
	 
	 WeightValue(1:3)=WeightValue(1:3)*sin(theta)
	 if(i==1) then
	   DDN(1)=DDN(1)+WeightValue(1)/2.0
	   DDN(2)=DDN(2)+WeightValue(2)/2.0
	   DDN(3)=DDN(3)+WeightValue(3)/2.0
	 else if(i==nint_X) then
	        DDN(1)=DDN(1)+WeightValue(1)/2.0
	        DDN(2)=DDN(2)+WeightValue(2)/2.0
	        DDN(3)=DDN(3)+WeightValue(3)/2.0
		  else
		    DDN(1)=DDN(1)+WeightValue(1)
	        DDN(2)=DDN(2)+WeightValue(2)
	        DDN(3)=DDN(3)+WeightValue(3)
	 end if
  end do



  DDN(1:3)=DDN(1:3)*BB3

  return

end subroutine D_local_n

