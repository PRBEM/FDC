! obtain the mirror latitude

subroutine Lambda_mirror_new(alpha,lambda_m)

	implicit none

	real*8,intent(in)::alpha
	real*8,intent(out)::lambda_m

	real*8::BB
  real*8::P(7)
  real*8::rtr(6)
  real*8::rti(6)
  real*8::root

  !real*8::AA

  integer::i

  !open(112,file='mirrorlatitude2.txt')

  BB=sin(alpha)


  root=0.0

  P(7)=1.0
  P(6)=0.0
  P(5)=0.0
  P(4)=0.0
  P(3)=0.0
  P(2)=3.0*(BB**4.0)
  P(1)=-4.0*(BB**4.0)


!write(*,*)P
  !call zrhqr(P,6,rtr,rti)
  call zrhqr1(P,rtr,rti)

!write(*,*)rtr

  do i=1,6
     if(rti(i)==0) then
	   if((rtr(i)>=0.0).and.(rtr(i)<=1.0)) then
         root=rtr(i)
       end if
     end if
  end do

 !write(*,*)root

  lambda_m=acos(sqrt(root))*0.999
  
  !write(*,*)lambda_m*180.0/3.1415926
 !write(*,*)alpha/PI*180.0


  return

end subroutine Lambda_mirror_new

