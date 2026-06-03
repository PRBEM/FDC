! find the resontant frequencies

subroutine Roots_finder(alpha,E,L,lambda,nres,theta,roots)

use constants
use variables
implicit none

  real*8,intent(in)::alpha,E,L,lambda,theta
  integer,intent(in)::nres
  real*8,intent(out)::roots(4)

  real*8::coeff(11)
  real*8::rtr(10)
  real*8::rti(10)

  real*8::coeff1(4)
  real*8::rtr1(3)
  real*8::rti1(3)


  real*8::a7
  real*8::a8
  real*8::a9
  real*8::a10
  real*8::a11
  real*8::a12
  real*8::a13
  real*8::a14
  real*8::a15
  real*8::a16



  real*8::b1
  real*8::b2



  real*8::gamma1
  real*8::betta
  real*8::alpha_ne
   




  real*8::X



  integer::i
  integer::NN
  real*8::Temp,Temp2

  betta=sqrt(E*(E+2.0))/(E+1.0)
  gamma1=E+1.0



  call Local_alpha(alpha,lambda,alpha_ne)

  X=tan(theta)

							

  !write(*,*)lambda

  a1=omega_i**2.0+omega_e**2.0+omega_pe**2.0+omega_pi**2.0
  a2=(omega_i*omega_e)**2.0+(omega_pe*omega_i)**2.0+(omega_pi*omega_e)**2.0
  a3=omega_i+real(sqe)*omega_e
  a4=omega_i*real(sqe)*omega_e-omega_pe**2.0-omega_pi**2.0
  !a5=-((omega_pi**2.0)*real(sqe)*omega_e+(omega_pe**2.0)*omega_i)
  a5=0.0
  a6=omega_pi**2.0+omega_pe**2.0
  a7=(X**2.0)*a1+omega_i**2.0+omega_e**2.0+a6
  a8=(X**2.0)*a2+(omega_i*omega_e)**2.0+a6*(omega_i**2.0+omega_e**2.0)
  a9=a6*((omega_i*omega_e)**2.0)
  a10=2.0*a4-a3**2.0-a6
  a11=a4**2.0-2.0*a3*a5-2.0*a4*a6+(a3**2.0)*a6
  a12=a5**2.0+(a4**2.0)*a6-2.0*a3*a5*a6
  a13=(a5**2.0)*a6
  !write(*,*)a13
  a14=(X**2.0)*(2.0*a4-a3**2.0)-(2.0+X**2.0)*(a1+a6)
  a15=(X**2.0)*(a4**2.0-2.0*a3*a5)+(2.0+X**2.0)*(a2+a1*a6)
  a16=(X**2.0)*(a5**2.0)+(2.0+X**2.0)*a2*a6

  b1=real(nres)*real(sqe)*omega_e/gamma1
  b2=1.0/(betta*cos(alpha_ne)*cos(theta))

  !write(*,*)"parameters=",a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,b1,b2
    
  coeff(11)=(1+X**2.0)*((1-b2**2.0)**2.0)
  coeff(10)=4.0*(1+X**2.0)*b1*(b2**2.0)*(1-b2**2.0)
  coeff(9)=(6.0*(1+X**2.0)*(b1**2.0)-a7)*(b2**4.0)-(a14+2.0*(1+X**2.0)*(b1**2.0))*(b2**2.0)+(1+X**2.0)*a10
  coeff(8)=(4.0*a7*b1-4.0*(1+X**2.0)*(b1**3.0))*(b2**4.0)+2.0*a14*b1*(b2**2.0)
  coeff(7)=(a8-6.0*a7*(b1**2.0)+(1+X**2.0)*(b1**4.0))*(b2**4.0)-(a15+a14*(b1**2.0))*(b2**2.0)+(1+X**2.0)*a11
  coeff(6)=(4.0*a7*(b1**3.0)-4.0*a8*b1)*(b2**4.0)+2.0*a15*b1*(b2**2.0)
  coeff(5)=(6.0*a8*(b1**2.0)-a7*(b1**4.0)-a9)*(b2**4.0)-(a15*(b1**2.0)-a16)*(b2**2.0)-(1+X**2.0)*a12
  coeff(4)=(4.0*a9*b1-4.0*a8*(b1**3.0))*(b2**4.0)-2.0*a16*b1*(b2**2.0)
  coeff(3)=(a8*(b1**4.0)-6.0*a9*(b1**2.0))*(b2**4.0)+a16*(b1**2.0)*(b2**2.0)+(1+X**2.0)*a13
  coeff(2)=4.0*a9*(b1**3.0)*(b2**4.0)
  coeff(1)=-a9*(b1**4.0)*(b2**4.0)

 !write(*,*)"Equation Coefficients=",coeff
  !write(*,*)b1,b2
  if(nres/=0) then
    !coeff(11)=coeff(11)/coeff(11)
    !coeff(10)=coeff(10)/coeff(11)
    !coeff(9)=coeff(9)/coeff(11)
    !coeff(8)=coeff(8)/coeff(11)
    !coeff(7)=coeff(7)/coeff(11)
    !coeff(6)=coeff(6)/coeff(11)
    !coeff(5)=coeff(5)/coeff(11)
    !coeff(4)=coeff(4)/coeff(11)
    !coeff(3)=coeff(3)/coeff(11)
    !coeff(2)=coeff(2)/coeff(11)
    !coeff(1)=coeff(1)/coeff(11)
    !write(*,*)coeff
    
	!write(*,*) " finding roots"
    !call zrhqr(coeff,10,rtr,rti)   ! use the subroutine from Numerical Recipes in Fortran (1992)

	call zrhqr2(coeff,rtr,rti)       ! use LAPACK library function DGEEV
  
    !write(*,*)nres,omega_lc,omega_uc,rtr
  
    roots(1)=0.0
    roots(2)=0.0
    roots(3)=0.0
	roots(4)=0.0

    NN=0

    do i=1,10
       if(rti(i)==0) then
	     if((rtr(i)>omega_lc).and.(rtr(i)<omega_uc)) then
	       Temp2=abs((rtr(i)-b1)*b2/rtr(i))
		   if(Temp2>1.0) then
	         NN=NN+1
		     !write(*,*)NN,Temp,omega_lc,omega_uc
		     roots(NN)=rtr(i)
		   end if
	     end if
	   end if
    end do

    !write(*,*)omega_lc
  
  else
    !coeff1(7)=coeff(11)
	!coeff1(6)=coeff(10)
	!coeff1(5)=coeff(9)
	coeff1(4)=coeff(11)
	coeff1(3)=coeff(9)
	coeff1(2)=coeff(7)
	coeff1(1)=coeff(5)

	!write(*,*)coeff1
	!write(*,*) " finding roots"

	!call zrhqr(coeff1,3,rtr1,rti1)
	call zrhqr3(coeff1,rtr1,rti1)        ! use LAPACK library function DGEEV

    roots(1)=0.0
    roots(2)=0.0
    roots(3)=0.0
	roots(4)=0.0

    NN=0

    do i=1,3
       if(rti1(i)==0) then
	     if((rtr1(i)>(omega_lc**2.0)).and.(rtr1(i)<(omega_uc**2.0)).and.(abs(b2)>1.0)) then
		 !if((rtr1(i)>(omega_lc**2.0)).and.(rtr1(i)<(omega_uc**2.0))) then
	       NN=NN+1
		   !write(*,*)NN,Temp,omega_lc,omega_uc
		   roots(NN)=sqrt(rtr1(i))
	     end if
	   end if
    end do
  
  end if

  !write(*,*)NN
  !write(*,*)" finished finding roots "
  !write(100,*)roots

  !close(100)

  return
  
end subroutine Roots_finder  

	   
     
