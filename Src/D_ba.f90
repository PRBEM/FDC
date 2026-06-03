! bounce-averaging

subroutine D_ba(alpha,E,L,nres,BavRate)
use constants
implicit none



  real*8,intent(in)::alpha,E,L
  integer,intent(in)::nres
  real*8,intent(out)::BavRate(3)

  real*8::B_m
  real*8::B_local
  real*8::lambda_tmp
  
  integer::nint_tmp


  real*8::lambda_min
  real*8::lambda_max


  real*8::lambda_m
  real*8::lambda_mdeg



  integer::i,j


  
  real*8::lambda
  real*8::D_integrad(3)

  real*8::A   ! 1/S(alpha), where S(alpha) is proportional to bounce period
  real*8::delta_lambda   ! step of magnetic latitude in rad

  
  A=1.0/(1.30-0.56*sin(alpha))          ! 1/S(alpha)

  call Lambda_mirror_new(alpha,lambda_m)    ! determine the mirror latitude
	call B_field(L,lambda_m,B_m)

  ! determine the range of magnetic latitude
  lambda_min=lam_min*RAD                 ! in rad
  lambda_mdeg=lambda_m/RAD                ! in degree

  if(lambda_mdeg>=lam_max) then                ! mirror latitude larger than maximum latitude
    lambda_m=lam_max*RAD
  else if(lambda_mdeg<=lam_min) then           ! mirror latitude smaller than minimum latitude
         lambda_m=lam_min*RAD
  else
	     lambda_m=lambda_m
  end if
  

  if(nres==0) then
  	nint_tmp=nint_zero_res
  else
        nint_tmp=nint
  end if


  call Integration_lambda(lambda_m,lambda_min,nint_tmp,B_m,alpha,E,L,nres,BavRate)

    BavRate(1)=BavRate(1)*A
    BavRate(2)=BavRate(2)*A
    BavRate(3)=BavRate(3)*A

  return

end subroutine D_ba
