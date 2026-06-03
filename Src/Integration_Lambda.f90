! Integration over latitude

subroutine Integration_lambda(lambda_m,lambda_min,nint_tmp,B_m,alpha,E,L,nres,BavRate)
use constants
implicit none

  real*8,intent(in)::alpha,E,L,B_m,lambda_min,lambda_m
  integer,intent(in)::nint_tmp,nres
  real*8,intent(out)::BavRate(3)


  real*8::psi(nint_tmp+1)
  real*8::Bfield(nint_tmp+2)
  real*8::dB
  real*8::B_local

  real*8::lambda_tmp
  real*8::delta_lambda   ! step of magnetic latitude in rad
  real*8::lambda

  integer::i,j



  real*8::D_integrad(3)


  delta_lambda=(lambda_m-lambda_min)/real(nint_tmp)     ! step of magnetic latitude in rad

  BavRate(1)=0.0
  BavRate(2)=0.0
  BavRate(3)=0.0
  if(delta_lambda.gt.0.0)then
    do i=1,nint_tmp+1

        lambda_tmp=lambda_min+delta_lambda*real(i-1)

        call b_field(L,lambda_tmp,B_local)
        Bfield(i)=B_local

        if(i==1) then
            psi(i)=PI/2.0
        else if (i==nint_tmp+1 .and. B_local>B_m) then
            psi(i)=0.0
        else
            psi(i)=asin(sqrt(1-B_local/B_m)/cos(alpha))
        end if

    end do
	call b_field(L,lambda_m+delta_lambda,B_local)
    Bfield(nint_tmp+2)=B_local

    do j=1,nint_tmp+1
       lambda=lambda_min+delta_lambda*real(j-1)              ! magnetic latitude

       call int_D_loc(alpha,E,L,nres,lambda,D_integrad)    ! the integrad for integration over lambda

       if(j==1) then
             BavRate(1)=BavRate(1)+D_integrad(1)/2.0
             BavRate(2)=BavRate(2)+D_integrad(2)/2.0
             BavRate(3)=BavRate(3)+D_integrad(3)*B_m/Bfield(j)*tan(alpha)**2

       else if(j==nint_tmp+1) then

              dB=(Bfield(j+1)-Bfield(j-1))/2
              BavRate(1)=BavRate(1)+D_integrad(1)/2.0
              BavRate(2)=BavRate(2)+D_integrad(2)/2.0
              BavRate(3)=BavRate(3)+D_integrad(3)*2*B_m/abs(dB)*cos(psi(j))*abs(psi(j-1)-psi(j))


      else
              dB=(Bfield(j+1)-Bfield(j-1))/2
              BavRate(1)=BavRate(1)+D_integrad(1)
              BavRate(2)=BavRate(2)+D_integrad(2)
              BavRate(3)=BavRate(3)+D_integrad(3)*2*B_m/abs(dB)*cos(psi(j))*abs(psi(j-1)-psi(j+1))

       end if
    end do

    BavRate(1)=BavRate(1)*delta_lambda
    BavRate(2)=BavRate(2)*delta_lambda
    BavRate(3)=BavRate(3)*delta_lambda/2
	endif

  return


end subroutine Integration_lambda
