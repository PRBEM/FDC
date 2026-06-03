! the function over lambda for subsequent integration to get b-av diffusion coefficients

subroutine int_D_loc(alpha,E,L,nres,lambda,D_integrad)

implicit none

  real*8,intent(in)::alpha,E,L,lambda
  real*8,intent(out)::D_integrad(3)


  real*8::Bw
  real*8::Theta_function(7)
  
  real*8::alpha_ne
  real*8::ds
  real*8::D_local(3)

  integer::nres
 
  call Local_alpha(alpha,lambda,alpha_ne)   ! the local pitch angle

  call Bw_wave(L,lambda,Bw)
  call Theta_wave(L,lambda,Theta_function)

  call D_local_n(alpha,E,L,nres,lambda,Bw,Theta_function,D_local)        ! the local diffusion rates

!  D_integrad(1)=D_local(1)*cos(alpha_ne)*((cos(lambda))**7.0)/((cos(alpha))**2.0)
!  D_integrad(2)=D_local(2)*((cos(lambda))**4.0)*((1.0+3.0*((sin(lambda))**2.0))**0.25)/cos(alpha)
!  D_integrad(3)=D_local(3)*cos(lambda)*((1.0+3.0*((sin(lambda))**2.0))**0.5)/cos(alpha_ne)


  D_integrad(1)=D_local(1)/cos(alpha_ne)*((tan(alpha)/tan(alpha_ne))**2.0)*cos(lambda)*(1.0+3.0*((sin(lambda))**2.0))**0.5
  D_integrad(2)=D_local(2)*tan(alpha)/sin(alpha_ne)*cos(lambda)*(1.0+3.0*((sin(lambda))**2.0))**0.5
  D_integrad(3)=D_local(3)*cos(alpha)*cos(lambda)*(1.0+3.0*((sin(lambda))**2.0))**0.5
  
  return

end subroutine int_D_loc
