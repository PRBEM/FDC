! get the local pitch angle from equatorial pitch-angle

subroutine Local_alpha(alpha,lambda,alpha_ne)

implicit none
  
  real*8,intent(in)::alpha,lambda
  real*8,intent(out)::alpha_ne
  real*8::f_value

  real*8::temp

  call f_function(lambda,f_value)
  
  temp=sqrt(f_value)*sin(alpha)

  alpha_ne=asin(temp)
 
  !write(*,*)"A=",temp,"asin(A)=",alpha_ne

  return

end subroutine Local_alpha
