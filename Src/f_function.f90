! the variation of magnetic field with magnetic latitude
! refer to equation (3.9) in book Basic Space Plasma Physics 

subroutine f_function(lambda,f_value)

implicit none

  real*8,intent(in)::lambda
  real*8,intent(out)::f_value

  f_value=sqrt(4.0-3.0*((cos(lambda))**2.0))/((cos(lambda))**6.0)

  return

end subroutine f_function
