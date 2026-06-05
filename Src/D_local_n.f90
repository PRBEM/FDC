! ==============================================================================
! SUBROUTINE: D_local_n
! ==============================================================================
!
! PURPOSE:
!   Computes the local quasi-linear diffusion coefficients (D_aa, D_ap, D_pp)
!   for a single resonance harmonic at a fixed magnetic latitude (lambda), by
!   integrating the diffusion coefficient integrand over wave normal angle theta.
!
!   The integration uses the trapezoidal rule over a uniform grid in theta,
!   from theta_min to theta_max (taken from Theta_function). The integrand at
!   each theta point is provided by WeightFunction.
!
!   The result is the local (non-bounce-averaged) diffusion rate. Bounce-
!   averaging over latitude and summation over resonance harmonics are
!   performed by the calling routine.
!
!
! INPUTS:
!   alpha          - Equatorial electron pitch angle [radians]
!   E              - Electron kinetic energy normalised to rest-mass energy
!                    (E = gamma - 1, dimensionless)
!   L              - Magnetic L-shell [Earth radii, dimensionless]
!   lambda         - Magnetic latitude [radians]
!   Bw             - RMS wave magnetic field amplitude [T]
!   Theta_function - Array(7) defining the wave normal angle distribution:
!                      (1) theta_m1  : mean angle of 1st Gaussian [degrees]
!                      (2) theta_w1  : width of 1st Gaussian [degrees]
!                      (3) theta_m2  : mean angle of 2nd Gaussian [degrees]
!                      (4) theta_w2  : width of 2nd Gaussian [degrees]
!                      (5) theta_min : lower cutoff of integration range [degrees]
!                      (6) theta_max : upper cutoff of integration range [degrees]
!                      (7) theta_q   : relative amplitude of 2nd Gaussian
!   nres           - Resonance harmonic number n (integer, e.g. -1, 0, 1, ...)
!
! OUTPUTS:
!   DDN(1)         - Local D_aa at this latitude [s^-1] (pitch-angle diffusion)
!   DDN(2)         - Local D_ap at this latitude [s^-1] (mixed diffusion)
!   DDN(3)         - Local D_pp at this latitude [s^-1] (momentum diffusion)
!
! CALLS:
!   WeightFunction - Diffusion coefficient integrand at a single theta value
!
! MODULE DEPENDENCIES:
!   constants      - Physical constants and unit conversion factors (RAD, PI, ...)
!   (nint_X is read from a module: number of theta integration points)
!
! INTEGRATION METHOD:
!   Trapezoidal rule on a uniform grid of nint_X points in theta [radians].
!   The integrand is weighted by sin(theta) before integration, consistent
!   with integration over wave normal angle in spherical geometry.
!
! NOTES:
!   - theta_min and theta_max are read in degrees from Theta_function(5:6)
!     and converted to radians internally.
!   - nint_X is a global parameter from a module rather than an argument;
!     this is a known issue to be addressed in future refactoring.

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

