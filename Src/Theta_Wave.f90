!!!For wave nornal angle
subroutine Theta_wave(L,lambda,Theta_function)
        use constants
        implicit none

  real*8,intent(in)::L,lambda
  real*8,intent(out)::Theta_function(7)


  real*8::lambda_deg
  real*8::theta_min
  real*8::theta_max
  real*8::theta_m1
  real*8::theta_w1
  real*8::theta_m2
  real*8::theta_w2
  real*8::theta_q
  real*8::coeff(4,4),coeff_bb(4),theta_q1,theta_q2
  real*8::coeff_a(2,6),coeff_b(6)
  integer::i,j 
 lambda_deg=lambda/RAD




SELECT CASE(ThetaType)
    CASE('CONSTANT') 
         Theta_function=Theta_Input
    CASE('OBLIQUEHISS') 
         if(lambda_deg < 5) then
      theta_min=0.0
      theta_max=20.0
      theta_m1=0.0
      theta_w1=15.0
 else if (lambda_deg >= 5 .and. lambda_deg < 10) then
      theta_min=0.0
      theta_max=30.0
      theta_m1=10.0
      theta_w1=20.0
 else if (lambda_deg >= 10 .and. lambda_deg < 15) then
      theta_min=0.0
      theta_max=40.0
      theta_m1=20.0
      theta_w1=25.0
 else if (lambda_deg >= 15 .and. lambda_deg < 20) then
      theta_min=10.0
      theta_max=50.0
      theta_m1=30.0
      theta_w1=25.0
 else if (lambda_deg >= 20 .and. lambda_deg < 25) then
      theta_min=20.0
      theta_max=60.0
      theta_m1=40.0
      theta_w1=35.0
 else if (lambda_deg >= 25 .and. lambda_deg < 30) then
      theta_min=30.0
      theta_max=70.0
      theta_m1=50.0
      theta_w1=40.0
 else if (lambda_deg >= 30 .and. lambda_deg < 35) then
      theta_min=40.0
      theta_max=80.0
      theta_m1=60.0
      theta_w1=50.0
 else if (lambda_deg >= 35 .and. lambda_deg < 40) then
      theta_min=50.0
      theta_max=80.0
      theta_m1=70.0
      theta_w1=60.0
 else
      theta_min=55.0
      theta_max=85.0
      theta_m1=80.0
      theta_w1=70.0
end if
      theta_m2=0.0
      theta_w2=1.0
      theta_q=0.0
    Theta_function(6)=theta_max
        Theta_function(1)=theta_m1
    Theta_function(2)=theta_w1
    Theta_function(3)=theta_m2
    Theta_function(4)=theta_w2
    Theta_function(5)=theta_min
    Theta_function(6)=theta_max
    Theta_function(7)=theta_q

    CASE('OA2015')
	if(mlt.ge.3.and.mlt.le.15)then !dayside
		if(L.le.5)then
			coeff(1,1)=-2.96
			coeff(1,2)=5.43e-3
			coeff(1,3)=1.25e-4
			coeff(1,4)=-1.98e-6
			coeff(2,1)=-2.25e-2
			coeff(2,2)=2.57e-2
			coeff(2,3)=-4.71e-4
			coeff(2,4)=1.4e-6
			coeff(3,1)=4.96e-3
			coeff(3,2)=-1.57e-3
			coeff(3,3)=2.61e-5
			coeff(3,4)=-5.46e-8
			coeff(4,1)=-7.44e-5
			coeff(4,2)=2.39e-5
			coeff(4,3)=-3.9e-7
			coeff(4,4)=6.85e-10
		else
			coeff(1,1)=-2.33942
			coeff(1,2)=0.148424
			coeff(1,3)=-0.00442111
			coeff(1,4)=2.93478e-5
			coeff(2,1)=0.397584
			coeff(2,2)=-0.0260653
			coeff(2,3)=0.000500087
			coeff(2,4)=-2.39098e-6
			coeff(3,1)=-0.0200059
			coeff(3,2)=0.00102028
			coeff(3,3)=-1.37305e-5
			coeff(3,4)=3.316e-8
			coeff(4,1)=0.000272183
			coeff(4,2)=-1.1708e-5
			coeff(4,3)=1.08928e-7
			coeff(4,4)=1.04256e-11
		end if

        else
		if(L.le.5)then
			coeff(1,1)=-2.36144
			coeff(1,2)=-0.0679
			coeff(1,3)=0.003815125
			coeff(1,4)=-0.000043285
			coeff(2,1)=0.05889
			coeff(2,2)=0.00836
			coeff(2,3)=-0.00022328
			coeff(2,4)=0.000001204
			coeff(3,1)=-0.00201295
			coeff(3,2)=-0.00041678
			coeff(3,3)=0.000006725
			coeff(3,4)=4.719e-9
			coeff(4,1)=5.096e-5
			coeff(4,2)=4.968e-6
			coeff(4,3)=-5.058e-8
			coeff(4,4)=-5.518e-10
		else
			coeff(1,1)=-0.78803
			coeff(1,2)=-0.02769
			coeff(1,3)=1.479e-3
			coeff(1,4)=-4.1015e-5
			coeff(2,1)=0.0640415
			coeff(2,2)=0.0108776
			coeff(2,3)=-7.891e-4
			coeff(2,4)=1.328e-5
			coeff(3,1)=-5.297e-3
			coeff(3,2)=-5.60918e-4
			coeff(3,3)=4.53587e-5
			coeff(3,4)=-7.31029e-7
			coeff(4,1)=9.1149e-5
			coeff(4,2)=8.60196e-6
			coeff(4,3)=-6.94066e-7
			coeff(4,4)=1.07675e-8
		end if
        end if
	theta_q=0.0
	do i=0,3
		do j=0,3
		theta_q=theta_q+coeff(i+1,j+1)*(KP*10)**j*(lambda_deg)**i
		end do
	end do
	theta_q=theta_q**2
	coeff_a=0.0
	coeff_a(1,1)=11.54
	coeff_a(1,2)=14.32
	coeff_a(1,3)=-8.18
	coeff_a(1,4)=1.22
	coeff_a(2,1)=66.0
	coeff_a(2,2)=1.0
	coeff_b(1)=5.68
	coeff_b(2)=4.62
	coeff_b(3)=3.05
	coeff_b(4)=-5.06
	coeff_b(5)=1.83
	coeff_b(6)=-0.203
	theta_m1=0.0
	theta_m2=0.0
	theta_w1=0.0
	do i=0,5
		theta_m1=theta_m1+coeff_a(1,i+1)*(lambda_deg/10.0)**i
		theta_m2=theta_m2+coeff_a(2,i+1)*(lambda_deg/10.0)**i
		theta_w1=theta_w1+coeff_b(i+1)*(lambda_deg/10.0)**i
	end do
	theta_w2=theta_w1
	theta_min=0.0
	theta_max=85.0
    Theta_function(1)=theta_m1
    Theta_function(2)=theta_w1
    Theta_function(3)=theta_m2
    Theta_function(4)=theta_w2
    Theta_function(5)=theta_min
    Theta_function(6)=theta_max
    Theta_function(7)=theta_q


	CASE('VITALIICHORUS')
           if(lambda_deg>20.0)then
             if(mlt.ge.3.and.mlt.le.15)then !dayside
		if(L.le.5)then
			coeff(1,1)=-2.96
			coeff(1,2)=5.43e-3
			coeff(1,3)=1.25e-4
			coeff(1,4)=-1.98e-6
			coeff(2,1)=-2.25e-2
			coeff(2,2)=2.57e-2
			coeff(2,3)=-4.71e-4
			coeff(2,4)=1.4e-6
			coeff(3,1)=4.96e-3
			coeff(3,2)=-1.57e-3
			coeff(3,3)=2.61e-5
			coeff(3,4)=-5.46e-8
			coeff(4,1)=-7.44e-5
			coeff(4,2)=2.39e-5
			coeff(4,3)=-3.9e-7
			coeff(4,4)=6.85e-10
		else
			coeff(1,1)=-2.33942
			coeff(1,2)=0.148424
			coeff(1,3)=-0.00442111
			coeff(1,4)=2.93478e-5
			coeff(2,1)=0.397584
			coeff(2,2)=-0.0260653
			coeff(2,3)=0.000500087
			coeff(2,4)=-2.39098e-6
			coeff(3,1)=-0.0200059
			coeff(3,2)=0.00102028
			coeff(3,3)=-1.37305e-5
			coeff(3,4)=3.316e-8
			coeff(4,1)=0.000272183
			coeff(4,2)=-1.1708e-5
			coeff(4,3)=1.08928e-7
			coeff(4,4)=1.04256e-11
		end if

           else
		if(L.le.5)then
			coeff(1,1)=-2.36144
			coeff(1,2)=-0.0679
			coeff(1,3)=0.003815125
			coeff(1,4)=-0.000043285
			coeff(2,1)=0.05889
			coeff(2,2)=0.00836
			coeff(2,3)=-0.00022328
			coeff(2,4)=0.000001204
			coeff(3,1)=-0.00201295
			coeff(3,2)=-0.00041678
			coeff(3,3)=0.000006725
			coeff(3,4)=4.719e-9
			coeff(4,1)=5.096e-5
			coeff(4,2)=4.968e-6
			coeff(4,3)=-5.058e-8
			coeff(4,4)=-5.518e-10
		else
			coeff(1,1)=-0.78803
			coeff(1,2)=-0.02769
			coeff(1,3)=1.479e-3
			coeff(1,4)=-4.1015e-5
			coeff(2,1)=0.0640415
			coeff(2,2)=0.0108776
			coeff(2,3)=-7.891e-4
			coeff(2,4)=1.328e-5
			coeff(3,1)=-5.297e-3
			coeff(3,2)=-5.60918e-4
			coeff(3,3)=4.53587e-5
			coeff(3,4)=-7.31029e-7
			coeff(4,1)=9.1149e-5
			coeff(4,2)=8.60196e-6
			coeff(4,3)=-6.94066e-7
			coeff(4,4)=1.07675e-8
		end if
        end if
	theta_q=0.0
	do i=0,3
		do j=0,3
		theta_q=theta_q+coeff(i+1,j+1)*(KP*10)**j*(lambda_deg)**i
		end do
	end do
	theta_q=theta_q**2
	coeff_a=0.0
	coeff_a(1,1)=11.54
	coeff_a(1,2)=14.32
	coeff_a(1,3)=-8.18
	coeff_a(1,4)=1.22
	coeff_a(2,1)=66.0
	coeff_a(2,2)=1.0
	coeff_b(1)=5.68
	coeff_b(2)=4.62
	coeff_b(3)=3.05
	coeff_b(4)=-5.06
	coeff_b(5)=1.83
	coeff_b(6)=-0.203
	theta_m1=0.0
	theta_m2=0.0
	theta_w1=0.0
	do i=0,5
		theta_m1=theta_m1+coeff_a(1,i+1)*(lambda_deg/10.0)**i
		theta_m2=theta_m2+coeff_a(2,i+1)*(lambda_deg/10.0)**i
		theta_w1=theta_w1+coeff_b(i+1)*(lambda_deg/10.0)**i
	end do
				

     else

	if(mlt.ge.3.and.mlt.le.15)then !dayside
			coeff(1,1)=1.613e-2
			coeff(1,2)=-1.343e-2
			coeff(1,3)=3.497e-3
			coeff(1,4)=-2.766e-4
			coeff(2,1)=3.051e-2
			coeff(2,2)=-3.199e-2
			coeff(2,3)=1.05e-2
			coeff(2,4)=-9.992e-4
			coeff(3,1)=-5.596e-2
			coeff(3,2)=5.255e-2
			coeff(3,3)=-1.531e-2
			coeff(3,4)=1.318e-3
			coeff(4,1)=2.028e-2
			coeff(4,2)=-1.794e-2
			coeff(4,3)=4.905e-3
			coeff(4,4)=-3.958e-4

			coeff_bb(1)=5.926
			coeff_bb(2)=1.643e1
			coeff_bb(3)=3.666
			coeff_bb(4)=-9.124e-1

        else
			coeff(1,1)=2.49e-3
			coeff(1,2)=-2.252e-3
			coeff(1,3)=7.417e-4
			coeff(1,4)=-5.463e-5
			coeff(2,1)=5.711e-2
			coeff(2,2)=-6.312e-2
			coeff(2,3)=2.142e-2
			coeff(2,4)=-2.152e-3
			coeff(3,1)=-1.171e-1
			coeff(3,2)=1.217e-1
			coeff(3,3)=-3.868e-2
			coeff(3,4)=3.694e-3
			coeff(4,1)=4.599e-2
			coeff(4,2)=-4.62e-2
			coeff(4,3)=1.424e-2
			coeff(4,4)=-1.328e-3

			coeff_bb(1)=3.455
			coeff_bb(2)=9.938
			coeff_bb(3)=2.147
			coeff_bb(4)=-4.175e-1
        end if
	theta_q1=0.0
	do i=0,3
		do j=0,3
		theta_q1=theta_q1+coeff(i+1,j+1)*(lambda_deg/10)**i*L**j
		end do
	end do
	theta_q2=0.0
	do i=0,3
		theta_q2=theta_q2+coeff_bb(i+1)*Kp**i
	enddo
	theta_q=theta_q1*theta_q2/(1-theta_q1*theta_q2)
	coeff_a=0.0
	coeff_a(1,1)=11.54
	coeff_a(1,2)=14.32
	coeff_a(1,3)=-8.18
	coeff_a(1,4)=1.22
	coeff_a(2,1)=66.0
	coeff_a(2,2)=1.0
	coeff_b(1)=5.68
	coeff_b(2)=4.62
	coeff_b(3)=3.05
	coeff_b(4)=-5.06
	coeff_b(5)=1.83
	coeff_b(6)=-0.203
	theta_m1=0.0
	theta_m2=0.0
	theta_w1=0.0
	do i=0,5
		theta_m1=theta_m1+coeff_a(1,i+1)*(lambda_deg/10.0)**i
		theta_m2=theta_m2+coeff_a(2,i+1)*(lambda_deg/10.0)**i
		theta_w1=theta_w1+coeff_b(i+1)*(lambda_deg/10.0)**i
	end do
   endif
	theta_w2=theta_w1
	theta_min=0.0
	theta_max=85.0
        Theta_function(1)=theta_m1
    Theta_function(2)=theta_w1
    Theta_function(3)=theta_m2
    Theta_function(4)=theta_w2
    Theta_function(5)=theta_min
    Theta_function(6)=theta_max
    Theta_function(7)=theta_q
	CASE('NBBLBC')
		if(lambda_deg.le.5)then
			theta_m1=0
			theta_w1=30
			theta_min=0
			theta_max=50
		else if(lambda.le.10)then
			theta_m1=15
                        theta_w1=30
                        theta_min=0
                        theta_max=50
		else
			theta_m1=30
                        theta_w1=30
                        theta_min=0
                        theta_max=50
		endif	
		theta_m2=0
		theta_w2=1
		theta_q=0
    Theta_function(1)=theta_m1
    Theta_function(2)=theta_w1
    Theta_function(3)=theta_m2
    Theta_function(4)=theta_w2
    Theta_function(5)=theta_min
    Theta_function(6)=theta_max
    Theta_function(7)=theta_q
	CASE('NBBUBC')
                if(lambda_deg.le.5)then
                        theta_m1=0
                        theta_w1=30
                        theta_min=0
                        theta_max=45
                else
                        theta_m1=30
                        theta_w1=30
                        theta_min=0
                        theta_max=45
                endif   
        theta_m2=0
        theta_w2=1
        theta_q=0
Theta_function(1)=theta_m1
    Theta_function(2)=theta_w1
    Theta_function(3)=theta_m2
    Theta_function(4)=theta_w2
    Theta_function(5)=theta_min
    Theta_function(6)=theta_max
    Theta_function(7)=theta_q
    CASE DEFAULT
        WRITE(*,*) 'Invalid Theta Type.'
        STOP
    END SELECT
  return


end subroutine Theta_wave
