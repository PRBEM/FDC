!!! the normalization factor N(omega)

subroutine int_N_omega(alpha,E,L,lambda,Theta_function,theta_max_eff,nres,omega,res_cone,N_of_omega)
use constants
use variables
implicit none

  real*8,intent(in)::alpha,E,L,lambda,Theta_function(7),theta_max_eff,omega
  integer,intent(in)::nres
  real*8,intent(out)::N_of_omega
  real*8::Ntmp
  
  character::res_cone

  real*8::theta_min
  real*8::theta_max
  real*8::theta_m1
  real*8::theta_w1
  real*8::theta_m2
  real*8::theta_w2
  real*8::theta_q

  
  real*8::X_w1
  real*8::X_m1
  real*8::X_w2
  real*8::X_m2




  real*8::gamma1
  real*8::betta
   


  real*8::theta
  real*8::X
  real*8::step

  real*8::a7
  real*8::a8

  real*8::a14
  real*8::a15
  real*8::a16

  real*8::k_vec
  real*8::k_vec1

  real*8::n_ref_ind2,n_ref_ind2_0,n_ref_ind2_1,n_ref_ind2_2,n_ref_ind2_1_temp,n_ref_ind2_3

  real*8::PP0
  real*8::PP1
  real*8::DD1
  real*8::DD2
  real*8::dDdk,dPdw,dSdw,dRdw,dLdw,dF1dw,dF2dw,dF3dw,dDdw
  real*8::dw_dk

  real*8::PD



  integer::nint_N_tmp
  integer::i

  betta=sqrt(E*(E+2.0))/(E+1.0)
  gamma1=E+1.0

  if(res_cone=='Y') then
     nint_N_tmp=nint_N_res_cone
     !write(*,*)Theta_function(2)
  else !if (res_cone=='N') then
     nint_N_tmp=nint_N
  end if

  theta_m1=Theta_function(1)
  theta_w1=Theta_function(2)
  theta_m2=Theta_function(3)
  theta_w2=Theta_function(4)
  theta_min=Theta_function(5)
  theta_max=theta_max_eff !Theta_function(2)
  theta_q=Theta_function(7)

  X_m1=tan(theta_m1/180.0*PI)
  X_w1=tan(theta_w1/180.0*PI)
  X_m2=tan(theta_m2/180.0*PI)
  X_w2=tan(theta_w2/180.0*PI)

  step=(theta_max-theta_min)/180.0*PI/(real(nint_N_tmp)-1.0)

  N_of_omega=0

  Ntmp=0

  do i=1,nint_N_tmp
     theta=theta_min/180.0*PI+step*(real(i)-1.0)
	 !write(*,*)theta*180.0/PI
	 X=tan(theta)


	 a7=(X**2.0)*a1+omega_i**2.0+omega_e**2.0+a6
     a8=(X**2.0)*a2+(omega_i*omega_e)**2.0+a6*(omega_i**2.0+omega_e**2.0)
    
     a14=(X**2.0)*(2.0*a4-a3**2.0)-(2.0+X**2.0)*(a1+a6)
     a15=(X**2.0)*(a4**2.0-2.0*a3*a5)+(2.0+X**2.0)*(a2+a1*a6)
     a16=(X**2.0)*(a5**2.0)+(2.0+X**2.0)*a2*a6

	 		
		
		! modification by xudong 20111024
		! getting the refrective index from the full dispersion relation
!		PD=D_stix!*P_stix
!
!		if (PD<=0) then
!			n_ref_ind2=((R_stix*L_stix*(X**2.0)+(2.0+X**2.0)*P_stix*S_stix)+&
!					  sqrt((R_stix*L_stix*(X**2.0)+(2.0+X**2.0)*P_stix*S_stix)**2.0-&
!					  4.0*(S_stix*(X**2.0)+P_stix)*P_stix*R_stix*L_stix*(1.0+X**2.0)))/(2.0*(S_stix*(X**2.0)+P_stix))
!		else !if (PD<=0) then
			n_ref_ind2=((R_stix*L_stix*(X**2.0)+(2.0+X**2.0)*P_stix*S_stix)-&
					  sqrt((R_stix*L_stix*(X**2.0)+(2.0+X**2.0)*P_stix*S_stix)**2.0-&
					  4.0*(S_stix*(X**2.0)+P_stix)*P_stix*R_stix*L_stix*(1.0+X**2.0)))/(2.0*(S_stix*(X**2.0)+P_stix))
!		end if
!
!		if (n_ref_ind2<0) then
!		write(*,*)n_ref_ind2
!		end if

		k_vec=sqrt(n_ref_ind2)*omega/c
        !write(*,*)k_vec

		if (n_ref_ind2<0) then
			!write(*,*)"the infractive from Glauert index is negative,WF"
			!stop
			Ntmp=0
		!end if	
		else
	
	   ! getting the dw/dk from the full dispersion relation
       
		!Revised by Xudong in 20111012
					!!!!!Start!!!!!!!
						
						!!dD/dk
						
						dDdk=(S_stix*(X**2.0)+P_stix)*(c**4.0)/(omega**4.0)*4.0*(k_vec**3.0)-&
							 (R_stix*L_stix*(X**2.0)+(2.0+X**2.0)*P_stix*S_stix)*(c**2.0)/(omega**2.0)*2.0*k_vec
						
							
						!!dD/dw

						dPdw=2.0*a6/(omega**3.0)

						dSdw=(2.0*(a1-omega_i**2.0-omega_e**2.0)*(omega**5.0)+&
							 4.0*((omega_i**2.0)*(omega_e**2.0)-a2)*(omega**3.0)+&
							 2.0*(a2*(omega_i**2.0)+a2*(omega_e**2.0)-a1*(omega_i**2.0)*(omega_e**2.0))*omega)/&
							 (((omega**2.0-omega_i**2.0)**2.0)*((omega**2.0-omega_e**2.0)**2.0))

						dRdw=(2.0*(omega_i*real(sqe)*omega_e-a4)*(omega**3.0)+&
							 (a3*omega_i*real(sqe)*omega_e-a4*omega_i-a4*real(sqe)*omega_e-3.0*a5)*(omega**2.0)-&
							 2.0*a5*(omega_i+real(sqe)*omega_e)*omega-&
							 a5*omega_i*real(sqe)*omega_e)/&
							 ((omega**2.0)*((omega+omega_i)**2.0)*((omega+real(sqe)*omega_e)**2.0))

						dLdw=(2.0*(omega_i*real(sqe)*omega_e-a4)*omega**3.0-&
							 (a3*omega_i*real(sqe)*omega_e-a4*omega_i-a4*real(sqe)*omega_e-3.0*a5)*(omega**2.0)-&
							 2.0*a5*(omega_i+real(sqe)*omega_e)*omega+&
							 a5*omega_i*real(sqe)*omega_e)/&
							 ((omega**2.0)*((omega-omega_i)**2.0)*((omega-real(sqe)*omega_e)**2.0))

						dF1dw=((k_vec*c)**4.0)/(omega**5.0)*(omega*(X**2.0)*dSdw+omega*dPdw-4.0*(s_stix*(X**2.0)+P_stix))
						dF2dw=((k_vec*c)**2.0)/(omega**3.0)*(omega*((X**2.0)*(L_stix*dRdw+R_stix*dLdw)+(2.0+X**2.0)*(P_stix*dSdw+S_stix*dPdw))-&
							  2.0*(R_stix*L_stix*(X**2.0)+(2.0+X**2.0)*P_stix*S_stix))
						dF3dw=(1.0+X**2.0)*(R_stix*L_stix*dPdw+P_stix*L_stix*dRdw+P_stix*R_stix*dLdw)

						dDdw=dF1dw-dF2dw+dF3dw
					!!!!!End!!!!!!!!!
						
	 dw_dk=-dDdk/dDdw
	
	 Ntmp=sin(theta)*(exp(-((X-X_m1)/X_w1)**2)+theta_q*exp(-((X-X_m2)/X_w2)**2))*n_ref_ind2*c/abs(dw_dk)

	 end if

	 if(i==1) then
	        N_of_omega=N_of_omega+Ntmp/2.0
	 else if(i==nint_N_tmp) then
	        N_of_omega=N_of_omega+Ntmp/2.0
		  else
		    N_of_omega=N_of_omega+Ntmp
	 end if
  end do

  N_of_omega=N_of_omega*step

  return

end subroutine int_N_omega




     

     
