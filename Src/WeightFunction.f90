subroutine WeightFunction(alpha,E,L,lambda,Bw,Theta_function,nres,theta,WeightValue)
use constants
use variables
implicit none
  real*8,intent(in)::alpha,E,L,lambda,Bw,Theta_function(7),theta
  integer,intent(in)::nres
  real*8,intent(out)::WeightValue(3)
  real*8::D_tmp(3)
  real*8::lambda_eq

  real*8::theta_min
  real*8::theta_max
  real*8::theta_m1
  real*8::theta_w1
  real*8::theta_m2
  real*8::theta_w2
  real*8::theta_q

  real*8::theta_cone
  real*8::theta_max_eff
  character::res_cone

  real*8::gamma1
  real*8::betta
  real*8::alpha_ne
   
  real*8::omega_e_eq

  real*8::B_local
  real*8::B_eq

  real*8::p
  real*8::vpar

  real*8::X

  real*8::roots(4)

  integer::i
  
  real*8::omega
  real*8::k_vec

  real*8::n_ref_ind2,n_ref_ind2_0,n_ref_ind2_1,n_ref_ind2_2,n_ref_ind2_1_temp,n_ref_ind2_3
  
  real*8::PP0
  real*8::PP1
  real*8::DD1
  real*8::DD2
  real*8::dDdk,dPdw,dSdw,dRdw,dLdw,dF1dw,dF2dw,dF3dw,dDdw
  real*8::dw_dk

  real*8::argJ

  real*8::BesselF1
  real*8::BesselF2
  real*8::BesselF3

  real*8::F1
  real*8::F2
  real*8::F3
  real*8::F4
  real*8::F5
  real*8::F6
  real*8::Fsq

  real*8::dw_dk_par

  real*8::Delta_n
  real*8::G1
  real*8::G2
  real*8::G3
  real*8::G4

  real*8::N_of_omega
  real*8::X_m1
  real*8::X_w1
  real*8::X_m2
  real*8::X_w2

  real*8::feq             ! ratio of equatorial plasma frequency to gyrofrequency

  real*8::dB_local
  real*8::PSD

!real*8::erf1,erf2,mag1,mag2

  real*8::PD
  real*8::Ne_local

  D_tmp(1)=0.0
  D_tmp(2)=0.0
  D_tmp(3)=0.0

  WeightValue(1)=0.0
  WeightValue(2)=0.0
  WeightValue(3)=0.0

  theta_m1=Theta_function(1)
  theta_w1=Theta_function(2)
  theta_m2=Theta_function(3)
  theta_w2=Theta_function(4)
  theta_min=Theta_function(5)
  theta_max=Theta_function(6)
  theta_q=Theta_function(7)

  !write(*,*)Theta_function(2)

  call Local_alpha(alpha,lambda,alpha_ne)
  call ElectronDensity(L,lambda,Ne_local)
  
  X=tan(theta)
  X_m1=tan(theta_m1/180.0*PI)
  X_w1=tan(theta_w1/180.0*PI)
  X_m2=tan(theta_m2/180.0*PI)
  X_w2=tan(theta_w2/180.0*PI)

  betta=sqrt(E*(E+2.0))/(E+1.0)
  gamma1=E+1.0

  lambda_eq=0.0
  call B_field(L,lambda,B_local)
  call B_field(L,lambda_eq,B_eq)

  omega_e=qe*B_local/me
  omega_e_eq=qe*B_eq/me
  omega_i=omega_e*me/mi
  omega_pe=sqrt(Ne_local*qe*qe/(C0*me))
  omega_pi=omega_pe*sqrt(me/mi)
  
  call Freq_Spectral(L,lambda,omega_e_eq,omega_m,d_omega,omega_lc,omega_uc)
  
!mag1=(omega_m-omega_lc)/d_omega
!  mag2=(omega_uc-omega_m)/d_omega
  
  !write(*,*)"mag1=",mag1,"mag2=",mag2

 ! call errorfunction(mag1,erf1)
 ! call errorfunction(mag2,erf2)



p=me*c*(sqrt((E+1.0)**2.0-1.0))
  vpar=p/gamma1/me*cos(alpha_ne)


  call Roots_finder(alpha,E,L,lambda,nres,theta,roots)   

  do i=1,4
     if(roots(i)==0.or.roots(i).le.sqrt((omega_pe/omega_e)**2+1.0)) then
           WeightValue(1)=WeightValue(1)
       WeightValue(2)=WeightValue(2)
           WeightValue(3)=WeightValue(3)
     else
           omega=roots(i)
           ! wave number
           k_vec=(omega-(real(nres)*real(sqe)*omega_e)/gamma1)/(c*betta*cos(alpha_ne)*cos(theta))
                        !!!!!test by xudong 20111020
                                if ((nres>0).and.(k_vec<0)) then
                                        write(*,*)"this root is error"
                                        write(*,*)"N=",nres,"Root=",omega
                                        stop
                                end if
                        !!!!!
           ! Stix parameter
           R_stix=1.0-omega_pe**2.0/(omega*(omega+real(sqe)*omega_e))-omega_pi**2.0/(omega*(omega+omega_i))
           L_stix=1.0-omega_pe**2.0/(omega*(omega-real(sqe)*omega_e))-omega_pi**2.0/(omega*(omega-omega_i))
           P_stix=1.0-(omega_pe/omega)**2.0-(omega_pi/omega)**2.0
           S_stix=0.5*(R_stix+L_stix)
           D_stix=0.5*(R_stix-L_stix)

   res_cone='N'
   theta_max_eff=theta_max ! degree

   !write(*,*)theta_max_eff

!      if(lambda*180.0/3.1415926>44.6 .and. theta*180.0/PI>43.0) then
!                write(*,*)theta*180.0/PI
!                !write(*,*)theta_max
!      end if

   if (P_stix/S_stix<0) then
      theta_cone=atan(sqrt(-P_stix/S_stix))*180.0/PI ! degree

      if(theta_cone*0.9<theta_max .and. theta_cone*0.99>theta_max) then       ! compare the resonance cone angle and X_max and choose the smaller one
          res_cone='Y'

      else if(theta_cone*0.99<=theta_max) then
          theta_max_eff=theta_cone*0.99
          !Theta_function(2)=theta_max_eff ! degree
          res_cone='Y'
          !write(*,*)theta_cone
          !write(*,*)Theta_function(2)
      end if

    end if

                ! modification by xudong 20111024
                !getting the refrective index from the full dispersion relation
!		PD=D_stix!*P_stix
!
!		if (PD<=0) then
!			n_ref_ind2=((R_stix*L_stix*(X**2.0)+(2.0+X**2.0)*P_stix*S_stix)+&
!					  sqrt((R_stix*L_stix*(X**2.0)+(2.0+X**2.0)*P_stix*S_stix)**2.0-&
!					  4.0*(S_stix*(X**2.0)+P_stix)*P_stix*R_stix*L_stix*(1.0+X**2.0)))/(2.0*(S_stix*(X**2.0)+P_stix))
!		else
                !if (PD<=0) then
                        n_ref_ind2=((R_stix*L_stix*(X**2.0)+(2.0+X**2.0)*P_stix*S_stix)-&
                                          sqrt((R_stix*L_stix*(X**2.0)+(2.0+X**2.0)*P_stix*S_stix)**2.0-&
                                          4.0*(S_stix*(X**2.0)+P_stix)*P_stix*R_stix*L_stix*(1.0+X**2.0)))/(2.0*(S_stix*(X**2.0)+P_stix))
                !else
                        !write(*,*) "the infractive index is error,WF"
                        !stop
!		end if

                !if (n_ref_ind2<0) then
                        !write(*,*)"the infractive from Glauert index is negative,WF"
                        !stop
                !end if

        !write(*,*)D_stix
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

                                                !dw_dk=-dDdk/dDdw
                                        !!!!!End!!!!!!!!!	
        
                dw_dk=-dDdk/dDdw
                                

           argJ=abs(k_vec)*c*(sqrt((E+1.0)**2.0-1.0))*sin(theta)*sin(alpha_ne)/(real(sqe)*omega_e)

           call Bessel_FirstKind(argJ,nres-1,BesselF1)                    
           call Bessel_Firstkind(argJ,nres,BesselF2)
           call Bessel_FirstKind(argJ,nres+1,BesselF3)

           F1=BesselF3*(n_ref_ind2-L_stix)/(n_ref_ind2-S_stix)+BesselF1*(n_ref_ind2-R_stix)/(n_ref_ind2-S_stix)
           F2=(n_ref_ind2*((sin(theta))**2.0)-P_stix)/(2.0*n_ref_ind2)
           F3=k_vec/abs(k_vec)*1.0/tan(alpha_ne)*sin(theta)*cos(theta)*BesselF2
           F4=(R_stix-L_stix)/(2.0*(n_ref_ind2-S_stix))
           F5=(P_stix-n_ref_ind2*((sin(theta))**2.0))/n_ref_ind2
           F6=P_stix*cos(theta)/n_ref_ind2
           Fsq=abs(((F1*F2+F3)**2.0)/((F4*F5)**2.0+F6**2.0))
 
           dw_dk_par=dw_dk/cos(theta)

           if(vpar==dw_dk_par .or. (theta*180.0/PI)  > theta_max_eff .or. theta_min > theta_max_eff) then
             D_tmp(1)=0.0
                 D_tmp(2)=0.0
                 D_tmp(3)=0.0

                 !write(*,*)theta_max!_eff !theta*180.0/PI

           else
             Delta_n=PI/2.0*1.0/cos(theta)/((vpar/c)**3.0)*Fsq*((real(nres)*real(sqe)*omega_e/(gamma1*omega)-&
             (sin(alpha_ne))**2.0)**2.0)/abs(1.0-dw_dk_par/vpar)

        call PowerSpectralDensity(L,Lambda,omega_e_eq,omega,PSD) 
                G1=omega_e*PSD

! G1=2.0*omega_e*exp(-((omega-omega_m)/d_omega)**2.0)/(sqrt(PI)*d_omega*(erf1+erf2))
		 call int_N_omega(alpha,E,L,lambda,Theta_function,theta_max_eff,nres,omega,res_cone,N_of_omega)
	
         G2=(exp(-((X-X_m1)/X_w1)**2.0)+theta_q*exp(-((X-X_m2)/X_w2)**2.0))/N_of_omega

		 dB_local=Bw
	     !give the special value
!					 if((Lambda*180.0/PI)<=10.0) then
!						dB_local=100 !pT
!					 else if ((Lambda*180.0/PI)>10.0.and.(Lambda*180.0/PI)<20 ) then
!						dB_local=50 !pT
!					 else if ((Lambda*180.0/PI)>=20.and.(Lambda*180.0/PI)<=45 ) then
!						dB_local=10.0**(0.75+0.04*(Lambda*180.0/PI))
!							if (dB_local>200) then
!								dB_local=200
!							end if
!							if (dB_local<50) then
!								dB_local=50
!							end if

!					 else
!						stop
!					 end if
				!if((Lambda*180.0/PI)<=20.0) then
				!	dB_local=10**(-0.015*(Lambda*180.0/PI)+1.98)
				!		if(dB_local<(10.0**(0.75+0.04*(Lambda*180.0/PI)))) dB_local=10.0**(0.75+0.04*(Lambda*180.0/PI))
					
				!else if ((Lambda*180.0/PI)>20.0.and.(Lambda*180.0/PI)<=45 ) then
				!	dB_local=10.0**(0.75+0.04*(Lambda*180.0/PI))
				!	if (dB_local<10**(-0.015*(Lambda*180.0/PI)+1.98)) dB_local=10**(-0.015*(Lambda*180.0/PI)+1.98)
				!	if (dB_local>150) dB_local=150
				!else
				!	stop
				!end if

			!dB_local=dB_local*1.0e-12

		 !write(*,*) dB_local*1e12, dB
          
		 G3=omega_e/(gamma1**2.0)*((dB_local/B_local)**2.0)

		 G4=sin(alpha_ne)*cos(alpha_ne)/(real(nres)*real(sqe)*omega_e/(gamma1*omega)-((sin(alpha_ne))**2.0))

		 D_tmp(1)=G3*Delta_n*G1*G2
		 D_tmp(2)=D_tmp(1)*G4
		 D_tmp(3)=D_tmp(2)*G4

       end if

	   !write(*,*)"G1=",G1,"Delta_n=",Delta_n,"G2=",G2,"G3=",G3
       
	   WeightValue(1)=WeightValue(1)+D_tmp(1)
       	   WeightValue(2)=WeightValue(2)+D_tmp(2)
	   WeightValue(3)=WeightValue(3)+D_tmp(3)

     end if

  end do

  return

end subroutine WeightFunction

	   


  
  


  
