!!!! Integration over latitude is in a separate function now Integration_latitude.f90
!!!! <Dpp> is calculated using an approach from Orlova and Shprits [2011]

!!!! The code computes matrices of the diffusion coefficients using realistic hiss waves parameters

!! Grids:
!!!! Distance L from 1 till 6 with step 0.1
!!!! Kinetic energy Ek from 100 keV till 10 MeV, 46 points
!!!! Equatorial pitch angle alpha from 0.6 till 89.6 with step 1, 89 points
!!!! Kp=4 (for other values of Kp the diffusion coefficients can be scaled using analytical function
!!!! Daa(Kp)=Daa(Kp=4)*10^{g(Kp)-g(Kp=4)}, where g(Kp)=0.2607*Kp - 0.01547*Kp*Kp)

!! Realistic hiss waves and plasma parameters:
!!!! Hiss waves parameters are taken from Orlova, Spasojevic and Shprits [GRL, 2014]:
!!!! Wave amplitude is dependent on L, magnetic latitude, and Kp;
!!!! Plasmaspheric density is increasing with latitude following Denton et al.[2004, 2006];
!!!! Gaussian spectrum parameters are delta_f=300.00 Hz, f_m=550.00 Hz,
!!!! f_lc=100.00 Hz, f_uc=2000.00 Hz;
!!!! Wave normal angle is increasing with latitude following Thorne et al.[2013];
!!!! The number of resonances is dependent on distance and energy following
!!!! Mourenas and Ripoll [2012];

!! Notes about code changes:
!!!! Added abs(k_vec) to the argument of Bessel function and sign (k_vec) to F3 in WeightFunction subroutine
!!!! Corrected bug in integration N(omega) in int_N_omega subroutine. It now computes k_vec from the
!!!! refractive index. Before it used k_vec from resonance condition which is wrong.
!!!! The code integrates 3 times. Now you can define the number of points for each integral.
!!!! Before the number of points for integration was the same for each integral.
!!!! Put the summation over resonances in the main file (to change the number of points in integration
!!!! over latitude for n=0)
!!!! Added Theta_wave(L,Kp,lambda,Bw) function to int_D_loc subroutine
!!!! Added checking the resonance cone. If theta is near to the resonance cone integral N(omega)
!!!! is computed with large number of points (to "solve" accurately singularity). If theta is larger than
!!!! the resonance cone we compute two integrals over theta till resonance cone*0.99
!!!! Added Bw_wave(L,Kp,lambda,Bw) function to int_D_loc subroutine. For lambda>25  we take
!!!! the same amplitude as for lambda=25
!!!! Added Denton et al., 2004; 2006 density model for plasmasphere which is dependent on latitude
!!!! The total number of resonances is now computed using equation from Ripoll and Mourenas, 2012
!!!! The number of points in integral over latitude is dependent on latitude.
!!!! The code uses LAPACK library Function
!!!! Diffusion coefficietns are computed using Albert [2007 JGR] method

subroutine sum_n(InputData2,local_rows2,columns,OutputData,noutput_rows,noutput_columns)
use constants
  implicit none

  integer,intent(in)::local_rows2
  integer,intent(in)::columns

  integer,intent(in)::noutput_rows, noutput_columns

  real*8,intent(in)::InputData2(local_rows2,columns)
  real*8,intent(out)::OutputData(noutput_rows,noutput_columns)


  real*8::BavDaa
  real*8::BavDap
  real*8::BavDaE
  real*8::BavDpp
  real*8::BavDEE

  
  real*8::BavRate(3)
    
  
  real*8::L
  real*8::E
  real*8::alpha_eq


  integer::i
  integer::resonmax,nres,I_nres


  do i=1,local_rows2
     L=InputData2(i,1)
     E=InputData2(i,2)*1.0e3*qe/(me*c*c)
     alpha_eq=InputData2(i,3)*RAD
     BavDaa=0.0
  BavDap=0.0
  BavDaE=0.0
  BavDpp=0.0
  BavDEE=0.0
     


     if(L<=55.0) then
   select case(NresType)
	case('CONSTANT')

	case('MR12')
	    call resonanceorder(E,L,resonmax)
            nres_min=-resonmax
	    nres_max= resonmax
        CASE DEFAULT
            WRITE(*,*) 'Invalid resonance Type.'
            STOP
    END SELECT
     do I_nres=nres_min, nres_max
       nres=I_nres
       call D_ba(alpha_eq,E,L,nres,BavRate)
       BavDaa=BavRate(1)+BavDaa
       BavDap=BavRate(2)+BavDap
       BavDpp=BavRate(3)+BavDpp
     enddo

   BavDaE=BavDap*(E+2.0)/(E+1.0)
   BavDEE=BavDpp*((E+2.0)**2.0)/((E+1.0)**2.0)
       
     ! now in unit day^{-1}
   BavDaa=BavDaa*24.0*3600.0
   BavDap=BavDap*24.0*3600.0
   BavDpp=BavDpp*24.0*3600.0
   BavDaE=BavDaE*24.0*3600.0
   BavDEE=BavDEE*24.0*3600.0
     else
        BavDaa=0.0
   BavDap=0.0
   BavDpp=0.0
   BavDaE=0.0
   BavDEE=0.0
     end if
     
   OutputData(i,1)=L
   OutputData(i,2)=InputData2(i,2)
   OutputData(i,3)=InputData2(i,3)
   OutputData(i,4)=BavDaa
   OutputData(i,5)=BavDpp
   OutputData(i,6)=BavDap      
  end do
           
end subroutine sum_n
