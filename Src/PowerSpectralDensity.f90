!Distribution of wave power with frequency
subroutine PowerSpectralDensity(L,Lambda,omega_e_eq,Freq,PSD)
use constants
use variables
implicit none
real*8,intent(in)::L,Lambda,Freq,omega_e_eq
real*8,intent(out)::PSD
real*8::lambda_deg,mag1,mag2,erf1,erf2
real*8::c1,c2,kappa,tau,Bw_norm
real*8::beta0,beta1,beta2,beta3,beta4,beta5,beta6,beta7
lambda_deg=lambda/rad

SELECT CASE(FreqDistribution)
        CASE('GAUSSIAN')
                mag1=(omega_m-omega_lc)/d_omega
                mag2=(omega_uc-omega_m)/d_omega
  
                call errorfunction(mag1,erf1)
                call errorfunction(mag2,erf2)
                PSD=2.0*exp(-((Freq-omega_m)/d_omega)**2.0)/(sqrt(PI)*d_omega*(erf1+erf2))

        CASE('2016HISS')
                if(L.le.2.5)then
                        kappa=0.495
                        tau=5.2e4
                        c1=7.1e-7
                        c2=2.7e-6
                        Bw_norm=9.2446
                else if(L.gt.2.5.and.L.le.3.5)then
                        kappa=0.915
                        tau=1.6e5
                        c1=4.6e-7
                        c2=5.6e-7
                        Bw_norm=9.5703
                else if(L.gt.3.5.and.L.le.4.5)then
                        kappa=0.624
                        tau=1.6e5
                        c1=4.8e-7
                        c2=4.5e-7
                        Bw_norm=9.2774
                else
                        kappa=0.574
                        tau=2e5
                        c1=5.5e-7
                        c2=3.8e-7
                        Bw_norm=9.6841
                endif
                if (omega_uc<Freq .or. omega_lc>Freq) then
                        PSD=0
                else if (Freq<=omega_m .and. Freq>=omega_lc) then
                        PSD=c1*exp(-d_omega/Freq)/2/PI/Bw_norm/Bw_norm*1e6
                else
                        PSD=c2*(1+(Freq/2/PI)**2/kappa/tau)**(-kappa-1)/2/PI/Bw_norm/Bw_norm*1e6
                end if

        CASE('2017CHORUS')
        if(MLT.ge.6.and.MLT.lt.18)then !!Dayside
                if(lambda_deg.le.7)then
                        beta0=1.2608e-6
                        beta1=-4.6384e-5
                        beta2=6.2732e-4
                        beta3=-3.5303e-3
                        beta4=9.3343e-3
                        beta5=-1.0391e-2
                        beta6=5.7168e-4
                        beta7=4.973e-3
                else
                        beta0=5.974e-7
                        beta1=-2.2599e-5
                        beta2=3.3958e-4
                        beta3=-2.226e-3
                        beta4=7.985e-3
                        beta5=-1.6349e-2
                        beta6=1.7687e-2
                        beta7=-7.7332e-3
                endif
        else   !!!!Nightside
                if(lambda_deg.le.7)then
                        beta0=-1.2233e-6
                        beta1=5.1851e-5
                        beta2=-8.3298e-4
                        beta3=7.2506e-3
                        beta4=-3.3065e-2
                        beta5=8.1568e-2
                        beta6=-1.0443e-1
                        beta7=5.4606e-2
                else
                        beta0=1.2777e-7
                        beta1=-9.7504e-7
                        beta2=-2.4561e-5
                        beta3=6.7819e-4
                        beta4=-4.5157e-3
                        beta5=1.3066e-2
                        beta6=-1.7702e-2
                        beta7=9.2489e-3
                endif
        endif
                PSD=beta0+beta1*(Freq/omega_e_eq)+beta2*(Freq/omega_e_eq)**2+beta3*(Freq/omega_e_eq)**3+beta4*(Freq/omega_e_eq)**4+beta5*(Freq/omega_e_eq)**5+beta6*(Freq/omega_e_eq)**6+beta7*(Freq/omega_e_eq)**7
                erf1=beta0*(omega_lc/omega_e_eq)+beta1/2*(omega_lc/omega_e_eq)**2+beta2/3.0*(omega_lc/omega_e_eq)**3+beta3/4.0*(omega_lc/omega_e_eq)**4+beta4/5.0*(omega_lc/omega_e_eq)**5+&
                        beta5/6.0*(omega_lc/omega_e_eq)**6+beta6/7.0*(omega_lc/omega_e_eq)**7+beta7/8.0*(omega_lc/omega_e_eq)**8
                erf2=beta0*(omega_uc/omega_e_eq)+beta1/2*(omega_uc/omega_e_eq)**2+beta2/3.0*(omega_uc/omega_e_eq)**3+beta3/4.0*(omega_uc/omega_e_eq)**4+beta4/5.0*(omega_uc/omega_e_eq)**5+&
                        beta5/6.0*(omega_uc/omega_e_eq)**6+beta6/7.0*(omega_uc/omega_e_eq)**7+beta7/8.0*(omega_uc/omega_e_eq)**8
                PSD=PSD/((erf2-erf1)*omega_e_eq)
        CASE('2026UBCHORUS')
        if(MLT.ge.6.and.MLT.lt.18)then!!Dayside only for lat lower than 10
                        beta0=0.0014947221225619
                        beta1=-0.015803770764091
                        beta2=0.0711096808978984
                        beta3=-0.176587151860897
                        beta4=0.261515137534286
                        beta5=-0.231080159701795
                        beta6=0.112859681459824
                        beta7=-0.0235129368141491
        else !!!Nightside only for lat lower than 10 degree
                        beta0=0.0161805052165801
                        beta1=-0.177099131141267
                        beta2=0.825782709253675
                        beta3=-2.12703533294286
                        beta4=3.26969853152609
                        beta5=-3.0005692071183
                        beta6=1.52254312687049
                        beta7=-0.329626122022566
        endif
                PSD=beta0+beta1*(Freq/omega_e_eq)+beta2*(Freq/omega_e_eq)**2+beta3*(Freq/omega_e_eq)**3+beta4*(Freq/omega_e_eq)**4+beta5*(Freq/omega_e_eq)**5+beta6*(Freq/omega_e_eq)**6+beta7*(Freq/omega_e_eq)**7
                erf1=beta0*(omega_lc/omega_e_eq)+beta1/2*(omega_lc/omega_e_eq)**2+beta2/3.0*(omega_lc/omega_e_eq)**3+beta3/4.0*(omega_lc/omega_e_eq)**4+beta4/5.0*(omega_lc/omega_e_eq)**5+&
                        beta5/6.0*(omega_lc/omega_e_eq)**6+beta6/7.0*(omega_lc/omega_e_eq)**7+beta7/8.0*(omega_lc/omega_e_eq)**8
                erf2=beta0*(omega_uc/omega_e_eq)+beta1/2*(omega_uc/omega_e_eq)**2+beta2/3.0*(omega_uc/omega_e_eq)**3+beta3/4.0*(omega_uc/omega_e_eq)**4+beta4/5.0*(omega_uc/omega_e_eq)**5+&
                        beta5/6.0*(omega_uc/omega_e_eq)**6+beta6/7.0*(omega_uc/omega_e_eq)**7+beta7/8.0*(omega_uc/omega_e_eq)**8
                PSD=PSD/((erf2-erf1)*omega_e_eq)
                if(PSD.lt.1.0e-20)then
                        PSD=1.0e-20
                endif        
        CASE('2026LBCHORUS')
        if(MLT.ge.6.and.MLT.lt.18)then!!Dayside
                if(lambda_deg.le.10)then!!Dayside lowlatitude
                        beta0=1.11336540917184E-07
                        beta1=-3.72912612422969E-07
                        beta2=-9.16958991615852E-06
                        beta3=0.000215723693886208
                        beta4=-0.00126806058032446
                        beta5=0.00359190154486275
                        beta6=-0.00535827750817583
                        beta7=0.00335331572494158
                else if(lambda_deg.gt.10.and.lambda_deg.le.20)then!!Dayside midlatitude
                        beta0=4.30576556905264E-07
                        beta1=-1.93877149456913E-05
                        beta2=0.00035456331476089
                        beta3=-0.00247357343147156
                        beta4=0.00907129101805928
                        beta5=-0.0189341440278942
                        beta6=0.0211954576442875
                        beta7=-0.00980281543417116
                else!!Dayside highlatitude
                        beta0=1.01790061236938E-06
                        beta1=-4.82797987996845E-05
                        beta2=0.000929651327598179
                        beta3=-0.00771658195124899
                        beta4=0.0338877002399438
                        beta5=-0.0825886472235263
                        beta6=0.105301327489477
                        beta7=-0.0546805009108624
                endif
        else  !!!Nightside
                if(lambda_deg.le.10)then!!Nightside lowlatitude
                        beta0=-2.97874170950016E-07
                        beta1=1.4322737443405E-05
                        beta2=-0.000255070014505957
                        beta3=0.00266338023050004
                        beta4=-0.0144207301553904
                        beta5=0.0418529434765739
                        beta6=-0.0620248242140477
                        beta7=0.0366919071518278
                else if(lambda_deg.gt.10.and.lambda_deg.le.20)then!!Nightside midlatitude
                        beta0=1.87395736698176E-07
                        beta1=-7.4583549236649E-06
                        beta2=0.000115119097024855
                        beta3=-0.000712856986278894
                        beta4=0.00219744709855392
                        beta5=-0.00358152843005299
                        beta6=0.00290037927467941
                        beta7=-0.000878285849663509
                else!!Nightside highlatitude
                        beta0=1.41774960798836E-07
                        beta1=-4.18554459149626E-06
                        beta2=9.27712392661608E-05
                        beta3=-0.000844912962502075
                        beta4=0.0038042235566409
                        beta5=-0.00912206447572251
                        beta6=0.0112105133295341
                        beta7=-0.005566214822435
                endif
        endif
                PSD=beta0+beta1*(Freq/omega_e_eq)+beta2*(Freq/omega_e_eq)**2+beta3*(Freq/omega_e_eq)**3+beta4*(Freq/omega_e_eq)**4+beta5*(Freq/omega_e_eq)**5+beta6*(Freq/omega_e_eq)**6+beta7*(Freq/omega_e_eq)**7
                erf1=beta0*(omega_lc/omega_e_eq)+beta1/2*(omega_lc/omega_e_eq)**2+beta2/3.0*(omega_lc/omega_e_eq)**3+beta3/4.0*(omega_lc/omega_e_eq)**4+beta4/5.0*(omega_lc/omega_e_eq)**5&
                        +beta5/6.0*(omega_lc/omega_e_eq)**6+beta6/7.0*(omega_lc/omega_e_eq)**7+beta7/8.0*(omega_lc/omega_e_eq)**8
                erf2=beta0*(omega_uc/omega_e_eq)+beta1/2*(omega_uc/omega_e_eq)**2+beta2/3.0*(omega_uc/omega_e_eq)**3+beta3/4.0*(omega_uc/omega_e_eq)**4+beta4/5.0*(omega_uc/omega_e_eq)**5&
                        +beta5/6.0*(omega_uc/omega_e_eq)**6+beta6/7.0*(omega_uc/omega_e_eq)**7+beta7/8.0*(omega_uc/omega_e_eq)**8
                PSD=PSD/((erf2-erf1)*omega_e_eq)
                if(PSD.lt.1.0e-20)then
                        PSD=1.0e-20
                endif
        CASE('2018UBCHORUS')
        if(MLT.ge.6.and.MLT.lt.18)then!!Dayside only for lat lower than 10
                        beta0=2.0543522894746959e-4
                        beta1=-2.1251727655271302e-3
                        beta2=9.2309482352891624e-3
                        beta3=-2.1851729739887232e-2
                        beta4=3.0496655461723709e-2
                        beta5=-2.5133977878939284e-2
                        beta6=1.1343438413049637e-2
                        beta7=-2.1655977959280860e-3
        else !!!Nightside only for lat lower than 10 degree
                        beta0=-1.4010228888340406e-4
                        beta1=8.2669571248420001e-4
                        beta2=-1.3655153785780715e-3
                        beta3=-1.1351830489998828e-3
                        beta4=6.6938370901288491e-3
                        beta5=-9.0877873665073227e-3
                        beta6=5.4834375943244367e-3
                        beta7=-1.2753845539271731e-3
        endif
                PSD=beta0+beta1*(Freq/omega_e_eq)+beta2*(Freq/omega_e_eq)**2+beta3*(Freq/omega_e_eq)**3+beta4*(Freq/omega_e_eq)**4+beta5*(Freq/omega_e_eq)**5+beta6*(Freq/omega_e_eq)**6+beta7*(Freq/omega_e_eq)**7
                erf1=beta0*(omega_lc/omega_e_eq)+beta1/2*(omega_lc/omega_e_eq)**2+beta2/3.0*(omega_lc/omega_e_eq)**3+beta3/4.0*(omega_lc/omega_e_eq)**4+beta4/5.0*(omega_lc/omega_e_eq)**5+&
                        beta5/6.0*(omega_lc/omega_e_eq)**6+beta6/7.0*(omega_lc/omega_e_eq)**7+beta7/8.0*(omega_lc/omega_e_eq)**8
                erf2=beta0*(omega_uc/omega_e_eq)+beta1/2*(omega_uc/omega_e_eq)**2+beta2/3.0*(omega_uc/omega_e_eq)**3+beta3/4.0*(omega_uc/omega_e_eq)**4+beta4/5.0*(omega_uc/omega_e_eq)**5+&
                        beta5/6.0*(omega_uc/omega_e_eq)**6+beta6/7.0*(omega_uc/omega_e_eq)**7+beta7/8.0*(omega_uc/omega_e_eq)**8
                PSD=PSD/((erf2-erf1)*omega_e_eq)
                if(PSD.lt.1.0e-20)then
                        PSD=1.0e-20
                endif        
        CASE('2018LBCHORUS')
        if(MLT.ge.6.and.MLT.lt.18)then!!Dayside
                if(lambda_deg.le.10)then!!Dayside lowlatitude
                        beta0=3.6057e-7
                        beta1=-1.3014e-5
                        beta2=2.1507e-4
                        beta3=-1.4922e-3
                        beta4=5.5155e-3
                        beta5=-1.1273e-2
                        beta6=1.1739e-2
                        beta7=-4.7206e-3
                else
                        beta0=4.4083e-7
                        beta1=-1.4193e-5
                        beta2=2.0794e-4
                        beta3=-1.0750e-3
                        beta4=2.4121e-3
                        beta5=-1.6367e-3
                        beta6=-2.1965e-3
                        beta7=3.0513e-3
                endif
        else  !!!Nightside
                if(lambda_deg.le.10)then!!Nightside lowlatitude
                        beta0=-5.3446e-7
                        beta1=2.4890e-5
                        beta2=-4.2571e-4
                        beta3=4.1014e-3
                        beta4=-2.0035e-2
                        beta5=5.2005e-2
                        beta6=-6.9501e-2
                        beta7=3.7751e-2
                else
                        beta0=-1.2997e-7
                        beta1=9.0475e-6
                        beta2=-1.8619e-4
                        beta3=1.9125e-3
                        beta4=-9.8212e-3
                        beta5=2.6168e-2
                        beta6=-3.4857e-2
                        beta7=1.8411e-2
                endif
        endif
                PSD=beta0+beta1*(Freq/omega_e_eq)+beta2*(Freq/omega_e_eq)**2+beta3*(Freq/omega_e_eq)**3+beta4*(Freq/omega_e_eq)**4+beta5*(Freq/omega_e_eq)**5+beta6*(Freq/omega_e_eq)**6+beta7*(Freq/omega_e_eq)**7
                erf1=beta0*(omega_lc/omega_e_eq)+beta1/2*(omega_lc/omega_e_eq)**2+beta2/3.0*(omega_lc/omega_e_eq)**3+beta3/4.0*(omega_lc/omega_e_eq)**4+beta4/5.0*(omega_lc/omega_e_eq)**5&
                        +beta5/6.0*(omega_lc/omega_e_eq)**6+beta6/7.0*(omega_lc/omega_e_eq)**7+beta7/8.0*(omega_lc/omega_e_eq)**8
                erf2=beta0*(omega_uc/omega_e_eq)+beta1/2*(omega_uc/omega_e_eq)**2+beta2/3.0*(omega_uc/omega_e_eq)**3+beta3/4.0*(omega_uc/omega_e_eq)**4+beta4/5.0*(omega_uc/omega_e_eq)**5&
                        +beta5/6.0*(omega_uc/omega_e_eq)**6+beta6/7.0*(omega_uc/omega_e_eq)**7+beta7/8.0*(omega_uc/omega_e_eq)**8
                PSD=PSD/((erf2-erf1)*omega_e_eq)
                if(PSD.lt.1.0e-20)then
                        PSD=1.0e-20
                endif
        CASE DEFAULT
                        WRITE(*,*) 'Invalid Frequency Distribution'
                        STOP
        END SELECT




end subroutine PowerSpectralDensity
