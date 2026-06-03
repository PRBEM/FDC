!Frequency range of wave
subroutine Freq_Spectral(L,Lambda,omega_e_eq,omega_m,d_omega,omega_lc,omega_uc)
        use constants
        implicit none 
	real*8,intent(IN)::L,Lambda,omega_e_eq
	real*8,intent(out)::omega_m,d_omega,omega_lc,omega_uc
	real*8::lambda_deg

        lambda_deg=lambda/rad
        SELECT CASE(FreqUnit)
                CASE('HZ') 
                        d_omega=d_omega_per*2.0*PI
                        omega_m=omega_m_per*2.0*PI
                        omega_uc=o_uc*2.0*PI
                        omega_lc=o_lc*2.0*PI
                CASE('RATIO') 
                        d_omega=d_omega_per*omega_e_eq
                        omega_m=omega_m_per*omega_e_eq
                        omega_uc=o_uc*omega_e_eq
                        omega_lc=o_lc*omega_e_eq
                CASE('2016HISS')
                        if(L.le.2.5)then
                                d_omega=290.0*2.0*PI
                                omega_m=282*2.0*PI
                                omega_uc=4000.0*2.0*PI
                                omega_lc=20*2.0*PI
                        else if(L.gt.2.5.and.L.le.3.5)then
                                d_omega=160.0*2.0*PI
                                omega_m=252*2.0*PI
                                omega_uc=4000.0*2.0*PI
                                omega_lc=20*2.0*PI
                        else if(L.gt.3.5.and.L.le.4.5)then
                                d_omega=110.0*2.0*PI
                                omega_m=178*2.0*PI
                                omega_uc=4000.0*2.0*PI
                                omega_lc=20*2.0*PI
                        else
                                d_omega=90.0*2.0*PI
                                omega_m=140*2.0*PI
                                omega_uc=4000.0*2.0*PI
                                omega_lc=20*2.0*PI

                        endif
                CASE('NBBLBC')
                        if(kp.eq.0)then
                                if(lambda_deg.le.5)then
                                        omega_m=0.28*omega_e_eq
                                        d_omega=0.07*omega_e_eq
                                else if(lambda_deg.le.10)then
                                        omega_m=0.25*omega_e_eq
                                        d_omega=0.13*omega_e_eq
                                else
                                        omega_m=0.22*omega_e_eq
                                        d_omega=0.1*omega_e_eq
                                endif
                        else if(kp.eq.1)then
                                if(lambda_deg.le.5)then
                                        omega_m=0.28*omega_e_eq
                                        d_omega=0.12*omega_e_eq
                                else if(lambda_deg.le.10)then
                                        omega_m=0.3*omega_e_eq
                                        d_omega=0.08*omega_e_eq
                                else
                                        omega_m=0.26*omega_e_eq
                                        d_omega=0.07*omega_e_eq
                                end if
                        else
                                if(lambda_deg.le.5)then
                                        omega_m=0.35*omega_e_eq
                                        d_omega=0.15*omega_e_eq
                                else if(lambda_deg.le.10)then
                                        omega_m=0.3*omega_e_eq
                                        d_omega=0.15*omega_e_eq
                                else
                                        omega_m=0.25*omega_e_eq
                                        d_omega=0.1*omega_e_eq
                                end if

                        end if
                        omega_lc=0.1*omega_e_eq
                        omega_uc=0.5*omega_e_eq
                CASE('NBBUBC')
                        if(kp.eq.0)then
                                if(lambda_deg.le.5)then
                                        omega_m=0.54*omega_e_eq
                                        d_omega=0.08*omega_e_eq
                                else
                                        omega_m=0.56*omega_e_eq
                                        d_omega=0.06*omega_e_eq
                                endif
                        else if(kp.eq.1)then
                                if(lambda_deg.le.5)then
                                        omega_m=0.54*omega_e_eq
                                        d_omega=0.1*omega_e_eq
                                else
                                        omega_m=0.57*omega_e_eq
                                        d_omega=0.04*omega_e_eq
                                endif
                        else
                                if(lambda_deg.le.5)then
                                        omega_m=0.53*omega_e_eq
                                        d_omega=0.1*omega_e_eq
                                else
                                        omega_m=0.58*omega_e_eq
                                        d_omega=0.04*omega_e_eq
                                endif

                        end if
                        omega_lc=0.5*omega_e_eq
                        omega_uc=0.65*omega_e_eq
                CASE('2017CHORUS')
                        if(MLT.ge.6.and.MLT.lt.18)then
                                omega_lc=0.05*omega_e_eq
                                omega_uc=0.5*omega_e_eq
                        else
                                omega_lc=0.05*omega_e_eq
                                omega_uc=0.5*omega_e_eq
                        endif
                CASE('18LBCHORUS')
                        if(MLT.ge.6.and.MLT.lt.18)then
                                omega_lc=0.05*omega_e_eq
                                omega_uc=0.5*omega_e_eq
                        else
                                omega_lc=0.05*omega_e_eq
                                omega_uc=0.5*omega_e_eq
                        endif
                CASE('18UBCHORUS')
                        if(MLT.ge.6.and.MLT.lt.18)then
                                omega_lc=0.5*omega_e_eq
                                omega_uc=1.0*omega_e_eq
                        else
                                omega_lc=0.5*omega_e_eq
                                omega_uc=1.0*omega_e_eq
                        endif
                CASE('2026LBCHORUS')
                        omega_lc=0.05*omega_e_eq
                        omega_uc=0.5*omega_e_eq
                CASE('2026UBCHORUS')
                        omega_lc=0.5*omega_e_eq
                        omega_uc=0.8*omega_e_eq
                CASE DEFAULT
                        WRITE(*,*) 'Invalid Frequency Unit.'
                        STOP
        END SELECT
return
end subroutine Freq_Spectral

