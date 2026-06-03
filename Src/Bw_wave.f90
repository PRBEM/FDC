!This subroutine is to caculate the wave amplitude of waves in different models
!The case 'OA2015' is refered to Oleksiy Agapitov's 2015 JGR paper

subroutine Bw_wave(L,lambda,Bw)
        use constants

        implicit none

	real*8,intent(in)::L,lambda
	real*8,intent(out)::Bw
	real*8::lambda_deg
	real*8::Ltemp
	real*8::log10_Bwsquared
        integer::i,j
        real*8::coeff(4,4)
	real*8::a1
	real*8::a2
	real*8::a3
	real*8::a4
	real*8::a5
	real*8::a6
	real*8::a7
	real*8::a8
	real*8::a9
        real*8::g_Kp, G0, b_0, b_1, b_2, b_3, b_4
        real*8::a_0, a_1, a_2
        real*8::aaa,bbb
	real*8::aa(4,4),bb(4),Bw_Input1,Bw_Input2
        integer::ii,jj
  
        lambda_deg=lambda/RAD


        SELECT CASE(BwType)
                CASE('TYPICALDAYSIDECHORUS') !This case is a object for comparison in Hui's 2017 paper
                        if(Kp.le.2)then
                                Bw_Input=10**(0.75+0.04*lambda_deg)*sqrt(2.0*10**(0.73+0.91*Kp))/57.6
                        else 
                                Bw_Input=10**(0.75+0.04*lambda_deg)*sqrt(2.0*10**(2.5+0.18*Kp))/57.6
                        endif
                CASE('TYPICALNIGHTSIDECHORUS') !This case is a object for comparison in Hui's 2017 paper
                        if(Kp.le.2)then
                                Bw_Input=50.0*sqrt(2.0*10**(0.73+0.91*Kp))/57.6
                        else 
                                Bw_Input=50.0*sqrt(2.0*10**(2.5+0.18*Kp))/57.6
                        endif
                CASE('TYPICALHISS') 
                        Bw_Input=80.0*Kp/4.0
                CASE('KS14DAYSIDEHISS') 
                        a1=2.5593
                        a2=0.26074
                        a3=0.55554
                        a4=-0.1263
                        a5=0.018152
                        a6=-0.015465
                        a7=-0.14001
                        a8=0.0030924
                        if(lambda_deg>25)then
                                lambda_deg=25
                        end if
                        log10_Bwsquared= a1 + a2*Kp+ a3*L+ a4*lambda_deg + a5*L*lambda_deg + a6*Kp*Kp + a7*L*L + a8*lambda_deg*lambda_deg
                        Bw_Input=sqrt(10**log10_Bwsquared)
                CASE('KS14NIGHTSIDEHISS') 

                        a1=2.4319
                        a2=-0.088243
                        a3=0.63916
                        a4=-0.13609
                        a5=-0.012665
                        a6=0.03026
                        a7=0.042481
                        a8=-0.15691
                        a9=0.0014566

                        if(lambda_deg>25)then
                                lambda_deg=25
                        end if
                        log10_Bwsquared= a1 + a2*Kp+ a3*L+ a4*lambda_deg + a5*Kp*lambda_deg + a6*L*lambda_deg + a7*Kp*Kp + a8*L*L + a9*lambda_deg*lambda_deg

                        Bw_Input=sqrt(10**log10_Bwsquared)
                CASE('2016DAYSIDECHORUS') 
                        if(lambda_deg>20.0)then
                                lambda_deg=20.0
                        end if
                        log10_Bwsquared=1.0561 + 0.4095*Kp + 0.0351*lambda_deg
                        Bw_Input=sqrt(10**log10_Bwsquared)
                CASE('2016NIGHTSIDECHORUS') 
                        log10_Bwsquared=1.8821 + 0.3637*Kp - 0.05015*lambda_deg
                        Bw_Input=sqrt(10**log10_Bwsquared)
                CASE('2016NIGHTSIDECHORUS2') 
                        log10_Bwsquared=1.0662 + 0.6677*Kp + 0.0223*lambda_deg - 0.0312*Kp*lambda_deg
                        Bw_Input=sqrt(10**log10_Bwsquared)
                CASE('2016HISS')
                        log10_Bwsquared=0.05894*L**3-0.7768*L**2+3.163*L-1.036
                        Bw_Input=sqrt(10**log10_Bwsquared)
                CASE('UPPERBANDNIGHTSIDECHORUS')
                        log10_Bwsquared=1.6957 + 0.17491*Kp - 0.15828*lambda_deg
                        Bw_Input=sqrt(10**log10_Bwsquared)
                CASE('UPPERBANDDAYSIDECHORUS')
                        log10_Bwsquared=1.0769 + 0.23987*Kp - 0.067579*lambda_deg
                        Bw_Input=sqrt(10**log10_Bwsquared)
                CASE('NBBLBC')
                        if(kp.eq.0)then
                                if(lambda_deg.le.5)then
                                        Bw_Input=10.0
                                else if(lambda_deg.le.10)then
                                        Bw_Input=10.0
                                else
                                        Bw_Input=10.0
                                end if
                                        
                        else if(kp.eq.1)then
                                if(lambda_deg.le.5)then
                                        Bw_Input=20.0
                                else if(lambda_deg.le.10) then
                                        Bw_Input=35.0
                                else
                                        Bw_Input=20.0
                                end if
                        else
                                if(lambda_deg.le.5)then
                                        Bw_Input=65.0
                                else if(lambda_deg.le.10) then
                                        Bw_Input=80.0
                                else
                                        Bw_Input=60.0
                                end if
                        end if
                 CASE('NBBUBC')
                        if(kp.eq.0)then
                                if(lambda_deg.le.5)then
                                        Bw_Input=5.0
                                else 
                                        Bw_Input=5.0
                                end if
                                        
                        else if(kp.eq.1)then
                                if(lambda_deg.le.5)then
                                        Bw_Input=15.0
                                else 
                                        Bw_Input=10.0
                                end if
                        else 
                                if(lambda_deg.le.5)then
                                        Bw_Input=45.0
                                else
                                        Bw_Input=20.0
                                end if
                        end if
                

                CASE('CONSTANT') 

                CASE('OA2015')
                   
                       if(MLT.ge.3.and.MLT.Lt.15)then !dayside
                           if(L.le.5)then
                                coeff(1,1)=-2.53242
                                coeff(1,2)=0.14836
                                coeff(1,3)=-0.01648
                                coeff(1,4)=0.00016
                                coeff(2,1)=0.11290
                                coeff(2,2)=-0.47100
                                coeff(2,3)=0.23260
                                coeff(2,4)=-0.02317
                                coeff(3,1)=-0.08181
                                coeff(3,2)=0.28823
                                coeff(3,3)=-0.15163
                                coeff(3,4)=0.01726
                                coeff(4,1)=0.01027
                                coeff(4,2)=-0.04109
                                coeff(4,3)=0.02380
                                coeff(4,4)=-0.00311
                           else
                                coeff(1,1)=-2.40708
                                coeff(1,2)=-0.02297
                                coeff(1,3)=0.01223
                                coeff(1,4)=-0.00021
                                coeff(2,1)=-0.29267
                                coeff(2,2)=0.45874
                                coeff(2,3)=-0.10822
                                coeff(2,4)=0.00812
                                coeff(3,1)=0.18808
                                coeff(3,2)=-0.22663
                                coeff(3,3)=0.05810
                                coeff(3,4)=-0.00446
                                coeff(4,1)=-0.02839
                                coeff(4,2)=0.03140
                                coeff(4,3)=-0.00859
                                coeff(4,4)=0.00066
                           end if
                       else
                          if(L.le.5)then
                                coeff(1,1)=-2.39346
                                coeff(1,2)=0.10988
                                coeff(1,3)=-0.03083
                                coeff(1,4)=0.00431
                                coeff(2,1)=-0.03083
                                coeff(2,2)=0.16417
                                coeff(2,3)=-0.01875
                                coeff(2,4)=-0.00346
                                coeff(3,1)=-0.00884
                                coeff(3,2)=-0.13052
                                coeff(3,3)=0.00501
                                coeff(3,4)=0.00680
                                coeff(4,1)=0.00334
                                coeff(4,2)=0.00089
                                coeff(4,3)=0.01714
                                coeff(4,4)=-0.00523
                          else
                                coeff(1,1)=-2.69548
                                coeff(1,2)=0.22568
                                coeff(1,3)=-0.04232
                                coeff(1,4)=0.00423
                                coeff(2,1)=-0.00652
                                coeff(2,2)=-0.09346
                                coeff(2,3)=0.05615
                                coeff(2,4)=-0.00424
                                coeff(3,1)=0.07774
                                coeff(3,2)=0.12785
                                coeff(3,3)=-0.06845
                                coeff(3,4)=0.00524
                                coeff(4,1)=-0.01756
                                coeff(4,2)=-0.02291
                                coeff(4,3)=0.01166
                                coeff(4,4)=-0.00080
                          end if
                      end if
                   Bw_input1=0.0
                   do i=0,3
                      do j=0,3
                        Bw_input1=Bw_input1+coeff(i+1,j+1)*Kp**j*(lambda_deg/10)**i
                      end do
                   end do
                   Bw_input2=10**Bw_input1
                   Bw_input=1000*Bw_input2

                CASE('VITALIICHORUS')
                        if(lambda_deg>20.0)then
                                lambda_deg=20.0
                        end if
                        if(MLT.ge.3.and.MLT.Lt.15)then
                                aa(1,1)=-5.387e-3
                                aa(1,2)=7.326e-3
                                aa(1,3)=-2.976e-3
                                aa(1,4)=4.128e-4
                                aa(2,1)=8.242e-2
                                aa(2,2)=-8.72e-2
                                aa(2,3)=2.861e-2
                                aa(2,4)=-2.777e-3
                                aa(3,1)=-1.302e-1
                                aa(3,2)=1.309e-1
                                aa(3,3)=-4.07e-2
                                aa(3,4)=3.759e-3
                                aa(4,1)=5.032e-2
                                aa(4,2)=-4.968e-2
                                aa(4,3)=1.519e-2
                                aa(4,4)=-1.389e-3

                                bb(1)=8.408e-1
                                bb(2)=-2.484e-1
                                bb(3)=5.387e-1
                                bb(4)=-6.574e-2
                        else
                                aa(1,1)=8.671e-3
                                aa(1,2)=-8.872e-3
                                aa(1,3)=2.848e-3
                                aa(1,4)=-2.689e-4
                                aa(2,1)=3.613e-3
                                aa(2,2)=-5.749e-3
                                aa(2,3)=2.599e-3
                                aa(2,4)=-3.138e-4
                                aa(3,1)=-1.691e-2
                                aa(3,2)=2.035e-2
                                aa(3,3)=-7.574e-3
                                aa(3,4)=8.172e-4
                                aa(4,1)=6.087e-3
                                aa(4,2)=-7.417e-3
                                aa(4,3)=2.781e-3
                                aa(4,4)=-3.013e-4

                                bb(1)=3.684
                                bb(2)=-2.46
                                bb(3)=2.828
                                bb(4)=-3.018e-1
                        end if
                        Bw_Input1=0.0
                        do ii=1,4
                                do jj=1,4
                                        Bw_Input1=Bw_Input1+aa(ii,jj)*(lambda_deg/10.0)**(ii-1)*L**(jj-1)
                                enddo
                        enddo
                        Bw_Input2=0.0
                        do ii=1,4
                                Bw_Input2=Bw_Input2+bb(ii)*Kp**(ii-1)
                        enddo
                        Bw_Input=Bw_Input1*Bw_Input2*1e3 !nT to pT
                CASE('2017CHORUS')
                        if(lambda_deg.ge.20.0)then
                                lambda_deg=20.0
                        endif
                        if(L.lt.3.5)then
                                ltemp=3.5
                        else if(L.gt.6)then
                                ltemp=6
                        else
                                ltemp=L
                        end if
                        aaa=cos(MLT*pi/12.0)
                        bbb=sin(MLT*pi/12.0)
                        Bw_Input=-3266.7 + 546.23*aaa - 670.18*bbb + 1150*Ltemp + 87.776*lambda_deg + 199.54*aaa*bbb +&
                        228.42*bbb*Ltemp - 77.367*Ltemp**2 - 54.412*aaa*lambda_deg - 10.025*bbb*lambda_deg - 18.864*Ltemp*lambda_deg
                        Bw_Input=Bw_Input*(0 + 75.532*Kp - 45.677*Kp**2 + 50.591*Kp**3 - 6.265*Kp**4 + 0.53362*Kp**5)/1478.1
                        if(Bw_Input.lt.1.0e-2)then
                                Bw_Input=1.0e-2
                        endif
                        Bw_Input=sqrt(Bw_Input)
                CASE('2017UPPERBANDCHORUS')
                      if(lambda_deg.ge.17.0)then
                               lambda_deg=17.0
                      endif
                      if(L.lt.3.5)then
                              ltemp=3.5
                      else if(L.gt.6)then
                              ltemp=6
                      else
                              ltemp=L
                      endif
                        aaa=cos(MLT*pi/12.0)
                        bbb=sin(MLT*pi/12.0)
                        Bw_Input=510.0363206 + 206.0115305*aaa + 294.2069939*bbb - 64.57977972*Ltemp - 48.76766076*lambda_deg + 27.89711269*aaa*bbb - 26.18987181*aaa*Ltemp& 
                                - 5.740739266*aaa*lambda_deg - 27.68773262*bbb*Ltemp - 10.6820023*bbb*Ltemp + 0.2699685419*Ltemp**2 + 4.326157077*Ltemp*lambda_deg + 0.9335134531*lambda_deg**2
                        Bw_Input=Bw_Input*(0 - 9.752027449*Kp + 26.46676219*Kp**2 - 5.243283335*Kp**3 + 0.3785342851*Kp**4)/103.3238506
                       if(Bw_Input.lt.1.0e-2)then
                          Bw_Input=1.0e-2
                       end if
                        Bw_Input=sqrt(Bw_Input)
                CASE('CHORUSLAT1545')
                        if(L.lt.3.5)then
                                ltemp=3.5
                        else if(L.gt.6)then
                                ltemp=6
                        else
                                ltemp=L
                        end if
                        if(ltemp.le.5)then
                                if(MLT.ge.13.and.MLT.Lt.13)then
                                        a_0 = 0.0104
                                        a_1 = -0.8551
                                        a_2 = 21.2939
                                        b_0 = -0.0192
                                        b_1 = 0.1750
                                        b_2 = -0.4159
                                        b_3 = 0.4302
                                        b_4 = 0.2282
                                else if(MLT.ge.13.and.MLT.Lt.22)then
                                        a_0 = 0.0104
                                        a_1 = -0.8551
                                        a_2 = 21.2939
                                        b_0 = -0.0034
                                        b_1 = -0.0136
                                        b_2 = 0.2681
                                        b_3 = -0.3412
                                        b_4 = 0.4098
                                else
                                        a_0 = -0.0070
                                        a_1 = 0.2295
                                        a_2 = 3.9511
                                        b_0 = 0.0241
                                        b_1 = -0.1892
                                        b_2 = 0.3553
                                        b_3 = 0.0087
                                        b_4 = 0.6821
                                end if
                        else!high L
                                if(MLT.ge.4.and.MLT.Lt.14)then
                                        a_0 = -0.0173
                                        a_1 = 0.8612
                                        a_2 = 4.1508
                                        b_0 = -0.0082
                                        b_1 = 0.0941
                                        b_2 = -0.3079
                                        b_3 = 0.4680
                                        b_4 = 0.3394
                                else if(MLT.ge.14.and.MLT.Lt.21)then
                                        a_0 = -0.0084
                                        a_1 = 0.3892
                                        a_2 = 6.1371
                                        b_0 = -0.0036
                                        b_1 = 0.0327
                                        b_2 = -0.1312
                                        b_3 = 0.2713
                                        b_4 = 0.9189
                                else
                                        a_0 = -0.0058
                                        a_1 = 0.4085
                                        a_2 = -2.6257
                                        b_0 = -0.0160
                                        b_1 = 0.1826
                                        b_2 = -0.6751
                                        b_3 = 0.8386
                                        b_4 = 0.8622
                                end if
                        end if
                        g_Kp=b_0*Kp**4+b_1*Kp**3+b_2*Kp**2+b_3*Kp+b_4
                        Bw_Input=a_0*lambda_deg**2+a_1*lambda_deg+a_2
                        Bw_Input=Bw_Input*g_Kp
                        if(Bw_Input.lt.1.0e-2)then
                                Bw_Input=1.0e-2
                        end if


                CASE('MIDLATCHORUS')
                        if(lambda_deg.ge.45)then
                                lambda_deg=45.0
                        end if
                        if(L.lt.3.5)then
                                ltemp=3.5
                        else if(L.gt.6)then
                                ltemp=6
                        else
                                ltemp=L
                        end if
                        if(ltemp.le.5)then
                                if(MLT.ge.4.and.MLT.Lt.13)then
                                        a_0=2.5452653e-2
                                        a_1=-1.8204229
                                        a_2=3.5917442e1
                                        b_0=-2.1785381e-2
                                        b_1=2.1934154e-1
                                        b_2=-6.5135800e-1
                                        b_3=8.1637208e-1
                                        b_4=1.6950750e-1
                                else if(MLT.ge.13.and.MLT.Lt.22)then
                                        a_0=5.8055846e-3
                                        a_1=-5.5876391e-1
                                        a_2=1.5572106e1
                                        b_0=-1.7678770e-2
                                        b_1=1.6252408e-1
                                        b_2=-3.9468309e-1
                                        b_3=4.0186281e-1
                                        b_4=3.3664550e-1
                                else
                                        a_0=3.0811025e-3
                                        a_1=-4.0999943e-1
                                        a_2=1.3541850e1
                                        b_0=1.0852508e-2
                                        b_1=-2.0289397e-2
                                        b_2=-3.3138473e-1
                                        b_3=9.3497384e-1
                                        b_4=4.5435250e-1
                                end if
                        else
                                if(MLT.ge.4.and.MLT.Lt.14)then
                                        a_0=-2.5864386e-2
                                        a_1=1.4217082
                                        a_2=-4.7289716
                                        b_0=-1.2138862e-2
                                        b_1=1.3788319e-1
                                        b_2=-4.4842546e-1
                                        b_3=5.8112700e-1
                                        b_4=3.3842536e-1
                                else if(MLT.ge.14.and.MLT.Lt.21)then
                                        a_0=-1.2138022e-2
                                        a_1=6.4269273e-1
                                        a_2=2.0069819
                                        b_0=-2.4564293e-3
                                        b_1=2.2152503e-2
                                        b_2=-9.7797716e-2
                                        b_3=2.0362217e-1
                                        b_4=1.0031774
                                else
                                        a_0=-1.1856737e-2
                                        a_1=8.1229858e-1
                                        a_2=-9.2016650
                                        b_0=-1.1385737e-2
                                        b_1=1.3457093e-1
                                        b_2=-5.3830828e-1
                                        b_3=7.5101515e-1
                                        b_4=8.9100649e-1
                                end if
                        end if
                        g_Kp=b_0*Kp**4+b_1*Kp**3+b_2*Kp**2+b_3*Kp+b_4
                        Bw_Input=a_0*lambda_deg**2+a_1*lambda_deg+a_2
                        Bw_Input=Bw_Input*g_Kp
                        if(Bw_Input.lt.1.0e-2)then
                                Bw_Input=1.0e-2
                        end if
                CASE('2018UBCHORUS')
                        if(L.lt.3.5)then
                                ltemp=3.5
                        else if(L.gt.6)then
                                ltemp=6
                        else
                                ltemp=L
                        end if
                        if(ltemp.le.5)then
                                if(MLT.ge.6.and.MLT.Lt.12)then
                                        G0 =8.218
                                        b_0=0
                                        b_1=0.189
                                        b_2=0.301
                                        b_3=0.348
                                        b_4=-0.057
                                else if(MLT.ge.12.and.MLT.Lt.22)then
                                        G0=1.392
                                        b_0=0
                                        b_1=0.011
                                        b_2=0.012
                                        b_3=-0.003
                                        b_4=0.005
                                else
                                        G0 =11.633
                                        b_0=0
                                        b_1=0.070
                                        b_2=0.202
                                        b_3=0.335
                                        b_4=-0.035
                                endif
                        else
                                if(MLT.ge.7.and.MLT.Lt.13)then
                                        G0 =5.878
                                        b_0=0
                                        b_1=0.771
                                        b_2=0.901
                                        b_3=-0.048
                                        b_4=-0.018
                                else if(MLT.ge.13.and.MLT.Lt.21)then
                                        G0 =2.099
                                        b_0=0
                                        b_1=0.018
                                        b_2=0.066
                                        b_3=0.135
                                        b_4=-0.023
                                else
                                        G0 =9.754
                                        b_0=0
                                        b_1=1.469
                                        b_2=2.086
                                        b_3=-0.534
                                        b_4=0.035
                                end if
                        end if
                        g_Kp=(b_0+b_1*Kp+b_2*Kp**2+b_3*Kp**3+b_4*Kp**4)/G0
                        if(g_Kp.lt.1.0e-2)then
                                g_Kp=1.0e-2
                        endif
                        aaa=cos(MLT*pi/12.0)
                        bbb=sin(MLT*pi/12.0)
                        Bw_Input=-0.896+23.216*aaa+15.102*bbb+6.580*ltemp+0.155*lambda_deg+0.797*aaa*bbb-2.943*aaa*ltemp-1.148*aaa*lambda_deg-1.225*bbb*ltemp-0.796*bbb*lambda_deg&
                                -0.808*ltemp**2-0.057*ltemp*lambda_deg-0.092*lambda_deg**2 
                        if(Bw_Input.lt.1.0e-2)then
                                Bw_Input=1.0e-2
                        endif
                        Bw_Input=Bw_Input*g_Kp

                CASE('2018LBCHORUS')
                        if(lambda_deg.ge.20.0)then
                                lambda_deg=20.0
                        endif
                        if(L.lt.3.5)then
                                ltemp=3.5
                        else if(L.gt.6)then
                                ltemp=6
                        else
                                ltemp=L
                        end if
                        if(ltemp.le.5)then
                                if(MLT.ge.4.and.MLT.Lt.13)then
                                        G0=18.308
                                        b_0=0
                                        b_1=0.064
                                        b_2=0.179
                                        b_3=0.297
                                        b_4=-0.001
                                else if(MLT.ge.13.and.MLT.Lt.22)then
                                        G0=4.101
                                        b_0=0
                                        b_1=0.036
                                        b_2=0.097
                                        b_3=0.167
                                        b_4=-0.024
                                else
                                        G0=23.141
                                        b_0=0
                                        b_1=0.083
                                        b_2=0.278
                                        b_3=0.556
                                        b_4=-0.040
                                endif
                        else
                                if(MLT.ge.4.and.MLT.Lt.14)then
                                        G0=27.622
                                        b_0=0
                                        b_1=0.194
                                        b_2=0.456
                                        b_3=0.724
                                        b_4=-0.068
                                else if(MLT.ge.14.and.MLT.Lt.21)then
                                        G0=13.783
                                        b_0=0
                                        b_1=0.241
                                        b_2=0.471
                                        b_3=0.710
                                        b_4=-0.119
                                else
                                        G0=30.181
                                        b_0=0
                                        b_1=1.461
                                        b_2=2.632
                                        b_3=0.184
                                        b_4=-0.061
                                end if
                        end if
                        g_Kp=(b_0+b_1*Kp+b_2*Kp**2+b_3*Kp**3+b_4*Kp**4)/G0
                        if(g_Kp.lt.1.0e-2)then
                                g_Kp=1.0e-2
                        endif
                        aaa=cos(MLT*pi/12.0)
                        bbb=sin(MLT*pi/12.0)
                        Bw_Input=-136.706+35.308*aaa+4.571*bbb+49.943*ltemp+4.626*lambda_deg+2.410*aaa*bbb-3.871*aaa*ltemp-1.780*aaa*lambda_deg+1.739*bbb*ltemp-0.462*bbb*lambda_deg&
                                -3.801*ltemp**2-0.848*ltemp*lambda_deg-0.017*lambda_deg**2
                        if(Bw_Input.lt.1.0e-2)then
                                Bw_Input=1.0e-2
                        endif
                        Bw_Input=Bw_Input*g_Kp
                        
                        
                CASE('2026LBCHORUS')
                        if(lambda_deg.gt.40.0)then
                                lambda_deg=40.0
                        end if
                        if(L.lt.3.5)then
                                ltemp=3.5
                        else if(L.gt.6.5)then
                                ltemp=6.5
                        else
                                ltemp=L
                        end if
                        G0 = 29.169
                        b_0=8.847
                        b_1=1.896
                        b_2=1.834
                        b_3=-0.247
                        b_4=0.014
                        g_Kp=(b_0+b_1*Kp+b_2*Kp**2+b_3*Kp**3+b_4*Kp**4)/G0
                        if(g_Kp.lt.1.0e-2)then
                                g_Kp=1.0e-2
                        endif
                        aaa=cos(MLT*pi/12.0)
                        bbb=sin(MLT*pi/12.0)
                        Bw_Input=-79.283+33.851*aaa+10.731*bbb+31.739*ltemp+1.985*lambda_deg+1.258*aaa*bbb-4.564*aaa*ltemp-1.253*aaa*lambda_deg+0.48*bbb*ltemp-0.549*bbb*lambda_deg&
                                -2.497*ltemp**2-0.405*ltemp*lambda_deg+0.013*lambda_deg**2
                        if(Bw_Input.lt.1.0e-2)then
                                Bw_Input=1.0e-2
                        endif
                        Bw_Input=Bw_Input*g_Kp

                CASE('2026UBCHORUS')
                        if(lambda_deg.gt.10.0)then
                                lambda_deg=10.0
                        end if
                        if(L.lt.4)then
                                ltemp=4
                        else if(L.gt.6.5)then
                                ltemp=6.5
                        else
                                ltemp=L
                        end if
                        G0 = 12.449
                        b_0=2.283
                        b_1=1.125
                        b_2=1.429
                        b_3=-0.319
                        b_4=0.023
                    
                        g_Kp=(b_0+b_1*Kp+b_2*Kp**2+b_3*Kp**3+b_4*Kp**4)/G0
                        if(g_Kp.lt.1.0e-2)then
                                g_Kp=1.0e-2
                        endif
                        aaa=cos(MLT*pi/12.0)
                        bbb=sin(MLT*pi/12.0)
                        Bw_Input=18.59+20.906*aaa+15.359*bbb-0.408*ltemp-1.055*lambda_deg+0.843*aaa*bbb-2.64*aaa*ltemp-1.043*aaa*lambda_deg-1.353*bbb*ltemp-0.749*bbb*lambda_deg&
                                -0.202*ltemp**2+0.13*ltemp*lambda_deg-0.062*lambda_deg**2 
                        if(Bw_Input.lt.1.0e-2)then
                                Bw_Input=1.0e-2
                        endif
                        Bw_Input=Bw_Input*g_Kp



                CASE DEFAULT
                        WRITE(*,*) 'Invalid Wave Amplitude.'
                        STOP
                END SELECT

        Bw=Bw_Input*1.0e-12  !pT to T


        return

end subroutine Bw_wave
