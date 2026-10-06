! ==============================================================================
! SUBROUTINE: Bw_wave
! ==============================================================================
!
! PURPOSE:
!   Returns the RMS wave magnetic field amplitude Bw [T] at a given L-shell
!   and magnetic latitude, for the wave amplitude model specified by BwType
!   in module constants.
!
!   All models compute the amplitude in pT internally. The final conversion
!   to Tesla is applied at the end: Bw = Bw_pT * 1e-12.
!
!   Supported models and their parameter dependencies:
!
!     Model                   L    lambda   MLT   Kp   Description
!     TYPICALDAYSIDECHORUS          X             X    Hui 2017 comparison
!     TYPICALNIGHTSIDECHORUS                      X    Hui 2017 comparison
!     TYPICALHISS                                 X    Simple Kp scaling
!     KS14DAYSIDEHISS         X      X       X         ?
!     KS14NIGHTSIDEHISS       X      X       X         ?
!     2016DAYSIDECHORUS              X       X         ?
!     2016NIGHTSIDECHORUS            X       X         ?
!     2016NIGHTSIDECHORUS2           X       X         ?
!     2016HISS                X                        L-dependent hiss amplitude, probably Orlova et al. (2016)
!     OA2015                  X      X       X         probably Agapitov et al. (2015)
!     VITALIICHORUS           X      X       X         Polynomial model, ?
!     2017CHORUS              X      X       X         ?
!     2017UPPERBANDCHORUS     X      X       X         ?
!     2018LBCHORUS            X      X       X         2018 lower band, Dedong Wang
!     2018UBCHORUS            X      X       X         2018 upper band, Dedong Wang
!     2026LBCHORUS            X      X       X         2026 lower band, Alwin Roy
!     2026UBCHORUS            X      X       X         2026 upper band, Alwin Roy
!     CHORUSLAT1545           X      X       X         Latitude 15-45 deg model, ?? -> has a bug
!     MIDLATCHORUS            X      X       X         Mid-latitude chorus, ?
!     NBBLBC                         X       X         ?
!     NBBUBC                         X       X         ?
!     CONSTANT                                         Fixed amplitude from input
!
! INPUTS:
!   L        - Magnetic L-shell [Earth radii, dimensionless]
!   lambda   - Magnetic latitude [radians]
!
! OUTPUTS:
!   Bw       - RMS wave magnetic field amplitude [T]
!
! MODULE DEPENDENCIES:
!   constants - BwType, Bw_Input, Kp, MLT, PI, RAD
!
! NOTES:
!   - All models compute an amplitude in pT in the local variable Bw_pT;
!     the conversion Bw = Bw_pT * 1e-12 is applied once at the end. The
!     module variable Bw_Input holds the value read from the input file and
!     is only read here (by the CONSTANT case), never overwritten.
!   - Several models clamp L to a valid range (e.g. [3.5, 6] or [4, 6.5])
!     using a local variable Ltemp. Similarly, lambda_deg is clamped to a
!     maximum value in some models. These clamps reflect the observational
!     data range used to fit the model.
!   - The CONSTANT case takes Bw_Input as set from the input file in module
!     constants; no model evaluation is performed.
!   - The local variables a1-a9 in the KS14 cases are regression
!     coefficients specific to those models. They are entirely unrelated to
!     the module variables a1-a6 in Module_Variables.f90.
!   - SUSPECTED TYPO: 2017UPPERBANDCHORUS contains the term bbb*Ltemp twice
!     (coefficients -27.68773262 and -10.6820023) and no bbb*lambda_deg term,
!     unlike the other models of this family. Left unchanged pending a check
!     against the original reference.
!   - BUG: The CHORUSLAT1545 case contains the condition
!       if(MLT.ge.13.and.MLT.Lt.13)
!     which can never be true (MLT cannot be simultaneously >= 13 and < 13).
!     The intended boundary is unknown — likely MLT.ge.3.and.MLT.Lt.13
!     by analogy with OA2015 and other models, but left unchanged pending
!     verification against the original reference.
!
! ==============================================================================

subroutine Bw_wave(L, lambda, Bw)
use constants
implicit none

    real*8, intent(in)  :: L, lambda
    real*8, intent(out) :: Bw

    real*8  :: lambda_deg, Ltemp
    real*8  :: log10_Bwsquared
    real*8  :: coeff(4,4)
    real*8  :: a1, a2, a3, a4, a5, a6, a7, a8, a9   ! KS14 hiss wave regression coefficients
    real*8  :: g_Kp, G0, b_0, b_1, b_2, b_3, b_4
    real*8  :: a_0, a_1, a_2
    real*8  :: aaa, bbb                              ! cos/sin of MLT angle
    real*8  :: aa(4,4), bb(4)
    real*8  :: Bw_Input1, Bw_Input2
    real*8  :: Bw_pT                                 ! wave amplitude [pT], local
    integer :: i, j, ii, jj

    lambda_deg = lambda / RAD

    select case(BwType)

        ! ----------------------------------------------------------------------
        ! TYPICALDAYSIDECHORUS: latitude and Kp dependent.
        ! Reference: Hui et al. [2017] comparison case.
        ! ----------------------------------------------------------------------
        case('TYPICALDAYSIDECHORUS')
            if (Kp <= 2.0) then
                Bw_pT = 10**(0.75 + 0.04*lambda_deg) &
                           * sqrt(2.0*10**(0.73 + 0.91*Kp)) / 57.6
            else
                Bw_pT = 10**(0.75 + 0.04*lambda_deg) &
                           * sqrt(2.0*10**(2.5 + 0.18*Kp)) / 57.6
            end if

        ! ----------------------------------------------------------------------
        ! TYPICALNIGHTSIDECHORUS: Kp dependent only, fixed latitude factor.
        ! Reference: Hui et al. [2017] comparison case.
        ! ----------------------------------------------------------------------
        case('TYPICALNIGHTSIDECHORUS')
            if (Kp <= 2.0) then
                Bw_pT = 50.0 * sqrt(2.0*10**(0.73 + 0.91*Kp)) / 57.6
            else
                Bw_pT = 50.0 * sqrt(2.0*10**(2.5 + 0.18*Kp)) / 57.6
            end if

        ! ----------------------------------------------------------------------
        ! TYPICALHISS: simple linear Kp scaling.
        ! ----------------------------------------------------------------------
        case('TYPICALHISS')
            Bw_pT = 80.0 * Kp / 4.0

        ! ----------------------------------------------------------------------
        ! KS14DAYSIDEHISS: polynomial in L, lambda, Kp.
        ! Lambda clamped to 25 degrees maximum.
        ! ----------------------------------------------------------------------
        case('KS14DAYSIDEHISS')
            a1 =  2.5593;  a2 =  0.26074; a3 =  0.55554
            a4 = -0.1263;  a5 =  0.018152; a6 = -0.015465
            a7 = -0.14001; a8 =  0.0030924

            if (lambda_deg > 25.0) lambda_deg = 25.0

            log10_Bwsquared = a1 + a2*Kp + a3*L + a4*lambda_deg &
                            + a5*L*lambda_deg + a6*Kp**2 &
                            + a7*L**2 + a8*lambda_deg**2
            Bw_pT = sqrt(10**log10_Bwsquared)

        ! ----------------------------------------------------------------------
        ! KS14NIGHTSIDEHISS: polynomial in L, lambda, Kp.
        ! Lambda clamped to 25 degrees maximum.
        ! ----------------------------------------------------------------------
        case('KS14NIGHTSIDEHISS')
            a1 =  2.4319;   a2 = -0.088243; a3 =  0.63916
            a4 = -0.13609;  a5 = -0.012665; a6 =  0.03026
            a7 =  0.042481; a8 = -0.15691;  a9 =  0.0014566

            if (lambda_deg > 25.0) lambda_deg = 25.0

            log10_Bwsquared = a1 + a2*Kp + a3*L + a4*lambda_deg &
                            + a5*Kp*lambda_deg + a6*L*lambda_deg &
                            + a7*Kp**2 + a8*L**2 + a9*lambda_deg**2
            Bw_pT = sqrt(10**log10_Bwsquared)

        ! ----------------------------------------------------------------------
        ! 2016DAYSIDECHORUS: lambda clamped to 20 degrees maximum.
        ! ----------------------------------------------------------------------
        case('2016DAYSIDECHORUS')
            if (lambda_deg > 20.0) lambda_deg = 20.0
            log10_Bwsquared = 1.0561 + 0.4095*Kp + 0.0351*lambda_deg
            Bw_pT = sqrt(10**log10_Bwsquared)

        ! ----------------------------------------------------------------------
        ! 2016NIGHTSIDECHORUS
        ! ----------------------------------------------------------------------
        case('2016NIGHTSIDECHORUS')
            log10_Bwsquared = 1.8821 + 0.3637*Kp - 0.05015*lambda_deg
            Bw_pT = sqrt(10**log10_Bwsquared)

        ! ----------------------------------------------------------------------
        ! 2016NIGHTSIDECHORUS2: includes Kp*lambda cross term.
        ! ----------------------------------------------------------------------
        case('2016NIGHTSIDECHORUS2')
            log10_Bwsquared = 1.0662 + 0.6677*Kp + 0.0223*lambda_deg &
                            - 0.0312*Kp*lambda_deg
            Bw_pT = sqrt(10**log10_Bwsquared)

        ! ----------------------------------------------------------------------
        ! 2016HISS: L-dependent amplitude only.
        ! ----------------------------------------------------------------------
        case('2016HISS')
            log10_Bwsquared = 0.05894*L**3 - 0.7768*L**2 + 3.163*L - 1.036
            Bw_pT = sqrt(10**log10_Bwsquared)

        ! ----------------------------------------------------------------------
        ! OA2015: 4x4 polynomial in (lambda/10, Kp), MLT and L dependent.
        ! Reference: Agapitov et al. (2015), JGR.
        ! ----------------------------------------------------------------------
        case('OA2015')
            if (MLT >= 3.0 .and. MLT < 15.0) then   ! Dayside
                if (L <= 5.0) then
                    coeff(1,:) = (/ -2.53242,  0.14836, -0.01648,  0.00016 /)
                    coeff(2,:) = (/  0.11290, -0.47100,  0.23260, -0.02317 /)
                    coeff(3,:) = (/ -0.08181,  0.28823, -0.15163,  0.01726 /)
                    coeff(4,:) = (/  0.01027, -0.04109,  0.02380, -0.00311 /)
                else
                    coeff(1,:) = (/ -2.40708, -0.02297,  0.01223, -0.00021 /)
                    coeff(2,:) = (/ -0.29267,  0.45874, -0.10822,  0.00812 /)
                    coeff(3,:) = (/  0.18808, -0.22663,  0.05810, -0.00446 /)
                    coeff(4,:) = (/ -0.02839,  0.03140, -0.00859,  0.00066 /)
                end if
            else                                      ! Nightside
                if (L <= 5.0) then
                    coeff(1,:) = (/ -2.39346,  0.10988, -0.03083,  0.00431 /)
                    coeff(2,:) = (/ -0.03083,  0.16417, -0.01875, -0.00346 /)
                    coeff(3,:) = (/ -0.00884, -0.13052,  0.00501,  0.00680 /)
                    coeff(4,:) = (/  0.00334,  0.00089,  0.01714, -0.00523 /)
                else
                    coeff(1,:) = (/ -2.69548,  0.22568, -0.04232,  0.00423 /)
                    coeff(2,:) = (/ -0.00652, -0.09346,  0.05615, -0.00424 /)
                    coeff(3,:) = (/  0.07774,  0.12785, -0.06845,  0.00524 /)
                    coeff(4,:) = (/ -0.01756, -0.02291,  0.01166, -0.00080 /)
                end if
            end if

            Bw_Input1 = 0.0
            do i = 0, 3
                do j = 0, 3
                    Bw_Input1 = Bw_Input1 &
                                + coeff(i+1,j+1) * Kp**j * (lambda_deg/10.0)**i
                end do
            end do
            Bw_Input2 = 10**Bw_Input1
            Bw_pT  = 1000.0 * Bw_Input2

        ! ----------------------------------------------------------------------
        ! VITALIICHORUS: 4x4 polynomial in (lambda/10, L) scaled by Kp polynomial.
        ! Lambda clamped to 20 degrees maximum.
        ! ----------------------------------------------------------------------
        case('VITALIICHORUS')
            if (lambda_deg > 20.0) lambda_deg = 20.0

            if (MLT >= 3.0 .and. MLT < 15.0) then   ! Dayside
                aa(1,:) = (/ -5.387e-3,  7.326e-3, -2.976e-3,  4.128e-4 /)
                aa(2,:) = (/  8.242e-2, -8.720e-2,  2.861e-2, -2.777e-3 /)
                aa(3,:) = (/ -1.302e-1,  1.309e-1, -4.070e-2,  3.759e-3 /)
                aa(4,:) = (/  5.032e-2, -4.968e-2,  1.519e-2, -1.389e-3 /)
                bb = (/  8.408e-1, -2.484e-1,  5.387e-1, -6.574e-2 /)
            else                                      ! Nightside
                aa(1,:) = (/  8.671e-3, -8.872e-3,  2.848e-3, -2.689e-4 /)
                aa(2,:) = (/  3.613e-3, -5.749e-3,  2.599e-3, -3.138e-4 /)
                aa(3,:) = (/ -1.691e-2,  2.035e-2, -7.574e-3,  8.172e-4 /)
                aa(4,:) = (/  6.087e-3, -7.417e-3,  2.781e-3, -3.013e-4 /)
                bb = (/ 3.684e0, -2.46e0, 2.828e0, -3.018e-1 /)
            end if

            Bw_Input1 = 0.0
            do ii = 1, 4
                do jj = 1, 4
                    Bw_Input1 = Bw_Input1 &
                                + aa(ii,jj) * (lambda_deg/10.0)**(ii-1) * L**(jj-1)
                end do
            end do

            Bw_Input2 = 0.0
            do ii = 1, 4
                Bw_Input2 = Bw_Input2 + bb(ii) * Kp**(ii-1)
            end do

            Bw_pT = Bw_Input1 * Bw_Input2 * 1e3   ! nT to pT

        ! ----------------------------------------------------------------------
        ! 2017CHORUS: polynomial in L, lambda, MLT (cos/sin), scaled by Kp.
        ! L clamped to [3.5, 6], lambda to 20 degrees maximum.
        ! ----------------------------------------------------------------------
        case('2017CHORUS')
            if (lambda_deg >= 20.0) lambda_deg = 20.0
            Ltemp = max(3.5d0, min(L, 6.0d0))

            aaa = cos(MLT * PI / 12.0)
            bbb = sin(MLT * PI / 12.0)

            Bw_pT = -3266.7 + 546.23*aaa - 670.18*bbb + 1150.0*Ltemp &
                     + 87.776*lambda_deg + 199.54*aaa*bbb &
                     + 228.42*bbb*Ltemp - 77.367*Ltemp**2 &
                     - 54.412*aaa*lambda_deg - 10.025*bbb*lambda_deg &
                     - 18.864*Ltemp*lambda_deg
            Bw_pT = Bw_pT &
                     * (75.532*Kp - 45.677*Kp**2 + 50.591*Kp**3 &
                        - 6.265*Kp**4 + 0.53362*Kp**5) / 1478.1
            if (Bw_pT < 1.0e-2) Bw_pT = 1.0e-2
            Bw_pT = sqrt(Bw_pT)

        ! ----------------------------------------------------------------------
        ! 2017UPPERBANDCHORUS: polynomial in L, lambda, MLT, scaled by Kp.
        ! L clamped to [3.5, 6], lambda to 17 degrees maximum.
        ! ----------------------------------------------------------------------
        case('2017UPPERBANDCHORUS')
            if (lambda_deg >= 17.0) lambda_deg = 17.0
            Ltemp = max(3.5d0, min(L, 6.0d0))

            aaa = cos(MLT * PI / 12.0)
            bbb = sin(MLT * PI / 12.0)

            Bw_pT = 510.0363206 + 206.0115305*aaa + 294.2069939*bbb &
                     - 64.57977972*Ltemp - 48.76766076*lambda_deg &
                     + 27.89711269*aaa*bbb - 26.18987181*aaa*Ltemp &
                     - 5.740739266*aaa*lambda_deg - 27.68773262*bbb*Ltemp &
                     - 10.6820023*bbb*Ltemp + 0.2699685419*Ltemp**2 &
                     + 4.326157077*Ltemp*lambda_deg + 0.9335134531*lambda_deg**2
            Bw_pT = Bw_pT &
                     * (-9.752027449*Kp + 26.46676219*Kp**2 &
                        - 5.243283335*Kp**3 + 0.3785342851*Kp**4) / 103.3238506
            if (Bw_pT < 1.0e-2) Bw_pT = 1.0e-2
            Bw_pT = sqrt(Bw_pT)

        ! ----------------------------------------------------------------------
        ! CHORUSLAT1545: latitude 15-45 degree chorus model.
        ! L clamped to [3.5, 6].
        ! BUG: the first MLT branch condition is MLT.ge.13.and.MLT.Lt.13
        ! which can never be true. 
        ! ----------------------------------------------------------------------
        case('CHORUSLAT1545')
            Ltemp = max(3.5d0, min(L, 6.0d0))

            if (Ltemp <= 5.0) then
                if (MLT >= 13.0 .and. MLT < 13.0) then   ! BUG: never true
                    a_0 =  0.0104; a_1 = -0.8551; a_2 = 21.2939
                    b_0 = -0.0192; b_1 =  0.1750; b_2 = -0.4159
                    b_3 =  0.4302; b_4 =  0.2282
                else if (MLT >= 13.0 .and. MLT < 22.0) then
                    a_0 =  0.0104; a_1 = -0.8551; a_2 = 21.2939
                    b_0 = -0.0034; b_1 = -0.0136; b_2 =  0.2681
                    b_3 = -0.3412; b_4 =  0.4098
                else
                    a_0 = -0.0070; a_1 =  0.2295; a_2 =  3.9511
                    b_0 =  0.0241; b_1 = -0.1892; b_2 =  0.3553
                    b_3 =  0.0087; b_4 =  0.6821
                end if
            else
                if (MLT >= 4.0 .and. MLT < 14.0) then
                    a_0 = -0.0173; a_1 =  0.8612; a_2 =  4.1508
                    b_0 = -0.0082; b_1 =  0.0941; b_2 = -0.3079
                    b_3 =  0.4680; b_4 =  0.3394
                else if (MLT >= 14.0 .and. MLT < 21.0) then
                    a_0 = -0.0084; a_1 =  0.3892; a_2 =  6.1371
                    b_0 = -0.0036; b_1 =  0.0327; b_2 = -0.1312
                    b_3 =  0.2713; b_4 =  0.9189
                else
                    a_0 = -0.0058; a_1 =  0.4085; a_2 = -2.6257
                    b_0 = -0.0160; b_1 =  0.1826; b_2 = -0.6751
                    b_3 =  0.8386; b_4 =  0.8622
                end if
            end if

            g_Kp     = b_0*Kp**4 + b_1*Kp**3 + b_2*Kp**2 + b_3*Kp + b_4
            Bw_pT = a_0*lambda_deg**2 + a_1*lambda_deg + a_2
            Bw_pT = Bw_pT * g_Kp
            if (Bw_pT < 1.0e-2) Bw_pT = 1.0e-2

        ! ----------------------------------------------------------------------
        ! MIDLATCHORUS: mid-latitude chorus model.
        ! L clamped to [3.5, 6], lambda to 45 degrees maximum.
        ! ----------------------------------------------------------------------
        case('MIDLATCHORUS')
            if (lambda_deg >= 45.0) lambda_deg = 45.0
            Ltemp = max(3.5d0, min(L, 6.0d0))

            if (Ltemp <= 5.0) then
                if (MLT >= 4.0 .and. MLT < 13.0) then
                    a_0 =  2.5452653e-2; a_1 = -1.8204229;   a_2 =  3.5917442e1
                    b_0 = -2.1785381e-2; b_1 =  2.1934154e-1; b_2 = -6.5135800e-1
                    b_3 =  8.1637208e-1; b_4 =  1.6950750e-1
                else if (MLT >= 13.0 .and. MLT < 22.0) then
                    a_0 =  5.8055846e-3; a_1 = -5.5876391e-1; a_2 =  1.5572106e1
                    b_0 = -1.7678770e-2; b_1 =  1.6252408e-1; b_2 = -3.9468309e-1
                    b_3 =  4.0186281e-1; b_4 =  3.3664550e-1
                else
                    a_0 =  3.0811025e-3; a_1 = -4.0999943e-1; a_2 =  1.3541850e1
                    b_0 =  1.0852508e-2; b_1 = -2.0289397e-2; b_2 = -3.3138473e-1
                    b_3 =  9.3497384e-1; b_4 =  4.5435250e-1
                end if
            else
                if (MLT >= 4.0 .and. MLT < 14.0) then
                    a_0 = -2.5864386e-2; a_1 =  1.4217082;    a_2 = -4.7289716
                    b_0 = -1.2138862e-2; b_1 =  1.3788319e-1; b_2 = -4.4842546e-1
                    b_3 =  5.8112700e-1; b_4 =  3.3842536e-1
                else if (MLT >= 14.0 .and. MLT < 21.0) then
                    a_0 = -1.2138022e-2; a_1 =  6.4269273e-1; a_2 =  2.0069819
                    b_0 = -2.4564293e-3; b_1 =  2.2152503e-2; b_2 = -9.7797716e-2
                    b_3 =  2.0362217e-1; b_4 =  1.0031774
                else
                    a_0 = -1.1856737e-2; a_1 =  8.1229858e-1; a_2 = -9.2016650
                    b_0 = -1.1385737e-2; b_1 =  1.3457093e-1; b_2 = -5.3830828e-1
                    b_3 =  7.5101515e-1; b_4 =  8.9100649e-1
                end if
            end if

            g_Kp     = b_0*Kp**4 + b_1*Kp**3 + b_2*Kp**2 + b_3*Kp + b_4
            Bw_pT = a_0*lambda_deg**2 + a_1*lambda_deg + a_2
            Bw_pT = Bw_pT * g_Kp
            if (Bw_pT < 1.0e-2) Bw_pT = 1.0e-2

        ! ----------------------------------------------------------------------
        ! 2018UBCHORUS: L-and-MLT-dependent G0/g_Kp scaling, polynomial
        ! in L and lambda. L clamped to [3.5, 6].
        ! ----------------------------------------------------------------------
        case('2018UBCHORUS')
            Ltemp = max(3.5d0, min(L, 6.0d0))

            if (Ltemp <= 5.0) then
                if (MLT >= 6.0 .and. MLT < 12.0) then
                    G0 = 8.218
                    b_0 = 0.0;  b_1 = 0.189; b_2 = 0.301; b_3 = 0.348; b_4 = -0.057
                else if (MLT >= 12.0 .and. MLT < 22.0) then
                    G0 = 1.392
                    b_0 = 0.0;  b_1 = 0.011; b_2 = 0.012; b_3 = -0.003; b_4 = 0.005
                else
                    G0 = 11.633
                    b_0 = 0.0;  b_1 = 0.070; b_2 = 0.202; b_3 = 0.335; b_4 = -0.035
                end if
            else
                if (MLT >= 7.0 .and. MLT < 13.0) then
                    G0 = 5.878
                    b_0 = 0.0;  b_1 = 0.771; b_2 = 0.901; b_3 = -0.048; b_4 = -0.018
                else if (MLT >= 13.0 .and. MLT < 21.0) then
                    G0 = 2.099
                    b_0 = 0.0;  b_1 = 0.018; b_2 = 0.066; b_3 = 0.135; b_4 = -0.023
                else
                    G0 = 9.754
                    b_0 = 0.0;  b_1 = 1.469; b_2 = 2.086; b_3 = -0.534; b_4 = 0.035
                end if
            end if

            g_Kp = (b_0 + b_1*Kp + b_2*Kp**2 + b_3*Kp**3 + b_4*Kp**4) / G0
            if (g_Kp < 1.0e-2) g_Kp = 1.0e-2

            aaa = cos(MLT * PI / 12.0)
            bbb = sin(MLT * PI / 12.0)
            Bw_pT = -0.896 + 23.216*aaa + 15.102*bbb + 6.580*Ltemp &
                     + 0.155*lambda_deg + 0.797*aaa*bbb - 2.943*aaa*Ltemp &
                     - 1.148*aaa*lambda_deg - 1.225*bbb*Ltemp &
                     - 0.796*bbb*lambda_deg - 0.808*Ltemp**2 &
                     - 0.057*Ltemp*lambda_deg - 0.092*lambda_deg**2
            if (Bw_pT < 1.0e-2) Bw_pT = 1.0e-2
            Bw_pT = Bw_pT * g_Kp

        ! ----------------------------------------------------------------------
        ! 2018LBCHORUS: L-and-MLT-dependent G0/g_Kp scaling, polynomial
        ! in L and lambda. L clamped to [3.5, 6], lambda to 20 degrees.
        ! ----------------------------------------------------------------------
        case('2018LBCHORUS')
            if (lambda_deg >= 20.0) lambda_deg = 20.0
            Ltemp = max(3.5d0, min(L, 6.0d0))

            if (Ltemp <= 5.0) then
                if (MLT >= 4.0 .and. MLT < 13.0) then
                    G0 = 18.308
                    b_0 = 0.0; b_1 = 0.064; b_2 = 0.179; b_3 = 0.297; b_4 = -0.001
                else if (MLT >= 13.0 .and. MLT < 22.0) then
                    G0 = 4.101
                    b_0 = 0.0; b_1 = 0.036; b_2 = 0.097; b_3 = 0.167; b_4 = -0.024
                else
                    G0 = 23.141
                    b_0 = 0.0; b_1 = 0.083; b_2 = 0.278; b_3 = 0.556; b_4 = -0.040
                end if
            else
                if (MLT >= 4.0 .and. MLT < 14.0) then
                    G0 = 27.622
                    b_0 = 0.0; b_1 = 0.194; b_2 = 0.456; b_3 = 0.724; b_4 = -0.068
                else if (MLT >= 14.0 .and. MLT < 21.0) then
                    G0 = 13.783
                    b_0 = 0.0; b_1 = 0.241; b_2 = 0.471; b_3 = 0.710; b_4 = -0.119
                else
                    G0 = 30.181
                    b_0 = 0.0; b_1 = 1.461; b_2 = 2.632; b_3 = 0.184; b_4 = -0.061
                end if
            end if

            g_Kp = (b_0 + b_1*Kp + b_2*Kp**2 + b_3*Kp**3 + b_4*Kp**4) / G0
            if (g_Kp < 1.0e-2) g_Kp = 1.0e-2

            aaa = cos(MLT * PI / 12.0)
            bbb = sin(MLT * PI / 12.0)
            Bw_pT = -136.706 + 35.308*aaa + 4.571*bbb + 49.943*Ltemp &
                     + 4.626*lambda_deg + 2.410*aaa*bbb - 3.871*aaa*Ltemp &
                     - 1.780*aaa*lambda_deg + 1.739*bbb*Ltemp &
                     - 0.462*bbb*lambda_deg - 3.801*Ltemp**2 &
                     - 0.848*Ltemp*lambda_deg - 0.017*lambda_deg**2
            if (Bw_pT < 1.0e-2) Bw_pT = 1.0e-2
            Bw_pT = Bw_pT * g_Kp

        ! ----------------------------------------------------------------------
        ! 2026LBCHORUS: L-and-MLT-dependent g_Kp scaling, polynomial
        ! in L and lambda. L clamped to [3.5, 6.5], lambda to 40 degrees.
        ! ----------------------------------------------------------------------
        case('2026LBCHORUS')
            if (lambda_deg > 40.0) lambda_deg = 40.0
            Ltemp = max(3.5d0, min(L, 6.5d0))

            G0  = 29.169
            b_0 =  8.847; b_1 = 1.896; b_2 = 1.834; b_3 = -0.247; b_4 = 0.014

            g_Kp = (b_0 + b_1*Kp + b_2*Kp**2 + b_3*Kp**3 + b_4*Kp**4) / G0
            if (g_Kp < 1.0e-2) g_Kp = 1.0e-2

            aaa = cos(MLT * PI / 12.0)
            bbb = sin(MLT * PI / 12.0)
            Bw_pT = -79.283 + 33.851*aaa + 10.731*bbb + 31.739*Ltemp &
                     + 1.985*lambda_deg + 1.258*aaa*bbb - 4.564*aaa*Ltemp &
                     - 1.253*aaa*lambda_deg + 0.48*bbb*Ltemp &
                     - 0.549*bbb*lambda_deg - 2.497*Ltemp**2 &
                     - 0.405*Ltemp*lambda_deg + 0.013*lambda_deg**2
            if (Bw_pT < 1.0e-2) Bw_pT = 1.0e-2
            Bw_pT = Bw_pT * g_Kp

        ! ----------------------------------------------------------------------
        ! 2026UBCHORUS: L-and-MLT-dependent g_Kp scaling, polynomial
        ! in L and lambda. L clamped to [4, 6.5], lambda to 10 degrees.
        ! ----------------------------------------------------------------------
        case('2026UBCHORUS')
            if (lambda_deg > 10.0) lambda_deg = 10.0
            Ltemp = max(4.0d0, min(L, 6.5d0))

            G0  = 12.449
            b_0 = 2.283; b_1 = 1.125; b_2 = 1.429; b_3 = -0.319; b_4 = 0.023

            g_Kp = (b_0 + b_1*Kp + b_2*Kp**2 + b_3*Kp**3 + b_4*Kp**4) / G0
            if (g_Kp < 1.0e-2) g_Kp = 1.0e-2

            aaa = cos(MLT * PI / 12.0)
            bbb = sin(MLT * PI / 12.0)
            Bw_pT = 18.59 + 20.906*aaa + 15.359*bbb - 0.408*Ltemp &
                     - 1.055*lambda_deg + 0.843*aaa*bbb - 2.64*aaa*Ltemp &
                     - 1.043*aaa*lambda_deg - 1.353*bbb*Ltemp &
                     - 0.749*bbb*lambda_deg - 0.202*Ltemp**2 &
                     + 0.13*Ltemp*lambda_deg - 0.062*lambda_deg**2
            if (Bw_pT < 1.0e-2) Bw_pT = 1.0e-2
            Bw_pT = Bw_pT * g_Kp

        ! ----------------------------------------------------------------------
        ! NBBLBC: lower-band chorus. Kp and latitude dependent.
        ! ----------------------------------------------------------------------
        case('NBBLBC')
            if (Kp == 0.0) then
                if     (lambda_deg <= 5.0)  then; Bw_pT = 10.0
                else if(lambda_deg <= 10.0) then; Bw_pT = 10.0
                else;                             Bw_pT = 10.0
                end if
            else if (Kp == 1.0) then
                if     (lambda_deg <= 5.0)  then; Bw_pT = 20.0
                else if(lambda_deg <= 10.0) then; Bw_pT = 35.0
                else;                             Bw_pT = 20.0
                end if
            else
                if     (lambda_deg <= 5.0)  then; Bw_pT = 65.0
                else if(lambda_deg <= 10.0) then; Bw_pT = 80.0
                else;                             Bw_pT = 60.0
                end if
            end if

        ! ----------------------------------------------------------------------
        ! NBBUBC: NBB upper-band chorus. Kp and latitude dependent.
        ! ----------------------------------------------------------------------
        case('NBBUBC')
            if (Kp == 0.0) then
                if (lambda_deg <= 5.0) then; Bw_pT = 5.0
                else;                        Bw_pT = 5.0
                end if
            else if (Kp == 1.0) then
                if (lambda_deg <= 5.0) then; Bw_pT = 15.0
                else;                        Bw_pT = 10.0
                end if
            else
                if (lambda_deg <= 5.0) then; Bw_pT = 45.0
                else;                        Bw_pT = 20.0
                end if
            end if

        ! ----------------------------------------------------------------------
        ! CONSTANT: fixed amplitude read from the input file.
        ! ----------------------------------------------------------------------
        case('CONSTANT')
            Bw_pT = Bw_Input          ! module variable, set from the input file

        case default
            write(*,*) 'ERROR: Invalid BwType = ', trim(BwType)
            stop

    end select

    ! Convert from pT to T
    Bw = Bw_pT * 1.0d-12

    return

end subroutine Bw_wave