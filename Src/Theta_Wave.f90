! ==============================================================================
! SUBROUTINE: Theta_wave
! ==============================================================================
!
! PURPOSE:
!   Returns the wave normal angle distribution parameters (Theta_function)
!   at a given L-shell and magnetic latitude, for the model specified by
!   ThetaType in module constants.
!
!   The seven output parameters define a two-Gaussian distribution in
!   tan(theta) with cutoff angles:
!     Theta_function(1) = theta_m1  : mean angle of 1st Gaussian [degrees]
!     Theta_function(2) = theta_w1  : width of 1st Gaussian [degrees]
!     Theta_function(3) = theta_m2  : mean angle of 2nd Gaussian [degrees]
!     Theta_function(4) = theta_w2  : width of 2nd Gaussian [degrees]
!     Theta_function(5) = theta_min : lower cutoff angle [degrees]
!     Theta_function(6) = theta_max : upper cutoff angle [degrees]
!     Theta_function(7) = theta_q   : relative amplitude of 2nd Gaussian
!
!   A model with theta_q = 0 reduces to a single Gaussian.
!
!   Supported models and their parameter dependencies:
!
!     Model           L    lambda   MLT   Kp   Description
!     CONSTANT                                 Fixed angles from input file
!     OBLIQUEHISS            X                 Latitude-stepping hiss model
!     OA2015          X      X       X    X    Agapitov et al. (2015)
!     VITALIICHORUS   X      X       X    X    Two-regime chorus model
!     NBBLBC                 X                 Fixed angles, latitude dependent
!     NBBUBC                 X                 Fixed angles, latitude dependent
!
! INPUTS:
!   L        - Magnetic L-shell [Earth radii, dimensionless]
!   lambda   - Magnetic latitude [radians]
!
! OUTPUTS:
!   Theta_function(7) - Wave normal angle distribution parameters [degrees]
!                       except theta_q which is dimensionless
!
! MODULE DEPENDENCIES:
!   constants - ThetaType, Theta_Input, Kp, MLT, RAD
!
! NOTES:
!   - OA2015 uses Kp*10 in the polynomial evaluation for theta_q. This
!     reflects the original fitting which used Kp on a 0-90 scale. Since
!     module constants stores Kp in its natural 0-9 scale, the *10 factor
!     is applied internally here.
!   - VITALIICHORUS uses two different approaches depending on latitude:
!       lambda_deg > 20: uses OA2015 polynomial coefficients for theta_q
!                        (code is duplicated from OA2015 case)
!       lambda_deg <= 20: uses separate coeff/coeff_bb arrays for theta_q
!     The theta_m1, theta_m2, theta_w1 polynomials (coeff_a, coeff_b) are
!     the same in both regimes and both models. This duplication is a known
!     issue to be addressed in a future refactoring step.
!   - OBLIQUEHISS has a duplicate assignment of Theta_function(6) = theta_max.
!     The duplicate has been removed in this version.
!   - coeff_a is declared as (2,6) but only indices up to (1,4) and (2,2)
!     are set; the remaining entries are zero by initialisation.
!
! ==============================================================================

subroutine Theta_wave(L, lambda, Theta_function)
use constants
implicit none

    real*8, intent(in)  :: L, lambda
    real*8, intent(out) :: Theta_function(7)

    real*8  :: lambda_deg
    real*8  :: theta_min, theta_max
    real*8  :: theta_m1, theta_w1, theta_m2, theta_w2, theta_q
    real*8  :: coeff(4,4), coeff_bb(4)
    real*8  :: coeff_a(2,6), coeff_b(6)
    real*8  :: theta_q1, theta_q2
    integer :: i, j

    lambda_deg = lambda / RAD

    select case(ThetaType)

        ! ----------------------------------------------------------------------
        ! CONSTANT: fixed angles from input file, same at all locations.
        ! ----------------------------------------------------------------------
        case('CONSTANT')
            Theta_function = Theta_Input

        ! ----------------------------------------------------------------------
        ! OBLIQUEHISS: latitude-stepping single-Gaussian model.
        ! Wave normal angle increases with latitude in discrete bands.
        ! theta_q = 0 (single Gaussian).
        ! ----------------------------------------------------------------------
        case('OBLIQUEHISS')
            if (lambda_deg < 5.0) then
                theta_min = 0.0;  theta_max = 20.0
                theta_m1  = 0.0;  theta_w1  = 15.0
            else if (lambda_deg < 10.0) then
                theta_min = 0.0;  theta_max = 30.0
                theta_m1  = 10.0; theta_w1  = 20.0
            else if (lambda_deg < 15.0) then
                theta_min = 0.0;  theta_max = 40.0
                theta_m1  = 20.0; theta_w1  = 25.0
            else if (lambda_deg < 20.0) then
                theta_min = 10.0; theta_max = 50.0
                theta_m1  = 30.0; theta_w1  = 25.0
            else if (lambda_deg < 25.0) then
                theta_min = 20.0; theta_max = 60.0
                theta_m1  = 40.0; theta_w1  = 35.0
            else if (lambda_deg < 30.0) then
                theta_min = 30.0; theta_max = 70.0
                theta_m1  = 50.0; theta_w1  = 40.0
            else if (lambda_deg < 35.0) then
                theta_min = 40.0; theta_max = 80.0
                theta_m1  = 60.0; theta_w1  = 50.0
            else if (lambda_deg < 40.0) then
                theta_min = 50.0; theta_max = 80.0
                theta_m1  = 70.0; theta_w1  = 60.0
            else
                theta_min = 55.0; theta_max = 85.0
                theta_m1  = 80.0; theta_w1  = 70.0
            end if

            theta_m2 = 0.0; theta_w2 = 1.0; theta_q = 0.0

            Theta_function(1) = theta_m1
            Theta_function(2) = theta_w1
            Theta_function(3) = theta_m2
            Theta_function(4) = theta_w2
            Theta_function(5) = theta_min
            Theta_function(6) = theta_max
            Theta_function(7) = theta_q

        ! ----------------------------------------------------------------------
        ! OA2015: two-Gaussian model from Agapitov et al. (2015).
        ! theta_q computed from 4x4 polynomial in (Kp*10, lambda_deg).
        ! theta_m1, theta_m2 from polynomial in (lambda_deg/10).
        ! theta_w1 = theta_w2 from polynomial in (lambda_deg/10).
        ! Note: Kp*10 is used in the theta_q polynomial (original fitting
        ! used Kp on a 0-90 scale).
        ! ----------------------------------------------------------------------
        case('OA2015')
            if (MLT >= 3.0 .and. MLT <= 15.0) then   ! Dayside
                if (L <= 5.0) then
                    coeff(1,:) = (/ -2.96,     5.43e-3,  1.25e-4,  -1.98e-6 /)
                    coeff(2,:) = (/ -2.25e-2,  2.57e-2, -4.71e-4,   1.4e-6  /)
                    coeff(3,:) = (/  4.96e-3, -1.57e-3,  2.61e-5,  -5.46e-8 /)
                    coeff(4,:) = (/ -7.44e-5,  2.39e-5, -3.9e-7,    6.85e-10/)
                else
                    coeff(1,:) = (/ -2.33942,   0.148424,  -0.00442111,  2.93478e-5  /)
                    coeff(2,:) = (/  0.397584,  -0.0260653,  0.000500087,-2.39098e-6  /)
                    coeff(3,:) = (/ -0.0200059,  0.00102028,-1.37305e-5,  3.316e-8   /)
                    coeff(4,:) = (/  0.000272183,-1.1708e-5, 1.08928e-7,  1.04256e-11/)
                end if
            else                                        ! Nightside
                if (L <= 5.0) then
                    coeff(1,:) = (/ -2.36144,  -0.0679,    0.003815125, -0.000043285/)
                    coeff(2,:) = (/  0.05889,   0.00836,  -0.00022328,   0.000001204/)
                    coeff(3,:) = (/ -0.00201295,-0.00041678,0.000006725,  4.719e-9   /)
                    coeff(4,:) = (/  5.096e-5,  4.968e-6, -5.058e-8,    -5.518e-10  /)
                else
                    coeff(1,:) = (/ -0.78803,  -0.02769,   1.479e-3,   -4.1015e-5  /)
                    coeff(2,:) = (/  0.0640415, 0.0108776, -7.891e-4,    1.328e-5   /)
                    coeff(3,:) = (/ -5.297e-3, -5.60918e-4, 4.53587e-5, -7.31029e-7 /)
                    coeff(4,:) = (/  9.1149e-5, 8.60196e-6,-6.94066e-7,  1.07675e-8 /)
                end if
            end if

            theta_q = 0.0
            do i = 0, 3
                do j = 0, 3
                    theta_q = theta_q + coeff(i+1,j+1) * (Kp*10)**j * lambda_deg**i
                end do
            end do
            theta_q = theta_q**2

            ! theta_m1, theta_m2, theta_w1 polynomials in lambda_deg/10
            coeff_a      = 0.0
            coeff_a(1,1) = 11.54; coeff_a(1,2) = 14.32
            coeff_a(1,3) = -8.18; coeff_a(1,4) =  1.22
            coeff_a(2,1) = 66.0;  coeff_a(2,2) =  1.0

            coeff_b = (/ 5.68, 4.62, 3.05, -5.06, 1.83, -0.203 /)

            theta_m1 = 0.0; theta_m2 = 0.0; theta_w1 = 0.0
            do i = 0, 5
                theta_m1 = theta_m1 + coeff_a(1,i+1) * (lambda_deg/10.0)**i
                theta_m2 = theta_m2 + coeff_a(2,i+1) * (lambda_deg/10.0)**i
                theta_w1 = theta_w1 + coeff_b(i+1)   * (lambda_deg/10.0)**i
            end do
            theta_w2 = theta_w1
            theta_min = 0.0; theta_max = 85.0

            Theta_function(1) = theta_m1
            Theta_function(2) = theta_w1
            Theta_function(3) = theta_m2
            Theta_function(4) = theta_w2
            Theta_function(5) = theta_min
            Theta_function(6) = theta_max
            Theta_function(7) = theta_q

        ! ----------------------------------------------------------------------
        ! VITALIICHORUS: two-regime model depending on latitude.
        !   lambda_deg > 20: uses OA2015 polynomial for theta_q
        !                    (same coeff blocks as OA2015 — code duplication)
        !   lambda_deg <= 20: uses separate coeff/coeff_bb for theta_q,
        !                     computed as theta_q1*theta_q2/(1-theta_q1*theta_q2)
        ! theta_m1, theta_m2, theta_w1 use the same coeff_a/coeff_b as OA2015
        ! in both regimes.
        ! ----------------------------------------------------------------------
        case('VITALIICHORUS')
            if (lambda_deg > 20.0) then

                ! High latitude: use OA2015 theta_q polynomial
                if (MLT >= 3.0 .and. MLT <= 15.0) then   ! Dayside
                    if (L <= 5.0) then
                        coeff(1,:) = (/ -2.96,     5.43e-3,  1.25e-4,  -1.98e-6 /)
                        coeff(2,:) = (/ -2.25e-2,  2.57e-2, -4.71e-4,   1.4e-6  /)
                        coeff(3,:) = (/  4.96e-3, -1.57e-3,  2.61e-5,  -5.46e-8 /)
                        coeff(4,:) = (/ -7.44e-5,  2.39e-5, -3.9e-7,    6.85e-10/)
                    else
                        coeff(1,:) = (/ -2.33942,   0.148424,  -0.00442111,  2.93478e-5  /)
                        coeff(2,:) = (/  0.397584,  -0.0260653,  0.000500087,-2.39098e-6  /)
                        coeff(3,:) = (/ -0.0200059,  0.00102028,-1.37305e-5,  3.316e-8   /)
                        coeff(4,:) = (/  0.000272183,-1.1708e-5, 1.08928e-7,  1.04256e-11/)
                    end if
                else                                        ! Nightside
                    if (L <= 5.0) then
                        coeff(1,:) = (/ -2.36144,  -0.0679,    0.003815125, -0.000043285/)
                        coeff(2,:) = (/  0.05889,   0.00836,  -0.00022328,   0.000001204/)
                        coeff(3,:) = (/ -0.00201295,-0.00041678,0.000006725,  4.719e-9   /)
                        coeff(4,:) = (/  5.096e-5,  4.968e-6, -5.058e-8,    -5.518e-10  /)
                    else
                        coeff(1,:) = (/ -0.78803,  -0.02769,   1.479e-3,   -4.1015e-5  /)
                        coeff(2,:) = (/  0.0640415, 0.0108776, -7.891e-4,    1.328e-5   /)
                        coeff(3,:) = (/ -5.297e-3, -5.60918e-4, 4.53587e-5, -7.31029e-7 /)
                        coeff(4,:) = (/  9.1149e-5, 8.60196e-6,-6.94066e-7,  1.07675e-8 /)
                    end if
                end if

                theta_q = 0.0
                do i = 0, 3
                    do j = 0, 3
                        theta_q = theta_q + coeff(i+1,j+1) * (Kp*10)**j * lambda_deg**i
                    end do
                end do
                theta_q = theta_q**2

            else

                ! Low latitude: use separate coeff/coeff_bb for theta_q
                if (MLT >= 3.0 .and. MLT <= 15.0) then   ! Dayside
                    coeff(1,:) = (/  1.613e-2, -1.343e-2,  3.497e-3, -2.766e-4 /)
                    coeff(2,:) = (/  3.051e-2, -3.199e-2,  1.05e-2,  -9.992e-4 /)
                    coeff(3,:) = (/ -5.596e-2,  5.255e-2, -1.531e-2,  1.318e-3 /)
                    coeff(4,:) = (/  2.028e-2, -1.794e-2,  4.905e-3, -3.958e-4 /)
                    coeff_bb   = (/  5.926, 1.643e1, 3.666, -9.124e-1 /)
                else                                        ! Nightside
                    coeff(1,:) = (/  2.49e-3,  -2.252e-3,  7.417e-4, -5.463e-5 /)
                    coeff(2,:) = (/  5.711e-2, -6.312e-2,  2.142e-2, -2.152e-3 /)
                    coeff(3,:) = (/ -1.171e-1,  1.217e-1, -3.868e-2,  3.694e-3 /)
                    coeff(4,:) = (/  4.599e-2, -4.62e-2,   1.424e-2, -1.328e-3 /)
                    coeff_bb   = (/  3.455, 9.938, 2.147, -4.175e-1 /)
                end if

                theta_q1 = 0.0
                do i = 0, 3
                    do j = 0, 3
                        theta_q1 = theta_q1 &
                                 + coeff(i+1,j+1) * (lambda_deg/10.0)**i * L**j
                    end do
                end do

                theta_q2 = 0.0
                do i = 0, 3
                    theta_q2 = theta_q2 + coeff_bb(i+1) * Kp**i
                end do

                theta_q = theta_q1 * theta_q2 / (1.0 - theta_q1*theta_q2)

            end if

            ! theta_m1, theta_m2, theta_w1 same polynomial for both regimes
            coeff_a      = 0.0
            coeff_a(1,1) = 11.54; coeff_a(1,2) = 14.32
            coeff_a(1,3) = -8.18; coeff_a(1,4) =  1.22
            coeff_a(2,1) = 66.0;  coeff_a(2,2) =  1.0
            coeff_b = (/ 5.68, 4.62, 3.05, -5.06, 1.83, -0.203 /)

            theta_m1 = 0.0; theta_m2 = 0.0; theta_w1 = 0.0
            do i = 0, 5
                theta_m1 = theta_m1 + coeff_a(1,i+1) * (lambda_deg/10.0)**i
                theta_m2 = theta_m2 + coeff_a(2,i+1) * (lambda_deg/10.0)**i
                theta_w1 = theta_w1 + coeff_b(i+1)   * (lambda_deg/10.0)**i
            end do
            theta_w2  = theta_w1
            theta_min = 0.0; theta_max = 85.0

            Theta_function(1) = theta_m1
            Theta_function(2) = theta_w1
            Theta_function(3) = theta_m2
            Theta_function(4) = theta_w2
            Theta_function(5) = theta_min
            Theta_function(6) = theta_max
            Theta_function(7) = theta_q

        ! ----------------------------------------------------------------------
        ! NBBLBC: simple latitude-stepping single-Gaussian model.
        ! ----------------------------------------------------------------------
        case('NBBLBC')
            if (lambda_deg <= 5.0) then
                theta_m1 = 0.0; theta_w1 = 30.0
                theta_min = 0.0; theta_max = 50.0
            else if (lambda_deg <= 10.0) then
                theta_m1 = 15.0; theta_w1 = 30.0
                theta_min = 0.0; theta_max = 50.0
            else
                theta_m1 = 30.0; theta_w1 = 30.0
                theta_min = 0.0; theta_max = 50.0
            end if
            theta_m2 = 0.0; theta_w2 = 1.0; theta_q = 0.0

            Theta_function(1) = theta_m1
            Theta_function(2) = theta_w1
            Theta_function(3) = theta_m2
            Theta_function(4) = theta_w2
            Theta_function(5) = theta_min
            Theta_function(6) = theta_max
            Theta_function(7) = theta_q

        ! ----------------------------------------------------------------------
        ! NBBUBC: simple latitude-stepping single-Gaussian model.
        ! ----------------------------------------------------------------------
        case('NBBUBC')
            if (lambda_deg <= 5.0) then
                theta_m1 = 0.0; theta_w1 = 30.0; theta_max = 45.0
            else
                theta_m1 = 30.0; theta_w1 = 30.0; theta_max = 45.0
            end if
            theta_m2 = 0.0; theta_w2 = 1.0
            theta_q  = 0.0; theta_min = 0.0

            Theta_function(1) = theta_m1
            Theta_function(2) = theta_w1
            Theta_function(3) = theta_m2
            Theta_function(4) = theta_w2
            Theta_function(5) = theta_min
            Theta_function(6) = theta_max
            Theta_function(7) = theta_q

        case default
            write(*,*) 'ERROR: Invalid ThetaType = ', trim(ThetaType)
            stop

    end select

    return

end subroutine Theta_wave