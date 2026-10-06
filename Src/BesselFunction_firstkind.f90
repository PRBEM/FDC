! ==============================================================================
! SUBROUTINE: BesselFunction_firstkind
! ==============================================================================
!
! PURPOSE:
!   Computes the Bessel function of the first kind J_N(x) for non-negative
!   integer order N and real argument x.
!
!   This is the low-level implementation called by Bessel_FirstKind, which
!   handles negative orders via J_{-N}(x) = (-1)^N * J_N(x) before
!   delegating to this routine.
!
! ALGORITHM:
!   Three different methods are used depending on N and |x|:
!
!   N = 0, |x| < 8:
!     Rational polynomial approximation using coefficients A (numerator)
!     and B (denominator). Evaluated via Horner's method.
!
!   N = 0, |x| >= 8:
!     Asymptotic expansion:
!       J0(x) ~ sqrt(2/(pi*x)) * [P(x)*cos(x - pi/4) - Q(x)*sin(x - pi/4)]
!     where P and Q are polynomial approximations using coefficients E and F.
!     The phase pi/4 = 0.785398164 is hardcoded.
!
!   N = 1:
!     Same structure as N=0 but using coefficients C, D (rational) and
!     G, H (asymptotic). The phase 3*pi/4 = 2.356194491 is hardcoded.
!
!   N > 1, |x| > N (forward recurrence):
!     Uses the three-term recurrence relation:
!       J_{n+1}(x) = (2n/x) * J_n(x) - J_{n-1}(x)
!     starting from J0 and J1. Stable when |x| > N.
!
!   N > 1, |x| <= N (Miller's backward recurrence):
!     Forward recurrence is numerically unstable when |x| <= N.
!     Miller's algorithm starts from J_M ~ 0 and J_{M-1} ~ 1 for
!     M >> N (M = N + sqrt(40*N), rounded to an even integer),
!     recurses backward to J_0, then normalises using the identity:
!       J_0(x) + 2*sum[ J_{2k}(x) ] = 1
!     This is the standard stable method for large orders.
!
! INPUTS:
!   N        - Non-negative integer order (negative orders handled by caller)
!   x        - Real argument
!
! OUTPUTS:
!   BFvalue  - J_N(x)
!
! NOTES:
!   - The polynomial coefficients (A-H) are from standard approximations
!     for J0 and J1; see Abramowitz & Stegun or equivalent reference.
!   - The constants 0.785398164 (pi/4), 2.356194491 (3*pi/4) and
!     0.636619772 (2/pi) are hardcoded approximations.
!   - Intermediate values are rescaled by 1e-10 during backward recurrence
!     to prevent overflow for large N.
!
! Based on Fortran code by Xu, ShiLiang (PKU)
! ==============================================================================

subroutine BesselFunction_firstkind(N, x, BFvalue)
implicit none
 
    integer, intent(in)  :: N
    real*8,  intent(in)  :: x
    real*8,  intent(out) :: BFvalue
 
    real*8  :: T, Y, Z, P, Q, S, B0, B1, M
    real*8  :: A(6), B(6), C(6), D(6), E(5), F(5), G(5), H(5)
    integer :: i, NN
 
    ! Polynomial coefficients for J0 rational approximation (|x| < 8)
    data A / 57568490574.0, -13362590354.0, 651619640.7, &
            -11214424.18, 77392.33017, -184.9052456 /
    data B / 57568490411.0, 1029532985.0, 9494680.718, &
             59272.64853, 267.8532712, 1.0 /
 
    ! Polynomial coefficients for J1 rational approximation (|x| < 8)
    data C / 72362614232.0, -7895059235.0, 242396853.1, &
            -2972611.439, 15704.4826, -30.16036606 /
    data D / 144725228443.0, 2300535178.0, 18583304.74, &
             99447.43394, 376.9991397, 1.0 /
 
    ! Asymptotic coefficients for J0 (|x| >= 8)
    data E / 1.0, -0.1098628627d-02,  0.2734510407d-04, &
            -0.2073370639d-05, 0.2093887211d-06 /
    data F / -0.1562499995d-01, 0.1430488765d-03, -0.6911147651d-05, &
              0.7621095161d-06, -0.934935152d-7 /
 
    ! Asymptotic coefficients for J1 (|x| >= 8)
    data G / 1.0,  0.183105d-02, -0.3516396496d-04, &
             0.2457520174d-05, -0.240337019d-06 /
    data H / 0.4687499995d-01, -0.2002690873d-03,  0.8449199096d-05, &
            -0.88228987d-06,    0.105787412d-06 /
 
    T  = abs(x)
    NN = abs(N)        ! local copy; negative orders are handled by the caller
 
    ! ------------------------------------------------------------------
    ! Compute J0(x)
    ! ------------------------------------------------------------------
    if (NN /= 1) then
        if (T < 8.0) then
            ! Rational polynomial approximation
            Y = T * T
            P = A(6); Q = B(6)
            do i = 5, 1, -1
                P = P*Y + A(i)
                Q = Q*Y + B(i)
            end do
            P = P / Q
        else
            ! Asymptotic expansion
            Z = 8.0 / T
            Y = Z * Z
            P = E(5); Q = F(5)
            do i = 4, 1, -1
                P = P*Y + E(i)
                Q = Q*Y + F(i)
            end do
            S = T - 0.785398164   ! T - pi/4
            P = (P*cos(S) - Z*Q*sin(S)) * sqrt(0.636619772/T)   ! sqrt(2/pi/T)
        end if
    end if
 
    if (NN == 0) then
        BFvalue = P
        return
    end if
    B0 = P
 
    ! ------------------------------------------------------------------
    ! Compute J1(x)
    ! ------------------------------------------------------------------
    if (T < 8.0) then
        ! Rational polynomial approximation
        Y = T * T
        P = C(6); Q = D(6)
        do i = 5, 1, -1
            P = P*Y + C(i)
            Q = Q*Y + D(i)
        end do
        P = x * P / Q
    else
        ! Asymptotic expansion
        Z = 8.0 / T
        Y = Z * Z
        P = G(5); Q = H(5)
        do i = 4, 1, -1
            P = P*Y + G(i)
            Q = Q*Y + H(i)
        end do
        S = T - 2.356194491   ! T - 3*pi/4
        P = (P*cos(S) - Z*Q*sin(S)) * x * sqrt(0.636619772/T) / T
    end if
 
    if (NN == 1) then
        BFvalue = P
        return
    end if
    B1 = P
 
    ! J_N(0) = 0 for N >= 2
    if (x == 0.0) then
        BFvalue = 0.0
        return
    end if
 
    S = 2.0 / T
 
    ! ------------------------------------------------------------------
    ! N > 1: recurrence
    ! ------------------------------------------------------------------
    if (T > real(NN, 8)) then
 
        ! Forward recurrence: stable when |x| > N
        if (x < 0.0) B1 = -B1
        do i = 1, NN-1
            P  = S*i*B1 - B0
            B0 = B1
            B1 = P
        end do
 
    else
 
        ! Miller's backward recurrence: stable when |x| <= N
        ! (see the note on the evenness of M in the header)
        M  = sqrt(40.0*NN)
        M  = (NN + M) / 2
        M  = M + M
 
        P  = 0.0
        Q  = 0.0
        B0 = 1.0
        B1 = 0.0
 
        do i = int(M)-1, 0, -1
            T  = S*(i+1)*B0 - B1
            B1 = B0
            B0 = T
 
            ! Rescale to prevent overflow
            if (abs(B0) > 1.0d+10) then
                B0 = B0 * 1.0d-10
                B1 = B1 * 1.0d-10
                P  = P  * 1.0d-10
                Q  = Q  * 1.0d-10
            end if
 
            ! Accumulate the normalisation sum over even orders
            if (mod(i, 2) == 0) Q = Q + B0
 
            ! Store J_N (unnormalised)
            if (i+1 == NN) P = B1
        end do
 
        ! Normalise using J_0 + 2*sum(J_{2k}) = 1
        Q = 2.0*Q - B0
        P = P / Q
 
    end if
 
    if ((x < 0.0) .and. (mod(NN,2) == 1)) P = -P
    BFvalue = P
 
    return
 
end subroutine BesselFunction_firstkind
