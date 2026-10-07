! ==============================================================================
! SUBROUTINE: Bessel_FirstKind
! ==============================================================================
!
! PURPOSE:
!   Returns the value of the Bessel function of the first kind J_n(x) for
!   integer order n (positive, negative, or zero) and real argument x.
!
!   For negative orders, the identity is used:
!     J_{-n}(x) = (-1)^n * J_n(x)
!
!   The actual evaluation is delegated to BesselFunction_firstkind, which
!   handles non-negative integer orders.
!
!   The argument is clamped to [-1e5, 1e5] before evaluation to avoid
!   numerical overflow in BesselFunction_firstkind. For radiation belt
!   applications the argument argJ = |k|*v_perp/omega_e is typically
!   well below this limit; the clamp is a safety guard.
!
! INPUTS:
!   argJ        - Argument of the Bessel function x = |k|*v_perp/omega_e [-]
!   n           - Integer order (may be negative, zero, or positive)
!
! OUTPUTS:
!   BesselValue - J_n(argJ) [-]
!
! CALLS:
!   BesselFunction_firstkind - Evaluates J_n(x) for non-negative integer n
!
! ==============================================================================

subroutine Bessel_FirstKind(argJ, n, BesselValue)
implicit none

    real*8,  intent(in)  :: argJ
    integer, intent(in)  :: n
    real*8,  intent(out) :: BesselValue

    real*8 :: x, temp

    ! Clamp argument to avoid overflow in BesselFunction_firstkind
    x = argJ
    if (abs(x) > 1.0d5) x = sign(1.0d5, x)

    call BesselFunction_firstkind(abs(n), x, temp)

    if (n >= 0 .or. mod(n, 2) == 0) then
        BesselValue = temp
    else
        ! J_{-n}(x) = (-1)^n * J_n(x); odd order flips the sign
        BesselValue = -temp
    end if

    return

end subroutine Bessel_FirstKind