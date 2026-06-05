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
! NOTES:
!   - argJ is declared intent(inout) because it is clamped in place if
!     |argJ| > 1e5. Callers should be aware their variable may be modified.
!     A future improvement would be to use a local copy and declare argJ
!     as intent(in), avoiding the side effect on the caller.
! ==============================================================================

subroutine Bessel_FirstKind(argJ, n, BesselValue)
implicit none
 
    real*8,  intent(inout) :: argJ        ! modified if |argJ| > 1e5
    integer, intent(in)    :: n
    real*8,  intent(out)   :: BesselValue
 
    real*8 :: temp
 
    ! Clamp argument to avoid overflow in BesselFunction_firstkind
    if (abs(argJ) > 1.0e5) then
        if (argJ > 0.0) then
            argJ =  1.0e5
        else
            argJ = -1.0e5
        end if
    end if
 
    if (n >= 0) then
        call BesselFunction_firstkind(n, argJ, temp)
        BesselValue = temp
    else
        ! J_{-n}(x) = (-1)^n * J_n(x)
        call BesselFunction_firstkind(abs(n), argJ, temp)
        BesselValue = ((-1.0)**real(abs(n))) * temp
    end if
 
    return
 
end subroutine Bessel_FirstKind
