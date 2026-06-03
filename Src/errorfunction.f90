
SUBROUTINE errorfunction(x, ERR)
use constants,only:pi
!       =========================================
!       Purpose: Compute error function erf(x)
!       Input:   x   --- Argument of erf(x)
!       Output:  ERR --- erf(x)
!       =========================================

REAL*8,intent(in):: x
REAL*8,intent(out):: ERR

REAL*8:: c0, eps, er, r, x2
INTEGER:: k

eps = 1.0D-15
x2 = x * x
IF (ABS(x) < 3.5) THEN
  er = 1.0
  r = 1.0
  DO  k = 1, 50
    r = r * x2 / (k+0.5)
    er = er + r
    IF (ABS(r) <= ABS(er)*eps) EXIT
  END DO
  c0 = 2.0 / SQRT(pi) * x * EXP(-x2)
  ERR = c0 * er
ELSE
  er = 1.0
  r = 1.0
  DO  k = 1, 12
    r = -r * (k-0.5) / x2
    er = er + r
  END DO
  c0 = EXP(-x2) / (ABS(x)*SQRT(pi))
  ERR = 1.0 - c0 * er
  IF (x < 0.0) ERR = -ERR
END IF
RETURN
END SUBROUTINE errorfunction
