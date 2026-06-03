! Get the value of Bessel Function of First Kind with integer order and real argument
! Please note here the integer order can be negative or zero, we have determine this order first

subroutine Bessel_FirstKind(argJ,n,BesselValue)

!use imsl

implicit none

  real*8::argJ
  integer::n
  real*8::BesselValue
  real*8::temp

  if(abs(argJ)>1.0e5) then
    if(argJ>0) then
	  argJ=1.0e5
    else
	  argJ=-1.0e5
	end if
  end if

  if(n>=0) then
	call BesselFunction_firstkind(n,argJ,temp)
	BesselValue=temp
  else
    call BesselFunction_firstkind(abs(n),argJ,temp)
    BesselValue=((-1.0)**real(abs(n)))*temp
  end if

  return

end subroutine Bessel_FirstKind