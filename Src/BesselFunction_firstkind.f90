! Based on Fortran Code provided by Xu, ShiLiang (PKU)

subroutine BesselFunction_firstkind(N, x, BFvalue)

implicit none

  integer::N
  real*8::x
  real*8::BFvalue
  real*8::T, Y, Z, P, Q, S, B0, B1, M
  real*8::A(6), B(6), C(6), D(6), E(5), F(5), G(5), H(5)

  integer::i
  
  Data A/57568490574.0, -13362590354.0, 651619640.7, -11214424.18, 77392.33017, -184.9052456/
  Data B/57568490411.0, 1029532985.0, 9494680.718, 59272.64853, 267.8532712, 1.0/
  Data C/72362614232.0, -7895059235.0, 242396853.1, -2972611.439, 15704.4826, -30.16036606/
  Data D/144725228443.0, 2300535178.0, 18583304.74, 99447.43394, 376.9991397, 1.0/
  Data E/1.0, -0.1098628627D-02, 0.2734510407D-04, -0.2073370639D-05, 0.2093887211D-06/
  Data F/-0.1562499995D-01, 0.1430488765D-03, -0.6911147651D-05, 0.7621095161D-06, -0.934935152D-7/
  Data G/1.0, 0.183105D-02, -0.3516396496D-04, 0.2457520174D-05, -0.240337019D-06/
  Data H/0.4687499995D-01, -0.2002690873D-03, 0.8449199096D-05, -0.88228987D-06, 0.105787412D-06/ 

  T=abs(X)

  if(N.LT.0) N=-N
  if(N.NE.1) then
    if(T.LT.8.0) then
	  Y=T*T
	  P=A(6)
	  Q=B(6)
	  do i=5,1,-1
	     P=P*Y+A(i)
		 Q=Q*Y+B(i)
      end do
	  P=P/Q
    else
	  Z=8.0/T
	  Y=Z*Z
	  P=E(5)
	  Q=F(5)
	  do i=4,1,-1
	     P=P*Y+E(i)
		 Q=Q*Y+F(i)
      end do
      S=T-0.785398164
	  P=P*cos(S)-Z*Q*sin(S)
	  P=P*sqrt(0.636619772/T)
	end if
  end if

  if(N.EQ.0) then
    BFvalue=P
	return
  end if

  B0=P
  
  if(T.LT.8.0) then
    Y=T*T
	P=C(6)
	Q=D(6)
	do i=5,1,-1
	   P=P*Y+C(i)
	   Q=Q*Y+D(i)
	end do
	P=X*P/Q
  else
    Z=8.0/T
	Y=Z*Z
	P=G(5)
	Q=H(5)
	do i=4,1,-1
	   P=P*Y+G(i)
	   Q=Q*Y+H(i)
	end do
	S=T-2.356194491
	P=P*cos(S)-Z*Q*sin(S)
    P=P*X*sqrt(0.636619772/T)/T
  end if

  if(N.EQ.1) then
    BFvalue=P
    return
  end if

  B1=P

  if(X.EQ.0.0) then
    BFvalue=0.0
	return
  end if
  
  S=2.0/T
  
  if(T.GT.(1.0*N)) then
    if(X.LT.0.0) B1=-B1
	do i=1,(N-1)
	   P=S*i*B1-B0
	   B0=B1
	   B1=P
	end do
  else
    M=sqrt(40.0*N)
	M=(N+M)/2
	M=M+M
	P=0.0
	Q=0.0
	B0=1.0
	B1=0.0
	do i=(int(M)-1),0,-1
	   T=S*(i+1)*B0-B1
	   B1=B0
	   B0=T
	   if(abs(B0).GT.(1.0D+10)) then
	     B0=B0*1.0D-10
		 B1=B1*1.0D-10
		 P=P*1.0D-10
		 Q=Q*1.0D-10
	   end if
	   if(mod(i+2,2).EQ.0) Q=Q+B0
	   if((i+1).EQ.N) P=B1
	end do
	Q=2.0*Q-B0
	P=P/Q
  end if
  
  if((X.LT.0.0).AND.(mod(N,2).EQ.1)) P=-P
  BFvalue=P
  return
 
end subroutine BesselFunction_firstkind 
