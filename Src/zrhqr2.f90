!!!! getting the eigenvalues by using the LAPACK library


subroutine zrhqr2(A,WR,WI)

implicit none

  integer,parameter::N2=10
  !integer,parameter::MAXM=50

  real*8,intent(in)::A(N2+1)
  real*8,intent(out)::WR(N2)
  real*8,intent(out)::WI(N2)

  integer::j,k

  real*8::hess(N2,N2),xr,xi


  integer,parameter::LDA=10 
  integer,parameter::LDVL=10 
  integer,parameter::LDVR=10
  integer,parameter::LWORK=45!40

  integer::info

  real*8::VL(LDVL,N2), VR(LDVR,N2)
  
  integer::WORK(LWORK)
  !vl(:)=0
  !vr(:)=0
  WORK(:)=0                   

  !if((m.gt.MAXM).or.(a(m+1).eq.0.))  pause 'bad args in zrhqr'

  !! construct the matrix

  hess(1,:)=-a(N2:1:-1)/a(N2+1)
  hess(2:N2,:)=0.0
  do k=1,N2-1
     hess(k+1,k)=1.0
  end do

  !write(*,*)hess

  !call balanc(hess,N2)
  !write(*,*)hess,m
  !call hqr(hess,m,rtr,rti)
  
  call dgeev('N', 'N', N2, hess, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO)

  !SUBROUTINE SGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR,
  !   $                  LDVR, WORK, LWORK, INFO )

  do j=2,N2
     xr=WR(j)
	 xi=WI(j)
	 do k=j-1,1,-1
	    if(WR(k).le.xr) goto 1
		 WR(k+1)=WR(k)
		 WI(k+1)=WI(k)
	 end do
     k=0
1	 WR(k+1)=xr
	 WI(k+1)=xi
  end do

  return

end subroutine zrhqr2
  
