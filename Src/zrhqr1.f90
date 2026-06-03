!!!! getting the eigenvalues by using the LAPACK library


subroutine zrhqr1(A,WR,WI)

implicit none

  integer,parameter::N2=6
  !integer,parameter::MAXM=50

  real*8,intent(in)::A(N2+1)
  real*8,intent(out)::WR(N2)
  real*8,intent(out)::WI(N2)

  integer::j,k

  real*8::hess(N2,N2),xr,xi


  integer,parameter::LDA=6
  integer,parameter::LDVL=6 
  integer,parameter::LDVR=6
  integer,parameter::LWORK=30!24

!  integer,parameter::ILO=1
!  integer,parameter::IHI=6
!  integer,parameter::LDZ=6
!   integer,parameter::LWORK1=66
!    integer::WORK1(LWORK1)

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

!!!!! New matrix changed by Orlova, 2014
!  hess(1:N2,:)=0.0
!  hess(:,N2)=-a(1:N2:1)/a(N2+1)
!  do k=1,N2-1
!     hess(k+1,k)=1.0
!  end do

!  write(*,*)hess(1,:)
!  write(*,*)hess(2,:)
!  write(*,*)hess(3,:)
!  write(*,*)hess(4,:)
!  write(*,*)hess(5,:)
!  write(*,*)hess(6,:)

  !call balanc(hess,N2)
  !write(*,*)hess,m
  !call hqr(hess,m,rtr,rti)
  
  call dgeev('N', 'N', N2, hess, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO)
  !call sgeev('N', 'N', N2, hess, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO)
  !call DHSEQR('E','N',N2,ILO,IHI,hess,WR,WI,LDZ, WORK1, LWORK1, INFO)
  !call sgeev('N', 'N', N2, hess, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO)

  !SUBROUTINE SGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR,
  !   $                  LDVR, WORK, LWORK, INFO )
!write(*,*)WR
!write(*,*)WI
!write(*,*)info

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

end subroutine zrhqr1
  
