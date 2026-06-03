subroutine resonanceorder (E,L,resonmax)
use constants
implicit none
  real*8,intent(in)::E,L
  integer,intent(out)::resonmax
  real*8::Ne_local

  real*8::p,B_eq,omega_p,omega_e_eq

  call B_field(L,0,B_eq)
  call ElectronDensity(L,0,Ne_local)

  omega_e_eq=qe*B_eq/me
  omega_p=sqrt(Ne_local*qe*qe/(C0*me))
  p=sqrt((E+1.0)**2.0-1.0)
  resonmax=int(p*sqrt(omega_m_per*omega_e_eq*omega_p**2/omega_e_eq**3)*sin(85.0*RAD))+1
  return
endsubroutine resonanceorder
