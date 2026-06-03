! Calculate B_local as a function of L and latitude

subroutine B_field(L,lambda,Blocal)
use constants,only:B0,PT
implicit none

  real*8,intent(in)::L,lambda
  real*8,intent(out)::Blocal

  real*8::f_value

  call f_function(lambda,f_value)
  select case(PT)
        case('JUPITER','EARTH','SATURN')
                Blocal=B0*f_value/(L**3.0)
        case('DIPOLE')
                Blocal=B0*f_value
        CASE DEFAULT
            WRITE(*,*) 'Invalid Planet.'
            STOP
    END SELECT
  

  return

end subroutine B_field
