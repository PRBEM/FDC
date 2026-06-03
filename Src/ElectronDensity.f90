!!This subroutine is to calculate the electron density in different L and Latitude utilizing different density models

subroutine ElectronDensity(L,lambda,Ne_local)
use constants
implicit none

  real*8,intent(in)::L,lambda
  real*8,intent(out)::Ne_local

SELECT CASE(DenModel)
     CASE('SHEELY01') 
         Ne_local=1.24e8*((3.0/L)**4.0)+3.6e7*((3.0/L)**3.5)*cos((MLT-(7.7*((3.0/L)**2.0)+12.0))*PI/12.0)   ! Sheeley et al. 2001 Equation 7, related with L and MLT
     CASE('SHEELYDIV1_1') 
         Ne_local=1.24e8*((3.0/L)**4.0)+3.6e7*((3.0/L)**3.5)*cos((MLT-(7.7*((3.0/L)**2.0)+12.0))*PI/12.0)/1.1   ! Sheeley et al. 2001 Equation 7, related with L and MLT
     CASE('SHEELYDIV1_2') 
         Ne_local=1.24e8*((3.0/L)**4.0)+3.6e7*((3.0/L)**3.5)*cos((MLT-(7.7*((3.0/L)**2.0)+12.0))*PI/12.0)/1.2   ! Sheeley et al. 2001 Equation 7, related with L and MLT
     CASE('SHEELYDIV1_3') 
         Ne_local=1.24e8*((3.0/L)**4.0)+3.6e7*((3.0/L)**3.5)*cos((MLT-(7.7*((3.0/L)**2.0)+12.0))*PI/12.0)/1.3   ! Sheeley et al. 2001 Equation 7, related with L and MLT
     CASE('SHEELYDIV1_4') 
         Ne_local=1.24e8*((3.0/L)**4.0)+3.6e7*((3.0/L)**3.5)*cos((MLT-(7.7*((3.0/L)**2.0)+12.0))*PI/12.0)/1.4   ! Sheeley et al. 2001 Equation 7, related with L and MLT
     CASE('SHEELYDIV1_5') 
         Ne_local=1.24e8*((3.0/L)**4.0)+3.6e7*((3.0/L)**3.5)*cos((MLT-(7.7*((3.0/L)**2.0)+12.0))*PI/12.0)/1.5   ! Sheeley et al. 2001 Equation 7, related with L and MLT
     CASE('SHEELYDIV1_6') 
         Ne_local=1.24e8*((3.0/L)**4.0)+3.6e7*((3.0/L)**3.5)*cos((MLT-(7.7*((3.0/L)**2.0)+12.0))*PI/12.0)/1.6   ! Sheeley et al. 2001 Equation 7, related with L and MLT
     CASE('SHEELYDIV1_7') 
         Ne_local=1.24e8*((3.0/L)**4.0)+3.6e7*((3.0/L)**3.5)*cos((MLT-(7.7*((3.0/L)**2.0)+12.0))*PI/12.0)/1.7   ! Sheeley et al. 2001 Equation 7, related with L and MLT
     CASE('SHEELYDIV1_8')     
         Ne_local=1.24e8*((3.0/L)**4.0)+3.6e7*((3.0/L)**3.5)*cos((MLT-(7.7*((3.0/L)**2.0)+12.0))*PI/12.0)/1.8   ! Sheeley et al. 2001 Equation 7, related with L and MLT  
     CASE('SHEELYDIV1_9') 
         Ne_local=1.24e8*((3.0/L)**4.0)+3.6e7*((3.0/L)**3.5)*cos((MLT-(7.7*((3.0/L)**2.0)+12.0))*PI/12.0)/1.9   ! Sheeley et al. 2001 Equation 7, related with L and MLT
     CASE('SHEELYDIV2') 
         Ne_local=(1.24e8*((3.0/L)**4.0)+3.6e7*((3.0/L)**3.5)*cos((MLT-(7.7*((3.0/L)**2.0)+12.0))*PI/12.0))*0.5   ! Sheeley et al. 2001 Equation 7, related with L and MLT
     CASE('SHEELYDIV2_5') 
         Ne_local=(1.24e8*((3.0/L)**4.0)+3.6e7*((3.0/L)**3.5)*cos((MLT-(7.7*((3.0/L)**2.0)+12.0))*PI/12.0))/2.5   ! Sheeley et al. 2001 Equation 7, related with L and MLT
     CASE('SHEELYDIV3') 
         Ne_local=(1.24e8*((3.0/L)**4.0)+3.6e7*((3.0/L)**3.5)*cos((MLT-(7.7*((3.0/L)**2.0)+12.0))*PI/12.0))/3.0   ! Sheeley et al. 2001 Equation 7, related with L and MLT
     CASE('SHEELYDIV4') 
         Ne_local=(1.24e8*((3.0/L)**4.0)+3.6e7*((3.0/L)**3.5)*cos((MLT-(7.7*((3.0/L)**2.0)+12.0))*PI/12.0))/4.0   ! Sheeley et al. 2001 Equation 7, related with L and MLT
     CASE('SHEELYDIV5') 
         Ne_local=(1.24e8*((3.0/L)**4.0)+3.6e7*((3.0/L)**3.5)*cos((MLT-(7.7*((3.0/L)**2.0)+12.0))*PI/12.0))*0.2   ! Sheeley et al. 2001 Equation 7, related with L and MLT
     CASE('SHEELYDIV6') 
         Ne_local=(1.24e8*((3.0/L)**4.0)+3.6e7*((3.0/L)**3.5)*cos((MLT-(7.7*((3.0/L)**2.0)+12.0))*PI/12.0))/6.0   ! Sheeley et al. 2001 Equation 7, related with L and MLT
     CASE('SHEELYDIV7') 
         Ne_local=(1.24e8*((3.0/L)**4.0)+3.6e7*((3.0/L)**3.5)*cos((MLT-(7.7*((3.0/L)**2.0)+12.0))*PI/12.0))/7.0   ! Sheeley et al. 2001 Equation 7, related with L and MLT
     CASE('SHEELYDIV8') 
         Ne_local=(1.24e8*((3.0/L)**4.0)+3.6e7*((3.0/L)**3.5)*cos((MLT-(7.7*((3.0/L)**2.0)+12.0))*PI/12.0))/8.0   ! Sheeley et al. 2001 Equation 7, related with L and MLT
     CASE('SHEELYDIV9') 
         Ne_local=(1.24e8*((3.0/L)**4.0)+3.6e7*((3.0/L)**3.5)*cos((MLT-(7.7*((3.0/L)**2.0)+12.0))*PI/12.0))/9.0   ! Sheeley et al. 2001 Equation 7, related with L and MLT
     CASE('SHEELYDIV1O')
         Ne_local=(1.24e8*((3.0/L)**4.0)+3.6e7*((3.0/L)**3.5)*cos((MLT-(7.7*((3.0/L)**2.0)+12.0))*PI/12.0))*0.1   ! Sheeley et al. 2001 Equation 7, related with L and MLT
     CASE('CA92') 
         Ne_local=1e6*10.0**(-0.3145*L+3.9043)  !Carpenter and Anderson 1992 Plasmaspheric Density Model, only related with L
     CASE('DENTON06') 
         Ne_local=1e6*10.0**(-0.324*L+3.78+(0.00127*78.0-0.0635)*exp((2.0-L)/1.5));  !  Denton et al. 2004,2006 Density Model, only related with L
     CASE('DENTONMUL2')
         Ne_local=1e6*10.0**(-0.324*L+3.78+(0.00127*78.0-0.0635)*exp((2.0-L)/1.5))*2  !  Denton et al. 2004,2006 Density Model, related with L and Latitude
     CASE('DENTONMUL3')
         Ne_local=1e6*10.0**(-0.324*L+3.78+(0.00127*78.0-0.0635)*exp((2.0-L)/1.5))*3  !  Denton et al. 2004,2006 Density Model, related with L and Latitude
     CASE('DENTONMUL4')
         Ne_local=1e6*10.0**(-0.324*L+3.78+(0.00127*78.0-0.0635)*exp((2.0-L)/1.5))*4  !  Denton et al. 2004,2006 Density Model, related with L and Latitude
     CASE('DENTONMUL5')
         Ne_local=1e6*10.0**(-0.324*L+3.78+(0.00127*78.0-0.0635)*exp((2.0-L)/1.5))*5  !  Denton et al. 2004,2006 Density Model, related with L and Latitude
     CASE('DENTONDIV2')
         Ne_local=1e6*10.0**(-0.324*L+3.78+(0.00127*78.0-0.0635)*exp((2.0-L)/1.5))/2  !  Denton et al. 2004,2006 Density Model, related with L and Latitude
     CASE('DENTONDIV3')
         Ne_local=1e6*10.0**(-0.324*L+ 3.78+(0.00127*78.0-0.0635)*exp((2.0-L)/1.5))/3  !  Denton et al. 2004,2006 Density Model, related with L and Latitude
     CASE('DENTONDIV4')
         Ne_local=1e6*10.0**(-0.324*L+ 3.78+(0.00127*78.0-0.0635)*exp((2.0-L)/1.5))/4  !  Denton et al. 2004,2006 Density Model, related with L and Latitude
     CASE('DENTONDIV5')
         Ne_local=1e6*10.0**(-0.324*L+ 3.78+(0.00127*78.0-0.0635)*exp((2.0-L)/1.5))/5  !  Denton et al. 2004,2006 Density Model, related with L and Latitude
     CASE('CONSTANT') 
         Ne_local=DenValue*1e6
     CASE('NBBLBUB')
        if(kp.eq.0)then
                Ne_local=16.9e6
        else if(kp.eq.1)then
                Ne_local=5.52e6
        else
                Ne_local=2.88e6
        endif
    CASE('XD_SATURN')
        Ne_local=60e6*exp(-1/3.0*L**2/0.735**2*(1.0-cos(lambda)**6))
    CASE DEFAULT
        WRITE(*,*) 'Invalid DENSITY MODEL.'
        STOP
    END SELECT
 

Ne_local=Ne_local*(cos(lambda))**DenLatIndex

  

  return

end subroutine ElectronDensity
