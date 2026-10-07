! ==============================================================================
! SUBROUTINE: ElectronDensity
! ==============================================================================
!
! PURPOSE:
!   Returns the local electron number density Ne at a given L-shell and
!   magnetic latitude, using the density model specified by DenModel in
!   module constants.
!
!   The base density (at the equator or averaged over latitude) is computed
!   first, then scaled by cos(lambda)^DenLatIndex to account for the
!   latitude dependence of the plasma density along the field line.
!
!   Supported models:
!     SHEELY01      - Sheeley et al. [2001], eq. 7: empirical plasmaspheric
!                     density as a function of L and MLT
!     SHEELYDIVx    - Sheeley et al. [2001] divided by scaling factor x,
!                     where x is 1.1, 1.2, ..., 1.9, 2, 2.5, 3, ..., 10
!     CA92          - Carpenter & Anderson [1992]: equatorial plasmaspheric
!                     density as a function of L only
!     DENTON06      - Denton et al. [2004, 2006]: equatorial plasmaspheric
!                     density as a function of L
!     DENTONMULx    - Denton et al. [2004, 2006] multiplied by x (x=2,3,4,5)
!     DENTONDIVx    - Denton et al. [2004, 2006] divided by x (x=2,3,4,5)
!     CONSTANT      - Uniform density DenValue [cm^-3] everywhere
!     NBBLBUB       - Kp-dependent fixed densities (Binbin Ni study)
!     XD_SATURN     - Saturn density model (Gu et al. 2012)
!
! INPUTS:
!   L        - Magnetic L-shell [Earth radii, dimensionless]
!   lambda   - Magnetic latitude [radians]
!
! OUTPUTS:
!   Ne_local - Local electron number density [m^-3]
!
! MODULE DEPENDENCIES:
!   constants - DenModel, DenValue, DenLatIndex, MLT, Kp, PI
!
! NOTES:
!   - All model formulas return density in m^-3. The Sheeley and Denton
!     formulas include a 1e6 factor to convert from cm^-3 to m^-3.
!   - The latitude dependence cos(lambda)^DenLatIndex is applied to all
!     models after the base density is computed, including CONSTANT.
!     DenLatIndex = 0 gives a latitude-independent density.
!   - The Denton formula contains the term (0.00127*F10.7 - 0.0635) with
!     the solar radio flux fixed at F10.7 = 78.
!
! REFERENCES:
!   Sheeley et al. [2001]
!   Carpenter & Anderson [1992]
!   Denton et al. [2004]
!   Denton et al. [2006]
!   Gu et al. [2012]
!
! ==============================================================================

subroutine ElectronDensity(L, lambda, Ne_local)
use constants
implicit none

    real*8, intent(in)  :: L, lambda
    real*8, intent(out) :: Ne_local

    real*8 :: Ne_base         ! model density before the latitude factor [m^-3]

    select case(DenModel)

        ! ------------------------------------------------------------------
        ! Sheeley et al. [2001], eq. 7: equatorial density from L and MLT,
        ! optionally divided by a fixed factor
        ! ------------------------------------------------------------------
        case('SHEELY01', 'SHEELYDIV1_1', 'SHEELYDIV1_2', 'SHEELYDIV1_3', &
             'SHEELYDIV1_4', 'SHEELYDIV1_5', 'SHEELYDIV1_6', 'SHEELYDIV1_7', &
             'SHEELYDIV1_8', 'SHEELYDIV1_9', 'SHEELYDIV2',   'SHEELYDIV2_5', &
             'SHEELYDIV3',   'SHEELYDIV4',   'SHEELYDIV5',   'SHEELYDIV6',   &
             'SHEELYDIV7',   'SHEELYDIV8',   'SHEELYDIV9',   'SHEELYDIV10')

            Ne_base = 1.24e8 * ((3.0/L)**4) &
                    + 3.6e7  * ((3.0/L)**3.5) &
                    * cos((MLT - (7.7*((3.0/L)**2) + 12.0)) * PI/12.0)

            select case(DenModel)
                case('SHEELYDIV1_1');  Ne_base = Ne_base / 1.1
                case('SHEELYDIV1_2');  Ne_base = Ne_base / 1.2
                case('SHEELYDIV1_3');  Ne_base = Ne_base / 1.3
                case('SHEELYDIV1_4');  Ne_base = Ne_base / 1.4
                case('SHEELYDIV1_5');  Ne_base = Ne_base / 1.5
                case('SHEELYDIV1_6');  Ne_base = Ne_base / 1.6
                case('SHEELYDIV1_7');  Ne_base = Ne_base / 1.7
                case('SHEELYDIV1_8');  Ne_base = Ne_base / 1.8
                case('SHEELYDIV1_9');  Ne_base = Ne_base / 1.9
                case('SHEELYDIV2');  Ne_base = Ne_base / 2.0
                case('SHEELYDIV2_5');  Ne_base = Ne_base / 2.5
                case('SHEELYDIV3');  Ne_base = Ne_base / 3.0
                case('SHEELYDIV4');  Ne_base = Ne_base / 4.0
                case('SHEELYDIV5');  Ne_base = Ne_base / 5.0
                case('SHEELYDIV6');  Ne_base = Ne_base / 6.0
                case('SHEELYDIV7');  Ne_base = Ne_base / 7.0
                case('SHEELYDIV8');  Ne_base = Ne_base / 8.0
                case('SHEELYDIV9');  Ne_base = Ne_base / 9.0
                case('SHEELYDIV10');  Ne_base = Ne_base / 10.0
            end select

        ! ------------------------------------------------------------------
        ! Carpenter & Anderson [1992]: L-dependent only
        ! ------------------------------------------------------------------
        case('CA92')
            Ne_base = 1e6 * 10.0**(-0.3145*L + 3.9043)

        ! ------------------------------------------------------------------
        ! Denton et al. [2004, 2006]: L-dependent equatorial density,
        ! optionally multiplied or divided by a fixed factor
        ! ------------------------------------------------------------------
        case('DENTON06', 'DENTONMUL2', 'DENTONMUL3', 'DENTONMUL4', 'DENTONMUL5', &
             'DENTONDIV2', 'DENTONDIV3', 'DENTONDIV4', 'DENTONDIV5')

            Ne_base = 1e6 * 10.0**(-0.324*L + 3.78 &
                      + (0.00127*78.0 - 0.0635)*exp((2.0-L)/1.5))

            select case(DenModel)
                case('DENTONMUL2');  Ne_base = Ne_base * 2.0
                case('DENTONMUL3');  Ne_base = Ne_base * 3.0
                case('DENTONMUL4');  Ne_base = Ne_base * 4.0
                case('DENTONMUL5');  Ne_base = Ne_base * 5.0
                case('DENTONDIV2');  Ne_base = Ne_base / 2.0
                case('DENTONDIV3');  Ne_base = Ne_base / 3.0
                case('DENTONDIV4');  Ne_base = Ne_base / 4.0
                case('DENTONDIV5');  Ne_base = Ne_base / 5.0
                ! DENTON06: unscaled
            end select

        ! ------------------------------------------------------------------
        ! Uniform density from the input file
        ! ------------------------------------------------------------------
        case('CONSTANT')
            Ne_base = DenValue * 1e6

        ! ------------------------------------------------------------------
        ! Kp-dependent fixed densities (Binbin Ni study)
        ! ------------------------------------------------------------------
        case('NBBLBUB')
            if (Kp == 0.0) then
                Ne_base = 16.9e6
            else if (Kp == 1.0) then
                Ne_base = 5.52e6
            else
                Ne_base = 2.88e6
            end if

        ! ------------------------------------------------------------------
        ! Saturn density model, Gu et al. [2012]
        ! ------------------------------------------------------------------
        case('XD_SATURN')
            Ne_base = 60e6 * exp(-1.0/3.0 * L**2/0.735**2 &
                      * (1.0 - cos(lambda)**6))

        case default
            write(*,*) 'ERROR: Invalid DenModel = ', trim(DenModel)
            stop

    end select

    ! Latitude dependence: Ne(lambda) = Ne_eq * cos(lambda)^DenLatIndex
    Ne_local = Ne_base * (cos(lambda))**DenLatIndex

    return

end subroutine ElectronDensity