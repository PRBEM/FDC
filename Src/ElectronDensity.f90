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
!     NBBLBUB       - Kp-dependent fixed densities (Saturn NBB study)
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
!     models after the base density is computed. DenLatIndex = 0 gives
!     a latitude-independent density.
!   - BUG FIXED: the original code had 'SHEELYDIV1O' (letter O) instead of
!     'SHEELYDIV10' (number zero) in the case statement, which would have
!     caused a silent fall-through to the DEFAULT error for that model.
!
! REFERENCES:
!   Sheeley et al. [2001] - J. Geophys. Res., 106, 25631-25641
!   Carpenter & Anderson [1992] - J. Geophys. Res., 97, 1097-1108
!   Denton et al. [2004] - J. Geophys. Res., 109, A09215
!   Denton et al. [2006] - J. Geophys. Res., 111, A04213
!
! ==============================================================================

subroutine ElectronDensity(L, lambda, Ne_local)
use constants
implicit none

    real*8, intent(in)  :: L, lambda
    real*8, intent(out) :: Ne_local

    ! Base Sheeley density before scaling (used by all SHEELY variants)
    real*8 :: Ne_sheeley_base
    real*8 :: sheeley_scale

    ! -------------------------------------------------------------------------
    ! Base Sheeley [2001] equatorial density formula (eq. 7)
    ! Depends on L and MLT; units m^-3
    ! -------------------------------------------------------------------------
    Ne_sheeley_base = 1.24e8 * ((3.0/L)**4.0) &
                    + 3.6e7  * ((3.0/L)**3.5) &
                    * cos((MLT - (7.7*((3.0/L)**2.0) + 12.0)) * PI/12.0)

    select case(DenModel)

        ! --- Sheeley variants: base density divided by a scaling factor ---
        case('SHEELY01');      sheeley_scale = 1.0
        case('SHEELYDIV1_1');  sheeley_scale = 1.1
        case('SHEELYDIV1_2');  sheeley_scale = 1.2
        case('SHEELYDIV1_3');  sheeley_scale = 1.3
        case('SHEELYDIV1_4');  sheeley_scale = 1.4
        case('SHEELYDIV1_5');  sheeley_scale = 1.5
        case('SHEELYDIV1_6');  sheeley_scale = 1.6
        case('SHEELYDIV1_7');  sheeley_scale = 1.7
        case('SHEELYDIV1_8');  sheeley_scale = 1.8
        case('SHEELYDIV1_9');  sheeley_scale = 1.9
        case('SHEELYDIV2');    sheeley_scale = 2.0
        case('SHEELYDIV2_5');  sheeley_scale = 2.5
        case('SHEELYDIV3');    sheeley_scale = 3.0
        case('SHEELYDIV4');    sheeley_scale = 4.0
        case('SHEELYDIV5');    sheeley_scale = 5.0
        case('SHEELYDIV6');    sheeley_scale = 6.0
        case('SHEELYDIV7');    sheeley_scale = 7.0
        case('SHEELYDIV8');    sheeley_scale = 8.0
        case('SHEELYDIV9');    sheeley_scale = 9.0
        case('SHEELYDIV10');   sheeley_scale = 10.0

        ! --- Other density models ---
        case('CA92')
            ! Carpenter & Anderson [1992]: L-dependent only [m^-3]
            Ne_local = 1e6 * 10.0**(-0.3145*L + 3.9043)

        case('DENTON06')
            ! Denton et al. [2004, 2006]: L-dependent equatorial density [m^-3]
            Ne_local = 1e6 * 10.0**(-0.324*L + 3.78 &
                       + (0.00127*78.0 - 0.0635)*exp((2.0-L)/1.5))

        case('DENTONMUL2')
            Ne_local = 1e6 * 10.0**(-0.324*L + 3.78 &
                       + (0.00127*78.0 - 0.0635)*exp((2.0-L)/1.5)) * 2.0

        case('DENTONMUL3')
            Ne_local = 1e6 * 10.0**(-0.324*L + 3.78 &
                       + (0.00127*78.0 - 0.0635)*exp((2.0-L)/1.5)) * 3.0

        case('DENTONMUL4')
            Ne_local = 1e6 * 10.0**(-0.324*L + 3.78 &
                       + (0.00127*78.0 - 0.0635)*exp((2.0-L)/1.5)) * 4.0

        case('DENTONMUL5')
            Ne_local = 1e6 * 10.0**(-0.324*L + 3.78 &
                       + (0.00127*78.0 - 0.0635)*exp((2.0-L)/1.5)) * 5.0

        case('DENTONDIV2')
            Ne_local = 1e6 * 10.0**(-0.324*L + 3.78 &
                       + (0.00127*78.0 - 0.0635)*exp((2.0-L)/1.5)) / 2.0

        case('DENTONDIV3')
            Ne_local = 1e6 * 10.0**(-0.324*L + 3.78 &
                       + (0.00127*78.0 - 0.0635)*exp((2.0-L)/1.5)) / 3.0

        case('DENTONDIV4')
            Ne_local = 1e6 * 10.0**(-0.324*L + 3.78 &
                       + (0.00127*78.0 - 0.0635)*exp((2.0-L)/1.5)) / 4.0

        case('DENTONDIV5')
            Ne_local = 1e6 * 10.0**(-0.324*L + 3.78 &
                       + (0.00127*78.0 - 0.0635)*exp((2.0-L)/1.5)) / 5.0

        case('CONSTANT')
            ! Uniform density from input file [m^-3]
            Ne_local = DenValue * 1e6

        case('NBBLBUB')
            ! Kp-dependent fixed densities [m^-3]
            if (Kp == 0.0) then
                Ne_local = 16.9e6
            else if (Kp == 1.0) then
                Ne_local = 5.52e6
            else
                Ne_local = 2.88e6
            end if

        case('XD_SATURN')
            ! Saturn density model, Gu et al. [2012]
            Ne_local = 60e6 * exp(-1.0/3.0 * L**2.0/0.735**2.0 &
                       * (1.0 - cos(lambda)**6.0))

        case default
            write(*,*) 'ERROR: Invalid DenModel = ', trim(DenModel)
            stop

    end select

    ! -------------------------------------------------------------------------
    ! Apply Sheeley scaling if a Sheeley variant was selected
    ! -------------------------------------------------------------------------
    select case(DenModel)
        case('SHEELY01', 'SHEELYDIV1_1', 'SHEELYDIV1_2', 'SHEELYDIV1_3', &
             'SHEELYDIV1_4', 'SHEELYDIV1_5', 'SHEELYDIV1_6', 'SHEELYDIV1_7', &
             'SHEELYDIV1_8', 'SHEELYDIV1_9', 'SHEELYDIV2',   'SHEELYDIV2_5', &
             'SHEELYDIV3',   'SHEELYDIV4',   'SHEELYDIV5',   'SHEELYDIV6',   &
             'SHEELYDIV7',   'SHEELYDIV8',   'SHEELYDIV9',   'SHEELYDIV10')
            Ne_local = Ne_sheeley_base / sheeley_scale
    end select

    ! -------------------------------------------------------------------------
    ! Apply latitude dependence: Ne(lambda) = Ne_eq * cos(lambda)^DenLatIndex
    ! DenLatIndex = 0 gives latitude-independent density.
    ! -------------------------------------------------------------------------
    Ne_local = Ne_local * (cos(lambda))**DenLatIndex

    return

end subroutine ElectronDensity