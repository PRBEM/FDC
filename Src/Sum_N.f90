! ==============================================================================
! SUBROUTINE: sum_n
! ==============================================================================
!
! PURPOSE:
!   Top-level driver routine. For each row of the input grid (L, E, alpha),
!   computes the fully bounce-averaged quasi-linear diffusion coefficients
!   by summing over resonance harmonics:
!
!     D_aa, D_ap, D_pp  (pitch-angle, mixed, momentum)
!
!   and derives the mixed energy and energy-energy diffusion coefficients
!   using the relativistic energy-momentum relation:
!
!     D_aE = D_ap * (E+2)/(E+1)
!     D_EE = D_pp * (E+2)^2/(E+1)^2
!
!   Results are converted from s^-1 to day^-1 before output.
!
!   The full call chain is:
!     sum_n
!       -> D_ba               (bounce-average for one harmonic)
!          -> Integration_lambda
!             -> int_D_loc
!                -> D_local_n
!                   -> WeightFunction  (integrand at one theta point)
!
!   The formulation follows:
!     Albert, J. M. (2007), Simple approximations of quasi-linear diffusion
!     coefficients, J. Geophys. Res., 112, A12202, doi:10.1029/2007JA012551
!   which builds on the PADIE framework of:
!     Glauert, S. A., and R. B. Horne (2005), Calculation of pitch angle and
!     energy diffusion coefficients with the PADIE code, J. Geophys. Res.,
!     110, A04206, doi:10.1029/2004JA010851
!
!   Wave and plasma parameters used:
!     - Hiss wave amplitudes from Orlova, Spasojevic & Shprits [GRL, 2014]
!     - Plasmaspheric density model from Denton et al. [2004, 2006]
!     - Wave normal angle model from Thorne et al. [2013]
!     - Resonance order from Mourenas & Ripoll [2012] (if NresType='MR12')
!     - D_pp bounce-averaging approach from Orlova & Shprits [2011]
!
! INPUTS:
!   InputData2(local_rows2, columns) - Input grid, one row per (L, E, alpha):
!     column 1: L-shell [Earth radii, dimensionless]
!     column 2: Kinetic energy [keV]
!     column 3: Equatorial pitch angle [degrees]
!   local_rows2    - Number of rows in InputData2
!   columns        - Number of columns in InputData2
!   noutput_rows   - Number of rows in OutputData
!   noutput_columns- Number of columns in OutputData
!
! OUTPUTS:
!   OutputData(noutput_rows, noutput_columns):
!     column 1: L-shell [Earth radii]
!     column 2: Kinetic energy [keV]
!     column 3: Equatorial pitch angle [degrees]
!     column 4: Bounce-averaged D_aa [day^-1]
!     column 5: Bounce-averaged D_pp [day^-1]
!     column 6: Bounce-averaged D_ap [day^-1]
!
! CALLS:
!   resonanceorder  - Computes maximum resonance harmonic (Mourenas & Ripoll 2012)
!   D_ba            - Bounce-averaged diffusion coefficients for one harmonic
!
! MODULE DEPENDENCIES:
!   constants  - Physical constants (qe, me, c, RAD, ...),
!                and control parameters:
!                  NresType    : resonance order mode ('CONSTANT' or 'MR12')
!                  nres_min    : minimum resonance harmonic (module default for
!                                'CONSTANT' mode; set locally for 'MR12')
!                  nres_max    : maximum resonance harmonic (same)
!
! NOTES:
!   - Energy is converted from keV to normalised units internally:
!       E = E_keV * 1e3 * qe / (me * c^2)   (= gamma - 1, dimensionless)
!   - Pitch angle is converted from degrees to radians internally.
!   - All diffusion coefficients are multiplied by 24*3600 to convert s^-1
!     to day^-1 before being written to OutputData.
!   - D_aE and D_EE are computed internally but not written to OutputData.
!     If they are needed downstream, OutputData should be extended to 7+
!     columns and the assignments below uncommented.
!
! ==============================================================================

subroutine sum_n(InputData2, local_rows2, columns, OutputData, noutput_rows, noutput_columns)
use constants
implicit none

    integer, intent(in) :: local_rows2, columns
    integer, intent(in) :: noutput_rows, noutput_columns

    real*8, intent(in)  :: InputData2(local_rows2, columns)
    real*8, intent(out) :: OutputData(noutput_rows, noutput_columns)

    real*8  :: BavDaa, BavDap, BavDaE, BavDpp, BavDEE
    real*8  :: BavRate(3)
    real*8  :: L, E, alpha_eq

    integer :: i, resonmax, nres, I_nres

    do i = 1, local_rows2

        L        = InputData2(i,1)
        E        = InputData2(i,2) * 1.0e3 * qe / (me*c*c)
        alpha_eq = InputData2(i,3) * RAD

        BavDaa = 0.0
        BavDap = 0.0
        BavDaE = 0.0
        BavDpp = 0.0
        BavDEE = 0.0

        if (L <= L_max) then

            select case(NresType)

                case('CONSTANT')
                    ! nres_min and nres_max are used as set in module constants
                    ! (via the input file).
                    nres_min = nres_min
                    nres_max = nres_max

                case('MR12')
                    call resonanceorder(E, L, resonmax)
                    nres_min = -resonmax
                    nres_max =  resonmax

                case default
                    write(*,*) 'ERROR: Invalid NresType = ', trim(NresType)
                    write(*,*) 'Valid options are: CONSTANT, MR12'
                    stop

            end select

            do I_nres = nres_min, nres_max
                nres = I_nres
                call D_ba(alpha_eq, E, L, nres, BavRate)
                BavDaa = BavDaa + BavRate(1)
                BavDap = BavDap + BavRate(2)
                BavDpp = BavDpp + BavRate(3)
            end do

            ! Derive energy-related diffusion coefficients from momentum ones
            BavDaE = BavDap * (E + 2.0) / (E + 1.0)
            BavDEE = BavDpp * ((E + 2.0)**2.0) / ((E + 1.0)**2.0)

            ! Convert s^-1 to day^-1
            BavDaa = BavDaa * 24.0 * 3600.0
            BavDap = BavDap * 24.0 * 3600.0
            BavDpp = BavDpp * 24.0 * 3600.0
            BavDaE = BavDaE * 24.0 * 3600.0
            BavDEE = BavDEE * 24.0 * 3600.0

        else

            BavDaa = 0.0
            BavDap = 0.0
            BavDpp = 0.0
            BavDaE = 0.0
            BavDEE = 0.0

        end if

        OutputData(i,1) = L
        OutputData(i,2) = InputData2(i,2)
        OutputData(i,3) = InputData2(i,3)
        OutputData(i,4) = BavDaa
        OutputData(i,5) = BavDpp
        OutputData(i,6) = BavDap
        ! D_aE and D_EE are computed above but not currently written to output.
        ! To include them, extend OutputData to 8 columns and uncomment:
        ! OutputData(i,7) = BavDaE
        ! OutputData(i,8) = BavDEE

    end do

end subroutine sum_n