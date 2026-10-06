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
!                  nres_min    : minimum resonance harmonic ('CONSTANT' mode)
!                  nres_max    : maximum resonance harmonic ('CONSTANT' mode)
!
! NOTES:
!   - In local mode (lam_min = lam_max = 0) the same columns hold the LOCAL
!     equatorial coefficients; D_ba returns them without bounce averaging.
!   - HOT/EMIC: the WHAMP dispersion curves are built here, once per lambda
!     grid point, and passed down so that every harmonic reuses them. They
!     are rebuilt on every call, i.e. for every (L, E, alpha) point, although
!     they depend only on (L, lambda): the same curves are recomputed for all
!     NE energies and all NPA pitch angles of a given L. Reusing them across
!     a whole (L, alpha) column would cut the cost of a HOT run by close to a
!     factor NE. Note that lambda_m, and therefore the lambda grid, does
!     depend on alpha.
!   - NresType = 'MR12' overwrites nres_min and nres_max locally, not in
!     module constants, so the values read from the input file survive.
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
use dispersion
implicit none

    integer, intent(in) :: local_rows2, columns
    integer, intent(in) :: noutput_rows, noutput_columns

    real*8, intent(in)  :: InputData2(local_rows2, columns)
    real*8, intent(out) :: OutputData(noutput_rows, noutput_columns)

    real*8  :: BavDaa, BavDap, BavDaE, BavDpp, BavDEE
    real*8  :: BavRate(3)
    real*8  :: L, E, alpha_eq

    integer :: i, resonmax, nres, I_nres
    integer :: nres_lo, nres_hi          ! harmonic range for this grid point

    type(WHAMPCurve), allocatable :: whamp_curves_shared(:)
    real*8, allocatable :: lambda_grid_shared(:)
    integer :: n_lambda_shared
    real*8 :: lambda_m_shared, lambda_min_shared, delta_lambda_shared
    real*8 :: fce_khz_s, fpe_khz_s, B_local_s, Ne_local_s
    real*8 :: Theta_function_s(7)
    integer :: jj

    do i = 1, local_rows2

        L        = InputData2(i,1)
        E        = InputData2(i,2) * 1.0e3 * qe / (me*c*c)
        alpha_eq = InputData2(i,3) * RAD

        if (DispersionType == 'HOT' .or. DispersionType == 'EMIC') then
            ! lambda grid for the shared curves: the same grid D_ba will use
            call Lambda_mirror_new(alpha_eq, lambda_m_shared)
            lambda_min_shared = lam_min * RAD
            if (lambda_m_shared/RAD >= lam_max) then
                lambda_m_shared = lam_max * RAD
            else if (lambda_m_shared/RAD <= lam_min) then
                lambda_m_shared = lam_min * RAD
            end if

            n_lambda_shared = nint + 1
            delta_lambda_shared = (lambda_m_shared - lambda_min_shared) / real(nint)

            if (allocated(whamp_curves_shared)) then
                do jj = 1, size(whamp_curves_shared)
                    deallocate(whamp_curves_shared(jj)%theta, whamp_curves_shared(jj)%k, &
                            whamp_curves_shared(jj)%k_whamp, whamp_curves_shared(jj)%omega, &
                            whamp_curves_shared(jj)%dw_dk, whamp_curves_shared(jj)%valid, &
                            whamp_curves_shared(jj)%EFL)
                end do
                deallocate(whamp_curves_shared)
            end if
            if (allocated(lambda_grid_shared)) deallocate(lambda_grid_shared)

            allocate(lambda_grid_shared(n_lambda_shared))
            allocate(whamp_curves_shared(n_lambda_shared))

            do jj = 1, n_lambda_shared
                lambda_grid_shared(jj) = lambda_min_shared + delta_lambda_shared * real(jj-1)
                call B_field(L, lambda_grid_shared(jj), B_local_s)
                call ElectronDensity(L, lambda_grid_shared(jj), Ne_local_s)
                fce_khz_s = qe * B_local_s / me / (2.0d0 * PI * 1.0d3)
                fpe_khz_s = sqrt(Ne_local_s * qe * qe / (C0 * me)) / (2.0d0 * PI * 1.0d3)
                call Theta_wave(L, lambda_grid_shared(jj), Theta_function_s)
                call get_whamp_dispersion(fce_khz_s, fpe_khz_s, &
                                        Theta_function_s(5)*PI/180.0d0, &
                                        Theta_function_s(6)*PI/180.0d0, &
                                        whamp_curves_shared(jj))
            end do

        else
            ! COLD case: allocate trivially so the array is always in a valid state
            if (allocated(whamp_curves_shared)) then
                do jj = 1, size(whamp_curves_shared)
                    if (allocated(whamp_curves_shared(jj)%theta)) then
                        deallocate(whamp_curves_shared(jj)%theta, whamp_curves_shared(jj)%k, &
                                whamp_curves_shared(jj)%k_whamp, whamp_curves_shared(jj)%omega, &
                                whamp_curves_shared(jj)%dw_dk, whamp_curves_shared(jj)%valid, &
                                whamp_curves_shared(jj)%EFL)
                    end if
                end do
                deallocate(whamp_curves_shared)
            end if
            if (allocated(lambda_grid_shared)) deallocate(lambda_grid_shared)

            n_lambda_shared = 1
            allocate(whamp_curves_shared(1))
            allocate(lambda_grid_shared(1))
            lambda_grid_shared(1) = 0.0d0
            whamp_curves_shared(1)%n_k     = 0
            whamp_curves_shared(1)%n_theta = 0
        end if

        BavDaa = 0.0
        BavDap = 0.0
        BavDaE = 0.0
        BavDpp = 0.0
        BavDEE = 0.0

        if (L <= L_max) then

            select case(NresType)

                case('CONSTANT')
                    ! range as read from the input file
                    nres_lo = nres_min
                    nres_hi = nres_max

                case('MR12')
                    ! range estimated per grid point (Mourenas & Ripoll 2012)
                    call resonanceorder(E, L, resonmax)
                    nres_lo = -resonmax
                    nres_hi =  resonmax

                case default
                    write(*,*) 'ERROR: Invalid NresType = ', trim(NresType)
                    write(*,*) 'Valid options are: CONSTANT, MR12'
                    stop

            end select

            do I_nres = nres_lo, nres_hi
                nres = I_nres
                call D_ba(alpha_eq, E, L, nres, BavRate, whamp_curves_shared, lambda_grid_shared, n_lambda_shared)
                BavDaa = BavDaa + BavRate(1)
                BavDap = BavDap + BavRate(2)
                BavDpp = BavDpp + BavRate(3)
            end do

            ! Derive energy-related diffusion coefficients from momentum ones
            BavDaE = BavDap * (E + 2.0) / (E + 1.0)
            BavDEE = BavDpp * ((E + 2.0)**2) / ((E + 1.0)**2)

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