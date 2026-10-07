! ==============================================================================
! PROGRAM: main
! ==============================================================================
!
! PURPOSE:
!   Entry point for the Full Diffusion Code (FDC). Reads all run parameters
!   from an input file, initialises module constants, and launches the MPI
!   parallel computation via mpi_parallel.
!
!   The program computes bounce-averaged quasi-linear diffusion coefficients
!   (D_aa, D_ap, D_pp) on a (L, E, alpha) grid for a specified
!   wave environment and plasma model.
!
!   The formulation follows:
!     Albert, J. M. (2007), Simple approximations of quasi-linear diffusion
!     coefficients, J. Geophys. Res., 112, A12202, doi:10.1029/2007JA012551
!     Glauert, S. A., and R. B. Horne (2005), Calculation of pitch angle and
!     energy diffusion coefficients with the PADIE code, J. Geophys. Res.,
!     110, A04206, doi:10.1029/2004JA010851
!
! INPUT FILE FORMAT:
!   The input file is read sequentially. Comment lines (starting with !)
!   are skipped by reading into char_tmp. The expected order is:
!
!   Line 12: Output filename
!   Line 14: Planet name (EARTH, JUPITER, SATURN, DIPOLE)
!   Line 16: Magnetic field strength B0 [nT] (used only for DIPOLE)
!   Line 18: Grid dimensions: NPA, NE, NL
!   Line 20: L range: L_min*10, L_max*10  (divided by 10 internally)
!   Line 22: Energy range: E_min, E_max [keV]
!   Line 24: Pitch angle range: PA_min, PA_max [degrees]
!   Line 26: Kp*10, MLT  (Kp divided by 10 internally)
!   Line 28: NresType, Nres_min, Nres_max
!   Line 30: DenModel, DenValue [cm^-3], DenLatIndex
!   Line 32: FreqDistribution
!   Line 34: FreqUnit, omega_m_per, d_omega_per, o_lc, o_uc
!   Line 36: nint, nint_X, nint_N
!   Line 38: lam_min, lam_max [degrees]
!   Line 40: BwType, Bw_Input [pT]
!   Line 42: ThetaType, Theta_Input(7)
!   Line 44: DispersionType (COLD, HOT or EMIC)
!   Line 46: DiffusionFormulation (LYONS or CUNNINGHAM) -- optional, see below
!   Then, only for DispersionType = HOT, the hot plasma species block:
!            k_min_exp, k_max_exp, n_species
!            followed by one line per species:
!            mass_number, density_fraction, T_par_keV, T_anis
!
!   The file is read POSITIONALLY: every parameter block consists of a
!   comment line followed by a value line, and a missing block shifts every
!   subsequent read. The only exceptions are the two blocks at the end:
!   DiffusionFormulation defaults to LYONS if absent (for backward
!   compatibility with older input files), and the hot plasma block is read
!   only for DispersionType = HOT, where it is mandatory.
!
! INPUT CONVENTIONS:
!   - L is given as L*10 (e.g. 30 = L=3.0) and divided by 10 on reading
!   - Kp is given as Kp*10 (e.g. 40 = Kp=4.0) and divided by 10 on reading
!   - Energy is in keV; converted to normalised units inside sum_n
!   - Pitch angle is in degrees; converted to radians inside sum_n
!   - nint_zero_res = 3*nint (finer grid for nres=0 resonance)
!   - nint_N_res_cone = 3*nint_N (finer grid near resonance cone)
!   - All model name strings are case-insensitive (converted to uppercase
!     internally via the to_upper function)
!
! OUTPUTS:
!   - Output file: Tecplot ASCII format with columns
!       L, E [keV], alpha [deg], Daa [day^-1], Dpp [day^-1], Dap [day^-1]
!   - Log file: run parameters and timing, named <output>_log.txt
!
! CALLS:
!   mpi_parallel - Distributes computation across MPI processes
!
! MODULE DEPENDENCIES:
!   constants - All run parameters are stored here after reading
!
! NOTES:
!   - B0 (equatorial surface field) for known planets is hardcoded:
!       EARTH:   2.997527836235779e-5 T
!       JUPITER: 4.2818e-4 T
!       SATURN:  0.2e-4 T
!       DIPOLE:  read from input file in nT, converted to T (* 1e-9)
!   - DispersionType selects the dispersion relation: COLD (cold plasma),
!     HOT (kinetic, via WHAMP; see Module_Dispersion.f90) or EMIC.
!   - The DIPOLE planet case in B_field deliberately omits the L^3 scaling:
!     the field read from the input file is then the equatorial field AT the
!     L-shell of interest, and L only labels the field line. Use this when
!     the field strength is prescribed directly, e.g. for code comparisons.
!   - LOCAL MODE: lam_min = lam_max = 0 switches off the bounce averaging and
!     returns the local equatorial coefficients (see D_ba). The log records
!     'LOCAL MODE' when it is active; any other lam_min = lam_max is rejected.
!
! COMPILE SWITCHES (all off by default; see build.sh):
!   -DPROFILE        per-rank timings for Integration_lambda, Roots_finder,
!                    WeightFunction and int_N_omega
!   -DDIAG           diagnostic output (accounting, N(omega) integrand,
!                    WHAMP surface dumps)
!   -DROOTS_POLY     cold resonant roots from the degree-10 polynomial
!                    instead of the scan in Roots_scan.f90
!   -DNO_THETA_SKIP  do not skip theta values that cannot resonate
!   -DINTN_ORIG      original N(omega) loop instead of int_N_cold_fast
!   -DOLD_BM         mirror field from the shortened mirror latitude
!   -DROOTCHECK      run both root finders and log differences
!   None of these change results except -DOLD_BM; the others are
!   performance or diagnostic switches.
!
! ==============================================================================


program main
use constants
implicit none

    character(len=100) :: Input_file, filename2, char_tmp
    integer :: i, ios_diff_formulation, ios_hot_params
    real*8  :: mag

    call cpu_time(time_start)
    call getarg(1, Input_file)

    ! -------------------------------------------------------------------------
    ! Open input file and skip comment header (first 11 lines)
    ! -------------------------------------------------------------------------
    open(unit=unit_input, file=Input_file, form='FORMATTED', status='OLD', &
         access='SEQUENTIAL', action='READ')
    do i = 1, 11
        read(unit_input,*) char_tmp
    end do

    ! -------------------------------------------------------------------------
    ! Output and log filenames
    ! -------------------------------------------------------------------------
    read(unit_input,*) filename2
    Output_file = trim(filename2)

    char_tmp = '_log.txt'
    log_file = Output_file(1:(len_trim(Output_file)-4)) // char_tmp

    open(unit=unit_log, file=log_file, form='FORMATTED', status='REPLACE', &
         access='SEQUENTIAL')
    write(unit_log,"('**********************************')")
    write(unit_log,"('Full Diffusion Code')")
    write(unit_log,"('**********************************')")
    write(unit_log,"(/)")
    write(unit_log,"('Input file name:' 4X,A30)") Input_file
    write(unit_log,"(/)")
    write(unit_log,"('**********Start initialization**********')")
    write(unit_log,"(/)")
    write(unit_log,"('Output file name:' 4X,A30)") Output_file
    write(unit_log,"('log file name:' 4X,A30)") log_file

    ! -------------------------------------------------------------------------
    ! Planet and magnetic field
    ! -------------------------------------------------------------------------
    read(unit_input,*) char_tmp
    read(unit_input,*) char_tmp
    PT = to_upper(trim(char_tmp))

    read(unit_input,*) char_tmp
    read(unit_input,*) mag

    select case(PT)
        case('JUPITER'); B0 = 4.2818e-4
        case('EARTH');   B0 = 2.997527836235779e-5
        case('SATURN');  B0 = 0.2e-4
        case('DIPOLE');  B0 = mag * 1e-9
        case default
            write(*,*) 'ERROR: Invalid Planet = ', trim(PT)
            stop
    end select
    write(unit_log,"('Planet:' 4X,A10)") PT

    ! -------------------------------------------------------------------------
    ! Grid dimensions
    ! -------------------------------------------------------------------------
    read(unit_input,*) char_tmp
    read(unit_input,*) NPA, NE, NL

    ! L range: input is L*10, divided by 10 here
    read(unit_input,*) char_tmp
    read(unit_input,*) L_min, L_max
    L_min = L_min / 10.0
    L_max = L_max / 10.0
    write(unit_log,"('Number of Lshell:' 4X, I4)") NL
    write(unit_log,"('Minimum of Lshell:' 4X, f5.2)") L_min
    write(unit_log,"('Maximum of Lshell:' 4X, f5.2)") L_max

    read(unit_input,*) char_tmp
    read(unit_input,*) E_min, E_max
    write(unit_log,"('Number of Energy:' 4X, I4)") NE
    write(unit_log,"('Minimum of Energy:' 4X, f5.2, 1X 'keV')") E_min
    write(unit_log,"('Maximum of Energy:' 4X, f10.1, 1X 'keV')") E_max

    read(unit_input,*) char_tmp
    read(unit_input,*) PA_min, PA_max
    write(unit_log,"('Number of Pitch Angle:'4X,I4)") NPA
    write(unit_log,"('Minimum of Pitch Angle:' 4X, f5.2, 1X 'Degree')") PA_min
    write(unit_log,"('Maximum of Pitch Angle:' 4X, f5.2, 1X 'Degree')") PA_max

    ! -------------------------------------------------------------------------
    ! Geomagnetic activity: Kp*10 in input, divided by 10 here
    ! -------------------------------------------------------------------------
    read(unit_input,*) char_tmp
    read(unit_input,*) Kp10, MLT
    write(unit_log,"('Kp10:',4X,f5.2)") Kp10
    Kp = Kp10 / 10.0
    write(unit_log,"('Magnetic Local Time:'4X,f5.2)") MLT

    ! -------------------------------------------------------------------------
    ! Resonance order
    ! -------------------------------------------------------------------------
    read(unit_input,*) char_tmp
    read(unit_input,*) NresType, Nres_min, Nres_max
    NresType = to_upper(trim(NresType))
    select case(NresType)
        case('CONSTANT')
            write(unit_log,"('resonance type:'4X,A10)") NresType
            write(unit_log,"('Minimum resonance order:'4X,I3)") Nres_min
            write(unit_log,"('Maximum resonance order:'4X,I3)") Nres_max
        case('MR12')
            write(unit_log,"('resonance type:'4X,A10)") NresType
            write(unit_log,"('Minimum and Maximum resonance order are " // &
                  "determined by equation 22 in MR 2012')")
        case default
            write(unit_log,"('Invalid resonance type:'4X,A10)") NresType
            stop
    end select

    ! -------------------------------------------------------------------------
    ! Electron density model
    ! -------------------------------------------------------------------------
    read(unit_input,*) char_tmp
    read(unit_input,*) DenModel, DenValue, DenLatIndex
    DenModel = to_upper(trim(DenModel))
    select case(DenModel)
        case('SHEELY01');     write(unit_log,"('electron density model:    Sheely et al. 2001')")
        case('SHEELYDIV1_1'); write(unit_log,"('electron density model:    Sheely div. by 1.1')")
        case('SHEELYDIV1_2'); write(unit_log,"('electron density model:    Sheely div. by 1.2')")
        case('SHEELYDIV1_3'); write(unit_log,"('electron density model:    Sheely div. by 1.3')")
        case('SHEELYDIV1_4'); write(unit_log,"('electron density model:    Sheely div. by 1.4')")
        case('SHEELYDIV1_5'); write(unit_log,"('electron density model:    Sheely div. by 1.5')")
        case('SHEELYDIV1_6'); write(unit_log,"('electron density model:    Sheely div. by 1.6')")
        case('SHEELYDIV1_7'); write(unit_log,"('electron density model:    Sheely div. by 1.7')")
        case('SHEELYDIV1_8'); write(unit_log,"('electron density model:    Sheely div. by 1.8')")
        case('SHEELYDIV1_9'); write(unit_log,"('electron density model:    Sheely div. by 1.9')")
        case('SHEELYDIV2');   write(unit_log,"('electron density model:    Sheely div. by 2')")
        case('SHEELYDIV2_5'); write(unit_log,"('electron density model:    Sheely div. by 2.5')")
        case('SHEELYDIV3');   write(unit_log,"('electron density model:    Sheely div. by 3')")
        case('SHEELYDIV4');   write(unit_log,"('electron density model:    Sheely div. by 4')")
        case('SHEELYDIV5');   write(unit_log,"('electron density model:    Sheely div. by 5')")
        case('SHEELYDIV6');   write(unit_log,"('electron density model:    Sheely div. by 6')")
        case('SHEELYDIV7');   write(unit_log,"('electron density model:    Sheely div. by 7')")
        case('SHEELYDIV8');   write(unit_log,"('electron density model:    Sheely div. by 8')")
        case('SHEELYDIV9');   write(unit_log,"('electron density model:    Sheely div. by 9')")
        case('SHEELYDIV10');  write(unit_log,"('electron density model:    Sheely div. by 10')")
        case('CA92');         write(unit_log,"('electron density model:    Carpenter and Anderson 1992')")
        case('DENTON06');     write(unit_log,"('electron density model:    Denton et al. 2004 2006')")
        case('DENTONMUL2');   write(unit_log,"('electron density model:    Denton et al. 2004 2006, multiplied by 2')")
        case('DENTONMUL3');   write(unit_log,"('electron density model:    Denton et al. 2004 2006, multiplied by 3')")
        case('DENTONMUL4');   write(unit_log,"('electron density model:    Denton et al. 2004 2006, multiplied by 4')")
        case('DENTONMUL5');   write(unit_log,"('electron density model:    Denton et al. 2004 2006, multiplied by 5')")
        case('DENTONDIV2');   write(unit_log,"('electron density model:    Denton et al. 2004 2006, divided by 2')")
        case('DENTONDIV3');   write(unit_log,"('electron density model:    Denton et al. 2004 2006, divided by 3')")
        case('DENTONDIV4');   write(unit_log,"('electron density model:    Denton et al. 2004 2006, divided by 4')")
        case('DENTONDIV5');   write(unit_log,"('electron density model:    Denton et al. 2004 2006, divided by 5')")
        case('CONSTANT');     write(unit_log,"('electron density:'4x,f8.4,1x'cm^-3')") DenValue
        case('NBBLBUB');      write(unit_log,"('Nbb PAPER')")
        case('XD_SATURN');    write(unit_log,"('see Gu et al. 2012 paper')")
        case default
            write(unit_log,"('Invalid density type'4x,A10)") DenModel
            stop
    end select
    write(unit_log,"('index of density latitudinal variation:'4X,f5.2)") DenLatIndex

    ! -------------------------------------------------------------------------
    ! Wave frequency distribution model
    ! -------------------------------------------------------------------------
    read(unit_input,*) char_tmp
    read(unit_input,*) FreqDistribution
    FreqDistribution = to_upper(trim(FreqDistribution))
    select case(FreqDistribution)
        case('GAUSSIAN', '2016HISS', '2017CHORUS', '2018UBCHORUS', &
             '2018LBCHORUS', '2026LBCHORUS', '2026UBCHORUS')
            write(unit_log,"('Frequency distribution type:'4x,A30)") FreqDistribution
        case default
            write(unit_log,"('Invalid frequency distribution type:'4X,A20)") FreqDistribution
            stop
    end select

    ! -------------------------------------------------------------------------
    ! Wave frequency band parameters
    ! -------------------------------------------------------------------------
    read(unit_input,*) char_tmp
    read(unit_input,*) FreqUnit, omega_m_per, d_omega_per, o_lc, o_uc
    FreqUnit = to_upper(trim(FreqUnit))
    select case(FreqUnit)
        case('HZ')
            write(unit_log,"('f_m=' f7.2,1X 'Hz,' 1X 'f_w=' f7.2,1X 'Hz,' 1X " // &
                  "'f_l=' f7.2,1X 'Hz,' 1X 'f_u=' f7.2,1X 'Hz')") &
                  omega_m_per, d_omega_per, o_lc, o_uc
        case('RATIO')
            write(unit_log,"('omega_m/Omega_e=' f7.2,1X ',' 1X " // &
                  "'delta_omega/Omega_e=' f7.2,1X ',' 1X " // &
                  "'omega_lc/Omega_e=' f7.2,1X ',' 1X " // &
                  "'omega_uc/Omega_e=' f7.2,1X)") omega_m_per, d_omega_per, o_lc, o_uc
        case('NBBLBC', 'NBBUBC', '2016HISS', '2017CHORUS', &
             '2018LBCHORUS', '2018UBCHORUS', '2026LBCHORUS', '2026UBCHORUS', &
             'EMIC_HBAND', 'EMIC_HEBAND', 'EMIC_OBAND')
            write(unit_log,"('Frequency unit type:'4X,A20)") FreqUnit
        case default
            write(unit_log,"('Invalid Frequency Unit.'4x,A10)") FreqUnit
            stop
    end select

    ! -------------------------------------------------------------------------
    ! Integration resolution
    ! -------------------------------------------------------------------------
    read(unit_input,*) char_tmp
    read(unit_input,*) nint, nint_X, nint_N
    nint_zero_res    = 3 * nint
    nint_N_res_cone  = 3 * nint_N
    write(unit_log,"('Points of latitude:'4X,I4)") nint
    write(unit_log,"('Points of normal angle distribution:' 4X,I4)") nint_X
    write(unit_log,"('Points of normalized factor:' 4X,I4)") nint_N

    read(unit_input,*) char_tmp
    read(unit_input,*) lam_min, lam_max
    write(unit_log,"('Minimum of latitude:' 4X, f5.2, 1X 'Degree')") lam_min
    write(unit_log,"('Maximum of latitude:' 4X, f5.1, 1X 'Degree')") lam_max

    ! Local mode: lam_min = lam_max gives local D at that latitude, with no
    ! bounce averaging. Only the equator is supported: at lambda /= 0 the
    ! local D refers to the local pitch angle, but the output is labelled
    ! with the equatorial one.
    local_mode = (abs(lam_max - lam_min) < 1.0d-12)
    if (local_mode) then
        if (abs(lam_min) > 1.0d-12) then
            write(unit_log,"('ERROR: local mode (lam_min = lam_max) ' // &
                  'requires lam_min = 0 (equator only).')")
            write(*,*) 'ERROR: local mode (lam_min = lam_max) requires lam_min = 0.'
            stop
        end if
        write(unit_log,"('LOCAL MODE: equatorial local D, no bounce averaging')")
    end if

    ! -------------------------------------------------------------------------
    ! Wave amplitude model
    ! -------------------------------------------------------------------------
    read(unit_input,*) char_tmp
    read(unit_input,*) BwType, Bw_Input
    BwType = to_upper(trim(BwType))
    select case(BwType)
        case('TYPICALDAYSIDECHORUS', 'TYPICALNIGHTSIDECHORUS', 'TYPICALHISS', &
             'KS14DAYSIDEHISS', 'KS14NIGHTSIDEHISS', '2016DAYSIDECHORUS', &
             '2016NIGHTSIDECHORUS', '2016NIGHTSIDECHORUS2')
            write(unit_log,"('wave amplitude type:'4x,A30)") BwType
        case('OA2015', 'VITALIICHORUS', '2016HISS', '2017CHORUS', &
             '2017UPPERBANDCHORUS', '2017LOWERBANDCHORUS', &
             '2018UBCHORUS', '2018LBCHORUS', 'MIDLATCHORUS')
            write(unit_log,"('wave amplitude type:'4x,A30)") BwType
        case('NBBLBC', 'NBBUBC', '2018CHORUSLAT1545', '2026LBCHORUS', '2026UBCHORUS')
            write(unit_log,"('wave amplitude type:'4x,A30)") BwType
        case('CONSTANT')
            write(unit_log,"('wave amplitude type:'4x,A30)") BwType
            write(unit_log,"('wave amplitude:' 4X, f8.3, 1x 'pT')") Bw_Input
        case default
            write(unit_log,"('Invalid Wave Amplitude Type:' 4x, A30)") BwType
            stop
    end select

    ! -------------------------------------------------------------------------
    ! Wave normal angle model
    ! -------------------------------------------------------------------------
    read(unit_input,*) char_tmp
    read(unit_input,*) ThetaType, Theta_Input
    ThetaType = to_upper(trim(ThetaType))
    select case(ThetaType)
        case('CONSTANT')
            write(unit_log,"('normal angle model type:'4x,A10)") ThetaType
            write(unit_log,"('theta_m=' f4.1,1X 'deg,' 1X, 'theta_w=' f4.1," // &
                  "1X 'deg,' 1X, 'theta_lc=' f4.1,1X 'deg,' 1X 'theta_uc=' f4.1,1X 'deg.')") &
                  Theta_Input
        case('OBLIQUEHISS', 'OA2015', 'NBBLBC', 'NBBUBC', 'VITALIICHORUS')
            write(unit_log,"('normal angle model type:'4x, A10)") ThetaType
        case default
            write(unit_log,"('Invalid Theta Type:' 4x, A10)")
            stop
    end select

    ! -------------------------------------------------------------------------
    ! Dispersion type
    ! -------------------------------------------------------------------------
    read(unit_input,*) char_tmp
    read(unit_input,*) DispersionType
    DispersionType = to_upper(trim(DispersionType))
    select case(DispersionType)
        case('COLD', 'HOT', 'EMIC')
            if (DispersionType == 'EMIC') then
                write(unit_log,"('WARNING: EMIC dispersion is work in progress, ' // &
                      'unvalidated and incomplete. Results should not be used ' // &
                      'without checking the implementation.')")
                write(*,*) 'WARNING: EMIC dispersion is work in progress and unvalidated.'
            end if
            write(unit_log,"('Dispersion type:'4x,A10)") DispersionType
        case default
            write(unit_log,"('Invalid Dispersion Type:' 4x, A10)") DispersionType
            stop
    end select

    ! -------------------------------------------------------------------------
    ! Diffusion coefficient formulation: 'LYONS' (default) or 'CUNNINGHAM'
    !
    ! This is read last and is OPTIONAL for backward compatibility: input
    ! files written before this option existed will hit end-of-file on the
    ! comment-line read below, which we detect via IOSTAT and use to fall
    ! back to 'LYONS' (the code's original, unchanged behaviour) without
    ! aborting. Input files that DO include this option must supply both
    ! the comment line and the keyword line, exactly like every other
    ! parameter block above.
    ! -------------------------------------------------------------------------
    read(unit_input, *, iostat=ios_diff_formulation) char_tmp
    if (ios_diff_formulation /= 0) then
        DiffusionFormulation = 'LYONS'
        write(unit_log,"('Diffusion formulation:'4x,A30," // &
              "' (defaulted -- not present in input file)')") DiffusionFormulation
    else
        read(unit_input, *, iostat=ios_diff_formulation) DiffusionFormulation
        if (ios_diff_formulation /= 0) then
            DiffusionFormulation = 'LYONS'
            write(unit_log,"('Diffusion formulation:'4x,A30," // &
                  "' (defaulted -- malformed/missing value)')") DiffusionFormulation
        else
            DiffusionFormulation = to_upper(trim(DiffusionFormulation))
            select case(DiffusionFormulation)
                case('LYONS', 'CUNNINGHAM')
                    write(unit_log,"('Diffusion formulation:'4x,A30)") DiffusionFormulation
                case default
                    write(unit_log,"('Invalid Diffusion Formulation:' 4x, A30)") DiffusionFormulation
                    stop
            end select
        end if
    end if

    ! Hot plasma (WHAMP) species parameters: required when DispersionType='HOT'.
    if (trim(DispersionType) == 'HOT') then
        read(unit_input, *, iostat=ios_hot_params) char_tmp
        if (ios_hot_params /= 0) then
            write(unit_log,"('ERROR: DispersionType=HOT requires a Hot plasma " // &
                  "species block in the input file, but none was found.')")
            write(*,*) 'ERROR: DispersionType=HOT requires a Hot plasma species ' // &
                        'block in the input file, but none was found.'
            stop
        end if

        read(unit_input, *, iostat=ios_hot_params) &
             k_min_exp_val_input, k_max_exp_val_input, n_hot_species
        if (ios_hot_params /= 0) then
            write(unit_log,"('ERROR: Malformed Hot plasma species header line.')")
            write(*,*) 'ERROR: Malformed Hot plasma species header line ' // &
                        '(expected: k_min_exp, k_max_exp, n_species).'
            stop
        end if

        if (n_hot_species < 1 .or. n_hot_species > max_hot_species) then
            write(unit_log,"('ERROR: n_hot_species=' I4 ' out of valid range [1,' I3 '].')") &
                  n_hot_species, max_hot_species
            write(*,*) 'ERROR: n_hot_species out of valid range.'
            stop
        end if

        write(unit_log,"('k_min_exp=' f6.2, ' k_max_exp=' f6.2, ' n_species=' I3)") &
              k_min_exp_val_input, k_max_exp_val_input, n_hot_species

        do i = 1, n_hot_species
            read(unit_input, *, iostat=ios_hot_params) char_tmp   ! comment line
            if (ios_hot_params == 0) then
                read(unit_input, *, iostat=ios_hot_params) &
                     hot_species_mass(i), hot_species_density_frac(i), &
                     hot_species_T_par_keV(i), hot_species_T_anis(i)
            end if
            if (ios_hot_params /= 0) then
                write(unit_log,"('ERROR: Malformed Hot plasma species block at species' I3)") i
                write(*,*) 'ERROR: Malformed Hot plasma species block at species', i
                stop
            end if
            write(unit_log,"('  Species' I2 ': mass=' f6.2 ' frac=' f8.4 " // &
                  "' T_par=' f8.3 'keV T_anis=' f6.2)") &
                  i, hot_species_mass(i), hot_species_density_frac(i), &
                  hot_species_T_par_keV(i), hot_species_T_anis(i)
        end do
    end if

    close(unit_input)
    write(unit_log,"('**********End initialization**********')")
    write(unit_log,"(/)")
    write(unit_log,"('**********Start calculation**********')")
    write(unit_log,"(/)")

    ! -------------------------------------------------------------------------
    ! Launch MPI parallel computation
    ! -------------------------------------------------------------------------
    call mpi_parallel()

    call cpu_time(time_finish)
    write(unit_log,"('**********End calculation**********')")
    write(unit_log,"(/)")
    write(unit_log, "('Calculation Time = ',f8.3,' seconds.')") time_finish - time_start
    write(unit_log,"(/)")

contains

    ! --------------------------------------------------------------------------
    ! FUNCTION: to_upper
    ! Converts a string to uppercase. Used to make model name matching in
    ! select case blocks case-insensitive.
    ! --------------------------------------------------------------------------
    function to_upper(str) result(out)
        character(len=*), intent(in) :: str
        character(len=len(str))      :: out
        integer :: i
        do i = 1, len(str)
            select case(str(i:i))
                case('a':'z')
                    out(i:i) = achar(iachar(str(i:i)) - 32)
                case default
                    out(i:i) = str(i:i)
            end select
        end do
    end function to_upper

end program main