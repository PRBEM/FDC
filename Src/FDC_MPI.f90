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
!   (D_aa, D_ap, D_pp, D_aE, D_EE) on a (L, E, alpha) grid for a specified
!   wave environment and plasma model.
!
!   The formulation follows:
!     Albert, J. M. (2007), Simple approximations of quasi-linear diffusion
!     coefficients, J. Geophys. Res., 112, A12202, doi:10.1029/2007JA012551
!   which builds on the PADIE framework of:
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
!   Line 44: DispersionType (COLD or HOT)
!   Line 46: DiffusionFormulation (LYONS or CUNNINGHAM)
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
!   - DispersionType is hardcoded to 'COLD' (cold plasma dispersion relation).
!     This is the only currently implemented option. When the hot plasma
!     (WHAMP) interface is implemented, DispersionType should be read from
!     the input file instead. See Module_Dispersion.f90 for the interface.
!   - The DIPOLE planet case in B_field does not apply the L^3 scaling,
!     making the field independent of L-shell. This may be intentional for
!     uniform-field test cases, but should be verified before using.
!
! ==============================================================================

program main
use constants
implicit none

    character(len=100) :: Input_file, filename2, char_tmp
    integer :: i, ios_diff_formulation
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
             '2018LBCHORUS', '2018UBCHORUS', '2026LBCHORUS', '2026UBCHORUS')
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
    DispersionType = trim(DispersionType)
    select case(DispersionType)
        case('COLD', 'HOT')
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