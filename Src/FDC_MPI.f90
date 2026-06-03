program main
use constants
implicit none
character(len=100)::Input_file,filename1,filename2,char_tmp
integer::i
real*8::mag
call cpu_time(time_start)



call getarg(1,Input_file)

OPEN(UNIT=unit_input,FILE=Input_file,FORM='FORMATTED',STATUS='OLD',ACCESS='SEQUENTIAL',ACTION='READ')
do i=1, 11
    read(unit_input,*)char_tmp
enddo

    read(unit_input,*)filename2
    Output_file=trim(filename2)


char_tmp='_log.txt'
log_file=Output_file(1:(LEN_TRIM(Output_file)-4))//char_tmp
OPEN(UNIT=unit_log,FILE=log_file,FORM='FORMATTED',STATUS='REPLACE',ACCESS='SEQUENTIAL')
WRITE(unit_log,"('**********************************')")
WRITE(unit_log,"('Full Diffusion Code')")
WRITE(unit_log,"('**********************************')")
WRITE(unit_log,"(/)")
WRITE(unit_log,"('Input file name:' 4X,A30)")Input_file
WRITE(unit_log,"(/)")
WRITE(unit_log,"('**********Start initialization**********')")
WRITE(unit_log,"(/)")
WRITE(unit_log,"('Output file name:' 4X,A30)")Output_file
WRITE(unit_log,"('log file name:' 4X,A30)")log_file

    read(unit_input,*)char_tmp
    read(unit_input,*)char_tmp
    PT=to_upper(trim(char_tmp))

        read(unit_input,*)char_tmp
        read(unit_input,*)mag
    select case(PT)
        case('JUPITER')
            B0=4.2818e-4   !for Jupiter
        case('EARTH')
            B0=2.997527836235779e-5 !for Earth
        case('SATURN')
            B0=0.2e-4  !for Saturn
        case('DIPOLE')
            B0=mag*1e-9
        CASE DEFAULT
            WRITE(*,*) 'Invalid Planet.'
            STOP
    END SELECT
WRITE(unit_log,"('Planet:' 4X,A10)")PT

    read(unit_input,*)char_tmp
    read(unit_input,*)NPA,NE,NL


    read(unit_input,*)char_tmp
    read(unit_input,*)L_min,L_max
        L_min=L_min/10.0
        L_max=L_max/10.0
WRITE(unit_log,"('Number of Lshell:' 4X, I4)")NL
WRITE(unit_log,"('Minimum of Lshell:' 4X, f5.2)")L_min
WRITE(unit_log,"('Maximum of Lshell:' 4X, f5.2)")L_max
    
        read(unit_input,*)char_tmp
    read(unit_input,*)E_min,E_max
WRITE(unit_log,"('Number of Energy:' 4X, I4)")NE
WRITE(unit_log,"('Minimum of Energy:' 4X, f5.2, 1X 'keV')")E_min
WRITE(unit_log,"('Maximum of Energy:' 4X, f10.1, 1X 'keV')")E_max

    read(unit_input,*)char_tmp
    read(unit_input,*)PA_min,PA_max
WRITE(unit_log,"('Number of Pitch Angle:'4X,I4)")NPA
WRITE(unit_log,"('Minimum of Pitch Angle:' 4X, f5.2, 1X 'Degree')")PA_min
WRITE(unit_log,"('Maximum of Pitch Angle:' 4X, f5.2, 1X 'Degree')")PA_max

    read(unit_input,*)char_tmp
    read(unit_input,*)Kp10,MLT
WRITE(unit_log,"('Kp10:',4X,f5.2)")Kp10
Kp=Kp10/10.0
WRITE(unit_log,"('Magnetic Local Time:'4X,f5.2)")MLT

    read(unit_input,*)char_tmp
    read(unit_input,*)NresType,Nres_min,Nres_max
    NresType=to_upper(trim(NresType))
   select case(NresType)
        case('CONSTANT')
                WRITE(unit_log,"('resonance type:'4X,A10)")NresType
                WRITE(unit_log,"('Minimum resonance order:'4X,I3)")Nres_min
                WRITE(unit_log,"('Maximum resonance order:'4X,I3)")Nres_max
        case('MR12')
                WRITE(unit_log,"('resonance type:'4X,A10)")NresType
                WRITE(unit_log,"('Minimum and Maximum resonance order are determined by equation 22 in MR 2012')")
        CASE DEFAULT
                WRITE(unit_log,"('Invalid resonance type:'4X,A10)")NresType
                STOP
    END SELECT
    read(unit_input,*)char_tmp
    read(unit_input,*)DenModel,DenValue,DenLatINdex
        DenModel=to_upper(trim(DenModel))
        SELECT CASE(DenModel)
        CASE('SHEELY01') 
                WRITE(unit_log,"('electron density model:    Sheely et al. 2001')")
        CASE('SHEELYDIV1_1') 
                WRITE(unit_log,"('electron density model:    Sheely div. by 1.1')")
        CASE('SHEELYDIV1_2') 
                WRITE(unit_log,"('electron density model:    Sheely div. by 1.2')")
        CASE('SHEELYDIV1_3') 
                WRITE(unit_log,"('electron density model:    Sheely div. by 1.3')")
        CASE('SHEELYDIV1_4') 
                WRITE(unit_log,"('electron density model:    Sheely div. by 1.4')")
        CASE('SHEELYDIV1_5') 
                WRITE(unit_log,"('electron density model:    Sheely div. by 1.5')")
        CASE('SHEELYDIV1_6') 
                WRITE(unit_log,"('electron density model:    Sheely div. by 1.6')")
        CASE('SHEELYDIV1_7') 
                WRITE(unit_log,"('electron density model:    Sheely div. by 1.7')")
        CASE('SHEELYDIV1_8') 
                WRITE(unit_log,"('electron density model:    Sheely div. by 1.8')")
        CASE('SHEELYDIV1_9') 
                WRITE(unit_log,"('electron density model:    Sheely div. by 1.9')")
        CASE('SHEELYDIV2') 
                WRITE(unit_log,"('electron density model:    Sheely div. by 2')")
        CASE('SHEELYDIV2_5')
                WRITE(unit_log,"('electron density model:    Sheely div. by 2.5')")
        CASE('SHEELYDIV3') 
                WRITE(unit_log,"('electron density model:    Sheely div. by 3')")
        CASE('SHEELYDIV4') 
                WRITE(unit_log,"('electron density model:    Sheely div. by 4')")
        CASE('SHEELYDIV5') 
                WRITE(unit_log,"('electron density model:    Sheely div. by 5')")
        CASE('SHEELYDIV6') 
                WRITE(unit_log,"('electron density model:    Sheely div. by 6')")
        CASE('SHEELYDIV7') 
                WRITE(unit_log,"('electron density model:    Sheely div. by 7')")
        CASE('SHEELYDIV8') 
                WRITE(unit_log,"('electron density model:    Sheely div. by 8')")
        CASE('SHEELYDIV9') 
                WRITE(unit_log,"('electron density model:    Sheely div. by 9')")
        CASE('SHEELYDIV10') 
                WRITE(unit_log,"('electron density model:    Sheely div. by 10')")
        CASE('CA92') 
                WRITE(unit_log,"('electron density model:    Carpenter and Anderson 1992')")
        CASE('DENTON06') 
                WRITE(unit_log,"('electron density model:    Denton et al. 2004 2006')")
        CASE('DENTONMUL2')
                WRITE(unit_log,"('electron density model:    Denton et al. 2004 2006, multiplied by 2')")
        CASE('DENTONMUL3')
                WRITE(unit_log,"('electron density model:    Denton et al. 2004 2006, multiplied by 3')")
        CASE('DENTONMUL4')
                WRITE(unit_log,"('electron density model:    Denton et al. 2004 2006, multiplied by 4')")
        CASE('DENTONMUL5')
                WRITE(unit_log,"('electron density model:    Denton et al. 2004 2006, multiplied by 5')")
        CASE('DENTONDIV2')
                WRITE(unit_log,"('electron density model:    Denton et al. 2004 2006, divided by 2')")
        CASE('DENTONDIV3')
                WRITE(unit_log,"('electron density model:    Denton et al. 2004 2006, divided by 3')")
        CASE('DENTONDIV4')
                WRITE(unit_log,"('electron density model:    Denton et al. 2004 2006, divided by 4')")
        CASE('DENTONDIV5')
                WRITE(unit_log,"('electron density model:    Denton et al. 2004 2006, divided by 5')")
        CASE('CONSTANT') 
                WRITE(unit_log,"('electron density:'4x,f8.4,1x'cm^-3')")DenValue
        CASE('NBBLBUB')
                WRITE(unit_log,"('Nbb PAPER')")
        CASE('XD_SATURN')
                write(unit_log,"('see Gu et al. 2012 paper')")
        CASE DEFAULT
                WRITE(unit_log,"('Invalid density type'4x,A10)")DenModel
                STOP
        END SELECT
        WRITE(unit_log,"('index of density latitudinal variation:'4X,f5.2)")DenLatINdex

    read(unit_input,*)char_tmp
    read(unit_input,*)FreqDistribution
    FreqDistribution=to_upper(trim(FreqDistribution))
        select case(FreqDistribution)
        case('GAUSSIAN', '2016HISS', '2017CHORUS', '2018UBCHORUS', '2018LBCHORUS', '2026LBCHORUS', '2026UBCHORUS')
                WRITE(unit_log,"('Frequency distribution type:'4x,A30)")FreqDistribution
        case default
                WRITE(unit_log,"('Invalid frequency distribution type:'4X,A20)")FreqDistribution
        stop
        end select 


    read(unit_input,*)char_tmp
    read(unit_input,*)FreqUnit,omega_m_per,d_omega_per,o_lc,o_uc
    FreqUnit=to_upper(trim(FreqUnit))
        SELECT CASE(FreqUnit)
        CASE('HZ') 
                WRITE(unit_log,"('f_m=' f7.2,1X 'Hz,' 1X 'f_w=' f7.2,1X 'Hz,' 1X 'f_l=' f7.2,1X 'Hz,' 1X 'f_u=' f7.2,1X 'Hz')")omega_m_per, d_omega_per, o_lc, o_uc
        CASE('RATIO') 
                WRITE(unit_log,"('omega_m/Omega_e=' f7.2,1X ',' 1X 'delta_omega/Omega_e=' f7.2,1X ',' 1X 'omega_lc/Omega_e=' f7.2,1X ',' 1X 'omega_uc/Omega_e=' f7.2,1X)")omega_m_per, d_omega_per,  o_lc, o_uc
        CASE('NBBLBC','NBBUBC', '2016HISS', '2017CHORUS', '2018LBCHORUS', '2026LBCHORUS', '2026UBCHORUS')
                WRITE(unit_log,"('Frequency unit type:'4X,A20)")FreqUnit
        CASE DEFAULT
                WRITE(unit_log,"('Invalid Frequency Unit.'4x,A10)") FreqUnit
                STOP
        END SELECT


    read(unit_input,*)char_tmp
    read(unit_input,*)nint,nint_X,nint_N
    nint_zero_res=3*nint
    nint_N_res_cone=3*nint_N
        WRITE(unit_log,"('Points of latitude:'4X,I4)")nint
        WRITE(unit_log,"('Points of normal angle distribution:' 4X,I4)")nint_X
        WRITE(unit_log,"('Points of normalized factor:' 4X,I4)")nint_N

    read(unit_input,*)char_tmp
    read(unit_input,*)lam_min,lam_max
        WRITE(unit_log,"('Minimum of latitude:' 4X, f5.2, 1X 'Degree')")lam_min
        WRITE(unit_log,"('Maximum of latitude:' 4X, f5.1, 1X 'Degree')")lam_max

    read(unit_input,*)char_tmp
    read(unit_input,*)BwType,Bw_Input
    BwType=(to_upper(trim(BwType)))
        
        SELECT CASE(BwType)
        CASE('TYPICALDAYSIDECHORUS','TYPICALNIGHTSIDECHORUS','TYPICALHISS','KS14DAYSIDEHISS','KS14NIGHTSIDEHISS','2016DAYSIDECHORUS','2016NIGHTSIDECHORUS','2016NIGHTSIDECHORUS2') 
                WRITE(unit_log,"('wave amplitude type:'4x,A30)")BwType
        CASE('OA2015', 'VITALIICHORUS', '2016HISS','2017CHORUS', '2017UPPERBANDCHORUS', '2017LOWERBANDCHORUS', '2018UBCHORUS', '2018LBCHORUS', 'MIDLATCHORUS') 
        WRITE(unit_log,"('wave amplitude type:'4x,A30)")BwType
        CASE('NBBLBC', 'NBBUBC', '2018CHORUSLAT1545', '2026LBCHORUS', '2026UBCHORUS')
                WRITE(unit_log,"('wave amplitude type:'4x,A30)")BwType
        CASE('CONSTANT')
                WRITE(unit_log,"('wave amplitude type:'4x,A30)")BwType
                WRITE(unit_log,"('wave amplitude:' 4X, f8.3, 1x 'pT')")Bw_Input
        CASE DEFAULT
                WRITE(unit_log,"('Invalid Wave Amplitude Type:' 4x, A30)") BwType
                STOP
        END SELECT

    read(unit_input,*)char_tmp
    read(unit_input,*)ThetaType,Theta_Input
    ThetaType=to_upper(trim(ThetaType))
        SELECT CASE(ThetaType)
        CASE('CONSTANT') 
                WRITE(unit_log,"('normal angle model type:'4x,A10)")ThetaType
                WRITE(unit_log,"('theta_m=' f4.1,1X 'deg,' 1X, 'theta_w=' f4.1,1X 'deg,' 1X, 'theta_lc=' f4.1,1X 'deg,' 1X 'theta_uc=' f4.1,1X 'deg.')")Theta_Input
        CASE('OBLIQUEHISS', 'OA2015', 'NBBLBC','NBBUBC', 'VITALIICHORUS') 
                WRITE(unit_log,"('normal angle model type:'4x, A10)")ThetaType
        CASE DEFAULT
                WRITE(unit_log,"('Invalid Theta Type:' 4x, A10)") 
                STOP
        END SELECT

close(unit_input)
WRITE(unit_log,"('**********End initialization**********')")
WRITE(unit_log,"(/)")
WRITE(unit_log,"('**********Start calculation**********')")
WRITE(unit_log,"(/)")

call mpi_parallel()

call cpu_time(time_finish)
WRITE(unit_log,"('**********End calculation**********')")
WRITE(unit_log,"(/)")
WRITE(unit_log, "('Calculation Time = ',f8.3,' seconds.')")time_finish-time_start
WRITE(unit_log,"(/)")

contains

function to_upper(str) result(out)
    character(len=*), intent(in) :: str
    character(len=len(str)) :: out
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
