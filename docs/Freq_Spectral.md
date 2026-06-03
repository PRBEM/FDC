## Subroutine `Freq_Spectral`

**Purpose**: Sets wave frequency range parameters (`ω_m, Δω, ω_lc, ω_uc`) based on L-shell, latitude, and model selection.

**Key feature**: Converts input units (Hz → rad/s or ratios → absolute) and applies L/MLT/Kp-dependent empirical ranges.

**Output**: All frequencies in **rad/s**, ready for spectral density calculations.

### Interface

```fortran
SUBROUTINE Freq_Spectral(L, Lambda, omega_e_eq, omega_m, d_omega, omega_lc, omega_uc)
    REAL*8, INTENT(IN)  :: L           ! L-shell
    REAL*8, INTENT(IN)  :: Lambda      ! Magnetic latitude [radians]
    REAL*8, INTENT(IN)  :: omega_e_eq  ! Equatorial electron gyrofrequency [rad/s]
    REAL*8, INTENT(OUT) :: omega_m     ! Peak/mean frequency [rad/s]
    REAL*8, INTENT(OUT) :: d_omega     ! Frequency width [rad/s]
    REAL*8, INTENT(OUT) :: omega_lc    ! Lower cutoff [rad/s]
    REAL*8, INTENT(OUT) :: omega_uc    ! Upper cutoff [rad/s]
```
Global inputs: `FreqUnit`, `omega_m_per`, `d_omega_per`, `o_lc`, `o_uc`, `Kp`, `MLT`

### Workflow
1. Convert latitude: `lambda_deg = Lambda/rad`
2. `SELECT CASE(FreqUnit)` → determine scaling:
   - `'HZ'`: `*2π` (convert Hz → rad/s)
   - `'ratio'`: `*ω_e_eq` (normalize to gyrofrequency)  
   - **Preset models**: Override with hardcoded L/Kp/lat tables
3. **Special cases** (NbbLBC/UBC, 2016Hiss): L-binned or Kp×lat tables
4. **Chorus presets**: Simple ω_ce fractions (0.05-1.0 ω_ce)
5. Return all 4 frequencies in rad/s

### Error Handling
```fortran
CASE DEFAULT
WRITE(,) 'Invalid Frequency Unit.'
STOP
```


### Output Summary
| Parameter | Description | Units | Range |
|-----------|-------------|-------|-------|
| `omega_m` | Peak/mean frequency | rad/s | Model-dependent |
| `d_omega` | Bandwidth/width | rad/s | Model-dependent |
| `omega_lc` | Lower cutoff | rad/s | Typically 0.05-0.1ω_ce |
| `omega_uc` | Upper cutoff | rad/s | Typically 0.5-1.0ω_ce |

