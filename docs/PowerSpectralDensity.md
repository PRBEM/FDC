## Subroutine `PowerSpectralDensity`

**Purpose**: Computes wave power spectral density `PSD` as a function of frequency for different empirical models.

**Key feature**: Normalized PSD distributions matching total integrated wave power across `[ω_lc, ω_uc]`.

**Units**: Frequency-normalized PSD (integrates to 1 over frequency range).

### Interface

```fortran
SUBROUTINE PowerSpectralDensity(L, Lambda, omega_e_eq, Freq, PSD)
    REAL*8, INTENT(IN)  :: L           ! L-shell
    REAL*8, INTENT(IN)  :: Lambda      ! Magnetic latitude [radians]
    REAL*8, INTENT(IN)  :: omega_e_eq  ! Equatorial electron gyrofrequency
    REAL*8, INTENT(IN)  :: Freq        ! Wave frequency [rad/s]
    REAL*8, INTENT(OUT) :: PSD         ! Power spectral density
```

Global dependencies (from `constants`/`variables`):
- `FreqDistribution`: model selector
- `omega_m`, `d_omega`, `omega_lc`, `omega_uc`: frequency parameters
- `MLT`: magnetic local time
- `rad`: radians-to-degrees constant

### Workflow
1. Convert latitude: `lambda_deg = Lambda/rad`
2. `SELECT CASE(FreqDistribution)` → pick model
3. Compute `PSD` based on model:
   - **Gaussian**: Normalized Gaussian centered at `ω_m`
   - **2016Hiss**: Kappa distribution with L-binned parameters
   - **2017/2018 Chorus**: 7th-order polynomials in `ω/ω_ce`
4. **Normalize** so ∫PSD dω = 1 over `[ω_lc, ω_uc]`
5. Apply minimum floor where needed

### Frequency Distribution Models

| Model | Type | Dependencies | Normalization |
|-------|------|--------------|---------------|
| `GAUSSIAN` | Analytical | `ω_m, d_ω, ω_lc, ω_uc` | Error function |
| `2016Hiss` | Kappa | L-shell bins | `Bw_norm`, explicit |
| `2017Chorus` | Polynomial | MLT (day/night), lat<7° | 8th-order integral |
| `2018UBChorus` | Polynomial | MLT (day/night) | 8th-order integral |
| `2018LBChorus` | Polynomial | MLT, lat≤10° bins | 8th-order integral |

### Safeguards
```fortran
CASE DEFAULT
WRITE(,) 'Invalid Frequency Distribution'
STOP
```

- **Invalid model** → program crash
- **2018 models**: `PSD ≥ 1e-20` (numerical floor)
- **Hiss**: `PSD=0` outside `[ω_lc, ω_uc]`

### Output
- `PSD`: **dimensionless**, normalized so ∫PSD df = 1
- **Ready for quasi-linear diffusion**: multiplies Bw² later
- **No units conversion** needed

**Dependencies**: `errorfunction()` external routine
