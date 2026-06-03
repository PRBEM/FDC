## Subroutine `Bw_wave` 

**Purpose**: Calculates wave amplitude `Bw` (in Tesla) for different empirical wave models based on L-shell, magnetic latitude, Kp index, and MLT.

**Units**: Input L-shell and latitude (radians) → Output Bw in Tesla (pT internally, converted to T).

**Key feature**: Massive `SELECT CASE(BwType)` covering 15+ published chorus/hiss models from 2014-2018 papers.

### Interface

**Inputs** (from `constants` module):
```fortran
SUBROUTINE Bw_wave(L, lambda, Bw)
REAL8, INTENT(IN) :: L ! L-shell
REAL8, INTENT(IN) :: lambda ! Magnetic latitude [radians]
REAL*8, INTENT(OUT) :: Bw ! Wave amplitude [Tesla]
```

**Global dependencies** (set in `main`):
- `BwType`     : string selecting model
- `Kp`         : geomagnetic index  
- `MLT`        : magnetic local time [hours]
- `RAD`        : radians-to-degrees constant

### Workflow
1. Convert `lambda` from radians → degrees (`lambda_deg = lambda/RAD`)
2. `SELECT CASE(BwType)` → pick empirical model
3. Each case computes `Bw_Input` [pT] using:
   - Polynomial fits to log₁₀(Bw²)
   - Tabulated coefficients 
   - Latitude/MLT/L-shell clipping
   - Kp scaling factors
4. `Bw = Bw_Input * 1e-12` (pT → T)
5. Return

### Wave Models Available

| Model | Type | Dependencies | Key Features |
|-------|------|--------------|--------------|
| `TypicalDaysideChorus` | Chorus | Kp, latitude | Hui 2017 baseline |
| `TypicalNightsideChorus` | Chorus | Kp | Hui 2017 baseline |
| `TypicalHiss` | Hiss | Kp | Simple linear |
| `KS14DaysideHiss` | Hiss | Kp, L, lat | 8-coeff polynomial |
| `2016DaysideChorus` | Chorus | Kp, lat | Simple 3-term |
| `OA2015` | Chorus | Kp, lat, L, MLT | 16-coeff matrix, dayside/nightside |
| `NbbLBC/UBC` | Chorus | Kp, lat | Ni et al. tabulated |
| `2017Chorus` | Chorus | Kp(¹), L, MLT, lat | 5th-order Kp poly |
| `2018UBChorus/LBChorus` | Chorus | Kp(⁴), L, MLT, lat | Complex MLT/L bins |

¹ Most clip latitude >17-45°, L to [3.5,6]

### Error Handling
```fortran
CASE DEFAULT
WRITE(,) 'Invalid Wave Amplitude.'
STOP
```

**Invalid `BwType`** → immediate program termination.

### Output
- `Bw` [Tesla] = `Bw_Input` [pT] × 1e-12
- All models return **physical values** (no NaN/inf)
- **Units consistent** across models after conversion

