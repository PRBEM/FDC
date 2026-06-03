## Subroutine `ElectronDensity`

**Purpose**: Computes local electron density `Ne_local` [cm⁻³] using empirical plasmaspheric models based on L-shell and magnetic latitude.

**Key feature**: Sheeley01 family (10 variants) + classic models (CA92, Denton06) with latitudinal scaling.

**Output**: Physical electron density in cm⁻³, ready for plasma frequency/gyrofrequency calculations.

### Interface

```fortran
SUBROUTINE ElectronDensity(L, lambda, Ne_local)
    REAL*8, INTENT(IN)  :: L        ! L-shell
    REAL*8, INTENT(IN)  :: lambda   ! Magnetic latitude [radians]
    REAL*8, INTENT(OUT) :: Ne_local ! Electron density [cm⁻³]
```

Global dependencies:
- `DenModel`: density model selector
- `MLT`: magnetic local time [hours]
- `DenValue`: constant density value
- `DenLatIndex`: latitude scaling exponent
- `Kp`: used by Nbb model

### Workflow
1. `SELECT CASE(DenModel)` → pick empirical formula
2. Compute base density from L(-shell), MLT, Kp (model-specific)
3. **Apply universal latitude scaling**:
$$Ne\_local = Ne\_base \times (cos(\lambda))^{DenLatIndex}$$
4. Return density in cm⁻³

### Density Models Available

| Model | Reference | Dependencies | Formula Type |
|-------|-----------|--------------|--------------|
| `Sheely01`+`Div*` (10 variants) | Sheeley 2001 | L, MLT | L⁻³.⁵⁻⁴ × cos(MLT phase) |
| `CA92` | Carpenter & Anderson 1992 | L only | 10⁽³.⁹⁰⁴³⁻⁰.³¹⁴⁵ᴸ⁾ |
| `DENTON06` | Denton 2004/2006 | L only | 10⁽³.⁷⁸⁻⁰.³²⁴ᴸ⁺…⁾ |
| `CONSTANT` | User | `DenValue` | Fixed ×10⁶ cm⁻³ |
| `NbbLBUB` | Ni et al. | Kp | Kp=0:16.9e6, Kp≥2:2.88e6 |
| `XD_Saturn` | Gu 2012 | L, λ | Saturn-specific Gaussian |

### Error Handling
```fortran
CASE DEFAULT
WRITE(,) 'Invalid DENSITY MODEL.'
STOP
```

### Output
- `Ne_local`: **cm⁻³**, physically realistic plasmaspheric values
- **Typical range**: 10²-10⁶ cm⁻³ (L=2→8)
- **Ready for**: plasma frequency `f_pe = 9√Ne`, refractive index

