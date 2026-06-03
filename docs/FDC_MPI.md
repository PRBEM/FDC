## Program `main` (FDC_MPI.f90)

This program is the main driver of the Full Diffusion Code (FDC). It reads a configuration file, initializes physical and numerical parameters, writes a log of the setup, and then calls an MPI-parallel routine to compute diffusion coefficients.

### Command-line interface

The program expects one command-line argument:

- `Input_file`: path to the text configuration file.

Example:
```bash
./FDC_MPI.x config.in
```

From this file, the code reads:

- The output file name (Output_file).
- All model and grid configuration parameters.

A log file is created automatically as:
- `log_file = <Output_file without last 4 characters>_log.txt`.


### Initialization workflow

1. Start CPU time measurement.
2. Read `Input_file` name from the command line.
3. Open `Input_file` and skip the first 11 header lines.
4. Read `Output_file` and construct `log_file`.
5. Open `log_file` and write:
   - A header banner,
   - The input and output file names,
   - A “Start initialization” marker.
6. Sequentially read blocks from the input file that define:
   - Planet and background magnetic field,
   - Phase-space grids (L-shell, energy, pitch angle),
   - Geomagnetic activity and local time,
   - Resonance-type settings,
   - Plasma density model,
   - Wave frequency distribution,
   - Numerical integration resolution,
   - Latitude range for integration,
   - Wave amplitude model,
   - Wave normal-angle model.
7. Echo key choices and basic values to `log_file`.
8. Close `Input_file` and write an “End initialization” marker to the log.

### Physical and model configuration (high level)

From the configuration file, the program selects:

- A planet or magnetic-field model and sets the reference field `B0`.
- Discrete grids in:
  - L-shell,
  - Particle energy,
  - Pitch angle.
- A geomagnetic activity level (via Kp) and magnetic local time.
- A resonance-order prescription (either fixed bounds or model-based).
- An electron density model (empirical, scaled, or constant).
- A wave frequency distribution (e.g., Gaussian or empirical presets).
- A wave amplitude model (typical, empirical, or constant amplitude).
- A wave normal-angle (θ) distribution model (e.g., constant or predefined).

The code uses these choices later (in other modules) to compute local and bounce-averaged diffusion coefficients.

### Numerical integration and parallel run

The configuration also specifies:

- Numbers of latitude points, wave normal-angle points, and normalization points.
- Latitude limits (minimum and maximum latitude) for integration.

From these, the program derives some internal resolution parameters and logs all counts and limits.

After successful initialization, the main program:

1. Writes a “Start calculation” marker to `log_file`.
2. Calls:
   ```fortran
   call mpi_parallel()
    ```
3. After the parallel calculation finishes, it records:
    - An “End calculation” marker,
    - The total CPU time elapsed.

### Error handling

At each configuration step, the program checks that the chosen option strings (e.g., planet, density model, frequency unit, wave amplitude type, theta type) match known keywords.

- If an unknown option is encountered, the code:
  - Writes an informative message (e.g., “Invalid Planet.” or similar),
  - Calls `STOP`, terminating the run.

This ensures the run fails fast when the configuration file contains unsupported or misspelled options.
