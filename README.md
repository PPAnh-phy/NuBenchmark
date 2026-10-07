# NuBenchmark

A Python pipeline for computing three-flavour neutrino oscillation probabilities in vacuum and constant-density matter, and comparing the fundamentally developed Hamiltonian calculations with [NuFast-LBL](https://github.com/PeterDenton/NuFast-LBL) and [NuOscProbExact](https://github.com/mbustama/NuOscProbExact).

## Calculation

For neutrinos, the flavour-basis Hamiltonian is

$$
H(E)=\frac{1}{2E}U\,\mathrm{diag}(0,\Delta m^2_{21},\Delta m^2_{31})U^\dagger+\mathrm{diag}(V_{CC},0,0).
$$

Here, $U$ is the PMNS matrix and $V_{CC}=\sqrt{2}G_F N_e$ is the charged-current matter potential. Vacuum calculations set the matter term to zero. For antineutrinos, the calculation uses $U^*$ and reverses the sign of the matter potential.

The evolution operator and probabilities are

$$
S(L,E)=\exp(-iHL),\qquad P_{\alpha\to\beta}=|S_{\beta\alpha}|^2.
$$
---
Two Hamiltonian methods are available for cross-checking:

- **`eigh`** (default mode) diagonalises the Hermitian Hamiltonian and constructs the evolution operator from its eigenvalues and eigenvectors.
- **`expm`** evaluates the matrix exponential with `scipy.linalg.expm`.

Both solve the same constant-density problem and return the same array layout. NuFast uses its specialised oscillation algorithm, while the three-flavour NuOscProbExact calculation uses an SU(3) formulation.

## Files

| File | Purpose |
|---|---|
| `PMNS.py` | PMNS matrix, parameter sets, and experiment settings |
| `hamiltonian.py` | Build the Hamiltonian and compute probabilities|
| `Benchmark.py` | Calls the three computators and calculates their differences |
| `Plot.py` | Visualise the results by plots |
| `pipeline.py` | Command-line interface and result saving |
| `nufast_bridge.cpp` | Connects the NuFast C++ calculation to Python |
| `build_nufast.py` | Builds the harmonised and native NuFast libraries |

## Installation

Install the Python dependencies in the environment used to run the pipeline:

```bash
python -m pip install numpy scipy matplotlib
```

Download the two external repositories of [NuFast-LBL](https://github.com/PeterDenton/NuFast-LBL) and [NuOscProbExact](https://github.com/mbustama/NuOscProbExact). The supplied build and import paths expect the following folders under the same parent directory:
- `benchmark_pipeline/`: the Python pipeline files and `nufast_bridge.cpp`.
- `NuFast-LBL-main/`: including `Benchmarks/src/NuFast_LBL.cpp` and `Benchmarks/include/`.
- `NuOscProbExact-main/`: including its `src/` directory.

Use a compatible source version containing the modules called by `Benchmark.py`. If you change the folder names, update the paths in `build_nufast.py` and `Benchmark.py`.

A C++ compiler available as `g++` is required for NuFast. On Windows, the MSYS2 UCRT64 compiler can be used; its `bin` directory must be available in `PATH`. Use a compiler architecture compatible with your Python installation.

From the pipeline directory, run:

```bash
python build_nufast.py
```

This builds both libraries:

| Platform | Harmonised library | Native library |
|---|---|---|
| Windows | `nufast_harmonised.dll` | `nufast_native.dll` |
| Linux | `libnufast_harmonised.so` | `libnufast_native.so` |
| macOS | `libnufast_harmonised.dylib` | `libnufast_native.dylib` |

Build the libraries on the platform where they will run. Rebuild after changing the bridge, its constants, or the NuFast source. The benchmark selects the appropriate library through `--mode`.

On Windows, `py` can replace `python` in these commands. On Linux or macOS, use `python3` if required by your installation.

## Benchmark modes

Both modes use the same **NuFIT 6.0 best-fit values** for the selected mass ordering, together with the same baseline, energy grid, density, and electron fraction.

| Setting | `harmonised` | `native` |
|---|---|---|
| Oscillation parameters | Shared reference values | Each package's input style |
| Natural-unit conversion | Hamiltonian conventions for all computators | Each package's own conventions |
| Matter coefficient | Hamiltonian coefficient for all calculators | Each package's own coefficient |

Harmonised mode compares the absolute differences under matching physical inputs and constants. Native mode also includes differences caused by the packages' numerical conversion and matter conventions. 

NuFast receives squared sines of mixing angles and a CP phase in radians. NuOscProbExact's native angle inputs are sines, obtained from the same squared-sine reference values. 

These different input sets are defined in `PMNS.py`.

## Experiment settings
The experimental parameters include:
- Baseline [km]
- Energy range [GeV] 
- Density [g/cm³] 
- Electron fraction 
- Default energy points 
- Selected flavour channel 

Change `PMNS_PARAMETERS` or `EXPERIMENTS` in `PMNS.py` to change the shared physics settings. The `alpha` and `beta` entries select the channel shown and saved by the benchmark.

## Running the benchmark

Run commands from the pipeline directory. For example:

```bash
python pipeline.py --mode harmonised --experiment JUNO --ordering NO --particle antinu --medium matter --hamiltonian-method eigh --show
```

```bash
python pipeline.py --mode native --experiment DUNE --ordering IO --particle nu --medium matter --hamiltonian-method expm --show
```

For a vacuum calculation, use `--medium vacuum`.

| Option | Values or meaning | Default |
|---|---|---|
| `--mode` | `harmonised`, `native` | `harmonised` |
| `--experiment` | `JUNO`, `DUNE` | `DUNE` |
| `--ordering` | `NO`, `IO` | `NO` |
| `--particle` | `nu`, `antinu` | `nu` |
| `--medium` | `vacuum`, `matter` | `matter` |
| `--hamiltonian-method` | `eigh`, `expm` | `eigh` |
| `--points` | Number of sampled energies | Experiment setting |
| `--n-newton` | NuFast Newton iteration setting | `2` |
| `--output` | Parent folder for results | `results` |
| `--show` | Display plots as well as saving them | Off |
| `--all` | Run every experiment/ordering/particle/medium combination | Off |

To run all 16 combinations for one mode and Hamiltonian method:

```bash
python pipeline.py --all --mode harmonised --hamiltonian-method eigh
```

Repeat with `--mode native` or `--hamiltonian-method expm` as needed.

## Results

Each run creates a new timestamped directory under the output folder. Repeating a calculation keeps earlier results. 

Each case produces a CSV file and a PNG plot. Filenames identify the experiment, ordering, particle, medium, mode, and Hamiltonian method, for example:

```text
JUNO_NO_antinu_matter_harmonised_eigh.csv
JUNO_NO_antinu_matter_harmonised_eigh.png
```

The CSV contains the energy, the selected-channel probability from each calculator, and three absolute differences:

- |NuFast − Hamiltonian|
- |NuOscProbExact − Hamiltonian|
- |NuFast − NuOscProbExact|

The upper plot shows probabilities. The lower plot shows the absolute differences on a logarithmic scale. Exactly zero differences cannot be displayed on a logarithmic axis.

The terminal summary also reports probability conservation, evaluated over all energies and initial flavours:

$$
\max_{E,\alpha}\left|\sum_\beta P_{\alpha\to\beta}(E)-1\right|.
$$

## Using the Hamiltonian independently
Place `hamiltonian.py` and `PMNS.py` in the same directory as your script.

Example script:
```python
import numpy as np

from PMNS import reference_parameters
from hamiltonian import oscillation_probabilities

# The NO/IO parameter set or the experimental parameters can be modified based on your interests.

# Select the NO/IO parameter set defined in PMNS.py.
parameters = reference_parameters("NO")

P = oscillation_probabilities(
    L=52.5,                              # Baseline [km]
    E=np.linspace(0.0018, 0.008, 1000),    # Energies [GeV]
    ordering="NO",                        # "NO" or "IO"
    parameters=parameters,
    rho=2.6,                              # Matter density [g/cm³]
    Ye=0.5,                               # Electron fraction
    matter=True,                         # False for vacuum
    antineutrino=True,                    # False for neutrinos
    method="eigh",                        # "eigh" or "expm"
)

# Indices: [energy, initial flavour, final flavour].
# Flavours: 0 = electron, 1 = muon, 2 = tau.
P_ee = P[:, 0, 0]       
P_mue = P[:, 1, 0]      

print(P[0])             # All nine probabilities at the first energy
```

The Hamiltonian calls `pmns_matrix()` from `PMNS.py` internally. The result has shape `(number_of_energies, 3, 3)`. A scalar energy returns shape `(1, 3, 3)`.

## References

- [NuFast-LBL source code](https://github.com/PeterDenton/NuFast-LBL)
- [NuFast-LBL paper](https://arxiv.org/abs/2405.02400)
- [NuOscProbExact source code](https://github.com/mbustama/NuOscProbExact)
- [NuOscProbExact paper](https://arxiv.org/abs/1904.12391)
- [NuFIT 6.0](https://arxiv.org/abs/2410.05380, https://www.nu-fit.org/)
- [Neutrino Formalism](https://arxiv.org/abs/1802.05781v2)

