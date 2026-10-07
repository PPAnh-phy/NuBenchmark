from pathlib import Path
import ctypes
import os
import shutil
import sys
import numpy as np
from PMNS import (
    experiment,
    nufast_parameters,
    nuexact_parameters,
)
from hamiltonian import (
    oscillation_probabilities,
    matter_potential as reference_matter_potential,
    KM_TO_INV_EV as REF_KM_TO_INV_EV,
)

# Hamiltonian adapter
def hamiltonian_probabilities(L, E, ordering="NO", rho=0.0, Ye=0.5, matter=True, antineutrino=False, method="eigh"):
    return oscillation_probabilities(L=L, E=E, ordering=ordering, rho=rho, Ye=Ye, matter=matter, antineutrino=antineutrino, method=method)

# Call NuOscProbExact
PROJECT_DIR = Path(__file__).resolve().parent
NUEXACT_SRC = PROJECT_DIR.parent / "NuOscProbExact-main" / "src"
if not NUEXACT_SRC.exists():
    raise FileNotFoundError(
        "NuOscProbExact source directory was not found:\n"
        f"{NUEXACT_SRC}\n\n"
        "Expected layout: Make sure 'NuOscProbExact-main/src/' exists in the same directory as the benchmark script."
    )
sys.path.insert(0, str(NUEXACT_SRC))

import earth as nuexact_earth
import globaldefs as nuexact_gd
import hamiltonians3nu
import oscprob3nu

# NuOscProbExact adapter
def nuexact_probabilities(L, E, ordering="NO", rho=0.0, Ye=0.5, matter=True, antineutrino=False, mode="harmonised"):
    mode = mode.lower()
    if mode not in {"harmonised", "native"}:
        raise ValueError("mode must be 'harmonised' or 'native'")

    E = np.atleast_1d(np.asarray(E, dtype=float))
    E_eV = E * 1.0e9
    p = nuexact_parameters(ordering, mode=mode, gd=nuexact_gd)
    delta = -p["delta"] if antineutrino else p["delta"]

    if mode == "harmonised":
        H0 = hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(
            p["s12"],
            p["s23"],
            p["s13"],
            delta,
            p["dm21"],
            p["dm31"],
            angles=p["angles"],
        )

        L_eVinv = float(L) * REF_KM_TO_INV_EV

        if matter:
            Vcc = reference_matter_potential(rho, Ye)
            if antineutrino:
                Vcc = -Vcc
            H = hamiltonians3nu.hamiltonian_3nu_matter(H0, E_eV, Vcc)
        else:
            H = H0 / E_eV[:, None, None]

    else:
        H0 = hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(
            p["s12"],
            p["s23"],
            p["s13"],
            delta,
            p["dm21"],
            p["dm31"],
            angles=p["angles"],
        )

        L_eVinv = float(L) * nuexact_gd.CONV_KM_TO_INV_EV

        if matter:
            Vcc = nuexact_earth.matter_potential(float(rho), float(Ye))
            if antineutrino:
                Vcc = -Vcc
            H = hamiltonians3nu.hamiltonian_3nu_matter(H0, E_eV, Vcc)
        else:
            H = H0 / E_eV[:, None, None]

    P = np.asarray(oscprob3nu.probabilities_3nu(H, L_eVinv))

    # Current NuOscProbExact returns (..., 9) for a batched scan.
    return P.reshape(-1, 3, 3)


# NuFast adapter
_NUFAST = {}
_DLL_DIR_HANDLES = []
_DOUBLE_PTR = ctypes.POINTER(ctypes.c_double)

def _library_name(mode):
    mode = mode.lower()
    if mode == "harmonised":
        base = "nufast_harmonised"
    elif mode == "native":
        base = "nufast_native"
    else:
        raise ValueError("mode must be 'harmonised' or 'native'")

    if sys.platform.startswith("win"):
        return f"{base}.dll"

    if sys.platform == "darwin":
        return f"lib{base}.dylib"

    return f"lib{base}.so"

def _add_windows_runtime_directory():
    """Help Windows locate MSYS2/UCRT64 runtime DLLs, if required."""
    if not sys.platform.startswith("win"):
        return

    candidates = []

    gpp = shutil.which("g++")
    if gpp:
        candidates.append(Path(gpp).resolve().parent)

    env_dir = os.environ.get("MSYS2_UCRT64_BIN")
    if env_dir:
        candidates.append(Path(env_dir))

    candidates.append(Path(r"C:\msys64\ucrt64\bin"))

    seen = set()

    for directory in candidates:
        key = str(directory).lower()
        if key in seen or not directory.exists():
            continue
        seen.add(key)
        try:
            # Keep the handle alive for the life of the Python process.
            _DLL_DIR_HANDLES.append(os.add_dll_directory(str(directory)))
        except (AttributeError, OSError):
            pass

def _load_nufast(mode):
    mode = mode.lower()
    if mode in _NUFAST:
        return _NUFAST[mode]
    path = Path(__file__).with_name(_library_name(mode))
    if not path.exists():
        raise FileNotFoundError(
            f"{path.name} is missing.\n\n"
            "Create both NuFast libraries first:\n"
            "  py build_nufast.py\n"
        )
    _add_windows_runtime_directory()

    lib = ctypes.CDLL(str(path))
    lib.nufast_matter.argtypes = [
        ctypes.c_double,  # s12sq
        ctypes.c_double,  # s13sq
        ctypes.c_double,  # s23sq
        ctypes.c_double,  # delta
        ctypes.c_double,  # dm21
        ctypes.c_double,  # dm31
        ctypes.c_double,  # L
        _DOUBLE_PTR,      # energies
        ctypes.c_int,     # nE
        ctypes.c_double,  # rhoYe
        ctypes.c_int,     # N_Newton
        _DOUBLE_PTR,      # output
    ]
    lib.nufast_matter.restype = None

    lib.nufast_vacuum.argtypes = [
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_double,
        _DOUBLE_PTR,
        ctypes.c_int,
        _DOUBLE_PTR,
    ]
    lib.nufast_vacuum.restype = None

    _NUFAST[mode] = lib
    return lib

def nufast_probabilities(L, E, ordering="NO", rho=0.0, Ye=0.5, matter=True, antineutrino=False, n_newton=2, mode="harmonised"):
    mode = mode.lower()
    lib = _load_nufast(mode)
    p = nufast_parameters(ordering, mode=mode)
    energies = np.ascontiguousarray(np.atleast_1d(np.asarray(E, dtype=float)), dtype=np.float64) 
    # NuFast convention: negative E denotes antineutrinos.
    if antineutrino:
        energies = -energies
    out = np.empty((energies.size, 3, 3), dtype=np.float64)

    common = (
        p["s12sq"],
        p["s13sq"],
        p["s23sq"],
        p["delta"],
        p["dm21"],
        p["dm31"],
        float(L),
        energies.ctypes.data_as(_DOUBLE_PTR),
        int(energies.size),
    )

    if matter:
        lib.nufast_matter(*common, float(rho * Ye), int(n_newton), out.ctypes.data_as(_DOUBLE_PTR))
    else:
        lib.nufast_vacuum(*common, out.ctypes.data_as(_DOUBLE_PTR))

    return out


# Benchmark interface
def run_case(
    experiment_name="DUNE",
    ordering="NO",
    particle="nu",
    medium="matter",
    n_newton=2,
    n_energy=None,
    mode="harmonised",
    hamiltonian_method="eigh",
):

    experiment_name = experiment_name.upper()
    ordering = ordering.upper()
    particle = particle.lower()
    medium = medium.lower()
    mode = mode.lower()

    hamiltonian_method = hamiltonian_method.lower()
    if hamiltonian_method not in {"eigh", "expm"}:
        raise ValueError("hamiltonian_method must be 'eigh' or 'expm'")

    if particle not in {"nu", "antinu"}:
        raise ValueError("particle must be 'nu' or 'antinu'")

    if medium not in {"vacuum", "matter"}:
        raise ValueError("medium must be 'vacuum' or 'matter'")

    if mode not in {"harmonised", "native"}:
        raise ValueError("mode must be 'harmonised' or 'native'")

    cfg = experiment(experiment_name, n_energy=n_energy)

    anti = particle == "antinu"
    use_matter = medium == "matter"

    L = cfg["L"]
    E = cfg["E"]
    rho = cfg["rho"]
    Ye = cfg["Ye"]
    alpha = cfg["alpha"]
    beta = cfg["beta"]

    P_ref = hamiltonian_probabilities(L, E, ordering=ordering, rho=rho, Ye=Ye, matter=use_matter, antineutrino=anti, method=hamiltonian_method)
    P_exact = nuexact_probabilities(L, E, ordering=ordering, rho=rho, Ye=Ye, matter=use_matter, antineutrino=anti, mode=mode)
    P_fast = nufast_probabilities(L, E, ordering=ordering, rho=rho, Ye=Ye, matter=use_matter, antineutrino=anti, n_newton=n_newton, mode=mode)

    ref = P_ref[:, alpha, beta]
    exact = P_exact[:, alpha, beta]
    fast = P_fast[:, alpha, beta]

    return {
        "experiment": experiment_name,
        "ordering": ordering,
        "particle": particle,
        "medium": medium,
        "mode": mode,
        "hamiltonian_method": hamiltonian_method,
        "L": L,
        "rho": rho,
        "Ye": Ye,
        "alpha": alpha,
        "beta": beta,
        "energy": E,
        "hamiltonian": ref,
        "NuOscProbExact": exact,
        "NuFast": fast,
        "diff_NuOscProbExact": np.abs(exact - ref),
        "diff_NuFast": np.abs(fast - ref),
        "diff_NuFast_NuOscProbExact": np.abs(fast - exact),
        "conservation_hamiltonian": np.max(np.abs(np.sum(P_ref, axis=2) - 1.0)),
        "conservation_NuOscProbExact": np.max(np.abs(np.sum(P_exact, axis=2) - 1.0)),
        "conservation_NuFast": np.max(np.abs(np.sum(P_fast, axis=2) - 1.0)),
    }
