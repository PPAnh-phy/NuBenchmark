"""
Define the PMNS matrix and experimental parameters for neutrino oscillation probability calculations.
The PMNS paramerters follow the NuFIT 6.0 best-fit values.
"""
import numpy as np

FLAVORS = {"e": 0, "mu": 1, "tau": 2}
FLAVOR_NAMES = ("e", "mu", "tau")

# Modifiable by users
PMNS_PARAMETERS = {
    "NO": {
        "s12sq": 0.308,
        "s13sq": 0.02215,
        "s23sq": 0.470,
        "delta_deg": 212.0,
        "dm21": 7.49e-5,
        "dm3l": 2.513e-3,
    },
    "IO": {
        "s12sq": 0.308,
        "s13sq": 0.02224,
        "s23sq": 0.562,
        "delta_deg": 285.0,
        "dm21": 7.49e-5,
        "dm32": -2.510e-3,  
    },
}

def dm31(ordering):
    """Return dm31 [eV^2] from NuFIT's dm3l convention."""
    ordering = ordering.upper()
    p = PMNS_PARAMETERS[ordering]

    if ordering == "NO":
        return p["dm3l"]

    if ordering == "IO":
        return p["dm32"] + p["dm21"]

    raise ValueError("ordering must be 'NO' or 'IO'")


def reference_parameters(ordering):
    """Shared parameter set for harmonised mode"""
    ordering = ordering.upper()
    p = PMNS_PARAMETERS[ordering]

    return {
        "s12sq": p["s12sq"],
        "s13sq": p["s13sq"],
        "s23sq": p["s23sq"],
        "delta": np.deg2rad(p["delta_deg"]),
        "dm21": p["dm21"],
        "dm31": dm31(ordering),
    }


NUFAST_NATIVE_NO = {
    "s12sq": 0.31,
    "s13sq": 0.02,
    "s23sq": 0.47,
    "delta": -0.8 * np.pi,
    "dm21": 7.5e-5,
    "dm31": 2.5e-3,
}

NUFAST_NATIVE_IO = {
    **NUFAST_NATIVE_NO,
    "dm31": -2.5e-3,
}

def nufast_parameters(ordering, mode="harmonised"):
    ordering = ordering.upper()
    mode = mode.lower()

    if ordering not in {"NO", "IO"}:
        raise ValueError("ordering must be 'NO' or 'IO'")

    if mode == "harmonised":
        return reference_parameters(ordering)

    if mode == "native":
        source = (
            NUFAST_NATIVE_NO
            if ordering == "NO"
            else NUFAST_NATIVE_IO
        )
        return dict(source)

    raise ValueError("mode must be 'harmonised' or 'native'")


def nuexact_parameters(ordering, mode="harmonised", gd=None):
    ordering = ordering.upper()
    mode = mode.lower()

    if ordering not in {"NO", "IO"}:
        raise ValueError("ordering must be 'NO' or 'IO'")

    if mode == "harmonised":
        p = reference_parameters(ordering)

        return {
            "s12": p["s12sq"],
            "s23": p["s23sq"],
            "s13": p["s13sq"],
            "angles": "sin2",
            "delta": p["delta"],
            "dm21": p["dm21"],
            "dm31": p["dm31"],
        }

    if mode == "native":
        if gd is None:
            raise ValueError(
                "Native mode requires NuOscProbExact's globaldefs module"
            )

        suffix = f"{ordering}_BF"

        return {
            "s12": getattr(gd, f"S12_{suffix}"),
            "s23": getattr(gd, f"S23_{suffix}"),
            "s13": getattr(gd, f"S13_{suffix}"),
            "angles": "sin",
            "delta": getattr(gd, f"DCP_{suffix}"),
            "dm21": getattr(gd, f"D21_{suffix}"),
            "dm31": getattr(gd, f"D31_{suffix}"),
        }

    raise ValueError("mode must be 'harmonised' or 'native'")


# Physical configurations
# Modifiable by users
EXPERIMENTS = {
    "JUNO": {
        "L": 52.5,
        "E_min": 0.0018,
        "E_max": 0.0080,
        "N_energy": 1000,
        "rho": 2.825,
        "Ye": 0.5,
        "alpha": FLAVORS["e"],
        "beta": FLAVORS["e"],
    },
    "DUNE": {
        "L": 1285.0,
        "E_min": 0.1,
        "E_max": 10.0,
        "N_energy": 1000,
        "rho": 2.825,
        "Ye": 0.5,
        "alpha": FLAVORS["mu"],
        "beta": FLAVORS["e"],
    },
}


def pmns_matrix(ordering="NO", parameters=None):
    """Return the standard three-flavour PMNS matrix."""
    p = reference_parameters(ordering) if parameters is None else parameters

    s12 = np.sqrt(p["s12sq"])
    s13 = np.sqrt(p["s13sq"])
    s23 = np.sqrt(p["s23sq"])

    c12 = np.sqrt(1.0 - p["s12sq"])
    c13 = np.sqrt(1.0 - p["s13sq"])
    c23 = np.sqrt(1.0 - p["s23sq"])

    d = p["delta"]

    return np.array([
        [c12*c13, s12*c13, s13*np.exp(-1j*d)],
        [-s12*c23 - c12*s23*s13*np.exp(1j*d), c12*c23 - s12*s23*s13*np.exp(1j*d), s23*c13],
        [s12*s23 - c12*c23*s13*np.exp(1j*d), -c12*s23 - s12*c23*s13*np.exp(1j*d), c23*c13]
        ], dtype=complex)


def experiment(name, n_energy=None):
    """Return one experiment configuration."""
    cfg = dict(EXPERIMENTS[name.upper()])
    n = cfg["N_energy"] if n_energy is None else int(n_energy)
    cfg["E"] = np.linspace(cfg["E_min"], cfg["E_max"], n)
    return cfg


def pmns_unitarity(ordering="NO"):
    """Frobenius norm of U^dagger U - I for the Hamiltonian's PMNS matrix."""
    U = pmns_matrix(ordering)
    return np.linalg.norm(U.conj().T @ U - np.eye(3), ord="fro")
