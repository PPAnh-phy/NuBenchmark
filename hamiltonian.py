"""
Three-flavour oscillations probability computation in vacuum and constant-density matter.

Inputs: E in GeV, L in km, rho in g/cm^3, delta in radians, dm in eV^2.
Outputs: Returns probabilities shaped (N_energy, 3, 3), including scalar energy.
"""
import numpy as np

from PMNS import pmns_matrix, reference_parameters

HBARC_EV_M = 197.3269804e-9
KM_TO_INV_EV = 1.0e3 / HBARC_EV_M
GF_GEV = 1.1663787e-5
NA = 6.02214076e23
DENSITY_CONVERSION = 7.645373e-33
VCC_PER_RHOYE = np.sqrt(2.0) * GF_GEV * NA * DENSITY_CONVERSION
YERHOE2A = 2.0e9 * VCC_PER_RHOYE


def _scalar(name, value):
    value = np.asarray(value, dtype=float)
    if value.ndim != 0 or not np.isfinite(value):
        raise ValueError(f'{name} must be a finite scalar')
    return float(value)


def matter_potential(rho, Ye=0.5):
    rho = _scalar('rho', rho)
    Ye = _scalar('Ye', Ye)
    if rho < 0.0 or not 0.0 <= Ye <= 1.0:
        raise ValueError('Require rho >= 0 and 0 <= Ye <= 1')
    return VCC_PER_RHOYE * rho * Ye


def build_hamiltonian(E, ordering="NO", parameters=None, *, rho=0.0, Ye=0.5, matter=True, antineutrino=False):
    """Return H [eV], shaped (N_energy, 3, 3)"""
    E = np.atleast_1d(np.asarray(E, dtype=float))
    if (
        E.ndim != 1
        or E.size == 0
        or not np.all(np.isfinite(E))
        or np.any(E <= 0.0)
    ):
        raise ValueError(
            "E must be a positive finite scalar or nonempty 1-D array"
        )

    source = (
        reference_parameters(ordering)
        if parameters is None
        else parameters
    )

    names = ("s12sq", "s13sq", "s23sq", "delta", "dm21", "dm31")
    p = {name: _scalar(name, source[name]) for name in names}

    for name in ("s12sq", "s13sq", "s23sq"):
        if not 0.0 <= p[name] <= 1.0:
            raise ValueError(f"{name} must lie in [0, 1]")

    U = pmns_matrix(ordering=ordering, parameters=p)

    if antineutrino:
        U = U.conj()

    M2 = np.diag([0.0, p["dm21"], p["dm31"]])
    H0 = U @ M2 @ U.conj().T / 2.0
    H = H0[None, :, :] / (E[:, None, None] * 1.0e9)

    if matter:
        Vcc = matter_potential(rho, Ye)
        H[:, 0, 0] += -Vcc if antineutrino else Vcc

    return H


def probabilities_from_hamiltonian(H, L, method='eigh'):
    """Propagate a constant Hermitian H [eV] over scalar L [km]. Accepts (3, 3) or (N_energy, 3, 3). Returns (N_energy, 3, 3)
    eigh: numerical Hermitian diagonalisation and spectral exponential.
    expm: scipy.linalg.expm applied separately to each matrix.
    """
    method = method.lower()
    if method not in {'eigh', 'expm'}:
        raise ValueError("method must be 'eigh' or 'expm'")
    
    L = _scalar('L', L)
    if L < 0:
        raise ValueError('L must be nonnegative')
    
    H = np.asarray(H, dtype=np.complex128)
    if H.shape == (3, 3):
        H = H[None, :, :]
    if H.ndim != 3 or H.shape[1:] != (3, 3) or H.shape[0] == 0:
        raise ValueError('H must have shape (3, 3) or (N_energy, 3, 3)')
    if not np.all(np.isfinite(H)):
        raise ValueError('H must contain only finite values')
    
    Hdag = H.conj().transpose(0, 2, 1)
    scale = np.max(np.abs(H), axis=(1, 2))
    defect = np.max(np.abs(H - Hdag), axis=(1, 2))
    
    if np.any(defect > 1.0e-12 * scale):
        raise ValueError('H must be Hermitian to relative tolerance 1e-12')
    
    if L == 0:
        return np.broadcast_to(np.eye(3), H.shape).copy()
    L_natural = L * KM_TO_INV_EV
    
    if method == 'eigh':
        eigenvalues, W = np.linalg.eigh(H)
        phase = np.exp(-1j * eigenvalues * L_natural)
        S = (W * phase[:, None, :]) @ W.conj().transpose(0, 2, 1)
    else:
        try:
            from scipy.linalg import expm
        except ImportError as exc:
            raise ImportError("method='expm' requires SciPy: python -m pip install scipy") from exc
        S = np.stack([expm(-1j * h * L_natural) for h in H])
    # S has [final, initial]; expose P as [initial, final].
    return (np.abs(S) ** 2).transpose(0, 2, 1)


def oscillation_probabilities(L, E, ordering="NO", parameters=None, *, rho=0.0, Ye=0.5, matter=True, antineutrino=False, method="eigh"):
    """Return all nine probabilities [energy, initial, final]."""
    H = build_hamiltonian(E, ordering=ordering, parameters=parameters, rho=rho, Ye=Ye, matter=matter, antineutrino=antineutrino)
    return probabilities_from_hamiltonian(H, L, method=method)
