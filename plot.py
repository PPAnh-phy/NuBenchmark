from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

from PMNS import FLAVOR_NAMES


def _for_log(values):
    """Hide exact zeros on a logarithmic residual axis."""
    values = np.asarray(values, dtype=float)
    return np.where(values > 0.0, values, np.nan)


def plot_result(result, output=None, show=False):
    E = result["energy"]
    if result["experiment"] == "JUNO":
        x = E * 1.0e3
        xlabel = "Energy [MeV]"
    else:
        x = E
        xlabel = "Energy [GeV]"

    alpha = FLAVOR_NAMES[result["alpha"]]
    beta = FLAVOR_NAMES[result["beta"]]

    particle_symbol = (
        r"\bar{\nu}"
        if result["particle"] == "antinu"
        else r"\nu"
    )

    reference_label = f'Hamiltonian ({result["hamiltonian_method"]})'

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 7), sharex=True, gridspec_kw={"height_ratios": [2.2, 1.0]})

    ax1.plot(x, result["hamiltonian"], label=reference_label)
    ax1.plot(x, result["NuOscProbExact"], "--", label="NuOscProbExact")
    ax1.plot(x, result["NuFast"], ":", label="NuFast")
    ax1.set_ylabel("Probability")
    ax1.set_ylim(0, 1.02)
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    ax1.set_title(
        f'{result["experiment"]} | {result["ordering"]} | {result["medium"]} | {result["mode"]}\n'
        rf'$P({particle_symbol}_{{{alpha}}}\rightarrow'
        rf'{particle_symbol}_{{{beta}}})$'
    )

    ax2.semilogy(x, _for_log(result["diff_NuOscProbExact"]), linestyle="-", alpha=0.75, label="|NuOscProbExact - Hamiltonian|")
    ax2.semilogy(x, _for_log(result["diff_NuFast"]), "-.", alpha=0.75, label="|NuFast - Hamiltonian|")
    ax2.semilogy(x, _for_log(result["diff_NuFast_NuOscProbExact"]), "--", alpha=0.75, label="|NuFast - NuOscProbExact|")
    ax2.set_ylabel("Absolute difference")
    ax2.grid(True, which="major", alpha=0.20, linewidth=0.6)
    ax2.grid(True, which="minor", alpha=0.01, linewidth=0.4)
    ax2.legend(loc="upper center", bbox_to_anchor=(0.5, -0.25), ncol=3, fontsize=8, frameon=False)

    fig.tight_layout()

    if output is not None:
        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, dpi=180, bbox_inches="tight")

    if show:
        plt.show()
    else:
        plt.close(fig)
