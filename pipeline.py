"""Command-line entry point for the neutrino probability benchmark.

General syntax
--------------
py pipeline.py 
    --mode <harmonised|native> --hamiltonian-method <eigh|expm>
    --experiment <JUNO|DUNE> --ordering <NO|IO>
    --particle <nu|antinu> --medium <vacuum|matter>

Enter the command on one line. The line breaks above are for readability.
On Linux/macOS, replace "py" with "python3".

Benchmark modes
---------------
harmonised: 
    Every method is benchmark with the same oscillation parameters and natural-unit conventions and matter-potential coefficient of the Hamiltonian, follow NuFIT 6.0 values.
    The parameter set and conventions expressed in the Hamiltonian prefer precision over computational efficiency.
native:
    The parameters also follow NuFIT 6.0 values, while each package's own input style and natural-unit conventions and matter-potential coefficient remain.
    We want to see the effect on precision of the results from each package's own parameters and conventions expression.

Hamiltonian methods
-------------------
eigh: Numerical Hermitian diagonalisation (default).
expm: Matrix exponential using scipy.linalg.expm.

Optional arguments
------------------
--points N       Number of energy samples; default: experiment setting.
--n-newton N     NuFast matter eigenvalue refinement; default: 2.
--output DIR     Output root directory; default: results.
--show           Display plots after a run.
--all            Run every combination of experiment, ordering, particle and medium for the selected mode and Hamiltonian method.

Examples
--------
py pipeline.py --mode harmonised --hamiltonian-method eigh --experiment JUNO --ordering NO --particle antinu --medium matter --show
py pipeline.py --mode native --hamiltonian-method expm --experiment JUNO --ordering IO --particle antinu --medium vacuum --show
py pipeline.py --mode harmonised --hamiltonian-method expm --experiment DUNE --ordering NO --particle nu --medium matter --show
py pipeline.py --mode native --hamiltonian-method eigh --experiment DUNE --ordering IO --particle antinu --medium matter --show
py pipeline.py --all --mode harmonised --hamiltonian-method expm
"""
from pathlib import Path
import argparse
from datetime import datetime
import numpy as np

from PMNS import FLAVOR_NAMES, pmns_unitarity
from benchmark import run_case
from plot import plot_result


def case_name(experiment, ordering, particle, medium, mode):
    """Return the output stem for one benchmark case."""
    return (f"{experiment}_{ordering}_{particle}_{medium}_{mode}")


def save_probabilities(result, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    data = np.column_stack([result["energy"],
                            result["hamiltonian"],
                            result["NuOscProbExact"],
                            result["NuFast"],
                            result["diff_NuOscProbExact"],
                            result["diff_NuFast"],
                            result["diff_NuFast_NuOscProbExact"]])

    np.savetxt(path, 
               data, 
               delimiter=",", 
               header=("energy_GeV,"
                       "hamiltonian,"
                       "NuOscProbExact,"
                       "NuFast,"
                       "absdiff_NuOscProbExact,"
                       "absdiff_NuFast,"
                       "absdiff_NuFast_NuOscProbExact"), comments="")


def print_summary(result):
    alpha = FLAVOR_NAMES[result["alpha"]]
    beta = FLAVOR_NAMES[result["beta"]]

    print(
        f'{result["experiment"]} | '
        f'{result["ordering"]} | '
        f'{result["particle"]} | '
        f'{result["medium"]} | '
        f'{result["mode"]} | '
        f'P({alpha}->{beta})'
    )
    print(f'Hamiltonian solver: {result["hamiltonian_method"]}')

    if result["mode"] == "harmonised":
        print(
            "  max numerical difference of probabilities |NuOscProbExact - Hamiltonian| = "
            f'{result["diff_NuOscProbExact"].max():.6e}'
        )
        print(
            "  max numerical difference of probabilities |NuFast - Hamiltonian|         = "
            f'{result["diff_NuFast"].max():.6e}'
        )
    else:
        print(
            "  max numerical difference of probabilities |NuOscProbExact - Hamiltonian| = "
            f'{result["diff_NuOscProbExact"].max():.6e}'
        )
        print(
            "  max numerical difference of probabilities |NuFast - Hamiltonian|         = "
            f'{result["diff_NuFast"].max():.6e}'
        )

    print(
        "  max numerical difference of probabilities |NuFast - NuOscProbExact|      = "
        f'{result["diff_NuFast_NuOscProbExact"].max():.6e}'
    )

    print(
        "Probability conservation: "
        f'H={result["conservation_hamiltonian"]:.3e}, '
        f'Exact={result["conservation_NuOscProbExact"]:.3e}, '
        f'Fast={result["conservation_NuFast"]:.3e}'
    )


def run_one(experiment, ordering, particle, medium, args):
    name = case_name(experiment, ordering, particle, medium, args.mode)
    name = f"{name}_{args.hamiltonian_method}"
    out = Path(args.output)
    result = run_case(
        experiment_name=experiment,
        ordering=ordering,
        particle=particle,
        medium=medium,
        hamiltonian_method=args.hamiltonian_method,
        n_newton=args.n_newton,
        n_energy=args.points,
        mode=args.mode,
    )
    print_summary(result)
    save_probabilities(result, out / f"{name}.csv")
    plot_result(result, out / f"{name}.png", show=args.show)
    print()


def main():
    parser = argparse.ArgumentParser(
        description=("Benchmark a three-flavour Hamiltonian solver against NuOscProbExact and NuFast-LBL.")
    )

    parser.add_argument(
        "--mode",
        choices=["harmonised", "native"],
        default="harmonised",
        help=(
            "harmonised: all methods use the Hamiltonian's parameter and unit conventions; "
            "native: same NuFIT 6.0 physics values, but each package keeps its own input/unit/matter convention"
        ),
    )

    parser.add_argument("--experiment", choices=["JUNO", "DUNE"], default="DUNE")
    parser.add_argument("--ordering", choices=["NO", "IO"], default="NO")
    parser.add_argument("--particle", choices=["nu", "antinu"], default="nu")
    parser.add_argument("--medium", choices=["vacuum", "matter"], default="matter")
    parser.add_argument("--hamiltonian-method", choices=["eigh", "expm"], default="eigh", help="Method used by the Hamiltonian reference solver")
    parser.add_argument("--n-newton", type=int, default=2)
    parser.add_argument("--points", type=int, default=None)
    parser.add_argument("--output", default="results")
    parser.add_argument("--show", action="store_true")
    parser.add_argument("--all", action="store_true")

    args = parser.parse_args()
    
    output_root = Path(args.output)
    output_root.mkdir(parents=True, exist_ok=True)

    while True:
        run_id = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
        run_dir = output_root / run_id
        try:
            run_dir.mkdir(exist_ok=False)
            break
        except FileExistsError:
            continue

    args.output = str(run_dir)

    print(f"Results folder: {run_dir.resolve()}\n")
    print(f"PMNS Unitarity (NO): {pmns_unitarity('NO'):.3e}")
    print(f"PMNS Unitarity (IO): {pmns_unitarity('IO'):.3e}")
    print()

    if args.all:
        for ex in ("JUNO", "DUNE"):
            for order in ("NO", "IO"):
                for particle in ("nu", "antinu"):
                    for medium in ("vacuum", "matter"):
                        run_one(ex, order, particle, medium, args)
    else:
        run_one(args.experiment, args.ordering, args.particle, args.medium, args)


if __name__ == "__main__":
    main()
