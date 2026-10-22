"""
Settings:
    --oscpars common : Use reference oscillation parameters in all programs.
    --oscpars native : Use each program's native oscillation parameters.
    --Vcc common     : Use the reference matter conversion in all programs.
    --Vcc native     : Use each program's native matter conversion.

Examples:
py pipeline.py --oscpars common --Vcc common --experiment JUNO --ordering NO --particle antinu --medium matter --hamiltonian-method eigh --show
py pipeline.py --oscpars common --Vcc native --experiment JUNO --ordering NO --particle antinu --medium matter --hamiltonian-method eigh --show
py pipeline.py --oscpars native --Vcc common --experiment JUNO --ordering NO --particle antinu --medium matter --hamiltonian-method eigh --show
py pipeline.py --oscpars native --Vcc native --experiment JUNO --ordering NO --particle antinu --medium matter --hamiltonian-method eigh --show

# Display all available options:
py pipeline.py --help

Use python or python3 instead of py depending on your system configuration.
"""
from pathlib import Path
import argparse
from datetime import datetime
import shutil
import numpy as np

from PMNS import FLAVOR_NAMES, pmns_unitarity
from benchmark import run_case
from plot import plot_result


def case_name(experiment, ordering, particle, medium, oscpars, vcc):
    """Return the output stem for one benchmark case."""
    return (f"{experiment}_{ordering}_{particle}_{medium}_oscpars-{oscpars}_Vcc-{vcc}")


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
        f'P({alpha}->{beta})'
    )
    print(f'Hamiltonian solver: {result["hamiltonian_method"]}')

    print(
        f'  Settings: oscpars={result["oscpars"]}, '
        f'Vcc={result["Vcc"]}'
    )

    print(
        "  max absolute probability difference |NuOscProbExact - Hamiltonian| = "
        f'{result["diff_NuOscProbExact"].max():.6e}'
    )

    print(
        "  max absolute probability difference |NuFast - Hamiltonian| = "
        f'{result["diff_NuFast"].max():.6e}'
    )

    print(
        "  max absolute probability difference |NuFast - NuOscProbExact| = "
        f'{result["diff_NuFast_NuOscProbExact"].max():.6e}'
    )

    print(
        "Probability conservation: "
        f'H={result["conservation_hamiltonian"]:.3e}, '
        f'Exact={result["conservation_NuOscProbExact"]:.3e}, '
        f'Fast={result["conservation_NuFast"]:.3e}'
    )


def run_one(experiment, ordering, particle, medium, args):
    name = case_name(experiment, ordering, particle, medium, args.oscpars, args.vcc)
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
        oscpars=args.oscpars,
        Vcc=args.vcc,
    )
    print_summary(result)

    # Both files in the timestamp folder.
    csv_path = out / f"{name}.csv"
    png_path = out / f"{name}.png"

    save_probabilities(result, csv_path)
    plot_result(result, png_path, show=args.show)

    # Seperate folders for all CSV and PNG results
    csv_dir = out.parent / "_csv"
    png_dir = out.parent / "_png"

    csv_dir.mkdir(parents=True, exist_ok=True)
    png_dir.mkdir(parents=True, exist_ok=True)

    run_id = out.name
    shutil.copy2(csv_path, csv_dir / f"{run_id}_{name}.csv")
    shutil.copy2(png_path, png_dir / f"{run_id}_{name}.png")

    print()


def main():
    parser = argparse.ArgumentParser(
        description=("Benchmark a three-flavour Hamiltonian solver against NuOscProbExact and NuFast-LBL.")
    )

    parser.add_argument("--oscpars", choices=["common", "native"], default="common", help="Oscillation-parameter selection.")
    parser.add_argument("--Vcc", dest="vcc", choices=["common", "native"], default="common", help="Matter-potential conversion selection.")
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
