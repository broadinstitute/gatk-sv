"""
CLI entry point for gatk-sv-ploidy.

Dispatches to subcommands: preprocess, infer, call, plot, eval.
"""

import sys


SUBCOMMANDS = {
    "preprocess": "gatk_sv_ploidy.preprocess",
    "infer": "gatk_sv_ploidy.infer",
    "call": "gatk_sv_ploidy.call",
    "plot": "gatk_sv_ploidy.plot",
    "eval": "gatk_sv_ploidy.eval",
}

DESCRIPTIONS = {
    "preprocess": "Normalise and filter depth data for aneuploidy inference",
    "infer": "Train Bayesian model and run discrete CN inference",
    "call": "Assign sex karyotype and aneuploidy type per sample",
    "plot": "Generate diagnostic and summary plots",
    "eval": "Evaluate predictions against a truth set",
}


def _print_usage():
    """Print top-level usage information."""
    prog = "gatk-sv-ploidy"
    print(f"Usage: {prog} <subcommand> [options]\n")
    print("Whole-genome aneuploidy detection from binned read counts.\n")
    print("Subcommands:")
    for name, desc in DESCRIPTIONS.items():
        print(f"  {name:12s}  {desc}")
    print(f"\nRun '{prog} <subcommand> --help' for subcommand-specific options.")


def main():
    """Main CLI dispatcher."""
    if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help"):
        _print_usage()
        sys.exit(0 if len(sys.argv) >= 2 else 1)

    subcommand = sys.argv[1]

    if subcommand not in SUBCOMMANDS:
        print(f"Error: unknown subcommand '{subcommand}'\n")
        _print_usage()
        sys.exit(1)

    # Rewrite sys.argv so the submodule's argparse sees the correct prog name
    sys.argv = [f"gatk-sv-ploidy {subcommand}"] + sys.argv[2:]

    # Lazy import to avoid loading heavy dependencies (torch, pyro, matplotlib)
    # for subcommands that don't need them.
    import importlib

    module = importlib.import_module(SUBCOMMANDS[subcommand])
    module.main()


if __name__ == "__main__":
    main()
