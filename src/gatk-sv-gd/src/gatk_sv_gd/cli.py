"""
CLI entry point for gatk-sv-gd.

Dispatches to subcommands: infer, call, plot, eval.
"""

import sys


SUBCOMMANDS = {
    "infer": "gatk_sv_gd.infer",
    "call": "gatk_sv_gd.call",
    "plot": "gatk_sv_gd.plot",
    "eval": "gatk_sv_gd.eval",
}

DESCRIPTIONS = {
    "infer": "Run Bayesian CNV inference at GD loci (Pyro model)",
    "call": "Call GD CNVs from model posterior probabilities",
    "plot": "Generate visualisation plots for GD CNV calls",
    "eval": "Evaluate GD CNV calls against a truth table",
}


def _print_usage():
    """Print top-level usage information."""
    prog = "gatk-sv-gd"
    print(f"Usage: {prog} <subcommand> [options]\n")
    print("Genomic Disorder CNV detection from binned read counts.\n")
    print("Subcommands:")
    for name, desc in DESCRIPTIONS.items():
        print(f"  {name:8s}  {desc}")
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
    sys.argv = [f"gatk-sv-gd {subcommand}"] + sys.argv[2:]

    # Lazy import to avoid loading heavy dependencies (torch, pyro, matplotlib)
    # for subcommands that don't need them.
    import importlib
    module = importlib.import_module(SUBCOMMANDS[subcommand])
    module.main()


if __name__ == "__main__":
    main()
