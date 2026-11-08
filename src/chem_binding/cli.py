"""Command line interface for ChemBinding molecular docking and resistance assay."""

from __future__ import annotations

import argparse
import sys
import numpy as np

from chem_binding.molecular.receptor import (
    create_egfr_kinase_pocket,
    create_gefitinib_ligand,
)
from chem_binding.optimization.docking_engine import MolecularDockingEngine
from chem_binding.clinical.pharmacophore_assay import PharmacophoreAssay
from chem_binding.visualizer import render_interaction_summary, format_energy_breakdown


def handle_dock(args: argparse.Namespace) -> int:
    """Run molecular docking simulation."""
    print("=" * 76)
    print("  ChemBinding: Differentiable Molecular Docking Simulation")
    print("=" * 76)

    receptor = create_egfr_kinase_pocket(t790m_mutant=args.mutant)
    ligand = create_gefitinib_ligand()

    print(f"Target: {receptor.name}")
    print(f"Ligand: {ligand.name} ({ligand.num_atoms} atoms)")
    print(f"Searching {args.cycles} Monte Carlo Basin-Hopping cycles...")

    engine = MolecularDockingEngine()
    result = engine.dock(ligand, receptor, num_cycles=args.cycles)

    print("\n" + format_energy_breakdown(result))
    print("\n" + render_interaction_summary(result, receptor))

    assay = PharmacophoreAssay()
    report = assay.calculate_affinity(result)
    print("\n" + "=" * 76)
    print("  Clinical Pharmacology & Binding Affinity Assessment")
    print("=" * 76)
    print(f"Predicted Kd:     {report.kd_nanomolar:.2f} nM ({report.kd_molar:.2e} M)")
    print(f"Estimated pIC50:  {report.pic50:.2f}")
    print(f"Potency Tier:     {report.potency_tier}")
    print(f"Clinical Note:    {report.clinical_interpretation}")

    return 0


def handle_screen_mutations(args: argparse.Namespace) -> int:
    """Compare docking against Wild-Type vs T790M mutant EGFR kinase."""
    print("=" * 76)
    print("  ChemBinding: Oncogenic Drug Resistance Assay (EGFR WT vs T790M)")
    print("=" * 76)

    ligand = create_gefitinib_ligand()
    rec_wt = create_egfr_kinase_pocket(t790m_mutant=False)
    rec_mut = create_egfr_kinase_pocket(t790m_mutant=True)

    engine = MolecularDockingEngine()
    print("1. Docking Gefitinib into Wild-Type Sensitive EGFR...")
    res_wt = engine.dock(ligand, rec_wt, num_cycles=args.cycles)

    print("2. Docking Gefitinib into T790M Gatekeeper Resistant EGFR...")
    res_mut = engine.dock(ligand, rec_mut, num_cycles=args.cycles)

    assay = PharmacophoreAssay()
    screen = assay.screen_mutation_resistance(res_wt, res_mut)

    header = f"{'Target Condition':<28} | {'Delta G (kcal/mol)':<20} | {'Estimated Kd':<16}"
    print("\n" + header)
    print("-" * len(header))
    kd_wt = assay.calculate_affinity(res_wt).kd_nanomolar
    kd_mut = assay.calculate_affinity(res_mut).kd_nanomolar
    print(f"{'EGFR Wild-Type':<28} | {res_wt.energy.delta_g_total:+8.2f} kcal/mol     | {kd_wt:8.2f} nM")
    print(f"{'EGFR T790M Gatekeeper Mutant':<28} | {res_mut.energy.delta_g_total:+8.2f} kcal/mol     | {kd_mut:8.2f} nM")

    print("\nResistance Classification:")
    print(f"  Delta-Delta G:  {screen.delta_delta_g:+.2f} kcal/mol")
    print(f"  Fold Shift:     {screen.fold_resistance:.1f}x affinity loss")
    print(f"  Status:         {screen.resistance_classification}")
    print("Oncology Clinical Recommendations:")
    for rec in screen.oncology_recommendations:
        print(f"  * {rec}")

    return 0


def main(argv: list[str] | None = None) -> int:
    """CLI main entry point."""
    parser = argparse.ArgumentParser(
        prog="chem-binding",
        description="Differentiable Molecular Docking and Free Energy of Binding",
    )
    subparsers = parser.add_subparsers(dest="command", help="Subcommand to run")

    dock_p = subparsers.add_parser("dock", help="Dock ligand into protein pocket")
    dock_p.add_argument("--mutant", action="store_true", help="Dock into T790M mutant")
    dock_p.add_argument("--cycles", type=int, default=15, help="Number of Basin-Hopping cycles")

    screen_p = subparsers.add_parser("screen-mutations", help="Screen wild-type vs resistance mutation")
    screen_p.add_argument("--cycles", type=int, default=15, help="Number of Basin-Hopping cycles")

    args = parser.parse_args(argv)
    if args.command == "dock":
        return handle_dock(args)
    elif args.command == "screen-mutations":
        return handle_screen_mutations(args)
    else:
        parser.print_help()
        return 0


if __name__ == "__main__":
    sys.exit(main())
