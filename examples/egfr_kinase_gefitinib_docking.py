"""Example: Molecular docking of Gefitinib into EGFR tyrosine kinase."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from chem_binding import (
    create_egfr_kinase_pocket,
    create_gefitinib_ligand,
    MolecularDockingEngine,
    PharmacophoreAssay,
    render_interaction_summary,
    format_energy_breakdown,
)


def main() -> None:
    print("Preparing EGFR Tyrosine Kinase binding pocket and Gefitinib...")
    rec = create_egfr_kinase_pocket(t790m_mutant=False)
    lig = create_gefitinib_ligand()

    engine = MolecularDockingEngine()
    print("Executing Monte Carlo Basin-Hopping docking search...")
    result = engine.dock(lig, rec, num_cycles=15)

    print("\n" + format_energy_breakdown(result))
    print("\n" + render_interaction_summary(result, rec))

    assay = PharmacophoreAssay()
    report = assay.calculate_affinity(result)
    print("\nPharmacology Summary:")
    print(f"Predicted Dissociation Constant (Kd): {report.kd_nanomolar:.2f} nM")
    print(f"Potency: {report.potency_tier}")


if __name__ == "__main__":
    main()
