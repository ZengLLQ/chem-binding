"""Example: Screen EGFR gatekeeper mutation (T790M) for Gefitinib resistance."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from chem_binding import (
    create_egfr_kinase_pocket,
    create_gefitinib_ligand,
    MolecularDockingEngine,
    PharmacophoreAssay,
)


def main() -> None:
    print("Initiating comparative docking screen: EGFR WT vs T790M Mutant...")
    lig = create_gefitinib_ligand()
    rec_wt = create_egfr_kinase_pocket(t790m_mutant=False)
    rec_mut = create_egfr_kinase_pocket(t790m_mutant=True)

    engine = MolecularDockingEngine()
    print("Docking into Wild-Type...")
    res_wt = engine.dock(lig, rec_wt, num_cycles=15)

    print("Docking into T790M Mutant...")
    res_mut = engine.dock(lig, rec_mut, num_cycles=15)

    assay = PharmacophoreAssay()
    screen = assay.screen_mutation_resistance(res_wt, res_mut)

    print("\nResults:")
    print(f"  Delta G (WT):      {res_wt.energy.delta_g_total:.2f} kcal/mol")
    print(f"  Delta G (T790M):   {res_mut.energy.delta_g_total:.2f} kcal/mol")
    print(f"  Delta-Delta G:     {screen.delta_delta_g:+.2f} kcal/mol")
    print(f"  Fold Shift:        {screen.fold_resistance:.1f}x")
    print(f"  Classification:    {screen.resistance_classification}")
    print("\nOncology Action Items:")
    for rec in screen.oncology_recommendations:
        print(f"  - {rec}")


if __name__ == "__main__":
    main()
