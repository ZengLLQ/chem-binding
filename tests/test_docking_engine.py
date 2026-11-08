"""Unit tests for molecular docking search."""

import numpy as np
from chem_binding.molecular.receptor import (
    create_egfr_kinase_pocket,
    create_gefitinib_ligand,
)
from chem_binding.optimization.docking_engine import MolecularDockingEngine


def test_docking_optimization() -> None:
    rec = create_egfr_kinase_pocket(t790m_mutant=False)
    lig = create_gefitinib_ligand()

    engine = MolecularDockingEngine(random_seed=42)
    result = engine.dock(lig, rec, num_cycles=5)

    assert result.docked_coords.shape == (lig.num_atoms, 3)
    assert not np.isnan(result.energy.delta_g_total)
    assert np.isclose(np.linalg.norm(result.quaternion), 1.0, atol=1e-3)
