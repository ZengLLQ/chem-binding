"""Unit tests for semi-empirical free energy of binding scoring."""

import numpy as np
from chem_binding.molecular.receptor import (
    create_egfr_kinase_pocket,
    create_gefitinib_ligand,
)
from chem_binding.scoring.free_energy import EmpiricalScoringFunction


def test_free_energy_evaluation() -> None:
    rec = create_egfr_kinase_pocket(t790m_mutant=False)
    lig = create_gefitinib_ligand()
    scorer = EmpiricalScoringFunction()

    # Place ligand directly in pocket
    coords = lig.base_coords.copy()
    energy = scorer.evaluate(coords, lig, rec)

    assert not np.isnan(energy.delta_g_total)
    assert not np.isinf(energy.delta_g_total)
    assert energy.vdw_energy != 0.0
    assert energy.electrostatic_energy != 0.0
