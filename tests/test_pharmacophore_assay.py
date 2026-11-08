"""Unit tests for clinical affinity calculations and resistance screening."""

import numpy as np
from chem_binding.molecular.receptor import (
    create_egfr_kinase_pocket,
    create_gefitinib_ligand,
)
from chem_binding.scoring.free_energy import BindingEnergyComponents
from chem_binding.optimization.docking_engine import DockingResult
from chem_binding.clinical.pharmacophore_assay import PharmacophoreAssay


def test_affinity_conversion() -> None:
    assay = PharmacophoreAssay()
    # A Delta G of ~ -10.9 kcal/mol gives Kd ~ 10 nM at 298.15 K
    comp = BindingEnergyComponents(
        delta_g_total=-10.9,
        vdw_energy=-6.0,
        electrostatic_energy=-2.0,
        hbond_energy=-3.0,
        desolvation_energy=-0.5,
        torsional_penalty=0.6,
    )
    res = DockingResult(
        ligand_name="Gefitinib",
        receptor_name="EGFR",
        docked_coords=np.zeros((8, 3)),
        translation=np.zeros(3),
        quaternion=np.array([1.0, 0.0, 0.0, 0.0]),
        dihedrals=np.zeros(0),
        energy=comp,
        iterations=5,
    )

    report = assay.calculate_affinity(res)
    assert report.kd_nanomolar < 100.0  # Nanomolar range
    assert report.pic50 > 7.0
    assert "Potency" in report.potency_tier


def test_mutation_resistance_screen() -> None:
    assay = PharmacophoreAssay()
    res_wt = DockingResult(
        ligand_name="Gefitinib",
        receptor_name="EGFR WT",
        docked_coords=np.zeros((8, 3)),
        translation=np.zeros(3),
        quaternion=np.array([1.0, 0.0, 0.0, 0.0]),
        dihedrals=np.zeros(0),
        energy=BindingEnergyComponents(-11.0, -6.0, -2.0, -3.5, -0.5, 1.0),
        iterations=5,
    )
    res_mut = DockingResult(
        ligand_name="Gefitinib",
        receptor_name="EGFR T790M",
        docked_coords=np.zeros((8, 3)),
        translation=np.zeros(3),
        quaternion=np.array([1.0, 0.0, 0.0, 0.0]),
        dihedrals=np.zeros(0),
        energy=BindingEnergyComponents(-7.5, -4.0, -1.0, -2.0, -0.5, 0.0),
        iterations=5,
    )

    screen = assay.screen_mutation_resistance(res_wt, res_mut)
    assert screen.delta_delta_g > 3.0  # Significant loss of affinity
    assert screen.fold_resistance > 100.0
    assert "High Clinical Drug Resistance" in screen.resistance_classification
