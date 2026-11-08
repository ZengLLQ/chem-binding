"""Clinical pharmacophore binding affinity assay and oncogenic resistance profiling."""

from __future__ import annotations

import math
from dataclasses import dataclass
import numpy as np

from chem_binding.optimization.docking_engine import DockingResult


@dataclass
class ClinicalAffinityReport:
    """Pharmacological potency and thermodynamic dissociation metrics."""

    ligand_name: str
    target_name: str
    delta_g_bind: float         # kcal/mol
    kd_molar: float             # Molar
    kd_nanomolar: float         # nM
    pic50: float                # -log10(Kd)
    potency_tier: str           # e.g. "Sub-nanomolar Potency"
    clinical_interpretation: str


@dataclass
class ResistanceMutationScreen:
    """Differential binding assay across wild-type and drug-resistant oncogenic target."""

    wild_type_result: DockingResult
    mutant_result: DockingResult
    delta_delta_g: float        # Delta G_mutant - Delta G_wt (kcal/mol)
    fold_resistance: float      # Ratio Kd_mutant / Kd_wt
    resistance_classification: str
    oncology_recommendations: list[str]


class PharmacophoreAssay:
    """Translates free energies of binding into clinical pharmacology and resistance metrics."""

    GAS_CONSTANT_R: float = 1.987204e-3  # kcal / (mol * K)
    TEMPERATURE_K: float = 298.15        # Physiological 25 C (298.15 K)

    def calculate_affinity(self, result: DockingResult) -> ClinicalAffinityReport:
        """Convert Delta G_bind to dissociation constant Kd and pIC50."""
        dg = result.energy.delta_g_total
        rt = self.GAS_CONSTANT_R * self.TEMPERATURE_K

        # Kd = exp(Delta G / RT)
        # Clamped exponent to prevent overflow
        exp_arg = max(-100.0, min(100.0, dg / rt))
        kd = float(math.exp(exp_arg))
        kd_nm = kd * 1e9
        pic50 = float(-math.log10(max(1e-18, kd)))

        if kd_nm < 1.0:
            tier = "Sub-nanomolar High Potency (< 1 nM)"
        elif kd_nm < 100.0:
            tier = "Nanomolar Potency (1 - 100 nM)"
        elif kd_nm < 10_000.0:
            tier = "Micromolar Moderate Affinity (0.1 - 10 uM)"
        else:
            tier = "Weak / Non-Specific Binder (> 10 uM)"

        interp = (
            f"Ligand '{result.ligand_name}' docked against '{result.receptor_name}'. "
            f"Predicted Delta G_bind = {dg:.2f} kcal/mol yields an estimated Kd = {kd_nm:.2f} nM (pIC50 = {pic50:.2f}). "
            f"Classified as {tier}."
        )

        return ClinicalAffinityReport(
            ligand_name=result.ligand_name,
            target_name=result.receptor_name,
            delta_g_bind=dg,
            kd_molar=kd,
            kd_nanomolar=kd_nm,
            pic50=pic50,
            potency_tier=tier,
            clinical_interpretation=interp,
        )

    def screen_mutation_resistance(
        self,
        wt_result: DockingResult,
        mutant_result: DockingResult,
    ) -> ResistanceMutationScreen:
        """Evaluate sensitivity shift and resistance fold-change for oncogenic gatekeeper mutants."""
        dg_wt = wt_result.energy.delta_g_total
        dg_mut = mutant_result.energy.delta_g_total

        ddg = dg_mut - dg_wt
        rt = self.GAS_CONSTANT_R * self.TEMPERATURE_K
        fold = float(math.exp(np.clip(ddg / rt, -50.0, 50.0)))

        if ddg >= 2.5:
            resist_class = "High Clinical Drug Resistance (> 50-fold affinity loss)"
            recs = [
                "Discontinue 1st/2nd generation reversible TKI due to gatekeeper steric blockage.",
                "Switch patient to 3rd-generation covalent inhibitor (e.g. Osimertinib / AZD9291).",
                "Sequence cell-free tumor DNA (ctDNA) for emergence of tertiary C797S resistance.",
            ]
        elif ddg >= 0.8:
            resist_class = "Moderate Resistance (4 to 50-fold affinity loss)"
            recs = [
                "Monitor for progressive disease or acquired clonal outgrowth.",
                "Consider combination therapy or dose optimization under TDM guidance.",
            ]
        else:
            resist_class = "Target Sensitive / Unaltered Binding"
            recs = ["Maintain current targeted therapy regimen.", "Routine RECIST 1.1 imaging."]

        return ResistanceMutationScreen(
            wild_type_result=wt_result,
            mutant_result=mutant_result,
            delta_delta_g=ddg,
            fold_resistance=fold,
            resistance_classification=resist_class,
            oncology_recommendations=recs,
        )
