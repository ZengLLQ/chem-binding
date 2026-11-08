"""ASCII visualizer for molecular docking poses and pharmacophore interaction maps."""

from __future__ import annotations

import numpy as np

from chem_binding.optimization.docking_engine import DockingResult
from chem_binding.molecular.receptor import ReceptorPocket
from chem_binding.clinical.pharmacophore_assay import (
    ClinicalAffinityReport,
    ResistanceMutationScreen,
)


def render_interaction_summary(result: DockingResult, receptor: ReceptorPocket) -> str:
    """Render ASCII table of closest contact residues and hydrogen bonds."""
    lines: list[str] = []
    lines.append("=== Protein-Ligand 3D Binding Pocket Interaction Map ===")

    # Find closest receptor residues for each ligand atom
    diffs = result.docked_coords[:, np.newaxis, :] - receptor.coords[np.newaxis, :, :]
    dists = np.sqrt(np.sum(diffs**2, axis=-1))

    header = f"{'Ligand Atom':<14} | {'Receptor Residue':<16} | {'Distance (A)':<14} | {'Interaction Type':<20}"
    lines.append(header)
    lines.append("-" * len(header))

    contacts_found = 0
    for i in range(len(result.docked_coords)):
        closest_rec_idx = int(np.argmin(dists[i]))
        d = float(dists[i, closest_rec_idx])
        if d <= 4.2:
            contacts_found += 1
            lig_name = f"Lig-{i}"
            rec_name = receptor.atoms[closest_rec_idx].name

            # Determine interaction type
            rec_atom = receptor.atoms[closest_rec_idx]
            if d <= 3.2 and (rec_atom.is_donor or rec_atom.is_acceptor):
                itype = "Hydrogen Bond"
            elif rec_atom.is_hydrophobic:
                itype = "Hydrophobic Contact"
            elif abs(rec_atom.charge) > 0.4:
                itype = "Electrostatic / Salt"
            else:
                itype = "Van der Waals"

            lines.append(f"{lig_name:<14} | {rec_name:<16} | {d:6.2f} A        | {itype:<20}")

    if contacts_found == 0:
        lines.append("No close contacts under 4.2 A detected.")

    return "\n".join(lines)


def format_energy_breakdown(res: DockingResult) -> str:
    """Format energy components into ASCII summary table."""
    e = res.energy
    lines: list[str] = []
    lines.append(f"Target Receptor:     {res.receptor_name}")
    lines.append(f"Ligand Drug:         {res.ligand_name}")
    lines.append(f"Delta G_bind Total:  {e.delta_g_total:+.2f} kcal/mol")
    lines.append("\nEnergy Decomposition:")
    lines.append(f"  - Van der Waals:       {e.vdw_energy:+6.2f} kcal/mol")
    lines.append(f"  - Electrostatic:       {e.electrostatic_energy:+6.2f} kcal/mol")
    lines.append(f"  - Hydrogen Bonding:    {e.hbond_energy:+6.2f} kcal/mol")
    lines.append(f"  - Desolvation:         {e.desolvation_energy:+6.2f} kcal/mol")
    lines.append(f"  - Torsional Penalty:   {e.torsional_penalty:+6.2f} kcal/mol")
    return "\n".join(lines)
