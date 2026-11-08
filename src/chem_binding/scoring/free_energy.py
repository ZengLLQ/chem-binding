"""Empirical and physical binding free energy (Delta G_bind) scoring functions."""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from chem_binding.molecular.ligand import SmallMoleculeLigand
from chem_binding.molecular.receptor import ReceptorPocket


@dataclass
class BindingEnergyComponents:
    """Decomposition of protein-ligand binding free energy (kcal/mol)."""

    delta_g_total: float
    vdw_energy: float
    electrostatic_energy: float
    hbond_energy: float
    desolvation_energy: float
    torsional_penalty: float


class EmpiricalScoringFunction:
    """Semi-empirical physical force field for molecular docking scoring.

    Parameters
    ----------
    coulomb_factor : float
        Electrostatic constant (default 332.0636 kcal*A / (mol*e^2)).
    hbond_depth : float
        Well-depth for ideal hydrogen bonds (default -1.8 kcal/mol).
    hydrophobic_scale : float
        Favorable contact free energy for hydrophobic pairs (default -0.08 kcal/mol).
    torsional_entropy_cost : float
        Free energy cost per frozen rotatable bond (default 0.315 kcal/mol).
    """

    def __init__(
        self,
        coulomb_factor: float = 332.0636,
        hbond_depth: float = -1.80,
        hydrophobic_scale: float = -0.08,
        torsional_entropy_cost: float = 0.315,
    ) -> None:
        self.coulomb_const = float(coulomb_factor)
        self.e_hb = float(hbond_depth)
        self.e_hydro = float(hydrophobic_scale)
        self.w_tors = float(torsional_entropy_cost)

    def evaluate(
        self,
        ligand_coords: np.ndarray,
        ligand: SmallMoleculeLigand,
        receptor: ReceptorPocket,
    ) -> BindingEnergyComponents:
        """Compute Delta G_bind for a given ligand 3D conformation.

        Parameters
        ----------
        ligand_coords : np.ndarray
            (N_lig, 3) Cartesian coordinates.
        """
        # Pairwise distance matrix (N_lig, N_rec)
        # diffs: (N_lig, N_rec, 3)
        diffs = ligand_coords[:, np.newaxis, :] - receptor.coords[np.newaxis, :, :]
        dist_sq = np.sum(diffs**2, axis=-1)
        dist = np.sqrt(np.maximum(1e-8, dist_sq))

        # 1. Van der Waals (Buffered 12-6 Lennard-Jones)
        # sigma_ij = r_i + r_j (arithmetic mean), eps_ij = sqrt(eps_i * eps_j)
        sigmas = ligand_coords[:, 0:1] * 0.0 + (
            np.array([a.vdw_radius for a in ligand.atoms])[:, np.newaxis] + receptor.vdw_radii[np.newaxis, :]
        )
        epsilons = np.sqrt(
            np.array([a.vdw_well_depth for a in ligand.atoms])[:, np.newaxis] * receptor.vdw_eps[np.newaxis, :]
        )

        # Soft-core buffering at r < 0.7 * sigma to avoid infinite repulsive explosions
        r_soft = 0.75 * sigmas
        effective_dist = np.maximum(dist, r_soft)
        inv_r = sigmas / effective_dist
        inv_r6 = inv_r**6
        inv_r12 = inv_r6**2
        vdw_matrix = 4.0 * epsilons * (inv_r12 - inv_r6)

        # Repulsive penalty for steric penetration into r < r_soft
        penetration = np.maximum(0.0, r_soft - dist)
        steric_clash_penalty = np.sum(25.0 * (penetration**2))
        e_vdw = float(np.sum(vdw_matrix) + steric_clash_penalty)

        # 2. Screened Electrostatics (Distance-dependent dielectric epsilon(r) = 4r)
        q_lig = np.array([a.charge for a in ligand.atoms])[:, np.newaxis]
        q_rec = receptor.charges[np.newaxis, :]
        # V_elec = 332.0636 * q1 * q2 / (4 * r^2)
        elec_matrix = (self.coulomb_const * (q_lig * q_rec)) / (4.0 * (effective_dist**2))
        e_elec = float(np.sum(elec_matrix))

        # 3. Hydrogen Bonding (Donor - Acceptor pairs within 1.8 - 3.2 A)
        lig_donor = np.array([a.is_donor for a in ligand.atoms])[:, np.newaxis]
        lig_acc = np.array([a.is_acceptor for a in ligand.atoms])[:, np.newaxis]
        hb_active = (lig_donor & receptor.is_acceptor[np.newaxis, :]) | (
            lig_acc & receptor.is_donor[np.newaxis, :]
        )

        # 12-10 potential peaked at r0 = 2.8 A
        r0 = 2.80
        inv_hb = r0 / effective_dist
        hb_pot = self.e_hb * (5.0 * (inv_hb**12) - 6.0 * (inv_hb**10))
        # Cutoff beyond 3.5 A
        hb_mask = hb_active & (dist <= 3.5)
        e_hbond = float(np.sum(hb_pot[hb_mask])) if np.any(hb_mask) else 0.0

        # 4. Hydrophobic Desolvation (Non-polar contacts within 4.5 A)
        lig_hydro = np.array([a.is_hydrophobic for a in ligand.atoms])[:, np.newaxis]
        hydro_active = lig_hydro & receptor.is_hydrophobic[np.newaxis, :] & (dist <= 4.5)
        e_desolv = float(np.sum(self.e_hydro * np.exp(-((dist[hydro_active] - 3.5) ** 2)))) if np.any(hydro_active) else 0.0

        # 5. Torsional Penalty
        e_tors = float(self.w_tors * ligand.num_rotatable_bonds)

        # Net free energy
        total_g = e_vdw + e_elec + e_hbond + e_desolv + e_tors

        return BindingEnergyComponents(
            delta_g_total=total_g,
            vdw_energy=e_vdw,
            electrostatic_energy=e_elec,
            hbond_energy=e_hbond,
            desolvation_energy=e_desolv,
            torsional_penalty=e_tors,
        )
