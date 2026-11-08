"""Differentiable local optimization and Monte Carlo Basin-Hopping docking search."""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy.optimize import minimize

from chem_binding.molecular.ligand import SmallMoleculeLigand
from chem_binding.molecular.receptor import ReceptorPocket
from chem_binding.scoring.free_energy import (
    EmpiricalScoringFunction,
    BindingEnergyComponents,
)


@dataclass
class DockingResult:
    """Optimal docked pose and corresponding binding thermodynamics."""

    ligand_name: str
    receptor_name: str
    docked_coords: np.ndarray  # (N_lig, 3)
    translation: np.ndarray    # (3,)
    quaternion: np.ndarray     # (4,)
    dihedrals: np.ndarray      # (K,)
    energy: BindingEnergyComponents
    iterations: int


class MolecularDockingEngine:
    """Global and local search engine finding lowest Delta G_bind binding conformations."""

    def __init__(
        self,
        scoring_function: EmpiricalScoringFunction | None = None,
        random_seed: int | None = 42,
    ) -> None:
        self.scorer = scoring_function or EmpiricalScoringFunction()
        self.rng = np.random.default_rng(random_seed)

    def dock(
        self,
        ligand: SmallMoleculeLigand,
        receptor: ReceptorPocket,
        num_cycles: int = 15,
        temperature_k: float = 300.0,
    ) -> DockingResult:
        """Run Monte Carlo Basin-Hopping search around receptor pocket.

        Parameters
        ----------
        ligand : SmallMoleculeLigand
            Drug candidate molecule.
        receptor : ReceptorPocket
            Target protein pocket.
        num_cycles : int
            Number of Basin-Hopping trial perturbations.
        """
        # Initial guess: translation centered at receptor centroid
        rec_center = np.mean(receptor.coords, axis=0)
        lig_center = np.mean(ligand.base_coords, axis=0)
        t0 = rec_center - lig_center
        q0 = np.array([1.0, 0.0, 0.0, 0.0])
        phi0 = np.zeros(ligand.num_rotatable_bonds)

        def pack_params(t: np.ndarray, q: np.ndarray, phi: np.ndarray) -> np.ndarray:
            return np.concatenate([t, q, phi])

        def unpack_params(p: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
            t = p[:3]
            q = p[3:7]
            phi = p[7:]
            return t, q, phi

        def objective(p: np.ndarray) -> float:
            t, q, phi = unpack_params(p)
            coords = ligand.compute_conformation(t, q, phi)
            res = self.scorer.evaluate(coords, ligand, receptor)
            # Add penalty if quaternion deviates from unit sphere
            q_norm_pen = 100.0 * (np.linalg.norm(q) - 1.0) ** 2
            return float(res.delta_g_total + q_norm_pen)

        # Baseline optimization
        best_p = pack_params(t0, q0, phi0)
        curr_p = best_p.copy()
        curr_score = objective(curr_p)
        best_score = curr_score

        kb_t = 0.001987 * temperature_k  # kcal/mol

        for cycle in range(num_cycles):
            # Basin hopping perturbation
            t_pert = curr_p[:3] + self.rng.normal(0.0, 0.8, size=3)
            q_raw = curr_p[3:7] + self.rng.normal(0.0, 0.2, size=4)
            q_pert = q_raw / max(1e-6, np.linalg.norm(q_raw))
            phi_pert = curr_p[7:] + self.rng.uniform(-np.pi / 4, np.pi / 4, size=len(phi0))

            trial_p0 = pack_params(t_pert, q_pert, phi_pert)

            # Local Nelder-Mead / Powell minimization
            opt_res = minimize(
                objective,
                trial_p0,
                method="Powell",
                options={"maxiter": 35, "ftol": 1e-3},
            )
            trial_opt_p = opt_res.x
            trial_score = float(opt_res.fun)

            # Metropolis-Hastings acceptance
            delta_e = trial_score - curr_score
            accept = False
            if delta_e < 0.0:
                accept = True
            else:
                prob = np.exp(-delta_e / kb_t)
                if self.rng.random() < prob:
                    accept = True

            if accept:
                curr_p = trial_opt_p
                curr_score = trial_score

            if trial_score < best_score:
                best_p = trial_opt_p
                best_score = trial_score

        # Unpack best pose
        best_t, best_q, best_phi = unpack_params(best_p)
        best_q = best_q / max(1e-6, np.linalg.norm(best_q))
        final_coords = ligand.compute_conformation(best_t, best_q, best_phi)
        final_energy = self.scorer.evaluate(final_coords, ligand, receptor)

        return DockingResult(
            ligand_name=ligand.name,
            receptor_name=receptor.name,
            docked_coords=final_coords,
            translation=best_t,
            quaternion=best_q,
            dihedrals=best_phi,
            energy=final_energy,
            iterations=num_cycles,
        )
