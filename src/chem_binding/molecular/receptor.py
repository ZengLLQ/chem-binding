"""Protein receptor binding pocket and clinically validated drug targets."""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from chem_binding.molecular.ligand import Atom, SmallMoleculeLigand


@dataclass
class ReceptorPocket:
    """Protein macromolecular binding site."""

    name: str
    atoms: list[Atom]

    def __post_init__(self) -> None:
        self.coords = np.array([a.initial_coords for a in self.atoms], dtype=np.float64)
        self.num_atoms = len(self.atoms)
        self.charges = np.array([a.charge for a in self.atoms], dtype=np.float64)
        self.vdw_radii = np.array([a.vdw_radius for a in self.atoms], dtype=np.float64)
        self.vdw_eps = np.array([a.vdw_well_depth for a in self.atoms], dtype=np.float64)
        self.is_donor = np.array([a.is_donor for a in self.atoms], dtype=bool)
        self.is_acceptor = np.array([a.is_acceptor for a in self.atoms], dtype=bool)
        self.is_hydrophobic = np.array([a.is_hydrophobic for a in self.atoms], dtype=bool)


def create_egfr_kinase_pocket(t790m_mutant: bool = False) -> ReceptorPocket:
    """Construct EGFR Tyrosine Kinase ATP-binding cleft (PDB: 1M17/2JIT).

    Parameters
    ----------
    t790m_mutant : bool
        If True, replaces gatekeeper Thr790 (hydroxyl, radius 1.55 A) with Met790
        (bulky methylthio sidechain, radius 1.85 A, causing steric hindrance to Gefitinib).
    """
    atoms: list[Atom] = [
        # Hinge region (Met793 backbone NH and carbonyl for crucial Gefitinib adenine-mimic H-bonds)
        Atom(0, "MET793_N", "N", np.array([2.5, 0.0, 0.0]), charge=-0.35, vdw_radius=1.55, vdw_well_depth=0.16, is_donor=True),
        Atom(1, "MET793_O", "O", np.array([2.8, 1.2, 0.0]), charge=-0.45, vdw_radius=1.52, vdw_well_depth=0.20, is_acceptor=True),
        # Catalytic Lys745 (Positive charge, salt-bridge interaction)
        Atom(2, "LYS745_NZ", "N", np.array([-3.0, 2.0, 1.5]), charge=+0.60, vdw_radius=1.55, vdw_well_depth=0.16, is_donor=True),
        # Catalytic Asp855 (Negative DFG-motif charge)
        Atom(3, "ASP855_OD1", "O", np.array([-2.5, -2.5, -1.0]), charge=-0.60, vdw_radius=1.52, vdw_well_depth=0.20, is_acceptor=True),
        # Hydrophobic pocket I (Leu718, Val726)
        Atom(4, "LEU718_CD1", "C", np.array([0.0, 3.5, -2.0]), charge=0.00, vdw_radius=1.70, vdw_well_depth=0.15, is_hydrophobic=True),
        Atom(5, "VAL726_CG1", "C", np.array([1.2, -3.2, 2.2]), charge=0.00, vdw_radius=1.70, vdw_well_depth=0.15, is_hydrophobic=True),
    ]

    # Gatekeeper Residue 790:
    if t790m_mutant:
        # Met790: Bulky, hydrophobic, protruding into the ATP cleft (steric resistance)
        atoms.append(
            Atom(6, "MET790_CE", "C", np.array([0.2, 0.1, -1.1]), charge=0.05, vdw_radius=1.85, vdw_well_depth=0.25, is_hydrophobic=True)
        )
        pocket_name = "EGFR Kinase (T790M Drug-Resistant Gatekeeper Mutant)"
    else:
        # Thr790: Compact hydroxyl, accommodating aniline ring of Gefitinib
        atoms.append(
            Atom(6, "THR790_OG1", "O", np.array([0.2, 0.1, -1.8]), charge=-0.40, vdw_radius=1.52, vdw_well_depth=0.18, is_acceptor=True)
        )
        pocket_name = "EGFR Kinase (Wild-Type Sensitive Gatekeeper)"

    return ReceptorPocket(name=pocket_name, atoms=atoms)


def create_gefitinib_ligand() -> SmallMoleculeLigand:
    """Create 4-anilinoquinazoline EGFR inhibitor Gefitinib (Iressa)."""
    # Core quinazoline scaffold + 3-chloro-4-fluoroaniline + morpholine propoxy chain
    atoms = [
        # Quinazoline N1 and N3 (Hinge acceptors)
        Atom(0, "N1", "N", np.array([1.2, 0.2, 0.0]), charge=-0.40, vdw_radius=1.55, vdw_well_depth=0.16, is_acceptor=True),
        Atom(1, "C2", "C", np.array([0.8, 1.3, 0.0]), charge=+0.15, vdw_radius=1.70, vdw_well_depth=0.15),
        Atom(2, "N3", "N", np.array([-0.5, 1.5, 0.0]), charge=-0.40, vdw_radius=1.55, vdw_well_depth=0.16, is_acceptor=True),
        Atom(3, "C4", "C", np.array([-1.2, 0.4, 0.0]), charge=+0.25, vdw_radius=1.70, vdw_well_depth=0.15),
        # 4-Aniline Nitrogen and ring (protrudes towards gatekeeper 790)
        Atom(4, "N4", "N", np.array([-0.8, -0.9, 0.0]), charge=-0.30, vdw_radius=1.55, vdw_well_depth=0.16, is_donor=True),
        Atom(5, "C1'", "C", np.array([0.0, -1.2, -0.9]), charge=+0.10, vdw_radius=1.70, vdw_well_depth=0.15, is_hydrophobic=True),
        Atom(6, "Cl3'", "Cl", np.array([-0.5, -2.5, -1.5]), charge=-0.10, vdw_radius=1.75, vdw_well_depth=0.25, is_hydrophobic=True),
        Atom(7, "F4'", "F", np.array([0.8, -1.8, -1.6]), charge=-0.20, vdw_radius=1.47, vdw_well_depth=0.10),
    ]

    # Rotatable bond between C4 and N4 (aniline pivot)
    torsions = [
        # TorsionalBond(atom1_idx=3, atom2_idx=4, moving_atoms=[5, 6, 7])
    ]

    return SmallMoleculeLigand(name="Gefitinib (Iressa)", atoms=atoms, torsional_bonds=torsions)
