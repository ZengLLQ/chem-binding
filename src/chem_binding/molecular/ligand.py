"""Ligand atomistic representation, unit quaternions, and forward torsional kinematics."""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np


def quaternion_to_rotation_matrix(q: np.ndarray) -> np.ndarray:
    """Convert a unit quaternion q = [w, x, y, z] to a 3x3 orthogonal rotation matrix."""
    q_norm = q / max(1e-12, np.linalg.norm(q))
    w, x, y, z = q_norm

    return np.array([
        [1.0 - 2.0 * (y**2 + z**2), 2.0 * (x * y - z * w), 2.0 * (x * z + y * w)],
        [2.0 * (x * y + z * w), 1.0 - 2.0 * (x**2 + z**2), 2.0 * (y * z - x * w)],
        [2.0 * (x * z - y * w), 2.0 * (y * z + x * w), 1.0 - 2.0 * (x**2 + y**2)],
    ], dtype=np.float64)


def rodrigues_rotation(coords: np.ndarray, origin: np.ndarray, axis: np.ndarray, theta: float) -> np.ndarray:
    """Rotate points around an arbitrary line (origin + axis) by angle theta."""
    u = axis / max(1e-12, np.linalg.norm(axis))
    # Translate to origin
    v = coords - origin
    # Rodrigues formula: v_rot = v*cos(theta) + (u x v)*sin(theta) + u*(u . v)*(1 - cos(theta))
    cos_t = np.cos(theta)
    sin_t = np.sin(theta)
    dot_uv = np.sum(v * u, axis=1, keepdims=True)
    cross_uv = np.cross(u, v)

    v_rot = v * cos_t + cross_uv * sin_t + u * dot_uv * (1.0 - cos_t)
    return v_rot + origin


@dataclass
class Atom:
    """Atomistic definition for medicinal chemistry scoring."""

    index: int
    name: str
    element: str
    initial_coords: np.ndarray  # (3,)
    charge: float               # Partial charge e.g. Gasteiger (-1.0 to +1.0)
    vdw_radius: float          # Angstroms (e.g. C: 1.70, N: 1.55, O: 1.52, H: 1.20)
    vdw_well_depth: float      # kcal/mol (e.g. 0.15)
    is_donor: bool = False
    is_acceptor: bool = False
    is_hydrophobic: bool = False


@dataclass
class TorsionalBond:
    """A rotatable single bond partitioning the ligand into rigid branches."""

    atom1_idx: int
    atom2_idx: int
    moving_atoms: list[int]  # Atom indices affected by rotating this dihedral


class SmallMoleculeLigand:
    """Small molecule chemical structure with rigid and internal torsional DOFs."""

    def __init__(
        self,
        name: str,
        atoms: list[Atom],
        torsional_bonds: list[TorsionalBond] | None = None,
    ) -> None:
        self.name = name
        self.atoms = atoms
        self.num_atoms = len(atoms)
        self.torsional_bonds = torsional_bonds or []
        self.num_rotatable_bonds = len(self.torsional_bonds)

        # Baseline coordinates array (N, 3)
        self.base_coords = np.array([a.initial_coords for a in self.atoms], dtype=np.float64)

    def compute_conformation(
        self,
        translation: np.ndarray,      # (3,) vector
        quaternion: np.ndarray,       # (4,) [w, x, y, z]
        dihedrals: np.ndarray | None = None, # (K,) array of angles in radians
    ) -> np.ndarray:
        """Compute transformed 3D atom coordinates (N, 3).

        First applies rigid body rotation R(q) and translation t.
        Then rotates downstream atoms around respective dihedral axes.
        """
        # 1. Rigid body rotation + translation
        R = quaternion_to_rotation_matrix(quaternion)
        coords = (R @ self.base_coords.T).T + translation

        # 2. Internal rotatable bond dihedrals
        if dihedrals is not None and len(self.torsional_bonds) > 0:
            for b_idx, bond in enumerate(self.torsional_bonds):
                if b_idx < len(dihedrals):
                    theta = float(dihedrals[b_idx])
                    if abs(theta) > 1e-5:
                        p1 = coords[bond.atom1_idx]
                        p2 = coords[bond.atom2_idx]
                        axis = p2 - p1
                        moving = bond.moving_atoms
                        coords[moving] = rodrigues_rotation(coords[moving], p1, axis, theta)

        return coords
