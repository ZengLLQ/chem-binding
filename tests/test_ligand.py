"""Unit tests for ligand kinematics and quaternion transformations."""

import numpy as np
from chem_binding.molecular.ligand import (
    Atom,
    SmallMoleculeLigand,
    quaternion_to_rotation_matrix,
)


def test_quaternion_identity() -> None:
    q_id = np.array([1.0, 0.0, 0.0, 0.0])
    R = quaternion_to_rotation_matrix(q_id)
    assert np.allclose(R, np.eye(3))


def test_quaternion_rotation_90z() -> None:
    # 90 degrees around Z axis: w = cos(45) = 1/sqrt(2), z = sin(45) = 1/sqrt(2)
    s = 1.0 / np.sqrt(2.0)
    q_z90 = np.array([s, 0.0, 0.0, s])
    R = quaternion_to_rotation_matrix(q_z90)

    # Point [1, 0, 0] rotated 90 deg around z becomes [0, 1, 0]
    p = np.array([1.0, 0.0, 0.0])
    p_rot = R @ p
    assert np.allclose(p_rot, [0.0, 1.0, 0.0], atol=1e-5)


def test_ligand_conformation_translation() -> None:
    atoms = [
        Atom(0, "C1", "C", np.array([0.0, 0.0, 0.0]), 0.0, 1.7, 0.15),
        Atom(1, "C2", "C", np.array([1.5, 0.0, 0.0]), 0.0, 1.7, 0.15),
    ]
    lig = SmallMoleculeLigand("TestLig", atoms)
    t = np.array([5.0, 2.0, -1.0])
    q = np.array([1.0, 0.0, 0.0, 0.0])

    coords = lig.compute_conformation(t, q)
    assert np.allclose(coords[0], [5.0, 2.0, -1.0])
    assert np.allclose(coords[1], [6.5, 2.0, -1.0])
