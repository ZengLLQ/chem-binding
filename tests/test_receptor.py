"""Unit tests for protein receptor target definitions."""

from chem_binding.molecular.receptor import (
    create_egfr_kinase_pocket,
    create_gefitinib_ligand,
)


def test_egfr_pockets() -> None:
    pocket_wt = create_egfr_kinase_pocket(t790m_mutant=False)
    pocket_mut = create_egfr_kinase_pocket(t790m_mutant=True)

    assert pocket_wt.num_atoms == 7
    assert pocket_mut.num_atoms == 7

    # Check gatekeeper residue in WT vs Mutant
    names_wt = [a.name for a in pocket_wt.atoms]
    names_mut = [a.name for a in pocket_mut.atoms]

    assert "THR790_OG1" in names_wt
    assert "MET790_CE" in names_mut


def test_gefitinib_structure() -> None:
    gef = create_gefitinib_ligand()
    assert gef.num_atoms == 8
    assert gef.name.startswith("Gefitinib")
