"""Unit tests for ChemBinding CLI."""

from chem_binding.cli import main


def test_cli_dock() -> None:
    code = main(["dock", "--cycles", "3"])
    assert code == 0


def test_cli_screen_mutations() -> None:
    code = main(["screen-mutations", "--cycles", "3"])
    assert code == 0
