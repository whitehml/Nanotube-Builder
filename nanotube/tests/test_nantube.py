import numpy as np
import mbuild as mb
import pytest

from nanotube.tests.base_test import BaseTest
import SWCNT as nt


class TestNanoTubeBuilder(BaseTest):
    """
    Unit Tests for Nanotube class functionality.
    """

    # example from porebuilder
    # def test_save_dry(self, GraphenePore):
    #     GraphenePore.save(filename='dry_pore.gro')