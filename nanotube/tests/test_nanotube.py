import numpy as np
import mbuild as mb
import pytest

from nanotube.tests.base_test import BaseTest

class TestNanoTubeBuilder(BaseTest):
    """
    Unit Tests for Nanotube class functionality.
    """

    def test_save_dry(self, Nanotube):
        Nanotube.save(filename='dry_nanotube.gro')

    def test_positions(self, Nanotube):
        for C in Nanotube.children:
            assert(C.pos[2] >= 0)
            assert(C.pos[2] <= 4)
            radius = np.sqrt(C.pos[0]*C.pos[0]+C.pos[1]*C.pos[1])
            assert(radius <= .72)
            assert(radius >= .68)