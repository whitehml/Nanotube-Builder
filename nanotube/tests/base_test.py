import pytest
import mbuild as mb
import os
from pkg_resources import resource_filename

TESTFILE_DIR = resource_filename('nanotube', 'tests/test_molecules')


class BaseTest:

    @pytest.fixture(autouse=True)
    def initdir(self, tmpdir):
        tmpdir.chdir()

    # example to follow
    # @pytest.fixture
    # def GraphenePore(self):
    #     from porebuilder.porebuilder import GraphenePore
    #     return GraphenePore(pore_depth=3, side_dim=3, n_sheets=3, pore_width=1)
