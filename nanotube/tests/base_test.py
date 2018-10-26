import pytest
import mbuild as mb
import os
from pkg_resources import resource_filename

TESTFILE_DIR = resource_filename('nanotube', 'tests/test_molecules')


class BaseTest:

    @pytest.fixture(autouse=True)
    def initdir(self, tmpdir):
        tmpdir.chdir()

    @pytest.fixture
    def GraphenePore(self):
        from nanotube.nanotube import SWCNT
        return Nanotube(length=4,radius=.7)
