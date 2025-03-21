"""
Global pytest fixtures
"""

# Use this file if you need to share any fixtures
# across multiple modules
# More information at
# https://docs.pytest.org/en/stable/how-to/fixtures.html#scope-sharing-fixtures-across-classes-modules-packages-or-session

import os
import pytest

@pytest.fixture
def ggmv():
    from ggmolvis import GGMolVis
    return GGMolVis()


def pytest_runtest_setup(item):
    mark = item.get_closest_marker("xslow")
    if mark is not None:
        try:
            v = int(os.environ.get('GGMOLVIS_XSLOW', '0'))
        except ValueError:
            v = False
        if not v:
            pytest.skip("very slow test; "
                        "set environment variable GGMOLVIS_XSLOW=1 to run it")
