from ggmolvis.ggmolvis import GGMolVis
from ggmolvis.tests.data import PSF, DCD

from MDAnalysis import Universe
from ggmolvis.utils import MOL_AVAILABLE_STYLES, AVAILABLE_MATERIALS

import pytest
from pytest import fixture, mark
from .utils import cleanup

import bpy

@fixture(scope='module')
def universe():
    return Universe(PSF, DCD)

@fixture(scope='module')
def atomgroup(universe):
    return universe.select_atoms('protein')

@cleanup
def test_render_scene(ggmv, atomgroup):
    mol = ggmv.trajectory(atomgroup)
    ggmv.render()

# add test for rendering
render_options = [
    {'frame': 0, },
    {'frame_range': (0, 4, 1)},
    {'track': True},
]

@cleanup
@mark.parametrize('render_options', render_options)
def test_render(ggmv, atomgroup, render_options):
    mol = ggmv.trajectory(atomgroup)
    ggmv.render(mol, **render_options,
                    resolution=(300, 200),
                    )

@cleanup
def test_render_error_frame_range(ggmv, atomgroup):
    mol = ggmv.trajectory(atomgroup)
    with pytest.raises(ValueError):
        ggmv.render(mol, frame=0, frame_range=(0, 4, 1))
    
    with pytest.raises(ValueError, match='frame_range must be'):
        ggmv.render(mol, frame_range=(0, 4))
