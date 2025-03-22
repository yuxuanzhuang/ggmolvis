from ggmolvis.ggmolvis import GGMolVis
from ggmolvis.tests.data import PSF, DCD

from MDAnalysis import Universe
from ggmolvis.utils import MOL_AVAILABLE_STYLES, AVAILABLE_MATERIALS

import pytest
from pytest import fixture, mark

import bpy

@fixture(scope='module')
def universe():
    return Universe(PSF, DCD)

@fixture(scope='module')
def atomgroup(universe):
    return universe.select_atoms('protein')

@fixture(scope='module')
def ggmolvis(universe):
    return GGMolVis()

def test_render_scene(ggmolvis, atomgroup):
    mol = ggmolvis.trajectory(atomgroup)
    ggmolvis.render()

# add test for rendering
render_options = [
    {'frame': 0, },
    {'frame_range': (0, 4, 1)},
    {'track': True},
]

@mark.parametrize('render_options', render_options)
def test_render(ggmolvis, atomgroup, render_options):
    mol = ggmolvis.trajectory(atomgroup)
    ggmolvis.render(mol, **render_options,
                    resolution=(300, 200),
                    )

def test_render_error_frame_range(ggmolvis, atomgroup):
    mol = ggmolvis.trajectory(atomgroup)
    with pytest.raises(ValueError):
        ggmolvis.render(mol, frame=0, frame_range=(0, 4, 1))
    
    with pytest.raises(ValueError, match='frame_range must be'):
        ggmolvis.render(mol, frame_range=(0, 4))