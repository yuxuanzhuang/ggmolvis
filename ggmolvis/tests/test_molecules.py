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
def updating_atomgroup(universe):
    return universe.select_atoms('resid 127 40', updating=True)

@fixture(scope='module')
def ggmolvis(universe):
    return GGMolVis()

def test_show_trajectory(ggmolvis, atomgroup):
    ggmolvis.trajectory(atomgroup)

def test_show_trajectory_with_universe(ggmolvis, universe):
    ggmolvis.trajectory(universe)

def test_show_trajectory_updating(ggmolvis, updating_atomgroup):
    mol = ggmolvis.trajectory(updating_atomgroup)
    assert mol.atomgroup.__class__.__name__ == 'UpdatingAtomGroup'

def test_trajectory_connected(ggmolvis, atomgroup):
    ggmolvis.trajectory(atomgroup)
    # set frame to 1
    bpy.context.scene.frame_set(5)
    assert atomgroup.universe.trajectory.frame == 5


@mark.parametrize('material', AVAILABLE_MATERIALS)
def test_show_trajectory_with_material(ggmolvis, atomgroup, material):
    ggmolvis.trajectory(atomgroup, material=material)


@pytest.mark.parametrize('style', MOL_AVAILABLE_STYLES)
def test_show_trajectory_with_style(ggmolvis, atomgroup, style):
    # Mark specific styles as expected to fail
    if style in ['oxdna', 'density_surface', 'density_wire']:
        pytest.xfail(f"Style '{style}' is known to fail.")

    # Test the function with the given style
    ggmolvis.trajectory(atomgroup, style=style)

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