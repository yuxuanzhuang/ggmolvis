from ggmolvis.ggmolvis import GGMolVis
from ggmolvis.tests.data import PSF, DCD

from MDAnalysis import Universe
from ggmolvis.utils import MOL_AVAILABLE_STYLES, AVAILABLE_MATERIALS

import pickle

import pytest
from pytest import fixture, mark
import bpy
from .utils import cleanup

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

@cleanup
def test_show_trajectory(ggmv, atomgroup):
    ggmv.trajectory(atomgroup)

@cleanup
def test_show_trajectory_with_universe(ggmv, universe):
    ggmv.trajectory(universe)

@cleanup
def test_show_trajectory_updating(ggmv, updating_atomgroup):
    mol = ggmv.trajectory(updating_atomgroup)
    assert mol.atomgroup.__class__.__name__ == 'UpdatingAtomGroup'

@cleanup
def test_trajectory_connected(ggmv, atomgroup):
    mol = ggmv.trajectory(atomgroup)
    # set frame to 1
    bpy.context.scene.frame_set(5)
    assert mol.frame == 5
    assert atomgroup.universe.trajectory.frame == 5
    assert mol.n_frames == atomgroup.universe.trajectory.n_frames

@cleanup
@mark.parametrize('material', AVAILABLE_MATERIALS)
def test_show_trajectory_with_material(ggmv, atomgroup, material):
    ggmv.trajectory(atomgroup, material=material)

@cleanup
@pytest.mark.parametrize('style', MOL_AVAILABLE_STYLES)
def test_show_trajectory_with_style(ggmv, atomgroup, style):
    # Mark specific styles as expected to fail
    if style in ['oxdna', 'density_surface', 'density_wire']:
        pytest.xfail(f"Style '{style}' is known to fail.")

    # Test the function with the given style
    ggmv.trajectory(atomgroup, style=style)

@cleanup
def test_persistence(tmpdir, ggmv, atomgroup):
    # Save the object to a file
    atm_gg = ggmv.trajectory(atomgroup)
    traj_uuid = atm_gg.trajectory.uuid
    
    filename = tmpdir.join('vis.pkl')
    with open (filename, 'wb') as f:
        pickle.dump(ggmv, f)
    
    # Load the object from the file
    with open(filename, 'rb') as f:
        vis2 = pickle.load(f)

    # Check if the objects are correctly reloaed
    assert vis2.trajectories[0].trajectory.uuid == traj_uuid

    filename = tmpdir.join('vis.blend').strpath
    bpy.ops.wm.save_as_mainfile(filepath=filename)
    atm2_gg = vis2.trajectory(atomgroup)
    assert len(vis2.trajectories) == 2
    bpy.ops.wm.open_mainfile(filepath=filename)
    assert vis2.trajectories[0].trajectory.uuid == traj_uuid
    assert len(vis2.trajectories) == 1