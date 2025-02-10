from ggmolvis.ggmolvis import GGMolVis
from ggmolvis.tests.data import PSF, DCD

from MDAnalysis import Universe
from ggmolvis.utils import MOL_AVAILABLE_STYLES, AVAILABLE_MATERIALS

import pickle

import pytest
from pytest import fixture, mark
import bpy
from .utils import clean_up

@fixture(scope='module')
def universe():
    return Universe(PSF, DCD)

@fixture(scope='module')
def atomgroup(universe):
    return universe.select_atoms('protein')

@fixture(scope='module')
def updating_atomgroup(universe):
    return universe.select_atoms('resid 127 40', updating=True)


def test_show_molecule(ggmv, atomgroup):
    ggmv.molecule(atomgroup)

def test_show_molecule_updating(ggmv, updating_atomgroup):
    ggmv.molecule(updating_atomgroup)

@mark.parametrize('material', AVAILABLE_MATERIALS)
def test_show_molecule_with_material(ggmv, atomgroup, material):
    ggmv.molecule(atomgroup, material=material)

@pytest.mark.parametrize('style', MOL_AVAILABLE_STYLES)
def test_show_molecule_with_style(ggmv, atomgroup, style):
    # Mark specific styles as expected to fail
    if style in ['oxdna', 'density_surface', 'density_wire']:
        pytest.xfail(f"Style '{style}' is known to fail.")

    # Test the function with the given style
    ggmv.molecule(atomgroup, style=style)

def test_persistence(tmpdir, ggmv, atomgroup):
    clean_up(ggmv)
    # Save the object to a file
    atm_gg = ggmv.molecule(atomgroup)
    traj_uuid = atm_gg.trajectory.uuid
    
    filename = tmpdir.join('vis.pkl')
    with open (filename, 'wb') as f:
        pickle.dump(ggmv, f)
    
    # Load the object from the file
    with open(filename, 'rb') as f:
        vis2 = pickle.load(f)

    # Check if the objects are correctly reloaed
    assert vis2.molecules[0].trajectory.uuid == traj_uuid

    filename = tmpdir.join('vis.blend').strpath
    bpy.ops.wm.save_as_mainfile(filepath=filename)
    atm2_gg = vis2.molecule(atomgroup)
    assert len(vis2.molecules) == 2
    bpy.ops.wm.open_mainfile(filepath=filename)
    assert vis2.molecules[0].trajectory.uuid == traj_uuid
    assert len(vis2.molecules) == 1

