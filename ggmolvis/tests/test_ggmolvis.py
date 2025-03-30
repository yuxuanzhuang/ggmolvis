import pickle
import os
from unittest.mock import patch
import pytest
import importlib

import ggmolvis
from ggmolvis import GGMolVis
from ggmolvis.world import World
from ggmolvis.camera import Camera

@pytest.mark.skip(reason="Reload the module leads to cleanup issues")
def test_create_blend_file():
    with patch("os.makedirs") as mock_createdir, \
        patch("uuid.uuid4") as mock_uuid:
        mock_uuid.return_value = 'test-ggmolvis-uuid'

        mock_createdir.side_effect = Exception("cannot access temp dir")
        importlib.reload(ggmolvis)
        # Test if the blend file was created successfully
        assert os.path.exists('test-ggmolvis-uuid.blend')
        # remove the created file for cleanup
        os.remove('test-ggmolvis-uuid.blend')


def test_initialization(ggmv):
    # Check if the camera and world are set up correctly
    assert isinstance(ggmv.camera, Camera)
    assert isinstance(ggmv.global_world, World)

    # Check if the artists dictionary is populated with default values
    for key in ['trajectories', 'shapes', 'texts', 'lights', 'worlds']:
        assert key in ggmv._artists_dict

    # Check if the scene is set up correctly
    assert ggmv._subframes == 0

def test_singleton_pattern(ggmv):
    # Create second instance
    vis2 = GGMolVis()
    assert ggmv is vis2

def test_persistence(tmpdir, ggmv):
    # Save the object to a file
    filename = tmpdir.join('vis.pkl')
    with open(filename, 'wb') as f:
        pickle.dump(ggmv, f)

    # Load the object from the file
    with open(filename, 'rb') as f:
        vis2 = pickle.load(f)

    # Check if the loaded object is equivalent to the original object
    assert type(ggmv) == type(vis2)  # Ensure they are the same class
    assert ggmv.__dict__ == vis2.__dict__  # Check attributes for equality