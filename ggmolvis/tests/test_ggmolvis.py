import pickle
from ggmolvis import GGMolVis
from ggmolvis.world import World
from ggmolvis.camera import Camera

def test_initialization(ggmv):
    # Check if the global camera and world are set up correctly
    assert isinstance(ggmv.global_camera, Camera)
    assert isinstance(ggmv.global_world, World)

    # Check if the artists dictionary is populated with default values
    for key in ['molecules', 'shapes', 'texts', 'cameras', 'lights', 'worlds']:
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
    with open (filename, 'wb') as f:
        pickle.dump(ggmv, f)
    
    # Load the object from the file
    with open(filename, 'rb') as f:
        vis2 = pickle.load(f)

    # Check if the loaded object is the same as the original object
    assert ggmv is vis2