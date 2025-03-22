import unittest
from ggmolvis import GGMolVis
from ggmolvis.world import World
from ggmolvis.camera import Camera

class TestGGMolVis(unittest.TestCase):

    def test_singleton_pattern(self):
        # Create first instance
        vis1 = GGMolVis()
        self.assertIsNotNone(vis1)
        self.assertIsInstance(vis1, GGMolVis)

        # Create second instance
        vis2 = GGMolVis()
        self.assertIs(vis1, vis2)

    def test_initialization(self):
        # Create instance
        vis = GGMolVis()

        # Check if the global camera and world are set up correctly
        self.assertIsInstance(vis.camera, Camera)
        self.assertIsInstance(vis.global_world, World)

        # Check if the artists dictionary is populated with default values
        self.assertIn('trajectories', vis._artists_dict)
        self.assertIn('shapes', vis._artists_dict)
        self.assertIn('texts', vis._artists_dict)
        self.assertIn('lights', vis._artists_dict)
        self.assertIn('worlds', vis._artists_dict)

        # Check if the scene is set up correctly
        self.assertEqual(vis._subframes, 0)

if __name__ == '__main__':
    unittest.main()