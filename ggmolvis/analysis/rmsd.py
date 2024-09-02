from .visualizer import Visualizer
from ..molecules import Molecule
from ..ggmolvis import GGMolVis

class RMSDVisualizer(Visualizer):
    """Visualizer for RMSD analysis"""
    def _visualize(self):
        analysis = self.analysis
        mobile_atoms = analysis.mobile_atoms
        ggmv = GGMolVis()

        # create a molecule object for the mobile atoms
        ggmv.molecule(mobile_atoms, style='spheres', name='mobile')

        # set the color of the mobile atoms based on the RMSD values
        # TODO: implement this