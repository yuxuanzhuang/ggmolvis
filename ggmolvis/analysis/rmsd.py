from .base import Visualizer
from ..ggmolvis import GGMolVis

class RMSDVisualizer(Visualizer):
    """Visualizer for RMSD analysis"""
    def _visualize(self):
        analysis = self.analysis
        mobile_atoms = analysis.mobile_atoms
        ggmv = GGMolVis()
        # create a molecule object for the mobile atoms
        mobile_molecule = ggmv.molecule(mobile_atoms,
                                        style='cartoon',
                                        color='default',
                                        name='rmsd')

        self._camera_world = mobile_molecule.camera_world

        if not analysis.results:
            analysis.run()
        rmsd_results = analysis.results['rmsd'].T[2]
        mobile_molecule.color.set_map(values=rmsd_results, cmap='coolwarm')

        return ggmv