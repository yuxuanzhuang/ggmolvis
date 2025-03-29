from .base import Visualizer
from ..ggmolvis import GGMolVis

class RMSDVisualizer(Visualizer):
    """Visualizer for RMSD analysis"""
    def _visualize(self):
        analysis = self.analysis
        mobile_atoms = analysis.mobile_atoms
        ggmv = GGMolVis()
        # create a trajectory object for the mobile atoms
        mobile_mol = ggmv.trajectory(mobile_atoms,
                                        style='cartoon',
                                        color='default',
                                        name='rmsd')
        self.mobile_mol = mobile_mol

        if not analysis.results:
            analysis.run()
        rmsd_results = analysis.results['rmsd'].T[2]
        mobile_mol.color.set_map(values=rmsd_results, cmap='coolwarm')

        return self
    
    @property
    def _view_location(self):
        """Get the view location for the camera"""
        return self.mobile_mol._view_location
    
    @property
    def _view_rotation(self):
        """Get the view rotation for the camera"""
        return self.mobile_mol._view_rotation
    
    @property
    def _camera_view_active(self):
        return self.mobile_mol._camera_view_active
    
    @_camera_view_active.setter
    def _camera_view_active(self, value):
        self.mobile_mol._camera_view_active = value
    
    def _set_camera_view(self):
        self.mobile_mol._set_camera_view()