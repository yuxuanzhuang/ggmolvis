import bpy
from MDAnalysis.analysis.base import AnalysisBase
from ..base import GGMolvisArtist

class Visualizer(GGMolvisArtist):
    def __init__(self, analysis: AnalysisBase, **kwargs):
        self.analysis = analysis
        super().__init__(**kwargs)

    def _visualize(self):
        pass

    @classmethod
    def from_analysis(cls, analysis: AnalysisBase, **kwargs):
        visualizer_name = AVAILABLE_VISUALIZERS.get(analysis.__class__.__name__)
        visualizer = get_visualizer(analysis, **kwargs)
        return visualizer

    def _update_frame(self, frame):
        # TODO: Implement this method
        pass

    def render(self, **kwargs):
        # set the render movie end frame to the number of frames in the trajectory
        bpy.context.scene.frame_end = self.analysis.atomgroup.universe.trajectory.n_frames
        self.ggmolvis.render(object=self, **kwargs)
    
    def _set_camera_view(self):
        raise NotImplementedError("This method is only available in the subclass")


AVAILABLE_VISUALIZERS = {}

def register_visualizer(analysis_name: str):
    def decorator(visualizer_class):
        AVAILABLE_VISUALIZERS[analysis_name] = visualizer_class
        return visualizer_class
    return decorator

@register_visualizer('RMSD')
def _register_rmsd_visualizer(analysis, **kwargs):
    from .rmsd import RMSDVisualizer
    return RMSDVisualizer(analysis, **kwargs)

def get_visualizer(analysis: AnalysisBase, **kwargs):
    visualizer_name = AVAILABLE_VISUALIZERS.get(analysis.__class__.__name__)
    visualizer = visualizer_name(analysis, **kwargs)
    return visualizer
