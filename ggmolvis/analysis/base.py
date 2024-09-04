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

    @property
    def camera(self):
        return self._camera
    
    def render(self):
        self.camera.render()


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
