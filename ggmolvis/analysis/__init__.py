from MDAnalysis.analysis.base import AnalysisBase

from .visualizer import Visualizer
from .rmsd import RMSDVisualizer

def visualize(self):
    # Create a visualizer object based on the analysis class
    # e.g. if self is RMSD, then the visualizer object will be RMSDVisualizer
    visualizer = Visualizer.from_analysis(self)
    visualizer._visualize()
    return visualizer

AnalysisBase.visualize = visualize