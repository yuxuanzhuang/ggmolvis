from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.core.groups import GroupBase

from .base import Visualizer
from .rmsd import RMSDVisualizer
from ..ggmolvis import GGMolVis

def _visualize_analysis(self,
                        render=True,
                        **kwargs):
    """
    Create a visualizer object based on the analysis class
    e.g. if self is RMSD, then the visualizer object will be RMSDVisualizer
    """
    ggmv = GGMolVis()
    visualizer = Visualizer.from_analysis(self, **kwargs)
    visualizer._visualize()
    if render:
        ggmv.render(visualizer)
    return visualizer

AnalysisBase.visualize = _visualize_analysis

def _visualize_atom(self,
                    render=True,
                    **kwargs):
    """
    Visualize atoms function
    """
    ggmv = GGMolVis()
    # create a molecule object for the mobile atoms
    mol_vis = ggmv.molecule(self.atoms, **kwargs)
    if render:
        ggmv.render(mol_vis)
    return mol_vis

GroupBase.visualize = _visualize_atom