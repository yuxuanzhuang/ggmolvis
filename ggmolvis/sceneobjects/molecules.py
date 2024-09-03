import bpy
import molecularnodes as mn
from molecularnodes.entities.trajectory import Trajectory
from molecularnodes.blender.nodes import styles_mapping
import MDAnalysis as mda
import numpy as np
from typing import Union

from .sceneobjects import SceneObject
from ..utils.node import extract_mn_node, set_selection, swap_style


class Molecule(SceneObject):
    """Class for a molecule."""
    def __init__(
             self,
             atomgroup: Union[mda.AtomGroup, mda.Universe],
             style: str = 'spheres',
             name: str = 'atoms',
             location=None,
             rotation=None,
             scale=None,
             color='black',
             material='MN Default'
             ):
        """Show the molecule."""

        self.universe = atomgroup if isinstance(atomgroup, mda.Universe) else atomgroup.universe
        self.atomgroup = atomgroup if isinstance(atomgroup, mda.AtomGroup) else self.universe.atoms

        self._style = style


        super().__init__(name=name,
                         location=location,
                         rotation=rotation,
                         scale=scale,
                         color=color,
                         material=material)


    def create_object(self):
        traj = Trajectory(self.universe)
        self.trajectory = traj
        self.trajectory.create_object(name=self.name, style=self.style, subframes=self.subframes)
        self.trajectory.add_selection_from_atomgroup(self.atomgroup, name=self.name)
        # only render the selection
        set_selection(self.object, self.name)
    
    def draw(self):
        pass

    @property
    def object(self):
        return self.trajectory.object
    
    @property
    def style(self):
        return self._style
    
    @style.setter
    def style(self, value):
        if value not in styles_mapping:
            raise ValueError(f"Style {value} is not supported."
                             f"Supported styles are {styles_mapping.keys()}")
        self._style = value
        swap_style(self.object, styles_mapping[value])
