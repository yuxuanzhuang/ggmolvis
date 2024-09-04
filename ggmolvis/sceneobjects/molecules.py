import bpy
import molecularnodes as mn
from molecularnodes.entities.trajectory import Trajectory
import MDAnalysis as mda
import numpy as np
from typing import Union

from .base import SceneObject
from ..utils.node import set_selection
from ..properties import Color, Material, MoleculeStyle



class Molecule(SceneObject):
    """Class for a molecule."""
    def __init__(
        self,
        atomgroup: Union[mda.AtomGroup, mda.Universe],
        name: str = "atoms",
        location=None,
        rotation=None,
        scale=None,
        color="black",
        material="MN Default",
        style: str = "spheres",
    ):
        """Show the molecule."""

        self.universe = atomgroup if isinstance(atomgroup, mda.Universe) else atomgroup.universe
        self.atomgroup = atomgroup if isinstance(atomgroup, mda.AtomGroup) else self.universe.atoms

        # molecules need style name to create object
        self._style_name = style

        super().__init__(
            name=name,
            location=location,
            rotation=rotation,
            scale=scale,
            color=color,
            material=material,
            style=style,
        )

    def _init_style(self, style="default"):
        self._style = MoleculeStyle(self, style)

    def _create_object(self):
        traj = Trajectory(self.universe)
        self.trajectory = traj
        self.trajectory.create_object(
            name=self.name, style=self._style_name, subframes=self.subframes
        )
        self.trajectory.add_selection_from_atomgroup(
            self.atomgroup, name=self.name
        )
        # only render the selection
        set_selection(self.object, self.name)
    
    def draw(self):
        pass

    @property
    def object(self):
        return self.trajectory.object
