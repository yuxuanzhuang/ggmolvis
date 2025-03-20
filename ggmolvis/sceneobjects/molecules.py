"""
This module contains the Molecule class, which is a subclass of SceneObject.

Classes
=======
.. autoclass:: Molecule
    :members:
"""
import bpy
import molecularnodes as mn
from molecularnodes.entities.trajectory import Trajectory
import MDAnalysis as mda
import numpy as np
from typing import Union

from .base import SceneObject
from ..utils.node import set_selection
from ..properties import MoleculeColor, MoleculeMaterial, MoleculeStyle


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
        material="default",
        style: str = "spheres",
    ):
        """Show the molecule."""
        if isinstance(atomgroup, mda.Universe):
            self.atomgroup = atomgroup.atoms
            self.universe = atomgroup
        elif isinstance(atomgroup, mda.AtomGroup):
            self.atomgroup = atomgroup
            self.universe = atomgroup.universe
        else:
            raise TypeError("atomgroup must be an AtomGroup or Universe")

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

    def _init_material(self, material="default", color=None):
        self._material = MoleculeMaterial(self, material, color=color)

    def _init_color(self, color="black"):
        self._color = MoleculeColor(self, color)
    
    def _create_object(self):
        traj = Trajectory(self.universe)
        self.trajectory = traj
        self.trajectory.create_object(
            name=self.name, style=self._style_name
        )
        self.trajectory.subframes = self.subframes
        self.trajectory.average = self.average
        self.trajectory.add_selection_from_atomgroup(
            self.atomgroup, name=self.name
        )
        # only render the selection
        set_selection(self.object, self.name)
        return self.object
    
    @property
    def object(self):
        return self.trajectory.object

    
