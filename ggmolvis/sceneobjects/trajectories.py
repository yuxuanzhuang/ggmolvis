"""
This module contains the Trajectory class, which is a subclass of SceneObject.

Classes
=======
.. autoclass:: Trajectory
    :members:
"""
import bpy
import molecularnodes as mn
from molecularnodes.entities.trajectory import Trajectory as mn_Trajectory
import MDAnalysis as mda
import numpy as np
from typing import Union

from .base import SceneObject
from ..utils.node import set_selection
from ..properties import TrajectoryColor, TrajectoryMaterial, TrajectoryStyle


class Trajectory(SceneObject):
    """Class for a trajectory.

    Displays a trajectory in the scene with the given properties.

    Properties
    ----------
    object : bpy.types.Object
        The object in the scene.
    atomgroup : MDAnalysis.AtomGroup
        The AtomGroup to display.
    name : str
        The name of the object.
    trajectory:
        The MN Trajectory object.
    universe:
        The MDAnalysis Universe object.
    """
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
        if isinstance(atomgroup, mda.Universe):
            self._atomgroup = atomgroup.atoms
            self._universe = atomgroup
        elif isinstance(atomgroup, mda.AtomGroup):
            self._atomgroup = atomgroup
            self._universe = atomgroup.universe
        else:
            raise TypeError("atomgroup must be an AtomGroup or Universe")

        # trajectory need style name to create object
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
        self._style = TrajectoryStyle(self, style)

    def _init_material(self, material="default"):
        self._material = TrajectoryMaterial(self, material)

    def _init_color(self, color="black"):
        self._color = TrajectoryColor(self, color)
    
    def _create_object(self):
        traj = mn_Trajectory(self.universe)
        self._trajectory = traj
        self._trajectory.create_object(
            name=self.name, style=self._style_name
        )
        self._trajectory.subframes = self.subframes
        self._trajectory.average = self.average
        self._trajectory.add_selection_from_atomgroup(
            self.atomgroup, name=self.name
        )
        # only render the selection
        set_selection(self.object, self.name)
        return self.object
    
    @property
    def object(self):
        return self._trajectory.object

    @property
    def universe(self):
        return self._universe
    
    @property
    def atomgroup(self):
        return self._atomgroup
    
    @property
    def trajectory(self):
        return self._trajectory