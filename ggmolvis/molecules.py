import bpy
import molecularnodes as mn
from molecularnodes.entities.trajectory import Trajectory
from molecularnodes.entities.trajectory.selections import Selection
import MDAnalysis as mda
import numpy as np
from typing import Union
from pydantic import BaseModel, Field, validator, ValidationError

from .base import GGMolvisArtist
from .sceneobjects import SceneObject


class Molecule(SceneObject):
    """Class for a molecule."""
    def __init__(
             self,
             atomgroup: Union[mda.AtomGroup, mda.Universe],
             style: str = 'spheres',
             subframes: int = 0,
             name: str = 'atoms',
             location=None,
             rotation=None,
             scale=None
             ):
        """Show the molecule."""

        self.universe = atomgroup if isinstance(atomgroup, mda.Universe) else atomgroup.universe
        self.atomgroup = atomgroup if isinstance(atomgroup, mda.AtomGroup) else self.universe.atoms

        self.style = style
        self.subframes = subframes

        traj = Trajectory(self.universe)
        self.trajectory = traj

        super().__init__(location=location, rotation=rotation, scale=scale)


    def draw(self):
        self.trajectory.create_object(name=self.name, style=self.style, subframes=self.subframes)
        self.trajectory.add_selection_from_atomgroup(self.atomgroup, name=self.name)


    @property
    def object(self):
        return self.trajectory.object