import bpy
import molecularnodes as mn
from molecularnodes.entities.trajectory import Trajectory
from molecularnodes.entities.trajectory.selections import Selection
import MDAnalysis as mda
import numpy as np
from typing import Union
from pydantic import BaseModel, Field, validator, ValidationError

from .base import GGMolvisArtist


class Property(GGMolvisArtist):
    """Class for the property of the visulizations."""
    timeseries: Union[np.ndarray, int] = 0

    def update_frame(self, frame):
        # TODO: Implement the update_frame method
        pass

class Color(Property):
    """Class for the color."""
    color: np.ndarray = np.array([1.0, 1.0, 1.0, 1.0])


class Material(Property):
    """Class for the material."""
    material: str = 'default'
