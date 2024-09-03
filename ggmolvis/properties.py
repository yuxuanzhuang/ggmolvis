import bpy
from abc import abstractmethod

import molecularnodes as mn
from molecularnodes.entities.trajectory import Trajectory
from molecularnodes.entities.trajectory.selections import Selection
from molecularnodes.blender.material import materials
import MDAnalysis as mda
import numpy as np
from typing import Union
from pydantic import BaseModel, Field, validator, ValidationError

from .base import GGMolvisArtist
from .utils.node import set_material


class Property(GGMolvisArtist):
    """Class for the property of the visulizations."""
    def __init__(self, name=None, **kwargs):
        super().__init__(name=name, **kwargs)
        self.set_property(self.name)

    @abstractmethod
    def set_property(self, **kwargs):
        """Set the property of the visualization."""
        raise NotImplementedError("This method must be implemented in the subclass")

    def apply_to(self, obj, frame: int = 0):
        raise NotImplementedError("This method must be implemented in the subclass")


    def update_frame(self, frame):
        # TODO: Implement the update_frame method
        pass

class Color(Property):
    """Class for the color."""
    def set_property(self, color_name="black"):
        if color_name is None:
            color_name = "black"
        self.color_name = color_name

    def set_color(self, color_name):
        self.set_property(color_name)

    def apply_to(self, obj, frame: int = 0):
        pass

AVAILABLE_MATERIALS = materials
AVAILABLE_MATERIALS.append("Backdrop")


class Material(Property):
    """Class for the material."""
    def set_property(self, material_name="Backdrop"):
        if material_name is None:
            material_name = "Backdrop"

        if material_name not in AVAILABLE_MATERIALS:
            raise ValueError(f"Material {material_name} is not available."
                           f"Available materials are {AVAILABLE_MATERIALS}")
        
        # create a material copy
        default_material = bpy.data.materials[material_name]
        self.material = default_material.copy()
        self.material.name = f"{material_name}_{self.name}"

        self.material_name = self.material.name

    def set_material(self, material_name):
        self.set_property(material_name)

    def apply_to(self, obj, frame: int = 0):
        """Apply material to the object."""

        set_material(obj, self.material_name)