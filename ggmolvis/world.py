import bpy
import molecularnodes as mn
import MDAnalysis as mda
import numpy as np
from pydantic import BaseModel, Field, validator, ValidationError
from typing import Tuple, List, Union

from .base import GGMolvisArtist
from . import SESSION
from .utils import convert_list_to_array

class WorldTransformation(BaseModel):
    coordinates: Union[Tuple[float, float, float], List[Tuple[float, float, float]]] = Field(
        ...,
        description="Static transformation (x, y, z) or a list/array of transformations for animation"
    )


    @validator('coordinates', pre=True)
    def convert_lists_to_array(cls, value):
        return convert_list_to_array(value)

    def apply_to(self, obj, frame: int = 0):
        raise NotImplementedError("This method must be implemented in the subclass")
    

    def get_transformation_for_frame(self, frame: int) -> Tuple[float, float, float]:
        """Retrieve the coordinates for a specific frame"""
        if self.coordinates.ndim == 2:
            index = frame % len(self.coordinates)
            return self.coordinates[index]
        elif self.coordinates.ndim == 1:
            return self.coordinates
        else:
            raise ValueError("Invalid transformation coordinates")
        
class Location(WorldTransformation):
    # The location can be either a static tuple or a list of tuples for animations
    coordinates: Union[Tuple[float, float, float],
                      List[Tuple[float, float, float]]] = Field(
        ...,
    alias="location",
    description="Zero point of the object either 3D coordinates"
                    "or a list/array of 3D coordinates for animation"
    )

    def apply_to(self, obj, frame: int = 0):
        """Apply location to the object, considering if it's static or animated"""
        if isinstance(self.coordinates, np.ndarray):
            # Static location
            obj.location = self.coordinates * self.world_scale
        elif isinstance(self.coordinates, list):
            # Animated location
            obj.location = self.get_transformation_for_frame(frame) * self.world_scale
        

class Rotation(WorldTransformation):
    # The rotation can be either a static tuple or a list of tuples for animations
    coordinates: Union[Tuple[float, float, float], List[Tuple[float, float, float]]] = Field(
        ...,
        alias = "rotation",
        description="Static rotation (roll, pitch, yaw) or a list/array of rotations for animation"
    )
    

    def apply_to(self, obj, frame: int = 0):
        """Apply rotation to the object, considering if it's static or animated"""
        if isinstance(self.coordinates, np.ndarray):
            # Static rotation
            obj.rotation_euler = self.coordinates
        elif isinstance(self.coordinates, list):
            # Animated rotation
            obj.rotation_euler = self.get_transformation_for_frame(frame)


class Scale(WorldTransformation):
    # The scale can be either a static tuple or a list of tuples for animations
    coordinates: Union[Tuple[float, float, float], List[Tuple[float, float, float]]] = Field(
        ..., 
        alias="scale",
        description="Static scale (x, y, z) or a list/array of scales for animation"
    )
    
    def apply_to(self, obj, frame: int = 0):
        """Apply scale to the object, considering if it's static or animated"""
        if isinstance(self.coordinates, np.ndarray):
            # Static scale
            obj.scale = self.coordinates
        elif isinstance(self.coordinates, list):
            # Animated scale
            obj.scale = self.get_transformation_for_frame(frame)


class World(GGMolvisArtist):
    location: Location
    rotation: Rotation
    scale: Scale
    
    def __init__(self,
                 location=None,
                 rotation=None,
                 scale=None):
        super().__init__()
        if location is None:
            location = (0.0, 0.0, 0.0)
        if rotation is None:
            rotation = (0.0, 0.0, 0.0)
        if scale is None:
            scale = (1.0, 1.0, 1.0)
        self.location = Location(location=location)
        self.rotation = Rotation(rotation=rotation)
        self.scale = Scale(scale=scale)
    
    def draw(self):
        """Implement the draw method"""
        # This method would typically apply transformations to objects in the session
        pass
    
    def apply_to(self, obj, frame: int = 0):
        """Apply the world transformations to the object"""
        self.location.apply_to(obj, frame)
        self.rotation.apply_to(obj, frame)
        self.scale.apply_to(obj, frame)