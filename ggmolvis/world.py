"""
World is the class that represents the world transformations of an object.
It contains the location, rotation, and scale of the object.

Classes
=======
.. autoclass:: WorldTransformation
    :members:
.. autoclass:: Location
    :members:
.. autoclass:: Rotation
    :members:
.. autoclass:: Scale
    :members:
.. autoclass:: World
    :members:
"""
import bpy
import molecularnodes as mn
import numpy as np
from pydantic import BaseModel, Field, field_validator
from typing import Tuple, List, Union

from .base import GGMolvisArtist
from .utils import convert_list_to_array, euler_to_quaternion

class WorldTransformation(BaseModel):
    coordinates: Union[List[float], List[List[float]]] = Field(
        ...,
        description="Static transformation (x, y, z) or a list/array of transformations for animation"
    )

    @field_validator('coordinates')
    def _convert_lists_to_array(cls, value):
        return convert_list_to_array(value)

    def _apply_to(self, obj, frame: int = 0):
        raise NotImplementedError("This method must be implemented in the subclass")
    

    def _get_transformation_for_frame(self, frame: int) -> np.ndarray:
        """Retrieve the coordinates for a specific frame"""
        if self.coordinates.ndim == 2:
            if frame >= len(self.coordinates):
                frame = -1
            return np.radians(self.coordinates[frame])
        elif self.coordinates.ndim == 1:
            return np.radians(self.coordinates)
        else:
            raise ValueError("Invalid transformation coordinates")

    def _set_coordinates(self, coordinates):
        self.coordinates = convert_list_to_array(coordinates)
        

class Location(WorldTransformation):
    # The location can be either a static tuple or a list of tuples for animations
    coordinates: Union[List[float], List[List[float]]] = Field(
        ...,
    alias="location",
    description="Zero point of the object either 3D coordinates"
                    "or a list/array of 3D coordinates for animation"
    )

    def _apply_to(self, obj, frame: int = 0):
        """Apply location to the object, considering if it's static or animated"""
        if isinstance(self.coordinates, np.ndarray):
            # Static location
            obj.location = self.coordinates
        elif isinstance(self.coordinates, list):
            # Animated location
            obj.location = self._get_transformation_for_frame(frame)
        

class Rotation(WorldTransformation):
    # The rotation can be either a static tuple or a list of tuples for animations
    coordinates: Union[List[float], List[List[float]]] = Field(
        ...,
        alias = "rotation",
        description="Static rotation (roll, pitch, yaw) or a list/array of rotations for animation"
    )
    
    def _apply_to(self, obj, frame: int = 0):
        """Apply rotation to the object, considering if it's static or animated"""
        if isinstance(self.coordinates, np.ndarray):
            # Static rotation to rad
            obj.rotation_euler = np.radians(self.coordinates)
        elif isinstance(self.coordinates, list):
            # Animated rotation
            obj.rotation_euler = self._get_transformation_for_frame(frame)
    
    @property
    def quaternion(self):
        """Return the quaternion representation of the rotation"""
        return euler_to_quaternion(self.coordinates)
    
    def from_quaternion(self, quaternion):
        """Set the rotation from a quaternion"""
        self.coordinates = quaternion.to_euler()


class Scale(WorldTransformation):
    # The scale can be either a static list or a list of list for animations
    coordinates: Union[List[float], List[List[float]]] = Field(
        ..., 
        alias="scale",
        description="Static scale (x, y, z) or a list/array of scales for animation"
    )
    
    def _apply_to(self, obj, frame: int = 0):
        """Apply scale to the object, considering if it's static or animated"""
        if isinstance(self.coordinates, np.ndarray):
            # Static scale
            obj.scale = self.coordinates
        elif isinstance(self.coordinates, list):
            # Animated scale
            obj.scale = self._get_transformation_for_frame(frame)


class World(GGMolvisArtist):
    def __init__(self,
                 location=None,
                 rotation=None,
                 scale=None):
        super().__init__()
        if location is None:
            location = [0.0, 0.0, 0.0]
        if rotation is None:
            rotation = [0.0, 0.0, 0.0]
        if scale is None:
            scale = [1.0, 1.0, 1.0]
        
        self._location = Location(location=location)
        self._rotation = Rotation(rotation=rotation)
        self._scale = Scale(scale=scale)

    def _update_frame(self, frame_number):
        """Not implemented in the World class"""
        # TODO: Is it necessary to implement this method?
    
    def _apply_to(self, obj, frame: int = 0):
        """Apply the world transformations to the object"""
        self.location._apply_to(obj, frame)
        self.rotation._apply_to(obj, frame)
        self.scale._apply_to(obj, frame)

    @property
    def matrix_world(self):
        """Return the 4D world matrix"""

        mat_world = np.zeros((4, 4))

        # convert rotation x, y, z to quaternion
        rot_mat = self.rotation.quaternion
        mat_world[:3, :3] = rot_mat.to_matrix()
        mat_world[:3, 3] = self.location.coordinates
        mat_world[-1, -1] = self.scale

    @property
    def location(self):
        return self._location

    @location.setter
    def location(self, value):
        self._location._set_coordinates(value)

    @property
    def rotation(self):
        return self._rotation
    
    @rotation.setter
    def rotation(self, value):
        self._rotation._set_coordinates(value)

    @property
    def scale(self):
        return self._scale
    
    @scale.setter
    def scale(self, value):
        self._scale._set_coordinates(value)