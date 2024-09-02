import bpy
from abc import ABC, abstractmethod

import molecularnodes as mn
from molecularnodes.entities.trajectory import Trajectory
from molecularnodes.entities.trajectory.selections import Selection
import MDAnalysis as mda
import numpy as np
from typing import Tuple, List, Union
from pydantic import BaseModel, Field, validator, ValidationError

from .base import GGMolvisArtist
from .world import World
from .camera import Camera
from .properties import Color, Material
from .utils import convert_list_to_array, look_at

class SceneObject(GGMolvisArtist):
    """Class for the scene object.
    This class is the parent class for all the objects in the scene.
    Access the blender object using the object property, `self.object`.
    The name might be different from the initial name, as Blender might
    append a number to the name if the name already exists.
    """

    def __init__(self, name=None, location=None, rotation=None, scale=None):
        self.world = World(location=location, rotation=rotation, scale=scale)
        super().__init__(name=name)

        self.create_object()
        self.set_color()
        self.set_material()
        self.set_camera()
        self.draw()
        self.update_frame(bpy.context.scene.frame_current)

    def set_color(self):
        self.color = Color()

    def set_material(self):
        self.material = Material()

    def set_camera(self):
        self.camera = Camera(name=f"{self.name}_camera")
        size_obj_xyz = np.array(self.object.dimensions)
        # center of the object
        center_xyz = np.zeros(3)
        bbox = self.object.bound_box
        for v in bbox:
            center_xyz += np.array(v)
        center_xyz /= 8
        camera_center = center_xyz.copy()
        camera_center[1] = camera_center[1] - size_obj_xyz[1] * 3
        camera_center[2] = camera_center[2] + size_obj_xyz[2] * 1.3

        rotation_camera = look_at(camera_position=camera_center,
                                  target_position=center_xyz)

        #self.camera.world.rotation.set_coordinates(rotation_camera)
        self.camera.object.matrix_world = rotation_camera
        self.camera.world.location.set_coordinates(camera_center)

    def update_frame(self, frame):
        object = self.object
        self.world.apply_to(object, frame)

        camera = self.camera
        camera.update_frame(frame)
    
    @abstractmethod
    def create_object(self):
        raise NotImplementedError("This method is only available in the subclass")
    
    @abstractmethod
    def draw(self):
        """Draw the object"""
        # TODO: What should be done here?
        raise NotImplementedError("This method is only available in the subclass")

    @property
    def object(self):
        return bpy.data.objects[self.name]
    
    def render(self):
        self.camera.render()

class Text(SceneObject):
    """Class for the text."""
    # TODO: Implement the text class
    text: str = 'text'

