import bpy
from abc import abstractmethod

import numpy as np
from typing import Tuple, List, Union

from ..base import GGMolvisArtist
from ..world import World
from ..camera import Camera
from ..properties import Color, Material, Style
from ..utils import look_at


class SceneObject(GGMolvisArtist):
    """Class for the scene object.
    This class is the parent class for all the objects in the scene.
    Access the blender object using the object property, `self.object`.
    The name might be different from the initial name, as Blender might
    append a number to the name if the name already exists.


    """

    def __init__(
        self,
        name=None,
        location=None,
        rotation=None,
        scale=None,
        color="black",
        material="Backdrop",
        style="default",
    ):
        self.world = World(location=location, rotation=rotation, scale=scale)
        super().__init__()
        self._name = self.__class__.__name__

        self._subframes = 0

        self._create_object()
        self._init_color(color)
        self._init_material(material)
        self._init_style(style)
        self._init_camera()
        self.draw()
        self._update_frame(bpy.context.scene.frame_current)

    def _init_color(self, color="black"):
        self._color = Color(self, color)

    def _init_material(self, material="backdrop"):
        self._material = Material(self, material)

    def _init_style(self, style="default"):
        self._style = Style(self, style)

    def _init_camera(self):
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

        rotation_camera = look_at(
            camera_position=camera_center, target_position=center_xyz
        )

        # self.camera.world.rotation.set_coordinates(rotation_camera)
        self.camera.object.matrix_world = rotation_camera
        self.camera.world.location.set_coordinates(camera_center)

    def _update_frame(self, frame):
        object = self.object
        self.material.apply_to(object, frame)
        self.color.apply_to(object, frame)
        self.world.apply_to(object, frame)

        camera = self.camera
        camera._update_frame(frame)

    @abstractmethod
    def _create_object(self):
        raise NotImplementedError(
            "This method is only available in the subclass"
        )

    @abstractmethod
    def draw(self):
        """Draw the object"""
        # TODO: What should be done here?
        raise NotImplementedError(
            "This method is only available in the subclass"
        )

    @property
    def object(self):
        return bpy.data.objects[self.name]

    def set_style(self, style):
        self.style.set_style(style)

    def set_material(self, material):
        self.material.set_material(material)

    def set_color(self, color):
        self.color.set_color(color)

    def render(self):
        self.camera.render()

    @property
    def color(self):
        return self._color

    @color.setter
    def color(self, value):
        raise AttributeError("Use `set_color` instead")

    @property
    def style(self):
        return self._style

    @style.setter
    def style(self, value):
        raise AttributeError("Use `set_style` instead")

    @property
    def material(self):
        return self._material

    @material.setter
    def material(self, value):
        raise AttributeError("Use `set_material` instead")