"""
This module contains the base classes for the scene objects.

Classes
=======
.. autoclass:: SceneObject
    :members:
.. autoclass:: SceneObjectCollection
    :members:
"""
import bpy
from abc import abstractmethod

import numpy as np
from typing import Tuple, List, Union

from ..base import GGMolvisArtist
from ..world import World
from ..camera import Camera
from ..properties import Color, Material, Style
from ..utils import look_at
from ..delegated_property import DelegatedProperty

from loguru import logger

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
        material="backdrop",
        style="default",
    ):
        self.world = World(location=location, rotation=rotation, scale=scale)
        super().__init__()
        self.name = name

        # Create the object
        obj = self._create_object()

        # set the name that Blender assigned to the object
        self.name = obj.name

        self._init_color(color)
        self._init_material(material)
        self._init_style(style)

        self.world._apply_to(self.object)

        self.camera_world = World()
        self._set_camera_view()
        self._camera_view_active = False

        self.draw()
        self._update_frame(bpy.context.scene.frame_current)

    def _init_color(self, color="black"):
        self._color = Color(self, color)

    def _init_material(self, material="backdrop"):
        self._material = Material(self, material)

    def _init_style(self, style="default"):
        self._style = Style(self, style)

    def _set_camera_view(self):
        """
        Set camera view based on the object.
        # TODO: this is way to janky, needs to be improved
        """
        # only run if the camera view is active
        # to avoid unnecessary calculations
        size_obj_xyz = np.array(self.object.dimensions)
        # center of the object
        center_xyz = np.zeros(3)
        bbox = self.object.bound_box
        for v in bbox:
            center_xyz += np.array(v)
        center_xyz /= 8

        # shift by the location of the object
        center_xyz += np.array(self.object.location)

        camera_center = center_xyz.copy()
        camera_center[1] = camera_center[1] - size_obj_xyz[1] * 3
        camera_center[2] = camera_center[2] + size_obj_xyz[2] * 1.3

        rot = look_at(
            camera_position=camera_center, target_position=center_xyz
        )

        camera_degree = np.rad2deg(list(rot.to_euler()))

        self.camera_world.location.coordinates = camera_center
        self.camera_world.rotation.coordinates = camera_degree

    def _update_frame(self, frame):
        object = self.object
        self.material._apply_to(object, frame)
        self.color._apply_to(object, frame)
        self.world._apply_to(object, frame)
        if self._camera_view_active:
            self._set_camera_view()

    @abstractmethod
    def _create_object(self):
        """Create the object

        Returns
        -------
        bpy.types.Object
            The created object
        """

        raise NotImplementedError(
            "This method is only available in the subclass"
        )

    def draw(self):
        """Draw the object"""
        self.style._apply_to(self.object)
        self.color._apply_to(self.object)
        self.material._apply_to(self.object)
    
    object = DelegatedProperty().delegates(
        getter=lambda self: bpy.data.objects[self.name],
        # no setter here, as we don't want to set the object directly
        setter=None,
        doc="The Blender object in the scene.",
    )

    visible = DelegatedProperty().delegates(
        getter=lambda self: not self.object.hide_render,
        setter=lambda self, value: setattr(self.object, "hide_render", not value),
        doc="Visibility of the object during rendering.",
        allowed_type=bool,
    )

    visible_in_viewport = DelegatedProperty().delegates(
        getter=lambda self: not self.object.hide_viewport,
        setter=lambda self, value: setattr(self.object, "hide_viewport", not value),
        doc="Visibility of the object in the viewport.",
        allowed_type=bool,
    )

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        self._name = value

    def set_style(self, style):
        self.style.set_style(style)
        self.style._apply_to(self.object)

    def set_material(self, material):
        self.material.set_material(material)
        self.material._apply_to(self.object)

    def set_color(self, color):
        self.color.set_color(color)
        self.color._apply_to(self.object)

    def render(self, **kwargs):
        """
        Render the object with the corresponding
        camera settings.
        """
        self.ggmolvis.render(
            object=self,
            **kwargs)

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


class SceneObjectCollection:
    """Class for the collection of scene objects"""
    #TODO: Implement the class
