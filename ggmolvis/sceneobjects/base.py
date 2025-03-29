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
from ..properties import Color, Material, Style
from ..utils import look_at
from ..delegated_property import DelegatedProperty
from ..properties.dynamic_property import DynamicProperty

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
        super().__init__()

        self.name = name
        # Create the object
        obj = self._create_object()
        # set the name that Blender assigned to the object
        self.name = obj.name

        self._dynamic_properties = {}

        self.location = location
        self.rotation = rotation
        self.scale = scale

        self._init_color(color)
        self._init_material(material)
        self._init_style(style)

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

        self._view_location = camera_center
        self._view_rotation = camera_degree

    def _update_frame(self, frame):
        object = self.object
        self.material._apply_to(object, frame)
        self.color._apply_to(object, frame)
        for prop_name, dynamic_prop in self._dynamic_properties.items():
            """Update dynamic properties for the object"""
            # Update the dynamic properties
            dynamic_prop._update_frame(frame)
            # Apply the updated value to the object
            # This is done in the setter of the property
            setattr(object, prop_name, dynamic_prop.value)

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

    location = DelegatedProperty().delegates(
        getter=lambda self: self.object.location,
        setter=lambda self, value: self._set_transformation('location', value),
        doc="Camera location",
    )
    rotation = DelegatedProperty().delegates(
        getter=lambda self: np.rad2deg(self.object.rotation_euler),
        setter=lambda self, value: self._set_transformation('rotation', value),
        doc="Camera rotation",
    )
    scale = DelegatedProperty().delegates(
        getter=lambda self: self.object.scale,
        setter=lambda self, value: self._set_transformation('scale', value),
        doc="Scale of the object",
        # Allowing for a tuple/list of 3 floats for scale
    )

    def _set_transformation(
            self,
            property_name: str,
            value: Union[List[float], Tuple[float, ...], np.ndarray]):
        """
        Set a transformation property of the object (e.g. location, rotation, scale).

        Parameters
        ----------
        property_name : str
            The name of the property to set (e.g. 'location', 'rotation_euler', 'scale').
        value : Union[List[float], Tuple[float, ...], np.ndarray]
            The value to assign. This can be a 1D sequence of 3 floats (a single coordinate)
            or a 2D array-like with shape (N, 3) (multiple coordinates).

        Note, if property_name is 'scale', a 1-element sequence is expanded to a 3D vector.
        """
        if value is None:
            return
        self._dynamic_properties[property_name] = None  # Clear any existing dynamic property for this name
        value = np.asarray(value)

        # Convert degrees to radians for rotation_euler.
        if property_name == 'rotation':
            value = np.deg2rad(value)
            property_name = 'rotation_euler'  # Blender uses 'rotation_euler' for rotation

        if value.ndim == 1:
            # Handle the 'scale' special case:
            if property_name == 'scale':
                if value.size == 1:
                    value = np.array([value[0]] * 3)
                elif value.size != 3:
                    raise ValueError("Scale must be a 3D vector (either a single value or 3 values).")
            else:
                if value.size != 3:
                    raise ValueError("Coordinates must be a 3D vector (3 elements).")
            setattr(self.object, property_name, value)
        elif value.ndim == 2:
            if value.shape[1] != 3:
                raise ValueError("Each coordinate must be 3D (the array must have 3 columns).")
            self._dynamic_properties[property_name] = DynamicProperty(input_list=value)
        else:
            raise ValueError("Coordinates must be either a 1D vector of size 3 or a 2D array with 3 columns.")

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
