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
        material="backdrop",
        style="default",
    ):
        self.world = World(location=location, rotation=rotation, scale=scale)
        super().__init__()
        self.name = name

        self._subframes = 0

        # Create the object
        obj = self._create_object()

        # set the name that Blender assigned to the object
        self.name = obj.name

        self._init_color(color)
        self._init_material(material)
        self._init_style(style)

        self.world._apply_to(self.object)
        self._init_camera()
        self._move_to_collection()

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

        # shift by the location of the object
        center_xyz += np.array(self.object.location)

        camera_center = center_xyz.copy()
        camera_center[1] = camera_center[1] - size_obj_xyz[1] * 3
        camera_center[2] = camera_center[2] + size_obj_xyz[2] * 1.3

        rot = look_at(
            camera_position=camera_center, target_position=center_xyz
        )

        self.camera.world.location._set_coordinates(camera_center)
        self.camera.world.rotation._set_coordinates(np.rad2deg(list(rot.to_euler())))

    def _move_to_collection(self):
        """Move the object to the collection with the same name"""
        mn_coll = bpy.data.collections.get('MolecularNodes')
        coll = mn_coll.children.get(self.name)
        if coll is None:
            coll = bpy.data.collections.new(self.name)
            mn_coll.children.link(coll)
        coll.objects.link(self.object)
        mn_coll.objects.unlink(self.object)
        self.camera._move_to_collection(self.name)


    def _update_frame(self, frame):
        object = self.object
        self.material._apply_to(object, frame)
        self.color._apply_to(object, frame)
        self.world._apply_to(object, frame)

        camera = self.camera
        camera._update_frame(frame)

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

    @property
    def object(self):
        return bpy.data.objects[self.name]
    
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
        self.camera.render(**kwargs)

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
    pass
    