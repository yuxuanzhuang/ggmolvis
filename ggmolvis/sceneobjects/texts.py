import math
from typing import Optional


import bpy
from .base import SceneObject


class Text(SceneObject):
    """Class for the text."""

    def __init__(self,
                 frame_counter: Optional[str] = None,
                 text_size: float = 0.1,
                 name: str = "text",
                 location=None,
                 rotation=None,
                 scale=None,
                 color="black",
                 material="default",
                 style: str = "default"):
        """
        Parameters
        ----------
        frame_counter: str
            When not `None`, prefix frame counter text at the bottom of the
            camera view with the provided string (i.e., "frame: ")
        """
        # TODO: text color is more annoying to control--need a material?
        self.frame_counter = frame_counter
        self.location = location
        self.rotation = rotation
        self.text_size = text_size

        super().__init__(
            name=name,
            location=location,
            rotation=rotation,
            scale=scale,
            color=color,
            material=material,
            style=style,
        )

    def _text_creator(self):
        bpy.ops.object.text_add()
        text_object = bpy.context.active_object
        text_object.data.size = self.text_size
        #text_object.data.extrude = 0.2
        # TODO: avoid hard coding camera of interest for the render;
        # probably by unifying to a single scene camera
        camera = bpy.data.objects["atoms_camera"]
        try:
            con = text_object.constraints["Copy Transforms"]
        except KeyError:
            con = text_object.constraints.new("COPY_TRANSFORMS")
        con.target = camera
        con.mix_mode = 'BEFORE_FULL'
        return text_object

    def _create_object(self):
        self.text_object = self._text_creator()
        return self.object

    @property
    def object(self):
        return self.text_object

    def _update_frame(self, frame):
        # TODO: with multiple cameras and this Text object
        # not knowing which camera will be used for rendering,
        # we want to avoid manually specifying the identity of that
        # camera (probably by unifying to a single camera for the entire
        # scene?)
        camera = bpy.data.objects["atoms_camera"]
        cam_data = camera.data
        focal_length = cam_data.lens
        sensor_width = cam_data.sensor_width
        fov = 2 * math.tan(math.radians(sensor_width / (2 * focal_length)))
        adj_factor = 70
        new_scale = fov * adj_factor
        self.text_object.scale = (new_scale, new_scale, new_scale)
        # at larger focal length we need a smaller y drop
        y_drop = -1.8 + (focal_length / 100)
        self.text_object.location = (-0.5, y_drop, -5)
        self.text_object.data.body = f"{self.frame_counter}{frame}"
