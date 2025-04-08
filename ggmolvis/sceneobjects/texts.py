import math
from typing import Optional
import bpy
from .base import SceneObject

class Text(SceneObject):
    """Class for the text."""
    def __init__(self,
                 text: str = "text",
                 text_size: float = 0.1,
                 name: str = "Text",
                 location=None,
                 rotation=None,
                 scale=None,
                 color="black",
                 material="backdrop",
                 style: str = "default"):
        """
        Parameters
        ----------
        text: str
            The static text to display.
        text_size: float
            The size of the text.
        """
        self._text = text
        # TODO: text color is more annoying to control--need a material?
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
        try:
            con = text_object.constraints["Copy Transforms"]
        except KeyError:
            con = text_object.constraints.new("COPY_TRANSFORMS")
        con.target = self.ggmolvis.camera.object
        con.mix_mode = 'BEFORE_FULL'
        return text_object

    def _create_object(self):
        self.text_object = self._text_creator()
        return self.object

    @property
    def object(self):
        return self.text_object

    @property
    def text(self):
        return self._text
    
    @text.setter
    def text(self, value):
        self._text = value

    def _update_frame(self, frame):
        camera = self.ggmolvis.camera.object
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
        self.text_object.data.body = f"{self.text}"


class FrameCounter(Text):
    """Class for the frame counter."""
    def __init__(self,
                 object=None,
                 **kwargs):
        """
        Parameters
        ----------
        object: The object from which the `.frame` attribute will be read.
        """
        if object is None:
            object = self.ggmolvis
        self.frame_counter = object
        super().__init__(**kwargs)

    @property
    def text(self):
        return f"Frame: {self.frame_counter.frame}"