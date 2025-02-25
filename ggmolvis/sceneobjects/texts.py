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
        bpy.ops.object.text_add(location=self.location)
        text_object = bpy.context.active_object
        text_object.data.size = self.text_size
        text_object.rotation_euler = self.rotation
        return text_object

    def _create_object(self):
        self.text_object = self._text_creator()
        return self.object

    @property
    def object(self):
        return self.text_object

    def _update_frame(self, frame):
        self.text_object.data.body = f"{self.frame_counter}{frame}"
