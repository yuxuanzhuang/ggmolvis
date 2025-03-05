from typing import Optional


import bpy
from mathutils import Vector
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
        scene = bpy.context.scene
        camera = scene.camera
        cam_location = camera.location
        text_object.data.extrude = 0.2
        text_object.location = cam_location + camera.matrix_world @ Vector((0, 0, 10))
        direction_to_camera = (cam_location - text_object.location).normalized()
        rot_quat = direction_to_camera.to_track_quat('-Z', 'Y')
        text_object.rotation_mode = 'QUATERNION'
        text_object.rotation_quaternion = rot_quat
        #constraint = text_object.constraints.new(type='TRACK_TO')
        #constraint.target = camera
        #constraint.track_axis = 'TRACK_Z'
        #constraint.up_axis = 'UP_Z'
        #bpy.context.view_layer.update()
        return text_object

    def _create_object(self):
        self.text_object = self._text_creator()
        return self.object

    @property
    def object(self):
        return self.text_object

    def _update_frame(self, frame):
        self.text_object.data.body = f"{self.frame_counter}{frame}"
