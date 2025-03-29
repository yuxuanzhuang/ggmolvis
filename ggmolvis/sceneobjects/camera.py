"""
This module contains the Camera class, which is used to represent the camera in the scene.

Classes
=======
.. autoclass:: Camera
    :members:
    
"""
import bpy
import molecularnodes as mn
from molecularnodes.blender import coll

import MDAnalysis as mda
import numpy as np
from typing import Tuple, List, Union
from contextlib import ExitStack, nullcontext

from .. import SESSION
from ..renderer import Renderer, MovieRenderer
from ..compositor import _set_compositor_bg
from ..delegated_property import DelegatedProperty
from ..properties.dynamic_property import DynamicProperty
from .base import SceneObject

from loguru import logger

class Camera(SceneObject):
    """Class for the camera."""
    def __init__(self,
                 name='Camera',
                 location=None,
                 rotation=None,
                 lens=24,
                 clip_start=0.1,
                 clip_end=1000.0):
        super(SceneObject, self).__init__()
        
        self.name = name
        camera_obj =self._create_object()
        self.name = camera_obj.name

        self._dynamic_properties = {}

        self.location = location
        self.rotation = rotation
        
        self.lens = lens
        self.clip_start = clip_start
        self.clip_end = clip_end

        bpy.context.scene.camera = self.object

    def _create_object(self):
        if self.name not in bpy.data.cameras:
            _ = bpy.data.cameras.new(self.name)
            camera_obj = bpy.data.objects.new(self.camera.name, self.camera)
        else:
            camera_obj = bpy.data.objects[self.name]
        return camera_obj


    camera = DelegatedProperty().delegates(
        getter=lambda self: bpy.data.cameras[self.name],
        setter=None,
        doc="Camera object",
    )
    lens = DelegatedProperty().delegates(
        getter=lambda self: self.camera.lens,
        setter=lambda self, value: setattr(self.camera, 'lens', value),
        doc="Camera lens",
    )
    clip_start = DelegatedProperty().delegates(
        getter=lambda self: self.camera.clip_start,
        setter=lambda self, value: setattr(self.camera, 'clip_start', value),
        doc="Camera clip start",
    )
    clip_end = DelegatedProperty().delegates(
        getter=lambda self: self.camera.clip_end,
        setter=lambda self, value: setattr(self.camera, 'clip_end', value),
        doc="Camera clip end",
    )

    def _update_frame(self, frame_number):
        """Update the camera's state for the given frame"""
        for property_name, dynamic_property in self._dynamic_properties.items():
            if isinstance(dynamic_property, DynamicProperty):
                dynamic_property._update_frame(frame_number)
                setattr(self.object, property_name, dynamic_property.value)

    def render(self,
               mode='image',
               lens=None,
               clip_start=None,
               clip_end=None,
               frame=None,
               filepath=None,
               resolution=None,
               color_mode='RGBA',
               composite_bg_rgba=None):
        """Render the scene with this camera"""
        with ExitStack() as stack:
            stack.enter_context(
                type(self).lens.temporary_set_property(self, lens) if lens is not None else nullcontext()
            )
            stack.enter_context(
                type(self).clip_start.temporary_set_property(self, clip_start) if clip_start is not None else nullcontext()
            )
            stack.enter_context(
                type(self).clip_end.temporary_set_property(self, clip_end) if clip_end is not None else nullcontext()
            )
            stack.enter_context(
                type(self.ggmolvis).ggmolvis.frame.temporary_set_property(self.ggmolvis, frame)
                if frame is not None else nullcontext()
            )
            stack.enter_context(
                type(self.ggmolvis).resolution.temporary_set_property(self.ggmolvis, resolution)
                if resolution is not None else nullcontext()
            )
            stack.enter_context(
                type(self.ggmolvis).color_mode.temporary_set_property(self.ggmolvis, color_mode)
                if color_mode is not None else nullcontext()
            )
           
            bpy.context.scene.camera = self.object

            if composite_bg_rgba is not None:
                _set_compositor_bg(composite_bg_rgba)

            if mode == 'image':
                stack.enter_context(
                    type(self.ggmolvis).render_file_format.temporary_set_property(self.ggmolvis, 'PNG')
                )
                renderer = Renderer(filepath=filepath)
            elif mode == 'movie':
                ffmpeg_codec = self.ggmolvis.ffmpeg_codec
                ffmpeg_constant_rate_factor = self.ggmolvis.ffmpeg_constant_rate_factor
                ffmpeg_format = self.ggmolvis.ffmpeg_format
                stack.enter_context(
                    type(self.ggmolvis).render_file_format.temporary_set_property(self.ggmolvis, 'FFMPEG')
                )
                # when using FFMPEG, these settings are automatically reconfigured
                self.ggmolvis.ffmpeg_codec = ffmpeg_codec
                self.ggmolvis.ffmpeg_constant_rate_factor = ffmpeg_constant_rate_factor
                self.ggmolvis.ffmpeg_format = ffmpeg_format
                
                renderer = MovieRenderer(filepath=filepath)
            renderer.render()
            renderer.display_in_notebook()