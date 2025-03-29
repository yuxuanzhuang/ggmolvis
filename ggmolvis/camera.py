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

from .base import GGMolvisArtist
from .world import World, Location, Rotation
from . import SESSION
from .renderer import Renderer, MovieRenderer
from .compositor import _set_compositor_bg
from .delegated_property import DelegatedProperty

class Camera(GGMolvisArtist):
    """Class for the camera."""
    def __init__(self,
                 name='Camera',
                 location=None,
                 rotation=None,
                 lens=24.0,
                 clip_start=0.1,
                 clip_end=1000.0):
        super().__init__()
        
        # Initialize the location, rotation, and camera-specific properties
        self.world = World(location=location, rotation=rotation)
        
        # get existing camera or create a new one
        if name not in bpy.data.cameras:
            _ = bpy.data.cameras.new(name if name else "Camera")
            camera_obj = bpy.data.objects.new(self.camera.name, self.camera)
            self.name = camera_obj.name
        else:
            self.name = self.object.name
        self.lens = self.lens
        self.clip_start = self.clip_start
        self.clip_end = self.clip_end

        bpy.context.scene.camera = self.object

    camera = DelegatedProperty().delegates(
        getter=lambda self: bpy.data.cameras[self.name],
        setter=None,
        doc="Camera object",
    )
    lens = DelegatedProperty().delegates(
        getter=lambda self: self.camera.lens,
        setter=lambda self, value: setattr(self.camera, 'lens', value),
        doc="Camera lens",
        allowed_type=float,
    )
    clip_start = DelegatedProperty().delegates(
        getter=lambda self: self.camera.clip_start,
        setter=lambda self, value: setattr(self.camera, 'clip_start', value),
        doc="Camera clip start",
        allowed_type=float,
    )
    clip_end = DelegatedProperty().delegates(
        getter=lambda self: self.camera.clip_end,
        setter=lambda self, value: setattr(self.camera, 'clip_end', value),
        doc="Camera clip end",
        allowed_type=float,
    )
    object = DelegatedProperty().delegates(
        getter=lambda self: bpy.data.objects[self.name],
        setter=None,
        doc="Camera object",
    )
    
    def _update_frame(self, frame_number):
        """Update the camera's state for the given frame"""
        self.world._apply_to(self.object, frame_number)
    
    def set_position(self, location, rotation):
        """Set the position of the camera"""
        self.world.location = location
        self.world.rotation = rotation
    
    def render(self,
               mode='image',
               lens=None,
               clip_start=None,
               clip_end=None,
               frame=None,
               filepath=None,
               resolution=None,
               composite_bg_rgba=None):
        """Render the scene with this camera"""
        bpy.context.scene.camera = self.object
        if lens is not None:
            old_lens = self.lens
            self.lens = lens
        if clip_start is not None:
            old_clip_start = self.clip_start
            self.clip_start = clip_start
        if clip_end is not None:
            old_clip_end = self.clip_end
            self.clip_end = clip_end
        if frame is not None:
            old_frame = self.ggmolvis.frame
            self.ggmolvis.frame_current = frame
        if resolution is not None:
            old_resolution = self.ggmolvis.resolution
            self.ggmolvis.resolution = resolution

        if composite_bg_rgba is not None:
           _set_compositor_bg(composite_bg_rgba)

        if mode == 'image':
            renderer = Renderer(filepath=filepath)
        elif mode == 'movie':
            renderer = MovieRenderer(filepath=filepath)
        renderer.render()
        renderer.display_in_notebook()

        if lens is not None:
            self.lens = old_lens
        if clip_start is not None:
            self.clip_start = old_clip_start
        if clip_end is not None:
            self.clip_end = old_clip_end
        if frame is not None:
            self.ggmolvis.frame = frame
        if resolution is not None:
            self.ggmolvis.resolution = resolution