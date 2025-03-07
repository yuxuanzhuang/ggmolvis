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
from .renderer import Renderer, MovieRenderer, frame_number_render_handler
from .compositor import _set_compositor_bg

class Camera(GGMolvisArtist):
    """Class for the camera."""
    def __init__(self,
                 name=None,
                 location=None,
                 rotation=None,
                 lens=24.0,
                 clip_start=0.1,
                 clip_end=1000.0):
        super().__init__()
        
        # Initialize the location, rotation, and camera-specific properties
        self.world = World(location=location, rotation=rotation)
        self.lens = lens
        self.clip_start = clip_start
        self.clip_end = clip_end
        
        self.camera = bpy.data.cameras.new(name if name else "Camera")
        camera_obj = bpy.data.objects.new(self.camera.name, self.camera)
        self.name = camera_obj.name
        self.camera.lens = self.lens
        self.camera.clip_start = self.clip_start
        self.camera.clip_end = self.clip_end

        self.set_view()
        
    @property
    def object(self):
        return bpy.data.objects[self.name]
    
    def _update_frame(self, frame_number):
        """Update the camera's state for the given frame"""
        self.world._apply_to(self.object, frame_number)
        self.frame_number = frame_number
    
    def set_view(self):
        """Set the current view to this camera"""
        bpy.context.scene.camera = self.object

    def set_position(self, location, rotation):
        """Set the position of the camera"""
        self.world.location = location
        self.world.rotation = rotation
    
    def _move_to_collection(self, name):
        """Move the object to the collection with the same name"""
        coll_obj = coll.mn().children.get(name)
        coll_obj.objects.link(self.object)

    def __getstate__(self):
        state = self.__dict__.copy()
        del state['camera']
        return state
    
    def __setstate__(self, state):
        self.__dict__.update(state)
        self.camera = bpy.data.cameras[state['name']]
        self._set_property()

    def render(self,
               mode='image',
               frame=None,
               filepath=None,
               resolution=(640, 360),
               composite_bg_rgba=None,
               add_frame_numbers=None):
        """Render the scene with this camera"""
        bpy.context.scene.camera = self.object
        if frame is not None:
            bpy.context.scene.frame_set(frame)

        if composite_bg_rgba is not None:
           _set_compositor_bg(composite_bg_rgba)

        if add_frame_numbers is not None:
            bpy.app.handlers.render_pre.clear()
            bpy.app.handlers.render_pre.append(frame_number_render_handler)

        if mode == 'image':        
            renderer = Renderer(resolution=resolution,
                                filepath=filepath)
        elif mode == 'movie':
            renderer = MovieRenderer(resolution=resolution,
                            filepath=filepath)
        renderer.render()
        renderer.display_in_notebook()
