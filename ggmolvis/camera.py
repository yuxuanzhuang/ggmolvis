import bpy
import molecularnodes as mn
import MDAnalysis as mda
import numpy as np
from pydantic import BaseModel, Field, validator, ValidationError
from typing import Tuple, List, Union

from .base import GGMolvisArtist
from .world import World, Location, Rotation
from . import SESSION
from .renderer import Renderer

class Camera(GGMolvisArtist):
    """Class for the camera."""
    def __init__(self,
                 name=None,
                 location=None,
                 rotation=None,
                 lens=35.0,
                 clip_start=0.1,
                 clip_end=1000.0):
        super().__init__()
        
        # Initialize the location, rotation, and camera-specific properties
        self.world = World(location=location, rotation=rotation)
        self.lens = lens
        self.clip_start = clip_start
        self.clip_end = clip_end
        
        self.camera = bpy.data.cameras.new(name if name else "Camera")
        self.camera_obj = bpy.data.objects.new(self.camera.name, self.camera)
        bpy.context.scene.collection.objects.link(self.camera_obj)
        self.camera.lens = self.lens
        self.camera.clip_start = self.clip_start
        self.camera.clip_end = self.clip_end

        self.set_view()
        
        
    def update_frame(self, frame_number):
        """Update the camera's state for the given frame"""
        self.world.apply_to(self.camera_obj, frame_number)
    
    
    def set_view(self):
        """Set the current view to this camera"""
        bpy.context.scene.camera = self.camera_obj

    @property
    def object(self):
        return self.camera_obj
    
    def render(self,
               frame=None,
               filepath=None,
               resolution=(1920, 1080)):
        """Render the scene with this camera"""
        if frame is not None:
            bpy.context.scene.frame_set(frame)

        renderer = Renderer(resolution=resolution,
                        filepath=filepath)
        
        renderer.render()
        
        renderer.display_in_notebook()
