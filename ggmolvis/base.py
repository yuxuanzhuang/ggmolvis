import bpy
from abc import ABC, abstractmethod
import molecularnodes as mn
from molecularnodes.entities.trajectory import Trajectory
from molecularnodes.entities.trajectory.selections import Selection
import MDAnalysis as mda
import numpy as np
from typing import Union
from pydantic import BaseModel, Field, validator, ValidationError

from . import SESSION

class GGMolvisArtist:
    """Abstract class for all visualizations. It contains the MNSession in which
       every object is linked to and frame mapping functionality."""
    def __init__(self,
                 session: mn.session.MNSession = SESSION,
                 world_scale: float = 0.01,
                 name: str = None,
                 visible: bool = True,
                 z_order: int = 0):

        self.session = session
        self.world_scale = world_scale
        if not hasattr(self.session, '_ggmolvis'):
            self.session._ggmolvis = set()

        self.session._ggmolvis.add(self)
        self._name = name if name else self.__class__.__name__
        self._visible = visible
        self._z_order = z_order

        self.draw()
        bpy.context.view_layer.update()


    @property
    def visible(self):
        return self._visible
    
    @property
    def name(self):
        return self._name
    

    @name.setter
    def name(self, value):
        self._name = value
    

    @property
    def z_order(self):
        return self._z_order
    

    @z_order.setter
    def z_order(self, value):
        self._z_order = value
    

    @abstractmethod
    def draw(self):
        """Abstract method to draw the object within the session"""
        pass
    

    def set_visible(self, visibility):
        self._visible = visibility
    

    @abstractmethod
    def update_frame(self, frame):
        """Abstract method to update the object's state for the given frame"""
        pass
    

    def remove(self):
        """Remove the object from the session"""
        self.session.remove(self)


class DynamicGGMolvisArtist(GGMolvisArtist):
    """Abstract class for all dynamic visualizations."""
    pass