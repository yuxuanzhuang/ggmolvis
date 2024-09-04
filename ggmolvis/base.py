import bpy
from abc import ABC, abstractmethod
import molecularnodes as mn
import MDAnalysis as mda
import numpy as np
from typing import Union

from . import SESSION

class GGMolvisArtist(ABC):
    """Abstract class for all visualizations. It contains the MNSession in which
       every object is linked to and frame mapping functionality."""
    def __init__(self):
        # only one MNSession will be used to keep everything in sync.
        self._session = SESSION

        # default world_scale value used in Molecular Node
        self._world_scale = 0.01
        self._visible = True
        self._z_order = 0
        self._name = self.__class__.__name__

        #
        if not hasattr(self.session, '_ggmolvis'):
            self.session._ggmolvis = set()

        self.session._ggmolvis.add(self)

        # when creating a new rendering object
        # always update the blender view layer
        bpy.context.view_layer.update()

    @property
    def name(self):
        return self._name
    
    @name.setter
    def name(self, value):
        self._name = value

    @property
    def world_scale(self):
        return self._world_scale
    
    @world_scale.setter
    def world_scale(self, value):
        self._world_scale = value
    
    @property
    def visible(self):
        return self._visible

    @visible.setter
    def visible(self, value):
        self._visible = value

    @property
    def z_order(self):
        return self._z_order

    @z_order.setter
    def z_order(self, value):
        self._z_order = value

    @property
    def session(self):
        return self._session    
    
    def set_visible(self):
        self._visible = True

    def set_invisible(self):
        self._visible = False

    @abstractmethod
    def _update_frame(self, frame):
        """Abstract method to update the object's state for the given frame"""
        pass

    @property
    def ggmolvis(self):
        return self.session.ggmolvis

    @property
    def subframes(self):
        return self.ggmolvis.subframes

    @subframes.setter
    def subframes(self, value):
        self.ggmolvis.subframes = value

    def _remove(self):
        """Remove the object from the session"""
        self.session.remove(self)
