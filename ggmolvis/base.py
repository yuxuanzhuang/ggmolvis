"""
Base classes for all visualizations in GGMolVis.

Classes
=======

..autoclass:: GGMolvisArtist
    :members:
"""
import bpy
from abc import ABC, abstractmethod
import molecularnodes as mn
import MDAnalysis as mda
import numpy as np
from typing import Union

from . import SESSION
from .delegated_property import DelegatedProperty

class GGMolvisArtist(ABC):
    """Abstract class for all visualizations in GGMolVis.
    It creates a link to the `MNSession` in Molecular Nodes.    
    Every `GGMolvisArtist` object will be stored inside the
    `MNSession._ggmolvis`.

    Properties:
    -----------
    name: str
        Name of the object. Default is the class name.
    world_scale: float
        The scale of the object in the world. It is different
        from the scale in Blender as it is directly applied to
        the coordinates of the object. Default is 0.01
    session: MNSession
        The linked `MNSession` object in molecularnodes package
    ggmolvis: set
        The set of all `GGMolvisArtist` objects in the session
    """
    def __init__(self):
        # only one MNSession will be used to keep everything in sync.
        self._session = SESSION

        # default world_scale value used in Molecular Node
        self._world_scale = 0.01
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
    def session(self):
        return self._session    
    
    @abstractmethod
    def _update_frame(self, frame):
        """Abstract method to update the object's state for the given frame"""
        raise NotImplementedError("The method _update_frame must be implemented in the subclass")

    @property
    def ggmolvis(self):
        return self.session.ggmolvis

    def _remove(self):
        """Remove the object from the session"""
        self.session._ggmolvis.remove(self)
    
    def __getstate__(self):
        state = self.__dict__.copy()
        del state['_session']

        # remove bpy objects related to the object
        if 'object' in state:
            del state['object']

        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        # reattach the session
        SESSION = mn.session.get_session()
        self._session = SESSION
