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
    visible: bool
        Visibility of the object. Default is True
    z_order: int
        Z-order of the object. The object with higher z-order
        will be rendered on top of the object with lower z-order. Default is 0
    session: MNSession
        The linked `MNSession` object in molecularnodes package
    ggmolvis: set
        The set of all `GGMolvisArtist` objects in the session
    subframes: int
        Number of subframes to render. It will be a global setting
        for all objects. Default is 0. For clarity, when subframes is set to `1`
        the total frame count will double, and when it is set to `2` the
        total frame count will triple.
    average: int
        Number of flanking frames to average over--this can help reduce
        "jittering" in movies. In contrast to `subframes`, no new frames
        are added. It will be a global setting for all objects. Default is 0.
    """
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
        """Set the object to be visible"""
        self._visible = True

    def set_invisible(self):
        """Set the object to be invisible"""
        self._visible = False

    @abstractmethod
    def _update_frame(self, frame):
        """Abstract method to update the object's state for the given frame"""
        raise NotImplementedError("The method _update_frame must be implemented in the subclass")

    @property
    def ggmolvis(self):
        return self.session.ggmolvis

    @property
    def subframes(self):
        return self.ggmolvis.subframes

    @subframes.setter
    def subframes(self, value):
        self.ggmolvis.subframes = value

    @property
    def average(self):
        return self.ggmolvis.average

    @average.setter
    def average(self, value):
        self.ggmolvis.average = value

    def _remove(self):
        """Remove the object from the session"""
        self.session._ggmolvis.remove(self)
