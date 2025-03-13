import bpy
from abc import abstractmethod
from typing import Union

from ..base import GGMolvisArtist


class Property(GGMolvisArtist):
    """Class for the property of the visulizations."""

    def __init__(self, scene_object, property_name, **kwargs):
        super().__init__()
        self._scene_object = scene_object
        self._property_name = property_name
        self._set_property()
        self.color = kwargs.get("color")

    @abstractmethod
    def _set_property(self):
        """Set the property of the visualization."""
        raise NotImplementedError(
            "This method must be implemented in the subclass"
        )

    def _apply_to(self, obj, frame: int = 0):
        raise NotImplementedError(
            "This method must be implemented in the subclass"
        )

    @property
    def scene_object(self):
        return self._scene_object

    @property
    def property_name(self):
        return self._property_name

    @property_name.setter
    def property_name(self, value):
        self._property_name = value

    def _update_frame(self, frame):
        # TODO: Implement the update_frame method
        pass
