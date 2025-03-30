from .base import GGMolvisArtist
from loguru import logger
from typing import List, Tuple, Union
import numpy as np

class DynamicProperty(GGMolvisArtist):
    """
    A class to update a dynamic property on a target object with each frame change.

    This class updates a property on the target by applying a transformation
    function that takes the current frame and the current property value and returns a new value.
    """
    def __init__(self,
                 input_list: list):
        """
        Parameters:
        ----------
        target_property : object
            The property of the target to be updated.
        input_list : list, optional
            A list of values to be used for the transformation.
        """
        super().__init__()

        self.input_list = input_list

    def _update_frame(self, frame):
        try:
            new_value = self.input_list[frame]
        except IndexError:
            new_value = self.input_list[-1]  # Use the last value if frame is out of range

        self._value = new_value

    @property
    def value(self):
        return self._value
    
    @value.setter
    def value(self, new_value):
        self._value = new_value