import matplotlib as mpl
import numpy as np

from ..utils.node import set_mn_color
from .base import Property


class Color(Property):
    """Class for the color."""

    # TODO: implement
    def _set_property(self):
        color_name = self.property_name

    def set_color(self, color_name):
        self.property_name = color_name
        self._set_property()

    def _apply_to(self, obj, frame: int = 0):
       # TODO: implement a generic method
       pass
    
    def set_map(self, values, cmap="viridis",
                vmin=None, vmax=None):
        # create a color map
        vmin = np.min(values) if vmin is None else vmin
        vmax = np.max(values) if vmax is None else vmax
        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        self._norm = norm
        self.colormap = mpl.cm.get_cmap(cmap)

        self.colors = self.colormap(norm(values))


class TrajectoryColor(Color):
    """Class for the trajectory color."""

    def _apply_to(self, obj, frame: int = 0):
        """Apply color to the MN object."""
        try:
            self.colors
        except AttributeError:
            return
        if frame >= len(self.colors):
            frame = len(self.colors) - 1
        color = self.colors[frame]

        # set the color
        set_mn_color(obj, color)