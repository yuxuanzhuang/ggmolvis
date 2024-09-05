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
        pass