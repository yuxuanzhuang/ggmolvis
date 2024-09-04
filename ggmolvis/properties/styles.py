from .base import Property
from molecularnodes.blender.nodes import styles_mapping as mol_styles_mapping
from ..utils.node import swap_style

MOL_AVAILABLE_STYLES = mol_styles_mapping.keys()


class Style(Property):
    """Class for the style."""

    AVAILABLE_STYLES = ["default"]

    def _set_property(self):
        style_name = self.property_name
        if style_name is None:
            style_name = "spheres"
        if style_name not in self.AVAILABLE_STYLES:
            raise ValueError(
                f"Style {style_name} is not available."
                f"Available styles are {self.AVAILABLE_STYLES}"
            )
        self.style_name = style_name

    def set_style(self, style_name):
        self.property_name = style_name
        self._set_property()

    def _apply_to(self, obj, frame: int = 0):
        # TODO: implement it
        pass


class MoleculeStyle(Style):
    """Style explicitly for Molecules"""

    AVAILABLE_STYLES = MOL_AVAILABLE_STYLES

    def _apply_to(self, obj, frame: int = 0):
        """Apply style to the object."""

        swap_style(obj, mol_styles_mapping[self.style_name])
