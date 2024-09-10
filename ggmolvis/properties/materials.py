import bpy
from typing import Union

from ..utils.node import set_mn_material
from .base import Property
from ..utils import materials_mapping, AVAILABLE_MATERIALS

class Material(Property):
    """Class for the material."""
#    _material_modifier = {}

    def _set_property(self):
        self._material_modifier = {}
        material_name = self.property_name
        if material_name is None:
            material_name = "backdrop"

        if material_name not in AVAILABLE_MATERIALS:
            raise ValueError(
                f"Material {material_name} is not available."
                f"Available materials are {AVAILABLE_MATERIALS}"
            )

        # create a material copy
        default_material = bpy.data.materials[materials_mapping[material_name]]
        self.material = default_material.copy()
        self.material.name = f"{self.name}"

        self.material_name = self.material.name

    def set_material(self, material_name):
        self.property_name = material_name
        self._set_property()

    def _apply_to(self, obj, frame: int = 0):
        obj.data.materials.append(self.material)
    
    def _modify_material(self, material, frame: int = 0):
        for key, values in self.material_modifier.items():
            # if values is a number
            if isinstance(values, (int, float)):
                value = values
            else:
                if frame >= len(values):
                    frame = len(values) - 1
                value = values[frame]

            material.node_tree.nodes.get('Principled BSDF').inputs[key].default_value = value

    @property
    def material_modifier(self):
        """Return the material modifier.
        It is a dictionary with the key as property name and value as the modifier.
        """
        return self._material_modifier

    @material_modifier.setter
    def material_modifier(self, value: dict):
        if not isinstance(value, dict):
            raise ValueError("Material modifier must be a dictionary "
            "with the key as property name and value as the modifier.")
        self._material_modifier.update(value)


class MoleculeMaterial(Material):
    def _apply_to(self, obj, frame: int = 0):
        """Apply material to the object."""

        set_mn_material(obj, self.material_name)
        self._modify_material(self.material, frame)