import bpy
from typing import Union

from ..utils.node import set_mn_material
from .base import Property
from ..utils import materials_mapping, AVAILABLE_MATERIALS

class Material(Property):
    """Class for the material."""

    def _set_property(self):
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
        self._set_property(material_name)

    def _apply_to(self, obj, frame: int = 0):
        obj.data.materials.append(self.material)

class MoleculeMaterial(Material):
    def _apply_to(self, obj, frame: int = 0):
        """Apply material to the object."""

        set_mn_material(obj, self.material_name)
