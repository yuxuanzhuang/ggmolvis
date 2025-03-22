import bpy
from typing import Union

from ..utils.node import set_mn_material
from .base import Property
from ..utils import materials_mapping, AVAILABLE_MATERIALS

class Material(Property):
    """Class for the material."""
    _material_modifier = {}

    def _set_property(self):
        """Sets the material property and ensures it's valid."""
        material_name = self.property_name or "backdrop"

        if material_name not in AVAILABLE_MATERIALS:
            raise ValueError(
                f"Material {material_name} is not available. "
                f"Available materials are {AVAILABLE_MATERIALS}"
            )

        # Create a material copy
        default_material = bpy.data.materials.get(materials_mapping.get(material_name))
        if default_material is None:
            raise ValueError(f"Material mapping for {material_name} not found.")

        material_copied = default_material.copy()
        material_copied.name = f"{self.name}"
        self.material_name = material_copied.name

    @property
    def material(self):
        """Return the material, handling missing cases gracefully."""
        # Attempt to retrieve from Blender if not stored in memory
        if hasattr(self, "material_name") and self.material_name in bpy.data.materials:
            return bpy.data.materials[self.material_name]

        return None  # Material is not set

    @material.setter
    def material(self, material: bpy.types.Material):
        """Set the material and update material_name."""
        if not isinstance(material, bpy.types.Material):
            raise TypeError("Assigned material must be a Blender Material.")
        self.material_name = material.name


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


class TrajectoryMaterial(Material):
    def _apply_to(self, obj, frame: int = 0):
        """Apply material to the object."""

        set_mn_material(obj, self.material_name)
        self._modify_material(self.material, frame)