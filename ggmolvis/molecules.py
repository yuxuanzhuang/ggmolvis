import bpy
import molecularnodes as mn
from molecularnodes.entities.trajectory import Trajectory
from molecularnodes.entities.trajectory.selections import Selection
import MDAnalysis as mda
import numpy as np
from typing import Union
from pydantic import BaseModel, Field, validator, ValidationError

from .base import GGMolvisArtist
from .sceneobjects import SceneObject


class Molecule(SceneObject):
    """Class for a molecule."""
    def __init__(
             self,
             atomgroup: Union[mda.AtomGroup, mda.Universe],
             style: str = 'spheres',
             subframes: int = 0,
             name: str = 'atoms',
             location=None,
             rotation=None,
             scale=None
             ):
        """Show the molecule."""

        self.universe = atomgroup if isinstance(atomgroup, mda.Universe) else atomgroup.universe
        self.atomgroup = atomgroup if isinstance(atomgroup, mda.AtomGroup) else self.universe.atoms

        self.style = style
        self.subframes = subframes


        super().__init__(name=name,
                         location=location,
                         rotation=rotation,
                         scale=scale)


    def create_object(self):
        traj = Trajectory(self.universe)
        self.trajectory = traj
        self.trajectory.create_object(name=self.name, style=self.style, subframes=self.subframes)
        self.trajectory.add_selection_from_atomgroup(self.atomgroup, name=self.name)
        # only render the selection
        nodes_mn = self.object.modifiers["MolecularNodes"].node_group
        nodes = nodes_mn.nodes
        links = nodes_mn.links
        name_atoms = self.name
        named_attr_node = nodes.new(type='GeometryNodeInputNamedAttribute')
        named_attr_node.name = name_atoms
        named_attr_node.inputs["Name"].default_value = name_atoms

        # get current style node
        style_lists = [s for s in nodes.keys()  if s.startswith("Style")]
        style_node = nodes[style_lists[0]]
        links.new(named_attr_node.outputs["Attribute"],
                        style_node.inputs["Selection"])
    
    def draw(self):
        pass

    @property
    def object(self):
        return self.trajectory.object