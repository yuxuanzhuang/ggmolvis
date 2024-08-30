import bpy
from abc import ABC, abstractmethod
import molecularnodes as mn
from molecularnodes.entities.trajectory import Trajectory
from molecularnodes.entities.trajectory.selections import Selection
import MDAnalysis as mda
import numpy as np
from typing import Union
from pydantic import BaseModel, Field, validator, ValidationError

from . import SESSION
from .base import GGMolvisArtist

from .world import World
from .camera import Camera
from .light import Light
from .properties import Color, Material
from .sceneobjects import SceneObject, Text, Shape, Line
from .molecules import Molecule


class GGMolVis(GGMolvisArtist):
    """Top level class that contains all the elements of the visualization."""
    def __init__(self, session=SESSION):
        super().__init__(session)
        self.artists = []
        self._artists_dict = {
            'molecules': [],
            'shapes': [],
            'texts': [],
            'global_world': World(),
            'cameras': [Camera()],
            'lights': [],
            'worlds': []
        }

    @property
    def molecules(self):
        return self._artists_dict['molecules']
    
    @property
    def shapes(self):
        return self._artists_dict['shapes']
    
    @property
    def texts(self):
        return self._artists_dict['texts']
    
    @property
    def global_world(self):
        return self._artists_dict['global_world']
    
    @property
    def cameras(self):
        return self._artists_dict['cameras']
    
    @property
    def lights(self):
        return self._artists_dict['lights']
    
    @property
    def worlds(self):
        return self._artists_dict['worlds']

    def molecule(self,
                 universe: Union[mda.AtomGroup, mda.Universe],
                 style: str = 'spheres',
                 subframes: int = 0,
                 name: str = 'atoms',
                 world: World = None):
        molecule = Molecule(universe, style, subframes, name)

        if world is not None:
            molecule.world = world
        else:
            molecule.world = self.global_world
        self.molecules.append(molecule)

    def line(self,
             start_points: np.ndarray,
             end_points: np.ndarray,
             location=None,
             rotation=None,
             scale=None,
             world: World = None):
        
        line = Line(start_points, end_points, location, rotation, scale)
        if world is not None:
            line.world = world
        else:
            line.world = self.global_world
        self.shapes.append(line)