"""
========
GGMolVis
========

GGmolVis is a high-level API for creating molecular visualizations in Blender.
It is built on top of the `MolecularNodes` Blender add-on.
It provides a simple and intuitive interface for creating complex molecular visualizations
with just a few lines of code.

Classes
-------
.. autoclass:: GGMolVis
    :members:
"""
import bpy
from abc import ABC, abstractmethod
import molecularnodes as mn
from molecularnodes.entities.trajectory import Trajectory
from molecularnodes.entities.trajectory.selections import Selection
import MDAnalysis as mda
from MDAnalysis.core.groups import Atom, AtomGroup
import numpy as np
from typing import Union
from pydantic import BaseModel, Field, validator, ValidationError

from . import SESSION
from .base import GGMolvisArtist

from .world import World
from .camera import Camera
from .light import Light
from .properties import Color, Material
from .sceneobjects import SceneObject, Text, Molecule, Shape, Line
from .utils import validate_properties

from loguru import logger


class GGMolVis(GGMolvisArtist):
    """Top level class that contains all the elements of the visualization.
    It is similar to a `Figure` in matplotlib. It contains all the
    `Molecule`, `Shape`, `Text`, `Camera`, `Light`, and `World` objects.
    It also contains the global settings for the visualization like
    `subframes`. It is a singleton class, so only one instance will be
    created in a session.

    During initialization, it creates a global camera and a global world for
    object transformation. The global camera is set to a default position
    and rotation. The global world transformation is set to no positional,
    rotational, or scaling transformation.

    The artists are stored in a dictionary with keys as the type of the
    artist and values as the list of artists of that type.

    Properties:
    -----------
    molecules: list
        List of all `Molecule` objects in the visualization
    shapes: list
        List of all `Shape` objects in the visualization
    texts: list
        List of all `Text` objects in the visualization
    cameras: list
        List of all `Camera` objects in the visualization
    lights: list
        List of all `Light` objects in the visualization
    worlds: list
        List of all `World` transformation objects in the visualization
    global_world: World
        The global world transformation object
    global_camera: Camera
        The global camera object
    subframes: int
        Number of subframes to render. It will be a global setting
        for all objects. Default is 1
    

    
    """
    def __new__(cls):
        if hasattr(SESSION, 'ggmolvis'):
            # If SESSION already has an instance, return that instance
            return SESSION.ggmolvis
        logger.debug("Creating new GGMolVis")
        # Otherwise, create a new instance
        instance = super().__new__(cls)
        SESSION.ggmolvis = instance  # Store the instance in SESSION
        instance._initialized = False
        return instance
    

    def __init__(self):
        if self._initialized:
            return
        self._initialized = True

        super().__init__()
        self.session.ggmolvis = self
        self._artists_dict = {
            'molecules': [],
            'shapes': [],
            'texts': [],
            'cameras': [Camera(name='global_camera')],
            'lights': [],
            'worlds': [World()]
        }
        self._global_world = self.worlds[0]
        self._global_camera = self.cameras[0]

        self._subframes = 0

        # pre-defined camera position
        bpy.data.collections.get('MolecularNodes').objects.link(self._global_camera.object)
        self._global_camera.world.location._set_coordinates((0, -4, 1.3))
        self._global_camera.world.rotation._set_coordinates((83, 0, 0))

        # set up the scene
        self._set_scene()

        self._update_frame(bpy.context.scene.frame_current)

    def _update_frame(self, frame_number):
        """Update the camera's state for the given frame"""
        for artist in self._artists:
            artist._update_frame(frame_number)

        self._global_camera.world._apply_to(self._global_camera.object, frame_number)

    @property
    def _artists(self):
        return [item for sublist in self._artists_dict.values() for item in sublist]

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
    def cameras(self):
        return self._artists_dict['cameras']
    
    @property
    def lights(self):
        return self._artists_dict['lights']
    
    @property
    def worlds(self):
        return self._artists_dict['worlds']

    @property
    def global_world(self):
        return self._global_world
    
    @property
    def global_camera(self):
        return self._global_camera
    
    @property
    def subframes(self):
        return self._subframes
    
    @subframes.setter
    def subframes(self, value):
        self._subframes = value
        for molecule in self.molecules:
            molecule.trajectory.subframes = value

    def _set_scene(self):
        """Set up the scene with transparent background and CYCLES rendering."""
        bpy.context.scene.render.engine = 'CYCLES'
        bpy.context.scene.render.film_transparent = True
        try:
            bpy.context.scene.cycles.device = "GPU"
        except:
            pass
    
    @validate_properties
    def molecule(self,
                 universe: Union[AtomGroup, mda.Universe],
                 style: str = 'spheres',
                 name: str = 'atoms',
                 location: Union[np.ndarray, list] = None,
                 rotation: Union[np.ndarray, list] = None,
                 scale: Union[np.ndarray, list] = None,
                 color='default',
                 material='default'):
        """Create a `Molecule` object and add it to the visualization.
        
        Parameters:
        -----------
        universe: MDAnalysis.AtomGroup or MDAnalysis.Universe
            The AtomGroup or Universe object containing the atoms
        style: str
            The style of the molecule. Default is 'spheres'
        name: str
            The name of the molecule. Default is 'atoms'
        location: np.ndarray or list
            The location of the molecule. Default is None
        rotation: np.ndarray or list
            The rotation of the molecule. Default is None
        scale: np.ndarray or list
            The scale of the molecule. Default is None
        color: str
            The color of the molecule. Default is 'default'
        material: str
            The material of the molecule. Default is 'default'

        Returns:
        --------
        molecule: Molecule
            The created `Molecule` object
        """
        molecule = Molecule(atomgroup=universe.atoms,
                            style=style,
                            name=name,
                            color=color,
                            location=location,
                            rotation=rotation,
                            scale=scale,
                            material=material)
        self.molecules.append(molecule)
        return molecule

    @validate_properties
    def distance(self,
                 atom1: Union[AtomGroup, Atom],
                 atom2: Union[AtomGroup, Atom],
                 name: str = 'distance',
                 location: Union[np.ndarray, list] = None,
                 rotation: Union[np.ndarray, list] = None,
                 scale: Union[np.ndarray, list] = None,
                 mol_color: str ='default',
                 mol_material: str ='ambient',
                 mol_style: str ='sphere',
                 line_color: str ='black',
                 line_material: str ='backdrop',
                 line_style: str ='default',
                 ):
        """Create a `Distance` object and add it to the visualization.
        """
        if atom1.universe != atom2.universe:
            raise ValueError("The atoms belong to different universes")
        mol_atoms = Molecule(atomgroup=AtomGroup(atom1 + atom2),
                            style=mol_style,
                            name=f'{name}_atoms',
                            color=mol_color,
                            location=location,
                            rotation=rotation,
                            scale=scale,
                            material=mol_material)
        
        start_points = np.zeros((atom1.universe.trajectory.n_frames, 3))
        end_points = np.zeros((atom1.universe.trajectory.n_frames, 3))

        for i, ts in enumerate(atom1.universe.trajectory):
            start_points[i] = atom1.center_of_mass()
            end_points[i] = atom2.center_of_mass()
        
        line = Line(start_points=start_points,
                    end_points=end_points,
                    name=f'{name}_distance',
                    location=location,
                    rotation=rotation,
                    scale=scale,
                    color=line_color,
                    material=line_material)
        self.shapes.append(line)
        self.molecules.append(mol_atoms)
        return line
        
        
    @validate_properties
    def line(self,
             start_points: np.ndarray,
             end_points: np.ndarray,
             name: str = 'line',
             location: Union[np.ndarray, list] = None,
             rotation: Union[np.ndarray, list] = None,
             scale: Union[np.ndarray, list] = None,
             color='black',
             material='backdrop'
             ):
        """Create a `Line` object and add it to the visualization.
        
        Parameters:
        -----------
        start_points: np.ndarray
            The start points of the line
        end_points: np.ndarray
            The end points of the line
        name: str
            The name of the line. Default is 'line'
        location: np.ndarray or list
            The location of the line. Default is None
        rotation: np.ndarray or list
            The rotation of the line. Default is None
        scale: np.ndarray or list
            The scale of the line. Default is None
        color: str
            The color of the line. Default is 'black'
        material: str
            The material of the line. Default is 'backdrop'

        Returns:
        --------
        line: Line
            The created `Line` object
        """
        line = Line(start_points=start_points,
                    end_points=end_points,
                    name=name,
                    location=location,
                    rotation=rotation,
                    scale=scale,
                    color=color,
                    material=material)
        self.shapes.append(line)
        return line
    
    @staticmethod
    def check_color(value):
        if value not in bpy.data.materials:
            raise ValueError(f"Material {value} not found in the Blender data")
        return value