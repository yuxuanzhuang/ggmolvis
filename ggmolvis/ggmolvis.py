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
import functools
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
from .sceneobjects import SceneObject, Text, Trajectory, Shape, Line
from .utils import validate_properties
from .delegated_property import DelegatedProperty

from loguru import logger


class GGMolVis(GGMolvisArtist):
    """Top level class that contains all the elements of the visualization.
    It is similar to a `Figure` in matplotlib. It contains all the
    `Trajectory`, `Shape`, `Text`, `Light`, and `World` objects.
    It also contains the global settings for the visualization like
    `subframes`, `average`. It is a singleton class, so only one instance will be
    created in a session.

    During initialization, it creates a camera and a global world for
    object transformation. The camera is set to a default position
    and rotation. The global world transformation is set to no positional,
    rotational, or scaling transformation.

    The artists are stored in a dictionary with keys as the type of the
    artist and values as the list of artists of that type.

    Properties:
    -----------
    trajectories: list
        List of all `Trajectory` objects in the visualization
    shapes: list
        List of all `Shape` objects in the visualization
    texts: list
        List of all `Text` objects in the visualization
    lights: list
        List of all `Light` objects in the visualization
    worlds: list
        List of all `World` transformation objects in the visualization
    global_world: World
        The global world transformation object
    camera: Camera
        The camera object
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
            'trajectories': [],
            'shapes': [],
            'texts': [],
            'lights': [],
            'worlds': [World()]
        }
        self._global_world = self.worlds[0]
        self._camera = Camera()

        # pre-defined camera position
        self._camera.world.location._set_coordinates((0, -4, 1.3))
        self._camera.world.rotation._set_coordinates((83, 0, 0))

        # set up the scene
        self._set_scene()

        self._update_frame(self.frame_current)

    def _update_frame(self, frame_number):
        """Update the camera's state for the given frame"""
        for artist in self._artists:
            artist._update_frame(frame_number)

        self._camera.world._apply_to(self._camera.object, frame_number)

    @property
    def _artists(self):
        return [item for sublist in self._artists_dict.values() for item in sublist]

    @property
    def trajectories(self):
        return self._artists_dict['trajectories']
    
    @property
    def shapes(self):
        return self._artists_dict['shapes']
    
    @property
    def texts(self):
        return self._artists_dict['texts']
    
    @property
    def lights(self):
        return self._artists_dict['lights']
    
    @property
    def worlds(self):
        return self._artists_dict['worlds']

    @property
    def global_world(self):
        if not hasattr(self, '_global_world'):
            self._global_world = World()
        return self._global_world
    
    @property
    def camera(self):
        if not hasattr(self, '_camera'):
            self._camera = Camera()
            self._camera.world.location._set_coordinates((0, -4, 1.3))
            self._camera.world.rotation._set_coordinates((83, 0, 0))
        return self._camera

    def _set_scene(self):
        """Set up the scene with transparent background and CYCLES rendering."""
        self.render_engine = 'CYCLES'
        self.transparent_background = True
        try:
            self.cycles_device = "GPU"
        except:
            pass
    
    # Rendering properties
    render_engine = DelegatedProperty().delegates(
        getter=lambda self: bpy.context.scene.render.engine,
        setter=lambda self, value: setattr(bpy.context.scene.render, 'engine', value),
        default='CYCLES',
        doc="Render engine to use for rendering.",
        allowed_type=str
    )
    transparent_background = DelegatedProperty().delegates(
        getter=lambda self: bpy.context.scene.render.film_transparent,
        setter=lambda self, value: setattr(bpy.context.scene.render, 'film_transparent', value),
        default=True,
        doc="Whether to use a transparent background for rendering.",
        allowed_type=bool
    )
    cycles_device = DelegatedProperty().delegates(
        getter=lambda self: bpy.context.scene.cycles.device,
        setter=lambda self, value: setattr(bpy.context.scene.cycles, 'device', value),
        default='GPU',
        doc="Render device to use for rendering.",
        allowed_type=str
    )
    resolution = DelegatedProperty().delegates(
        getter=lambda self: (bpy.context.scene.render.resolution_x,
                             bpy.context.scene.render.resolution_y),
        setter=lambda self, value: (
            setattr(bpy.context.scene.render, 'resolution_x', value[0]),
            setattr(bpy.context.scene.render, 'resolution_y', value[1])
        ),
        default=(1920, 1080),
        doc="Resolution of the scene in pixels (Width, Height).",
        allowed_type=(tuple, list),
    )
    resolution_scale = DelegatedProperty().delegates(
        getter=lambda self: bpy.context.scene.render.resolution_percentage,
        setter=lambda self, value: setattr(bpy.context.scene.render, 'resolution_percentage', value),
        default=100,
        doc="Resolution scale of the scene in percentage.",
        allowed_type=int
    )
    frame_start = DelegatedProperty().delegates(
        getter=lambda self: bpy.context.scene.frame_start,
        setter=lambda self, value: setattr(bpy.context.scene, 'frame_start', value),
        default=0,
        doc="The first frame of the scene. This is the frame that will be rendered for a movie.",
        allowed_type=int
    )
    frame_end = DelegatedProperty().delegates(
        getter=lambda self: bpy.context.scene.frame_end,
        setter=lambda self, value: setattr(bpy.context.scene, 'frame_end', value),
        default=0,
        doc="The last frame of the scene. This is the frame that will be rendered for a movie.",
        allowed_type=int
    )
    frame_step = DelegatedProperty().delegates(
        getter=lambda self: bpy.context.scene.frame_step,
        setter=lambda self, value: setattr(bpy.context.scene, 'frame_step', value),
        default=1,
        doc="The number of frames to skip between renders. "
            "For example, if frame_step is set to 2, every second frame will be rendered.",
        allowed_type=int
    )
    frame_current = DelegatedProperty().delegates(
        getter=lambda self: bpy.context.scene.frame_current,
        setter=lambda self, value: setattr(bpy.context.scene, 'frame_current', value),
        default=0,
        doc="The current frame of the scene. This is the frame that will be rendered for an image.",
        allowed_type=int
    )
    frame_rate = DelegatedProperty().delegates(
        getter=lambda self: bpy.context.scene.render.fps,
        setter=lambda self, value: setattr(bpy.context.scene.render, 'fps', value),
        default=24,
        doc="The frame rate of the scene in frames per second.",
        allowed_type=int
    )
    render_file_format = DelegatedProperty().delegates(
        getter=lambda self: bpy.context.scene.render.image_settings.file_format,
        setter=lambda self, value: setattr(bpy.context.scene.render.image_settings, 'file_format', value),
        default='PNG',
        doc="The file format for the rendered image.",
        allowed_type=str
    )
    ffmpeg_format = DelegatedProperty().delegates(
        getter=lambda self: bpy.context.scene.render.ffmpeg.format,
        setter=lambda self, value: setattr(bpy.context.scene.render.ffmpeg, 'format', value),
        default='MPEG4',
        doc="The file format for the rendered movie.",
        allowed_type=str
    )
    ffmpeg_codec = DelegatedProperty().delegates(
        getter=lambda self: bpy.context.scene.render.ffmpeg.codec,
        setter=lambda self, value: setattr(bpy.context.scene.render.ffmpeg, 'codec', value),
        default='H264',
        doc="The codec for the rendered movie.",
        allowed_type=str
    )
    ffmpeg_constant_rate_factor = DelegatedProperty().delegates(
        getter=lambda self: bpy.context.scene.render.ffmpeg.constant_rate_factor,
        setter=lambda self, value: setattr(bpy.context.scene.render.ffmpeg, 'constant_rate_factor', value),
        default='HIGH',
        doc="The constant rate factor for the rendered movie.",
        allowed_type=str
    )

    # Trajectory properties
    subframes = DelegatedProperty().delegates(
        getter=lambda self: self.trajectories[0].trajectory.subframes if self.trajectories else 0,
        setter=lambda self, value: (
            setattr(self, '_subframes', value),
            [setattr(trajectory.trajectory, 'subframes', value) for trajectory in self.trajectories]
        ),
        default=0,
        doc="Number of subframes to render. It will be a global setting "
            "for all objects. Default is 0. For clarity, when subframes is set to `1` "
            "the total frame count will double, and when it is set to `2` the "
            "total frame count will triple. "
            "The first trajectory in the list will be used to get the value.",
        allowed_type=int
    )
    average = DelegatedProperty().delegates(
        getter=lambda self: self.trajectories[0].trajectory.average if self.trajectories else 0,
        setter=lambda self, value: (
            setattr(self, '_average', value),
            [setattr(trajectory.trajectory, 'average', value) for trajectory in self.trajectories]
        ),
        default=0,
        doc="Number of flanking frames to average over--this can help reduce "
            "jittering in movies. In contrast to `subframes`, no new frames "
            "are added. It will be a global setting for all objects. Default is 0. "
            "The first trajectory in the list will be used to get the value.",
        allowed_type=int
    )

    def render(self,
               object: SceneObject = None,
               track: bool = False,
               frame: int = None,
               frame_range: tuple = None,
               **kwargs):
        """
        Render the current scene.
        """
        if frame is not None and frame_range is not None:
            raise ValueError("Both frame and frame_range cannot be set")
        if frame is not None:
            render_mode = 'image'
            if kwargs.get('mode', None) == 'movie':
                logger.warning("mode is set to 'movie' but frame is set. "
                        "Changing mode to 'image'")
            kwargs['mode'] = 'image'
            self.frame_current = frame
        elif frame_range is not None:
            render_mode = 'movie'
            if kwargs.get('mode', None) == 'image':
                logger.warning("mode is set to 'image' but frame_range is set. "
                        "Changing mode to 'movie'")
            kwargs['mode'] = 'movie'
            if len(frame_range) != 3:
                raise ValueError("frame_range must be a tuple of 3 integers (start, end, step)")
            start, end, step = frame_range
            old_start = self.frame_start
            old_end = self.frame_end
            old_step = self.frame_step
            self.frame_start = start
            self.frame_end = end
            self.frame_step = step
        else:
            render_mode = kwargs.pop('mode', 'image')
        kwargs['mode'] = render_mode

        if object is not None:
            current_world = self.camera.world
            if track:
                object._camera_view_active = True
            object._set_camera_view()
            self.camera.world = object.camera_world
            self.camera.render(**kwargs)
            object._camera_view_active = False
            self.camera.world = current_world
        else:
            self.camera.render(**kwargs)
        
        if frame_range is not None:
            self.frame_start = old_start
            self.frame_end = old_end
            self.frame_step = old_step
    
    @validate_properties
    def trajectory(self,
                 universe: Union[AtomGroup, mda.Universe],
                 style: str = 'spheres',
                 name: str = 'atoms',
                 location: Union[np.ndarray, list] = None,
                 rotation: Union[np.ndarray, list] = None,
                 scale: Union[np.ndarray, list] = None,
                 color='default',
                 material='default',
                 ):
        """Create a `Trajectory` object and add it to the visualization.
        
        Parameters:
        -----------
        universe: MDAnalysis.AtomGroup or MDAnalysis.Universe
            The AtomGroup or Universe object containing the atoms
        style: str
            The style of the trajectory. Default is 'spheres'
        name: str
            The name of the trajectory. Default is 'atoms'
        location: np.ndarray or list
            The location of the trajectory. Default is None
        rotation: np.ndarray or list
            The rotation of the trajectory. Default is None
        scale: np.ndarray or list
            The scale of the trajectory. Default is None
        color: str
            The color of the trajectory. Default is 'default'
        material: str
            The material of the trajectory. Default is 'default'

        Returns:
        --------
        trajectory: Trajectory
            The created `Trajectory` object
        """
        trajectory = Trajectory(
                            atomgroup=universe,
                            style=style,
                            name=name,
                            color=color,
                            location=location,
                            rotation=rotation,
                            scale=scale,
                            material=material,
                            )
        self.trajectories.append(trajectory)
        return trajectory

    @functools.wraps(trajectory)
    def molecule(self, *args, **kwargs):
        logger.warning("molecule() is deprecated. Use trajectory() instead.")
        return self.trajectory(*args, **kwargs)


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
        mol_atoms = Trajectory(atomgroup=AtomGroup(atom1 + atom2),
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
        self.trajectories.append(mol_atoms)
        return mol_atoms, line
        
        
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