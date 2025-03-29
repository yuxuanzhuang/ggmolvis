"""
This module contains the classes for the different shapes that can be created in the scene.

Classes
=======
.. autoclass:: Shape
    :members:
.. autoclass:: Line
    :members:
"""

import bpy
from abc import ABC, abstractmethod

from molecularnodes.blender import coll
import numpy as np
from typing import Tuple, List, Union

from .base import SceneObject
from ..properties import Color, Material
from ..utils import convert_list_to_array, look_at, lerp


class Shape(SceneObject):
    def __init__(
        self,
        shape_type,
        name=None,
        location=None,
        rotation=None,
        scale=None,
        color="black",
        material="backdrop",
    ):
        self.shape_type = shape_type
        super().__init__(
            name=name,
            location=location,
            rotation=rotation,
            scale=scale,
            color=color,
            material=material,
        )


class Line(Shape):
    def __init__(
        self,
        start_points: Union[List[Tuple[float, float, float]], np.ndarray],
        end_points: Union[List[Tuple[float, float, float]], np.ndarray],
        name=None,
        location=None,
        rotation=None,
        scale=None,
        color="black",
        material="backdrop",
    ):
        self.start_points = convert_list_to_array(start_points)
        self.end_points = convert_list_to_array(end_points)
        super().__init__(
            shape_type="line",
            name=name,
            location=location,
            rotation=rotation,
            scale=scale,
            color=color,
            material=material,
        )

    def _create_object(self):
        line_data = bpy.data.curves.new(name=self.name, type="CURVE")
        line_data.dimensions = "3D"
        self.line_object = bpy.data.objects.new(self.name, line_data)
        coll.mn().objects.link(self.line_object)

        line = line_data.splines.new("POLY")
        self.line = line
        line.points.add(1)
        line.resolution_u = 4
        line.use_cyclic_u = False
        line.use_endpoint_u = True
        line.use_endpoint_v = True
        line.use_smooth = False

        line_data.bevel_depth = 0.004
        line_data.bevel_resolution = 10

        self._update_frame(bpy.context.scene.frame_current)
        return self.line_object

    def _update_frame(self, frame):
        object = self.object
        start_point, end_point = self._get_points_for_frame(frame)
        object.data.splines[0].points[0].co = (
            start_point[0],
            start_point[1],
            start_point[2],
            1.0,
        )
        object.data.splines[0].points[1].co = (
            end_point[0],
            end_point[1],
            end_point[2],
            1.0,
        )

    def _get_points_for_frame(self, frame: int) -> Tuple[float, float, float]:
        """Retrieve the coordinates for a specific frame"""

        if self.subframes == 0:
            frame_a = frame
        else:
            frame_a = int(frame / (self.subframes + 1))

        # get the next frame
        frame_b = frame_a + 1
        if frame_b >= self.start_points.shape[0]:
            return (
                self.start_points[-1] * self.world_scale,
                self.end_points[-1] * self.world_scale,
            )

        locations_a = []
        locations_b = []
        for points in [self.start_points, self.end_points]:
            if points.ndim == 2:
                locations_a.append(points[frame_a] * self.world_scale)
                locations_b.append(points[frame_b] * self.world_scale)
            elif points.ndim == 1:
                locations_a.append(points * self.world_scale)
                locations_b.append(points * self.world_scale)
            else:
                raise ValueError("Invalid transformation coordinates")

        if self.subframes > 0:
            fraction = frame % (self.subframes + 1) / (self.subframes + 1)

            # interpolate between the two sets of positions
            locations = lerp(locations_a, locations_b, t=fraction)
        else:
            locations = locations_a

        return locations
