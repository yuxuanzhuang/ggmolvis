"""
Utility functions for working with Blender objects and data.

Functions:
==========
.. autosummary::
    convert_list_to_array
    look_at
    quaternion_to_euler
    euler_to_quaternion
    validate_properties
    suppress_blender_output
"""

import numpy as np
from functools import wraps

import bpy
import mathutils
from molecularnodes.blender.nodes import styles_mapping as mol_styles_mapping
import os
import sys
from contextlib import contextmanager


def lerp(a, b, t) -> np.ndarray:
    """
    Linearly interpolates between two arrays a and b by a factor of t.
    Parameters:
    a (np.ndarray): The starting array.
    b (np.ndarray): The ending array.
    t (float or np.ndarray): The interpolation factor, where 0 <= t <= 1.
                             If t is an array, it should be broadcastable to the shape of a and b.
    Returns:
    np.ndarray: The interpolated array.
    """

    return np.add(a, np.multiply(np.subtract(b, a), t))


def convert_list_to_array(value):
    if isinstance(value, np.ndarray):
        # Check if it's a 1D array with size 3, convert to 2D array with shape (1, 3)
        if value.ndim == 1 and value.size == 3:
            return value
        # Check if it's a 2D array with shape (N, 3)
        elif value.ndim == 2 and value.shape[1] == 3:
            return value
        else:
            raise ValueError("NumPy array must have shape (3,) or (N, 3)")
    
    elif isinstance(value, list) or isinstance(value, tuple):
        # Convert list to NumPy array
        value = np.asarray(value)
        if value.ndim == 1 and value.size == 3:
            return value
        elif value.ndim == 2 and value.shape[1] == 3:
            return value
        else:
            raise ValueError("List of coordinates must have shape (N, 3)")

    raise ValueError("Coordinates must be a NumPy array or a list of tuples/lists with shape (N, 3)")

def look_at(camera_position, target_position):
    """
    Calculate the rotation matrix for a camera to point towards a target position.

    Parameters:
    camera_position (mathutils.Vector): The current position of the camera.
    target_position (mathutils.Vector): The position that the camera should look at.

    Returns:
    mathutils.Matrix: A 4x4 rotation matrix for the camera.
    """
    # Calculate the forward direction (from camera to target)
    target_position = mathutils.Vector(target_position)
    camera_position = mathutils.Vector(camera_position)
    forward = (target_position - camera_position).normalized()

    # Define the up vector (typically the global Z axis in Blender)
    up = mathutils.Vector((0.0, 0.0, 1.0))

    # Calculate the right vector (perpendicular to forward and up)
    right = forward.cross(up).normalized()

    # Recalculate the up vector to ensure it's orthogonal to the forward and right vectors
    up = right.cross(forward).normalized()

    # Create a 3x3 rotation matrix from the right, up, and forward vectors
    rotation_matrix = mathutils.Matrix((
        right,
        up,
        -forward  # Blender's camera looks along the -Z axis, so we negate the forward vector
    )).transposed()  # Transpose because Blender expects column-major order

    # Convert the 3x3 rotation matrix to a 4x4 matrix (no translation, just rotation)
    rotation_matrix_4x4 = rotation_matrix.to_4x4()

    _, rot, _ = rotation_matrix_4x4.decompose()

    return rot

def quaternion_to_euler(quaternion):
    """
    Convert quaternion to Euler angles.

    Parameters:
    quaternion (mathutils.Quaternion): A quaternion representing the rotation.

    Returns:
    tuple: A tuple of three floats representing the Euler angles (in radians).
    """
    # Create a mathutils.Euler object from the input quaternion
    euler = quaternion.to_euler('XYZ')
    
    # Convert the Euler angles to a tuple
    euler_angles = (euler.x, euler.y, euler.z)
    
    return euler_angles


def euler_to_quaternion(euler_angles):
    """
    Convert Euler angles to quaternion.

    Parameters:
    euler_angles (tuple): A tuple of three floats representing the Euler angles (in radians).

    Returns:
    mathutils.Quaternion: A quaternion representing the same rotation as the Euler angles.
    """
    # Create a mathutils.Euler object from the input angles
    euler = mathutils.Euler(euler_angles, 'XYZ')
    
    # Convert the Euler angles to a quaternion
    quaternion = euler.to_quaternion()
    
    return quaternion



materials_mapping = {
    "default": "MN Default",
    "flat": "MN Flat Outline",
    "squishy": "MN Squishy",
    "transparent": "MN Transparent Outline",
    "ambient": "MN Ambient Occlusion",
    # backdrop for all other shapes
    "backdrop": "Backdrop",
    # user-contributed add-in materials
    "metal": "ggmolvis_metal",
    "matte": "ggmolvis_matte",
}


AVAILABLE_MATERIALS = materials_mapping.keys()

AVAILABLE_STYLES = ["default"]
MOL_AVAILABLE_STYLES = mol_styles_mapping.keys()
AVAILABLE_STYLES.extend(MOL_AVAILABLE_STYLES)


def validate_properties(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        # Extract keyword arguments
        # get all style kwargs
        style_kwargs = {k: v for k, v in kwargs.items() if 'style' in k}
        color_kwargs = {k: v for k, v in kwargs.items() if 'color' in k}
        material_kwargs = {k: v for k, v in kwargs.items() if 'material' in k}
        
        # Validate style
        for style in style_kwargs.values():
            if style not in AVAILABLE_STYLES:
                raise ValueError(f"Invalid style '{style}'. "
                                 f"Valid options are: {AVAILABLE_STYLES}")     
        # Validate color
        #for color in color_kwargs.values():
        #    if color not in AVAILABLE_COLORS:
        #        raise ValueError(f"Invalid color '{color}'. Valid options are: {AVAILABLE_COLORS}")
        
        # Validate material
        for material in material_kwargs.values():
            if material not in AVAILABLE_MATERIALS:
                raise ValueError(f"Invalid material '{material}'. "
                                 f"Valid options are: {AVAILABLE_MATERIALS}")
        
        # Call the original function
        return func(*args, **kwargs)
    
    return wrapper



# Context manager to suppress stdout and stderr
@contextmanager
def suppress_blender_output():
    """Context manager to suppress Blender's low-level stdout and stderr output."""
    # Open /dev/null for discarding outputs
    null_fds = [os.open(os.devnull, os.O_RDWR) for _ in range(2)]
    # Save the original stdout and stderr file descriptors
    saved_fds = [os.dup(1), os.dup(2)]
    try:
        # Redirect stdout and stderr to /dev/null
        os.dup2(null_fds[0], 1)  # Redirect stdout (fd=1)
        os.dup2(null_fds[1], 2)  # Redirect stderr (fd=2)
        yield  # Execute code within the context
    finally:
        # Restore original stdout and stderr
        os.dup2(saved_fds[0], 1)
        os.dup2(saved_fds[1], 2)
        # Close file descriptors
        for fd in null_fds + saved_fds:
            os.close(fd)
