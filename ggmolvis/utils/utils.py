import numpy as np
from functools import wraps

import bpy
import mathutils
from molecularnodes.blender.nodes import styles_mapping as mol_styles_mapping


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
    # https://blender.stackexchange.com/questions/68834/recreate-to-track-quat-with-two-vectors-using-python
    camera_direction = camera_position - target_position
    camera_direction = camera_direction / np.linalg.norm(camera_direction)
    camera_right = np.cross(np.array([0.0, 0.0, 1.0]), camera_direction)
    camera_right = camera_right / np.linalg.norm(camera_right)
    camera_up = np.cross(camera_direction, camera_right)
    camera_up = camera_up / np.linalg.norm(camera_up)
    rotation_transform = np.zeros((4, 4))
    rotation_transform[0, :3] = camera_right
    rotation_transform[1, :3] = camera_up
    rotation_transform[2, :3] = camera_direction
    rotation_transform[-1, -1] = 1
    translation_transform = np.eye(4)
    translation_transform[:3, -1] = - camera_position
    look_at_transform = np.matmul(rotation_transform, translation_transform)
    return look_at_transform

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
}

AVAILABLE_MATERIALS = materials_mapping.keys()

AVAILABLE_STYLES = ["default"]
MOL_AVAILABLE_STYLES = mol_styles_mapping.keys()


def validate_properties(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        # Extract keyword arguments
        style = kwargs.get('style', 'spheres')
        color = kwargs.get('color', 'default')
        material = kwargs.get('material', 'default')
        
        # Validate style
        if style not in MOL_AVAILABLE_STYLES:
            raise ValueError(f"Invalid style '{style}'. Valid options are: {MOL_AVAILABLE_STYLES}")
                
        # Validate color
        #if color not in VALID_COLORS:
        #    raise ValueError(f"Invalid color '{color}'. Valid options are: {VALID_COLORS}")
        
        # Validate material
        if material not in AVAILABLE_MATERIALS:
            raise ValueError(f"Invalid material '{material}'. Valid options are: {AVAILABLE_MATERIALS}")
        
        # Call the original function
        return func(*args, **kwargs)
    
    return wrapper