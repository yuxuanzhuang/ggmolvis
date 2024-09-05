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