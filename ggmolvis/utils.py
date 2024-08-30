import numpy as np

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
