import pytest
import bpy
from functools import wraps

def cleanup(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        try:
            result = func(*args, **kwargs)  # Call the original test function
            return result
        finally:
            ggmv = kwargs.get("ggmv") or args[1]  # Extract ggmv safely
            _clean_up(ggmv)
    return wrapper

def _clean_up(ggmv):
    """Remove all objects after a test."""
    for name, entity in ggmv._artists_dict.items():
        for artist in entity:
            try:
                bpy.data.objects.remove(artist.object)
            except Exception:
                pass
        ggmv._artists_dict[name] = []