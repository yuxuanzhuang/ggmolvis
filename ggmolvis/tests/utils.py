import bpy
from functools import wraps
from ggmolvis.ggmolvis import GGMolVis

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
            # remove blender object
            try:
                bpy.data.objects.remove(artist.object)
            except Exception:
                pass
            # remove MN object
            try:
                ggmv._session.remove(artist.trajectory.uuid)
            except Exception:
                pass
    bpy.ops.wm.read_homefile(app_template="")
    ggmv._session.prune()
    ggmv._session._ggmolvis = set()
    ggmv._initialized = False
    ggmv = GGMolVis()
