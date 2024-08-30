"""
ggmolvis
Molecular visualization with Blender
"""

# Add imports here
from importlib.metadata import version
import bpy
from bpy.app.handlers import frame_change_post, load_post, save_post
from bpy.app.handlers import persistent

import molecularnodes as mn

__version__ = version("ggmolvis")

# sync to the molecularnodes session
SESSION = mn.session.get_session()

@persistent
def update_frame(scene):
    for artist in SESSION._ggmolvis:
        artist.update_frame(scene.frame_current)
frame_change_post.append(update_frame)