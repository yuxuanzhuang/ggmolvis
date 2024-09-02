"""
ggmolvis
Molecular visualization with Blender
"""

import os
from importlib.metadata import version
import bpy
from bpy.app.handlers import frame_change_post
from bpy.app.handlers import persistent
import molecularnodes as mn

__version__ = version("ggmolvis")

# sync to the molecularnodes session
try:
    SESSION = mn.session.get_session()
# if the module is started without Blender running
except AttributeError:
    from molecularnodes import register
    register()
    SESSION = mn.session.get_session()

mn_template_file = os.path.join(
        os.path.abspath(mn.utils.ADDON_DIR),
        "assets", "template", "startup.blend"
    )
bpy.ops.wm.open_mainfile(filepath=mn_template_file)


@persistent
def update_frame(scene):
    for artist in SESSION._ggmolvis:
        artist.update_frame(scene.frame_current)
#        try:
#            artist.update_frame(scene.frame_current)
#        except:
            # if the object under the artist is deleted,
            # remove the artist
#            SESSION._ggmolvis.remove(artist)
frame_change_post.append(update_frame)

# add visualize function to AnalysisBase
from .analysis import Visualizer
