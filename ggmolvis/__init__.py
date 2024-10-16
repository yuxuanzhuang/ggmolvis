"""
ggmolvis
Molecular visualization with Blender
"""

import os
import shutil
import tempfile
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
base_name = 'ggmolvis.blend'
name, ext = os.path.splitext(base_name)
dest_path = os.path.join('.', base_name)
dest_dir = '.'
if not os.access(dest_dir, os.W_OK):
    print(f"Directory {dest_dir} is not writable. Using a temporary directory.")
    dest_dir = tempfile.gettempdir()

count = 1
while os.path.exists(dest_path):
    dest_path = os.path.join(dest_dir, f"{name}_{count}{ext}")
    count += 1

shutil.copy(mn_template_file, dest_path)


bpy.ops.wm.open_mainfile(filepath=dest_path)


@persistent
def update_frame(scene):
    for artist in SESSION._ggmolvis:
        artist._update_frame(scene.frame_current)
#        try:
#            artist._update_frame(scene.frame_current)
#        except:
            # if the object under the artist is deleted,
            # remove the artist
#            SESSION._ggmolvis.remove(artist)
frame_change_post.append(update_frame)

# add visualize function to AnalysisBase
from .analysis import Visualizer

from .ggmolvis import GGMolVis
GGMOLVIS = GGMolVis()