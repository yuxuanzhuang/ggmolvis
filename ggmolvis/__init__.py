"""
ggmolvis
Molecular visualization with Blender
"""

import os
import shutil
import tempfile
import uuid
from importlib.metadata import version
import bpy
from bpy.app.handlers import frame_change_pre
from bpy.app.handlers import persistent
import molecularnodes as mn
from ggmolvis.utils import suppress_blender_output
import atexit

from loguru import logger

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

# we will save the current session to the temporary directory
try:
    dest_dir = f"{tempfile.gettempdir()}/{uuid.uuid4()}"
    os.makedirs(dest_dir, exist_ok=True)
    dest_path = os.path.join(dest_dir, base_name)
    shutil.copy(mn_template_file, dest_path)
except Exception as e:
    # fallback to the current directory if tempdir fails
    # create a unique filename in the current directory
    base_name = f'{uuid.uuid4()}.blend'
    dest_dir = '.'
    dest_path = os.path.join(dest_dir, base_name)
    shutil.copy(mn_template_file, dest_path)

logger.debug(f"Blend file is saved to: {dest_path}")

with suppress_blender_output():
    bpy.ops.wm.open_mainfile(filepath=dest_path)

# update every artist when the frame changes
@persistent
def update_frame(scene):
    # This is necessary to e.g. get the correct object
    # dimensions after frame change during rendering
    # see: https://blender.stackexchange.com/questions/61635/object-dimensions-changing-but-not-changing
    bpy.context.view_layer.update() 
    
    for artist in SESSION._ggmolvis:
        artist._update_frame(scene.frame_current)

frame_change_pre.append(update_frame)

# add visualize function to AnalysisBase and GroupBase so that
# AtomGroup and ResidueGroup can be visualized with `.visualize()`
from ggmolvis.analysis import Visualizer

from .ggmolvis import GGMolVis
GGMOLVIS = GGMolVis()

def cleanup_function():
#    print("Saving the current session to", dest_path)
#    bpy.ops.wm.save_as_mainfile(filepath=dest_path)
    from molecularnodes import unregister
    frame_change_pre.remove(update_frame)
    try:
        unregister()
    except Exception as e:
        print(e)

atexit.register(cleanup_function)
