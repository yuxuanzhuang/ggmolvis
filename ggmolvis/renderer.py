"""
This module provides a simple renderer class that can be used to render images and movies in Blender.

The `Renderer` class provides a simple interface to render images in Blender.
It allows you to set the resolution of the render and the output file path.
You can then call the `render` method to render the image and the
`display_in_notebook` method to display the rendered image in a Jupyter notebook.

Normally, you can use `.render()` function of the `Sceneobject` to render the scene instead of using this class.

Classes
=======

.. autoclass:: Renderer
    :members:

.. autoclass:: MovieRenderer
    :members:
"""
import bpy
import tempfile
from PIL import Image
from IPython.display import display, Video
from loguru import logger

class Renderer:
    def __init__(self, filepath=None):
        bpy.context.scene.render.filepath = filepath or tempfile.NamedTemporaryFile(suffix=".PNG").name
        logger.info(f'Rendering to: {bpy.context.scene.render.filepath}')
    
    def render(self):
        # Render the scene
        bpy.ops.render.render(write_still=True)
    
    def display_in_notebook(self):
        # Load the image from the filepath
        image = Image.open(bpy.context.scene.render.filepath)
        
        # Display the image in Jupyter notebook
        display(image)


class MovieRenderer(Renderer):
    def __init__(self, filepath=None):
        bpy.context.scene.render.filepath = filepath or tempfile.NamedTemporaryFile(suffix=".mp4").name
        logger.info(f'Rendering to: {bpy.context.scene.render.filepath}')
    
    def render(self):
        # Render the movie
        bpy.ops.render.render(animation=True)

    def display_in_notebook(self):
        return bpy.context.scene.render.filepath
