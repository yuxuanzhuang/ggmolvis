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

class Renderer:
    def __init__(self, resolution=(1280, 720), filepath=None):
        # Set the render resolution
        bpy.context.scene.render.resolution_x = resolution[0]
        bpy.context.scene.render.resolution_y = resolution[1]
        bpy.context.scene.render.filepath = filepath or tempfile.NamedTemporaryFile(suffix=".PNG").name
        print('Rendering to:', bpy.context.scene.render.filepath)
        bpy.context.scene.render.image_settings.file_format = 'PNG'
    
    def render(self):
        # Render the scene
        bpy.ops.render.render(write_still=True)
    
    def display_in_notebook(self):
        # Load the image from the filepath
        image = Image.open(bpy.context.scene.render.filepath)
        
        # Display the image in Jupyter notebook
        display(image)


class MovieRenderer(Renderer):
    def __init__(self, resolution=(640, 360), filepath=None):
        bpy.context.scene.render.resolution_x = resolution[0]
        bpy.context.scene.render.resolution_y = resolution[1]
        bpy.context.scene.render.filepath = filepath or tempfile.NamedTemporaryFile(suffix=".mp4").name
        print('Rendering to:', bpy.context.scene.render.filepath)
        bpy.context.scene.render.image_settings.file_format = 'FFMPEG'
        bpy.context.scene.render.ffmpeg.format = 'MPEG4'
        bpy.context.scene.render.ffmpeg.codec = 'H264'
        bpy.context.scene.render.ffmpeg.constant_rate_factor = 'HIGH'
    
    def render(self):
        # Render the movie
        bpy.ops.render.render(animation=True)

    def display_in_notebook(self):
        return bpy.context.scene.render.filepath
