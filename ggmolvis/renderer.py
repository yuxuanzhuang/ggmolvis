import bpy
import tempfile
from PIL import Image
from IPython.display import display

class Renderer:
    def __init__(self, resolution=(640, 360), filepath=None):
        # Set the render resolution
        bpy.context.scene.render.resolution_x = resolution[0]
        bpy.context.scene.render.resolution_y = resolution[1]
        bpy.context.scene.render.filepath = filepath or tempfile.NamedTemporaryFile(suffix=".png").name
        bpy.context.scene.render.image_settings.file_format = 'PNG'
    
    def render(self):
        # Render the scene
        bpy.ops.render.render(write_still=True)
    
    def display_in_notebook(self):
        # Load the image from the filepath
        image = Image.open(bpy.context.scene.render.filepath)
        
        # Display the image in Jupyter notebook
        display(image)
