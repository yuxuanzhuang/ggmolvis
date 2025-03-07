import os


import bpy
from PIL import Image, ImageDraw, ImageFont


def _set_compositor_bg(rgba: tuple[float, float, float, float]):
    """
    Set a solid output background color using the compositor.

    Parameters
    ----------
    rgba: tuple
        The red, green, blue, and alpha float values to use
        for the solid background color mixed in by the
        compositor.
    """
    scene = bpy.context.scene
    scene.render.use_compositing = True
    scene.use_nodes = True
    nodes = scene.node_tree.nodes
    links = scene.node_tree.links
    nodes.clear()
    render_layers = nodes.new(type="CompositorNodeRLayers")
    composite = nodes.new(type="CompositorNodeComposite")
    color_node = nodes.new(type="CompositorNodeValue")
    mix_node = nodes.new(type="CompositorNodeMixRGB")
    color_node.outputs[0].default_value = 1.0
    mix_node.blend_type = 'MIX'
    mix_node.inputs[0].default_value = 1.0
    mix_node.inputs[1].default_value = rgba
    links.new(render_layers.outputs['Alpha'], mix_node.inputs[0])
    links.new(render_layers.outputs['Image'], mix_node.inputs[2])
    links.new(mix_node.outputs[0], composite.inputs[0])


def _create_frame_image(frame_number: int,
                        width: int,
                        height: int,
                        text: str):
    img = Image.new('RGBA', (width, height), (0, 0, 0, 0))
    draw = ImageDraw.Draw(img)
    font = ImageFont.truetype("/System/Library/Fonts/Geneva.ttf", 80)
    text_bbox = draw.textbbox((0, 0), text, font=font)
    text_width = text_bbox[2] - text_bbox[0]
    text_height = text_bbox[3] - text_bbox[1]
    x = (width - text_width) // 2
    y = height - text_height - 40
    draw.text((x, y), text, font=font, fill=(255, 255, 255, 255))
    os.makedirs("frame_overlays", exist_ok=True)
    out_path = os.path.join("frame_overlays", f"frame_{frame_number}.png")
    img.save(out_path, "PNG")
    return out_path

def _composit_frame_label(img_path):
    scene = bpy.context.scene
    scene.use_nodes = True
    tree = scene.node_tree
    links = tree.links
    tree.nodes.clear()
    render_layers = tree.nodes.new('CompositorNodeRLayers')
    image_node = tree.nodes.new('CompositorNodeImage')
    alpha_over = tree.nodes.new('CompositorNodeAlphaOver')
    composite = tree.nodes.new('CompositorNodeComposite')
    image_node.image = bpy.data.images.load(img_path)
    links.new(render_layers.outputs['Image'], alpha_over.inputs[2])
    links.new(image_node.outputs['Image'], alpha_over.inputs[1])
    links.new(alpha_over.outputs['Image'], composite.inputs['Image'])
