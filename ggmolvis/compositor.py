import bpy


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
