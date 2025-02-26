import bpy
from ggmolvis import compositor


import pytest
from numpy.testing import assert_allclose



@pytest.mark.parametrize("rgba", [
    (0.0, 0.0, 0.0, 1.0), # black
    (0.3, 0.1, 0.9, 1.0), # mixture color
    # compositor allows > 1.0 intensities
    (5.0, 0.0, 0.0, 1.0),
])
def test_basic_background_setting(rgba):
    # verify that a simple scene has its composited
    # background "color" properly set
    bpy.ops.scene.new(type='NEW')
    bpy.ops.mesh.primitive_cube_add(size=2, location=(0, 0, 0))
    compositor._set_compositor_bg(rgba)
    nodes = bpy.context.scene.node_tree.nodes
    for n in nodes:
        if "MixRGB" in n.bl_idname:
            actual = n.inputs[1].default_value
            assert_allclose(actual, rgba)
    bpy.ops.scene.delete()
