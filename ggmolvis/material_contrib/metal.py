"""
ggmolvis_metal: A metallic, shiny material.
"""


import bpy


def material_generator():
    material = bpy.data.materials.new(name="ggmolvis_metal")
    material.use_nodes = True
    material.node_tree.nodes.clear()
    principled = material.node_tree.nodes.new("ShaderNodeBsdfPrincipled")
    principled.inputs["Base Color"].default_value = (0.8, 0.8, 0.8, 1.0)  # Silvery color
    principled.inputs["Metallic"].default_value = 0.95
    principled.inputs["Roughness"].default_value = 0.1
    output = material.node_tree.nodes.new("ShaderNodeOutputMaterial")
    material.node_tree.links.new(principled.outputs["BSDF"], output.inputs["Surface"])
