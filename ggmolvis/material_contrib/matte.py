"""
ggmolvis_matte: A matte material.
"""


import bpy


def material_generator():
    material = bpy.data.materials.new(name="ggmolvis_matte")
    material.use_nodes = True
    material.node_tree.nodes.clear()
    principled = material.node_tree.nodes.new("ShaderNodeBsdfPrincipled")
    principled.inputs["Base Color"].default_value = (0.8, 0.8, 0.8, 1.0)
    principled.inputs["Metallic"].default_value = 0.0
    principled.inputs["Roughness"].default_value = 0.9
    try:
        principled.inputs["Specular"].default_value = 1.0
    except KeyError:
        principled.inputs["Specular IOR Level"].default_value = 1.0
    for ele in principled.inputs:
        print(ele)
    principled.inputs["Subsurface Scale"].default_value = 5.5
    principled.inputs["Subsurface Radius"].default_value = (5.5, 5.5, 5.5)
    principled.inputs["Diffuse Roughness"].default_value = 1.0
    attribute = material.node_tree.nodes.new("ShaderNodeAttribute")
    attribute.attribute_name = "Color"
    output = material.node_tree.nodes.new("ShaderNodeOutputMaterial")
    material.node_tree.links.new(attribute.outputs["Color"], principled.inputs["Base Color"])
    material.node_tree.links.new(principled.outputs["BSDF"], output.inputs["Surface"])
