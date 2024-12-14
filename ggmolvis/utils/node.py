import bpy
import molecularnodes as mn
from molecularnodes.blender.nodes import swap, assign_material, add_custom

def extract_mn_node(object: bpy.types.Object):
    """Extract the node group from the MN object."""
    nodes_mn = object.modifiers["MolecularNodes"].node_group
    nodes = nodes_mn.nodes
    links = nodes_mn.links
    return nodes, links


def set_selection(object, selection_name):
    """Set the selection."""
    nodes, links = extract_mn_node(object)
    name_atoms = selection_name
    named_attr_node = nodes.new(type='GeometryNodeInputNamedAttribute')
    named_attr_node.name = name_atoms
    named_attr_node.inputs["Name"].default_value = name_atoms

    # get current style node
    style_lists = [s for s in nodes.keys()  if s.startswith("Style")]
    style_node = nodes[style_lists[0]]
    links.new(named_attr_node.outputs["Attribute"],
              style_node.inputs["Selection"])


def swap_style(object, style):
    """Swap the style."""
    nodes, links = extract_mn_node(object)
    style_lists = [s for s in nodes.keys()  if s.startswith("Style")]
    style_node = nodes[style_lists[0]]
    swap(style_node, style)
    
    # add DSSP node before that if cartoon style
    # kinda hacky
    if style == 'Style Cartoon' and 'Topology DSSP' not in nodes.keys():
        nodes_mn = object.modifiers["MolecularNodes"].node_group
        add_custom(nodes_mn, name='Topology DSSP')
        dssp_node = nodes_mn.nodes['Topology DSSP']

        # link dssp before style node
        links.new(dssp_node.outputs["Atoms"], style_node.inputs["Atoms"])
        # after Set Color
        links.new(dssp_node.inputs["Atoms"], nodes['Set Color'].outputs["Atoms"])

def set_mn_color(object, color):
    nodes, links = extract_mn_node(object)
    nodes_mn = object.modifiers["MolecularNodes"].node_group

    color_node = nodes_mn.nodes['Set Color']
    
    # cannot remove link with
    if color_node.inputs["Color"].links:
        nodes_mn.links.remove(color_node.inputs["Color"].links[0])
    color_node.inputs['Color'].default_value = color

def set_mn_material(object, material):
    """Set the material."""
    nodes, links = extract_mn_node(object)
    style_lists = [s for s in nodes.keys()  if s.startswith("Style")]
    style_node = nodes[style_lists[0]]
    assign_material(style_node, material)
