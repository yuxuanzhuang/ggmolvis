=====================================
Getting Started with GGMolVis
=====================================

GGMolVis is a tool that can be used in two different ways:

1. **Within Blender**: Integrate GGMolVis directly into Blender's environment.

2. **As a Standalone Python Module**: Use GGMolVis without any graphical user interface, suitable for computational workflows and integration with Python-based data processing pipelines. Note right now this approach is pretty noisy during rendering.


New to Blender?
===============

`Blender <https://www.blender.org/>`_ is an industry-grade 3D modeling and animation software. 
Check out the `Molecular Nodes tutorials <https://bradyajohnston.github.io/MolecularNodes/tutorials/>`_ to get started with
using Molecular Nodes Blender. You can also find a wealth of tutorials on YouTube and other online resources.

Alternatively, you can also use GGMolVis as a Python module to create visualizations without the need for Blender's GUI.


Using GGMolVis with Blender
===========================

If you plan to leverage Blender’s 3D visualization capabilities, 
follow the **Blender Installation Guide**:

- **Features**:

  - Full 3D visualization environment

  - Interactive molecule manipulation

  - Integration with other Blender add-ons.

- **Requirements**:

  - Blender 4.3 or later

  - Installation of MolecularNodes and BNotebooks

Check out the `GGMolVis Installation Guide with Blender <installation_blender.html>`_ for step-by-step instructions.


Using GGMolVis as a Python Module
=================================

For workflows that do not require Blender’s GUI, or where you prefer 
headless computation and scripting, install GGMolVis as a Python module.

- **Features**:

  - No GUI required

  - Easily integrate with Python-based tooling, Jupyter notebooks, and more

  - Ideal for batch processing or computational analysis

- **Requirements**:

  - Python 3.11

Refer to the `GGMolVis Installation Guide as a Python Module <installation_module.html>`_ for details on setting up your environment.


Installation Guides
===================

Choose the installation method that best fits your workflow and requirements.

.. toctree::
   :maxdepth: 2

   installation_blender
   installation_module
   

Usage Examples
==============

Once you have GGMolVis installed, you can start using it to show off your favorite molecules and trajectories!

See the example notebooks in the `examples` directory of the GGMolVis repository for inspiration.

.. toctree::
   :maxdepth: 2

   examples/ggmolvis_basic.ipynb


Getting Help
============

If you encounter issues or need help:

- Check the `GitHub Issues <https://github.com/yuxuanzhuang/ggmolvis/issues>`_ page for known problems or to report new ones.

- Refer to the `FAQ` page if available.

- Ask questions in the GGMolVis community channels or forums.
