=====================================
Getting Started with GGMolVis
=====================================

As molecular simulations and other molecular entities continue to grow in scale and complexity, 
there is an increasing demand for flexible, high-quality visualization tools capable of handling 
these complex datasets efficiently. Molecular Nodes (About – Molecular Nodes ) enables quick import 
and visualization of structural biology data inside of Blender, a free and open-source 3D computer 
graphics software tool. Blender has proven to be an innovative platform, providing versatile and 
industry-quality renderings of static molecules, molecular dynamics (MD) simulations (with the support 
of MDAnalysis), cryo-EM density maps, and EM tomography. It has become the go-to package for 
both experimentalists and computational scientists, widely used to create cover images and 
striking illustrations for academic and outreach purposes, including recent high-profile projects 
like the AlphaFold 3 release video by Google DeepMind.

Molecular Nodes is currently designed with a graphic user interface (GUI) workflow being the main input method.
This limits the ability and uptake by computations chemists and biophysicists who quite
often work in headless computational environments. Reusability and reproducibility
are also concerns in an effort to include Molecular Nodes in more computational pipelines
with FAIR principles in mind.

GGMolVis is designed to offer a streamlined, object-oriented API for
molecular visualization, inspired by the design philosophy of Matplotlib and ggplot2. 
Such an API would enable users to achieve automated, customizable visualizations with just a 
few lines of code, and allow platforms like Molecular Dynamics Data Bank (MDDB) and Protein Data 
Bank (PDB) to produce large-scale renderings of their molecular archive, leveraging the most 
up-to-date rendering technologies that Blender provides.

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
