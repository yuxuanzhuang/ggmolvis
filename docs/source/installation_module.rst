==============================================
GGMolVis Installation Guide as a Python Module
==============================================

This guide provides instructions for installing GGMolVis as a Python module, along with its dependencies.
Please note that this method of installation does not involve running Blender's GUI; instead, you will interact with GGMolVis through Python.

Installation Steps
==================

1. **Install GGMolVis**

   First, clone the GGMolVis repository from GitHub:

   .. code-block:: bash

      git clone git@github.com:yuxuanzhuang/ggmolvis.git
      cd ggmolvis

   Next, create and activate the conda environment specified in the repository:

   .. code-block:: bash

      conda env create -f devtools/conda-envs/test_env.yaml
      conda activate ggmolvis-test

   Then, install GGMolVis in editable mode:

   .. code-block:: bash

      python -m pip install -e .

   Verify the installation by checking the version number:

   .. code-block:: bash

      python -c "import ggmolvis; print(ggmolvis.__version__)"

   .. note::
      Currently, only Python 3.11 is supported.

2. **(Optional) Install a Different Version of MolecularNodes**

   If you need a development or alternative version of MolecularNodes, you can install it as follows:

   .. code-block:: bash

      git clone git@github.com:BradyAJohnston/MolecularNodes.git
      cd MolecularNodes
      python -m pip install -e .

Verification
============

1. Start JupyterLab using the ggmolvis-test environment:

   .. code-block:: bash

      jupyter lab

2. When JupyterLab opens, select the `ggmolvis-test` kernel in the Launcher or from the **Kernel** menu. No Blender GUI will appear during usage.

3. In a Jupyter notebook cell, run the following commands to confirm that GGMolVis is installed correctly:

   .. code-block:: python

      import ggmolvis
      print(ggmolvis.__version__)
      print(ggmolvis.__file__)

If GGMolVis is properly installed, these commands will print the currently installed version and the file path for the GGMolVis module.

You are now ready to use GGMolVis as a Python module!
