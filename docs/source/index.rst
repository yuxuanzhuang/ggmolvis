.. ggmolvis documentation master file, created by
   sphinx-quickstart on Thu Mar 15 13:55:56 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to ggmolvis's documentation!
=========================================================

As molecular simulations and other molecular entities continue to grow in scale and complexity, 
there is an increasing demand for flexible, high-quality visualization tools capable of handling 
these complex datasets efficiently. Molecular Nodes enables quick import 
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

GGMolVis provides a streamlined, object-oriented API for
molecular visualization, inspired by the design philosophy of Matplotlib and ggplot2. 
Such an API would enable users to achieve automated, customizable visualizations with just a 
few lines of code, and allow platforms like Molecular Dynamics Data Bank (MDDB) and Protein Data 
Bank (PDB) to produce large-scale renderings of their molecular archive, leveraging the most 
up-to-date rendering technologies that Blender provides.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   getting_started
   api

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`