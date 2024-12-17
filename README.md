ggmolvis
==============================
[//]: # (Badges)

| **Latest release** | [![Last release tag][badge_release]][url_latest_release] ![GitHub commits since latest release (by date) for a branch][badge_commits_since]  [![Documentation Status][badge_docs]][url_docs]|
| :----------------- | :------- |
| **Status**         | [![GH Actions Status][badge_actions]][url_actions] [![codecov][badge_codecov]][url_codecov] |
| **Community**      | [![License: GPL v2][badge_license]][url_license]  [![Powered by MDAnalysis][badge_mda]][url_mda]|

[badge_actions]: https://github.com/yuxuanzhuang/ggmolvis/actions/workflows/gh-ci.yaml/badge.svg
[badge_codecov]: https://codecov.io/gh/yuxuanzhuang/ggmolvis/branch/main/graph/badge.svg
[badge_commits_since]: https://img.shields.io/github/commits-since/yuxuanzhuang/ggmolvis/latest
[badge_docs]: https://readthedocs.org/projects/ggmolvis/badge/?version=latest
[badge_license]: https://img.shields.io/badge/License-GPLv2-blue.svg
[badge_mda]: https://img.shields.io/badge/powered%20by-MDAnalysis-orange.svg?logoWidth=16&logo=data:image/x-icon;base64,AAABAAEAEBAAAAEAIAAoBAAAFgAAACgAAAAQAAAAIAAAAAEAIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJD+XwCY/fEAkf3uAJf97wGT/a+HfHaoiIWE7n9/f+6Hh4fvgICAjwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACT/yYAlP//AJ///wCg//8JjvOchXly1oaGhv+Ghob/j4+P/39/f3IAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJH8aQCY/8wAkv2kfY+elJ6al/yVlZX7iIiI8H9/f7h/f38UAAAAAAAAAAAAAAAAAAAAAAAAAAB/f38egYF/noqAebF8gYaagnx3oFpUUtZpaWr/WFhY8zo6OmT///8BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgICAn46Ojv+Hh4b/jouJ/4iGhfcAAADnAAAA/wAAAP8AAADIAAAAAwCj/zIAnf2VAJD/PAAAAAAAAAAAAAAAAICAgNGHh4f/gICA/4SEhP+Xl5f/AwMD/wAAAP8AAAD/AAAA/wAAAB8Aov9/ALr//wCS/Z0AAAAAAAAAAAAAAACBgYGOjo6O/4mJif+Pj4//iYmJ/wAAAOAAAAD+AAAA/wAAAP8AAABhAP7+FgCi/38Axf4fAAAAAAAAAAAAAAAAiIiID4GBgYKCgoKogoB+fYSEgZhgYGDZXl5e/m9vb/9ISEjpEBAQxw8AAFQAAAAAAAAANQAAADcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAjo6Mb5iYmP+cnJz/jY2N95CQkO4pKSn/AAAA7gAAAP0AAAD7AAAAhgAAAAEAAAAAAAAAAACL/gsAkv2uAJX/QQAAAAB9fX3egoKC/4CAgP+NjY3/c3Nz+wAAAP8AAAD/AAAA/wAAAPUAAAAcAAAAAAAAAAAAnP4NAJL9rgCR/0YAAAAAfX19w4ODg/98fHz/i4uL/4qKivwAAAD/AAAA/wAAAP8AAAD1AAAAGwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAALGxsVyqqqr/mpqa/6mpqf9KSUn/AAAA5QAAAPkAAAD5AAAAhQAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADkUFBSuZ2dn/3V1df8uLi7bAAAATgBGfyQAAAA2AAAAMwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAB0AAADoAAAA/wAAAP8AAAD/AAAAWgC3/2AAnv3eAJ/+dgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA9AAAA/wAAAP8AAAD/AAAA/wAKDzEAnP3WAKn//wCS/OgAf/8MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIQAAANwAAADtAAAA7QAAAMAAABUMAJn9gwCe/e0Aj/2LAP//AQAAAAAAAAAA
[badge_release]: https://img.shields.io/github/release-pre/yuxuanzhuang/ggmolvis.svg
[url_actions]: https://github.com/yuxuanzhuang/ggmolvis/actions?query=branch%3Amain+workflow%3Agh-ci
[url_codecov]: https://codecov.io/gh/yuxuanzhuang/ggmolvis/branch/main
[url_docs]: https://ggmolvis.readthedocs.io/en/latest/?badge=latest
[url_latest_release]: https://github.com/yuxuanzhuang/ggmolvis/releases
[url_license]: https://www.gnu.org/licenses/gpl-2.0
[url_mda]: https://www.mdanalysis.org

**GGMolVis** is a Python package that provides a high-level interface to [MolecularNodes](https://github.com/BradyAJohnston/MolecularNodes) in
  [Blender](https://www.blender.org/) for molecular visualization.

It is inspired by the design patterns of ggplot2 and matplotlib. The goal is to create a stable Python API that enables users to visualize molecular trajectories (and potentially other entities in MN) with both automation and customization. Everything is designed to be dynamic during frame changes. The package also includes capabilities to visualize protein features and analysis results---such as distances, angles, and dihedrals---with texts, lines, and shapes.

![example](./rendering/example_ggmolvis.png)

### Features

- API Development: We will create a stable API for Molecular Nodes, empowering users to automate molecular rendering with minimal effort.

- Interactive Jupyter Integration: A Jupyter widget will be built to integrate with MDAnalysis, providing an interactive environment for controlling and rendering molecular objects directly within notebooks via Blender.

- Advanced Visualization Tools: We will develop tools for visualizing basic geometric features and even complex analysis results from MDAnalysis.

ggmolvis is bound by a [Code of Conduct](https://github.com/yuxuanzhuang/ggmolvis/blob/main/CODE_OF_CONDUCT.md).

### Installation

Follow the installation instructions in the [documentation](https://ggmolvis.readthedocs.io/en/latest/).

### Quickstart

```python
import MDAnalysis as mda
from MDAnalysis.tests.datafiles import PSF, DCD
from MDAnalysis.analysis.rms import RMSD

from ggmolvis.ggmolvis import GGMolVis

ggmv = GGMolVis()

u = mda.Universe(PSF, DCD)

# Trajectory visualization
residues_ag = u.select_atoms("resid 127 40")
residues_ag.visualize()

# Feature visualization
res_1 = residues_ag.residues[0].atoms
res_2 = residues_ag.residues[1].atoms
line = ggmv.distance(res_1, res_2, location=(5,0,0))
line.render()

# Analysis result visualization
rmsd = RMSD(u.select_atoms('name CA'))
rmsd.run()
vis = rmsd.visualize()
vis.render(mode='movie')

```

### Copyright

The ggmolvis source code is hosted at https://github.com/yuxuanzhuang/ggmolvis
and is available under the GNU General Public License, version 2 (see the file [LICENSE](https://github.com/yuxuanzhuang/ggmolvis/blob/main/LICENSE)).

Copyright (c) 2024, Yuxuan Zhuang


#### Acknowledgements
 
Project based on the 
[MDAnalysis Cookiecutter](https://github.com/MDAnalysis/cookiecutter-mda) version 0.1.
Please cite [MDAnalysis](https://github.com/MDAnalysis/mdanalysis#citation) when using ggmolvis in published work.
