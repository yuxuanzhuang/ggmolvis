from ggmolvis.ggmolvis import GGMolVis


import MDAnalysis as mda
from MDAnalysisTests.datafiles import GRO
import bpy
from PIL import Image
import numpy as np
import pytest


@pytest.mark.xslow
def test_render_basic_bg(tmpdir):
    # confirm correct propagation of background color
    # and image size specified
    ggmv = GGMolVis()
    u = mda.Universe(GRO)
    mol = ggmv.molecule(u.atoms)
    with tmpdir.as_cwd():
        ggmv.render(object=mol,
                    resolution=(50, 50),
                    filepath=f"test_red.png",
                    mode="image",
                    # red background
                    composite_bg_rgba=(0.9, 0, 0, 1.0),
                    lens=35)
        with Image.open("test_red.png") as img:
            img_array = np.array(img)
            red = img_array[..., 0]
            green = img_array[..., 1]
            blue = img_array[..., 2]
    # red should be the dominant channel
    assert red.sum() > (green.sum() + blue.sum())
    # image size should be 50 x 50 x 4 channels
    assert img_array.size == 10_000
