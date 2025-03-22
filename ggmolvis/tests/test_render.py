from ggmolvis.ggmolvis import GGMolVis
from ggmolvis.tests.data import PSF, DCD

from MDAnalysis import Universe
from ggmolvis.utils import MOL_AVAILABLE_STYLES, AVAILABLE_MATERIALS

import pytest
from pytest import fixture, mark
from .utils import cleanup

import bpy
from loguru import logger

@pytest.fixture
def render_png(request, tmp_path_factory):
    module_tmp_path = tmp_path_factory.mktemp("render")
    return str(module_tmp_path / f"{request.node.name}.png")

@pytest.fixture
def render_mp4(request, tmp_path_factory):
    module_tmp_path = tmp_path_factory.mktemp("render")
    logger.info(f"MP4: {module_tmp_path}")
    return str(module_tmp_path / f"{request.node.name}.mp4")

@fixture(scope='module')
def universe():
    return Universe(PSF, DCD)

@fixture(scope='module')
def atomgroup(universe):
    return universe.select_atoms('resid 40 127')

@cleanup
def test_render_scene(ggmv, atomgroup, render_png):
    mol = ggmv.trajectory(atomgroup)
    logger.info(f"PNG: {render_png}")
    ggmv.render(filepath=render_png, resolution=(300, 200))

# add test for rendering
render_image_options = [
    {'frame': 0},
    {'frame': 60},
]
@cleanup
@mark.parametrize('render_image_options', render_image_options)
def test_render_png(ggmv, atomgroup, render_image_options, render_png):
    mol = ggmv.trajectory(atomgroup)
    logger.info(f"PNG: {render_png}")
    ggmv.render(mol, **render_image_options,
                    resolution=(300, 200),
                    filepath=render_png)

render_mp4_options = [
    {'frame_range': (0, 60, 10)},
    {'frame_range': (0, 60, 10), 'track': True},
]
@cleanup
@mark.parametrize('render_mp4_options', render_mp4_options)
def test_render_mp4(ggmv, atomgroup, render_mp4_options, render_mp4):
    mol = ggmv.trajectory(atomgroup)
    logger.info(f"MP4: {render_mp4}")
    ggmv.render(mol, **render_mp4_options,
                resolution=(300, 200),
                filepath=render_mp4)


@cleanup
def test_render_error_frame_range(ggmv, atomgroup, render_mp4):
    mol = ggmv.trajectory(atomgroup)
    with pytest.raises(ValueError):
        ggmv.render(mol, frame=0, frame_range=(0, 4, 1),
                        resolution=(300, 200),
                        filepath=render_mp4)
    
    with pytest.raises(ValueError, match='frame_range must be'):
        ggmv.render(mol, frame_range=(0, 4),
                        resolution=(300, 200),
                        filepath=render_mp4)
