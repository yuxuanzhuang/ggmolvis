from importlib import resources

_base_ref = resources.files('ggmolvis.tests.data')

PSF = _base_ref / 'adk.psf'
DCD = _base_ref / 'adk_dims.dcd'