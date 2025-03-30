import bpy
import molecularnodes as mn
import MDAnalysis as mda
import numpy as np
from pydantic import BaseModel, Field, validator, ValidationError
from typing import Tuple, List, Union

from .base import SceneObject

class Light(SceneObject):
    """Class for the light."""