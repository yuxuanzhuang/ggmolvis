from loguru import logger
from .trajectories import *

logger.warn(
    "molecules.py is deprecated and will be removed in a future version. "
    "Use trajectory.py instead.", stacklevel=2
)