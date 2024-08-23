from typing import Tuple

import MDAnalysis as mda
from MDAnalysis.coordinates.memory import MemoryReader
import numpy as np


def make_Universe(
    extras: Tuple[str] = tuple(),
    size: Tuple[int, int, int] = (125, 25, 5),
    n_frames: int = 0,
    velocities: bool = False,
    forces: bool = False
) -> mda.Universe:
    """Make a dummy reference Universe

    Allows the construction of arbitrary-sized Universes. Suitable for
    the generation of structures for output.

    Preferable for testing core components because:
      * minimises dependencies within the package
      * very fast compared to a "real" Universe

    Parameters
    ----------
    extras : tuple of strings, optional
      extra attributes to add to Universe:
      u = make_Universe(('masses', 'charges'))
      Creates a lightweight Universe with only masses and charges.
    size : tuple of int, optional
      number of elements of the Universe (n_atoms, n_residues, n_segments)
    n_frames : int
      If positive, create a fake Reader object attached to Universe
    velocities : bool, optional
      if the fake Reader provides velocities
    force : bool, optional
      if the fake Reader provides forces

    Returns
    -------
    MDAnalysis.core.universe.Universe object

    """

    n_atoms, n_residues, n_segments = size
    trajectory = n_frames > 0
    u = mda.Universe.empty(
        # topology things
        n_atoms=n_atoms,
        n_residues=n_residues,
        n_segments=n_segments,
        atom_resindex=np.repeat(
            np.arange(n_residues), n_atoms // n_residues),
        residue_segindex=np.repeat(
            np.arange(n_segments), n_residues // n_segments),
        # trajectory things
        trajectory=trajectory,
        velocities=velocities,
        forces=forces,
    )
    if extras is None:
        extras = []
    for ex in extras:
        u.add_TopologyAttr(ex)

    if trajectory:
        pos = np.arange(3 * n_atoms * n_frames).reshape(n_frames, n_atoms, 3)
        vel = pos + 100 if velocities else None
        fcs = pos + 10000 if forces else None
        reader = MemoryReader(
            pos,
            velocities=vel,
            forces=fcs,
        )
        u.trajectory = reader

    return u
