from ggmolvis.utils import lerp
import numpy as np


def test_lerp():
    assert lerp(1, 2, 0.5) == 1.5
    assert lerp(1, 2, 0.3) == 1.3
    assert lerp(3, 7, 0.2) == 3.8
    assert np.allclose(lerp([1, 2, 3], [4, 5, 6], 0.5), [2.5, 3.5, 4.5])
    assert np.allclose(
        lerp(
            np.array([[1, 2, 3], [4, 5, 6]]), np.array([[7, 8, 9], [10, 11, 12]]), 0.5
        ),
        np.array(
            [
                [4.0, 5.0, 6.0],
                [7.0, 8.0, 9.0],
            ]
        ),
    )
