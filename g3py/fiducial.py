import numpy as np


def sign(p1, p2, p3):
    return (p1[0] - p3[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p3[1])


def check_within_fiducial(CoreX, CoreY):
    x = np.array(CoreX)
    y = np.array(CoreY)
    v0, v1, v2, v3, v4, v5 = [
        [-91.5, 16],
        [-64.0, 64.5],
        [40.8, 64.5],
        [64.0, 23.5],
        [16.0, -59.0],
        [-51.0, -59.0],
    ]

    b1 = sign((CoreX, CoreY), v0, v1) < 0.0
    b2 = sign((CoreX, CoreY), v1, v2) < 0.0
    b3 = sign((CoreX, CoreY), v2, v3) < 0.0
    b4 = sign((CoreX, CoreY), v3, v4) < 0.0
    b5 = sign((CoreX, CoreY), v4, v5) < 0.0
    b6 = sign((CoreX, CoreY), v5, v0) < 0.0

    return np.all([b1, b2, b3, b4, b5, b6], axis=0)
