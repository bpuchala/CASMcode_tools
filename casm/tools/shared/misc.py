import numpy as np


def pretty(x):
    y = x.copy()
    y[np.abs(y) < 1e-8] = 0.0
    return y
