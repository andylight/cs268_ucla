"""
A few utility functions.
"""
import numpy as np

def to_rgb(I, do_cpy=True):
    """ Convert input image into RGB (color). """
    if len(I.shape) == 3:
        if do_cpy:
            return I.copy()
        else:
            return I
    out = np.zeros([I.shape[0], I.shape[1], 3])
    out[:,:,0] = I
    out[:,:,1] = I
    out[:,:,2] = I
    return out

def intrnd(x):
    return int(round(x))
