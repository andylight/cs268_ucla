"""
Various utility functions relating to camera/geometry.
"""

def get_intrinsics(K):
    """ Return the camera intrinsic parameters from K.
    Input:
        nparray K
    Output:
        (float fx, float fy, principle_point)
    """
    return (K[0,0], K[1, 1], (K[0,2], K[1,2]))

def compute_x(line, y):
    """ ax + by + c = 0 => x = (-by - c) / a """
    return (-line[1]*y - line[2]) / line[0]
def compute_y(line, x):
    """ ax + by + c = 0 => y = (-ax - c) / b """
    return (-line[0]*x - line[2]) / line[1]
