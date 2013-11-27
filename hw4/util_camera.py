"""
Various utility functions relating to camera/geometry.
"""

import numpy as np, cv, cv2

from util import intrnd

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

def pt2homo(pt):
    """ Image point to homogeneous coords. """
    if len(pt) == 3:
        return pt
    if len(pt) != 2:
        raise TypeError("Must pass image coord to pt2homo: {0}".format(pt))
    return np.hstack((pt, np.array([1.0])))

def homo2pt(pt_h):
    """ Homogeneous coord to image (pixel) coord. """
    if len(pt_h) != 3:
        raise TypeError("Must pass 3-dim vector to homo2pt: {0}".format(pt_h))
    pt = pt_h[0:2]
    return pt / pt_h[2]

def draw_line(Irgb, line, color=(0, 255, 0)):
    """ Draws a line on an input image. If input image is not a three
    channeled image, then this will convert it.
    Input:
        nparray Irgb: (H x W x 1) or (H x W x 3)
        nparray line: [float a, float b, float c]
        tuple color: (int B, int G, int R)
    Output:
        nparray Irgb: (H x W x 3)
    """
    if len(Irgb.shape) != 3:
        Irgb = cv2.cvtColor(Irgb, cv.CV_GRAY2BGR)
        
    Irgb = Irgb.copy()
    h, w = Irgb.shape[0:2]
    pts = []
    for x in xrange(w):
        y = compute_line_y(line, x)
        if y > 0 and y < h:
            pts.append((x,y))
    cv.Line(cv.fromarray(Irgb), tuple(intrnd(*pts[0])), tuple(intrnd(*pts[-1])), color)
    return Irgb

def draw_points(Irgb, pts, color=(255, 0, 0)):
    """ Overlays points onto the image.
    Input:
        nparray Irgb: H x W x 3
        nparray pts: N x 2
            Rows of pixel coords (x,y)
    """
    for (x, y) in pts:
        Irgb[y, x, 0] = color[0]
        Irgb[y, x, 1] = color[1]
        Irgb[y, x, 2] = color[2]
    return Irgb

def find_line_segment(line, w, h):
    """ Computes good start/end points to display the line on the img
    with dimensions (w,h)
    NOTE: This function is just wrong and terrible. Use draw_line() instead.
    """
    ## Need at least 1 positive intercept.
    pt_xintercept = (compute_line_x(line, 0), 0)
    pt_yintercept = (0, compute_line_y(line, 0))
    is_xint_pos = (pt_xintercept[0] >= 0 and pt_xintercept[1] >= 0)
    is_yint_pos = (pt_yintercept[0] >= 0 and pt_yintercept[1] >= 0)
    if not any((is_xint_pos, is_yint_pos)):
        raise TypeError("Line does not occur in image")
    if is_xint_pos:
        ptA = pt_xintercept
    else:
        ptA = pt_yintercept
    def get_valid(pt1, pt2):
        pt1 = pt1 / pt1[2]
        pt2 = pt2 / pt2[2]
        if pt1[0] >= 0 and pt1[1] >= 0 and pt1[0] < w:
            return pt1
        elif pt2[0] >= 0 and pt2[0] >= 0 and pt2[1] < h:
            return pt2
        else:
            return None
    ptB = get_valid(np.cross(line, [0, 1, -(h-1)]),
                    np.cross(line, [1, 0, -(w-1)]))
    if ptB == None:
        raise Exception("WAT")
    return (ptA, tuple(ptB[0:2]))

def compute_line_x(line, y):
    """ Computes x coord of line at y:
        ax + by + c = 0
        x = (-by - c) / a
    """
    return (-line[1]*y - line[2]) / line[0]
def compute_line_y(line, x):
    """ Computes x coord of line at y:
        ax + by + c = 0
        y = (-ax - c) / b
    """
    return (-line[0]*x - line[2]) / line[1]

def make_crossprod_mat(v):
    """ Constructs the 3x3 skew-symmetric matrix representing the
    cross-product with v:
        v x w = v_hat * w
    Where 'x' denotes cross product, '*' denotes matrix mult, and
    v_hat is the output of: make_crossprod_mat(v)
    Input:
        nparray v: (2 x 1) or (3 x 1)
            If v is (2x1), then this will 'upgrade' it to homogenous
            coordinate prior to constructing v_hat.
    Output:
        nparray v_hat: (3x3)
    """
    if len(v) > 3:
        raise Exception("Must pass 2x1 or 3x1 vector to make_crossprod_mat \
(Received: {0})".format(v.shape))
    if len(v) == 2:
        v = np.hstack((v, np.array([1.0])))
    
    v_hat = np.array([[0, -v[2], v[1]],
                      [v[2], 0, -v[0]],
                      [-v[1], v[0], 0]], dtype=v.dtype)
    return v_hat

def normalize_det(A):
    """ Normalize A by a positive scalar factor c such that we get a
    matrix with determinant = +1.
    Recall: if A is nxn matrix: 
        We want: det(A*c) = 1 for some scaling factor c
        det(A*c) = (c^n)*det(A) = 1
        => c^n = (1 / det(A))
        => c = (1 / det(A)) ^ (1/n)
    Input:
        nparray A: N x N
    Output:
        nparray Anorm
            Anorm will have determinant +1.
    """
    det_A = np.linalg.det(A)
    if det_A > 0:
        c = np.power((1 / det_A), 1 / 3.0)
    else:
        c = -np.power((1 / -det_A), 1/ 3.0)
    return A * c
