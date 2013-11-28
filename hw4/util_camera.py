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

def normalize_coords(pts, K):
    """ Normalize pixel coordinates into image coordinates with
    intrinsic matrix K.
    Input:
        nparray pts: N x 2
    Output:
        nparray pts_norm: N x 2
    """
    out = np.zeros(pts.shape, dtype=pts.dtype)
    Kinv = np.linalg.inv(K)
    for i, pt in enumerate(pts):
        pt_h = np.hstack((pt, 1.0))
        pt_norm_h = np.dot(Kinv, pt_h)
        pt_norm = pt_norm_h[0:2]
        out[i, :] = pt_norm
    return out

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

def decompose_H(H, pts1=None, pts2=None):
    """ Decomposes homography H into (R, (1/d)*T, N), where (R,T) is
    the relative camera motion between the planes, and N is the normal
    vector of the plane w.r.t. the first camera view.
    Input:
        nparray H: 3x3
        nparray pts1, pts2: Nx2
            If given, then this function will only output the
            physically possible (R,T,N) by enforcing the positive
            depth constraint.
    Output:
        (nparray R, nparray Ts, nparray N)
    Where R is 3x3, N is a 3x1 column vector, and Ts is a 3x1 column
    vector defined UP TO an unknown scale (1/d).
    """
    U, S, Vt = np.linalg.svd(np.dot(H.T, H))
    if np.linalg.det(U) < 0:
        # We require that U, Vt have det=+1
        U = -U
        Vt = -Vt
    V = Vt.T    # Recall: svd outputs V as V.T
    v1 = V[:, 0]
    v2 = V[:, 1]
    v3 = V[:, 2]
    norm_ = np.sqrt(S[0]**2.0 - S[2]**2.0)
    u1 = ((np.sqrt(1 - S[2]**2.0)*v1) + (np.sqrt(S[0]**2.0 - 1)*v3)) / norm_
    u2 = ((np.sqrt(1 - S[2]**2.0)*v1) - (np.sqrt(S[0]**2.0 - 1)*v3)) / norm_
    def make_U1():
        U1 = np.zeros((3, 3))
        U1[:, 0] = v2
        U1[:, 1] = u1
        U1[:, 2] = np.dot(make_crossprod_mat(v2), u1)
        return U1
    def make_U2():
        U2 = np.zeros((3, 3))
        U2[:, 0] = v2
        U2[:, 1] = u2
        U2[:, 2] = np.dot(make_crossprod_mat(v2), u2)
        return U2
    def make_W1():
        W1 = np.zeros((3,3))
        W1[:, 0] = np.dot(H, v2)
        W1[:, 1] = np.dot(H, u1)
        W1[:, 2] = np.dot(make_crossprod_mat(np.dot(H, v2)), np.dot(H, u1))
        return W1
    def make_W2():
        W2 = np.zeros((3, 3))
        W2[:, 0] = np.dot(H, v2)
        W2[:, 1] = np.dot(H, u2)
        W2[:, 2] = np.dot(make_crossprod_mat(np.dot(H, v2)), np.dot(H, u2))
        return W2
    U1 = make_U1()
    U2 = make_U2()
    W1 = make_W1()
    W2 = make_W2()
    
    # Generate 4 possible solutions
    R1 = np.dot(W1, U1.T)
    N1 = np.dot(make_crossprod_mat(v2), u1)
    Ts1 = np.dot((H - R1), N1)
    
    R2 = np.dot(W2, U2.T)
    N2 = np.dot(make_crossprod_mat(v2), u2)
    Ts2 = np.dot((H - R2), N2)
    
    R3 = R1
    N3 = -N1
    Ts3 = -Ts1
    
    R4 = R2
    N4 = -N2
    Ts4 = -Ts2
    ## Remove physically impossible decomps: n3 < 0
    decomps = []
    if N1[2] < 0:
        # N1 is impossible, N3 must be correct
        R3 = normalize_det(R3)
        decomps.append((R3, Ts3, N3))
    else:
        R1 = normalize_det(R1)
        decomps.append((R1, Ts1, N1))
    if N2[2] < 0:
        # N2 is impossible, N4 must be correct
        R4 = normalize_det(R4)
        decomps.append((R4, Ts4, N4))
    else:
        R2 = normalize_det(R2)
        decomps.append((R2, Ts2, N2))
    return decomps
