import pdb
import util_camera, util
import numpy as np, numpy.linalg, cv2

from util import intrnd
from util_camera import compute_x, compute_y

def estimate_planar_homography(I, line1, line2, K, win1, win2, lane_width):
    """ Estimates the planar homography H between the camera image
    plane, and the World (ground) plane.
    The World reference frame is directly below the camera, with:
        - Z-axis parallel to the road
        - X-axis on the road plane
        - Y-axis pointing upward from the road plane
        - origin is on the road plane (Y=0), and halfway between the
          two lanes boundaries (center of road).
    Input:
        nparray I
        nparray line1, line2: [a, b, c]
        nparray K
            The 3x3 camera intrinsic matrix.
        tuple win1, win2: (float x, float y, float w, float h)
        float lane_width
    Output:
        nparray H
    where H is a 3x3 homography.
    """
    h, w = I.shape[0:2]
    x_win1 = intrnd(w * win1[0])
    y_win1 = intrnd(h * win1[1])
    x_win2 = intrnd(w * win2[0])
    y_win2 = intrnd(h * win2[1])
    
    pts = []
    NUM = 10
    h_win = intrnd(h*win1[3]) # Assume window heights same btwn left/right
    for i in xrange(NUM):
        frac = i / float(NUM)
        y_cur = intrnd((y_win1-(h_win/2) + frac*h_win))
        pt_i = (compute_x(line1, h - y_cur), y_cur)
        pt_j = (compute_x(line2, h - y_cur), y_cur)
        pts.append((pt_i, pt_j))
        
    #y_low = intrnd(y_win1 + 0.25*(h*win1[3]))
    #y_high = intrnd(y_win1 - 0.25*(h*win1[3]))
    ## First pair of points across from the lane
    #pt1 = (compute_x(line1, h - y_low), y_low)
    #pt2 = (compute_x(line2, h - y_low), y_low)
    ## Second pair of points across from the lane
    #pt3 = (compute_x(line1, h - y_high), y_high)
    #pt4 = (compute_x(line2, h - y_high), y_high)
    #pts = ((pt1, pt2), (pt3, pt4))
    r1 = solve_for_r1(pts,
                      K,
                      lane_width)
    vanishing_pt = np.cross(line1, line2)
    vanishing_pt = vanishing_pt / vanishing_pt[2]
    vanishing_pt[1] = h - vanishing_pt[1]
    ## DEBUG Plot points on image, save to file for visual verification
    Irgb = util.to_rgb(I)
    COLOURS = [(255, 0, 0), (0, 255, 0)]
    for i, (pt_i, pt_j) in enumerate(pts):
        clr = COLOURS[i % 2]
        cv2.circle(Irgb, tuple(map(intrnd, pt_i)), 5, clr)
        cv2.circle(Irgb, tuple(map(intrnd, pt_j)), 5, clr)
    cv2.circle(Irgb, (intrnd(vanishing_pt[0]), intrnd(vanishing_pt[1])), 5, (0, 0, 255))
    cv2.imwrite("_Irgb.png", Irgb)

    r3 = solve_for_r3(vanishing_pt, line1, line2, K)
    T = solve_for_t(pts, K, r1, r3, lane_width)
    
    H = np.zeros([3,3])
    H[:, 0] = r1
    H[:, 1] = r3
    H[:, 2] = T
    return np.dot(K, H)

def solve_for_r1(pts, K, lane_width):
    """ Solve for first column of the rotation matrix, utilizing the
    fact that we know the lane width. We require two pairs of points
    (p1,p2), (p3,p4), such that p1 is directly across from p2 (same
    for p3, p4).
    Input:
        tuple pts: ((p1, p2), (p3, p4))
            where each point is a pixel coord: (float x, float y)
        nparray K
            The 3x3 camera intrinsic matrix.
        float lane_width
            The width of the lane (e.g., 3.66 meters).
    Output:
        nparray r1
    A 3x1 column vector consisting of the first column of R.
    """
    (fx, fy, (cx, cy)) = util_camera.get_intrinsics(K)
    # Construct data matrix A
    A = np.zeros([len(pts) * 2, 4])
    i = 0 # Actual index into A
    for ii, (pi, pj) in enumerate(pts):
        xi, yi = pi
        xj, yj = pj
        A[i, :]   = (0, -fy, (-cy + yj), (yi - yj))
        A[i+1, :] = (fx, 0, (cx - xj), (-xj + xi))
        #A[i+2, :] = (-yj*fx, xj*fy, -yj*cx + xj*cy, yj*xi - xj*yi)
        i += 2
    U, S, V = numpy.linalg.svd(A)
    v = V[-1, :]
    residual = numpy.linalg.norm(np.dot(A, v.T))
    print "== compute_r1"
    print "Residual: {0}".format(residual)
    print "    Rank(A):", np.linalg.matrix_rank(A)
    gamma = v[-1]
    v_norm = v / gamma
    r11, r21, r31, _ = v_norm
    return np.array([r11, r21, r31])

def solve_for_r3(vanishing_pt, line1, line2, K):
    """ Solve for the third column r3 of the rotation matrix,
    utilizing the vanishing point of the lanes.
    Input:
        nparray vanishing_pt: [x, y, 1]
            Pixel image location of the vanishing point defined by the
            lanes.
        nparray line1, line2: [a, b, c]
            Lines vectors of the left/right lanes, in the form:
                [a, b, c]
            such that:
                ax + by + c = 0
        nparray K
            Camera intrinsic matrix.
    Output:
        nparray r3
            The third column of the rotation matrix R.
    """
    Kinv = numpy.linalg.inv(K)
    r3 = np.dot(Kinv, vanishing_pt)
    r3_norm = r3 / r3[2]
    return r3_norm

def solve_for_t(pts, K, r1, r3, lane_width):
    """ Recover the translation vector T, using the computed r1, r3.
    The input points pairs must be directly across from the lanes.
    Input:
        tuple pts: ((pt1, pt2), ...)
            Each point pair (pt_i, pt_j) must be directly across the
            lanes i.e. (X_j - X_i) = 3.66 meters, and:
                X_i = -1.83 meters
                X_j = +1.83 meters
        nparray K
            3x3 camera intrinsic matrix.
        nparray r1, r3
            The 3x1 column vectors comprising the first/third columns
            of the rotation matrix R.
        float lane_width
            Width of the lane (in meters).
    Output:
        nparray T
            The translation vector T as a 3x1 column vector.
    """
    (fx, fy, (cx, cy)) = util_camera.get_intrinsics(K)
    Kinv = numpy.linalg.inv(K)
    r11, r21, r31 = r1
    r13, r23, r33 = r3
    ww = lane_width / 2
    # Construct data matrix A
    A = np.zeros([len(pts) * 6, 5])
    i = 0 # Actual index into A
    for ii, (pi, pj) in enumerate(pts):
        xi, yi = pi
        xj, yj = pj
        bi = np.array([(xi - cx) / fx,
                       (yi - cy) / fy,
                       1]).T
        bj = np.array([(xj - cx) / fx,
                       (yj - cy) / fy,
                       1]).T
        A[i, :]   = [r13, 1, 0, 0, -ww*r11 - bi[0]]
        A[i+1, :] = [r23, 0, 1, 0, -ww*r21 - bi[1]]
        A[i+2, :] = [r33, 0, 0, 1, -ww*r31 - bi[2]]

        A[i+3, :] = [r13, 1, 0, 0, -ww*r11 - bj[0]]
        A[i+4, :] = [r23, 0, 1, 0, -ww*r21 - bj[1]]
        A[i+5, :] = [r33, 0, 0, 1, -ww*r31 - bj[2]]
        i += 6
    print "== compute_t"
    print "    Rank(A):", np.linalg.matrix_rank(A)
    U, S, V = numpy.linalg.svd(A)
    v = V[-1, :]
    gamma = v[-1]
    v_norm = v / gamma
    Z, tx, ty, tz, _ = v_norm
    return np.array([tx, ty, tz]).T

def main():
    # K matrix given by the Caltech Lanes dataset (CameraInfo.txt)
    K = np.array([[309.4362,     0,        317.9034],
                  [0,         344.2161,    256.5352],
                  [0,            0,            1   ]])
    line1 = np.array([  -1.35431663,    1.,          124.12336564])
    line2 = np.array([   1.23775052,    1.,         -695.71784631])
    win1 = (0.4, 0.60, 0.2, 0.25)
    win2 = (0.62, 0.60, 0.2, 0.25)
    lane_width = 3.66 # 3.66 meters
    imgpath = 'imgs_sample/f00001.png'
    I = cv2.imread(imgpath, cv2.CV_LOAD_IMAGE_GRAYSCALE)
    H = estimate_planar_homography(I, line1, line2, K, win1, win2, lane_width)
    print H
    print "    rank(H): {0}".format(np.linalg.matrix_rank(H))
    print "The following should be identity (inv(H) * H):"
    print np.dot(numpy.linalg.inv(H), H)

    print "(Evaluating a few world points to see where they lie on the image)"
    for pt in [(-1.83, 2, 1), (-1.83, 5, 1), (0, 0, 1), (-12347, 1281, 1)]:
        pt_np = np.array(pt)
        pt_img = np.dot(H, pt_np)
        pt_img = pt_img / pt_img[2]
        print pt_img
    
    print "Done."

if __name__ == '__main__':
    main()
