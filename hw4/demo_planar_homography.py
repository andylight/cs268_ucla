import sys, os, pdb, argparse
import cv2, cv, numpy as np, scipy.misc

import calibrate_camera, util, util_camera
from util import intrnd
from util_camera import pt2homo, homo2pt

"""
USAGE:

    $ python demo_planar_homography.py [-h] [--show_epipolar]

A quick demo about planar scenes and homographies. The test images
are photos of a poster taken from a variety of viewpoints.

test_koopa() takes two image views of the poster, and finds the 
homography mapping points from one image to another. Displays the
epipolar lines if --show_epipolar is passed.

"""

IMGSDIR_CALIB_MED = 'calibrate_ek_med/'
IMGSDIR_CALIB_SMALL = 'calibrate_ek_small/'
IMGSDIR_KOOPA_MED = 'planar_koopa_med/'
IMGSDIR_KOOPA_SMALL = 'planar_koopa_small/'

def test_koopa(SHOW_EPIPOLAR=False):
    """ Estimate the homography H between two image views of the planar
    poster. Also decomposes H into {R, (1/d)T, N}.
    """
    imgpath1 = os.path.join(IMGSDIR_KOOPA_SMALL, 'DSCN0643.png')
    imgpath2 = os.path.join(IMGSDIR_KOOPA_SMALL, 'DSCN0648.png')
    img1 = cv2.imread(imgpath1)
    img2 = cv2.imread(imgpath2)
    pts1_ = ((150.0, 39.0),    # Upperleft corner redbox
             (220.0, 47.0),    # Koopa's left-eye
             (191.0, 55.0),    # Tip of Koopa's pencil
             (122.0, 128.0),   # Lowerleft corner of bluebox
             (144.0, 100.0),   # Tip of glue bottle
             (278.0, 129.0),   # Lowerright corner of greenbox
             )
    pts2_ = ((276.0, 34.0),    # Upperleft corner redbox
             (327.0, 61.0),    # Koopa's left-eye
             (286.0, 57.0),    # Tip of Koopa's pencil
             (141.0, 86.0),    # Lowerleft corner of bluebox
             (188.0, 75.0),    # Tip of glue bottle
             (237.0, 154.0),   # Lowerright corner of greenbox
             )
    calib_imgpaths = util.get_imgpaths(IMGSDIR_CALIB_SMALL)
    print "(Calibrating camera...)"
    K = calibrate_camera.calibrate_camera(calib_imgpaths, 9, 6, 0.023)
    print "Finished calibrating camera, K is:"
    print K
    print "(Estimating homography...)"
    pts1 = tup2nparray(pts1_)
    pts2 = tup2nparray(pts2_)
    pts1_norm = normalize_coords(pts1, K)
    pts2_norm = normalize_coords(pts2, K)
    
    # H goes from img1 -> img2
    H_ = estimate_planar_homography(pts1_norm, pts2_norm)
    print "Estimated homography, H is:",
    print H_
    print "rank(H):", np.linalg.matrix_rank(H_)
    if np.linalg.matrix_rank(H_) != 3:
        print "    Oh no! H is not full rank!"

    #### Normalize H := H / sigma_2(H_)
    U, S, V = np.linalg.svd(H_)
    H = H_ / S[1]

    #### Enforce positive depth constraint
    flag_flipped = False
    for i, pt1 in enumerate(pts1_norm):
        pt2 = pts2_norm[i]
        val = np.dot(np.hstack((pt2, [1.0])), np.dot(H, np.hstack((pt1, [1.0]))))
        if val < 0:
            if flag_flipped:
                print "WOAH, flipping twice?!"
                pdb.set_trace()
            print "FLIP IT! val_{0} was: {1}".format(i, val)
            H = H * -1
            flag_flipped = True
    # (Sanity check positive depth constraint)
    for i, pt1 in enumerate(pts1_norm):
        pt2 = pts2_norm[i]
        val = np.dot(np.hstack((pt2, [1.0])), np.dot(H, np.hstack((pt1, [1.0]))))
        if val < 0:
            print "WOAH, positive depth constraint violated!?"
            pdb.set_trace()
            
    #### Check projection error from I1 -> I2.
    errs = []
    for i, pt1 in enumerate(pts1_norm):
        pt1_h = np.hstack((pt1, np.array([1.0])))
        pt2_h = np.hstack((pts2_norm[i], np.array([1.0])))
        pt2_pred = np.dot(H, pt1_h)
        pt2_pred = pt2_pred / pt2_pred[2]
        errs.append(np.linalg.norm(pt2_h - pt2_pred))
    print "Projection error: {0} (Mean: {1} std: {2})".format(sum(errs), np.mean(errs), np.std(errs))

    #### Check if planar epipolar constraint is satisfied:
    ####     x2_hat * H * x1 = 0.0
    errs = []
    for i, pt1 in enumerate(pts1_norm):
        pt1_h = np.hstack((pt1, [1.0]))
        pt2_hat = util_camera.make_crossprod_mat(pts2_norm[i])
        errs.append(np.linalg.norm(np.dot(pt2_hat, np.dot(H, pt1_h))))
    print "Epipolar constraint error: {0} (Mean: {1} std: {2})".format(sum(errs), np.mean(errs), np.std(errs))

    #### Draw epipolar lines
    Irgb1 = cv2.imread(imgpath1, cv2.CV_LOAD_IMAGE_COLOR)
    Irgb2 = cv2.imread(imgpath2, cv2.CV_LOAD_IMAGE_COLOR)
    if SHOW_EPIPOLAR:
        draw_epipolar_lines(Irgb1, Irgb2, pts1, pts2, H)

    decomps = decompose_H(H)
    print
    for i, (R, Ts, N) in enumerate(decomps):
        print "==== Decomposition {0} ====".format(i)
        print "    R is:", R
        print "    Ts is:", Ts
        print "    N is:", N
        H_redone = (R + np.dot(np.array([Ts]).T, np.array([N])))
        print "norm((R+Ts*N.T) - H, 'fro'):", np.linalg.norm(H - H_redone, 'fro')
        print "    det(R)={0} rank(R)={1}".format(np.linalg.det(R), np.linalg.matrix_rank(R))
        #### Sanity check that H_redone still maps I1 to I2
        errs = []
        for ii, pt1 in enumerate(pts1_norm):
            pt1_h = np.hstack(([pt1, [1.0]]))
            pt2_h = np.hstack((pts2_norm[ii], [1.0]))
            pt1_proj = np.dot(H_redone, pt1_h)
            pt1_proj /= pt1_proj[2]
            print pt1_proj
            print pt2_h
            print
            errs.append(np.linalg.norm(pt1_proj - pt2_h))
        print "Reprojection error: {0} (mean={1}, std={2})".format(sum(errs), np.mean(errs), np.std(errs))

def estimate_planar_homography(pts1, pts2):
    """ Estimates the planar homography between two images of a planar
    scene, given corresponding points.
    Input:
        nparray pts1, pts2: N x 2
        nparray K
    output:
        nparray H
    """
    pts1 = pts1.astype('float32')
    pts2 = pts2.astype('float32')
    retval, mask = cv2.findHomography(pts1, pts2)
    return retval

def tup2nparray(pts):
    """ Converts a tuple/list of points to a nparray of points, of
    dimension Nx2/Nx3.
    """
    out = np.zeros([len(pts), len(pts[0])])
    for i, pt in enumerate(pts):
        out[i, :] = pt
    return out

def get_worldplane_coords():
    # See: planar_koopa_med/reference_distances.png
    # Origin is at upper-left corner of red box, with X axis along
    # width of paper, and Y axis down height of paper. In other words,
    # X axis goes from one side of the red box to the other.
    # Y axis goes down from red box to the blue box.
    # Below coords are in centimeters (cm), but we will output in meters
    worldpts = np.array([[0.0, 0.0],       # A (upper-left red corner)
                         [17.5, 0.0],      # B (upper-right red corner)
                         [17.5, 13.4],     # C (lower-right red corner)
                         [9.0, 21.8],      # D (lower-left green corner)
                         [8.5, 21.8],      # E (lower-right blue corner)
                         [0.0, 21.8],      # F (lower-left blue corner)
                         [0.0, 14.3],      # G (upper-left blue corner)
                         [0.0, 13.4],      # H (lower-left red corner)
                         [6.7, 5.0],       # K (tip of pen)
                         [9.0, 14.3],      # L (upper-left of green corner)
                         [1.8, 16.2],      # M (tip of glue bottle)
                         ])
    return worldpts * 1e-2 # convert from cm -> meters
    
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
        U1[:, 2] = np.dot(util_camera.make_crossprod_mat(v2), u1)
        return U1
    def make_U2():
        U2 = np.zeros((3, 3))
        U2[:, 0] = v2
        U2[:, 1] = u2
        U2[:, 2] = np.dot(util_camera.make_crossprod_mat(v2), u2)
        return U2
    def make_W1():
        W1 = np.zeros((3,3))
        W1[:, 0] = np.dot(H, v2)
        W1[:, 1] = np.dot(H, u1)
        W1[:, 2] = np.dot(util_camera.make_crossprod_mat(np.dot(H, v2)), np.dot(H, u1))
        return W1
    def make_W2():
        W2 = np.zeros((3, 3))
        W2[:, 0] = np.dot(H, v2)
        W2[:, 1] = np.dot(H, u2)
        W2[:, 2] = np.dot(util_camera.make_crossprod_mat(np.dot(H, v2)), np.dot(H, u2))
        return W2
    U1 = make_U1()
    U2 = make_U2()
    W1 = make_W1()
    W2 = make_W2()
    
    # Generate 4 possible solutions
    R1 = np.dot(W1, U1.T)
    N1 = np.dot(util_camera.make_crossprod_mat(v2), u1)
    Ts1 = np.dot((H - R1), N1)
    
    R2 = np.dot(W2, U2.T)
    N2 = np.dot(util_camera.make_crossprod_mat(v2), u2)
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
        R3 = util_camera.normalize_det(R3)
        decomps.append((R3, Ts3, N3))
    else:
        R1 = util_camera.normalize_det(R1)
        decomps.append((R1, Ts1, N1))
    if N2[2] < 0:
        # N2 is impossible, N4 must be correct
        R4 = util_camera.normalize_det(R4)
        decomps.append((R4, Ts4, N4))
    else:
        R2 = util_camera.normalize_det(R2)
        decomps.append((R2, Ts2, N2))
    return decomps
    
def draw_epipolar_lines(Irgb1, Irgb2, pts1, pts2, H):
    """ Draws epipolar lines, and displays them to the user in an
    interactive manner.
    Input:
        nparray Irgb1, Irgb2: H x W x 3
        nparray pts1, pts2: N x 2
        nparray H: 3x3
    """
    h, w = Irgb1.shape[0:2]
    for i, pt1 in enumerate(pts1):
        Irgb1_ = Irgb1.copy()
        Irgb2_ = Irgb2.copy()
        pt2 = pts2[i]
        pt1_h = np.hstack((pt1, np.array([1.0])))
        pt2_h = np.hstack((pt2, np.array([1.0])))
        
        epiline2 = np.dot(util_camera.make_crossprod_mat(pt2_h), np.dot(H, pt1_h))
        epiline1 = np.dot(H.T, epiline2)
        epiline1 = epiline1 / epiline1[2]
        epiline2 = epiline2 / epiline2[2]
        print "Epiline1 is: slope={0} y-int={1}".format(-epiline1[0] / epiline1[1],
                                                         -epiline1[2] / epiline1[1])
        print "Epiline2 is: slope={0} y-int={1}".format(-epiline2[0] / epiline2[1],
                                                         -epiline2[2] / epiline2[1])
        Irgb1_ = util_camera.draw_line(Irgb1_, epiline1)
        Irgb2_ = util_camera.draw_line(Irgb2_, epiline2)
        cv.Circle(cv.fromarray(Irgb1_), tuple(intrnd(*pts1[i])), 3, (255, 0, 0))
        cv.Circle(cv.fromarray(Irgb2_), tuple(intrnd(*pts2[i])), 3, (255, 0, 0))
        print "(Displaying epipolar lines from img1 to img2. Press <enter> to continue.)"
        cv2.namedWindow('display1')
        cv2.imshow('display1', Irgb1_)
        cv2.namedWindow('display2')
        cv2.imshow('display2', Irgb2_)
        cv2.waitKey(0)
    
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

def parse_args():
    DESCRIPTION = """This is a demo about homographies relating a view(s)
to planar scenes."""
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument("--show_epipolar", action='store_true',
                        help="Interactively display epipolar lines.",
                        default=False)
    return parser.parse_args()

def main():
    args = parse_args()
    print "======== (1) test_kooopa: two view homography ========"
    test_koopa(SHOW_EPIPOLAR=args.show_epipolar)

if __name__ == '__main__':
    main()
