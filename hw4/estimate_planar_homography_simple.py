import sys, os, pdb, argparse
import cv2, cv, numpy as np, scipy.misc

import calibrate_camera, util, util_camera
from util import intrnd
from util_camera import pt2homo, homo2pt

"""
USAGE:

    $ python estimate_planar_homography_simple.py [-h] [--show_epipolar]

Two quick demos about planar scenes and homographies. The test images
are photos of a poster taken from a variety of viewpoints.

test_koopa() takes two image views of the poster, and finds the 
homography mapping points from one image to another.

test_koopa_singleimg() takes one image view, and finds the homography
mapping points from the image to points on the world frame (i.e. poster
plane). A perspective correction is then performed and displayed to the
user.
"""

IMGSDIR_CALIB_MED = 'calibrate_ek_med/'
IMGSDIR_CALIB_SMALL = 'calibrate_ek_small/'
IMGSDIR_KOOPA_MED = 'planar_koopa_med/'
IMGSDIR_KOOPA_SMALL = 'planar_koopa_small/'

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
    #return cv2.getPerspectiveTransform(pts1[0:4,:], pts2[0:4,:])
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

def test_koopa_singleimage(SHOW_EPIPOLAR=False):
    """ Here, I want to do "undo" the perspective distortion in the
    image, based on the fact that I know:
        - The scene is planar
        - The location of four points in the image that create a 
          rectangle in the world plane with known metric lengths
          Note: I don't have to necessarily know the metric lengths,
                but I should know the ratios to construct a corrected
                image with lengths preserved.
    In the literature, this is also known as: perspective rectification.
    """
    imgpath = os.path.join(IMGSDIR_KOOPA_MED, 'DSCN0643.png')
    pts1 = np.array([[373.0, 97.0],        # A (upper-left red corner)
                     [650.0, 95.0],        # B (upper-right red corner)
                     [674.0, 215.0],       # C (lower-right red corner)
                     [503.0, 322.0],       # D (lower-left green corner)
                     [494.0, 322.0],       # E (lower-right blue corner)
                     [303.0, 321.0],       # F (lower-left blue corner)
                     [332.0, 225.0],       # G (upper-left blue corner)
                     [335.0, 215.0],       # H (lower-left red corner)
                     [475.0, 135.0],       # K (tip of pen)
                     [507.0, 225.0],       # L (upper-left green corner)
                     [362.0, 248.0],       # M (tip of glue bottle)
                     ])
    worldpts = get_worldplane_coords()
    assert len(pts1) == len(worldpts)
    calib_imgpaths = util.get_imgpaths(IMGSDIR_CALIB_MED)
    print "(Calibrating camera...)"
    
    if True:
        print "(Using pre-computed camera matrix K!)"
        K = np.array([[ 158.23796519,    0.0,          482.07814366],
                      [   0.,           28.53758493,  333.32239125],
                      [   0.,            0.,            1.        ]])
    else:
        K = calibrate_camera.calibrate_camera(calib_imgpaths, 9, 6, 0.023)
    print "Finished calibrating camera, K is:"
    print K
    
    pts1_norm = normalize_coords(pts1, K)
    HL = estimate_planar_homography(pts1_norm, worldpts)

    print "Estimated homography, H is:"
    print HL
    print "rank(H):", np.linalg.matrix_rank(HL)
    if np.linalg.matrix_rank(HL) != 3:
        print "    Oh no! H is not full rank!"

    #### Normalize H := H / sigma_2(H_)
    U, S, V = np.linalg.svd(HL)
    H = HL / S[1]

    #### Enforce positive depth constraint
    flag_flipped = False
    for i, pt1 in enumerate(pts1_norm):
        pt2 = worldpts[i]
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
        pt2 = worldpts[i]
        val = np.dot(np.hstack((pt2, [1.0])), np.dot(H, np.hstack((pt1, [1.0]))))
        if val < 0:
            print "WOAH, positive depth constraint violated!?"
            pdb.set_trace()
    #### Check projection error from I1 -> I2.
    errs = []
    for i, pt1 in enumerate(pts1_norm):
        pt1_h = np.hstack((pt1, np.array([1.0])))
        pt2_h = np.hstack((worldpts[i], np.array([1.0])))
        pt2_pred = np.dot(H, pt1_h)
        pt2_pred = pt2_pred / pt2_pred[2]
        errs.append(np.linalg.norm(pt2_h - pt2_pred))
    print "Projection error: {0} (Mean: {1} std: {2})".format(sum(errs), np.mean(errs), np.std(errs))

    #### Check if planar epipolar constraint is satisfied:
    ####     x2_hat * H * x1 = 0.0
    errs = []
    for i, pt1 in enumerate(pts1_norm):
        pt1_h = np.hstack((pt1, [1.0]))
        pt2_hat = util_camera.make_crossprod_mat(worldpts[i])
        errs.append(np.linalg.norm(np.dot(pt2_hat, np.dot(H, pt1_h))))
    print "Epipolar constraint error: {0} (Mean: {1} std: {2})".format(sum(errs), np.mean(errs), np.std(errs))

    decomps = decompose_H(H)
    for i, (R, Ts, N) in enumerate(decomps):
        print "\n==== Decomposition {0} ====".format(i)
        print "    R is:"
        print R
        print "    Ts is:", Ts
        print "    N is:", N
        H_redone = (R + np.dot(np.array([Ts]).T, np.array([N])))
        print "norm((R+Ts*N.T) - H, 'fro'):", np.linalg.norm(H - H_redone, 'fro')
        print "    det(R)={0} rank(R)={1}".format(np.linalg.det(R), np.linalg.matrix_rank(R))
        #### Sanity check that H_redone still maps I1 to I2
        errs = []
        for ii, pt1 in enumerate(pts1_norm):
            pt1_h = np.hstack(([pt1, [1.0]]))
            pt2_h = np.hstack((worldpts[ii], [1.0]))
            pt1_proj = np.dot(H_redone, pt1_h)
            pt1_proj /= pt1_proj[2]
            #print pt1_proj
            #print pt2_h
            #print
            errs.append(np.linalg.norm(pt1_proj - pt2_h))
        print "Reprojection error: {0} (mean={1}, std={2})".format(sum(errs), np.mean(errs), np.std(errs))

    #### Sanity-check that world points project to pixel points
    Hinv = np.linalg.inv(H)
    errs = []
    for i, worldpt in enumerate(worldpts):
        pt_img = np.hstack((pts1[i], [1.0]))
        p = np.dot(Hinv, np.hstack((worldpt, [1.0])))
        p /= p[2]
        p = np.dot(K, p)
        errs.append(np.linalg.norm(pt_img - p))
    print "world2img errors (in pixels): {0} (mean={1} std={2})".format(sum(errs), np.mean(errs), np.std(errs))
        
    #### Perform Inverse Perspective Mapping (undo perspective effects)
    p_ul = (372, 97)    # Upperleft of redbox (x,y)
    p_ur = (651, 95)    # Upperright of redbox (x,y)
    p_ll = (303, 321)   # Lowerleft of bluebox (x,y)
    p_lr = (695, 324)   # Lowerright of greenbox (x,y)
    Hinv = np.linalg.inv(H)
    def compute_plane2I(pp_ul, pp_ur, pp_ll, pp_lr):
        """ Computes H_RI, i.e. road->orthogonal-image. Input points 
        are image pixel coords.
        """
        pts_img = np.array([pp_ul, pp_ur, pp_ll, pp_lr], dtype='float32')
        pts_out = np.array([[0.0, 0.0],    # upperleft
                           [175, 0.0],  # upperright
                           [0, 134],    # lowerleft
                           [175, 134], # lowerright
                           ], dtype='float32')
        HRI = cv2.getPerspectiveTransform(pts_img, pts_out)
        return HRI
        Kinv = np.linalg.inv(K)
        p_ul = util_camera.homo2pt(np.dot(Kinv, util_camera.pt2homo(pp_ul)))
        p_ur = util_camera.homo2pt(np.dot(Kinv, util_camera.pt2homo(pp_ur)))
        p_ll = util_camera.homo2pt(np.dot(Kinv, util_camera.pt2homo(pp_ll)))
        p_lr = util_camera.homo2pt(np.dot(Kinv, util_camera.pt2homo(pp_lr)))
        # Now, p_ul,p_ur,p_ll,p_lr are normalized coords
        w = p_lr[0] - p_ul[0]
        h = p_lr[1] - p_ul[1]
        pts_I = np.array([[0.0, 0.0],    # upperleft
                          [0.2, 0.0],      # upperright
                          [0.0, 0.4],      # lowerleft
                          [0.2, 0.4]])       # lowerright
        p1_w = homo2pt(np.dot(H, np.hstack((p_ul, [1.0]))))
        p2_w = homo2pt(np.dot(H, np.hstack((p_ur, [1.0]))))
        p3_w = homo2pt(np.dot(H, np.hstack((p_ll, [1.0]))))
        p4_w = homo2pt(np.dot(H, np.hstack((p_lr, [1.0]))))
        pts_w = np.array([p1_w,
                          p2_w,
                          p3_w,
                          p4_w])
        pts_w = pts_w.astype('float32')
        pts_I = pts_I.astype('float32')
        HRI = cv2.getPerspectiveTransform(pts_w, pts_I)
        print homo2pt(np.dot(HRI, pt2homo(p3_w)))
        print pts_I[3]
        return HRI

    HRI = compute_plane2I(p_ul, p_ur, p_ll, p_lr)
    #T = np.dot(HRI, Hinv)
    T = HRI
    print 'T is:'
    print T
    I = cv2.imread(imgpath, cv.CV_LOAD_IMAGE_COLOR).astype('float64')
    I = (2.0*I) + 40    # Up that contrast! Orig. images are super dark.
    print "Imin, Imax:", np.min(I), np.max(I)
    #Ipatch = I[p_ul[1]-50:p_lr[1]+50, p_ul[0]-40:p_lr[0]+40]
    Ipatch = I
    Iwarp = cv2.warpPerspective(Ipatch, T, None)
    print "Iwarpmin, Iwarpmax:", np.min(Iwarp), np.max(Iwarp)
    print "(Displaying before/after images, press <any key> to exit.)"
    cv2.namedWindow('original')
    cv2.imshow('original', I.astype('uint8'))
    cv2.namedWindow('corrected')
    cv2.imshow('corrected', Iwarp.astype('uint8'))
    cv2.waitKey(0)

def test_koopa(SHOW_EPIPOLAR=False):
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
            
        # This doesn't work, since I don't have T: I have (1/d)*T
        '''
        for ii, pt1 in enumerate(pts1_norm):
            pt2 = pts2_norm[ii]
            pt1_h = np.hstack((pt1, [1.0]))
            pt1_out = np.dot(R, pt1_h) + Ts
            pt1_out /= pt1_out[2]

            print pt1
            print pt1_out
            print pt2
            print
        '''    
    
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
    #V = Vt
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
        util_camera.draw_line(Irgb1_, epiline1)
        util_camera.draw_line(Irgb2_, epiline2)
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
    print "======== (2) test_koopa_singleimage: single view homography ========"
    test_koopa_singleimage(SHOW_EPIPOLAR=args.show_epipolar)

if __name__ == '__main__':
    main()
