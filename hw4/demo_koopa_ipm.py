import sys, os, pdb, argparse
import cv2, cv, numpy as np, scipy.misc

import calibrate_camera, util, util_camera
from util import intrnd
from util_camera import pt2homo, homo2pt

"""
USAGE:

    $ python demo_koopa_ipm.py [-h] [--show_epipolar]

Two quick demos about planar scenes and homographies. The test images
are photos of a poster taken from a variety of viewpoints.

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
        - The camera matrix K
        - The scene is planar
        - The location of four points in the image that create a 
          rectangle in the world plane with known metric lengths
          Note: I don't have to necessarily know the metric lengths,
                but at least I need to know length ratios to construct
                a corrected image with length ratios preserved.
    In the literature, this is also known as: perspective rectification,
    perspective correction.
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
    Kinv = np.linalg.inv(K)

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
    def compute_IPM(pp_ul, pp_ur, pp_ll, pp_lr):
        """ Computes homography that removes the perspective effects
        from a portion of the image. I believe some call this the 
        'Inverse Perspective Map' (IPM). We require that we know the
        world-shape (i.e. a box) of the image region. In this case, we
        know that the pixel coordinates describe a box with known
        metric lengths.
        Input points are image pixel coords.

        Note: This doesn't use the computed H. Why? Well, I hardcoded
        the worldpoints in pts_out: in general, you can use H, K to 
        compute the pts_out, i.e.:
            # pp_ur is pixel coord of upper-right corner of redbox
            p_ur = np.dot(Kinv, pp_ur)  # p_ur is normalized coord
            p_world_ur = np.dot(H, p_ur)
            p_world_ur /= p_world_ur[2] # world plane coord of redbox upper-right corner
        Now, p_world_ur is [0.175, 0.0, 1.0], as expected.
        """
        pts_img = np.array([pp_ul, pp_ur, pp_ll, pp_lr], dtype='float32')
        pts_out = np.array([[0.0, 0.0],    # upperleft
                           [175, 0.0],     # upperright, 17.5 cm * 10
                           [0.0, 134],       # lowerleft   13.4 cm * 10
                           [175, 134],     # lowerright
                           ], dtype='float32')
        #pts_out[:,0] += 100 # Translate points forward 100 units to show more of the left poster
        #pts_out[:,1] += 100 # Translate points down 100 units to show more of the top of the poster
        IPM = cv2.getPerspectiveTransform(pts_img, pts_out)
        return IPM

    IPM = compute_IPM(p_ul, p_ur, p_ll, p_lr)
    print 'IPM is:'
    print IPM
    I = cv2.imread(imgpath, cv.CV_LOAD_IMAGE_COLOR).astype('float64')
    I = (2.0*I) + 40    # Up that contrast! Orig. images are super dark.
    # Ipatch = I[p_ul[1]:p_lr[1], p_ul[0]:p_lr[1]]
    # Let's warp the entire image I with the homography, rather than
    # just the Ipatch.
    Iwarp = cv2.warpPerspective(I, IPM, None)
    print "(Displaying before/after images, press <any key> to exit.)"
    cv2.namedWindow('original')
    cv2.imshow('original', I.astype('uint8'))
    cv2.namedWindow('corrected')
    cv2.imshow('corrected', Iwarp.astype('uint8'))
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
    return parser.parse_args()

def main():
    args = parse_args()
    test_koopa_singleimage()

if __name__ == '__main__':
    main()
