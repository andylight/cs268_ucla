import sys, os, pdb, argparse
import cv2, cv, numpy as np, scipy.misc

import calibrate_camera, util, util_camera, transform_image
from util import intrnd, tupint
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
    ## Pixel coordinates of a region of the image known to:
    ##     a.) Lie on the plane
    ##     b.) Has known metric lengths
    Hinv = np.linalg.inv(H)
    # Define hardcoded pixel/world coordinates on the poster plane
    # TODO: We could try generalizing this to new images by:
    #     a.) Asking the user to choose four image points, and enter
    #         in metric length information (simplest)
    #     b.) Automatically detect four corner points that correspond
    #         to a world rectangle (somehow). If there are no 
    #         reference lengths, then I believe the best thing we can
    #         do is do a perspective-correction up to an affine 
    #         transformation (since we don't know the correct X/Y
    #         lengths).
    pts_pix = np.array([[372.0, 97.0],     # Upperleft of redbox (x,y)
                        [651.0, 95.0],     # Upperright of redbox (x,y)
                        [303.0, 321.0],    # Lowerleft of bluebox (x,y)
                        [695.0, 324.0]],   # Lowerright of greenbox (x,y)
                       )
    pts_world = np.zeros([0, 2])
    # Populate pts_world via H (maps image plane -> world plane)
    for pt_pix in pts_pix:
        pt_world = homo2pt(np.dot(H, np.dot(Kinv, pt2homo(pt_pix))))
        pts_world = np.vstack((pts_world, pt_world))
    pts_world *= 1000.0 # Let's get a reasonable-sized image please!
    
    # These are hardcoded world coords, but, we can auto-gen them via
    # H, Kinv as above. Isn't geometry nice?
    #pts_world = np.array([[0.0, 0.0],
    #                      [175.0, 0.0],
    #                      [0.0, 134.0],
    #                      [175.0, 134.0]])
    IPM0, IPM1 = compute_IPM(pts_pix, pts_world)
    print 'IPM0 is:'
    print IPM0
    print 'IPM1 is:'
    print IPM1
    I = cv2.imread(imgpath, cv.CV_LOAD_IMAGE_COLOR).astype('float64')
    I = (2.0*I) + 40    # Up that contrast! Orig. images are super dark.
    h, w = I.shape[0:2]
    # Determine bounding box of IPM-corrected image
    a = homo2pt(np.dot(IPM1, pt2homo([0.0, 0.0])))  # upper-left corner
    b = homo2pt(np.dot(IPM1, pt2homo([w-1, 0.0])))  # upper-right corner
    c = homo2pt(np.dot(IPM1, pt2homo([0.0, h-1])))  # lower-left corner
    d = homo2pt(np.dot(IPM1, pt2homo([w-1, h-1])))  # lower-right corner
    w_out = intrnd(max(b[0] - a[0], d[0] - c[0]))
    h_out = intrnd(max(c[1] - a[1], d[1] - b[1]))
    print "New image dimensions: ({0} x {1}) (Orig: {2} x {3})".format(w_out, h_out, w, h)
    # Warp the entire image I with the homography IPM1
    Iwarp = cv2.warpPerspective(I, IPM1, (w_out, h_out))
    # Draw selected points on I
    for pt in pts_pix:
        cv2.circle(I, tupint(pt), 5, (0, 255, 0))
    # Draw rectified-rectangle in Iwarp
    a = homo2pt(np.dot(IPM1, pt2homo(pts_pix[0])))
    b = homo2pt(np.dot(IPM1, pt2homo(pts_pix[3])))
    cv2.rectangle(Iwarp, tupint(a), tupint(b), (0, 255, 0))
    #Iwarp = transform_image.transform_image_forward(I, IPM)
    print "(Displaying before/after images, press <any key> to exit.)"
    cv2.namedWindow('original')
    cv2.imshow('original', I.astype('uint8'))
    cv2.namedWindow('corrected')
    cv2.imshow('corrected', Iwarp.astype('uint8'))
    cv2.waitKey(0)

def compute_IPM(pts_pix, pts_world):
    """ Computes the homography that maps four points in pixel coords
    to four points in world coords.
    This is known as the 'Inverse Perspective Map' (IPM), as it can
    be used to undo perspective effects in an image if you can
    identify four points in the image that correspond to a rectangle
    with known metric lengths.
    Note: for pts_world, you don't have to input exactly the metrix
          lengths: instead, you should scale the dimensions to achieve
          an output image of desired size. For instance, if you input:
          # A rectangle that is (1.2 meters x 1.4 meters):
              pts_world = [[0.0, 1.2],
                           [1.4, 0.0],
                           [0.0, 1.2],
                           [1.4, 1.2]]
          You'll get an output image where the selected region is of
          size (1.2 x 1.4) pixels, which is probably not what you
          wanted! Instead, you can scale the dimensions:
              pts_world *= 100     # Now, region is (120 x 140) pix!
    Input:
        nparray pts_pix: N x 2
        nparray pts_world: N x 2
    Output:
        (nparray IPM0, nparray IPM1)
    Here, IPM0 is the computed that results in the pts_pix as the new
    image origin. This will result in the rest of the image being
    omitted. If you wish to apply the IPM such that the entire image
    is displayed, use IPM1.
    """
    pts_pix = pts_pix.copy().astype('float32')
    pts_world = pts_world.copy().astype('float32')
    IPM0 = cv2.getPerspectiveTransform(pts_pix, pts_world)
    # Determine where (0,0) lands in IPM0's image, and do an offset
    # to ensure that the entire image is displayed with IPM1
    origin_new = homo2pt(np.dot(IPM0, pt2homo([0.0, 0.0])))
    if origin_new[0] < 0:
        pts_world[:, 0] += abs(origin_new[0])
    if origin_new[1] < 0:
        pts_world[:, 1] += abs(origin_new[1])
    IPM1 = cv2.getPerspectiveTransform(pts_pix, pts_world)
    return IPM0, IPM1

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
