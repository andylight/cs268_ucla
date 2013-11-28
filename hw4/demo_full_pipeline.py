import sys, os, pdb, argparse, time
import util_camera, util
import numpy as np, cv2, cv

import calibrate_camera, detect_lanes
from util import intrnd
from util_camera import compute_x, compute_y, pt2homo, homo2pt

"""
USAGE:

    $ python demo_full_pipeline.py

Demos a simple lane-departure warning system on the test images in:
    LDWS_test_short/

The script first determines the camera calibration parameters by
processing the planar calibration object in:
    LDWS_calibrate_short/

Then, each test image is processed, and after the lanes+extrinsic
parameters are estimated, the following output appears:
    1.) On stdout, the estimated camera position w.r.t. the middle
        of the lane (along with a warning if the camera position is
        too far away from the center).
    2.) Two image windows display. One with the detected lane positions,
        and another with a top-down ("birds-eye") view of the road,
        obtained by undo-ing the perspective distortion.
"""

IMGSDIR_CALIB = 'LDWS_calibrate_short'
IMGSDIR_TEST  = 'LDWS_test_short'

WIN_LEFT = (0.22, 0.55, 0.25, 0.1)
WIN_RIGHT = (0.47, 0.55, 0.25, 0.1)

LANE_W = 3.66 # 3.66 meters

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--imgsdir", default=IMGSDIR_TEST,
                        help="Directory of test images. Can be a single \
image path.")
    parser.add_argument("--reuse_calib", action='store_true', default=False,
                        help="Use a precomputed camera calibration matrix, \
rather than re-computing it.")
    return parser.parse_args()

def main():
    args = parse_args()
    if not os.path.exists(IMGSDIR_CALIB):
        print "(ERROR) Calibration images not found. Please place \
the LDWS_calibrate_short/ images in the current directory."
        exit(1)
    if not os.path.exists(IMGSDIR_TEST):
        print "(ERROR) Test images not found. Please place the \
LDWS_test_short/ images in the current directory."
        exit(1)

    imgpaths_calib = util.get_imgpaths(IMGSDIR_CALIB)
    if not os.path.isdir(args.imgsdir):
        imgpaths_test = [args.imgsdir]
    else:
        imgpaths_test  = util.get_imgpaths(args.imgsdir)
    if args.reuse_calib:
        print "(Reusing camera calibration matrix)"
        K = np.array([[ 674.07224154,    0.,          262.77722917],
                      [   0.,          670.26875783,  330.21546389],
                      [   0.,            0.,            1.        ]])
    else:
        print "(Estimating camera matrix...)"
        t = time.time()
        K = calibrate_camera.calibrate_camera(imgpaths_calib, 8, 8, 0.048)
        dur = time.time() - t
        print "(Finished. {0:.4f})".format(dur)

    for i, imgpath in enumerate(imgpaths_test):
        print "\n==== ({0}/{1}) Detecting lanes... [{2}]====".format(i+1, len(imgpaths_test), os.path.split(imgpath)[1])
        I = cv2.imread(imgpath, cv2.CV_LOAD_IMAGE_GRAYSCALE)
        h, w = I.shape[0:2]
        t = time.time()
        line1, line2 = detect_lanes.detect_lanes(I, win1=WIN_LEFT, win2=WIN_RIGHT,
                                                 threshold1=110, threshold2=220,
                                                 apertureSize=3,
                                                 show_edges=False)
        dur = time.time() - t
        print "    Finished detecting lanes ({0:.4f}s)".format(dur)
        if line1 == None or line2 == None:
            print "({0}/{1}) Error: Couldn't find lanes.".format(i+1, len(imgpaths_test))
            continue

        # Choose 4 points on the lanes to estimate the planar homography
        y1 = intrnd(0.45 * h)
        y2 = intrnd(0.65 * h)
        pts = np.array([[compute_x(line1, y1), y1],    # Left lane, far
                        [compute_x(line2, y1), y1],    # Right lane, far
                        [compute_x(line1, y2), y2],    # Left lane, close
                        [compute_x(line2, y2), y2]])   # Right lane, close
        # These world points have the origin at the middle of the lane,
        # directly below the camera.
        # Depths 5.0, 3.0 are chosen arbitrarily, as we don't have depth
        # information. Thus, the computed H will only be accurate in the
        # X axis.
        pts_worldm = np.array([[-LANE_W/2.0, 5.0],           # Left lane, far
                               [LANE_W/2.0, 5.0],            # Right lane, far
                               [-LANE_W/2.0, 3.0],           # Left lane, close
                               [LANE_W/2.0, 3.0]])           # Right lane, close
        # These world points are scaled+translated to generate a 
        # reasonably-sized image w/ perspective effects removed.
        # Note: the origin here is not the same origin as in pts_worldm!
        pts_world = np.array([[0.0, 0.0],
                              [LANE_W*1e2, 0.0],
                              [0.0, 500],
                              [LANE_W*1e2, 500.0]])
        pts_world[:,0] += LANE_W*1e2 # Let's see more of the horizontal image
        pts_world[:,1] += 200 # Let's see more down the road
        # H_metric := Homography mapping the image (pixel coords) to the 
        #             world plane defined by pts_worldm
        H_metric = cv2.getPerspectiveTransform(pts.astype('float32'), pts_worldm.astype('float32'))
        H = cv2.getPerspectiveTransform(pts.astype('float32'), pts_world.astype('float32'))
        
        ## Estimate where the camera is w.r.t. the world ref. frame
        R, T = estimate_extrinsic_parameters(H_metric)
        xdist = T[0] - (LANE_W / 2.0)
        print "    Distance from center of lane: X={0:.2f} meters".format(xdist)
        LEFT_THRESH = -1.0    # Stay within 1.0 meters of the center of the lane
        RIGHT_THRESH = 1.0
        if xdist >= RIGHT_THRESH:
            print "        WARNING: Camera center is awfully close to the \
RIGHT side of the lane!"
        elif xdist <= LEFT_THRESH:
            print "        WARNING: Camera center is awfully close to the \
LEFT side of the lane!"

        Irgb = cv2.imread(imgpath, cv2.CV_LOAD_IMAGE_COLOR)
        Iipm = cv2.warpPerspective(Irgb.astype('float64'), H, (1000, 700))
        cv2.namedWindow("win2: Perspective-rectified image")
        cv2.imshow("win2: Perspective-rectified image", Iipm.astype('uint8'))

        Irgb = detect_lanes.draw_subwindow(Irgb, WIN_LEFT)
        Irgb = detect_lanes.draw_subwindow(Irgb, WIN_RIGHT)
        for _pt in pts:
            _pt = tuple([intrnd(x) for x in _pt])
            cv2.circle(Irgb, _pt, 3, (0, 0, 255))

        print "    ({0}/{1}) Displaying detected lanes.".format(i+1, len(imgpaths_test))
        show_lanes(Irgb, line1, line2)

    print "Done."

def estimate_extrinsic_parameters(H):
    """ Outputs the (R,T) relative to the world reference frame. """
    #### Normalize H := H / sigma_2(H_)
    U, S, V = np.linalg.svd(H)
    H = H / S[1]

    decomps = util_camera.decompose_H(H)
    R = decomps[0][0]
    T = decomps[0][1]
    return R, T / T[2]

def show_lanes(Irgb, line1, line2):
    Irgb = util_camera.draw_line(Irgb, line1)
    Irgb = util_camera.draw_line(Irgb, line2)
    cv2.namedWindow("win1: Detected Lanes")
    cv2.imshow("win1: Detected Lanes", Irgb)
    print "(press <enter> to continue, with the imdisplay window active)"
    cv2.waitKey(0)

if __name__ == '__main__':
    main()
