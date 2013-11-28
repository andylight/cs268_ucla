import sys, os, pdb, argparse
import util_camera, util
import numpy as np, cv2, cv

import calibrate_camera, detect_lanes
from util import intrnd
from util_camera import compute_x, compute_y, pt2homo, homo2pt

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
    return parser.parse_args()

def main():
    args = parse_args()
    imgpaths_calib = util.get_imgpaths(IMGSDIR_CALIB)
    if not os.path.isdir(args.imgsdir):
        imgpaths_test = [args.imgsdir]
    else:
        imgpaths_test  = util.get_imgpaths(args.imgsdir)
    if True:
        K = np.array([[ 674.07224154,    0.,          262.77722917],
                      [   0.,          670.26875783,  330.21546389],
                      [   0.,            0.,            1.        ]])
    else:
        K = calibrate_camera.calibrate_camera(imgpaths_calib, 8, 8, 0.048)

    for i, imgpath in enumerate(imgpaths_test):
        print "({0}/{1}) Detecting lanes...".format(i+1, len(imgpaths_test))
        I = cv2.imread(imgpath, cv2.CV_LOAD_IMAGE_GRAYSCALE)
        h, w = I.shape[0:2]
        line1, line2 = detect_lanes.detect_lanes(I, win1=WIN_LEFT, win2=WIN_RIGHT)
        if line1 == None or line2 == None:
            print "({0}/{1}) Error: Couldn't find lanes.".format(i+1, len(imgpaths_test))
            continue

        # Choose 4 points on the lanes
        y1 = intrnd(0.45 * h)
        y2 = intrnd(0.65 * h)
        pts = np.array([[compute_x(line1, y1), y1],    # Left lane, far
                        [compute_x(line2, y1), y1],    # Right lane, far
                        [compute_x(line1, y2), y2],    # Left lane, close
                        [compute_x(line2, y2), y2]])   # Right lane, close
        # Depths 5.0, 3.0 are chosen arbitrarily.
        pts_worldm = np.array([[-1.83, 5.0],
                               [1.83, 5.0],
                               [-1.83, 3.0],
                               [1.83, 3.0]])
        pts_world = np.array([[0.0, 0.0],
                              [366, 0.0],
                              [0.0, 500],
                              [366, 500.0]])
        #pts_world[:,1] += 200
        #H, mask = cv2.findHomography(pts_world, pts)
        H_metric = cv2.getPerspectiveTransform(pts.astype('float32'), pts_worldm.astype('float32'))
        H = cv2.getPerspectiveTransform(pts.astype('float32'), pts_world.astype('float32'))
        
        R, T = estimate_extrinsic_parameters(H_metric)
        print "Rotation:"
        print R
        print "Translation:"
        print T

        Irgb = cv2.imread(imgpath, cv2.CV_LOAD_IMAGE_COLOR)
        Iipm = cv2.warpPerspective(Irgb.astype('float64'), H, None)
        cv2.namedWindow("win2: Perspective-rectified image")
        cv2.imshow("win2: Perspective-rectified image", Iipm.astype('uint8'))

        Irgb = detect_lanes.draw_subwindow(Irgb, WIN_LEFT)
        Irgb = detect_lanes.draw_subwindow(Irgb, WIN_RIGHT)
        for _pt in pts:
            _pt = tuple([intrnd(x) for x in _pt])
            cv2.circle(Irgb, _pt, 3, (0, 0, 255))

        print "({0}/{1}) Displaying detected lanes.".format(i+1, len(imgpaths_test))
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
    cv2.waitKey(0)

if __name__ == '__main__':
    main()
