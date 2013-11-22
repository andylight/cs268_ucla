import sys, os, pdb, argparse

import util
import numpy as np, cv, cv2

from util import intrnd

def calibrate_camera(imgpaths, rows, cols, boxdim, SHOW_CB=False):
    """ Determines intrinsic camera matrix K from a set of images of a
    calibration pattern (checkerboard pattern).
    Input:
        tuple imgpaths
        int rows, cols
            The number of inner-corners the checkerboard has.
        float boxdim
            The world length of a box side (in metric units, i.e., mm)
        boolean SHOW_CB
            If True, then we show the checkerboard results in an
            interactive manner.
    Output:
        nparray K
    """
    object_pts = []
    image_pts = []
    w_img, h_img = None, None
    for i, imgpath in enumerate(imgpaths):
        I = cv2.imread(imgpath, cv2.CV_LOAD_IMAGE_GRAYSCALE)
        if w_img == None:
            w_img = I.shape[1]
            h_img = I.shape[0]
        # nparray corners: N x 1 x 2 matrix of coords for all N corners
        retval, corners = cv2.findChessboardCorners(I, (rows, cols),
                                                    flags=cv.CV_CALIB_CB_ADAPTIVE_THRESH | 
                                                          cv.CV_CALIB_CB_NORMALIZE_IMAGE)
        if retval == 0:
            if corners == None:
                print "(i={0}) Warning: could not find any corners (0/{1})".format(i, rows*cols)
            else:
                print "(i={0}) Warning: could not find all corners ({1}/{2})".format(i, len(corners), rows*cols)
        else:
            pass
            #print "(i={0}) Found all corners! ({1}/{1})".format(i, rows*cols)
        if SHOW_CB:
            Irgb = cv2.cvtColor(I, cv.CV_GRAY2RGB)
            cv2.drawChessboardCorners(Irgb, (rows, cols), corners, retval)
            cv2.namedWindow('display')
            cv2.imshow('display', Irgb)
            cv2.waitKey(0)
        if retval == 0:
            continue
        cb_pts = compute_cb_pts(corners, rows, cols, boxdim)
        object_pts.append(cb_pts)
        image_pts.append(corners)
    print "Found all corners in {0}/{1} views".format(len(object_pts), len(imgpaths))
    # object_mat: M x N x 1 x 3 (M is # of views, N is # of corners)
    object_mat = np.zeros([len(object_pts), rows*cols, 1, 3])
    # image_mat: M x N x 1 x 2
    image_mat = np.zeros([len(image_pts), rows*cols, 1, 2])
    for i, objpts in enumerate(object_pts):
        object_mat[i,:,:,:] = objpts
    for i, imgpts in enumerate(image_pts):
        image_mat[i,:,:,:] = imgpts
    print "(Info) Calling cv2.calibrateCamera..."
    retval_calib, K, distCoeffs, rvecs, tvecs = cv2.calibrateCamera(object_mat.astype('float32'),
                                                                    image_mat.astype('float32'),
                                                                    (w_img, h_img))
    print "retval_calib:", retval_calib
    return K
        
def compute_cb_pts(corners, rows, cols, boxdim):
    """ Given the pixel coords of the corners, output world coords,
    i.e. with Z = 0 (where we set the world frame onto the calibration
    plane, with the origin on the upper-left-corner-most corner.
    Input:
        nparray corners: N x 1 x 2
        float boxdim
    Output:
        nparray corners_world: N x 1 x 3
    """
    out = np.zeros([corners.shape[0], 1, 3])
    for i in xrange(rows):
        # Y axis points along columns
        Y = i * boxdim
        for j in xrange(cols):
            # X axis points along rows
            X = j * boxdim
            idx = j + i*cols
            imgpt = corners[idx]
            out[idx, 0, 0] = X
            out[idx, 0, 1] = Y
            out[idx, 0, 2] = 0 # Z = 0
    return out

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("imgsdir", 
                        help="Directory of calibration images.")
    parser.add_argument("--patternsize", nargs=2, type=int,
                        metavar=("rows", "cols"),
                        help="Number of inner corners in the checkerboard.",
                        default=(8, 8))
    parser.add_argument("--boxdim", type=float,
                        help="Dimension of a box on the checkerboard, in \
metric units (in m)",
                        default=0.048)
    parser.add_argument("--show_cb", action='store_true', default=False,
                        help="Interactively display checkerboard.")
    return parser.parse_args()

def main():
    args = parse_args()
    imgsdir = args.imgsdir
    imgpaths = util.get_imgpaths(imgsdir)
    print "(Info) Processing {0} images".format(len(imgpaths))
    rows, cols = args.patternsize
    K = calibrate_camera(imgpaths, rows, cols, args.boxdim, SHOW_CB=args.show_cb)
    print "Computed K:"
    print K
    print "Done."

if __name__ == '__main__':
    main()
