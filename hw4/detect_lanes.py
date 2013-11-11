import sys, os, time, pdb, argparse
import numpy as np, cv2

from estimate_line import estimate_line
from util import intrnd

def detect_lanes(I, win1=(0.2, 0.1, 0.4, 0.55), win2=(0.2, 0.1, 0.6, 0.55),
                 threshold1=50, threshold2=100, apertureSize=3):
    """ Given a street image I, detect the (parallel) road lanes
    in image coordinates.
    Input:
        nparray I:
        tuple win1, win2: (float width, float height, float x, float y)
            Location+sizes of window-left and window-right used to try
            searching for lanes. Dimensions/positions are gived in 
            percentages of image size.
        float threshold1, threshold2
            threshold1 is low-threshold for hysterisis procedure of
            Canny. threshold2 is high-threshold.
        int apertureSize
            One of (1,3,5,7). Size of the Sobel filter.
    Output:
        (line1, line2)
    Where line1 = (a1, b1,c1) such that:
        a1x + b1y + c1 = 0
    Similarly, line2 = (a2, b2, c2).
    """
    h, w = np.shape(I)[0:2]
    
    edgemap = cv2.Canny(I, threshold1, threshold2, apertureSize=apertureSize)

    # Compute left window dimensions
    x_left = intrnd(win1[2]*w)
    y_left = intrnd(win1[3]*h)
    w_left = intrnd(win1[0]*w)
    h_left = intrnd(win1[1]*h)
    if w_left % 2 == 0:
        w_left += 1
    if h_left % 2 == 0:
        h_left += 1
    # Compute right window dimensions
    x_right = intrnd(win2[2]*w)
    y_right = intrnd(win2[3]*h)
    w_right = intrnd(win2[0]*w)
    h_right = intrnd(win2[1]*h)
    if w_right % 2 == 0:
        w_right += 1
    if h_right % 2 == 0:
        h_right += 1
    
    edges_left = edgemap[(y_left-(h_left/2)):(y_left+(h_left/2)),
                         (x_left-(w_left/2)):(x_left+(w_left/2))]
    edges_right = edgemap[(y_right-(h_right/2)):(y_right+(h_right/2)),
                          (x_right-(w_right/2)):(x_right+(w_right/2))]

    cv2.imwrite('edge_left.png', edges_left)
    cv2.imwrite('edge_right.png', edges_right)

    # Find dominant line in each window
    line1, inliers1 = estimate_line(edges_left)
    line2, inliers2 = estimate_line(edges_right)
    return line1, line2

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("imgsdir", help="Directory of street images (or a single imagepath).")
    parser.add_argument("--Tlow", type=int, help="LowThreshold (Canny).",
                        default=100)
    parser.add_argument("--ksize", type=int, help="Size of Sobel filter (Canny).",
                        default=3)
    return parser.parse_args()

def isimgext(path):
    p = path.lower()
    return p.endswith('.png') or p.endswith('.jpeg') or p.endswith('.jpg')

def main():
    args = parse_args()
    threshold1 = args.Tlow
    threshold2 = 2 * args.Tlow    # Canny recommends a ratio of 1:2
    imgsdir = args.imgsdir
    if not os.path.isdir(imgsdir):
        imgpaths = [imgsdir]
    else:
        imgpaths = []
        for dirpath, dirnames, filenames in os.walk(imgsdir):
            for f in [f for f in filenames if isimgext(f)]:
                imgpaths.append(os.path.join(dirpath, f))
    for i, imgpath in enumerate(imgpaths):
        print("({0}/{1}): Image={2}".format(i+1, len(imgpaths), imgpath))
        I = cv2.imread(imgpath, cv2.CV_LOAD_IMAGE_GRAYSCALE)
        line1, line2 = detect_lanes(I, threshold1=threshold1, threshold2=threshold2, apertureSize=args.ksize)
        if line1 == None and line2 == None:
            print("    Error: Couldn't find lanes.")
            continue
        if line1 == None:
            print("    Error: Couldn't find left lane")
            continue
        if line2 == None:
            print("    Error: Couldn't find right lane.")
            continue
        print("Found both lanes!")
    print("Done.")

if __name__ == '__main__':
    main()
