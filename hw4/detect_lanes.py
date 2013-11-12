import sys, os, time, pdb, argparse
import numpy as np, cv2

import util

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

    # Find dominant line in each window
    line1, inliers1 = estimate_line(edges_left)
    line2, inliers2 = estimate_line(edges_right)
    if line1[1] != 0:
        line1_norm = np.array([line1[0] / line1[1], 1, line1[2] / line1[1]])
    else:
        line1_norm = line1
    if line2[1] != 0:
        line2_norm = np.array([line2[0] / line2[1], 1, line2[2] / line2[1]])
    else:
        line2_norm = line2
    # Fix line to be in image coordinate system (not window coord sys)
    a1, b1, c1 = line1_norm
    c1_out = -a1*(x_left-(w_left/2)) - b1*((h-1-y_left)-(h_left/2)) + c1
    
    a2, b2, c2 = line2_norm
    c2_out = -a2*(x_right-(w_right/2)) - b2*((h-1-y_right)-(h_right/2)) + c2
    return np.array([a1, b1, c1_out]), np.array([a2, b2, c2_out])

def plot_lines(I, line1, line2):
    """ Plots lines on the input image.
    Input:
        nparray I
        tuple line1, line2: (float a, float b, float c)
    Output:
        nparray Iout
    """
    Irgb = util.to_rgb(I)
    if line1 != None:
        Irgb = overlay_line(Irgb, line1, colour=(255, 0, 0))
    if line2 != None:
        Irgb = overlay_line(Irgb, line2, colour=(0, 255, 0))
    return Irgb
    
def overlay_line(I, line, colour=(255, 0, 0), thick=2):
    """ Adds a colored line to the input image.
    Input:
        nparray I: 
            If I is not RGB, this will convert it.
        tuple line: (float a, float b, float c)
        tuple clr: (int R, int G, int B)
    Output:
        nparray I_rgb
    """
    a, b, c = line # ax + by + c = 0 => y = (-ax - c) / b
    Irgb = util.to_rgb(I)
    h, w = Irgb.shape[0:2]
    for x in xrange(w):
        y = (-a*x - c) / b
        y = (h-1) - y # Remember: lines are in 'normal' coord system,
                      # not image coord system (where origin is UL-corner)
        if y < 0:
            continue
        y_low = max(0, y - thick)
        y_high = min(h-1, y + thick)
        x_low = max(0, x - thick)
        x_high = min(w-1, x + thick)
        Irgb[y_low:y_high, x_low:x_high, 0] = colour[0]
        Irgb[y_low:y_high, x_low:x_high, 1] = colour[1]
        Irgb[y_low:y_high, x_low:x_high, 2] = colour[2]

    return Irgb

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
        if line2 == None:
            print("    Error: Couldn't find right lane.")
        Irgb = plot_lines(I, line1, line2)
        cv2.imwrite('{0}_lines.png'.format(get_filename(imgpath)), Irgb)
    print("Done.")

def get_filename(fpath):
    return os.path.splitext(os.path.split(fpath)[1])[0]

if __name__ == '__main__':
    main()
