import sys, os, time, pdb, argparse
import numpy as np, cv2

import util

from estimate_line import estimate_line
from util import intrnd

def detect_lanes(I, win1=(0.4, 0.55, 0.2, 0.1), win2=(0.6, 0.55, 0.2, 0.1),
                 threshold1=50, threshold2=100, apertureSize=3):
    """ Given a street image I, detect the (parallel) road lanes
    in image coordinates.
    Input:
        nparray I:
        tuple win1, win2: (float x, float y, float width, float height)
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
    x_left = intrnd(win1[0]*w)
    y_left = intrnd(win1[1]*h)
    w_left = intrnd(win1[2]*w)
    h_left = intrnd(win1[3]*h)
    if w_left % 2 == 0:
        w_left += 1
    if h_left % 2 == 0:
        h_left += 1
    # Compute right window dimensions
    x_right = intrnd(win2[0]*w)
    y_right = intrnd(win2[1]*h)
    w_right = intrnd(win2[2]*w)
    h_right = intrnd(win2[3]*h)
    if w_right % 2 == 0:
        w_right += 1
    if h_right % 2 == 0:
        h_right += 1
    
    edges_left = edgemap[(y_left-(h_left/2)):(y_left+(h_left/2)),
                         (x_left-(w_left/2)):(x_left+(w_left/2))]
    edges_right = edgemap[(y_right-(h_right/2)):(y_right+(h_right/2)),
                          (x_right-(w_right/2)):(x_right+(w_right/2))]

    # Find dominant line in each window
    res1 = estimate_line(edges_left)
    res2 = estimate_line(edges_right)
    if res1 == None:
        line1, inliers1 = None, None
    else:
        line1, inliers1 = res1
    if res2 == None:
        line2, inliers2 = None, None
    else:
        line2, inliers2 = res2

    if line1 != None and line1[1] != 0:
        line1_norm = np.array([line1[0] / line1[1], 1, line1[2] / line1[1]])
    else:
        line1_norm = line1
    if line2 != None and line2[1] != 0:
        line2_norm = np.array([line2[0] / line2[1], 1, line2[2] / line2[1]])
    else:
        line2_norm = line2
    # Fix line to be in image coordinate system (not window coord sys)
    if line1_norm != None:
        a1, b1, c1 = line1_norm
        c1_out = -a1*(x_left-(w_left/2)) - b1*((h-1-y_left)-(h_left/2)) + c1
        line1_out = np.array([a1, b1, c1_out])
    else:
        line1_out = None
    if line2_norm != None:
        a2, b2, c2 = line2_norm
        c2_out = -a2*(x_right-(w_right/2)) - b2*((h-1-y_right)-(h_right/2)) + c2
        line2_out = np.array([a2, b2, c2_out])
    else:
        line2_out = None
    return line1_out, line2_out

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
    
def overlay_line(I, line, colour=(255, 0, 0), thick=1):
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

def draw_subwindow(Irgb, win, colour=(125, 125, 0)):
    """ Draws subwindow on Irgb.
    Input:
        nparray Irgb
        tuple win: (float x, float y, float w, float h)
    Output:
        nparray Irgb_out
    """
    xf, yf, wf, hf = win
    Irgb_out = Irgb.copy()
    h, w = Irgb_out.shape[0:2]
    win_w = intrnd(w*wf)
    win_h = intrnd(h*hf)
    if win_w % 2 == 0:
        win_w += 1
    if win_h % 2 == 0:
        win_h += 1
    pt1 = (intrnd(w*xf - (win_w/2)), intrnd(h*yf - (win_h/2)))
    pt2 = (intrnd(w*xf + (win_w/2)), intrnd(h*yf + (win_h/2)))
    cv2.rectangle(Irgb_out, pt1, pt2, colour, thickness=1)
    return Irgb_out

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("imgsdir", help="Directory of street images (or a single imagepath).")
    parser.add_argument("--Tlow", type=int, help="LowThreshold (Canny).",
                        default=100)
    parser.add_argument("--ksize", type=int, help="Size of Sobel filter (Canny).",
                        default=3)
    parser.add_argument("--win1", nargs=4, type=float, metavar=("X", "Y", "W", "H"),
                        help="Left subwindow (X, Y, WIDTH, HEIGHT). Each \
value is in range [0.0,1.0]",
                        default=(0.4, 0.55, 0.2, 0.1))
    parser.add_argument("--win2", nargs=4, type=float, metavar=("X", "Y", "W", "H"),
                        help="Right subwindow. (See --win1)",
                        default=(0.6, 0.55, 0.2, 0.1))
    parser.add_argument("--n", type=int, help="Number of images to process.")
    return parser.parse_args()

def main():
    args = parse_args()
    threshold1 = args.Tlow
    threshold2 = 2 * args.Tlow    # Canny recommends a ratio of 1:2
    win1 = args.win1
    win2 = args.win2
    imgsdir = args.imgsdir
    if not os.path.isdir(imgsdir):
        imgpaths = [imgsdir]
    else:
        imgpaths = util.get_imgpaths(imgsdir, n=args.n)
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
        # Draw subwindows on image
        Irgb = draw_subwindow(Irgb, win1, colour=(255, 0, 0))
        Irgb = draw_subwindow(Irgb, win2, colour=(0, 255, 0))
        cv2.imwrite('{0}_lines.png'.format(util.get_filename(imgpath)), Irgb)
    print("Done.")

if __name__ == '__main__':
    main()
