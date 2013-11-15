import sys, os, pdb, time, argparse, random
import numpy as np, numpy.linalg as linalg
import cv2

def estimate_line(edgemap, MAX_ITERS=400, T=3.0, ALPHA=8):
    """ Given an edgemap, robustly determine the most dominant line.
    Input:
        nparray edgemap
        int MAX_ITERS
        float T
            Distance threshold between a candidate point and a line
        int ALPHA
            Min. number of inliers required for a model to be considered
    Output:
        tuple line := (float a, float b, float c)
    Where each a,b,c satisfies: ax + by + c = 0
    """
    best_nb_inliers = -np.inf
    best_line = None
    best_inliers = None

    edge_idxs = np.where(edgemap == 255) # (Ys, Xs)
    w = edgemap.shape[1]
    h = edgemap.shape[0]
    nb_active = len(edge_idxs[0])
    if nb_active == 0:
        return None # Couldn't detect any edges!
    
    cnt_iter = 0
    while cnt_iter < MAX_ITERS:
        idx1 = random.randint(0, nb_active - 1)
        idx2 = random.randint(0, nb_active - 1)
        if idx1 == idx2:
            cnt_iter += 1
            continue    # Degenerate case
        pt1 = (edge_idxs[1][idx1], (h-1) - edge_idxs[0][idx1]) # (x, y)
        pt2 = (edge_idxs[1][idx2], (h-1) - edge_idxs[0][idx2]) # subtract from (h-1) to get 'real' coords
        line_init, residual = fit_line((pt1, pt2))
        a, b, c = line_init
        inliers = [idx1, idx2] # Current set of inliers
        # Consider all other unchosen points
        for i in xrange(nb_active):
            if i == idx1 or i == idx2:
                continue
            pt_i = (edge_idxs[1][i], (h-1) - edge_idxs[0][i])
            d = np.abs(a*pt_i[0] + b*pt_i[1] + c) / np.sqrt(a**2 + b**2)
            if d <= T:
                inliers.append(i) # This point is close enough to our fitted line
        # We have a set of candidate inliers
        if len(inliers) < ALPHA:
            cnt_iter += 1
            continue # This model is probably junk
        elif len(inliers) > best_nb_inliers:
            # This is the best model so far!
            pts = []
            for idx in inliers:
                pt = (edge_idxs[1][idx], (h-1) - edge_idxs[0][idx])
                pts.append(pt)
            line, residual = fit_line(pts)
            best_nb_inliers = len(inliers)
            best_line = line
            best_inliers = inliers
        cnt_iter += 1
    return best_line, best_inliers

def fit_line(pts):
    """ Fits a 2D line through 2D points in pts.
    Input:
        tuple pts: ((int x_i, int y_i), ...)
    Output:
        tuple line: (float a, float b, float c)
    """
    if not pts or len(pts) <= 1:
        raise Exception("Can't fit line with less than 2 points!")
    if len(pts) == 2:
        # Solve analytically
        pt1, pt2 = pts
        a = -(pt2[1] - pt1[1])
        b = pt2[0] - pt1[0]
        c = -a*pt1[0] - b*pt1[1]
        return np.array([a, b, c]), 0.0
    else:
        # Solve min || Ap || via SVD of A
        n = len(pts)
        A = np.zeros((n, 3))
        for i, pt in enumerate(pts):
            A[i, :] = (pt[0], pt[1], 1)
        U, S, V = linalg.svd(A) 
        v = V[2, :]
        residual = linalg.norm(np.dot(A, v.T))
        return v, residual

def test_fitline():
    cases = [((1, 1), (2, 2)),
             ((0, 1), (1, 2)),
             ((1, 1), (2, 2), (3, 3)),
             ((1, 1), (2, 2), (3.1, 3.02)),
             ((1, 1), (2, 2), (3.1, 3.02), (3.98, 4.1)),
             ((0, 0), (1, 1), (9, 2), (3, 4))]

    for i, pts in enumerate(cases):
        line, residual = fit_line(pts)
        print("({0}/{1}) Line: {2}".format(i+1, len(cases), line))
        print("    Residual: {0}".format(residual))

def main():
    test_fitline()
    print("Done.")

if __name__ == '__main__':
    main()
