import os, sys, pdb

def parse_img_data(lines):
    lines = lines[1:] # Skip header
    xs, ys = [], []
    k = 1               # lines chomped
    for i, line in enumerate(lines):
        if line.startswith('x'):
            break
        elif i % 2 == 0:
            xs.extend([int(x) for x in line.strip().split(' ')])
        else:
            ys.extend([int(x) for x in line.strip().split(' ')])
        k += 1
    return xs, ys, k
    
def main():
    args = sys.argv[1:]
    fpath = args[0]
    f = open(fpath, 'r')
    lines = f.readlines()
    f.close()
    xs_all, ys_all = [], []
    k_cur = 0
    while True:
        xs, ys, k = parse_img_data(lines[k_cur:])
        xs_all.append(xs)
        ys_all.append(ys)
        k_cur += k
        if k_cur >= len(lines):
            break
    
    print "Done getting xs,ys for {0} images.".format(len(xs_all))
    print "    Saving images to separate .dat files"
    PTS_ROOTDIR = 'ptsdata'
    if not os.path.exists(PTS_ROOTDIR):
        os.makedirs(PTS_ROOTDIR)
    for imgid, xs in enumerate(xs_all):
        ys = ys_all[imgid]
        pts = zip(xs, ys)
        f = open(os.path.join(PTS_ROOTDIR, 'I_{0:02d}.dat'.format(imgid+1)), 'w')
        for (x,y) in pts:
            f.write("{0} {1}\n".format(x, y))
        f.close()
    print "Done."

if __name__ == '__main__':
    main()
