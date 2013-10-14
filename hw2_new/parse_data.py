import os, sys, pdb, argparse

def parse_img_data(lines):
    lines = lines[1:] # Skip header
    xs, ys = [], []
    k = 1               # lines chomped
    for i, line in enumerate(lines):
        if line.startswith('x'):
            break
        elif i % 2 == 0:
            xs.extend([int(round(float(x))) for x in line.strip().split(' ')])
        else:
            ys.extend([int(round(float(x))) for x in line.strip().split(' ')])
        k += 1
    return xs, ys, k
    
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('rawdata')
    parser.add_argument('outdir')
    return parser.parse_args()
    
def main():
    args = parse_args()
    f = open(args.rawdata, 'r')
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
    print "    Saving images to separate .pts files to: {0}".format(args.outdir)
    outdir = args.outdir
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    for imgid, xs in enumerate(xs_all):
        ys = ys_all[imgid]
        pts = zip(xs, ys)
        f = open(os.path.join(outdir, 'I_{0:02d}.pts'.format(imgid+1)), 'w')
        for (x,y) in pts:
            f.write("{0} {1}\n".format(x, y))
        f.close()
    print "Done."

if __name__ == '__main__':
    main()
