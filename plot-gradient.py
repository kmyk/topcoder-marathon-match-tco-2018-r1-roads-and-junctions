# Python Version: 3.x
from mpl_toolkits.mplot3d import Axes3D  # required for fig.add_subplot(111, projection='3d')
import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib.colors import Normalize
import matplotlib.tri as mtri
import numpy as np

# https://stackoverflow.com/a/20146989
class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

def main():
    import argparse
    parser = argparse.ArgumentParser(usage='''
    $ $CXX $CXXFLAGS plot-gradient.cpp -o plot-gradient.bin
    $ file test/$seed.in
    $ N=1000 time ./plot-gradient.bin < test/$seed.in | python3 plot-gradient.py test/$seed.in /dev/stdin --save $seed.png --no-frame
''')
    parser.add_argument('input')
    parser.add_argument('data')
    parser.add_argument('output', nargs='?')
    parser.add_argument('--3d', action='store_true', dest='three')
    parser.add_argument('--no-frame', action='store_true')
    parser.add_argument('--save')
    args = parser.parse_args()

    with open(args.input) as fh:
        S = int(fh.readline())
        NC = int(fh.readline()) // 2
        cities = []
        for _ in range(NC):
            x = int(fh.readline())
            y = int(fh.readline())
            cities += [ ( x, y ) ]
        _ = float(fh.readline())
        _ = float(fh.readline())

    with open(args.data) as fh:
        reference_value = float(fh.readline())
        xs, ys, zs = [], [], []
        for line in fh:
            x, y, z = line.split()
            xs += [ int(x) ]
            ys += [ int(y) ]
            zs += [ reference_value - float(z) ]

    if args.output is not None:
        with open(args.output) as fh:
            NJ = int(fh.readline()) // 2
            junctions = []
            for _ in range(NJ):
                x = int(fh.readline())
                y = int(fh.readline())
                junctions += [ ( x, y ) ]

    # prepare
    if args.three:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.zaxis.set_major_locator(LinearLocator(10))
        ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    else:
        fig, ax = plt.subplots()
    plt.xlim(0, S)
    plt.ylim(0, S)
    ax.invert_yaxis()
    ax.set_aspect('equal')

    if args.three:  # 3D
        # plot function
        triang = mtri.Triangulation(xs, ys)
        surf = ax.plot_trisurf(triang, zs, linewidth=0.1, cmap=plt.cm.Spectral)
        if not args.no_frame:
            fig.colorbar(surf, shrink=0.5, aspect=5)

        # NOTE: https://stackoverflow.com/questions/13932150/matplotlib-wrong-overlapping-when-plotting-two-3d-surfaces-on-the-same-axes
        # plot cities
        # for x, y in cities:
        #     z = 30
        #     ax.scatter(x, y, 0, s=40, c='black')
        #     ax.plot([ x, x ], [ y, y ], [ - z, z ], c='black')

    else:  # 2D
        # plot function
        xs1 = list(sorted(set(xs)))
        ys1 = list(sorted(set(ys)))
        zs1 = [ [ None for _ in xs1 ] for _ in ys1 ]
        for x, y, z in zip(xs, ys, zs):
            zs1[ys1.index(y)][xs1.index(x)] = z
        plt.pcolor(xs1, ys1, zs1, cmap=plt.cm.seismic, norm=MidpointNormalize(midpoint=0))
        if not args.no_frame:
            plt.colorbar()

        # plot cities
        for x, y in cities:
            plt.scatter(x, y, c='black')

        # plot junctions
        for x, y in junctions:
            plt.scatter(x, y, c='white')

    if args.save:
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        fig.savefig(args.save, bbox_inches='tight', pad_inches=0, dpi=2 * 80)
    else:
        plt.show()



if __name__ == '__main__':
    main()
