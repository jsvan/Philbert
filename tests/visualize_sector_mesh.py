import unittest

import numpy as np
import math
import test_base
from matplotlib import pyplot as plt
from misc import tools, euclidean
from hilbert import polygon, line, geometry
from misc.dual import Polar, Sectors
from misc.euclidean import X, Y, BasicLine, Point
#from operator import xor

three_colors = ['blue', 'orange', 'green']

class MyTestCase(unittest.TestCase):

    def test_loss_space(self):
        while True:
            try:
                self.visualize_loss_space()
            except Exception as e:
                print(e, e.__doc__)



    def visualize_loss_space(self):
        # For some reason, radius cant be some kind of multiple of width.
        width = 0.34
        # tointmultiple is the number you need to multiply your width by to get to an integer instead of a float.
        tointmultiple = 100
        mapradius = 13
        allXpoints = tools.linspace(-mapradius, mapradius, width)
        gridlength = len(allXpoints)
        i2fkeys = {i:round(f * tointmultiple) for i, f in enumerate(allXpoints)}
        f2ikeys = {round(f * tointmultiple):i for i, f in enumerate(allXpoints)}
        def gridmap(z):
            if isinstance(z, float):
                return f2ikeys[round(z * tointmultiple)]
            return i2fkeys[z] / tointmultiple
        grid = np.ndarray((gridlength, gridlength), dtype=object)
        for i in range(gridlength):
            for j in range(gridlength):
                grid[i][j] = []

        onaught = test_base.omega_for_dual
        j, g = tools.rand_points(2)
        o = test_base.primal_omega_weird
        mainline = line.Line(j, g, o)
        ball = mainline.hilbert_ball_about_line(1.1)
        p, q, r = ball.sample_from_ball_border(3, mainline)
        #plt.subplot(1, 2, 1)
        pnaught, qnaught, rnaught = Polar.to_line(p), Polar.to_line(q), Polar.to_line(r)
        psectors = Sectors(pnaught, onaught)
        qsectors = Sectors(qnaught, onaught)
        rsectors = Sectors(rnaught, onaught)
        lines = [psectors, qsectors, rsectors]
        minpoint = None
        for l in lines:
            regions = l.wedges_to_mosaic(mapradius)
            for region in regions:
                dist = None
                sect = region.sector
                points = region.polygon.gridspace(width)
                for point in points:
                    primal_line = Polar.to_line(point)
                    pside = euclidean.point_below_line(p, primal_line)
                    qside = euclidean.point_below_line(q, primal_line)
                    rside = euclidean.point_below_line(r, primal_line)
                    # Check here that (p, q) are opposites as well as (r, q)
                    # Also that (p, r) are same.
                    if (pside == qside or rside == qside) or (pside != rside):
                        continue
                    if point[X] > mapradius or point[Y] > mapradius:
                        print("skipping", point)
                        continue

                    xidx, yidx = gridmap(point[X]), gridmap(point[Y])
                    measure_point = sect.oXp
                    dist = Polar.dist(BasicLine(sect.wedge.aXp, sect.wedge.bXp),  # the line of our point-naught
                                      BasicLine(point, measure_point),            # the line of our point in space
                                      sect.o_edge,                                # one line from omega
                                      BasicLine(sect.o_tangent, measure_point))  # another line from omega
                    # Color of point a's derivative
                    # doesnt work
                    # plt.scatter(point[X], point[Y], c=str(min(1, Polar.derivative(BasicLine(point, measure_point),
                    #                                                               sect.o_edge,
                    #                                                               BasicLine(sect.o_tangent, measure_point)
                    #                                                              )
                    #                                             )
                    #                                       )
                    #             )

                    grid[yidx][xidx].append(dist)

        def loss_func(ds):
            return min(ds)

        loss = 0
        maxloss = 0
        #maxdist = 0
        for i in range(gridlength):
            for j in range(gridlength):
                three_distances = grid[i][j]
                if len(three_distances) == 3:
                    loss = loss_func(three_distances)
                    #maxdist = max(three_distances[0], maxdist)
                    if loss is not None:

                        # if you want a closest color map:
                        # y, x = gridmap(i), gridmap(j)
                        # plt.scatter(x, y, c=colors[np.argmin(three_distances)])

                        if loss > maxloss:
                            maxloss = loss
                            minpoint = Point((gridmap(j), gridmap(i)))
                    grid[i][j] = loss
                else:
                    grid[i][j] = None

        def avg_not_none(i, j):
            c, t = 0, 0
            for ii, jj in zip([-1, 0, 1, 0], [0, -1, 0, 1]):
                try:
                    t += abs(grid[i + ii][j + jj] - grid[i][j])
                    c += 1
                except:
                    continue
            print(f"returning avg deriv {t/c}")
            return t / c

        def max_not_none(i, j):
            t = 0
            if grid[i][j] is None:
                return None
            for ii, jj in zip([-1,0,1,0], [0,-1,0,1]):
                other = grid[i + ii][j + jj]
                if other is None: continue
                t = max(t, abs(other - grid[i][j]))
            return t

        def dir_not_none(i, j):
            t = 0
            direction = None
            if grid[i][j] is None:
                return None
            for ii, jj in zip([-1, 0, 1, 0], [0, -1, 0, 1]):
                other = grid[i + ii][j + jj]
                if other is None: continue
                if abs(other - grid[i][j]) > t:
                    pass

            return t


        maxD = 0
        derivs = dict()
        plt.subplot(1,2,1)
        for i in range(gridlength):
            for j in range(gridlength):
                if grid[i][j] is None:
                    continue
                try:
                    #loss = grid[i][j] / maxloss
                    deriv = avg_not_none(i, j)
                    #print(f"loss is {loss}, maxloss is {maxloss}, div is {loss/maxloss}")
                    y, x = gridmap(i), gridmap(j)
                    maxD = max(maxD, deriv)
                    #plt.scatter(x, y, c=str(loss)[:3])
                    derivs[(x, y)] = deriv
                except Exception as e:
                    print(e, e.__doc__)
                    continue
        print('maxD', maxD)
        for (x, y), deriv in derivs.items():
            print(f" color {str((deriv / maxD) + 0.05)[:3]}")
            plt.scatter(x, y, c=str((deriv / maxD) + 0.05)[:3])



        for l, c in zip(lines, three_colors):
            sectors = l.sektoren
            for sector in sectors:
                wedge = sector.wedge
                tools.plot_line(plt, *wedge.line_a, color=c, axline=True, linewidth=1, linestyle='dotted')
                tools.plot_line(plt, *wedge.line_b, color=c, axline=True, linewidth=1, linestyle='dotted')
                tools.plot_line(plt, *wedge.line_a, color='white', axline=False, linewidth=2)
                tools.plot_line(plt, *wedge.line_b, color='white', axline=False, linewidth=2)


        tools.plot_congruent(plt, onaught.vertices, color='gray')
        tools.scatter(plt, [Polar.to_point(mainline.l)], colors=['magenta'], size=4)
        tools.scatter(plt, [minpoint], colors='c', size=7)
        for l, c in zip([pnaught, qnaught, rnaught], three_colors):
            tools.plot_line(plt, *l, color='white', axline=True, linewidth=2.5)
            tools.plot_line(plt, *l, color=c, axline=True, linewidth=1.5)

        plt.xlim(-mapradius, mapradius)
        plt.ylim(-mapradius, mapradius)



        plt.subplot(1,2,2)
        tools.plot_congruent(plt, o.vertices, color='gray')

        tools.plot_line(plt, *Polar.to_line(minpoint), color='c', axline=True)
        tools.plot_line(plt, *mainline.l, color='magenta', axline=True)

        tools.scatter(plt, [p,q,r], colors=three_colors)
        tools.annotate(plt, 'brb', [p,q,r])
        for point, col in zip([p,q,r], three_colors):
            vers = o.hilbert_ball_around_point(point, maxloss).vertices
            tools.plot_congruent(plt, vers, color=col)

        plt.show()















if __name__ == '__main__':
    unittest.main()
