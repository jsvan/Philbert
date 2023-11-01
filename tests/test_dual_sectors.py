import unittest
from hilbert import geometry, omega, polygon, line
from matplotlib import pyplot as plt
from misc import tools, euclidean
from misc.dual import Polar, Sectors
import numpy as np

vertices=[(1.4, -0.1),
          (1.1, 0.8),
          (0.6, 1.1),
          (-0.2, 1.2),
          (-.5, 0.8),
          (-.8, 0.3),
          (-.7, -.2)]
polarO = omega.Omega(Polar.v2v(polygon.Polygon(vertices).vertices))

class PolarTests(unittest.TestCase):

    def test_wedges(self):
        o = polarO
        p = tools.rand_points(1)[0]
        s = Sectors(Polar.to_line(p), o)

        Polar.plot_congruent(plt, o.)
        tools.scatter(plt, s.intersections)
        for i, w in enumerate(s.wedges):
            tools.plot_line(plt, w[0][0]+0.01*i, w[0][1]+0.01*i, color=tools.COLORS[i], axline=True)
            tools.plot_line(plt, w[1][0]+0.01*i, w[1][1]+0.01*i, color=tools.COLORS[i], axline=True)

        plt.show()


    def test_wedge_in_std(self):
        p = tools.rand_points(1)[0]
        s = Sectors(Polar.to_line(p), polarO)
        o = omega.Omega(vertices)
        tools.plot_congruent(plt, o.vertices)
        #plt.axline(*o..edge(0))
        tools.scatter(plt, [p, Polar.to_point(s.wedges[0][0])], colors=[tools.COLORS[0]])
        tangents = o..point_tangents(Polar.to_point(s.wedges[0][0]))
        #print(Polar.to_point(s.wedges[0][0]), tangents)
        #plt.axline(Polar.to_point(s.wedges[0][0]), tangents[0].v)
        #plt.axline(Polar.to_point(s.wedges[0][0]), tangents[1].v)

        for xxx in s.wedgetangentidxs:
            print(xxx, o..edge(s.wedgetangentidxs[xxx]))
        # The vertices of polar omega correspond to the vertices i and i+1
        Aidx = 0
        Bidx = s.wedgetangentidxs[0]
        As, Bs = o..edge(Aidx), o..edge(Bidx)
        line0 = [geometry.hdist_to_euc(p, As[0], Bs[0], 1), geometry.hdist_to_euc(p, As[1], Bs[1], 1)]
        line1 = [geometry.hdist_to_euc(p, As[0], Bs[0], -1), geometry.hdist_to_euc(p, As[1], Bs[1], -1)]
        tools.scatter(plt, line0, colors=['magenta'])
        tools.scatter(plt, line1, colors=['lime'])
        plt.axline(*line0, color='gray')
        plt.axline(*line1, color='gray')

        plt.axline(*o..edge(s.wedgetangentidxs[0]))
        plt.show()

    def test_wedges_with_first_distance_point(self):
        o = omega.Omega(vertices)
        p = tools.rand_points(1)[0]
        s = Sectors(Polar.to_line(p), polarO)
        wpt0, wpt1 = s.wedges[0][0]
        #l = line.Line(p, Polar.to_point(wedgeline), o)
        #A, B = l.get_boundary_intersections()
        Aidx = 0
        Bidx = s.wedgetangentidxs[0]
        As, Bs = o..edge(Aidx), o.edge(Bidx)
        #line0 = [geometry.hdist_to_euc(p, As[0], Bs[0], 1), geometry.hdist_to_euc(p, As[1], Bs[1], 1)]
        #line1 = [geometry.hdist_to_euc(p, As[0], Bs[0], -1), geometry.hdist_to_euc(p, As[1], Bs[1], -1)]
        pinter = euclidean.intersect((wpt0, wpt1), s.p)
        rsamples = [geometry.hdist_to_euc(pinter, wpt0, wpt1, 1), geometry.hdist_to_euc(pinter, wpt0, wpt1, -1)]# , geometry.hdist_to_euc(pinter, wpt1, wpt0, 5), geometry.hdist_to_euc(pinter, wpt1, wpt0, -5)]
        #vanishingpt = euclidean.intersect(*l.get_boundaries())
        #polarwedgept = Polar.to_point(stdwedgept, p)

        Polar.plot_congruent(plt, polarO.vertices)
        #tools.scatter(plt, [Polar.to_point(line0), Polar.to_point(line1)], colors=['black', 'black'])
        #print("POLAR POINTS: ", Polar.to_point(line0), Polar.to_point(line1))
        for i, w in enumerate(s.wedges):
            tools.plot_line(plt, w[0][0]+0.01*i, w[0][1]+0.01*i, color=tools.COLORS[i], axline=True)
            tools.plot_line(plt, w[1][0]+0.01*i, w[1][1]+0.01*i, color=tools.COLORS[i], axline=True)

        tools.scatter(plt, rsamples, colors=['lime', 'magenta', 'black', 'gray'])
        tools.scatter(plt, s.intersections)

        #tools.plot_line(plt, *Polar.to_line(stdwedgept), color='black')
        plt.xlim(-3, 3)
        plt.ylim(-7, 4)
        plt.show()







    def test_what_polar_points_on_wedge_look_like_in_std(self):
        o = omega.Omega(vertices)
        p = tools.rand_points(1)[0]
        s = Sectors(Polar.to_line(p), polarO)
        s.pseudohyperbola(1, plt)


