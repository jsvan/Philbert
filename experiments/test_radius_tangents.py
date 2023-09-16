import unittest

from hilbert.polygon import Polygon
from hilbert.omega import Omega
from hilbert import geometry, line
from matplotlib import pyplot as plt
from numpy.random import dirichlet
import numpy as np
from misc import tools


class TangentsAsRadiusGrows(unittest.TestCase):
    """
    RESULTS OF TESTS:
    SPOKES CHANGE AS RADIUSES GROW, REGARDLESS OF WHETHER THEY INTERSECT OR NOT.
    """
    def test_random(self):
        p = dirichlet([1, 1, 1], 1)[0][:2]
        q = dirichlet([1, 1, 1], 1)[0][:2]
        self.test_tangents_v_radius_weird(p, q)

    def test_tangents_v_radius_simplex(self, p, q):
        o = Omega()

        print(p, q)
        for radius in [0.1,0.5, 1,1.4,1.8,2,2.2,2.4,2.5,22.6,2.7,2.8,2.9,3]:
            b1 = o.hilbert_ball_around_point(p, radius)
            b2 = o.hilbert_ball_around_point(q, radius)
            tangents = o.all_tangent_lines(b1, b2)
            print(tangents)
            tools.plot_congruent(plt, o.vertices, color='gray')
            tools.plot_congruent(plt, b1.vertices, color='black')
            tools.plot_congruent(plt, b2.vertices, color='black')
            tools.annotate(plt, 'pq', [p, q])
            colors = ["red", "orange", "blue", "green"]
            for i, t in enumerate(tangents):
                if t[0]['v'] is not None and t[1]['v'] is not None:
                    print(i, t)
                    plt.axline(t[0]['v'], t[1]['v'], color=colors[i])

            annots = tools.annotate(plt, [t[0]['i'] for t in tangents], [t[0]['v'] for t in tangents])
            tools.annotate(plt, [t[1]['i'] for t in tangents], [t[1]['v'] for t in tangents], annots)

            #plt.title(' -- '.join([str(x) for x in [tangentA[0]['i'], tangentA[1]['i'], tangentB[0]['i'], tangentB[1]['i']]]))
            plt.show()

    def test_tangents_v_radius_weird(self, p, q):
        o = Omega(vertices=[(1.4, 0.0),
                              (1.1, 0.8),
                              (0.6, 1.1),
                              (0.0, 1.2),
                              (-.5, 0.8),
                              (-.8, 0.3),
                              (-.7, -.2)])

        print(p, q)
        for radius in [0.1,0.5, 1,1.4,1.8,2,2.2,2.4,2.5,22.6,2.7,2.8,2.9,3]:
            b1 = o.hilbert_ball_around_point(p, radius)
            b2 = o.hilbert_ball_around_point(q, radius)
            tangents = o.all_tangent_lines(b1, b2)
            print(tangents)
            tools.plot_congruent(plt, o.vertices, color='gray')
            tools.plot_congruent(plt, b1.vertices, color='black')
            tools.plot_congruent(plt, b2.vertices, color='black')
            tools.annotate(plt, 'pq', [p, q])
            colors = ["red", "orange", "blue", "green"]
            for i, t in enumerate(tangents):
                if t[0]['v'] is not None and t[1]['v'] is not None:
                    print(i, t)
                    plt.axline(t[0]['v'], t[1]['v'], color=colors[i])

            annots = tools.annotate(plt, [t[0]['i'] for t in tangents], [t[0]['v'] for t in tangents])
            tools.annotate(plt, [t[1]['i'] for t in tangents], [t[1]['v'] for t in tangents], annots)

            #plt.title(' -- '.join([str(x) for x in [tangentA[0]['i'], tangentA[1]['i'], tangentB[0]['i'], tangentB[1]['i']]]))
            plt.show()





class TangentsAndBalls(unittest.TestCase):

    def test_random_simplex(self):
        p = dirichlet([1, 1, 1], 1)[0][:2]
        q = dirichlet([1, 1, 1], 1)[0][:2]
        self.test_simplex(p, q)

    def test_random_weird(self):
        p = dirichlet([1, 1, 1], 1)[0][:2]
        q = dirichlet([1, 1, 1], 1)[0][:2]
        self.test_weird(p, q)

    def test_simplex(self, p, q):
        o = Omega()
        l = line.Line(p, q, o)
        print(p, q)
        for radius in [0.1,0.5, 1,1.4,1.8,2,2.2,2.4,2.5,22.6,2.7,2.8,2.9,3]:
            b1 = o.hilbert_ball_around_point(p, radius)
            b2 = o.hilbert_ball_around_point(q, radius)
            tangents = o.all_tangent_lines(b1, b2)
            print(tangents)
            tools.plot_congruent(plt, o.vertices, color='gray')
            tools.plot_congruent(plt, b1.vertices, color='black')
            tools.plot_congruent(plt, b2.vertices, color='black')
            tools.plot_congruent(plt, l.hilbert_ball_about_line(radius).vertices, color="pink")
            tools.annotate(plt, 'pq', [p, q])
            colors = ["red", "orange", "blue", "green"]
            for i, t in enumerate(tangents):
                if t[0]['v'] is not None and t[1]['v'] is not None:
                    print(i, t)
                    plt.axline(t[0]['v'], t[1]['v'], color=colors[i])

            annots = tools.annotate(plt, [t[0]['i'] for t in tangents], [t[0]['v'] for t in tangents])
            tools.annotate(plt, [t[1]['i'] for t in tangents], [t[1]['v'] for t in tangents], annots)

            #plt.title(' -- '.join([str(x) for x in [tangentA[0]['i'], tangentA[1]['i'], tangentB[0]['i'], tangentB[1]['i']]]))
            plt.show()

    def test_weird(self, p, q):
        o = Omega(vertices=[(1.4, 0.0),
                              (1.1, 0.8),
                              (0.6, 1.1),
                              (0.0, 1.2),
                              (-.5, 0.8),
                              (-.8, 0.3),
                              (-.7, -.2)])
        l = line.Line(p, q, o)
        print(p, q)
        for radius in [0.1,0.5, 1,1.4,1.8,2,2.2,2.4,2.5,22.6,2.7,2.8,2.9,3]:
            b1 = o.hilbert_ball_around_point(p, radius)
            b2 = o.hilbert_ball_around_point(q, radius)
            tangents = o.all_tangent_lines(b1, b2)
            print(tangents)
            tools.plot_congruent(plt, o.vertices, color='gray')
            tools.plot_congruent(plt, b1.vertices, color='black')
            tools.plot_congruent(plt, b2.vertices, color='black')
            tools.plot_congruent(plt, l.hilbert_ball_about_line(radius).vertices, color="pink")
            tools.annotate(plt, 'pq', [p, q])
            colors = ["red", "orange", "blue", "green"]
            for i, t in enumerate(tangents):
                if t[0]['v'] is not None and t[1]['v'] is not None:
                    print(i, t)
                    plt.axline(t[0]['v'], t[1]['v'], color=colors[i])

            annots = tools.annotate(plt, [t[0]['i'] for t in tangents], [t[0]['v'] for t in tangents])
            tools.annotate(plt, [t[1]['i'] for t in tangents], [t[1]['v'] for t in tangents], annots)

            #plt.title(' -- '.join([str(x) for x in [tangentA[0]['i'], tangentA[1]['i'], tangentB[0]['i'], tangentB[1]['i']]]))
            plt.show()



if __name__ == '__main__':
    unittest.main()
