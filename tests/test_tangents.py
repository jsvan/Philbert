import unittest

from hilbert.polygon import Polygon
from hilbert.omega import Omega
from hilbert import geometry, line
from matplotlib import pyplot as plt
from numpy.random import dirichlet
import numpy as np
from misc import tools

class TestTangents(unittest.TestCase):

    def test_hilbertball_tangents(self, p, q, fix_q=False):
        o = Omega()
        b1 = o.hilbert_ball_around_point(p, 1)
        l = line.Line(p, q, o)
        if fix_q and l.get_hdist() <= 2:
            q = geometry.hdist_to_euc(p, l.A, l.B, 2.01)
            #q = dirichlet([1, 1, 1], 1)[0][:2]

        print(p, q)
        b2 = o.hilbert_ball_around_point(q, 1)
        try:
            tangents = o.all_tangent_lines(b1, b2)
        except Exception as e:
            print("Failed", e, e.__doc__)
            tangents = {"broken":(p, q)}

        for k, v in tangents.items():
            tools.plot_congruent(plt, o.vertices, color='gray')
            tools.plot_congruent(plt, b1.vertices, color='blue')
            tools.plot_congruent(plt, b2.vertices, color='green')
            tools.annotate(plt, 'pq', [p, q])
            plt.axline(v[0], v[1])
            plt.title(k)
            plt.show()


    def test_broken_a(self):
        p, q = [0.71075363, 0.20640696], [0.29698898, 0.46192805]
        self.test_hilbertball_tangents(p, q, True)

    def test_broken_b(self):
        p, q = [0.22620338, 0.29044416], [0.27322547, 0.30064652]
        self.test_hilbertball_tangents(p, q, True)

    def test_random(self):
        p = dirichlet([1, 1, 1], 1)[0][:2]
        q = dirichlet([1, 1, 1], 1)[0][:2]
        self.test_hilbertball_tangents(p, q, True)















if __name__ == '__main__':
    unittest.main()
