import unittest

import hilbert.omega
from hilbert.polygon import Polygon
from hilbert.omega import Omega
from hilbert import geometry, line
from matplotlib import pyplot as plt
from numpy.random import dirichlet
import numpy as np
from misc import tools

class MyTestCase(unittest.TestCase):

    def test_contains(self):
        p = np.array([0.25, 0.25])
        o = Omega()
        ball = o.hilbert_ball_around_point(p, 2)
        i = 8
        while i > 0:
            i -= 1
            q = dirichlet([1, 1, 1], 1)[0][:2]
            tools.plot_congruent(plt, o.vertices, color='gray')
            tools.plot_congruent(plt, ball.vertices, color='orange')
            tools.scatter(plt, [q])
            tools.annotate(plt, 'qp', [q, p])
            plt.title(f"{q} inside? {ball.contains(q)}")
            plt.show()


    def test_maxwards_tangent_linear(self):
        p = np.array([0.25, 0.25])
        o = Omega()
        ball = o.hilbert_ball_around_point(p, 1)
        i = 8
        while i > 0:
            i -= 1
            q = dirichlet([1, 1, 1], 1)[0][:2]
            while ball.contains(q):
                q = dirichlet([1, 1, 1], 1)[0][:2]
            try:
                t = ball.point_tangent_linear(q, maxwards=True)
            except Exception as e:
                print(e, e.__doc__)
                t = None
            tools.plot_line(plt, q, t)
            tools.plot_congruent(plt, o.vertices, color='gray')
            tools.plot_congruent(plt, ball.vertices, color='orange')
            tools.scatter(plt, [q, t])
            tools.annotate(plt, 'qpt', [q, p, t])
            plt.title(f"Q: {q}, T: {t}")
            plt.show()


    def test_minwards_tangent_linear(self):
        p = np.array([0.25, 0.25])
        o = Omega()
        ball = o.hilbert_ball_around_point(p, 1)
        i = 8
        while i > 0:
            i -= 1
            q = dirichlet([1, 1, 1], 1)[0][:2]
            while ball.contains(q):
                q = dirichlet([1, 1, 1], 1)[0][:2]
            try:
                t = ball.point_tangent_linear(q, maxwards=False)
            except Exception as e:
                print(e, e.__doc__)
                t = None
            tools.plot_line(plt, q, t)
            tools.plot_congruent(plt, o.vertices, color='gray')
            tools.plot_congruent(plt, ball.vertices, color='orange')
            tools.scatter(plt, [q, t])
            tools.annotate(plt, 'qpt', [q, p, t])
            plt.title(f"Q: {q}, T: {t}")
            plt.show()

    def test_maxwards_tangent_binary(self):
        p = np.array([0.25, 0.25])
        o = Omega()
        ball = o.hilbert_ball_around_point(p, 1)
        i = 8
        while i > 0:
            i -= 1
            q = dirichlet([1, 1, 1], 1)[0][:2]
            while ball.contains(q):
                q = dirichlet([1, 1, 1], 1)[0][:2]
            try:
                t = ball.point_tangent_binary(q, maxwards=True)
                plt.title(f"Broke")
            except Exception as e:
                print(e, e.__doc__)
                t = ball.point_tangent_linear(q, maxwards=True)
            tools.plot_line(plt, q, t)
            tools.plot_congruent(plt, o.vertices, color='gray')
            tools.plot_congruent(plt, ball.vertices, color='orange')
            tools.scatter(plt, [q, t])
            tools.annotate(plt, 'qpt', [q, p, t])
            tools.annotate(plt, range(len(ball.vertices)), ball.vertices)
            plt.show()


if __name__ == '__main__':
    unittest.main()
