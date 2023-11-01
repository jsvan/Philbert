import unittest

from hilbert.omega import Omega
from misc.euclidean import Point
from matplotlib import pyplot as plt
from misc import tools
from misc.dual import Polar, Sectors
import test_base

class MyTestCase(unittest.TestCase):

    def test_contains(self):
        p = Point([0.25, 0.25])
        o = Omega()
        ball = o.hilbert_ball_around_point(p, 2)
        i = 8
        while i > 0:
            i -= 1
            q = tools.rand_points(1)[0]
            tools.plot_congruent(plt, o.vertices, color='gray')
            tools.plot_congruent(plt, ball.vertices, color='orange')
            tools.scatter(plt, [q])
            tools.annotate(plt, 'qp', [q, p])
            plt.title(f"{q} inside? {ball.contains(q)}")
            plt.show()


    def test_maxwards_tangent_linear(self):
        p = Point([0.25, 0.25])
        o = Omega()
        ball = o.hilbert_ball_around_point(p, 1)
        i = 8
        while i > 0:
            i -= 1
            q = tools.rand_points(1)[0]
            while ball.contains(q):
                q = tools.rand_points(1)[0]
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
        p = Point([0.25, 0.25])
        o = Omega()
        ball = o.hilbert_ball_around_point(p, 1)
        i = 8
        while i > 0:
            i -= 1
            q = tools.rand_points(1)[0]
            while ball.contains(q):
                q = tools.rand_points(1)[0]
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
        p = Point([0.25, 0.25])
        o = Omega()
        ball = o.hilbert_ball_around_point(p, 1)
        i = 8
        while i > 0:
            i -= 1
            q = tools.rand_points(1)[0]
            while ball.contains(q):
                q = tools.rand_points(1)[0]
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

    def test_half_polygon_sorting(self):
        onaught = test_base.omega_for_dual
        naughts = [Polar.to_line(x) for x in tools.rand_points(1)]
        sectors = [Sectors(x, onaught) for x in naughts]
        # intersections = []
        for radius in [0.3, 0.5, 0.8, 1.2, 1.8, 2.5]:
            pseudohyperbolas = [x.pseudohyperbola(radius) for x in sectors]
            for i, hy in enumerate(pseudohyperbolas):
                hy.plot_congruent(plt, color=tools.COLORS[i])
                tools.plot_line(plt, *naughts[i], color=tools.COLORS[i], axline=True)
                tools.annotate(plt, range(len(hy.lowerhalf)),hy.lowerhalf)
                tools.annotate(plt, range(len(hy.upperhalf)),hy.upperhalf)

            tools.plot_congruent(plt, onaught.vertices, color='gray')
            plt.show()


    def test_gridspace(self):
        onaught = test_base.omega_for_dual
        onaught.gridspace(0.5)
if __name__ == '__main__':
    unittest.main()
