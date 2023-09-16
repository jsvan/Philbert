import unittest
from misc import tools, euclidean
from hilbert import geometry, line, omega, polygon
from matplotlib import pyplot as plt


class MyTestCase(unittest.TestCase):
    def test_something(self):
        o = omega.Omega(vertices=[(1.4, 0.0),
                            (1.1, 0.8),
                            (0.6, 1.1),
                            (0.0, 1.2),
                            (-.5, 0.8),
                            (-.8, 0.3),
                            (-.7, -.2)])
        best_dividing_line = line.Line(*tools.rand_points(2), o)
        bdl_ball = best_dividing_line.hilbert_ball_about_line(1)
        tools.plot_congruent(plt, o.vertices, color='gray')
        tools.plot_congruent(plt, [x for x in bdl_ball.vertices if euclidean.point_below_line(x, best_dividing_line.l)],
                             color="blue", connect_ends=False)
        plt.axline(best_dividing_line.p, best_dividing_line.q, color='orange')
        plt.show()











if __name__ == '__main__':
    unittest.main()
