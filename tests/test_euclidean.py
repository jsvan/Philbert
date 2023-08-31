from unittest import TestCase
from src.hilbert import omega, line
from src.misc import tools, euclidean
import numpy as np
from matplotlib import pyplot as plt

class Test_Point_Below_Line(TestCase):

    def test_horizontal(self):
        vs = [[0, 0], [0, 1], [1, 0]]
        o = omega.Omega(vs)
        l = line.Line(np.array([0.1, 0.4]), np.array([0.3, 0.4]), o)
        self.assertTrue(euclidean.point_below_line(vs[0], l.l))
        self.assertFalse(euclidean.point_below_line(vs[1], l.l))
        self.assertTrue(euclidean.point_below_line(vs[2], l.l))


    def test_vertical(self):
        """
        1 .
        | . \
        0-.--2
        """
        vs = [[0, 0], [0, 1], [1, 0]]
        o = omega.Omega(vs)
        l = line.Line(np.array([0.4, 0.1]), np.array([0.4, 0.3]), o)
        self.assertFalse(euclidean.point_below_line(vs[0], l.l))
        self.assertFalse(euclidean.point_below_line(vs[1], l.l))
        self.assertTrue(euclidean.point_below_line(vs[2], l.l))


class Test_Intersect(TestCase):

    def test_a(self):
        sa = np.array([1.4, 0.0])
        sb = np.array([-.7, -.2])
        la, lb = np.array([0.1, 0.4]), np.array([0.3, 0.4])

        tools.plot_line(plt, sa, sb)
        plt.axline(la, lb)
        z = euclidean.intersect([sa, sb], [la, lb])
        tools.plot_line(plt, z, sa)
        plt.scatter([z[0], sa[0], sb[0], la[0], lb[0]], [z[1],sa[1], sb[1], la[1], lb[1]])
        plt.show()