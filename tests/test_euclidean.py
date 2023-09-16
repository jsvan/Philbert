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

class Test_Point_On_Line(TestCase):

    def test_online(self):
        p = np.array([0.2, 0.4])
        l = np.array([0.1, 0.4]), np.array([0.3, 0.4])
        self.assertTrue(euclidean.point_on_line(p, l, segment=True))

    def test_toofar(self):
        p = np.array([0.4, 0.4])
        l = np.array([0.1, 0.4]), np.array([0.3, 0.4])
        self.assertFalse(euclidean.point_on_line(p, l, segment=True))

    def test_offline(self):
        p = np.array([0.2, 0.5])
        l = np.array([0.1, 0.4]), np.array([0.3, 0.4])
        self.assertFalse(euclidean.point_on_line(p, l, segment=True))

    def test_toofar_inagoodway(self):
        p = np.array([0.4, 0.4])
        l = np.array([0.1, 0.4]), np.array([0.3, 0.4])
        self.assertTrue(euclidean.point_on_line(p, l, segment=False))

    def test_vert(self):
        p = np.array([0., 0.4])
        l = np.array([0, 0]), np.array([0, 1])
        self.assertTrue(euclidean.point_on_line(p, l, segment=True))

    def test_diag(self):
        p = np.array([0.6, 0.4])
        l = np.array([0, 1]), np.array([1, 0])
        self.assertTrue(euclidean.point_on_line(p, l, segment=True))

    def test_get_point_on_line(self):
        p = np.array([0.8, 0.2])
        q = np.array([0.8, 0.4])
        dist = 1
        r = euclidean.get_point_on_line(p, q, dist)
        tools.annotate(plt, 'pqr', [p, q, r])
        tools.scatter(plt, [p, q, r])
        tools.plot_congruent(plt, [(0,0), (0,1), (1,1), (1, 0)], color='gray')
        plt.show()


    def test_uniform_segment_sampling(self):
        o = omega.Omega(vertices=[(1.4, 0.0),
                                    (1.1, 0.8),
                                    (0.6, 1.1),
                                    (0.0, 1.2),
                                    (-.5, 0.8),
                                    (-.8, 0.3),
                                    (-.7, -.2)])
        point, dist = euclidean.uniform_sample_from_line_segments(o.vertices)
        print(point, dist)

