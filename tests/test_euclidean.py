from unittest import TestCase
from hilbert import omega, line
from misc import tools, euclidean
import numpy as np
from matplotlib import pyplot as plt

class Test_Point_Below_Line(TestCase):

    def test_orient(self):
        while True:
            points = tools.rand_points(3)
            tools.scatter(plt, points)
            width = 0.0005
            for i, ii in tools.i_ii(3):
                plt.arrow(*points[i], *(points[ii] - points[i] ), width=width, head_width=20*width, length_includes_head=True)
            tools.annotate(plt, 'pqr', points)
            plt.title(euclidean.dirname(euclidean.orient(*points)))
            plt.show()
            plt.clf()

    def test_line_intersect(self):
        while True:
            points = tools.rand_points(4)
            tools.scatter(plt, points)
            tools.plot_line(plt, points[0], points[1])
            tools.plot_line(plt, points[2], points[3])
            tools.annotate(plt, 'aabb', points)
            plt.title(str(euclidean.line_intersects_on_segments(points[0:2], points[2:])))
            plt.show()
            plt.clf()


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

    def test_same(self):
        p, q = tools.rand_points(2)
        l1 = euclidean.BasicLine(p, q)
        l2 = euclidean.BasicLine(
            euclidean.get_point_on_line(p,q,3),
            euclidean.get_point_on_line(p,q,1)
            )
        print(l1)
        print(l2)
        plt.axline(*l1)
        plt.axline(*l2)
        z = euclidean.intersect(l1, l2)
        print(z)
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


    def test_spin(self):
        p = np.array([1,1])
        for i in range(50):
            q = euclidean.circle_point(i, 50)
            plt.axline(p, q, color='orange')

        plt.show()

class Test_Line(TestCase):
    def test_avg_line(self):
        p, pp, q, qq = tools.rand_points(4)
        la = euclidean.BasicLine(p, pp)
        lb = euclidean.BasicLine(q, qq)
        avgl = euclidean.avg_line(la, lb)
        tools.plot_line(plt, *la, color='blue', axline=True)
        tools.plot_line(plt, *lb, color='orange', axline=True)
        tools.plot_line(plt, *avgl, color='green', axline=True)
        plt.show()
