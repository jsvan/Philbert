import unittest
from hilbert import omega, line, geometry
from misc import tools
import numpy as np
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
        p1, p2 = np.array([0.5, 0.5]), np.array([0.6, 0.5])
        l = line.Line(p1, p2, o)
        q = tools.rand_points(1)[0]
        point, dist, idx = l.nearest_point(q, plt)



        print("Point, dist, idx", point,dist, idx)


        intersect_line = line.Line(q, point, o)
        A, B = intersect_line.get_boundary_intersections()
        tools.plot_congruent(plt, o.vertices, color='gray')
        plt.axline(p1, p2)
        ppp = geometry.hdist_to_euc(q,A,B,dist)
        tools.annotate(plt, 'xp', [ppp, q])
        tools.scatter(plt, [ppp, geometry.hdist_to_euc(q,A,B,-dist)])
        tools.plot_line(plt, q, point)
        plt.show()



if __name__ == '__main__':
    unittest.main()
