import unittest
from hilbert import geometry
from numpy import array as na
from hilbert import omega, line
from misc import tools

class TestDist(unittest.TestCase):
    def test_vert(self):
        print(geometry.hdist_to_euc(na([1, 3]), na([1,5]), na([1,1]), 0.5))
    def test_horiz(self):
        print(geometry.hdist_to_euc(na([3, 1]), na([5,1]), na([1,1]), 0.5))


    def test_dist(self):
        """
        You can mix up A and B, dist is still the same.
        """
        o = omega.Omega(vertices=[(1.4, -0.1),
                                  (1.1, 0.8),
                                  (0.6, 1.1),
                                  (-0.2, 1.2),
                                  (-.5, 0.8),
                                  (-.8, 0.3),
                                  (-.7, -.2)])

        for _ in range(1):
            p, q = tools.rand_points(2)
            l = line.Line(p, q, o)
            A, B = l.get_boundary_intersections()
            print(f"Standard dist: {l.get_hdist()}, opposite: {geometry.dist(p, q, B, A)}, other opposite: {geometry.dist(q, p, A, B)}")



if __name__ == '__main__':
    unittest.main()
