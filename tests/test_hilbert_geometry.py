import unittest

import matplotlib.pyplot as plt

from hilbert import geometry
import numpy as np
from numpy import array as na
from hilbert import omega, line
from misc import tools, euclidean
import test_base

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


    def test_vis_dist(self):
        p, q = tools.rand_points(2)
        o = test_base.primal_omega_simplex
        l = line.Line(p, q, o)
        A, B = l.get_boundary_intersections()
        xneg, xpos, yneg, ypos = [], [], [], []
        derivpos, derivneg = [], []
        #accelpos, accelneg = [], []
        minderiv,minidx = 100, 0

        for d in range(100):
            hdist = d / 10
            where_plus = geometry.hdist_to_euc(p, A, B, hdist)
            where_minus = geometry.hdist_to_euc(p, A, B, -hdist)
            xpos.append(np.linalg.norm(where_plus.v - p.v))
            xneg.append(- np.linalg.norm(p.v - where_minus.v))
            ypos.append(hdist)
            yneg.append(-hdist)

            derivpos.append(geometry.derivative(where_plus, A, B))
            derivneg.append(geometry.derivative(where_minus, A, B))
            #accelpos.append(geometry.acceleration(where_plus, A, B))
            #accelneg.append(geometry.acceleration(where_minus, A, B))
            if derivpos[-1] < minderiv: minderiv = derivpos[-1]; minidx = len(derivpos)-1
            if derivneg[-1] < minderiv: minderiv = derivneg[-1]; minidx = -(len(derivpos)-1)


        for i, ii in tools.i_ii(len(xpos), connect_ends=False):
            pa = (xpos[i], ypos[i])
            pb = (xpos[ii], ypos[ii])
            plt.scatter(xpos[i], euclidean.slope((pa, pb)), color='orange')
            #plt.scatter(xpos[i], accelpos[i], color='green', s=4)
            plt.scatter(xpos[i], derivpos[i], color='magenta', s=4, zorder=1000)

        for i, ii in tools.i_ii(len(xneg), connect_ends=False):
            pa = (xneg[i], yneg[i])
            pb = (xneg[ii], yneg[ii])
            plt.scatter(xneg[i], euclidean.slope((pa, pb)), color='orange')
            #plt.scatter(xneg[i], accelneg[i], color='green', s=4)
            plt.scatter(xneg[i], derivneg[i], color='magenta', s=4, zorder=1000)

        print(f"min derivative is {minderiv} \ndistance A to B is {np.linalg.norm(A-B)}")
        ycept = ypos[minidx] if minidx >=0 else yneg[-minidx]
        xcept = xpos[minidx] if minidx >= 0 else xneg[-minidx]
        plt.axline((xcept, ycept), (0.1+xcept, (0.1*minderiv) + ycept), color='magenta', zorder=1001)
        plt.scatter(xneg+xpos, yneg+ypos, zorder=1000)
        plt.ylim(-10, 10)
        plt.title("Distances from p along cord running through omega" +
                  f"\nminimum derivative is {str(minderiv)[:7]} \ndistance AB is {str(np.linalg.norm(A-B))[:7]}")
        plt.xlabel("Linear distance from p")
        plt.ylabel("Hilbert distance from p")
        plt.show()

    def test_dist_march_simplex(self):
        d = geometry.nielson_dist
        e = 0.0001
        eP = euclidean.Point
        p = euclidean.Point((1-e, e))
        qs = [eP((e, 1-e)), eP((.1, .9)), eP((.2, .8)), eP((.3, .7)), eP((.4, .6))
              , eP((.5, .5)), eP((.6, .4)), eP((.7, .3)), eP((.8, .2)), eP((.9, .1))
              , eP((.95, .05)), eP((.98, .02)), eP((.99, .01)), eP((1-e, e))]
        for q in qs:
            plt.scatter([np.linalg.norm(p.v-q.v)], [d(p, q)])
            print(q, d(p, q))
        plt.show()

if __name__ == '__main__':
    unittest.main()
