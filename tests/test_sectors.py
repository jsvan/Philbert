import unittest
import test_base
from hilbert import geometry, omega, polygon, line
from matplotlib import pyplot as plt
from misc import tools, euclidean
from misc.dual import Polar, Sectors
import numpy as np

class MyTestCase(unittest.TestCase):
    def test_wedges(self):
        onaught = test_base.omega_for_dual
        p = tools.rand_points(1)[0]
        pnaught = Polar.to_line(p)
        sectors = Sectors(pnaught, onaught)

        tools.plot_congruent(plt, onaught.vertices, color='gray')
        for i, s in enumerate(sectors.sektoren):
            l = s.wedge.line_a
            tools.plot_line(plt, *l, axline=True, color=tools.COLORS[s.o_i])
            tools.scatter(plt, [s.wedge.aXp], colors=[tools.COLORS[s.o_i]])
            tools.annotate(plt, [i], [s.wedge.aXp])
            tools.plot_line(plt, *s.wedge.line_b, axline=True, color=tools.COLORS[s.o_i])

        tools.plot_line(plt, *pnaught, color='black', axline=True)
        plt.show()

    def test_pseudohyperbola(self):
        onaught = test_base.omega_for_dual
        p = tools.rand_points(1)[0]
        pnaught = Polar.to_line(p)
        sectors = Sectors(pnaught, onaught)
        pos_hyper, neg_hyper = sectors.pseudohyperbola(1)

        tools.plot_congruent(plt, onaught.vertices, color='gray')
        for i, s in enumerate(sectors.sektoren):
            l = s.wedge.line_a
            tools.plot_line(plt, *l, axline=True, color=tools.COLORS[s.o_i])
            #tools.scatter(plt, [s.wedge.aXp], colors=[tools.COLORS[s.o_i]])
            #tools.annotate(plt, [i], [s.wedge.aXp])
            #tools.plot_line(plt, *s.wedge.line_b, axline=True, color=tools.COLORS[s.o_i])
        tools.plot_congruent(plt, pos_hyper, color='blue', connect_ends=False)
        tools.plot_congruent(plt, neg_hyper, color='orange', connect_ends=False)

        tools.scatter(plt, pos_hyper, colors=['blue'])
        tools.scatter(plt, neg_hyper, colors=['orange'])

        nodenames = dict()
        name = 0
        print("Negatives")
        for v in neg_hyper:
            k = tuple(v.v)
            if k not in nodenames:
                nodenames[k] = name
                print(name, v)
                name += 1
            else:
                print(nodenames[k], v)

        print("Positives")
        for v in pos_hyper:
            k = tuple(v.v)
            if k not in nodenames:
                nodenames[k] = name
                print(name, v)
                name += 1
            else:
                print(nodenames[k], v)


        tools.plot_line(plt, *pnaught, color='black', axline=True)
        plt.xlim(-10, 10)
        plt.ylim(-10, 10)
        plt.show()










if __name__ == '__main__':
    unittest.main()
