import unittest
import test_base
from matplotlib import pyplot as plt
from misc import tools, euclidean
from hilbert import polygon, line
from misc.dual import Polar, Sectors

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


    def test_three_pseudohyperbola(self):
        onaught = test_base.omega_for_dual
        naughts = list(map(Polar.to_line, tools.rand_points(3)))
        sectors = list(map(lambda x: Sectors(x, onaught), naughts))

        for ii in range(1, 400):
            print('\r', ii, '/ 400', flush=True)
            radius = ii / 200
            pseudohyperbolas = map(lambda x: x.pseudohyperbola(radius), sectors)
            tools.plot_congruent(plt, onaught.vertices, color='gray')
            for i, hy in enumerate(pseudohyperbolas):
                tools.scatter(plt, hy[0], colors=[tools.COLORS[i]])
                tools.scatter(plt, hy[1], colors=[tools.COLORS[i]])
                tools.plot_line(plt, *naughts[i], color=tools.COLORS[i], axline=True)
            plt.xlim(-25, 25)
            plt.ylim(-25, 25)
            #tools.save_movie_frame(plt, ii)
            plt.show()
        #tools.finish_movie("growing_hyperboloids")


    def test_three_hyperbola_intersections(self):
        import time
        lasttime = time.time()
        onaught = test_base.omega_for_dual
        naughts = [Polar.to_line(x) for x in tools.rand_points(3)]
        sectors = [Sectors(x, onaught) for x in naughts]
        #intersections = []
        halvedpolygons = [polygon.HalvedPolygon(None, None), polygon.HalvedPolygon(None, None), polygon.HalvedPolygon(None, None)]
        for ii in range(1, 400):
            print('\r', ii, '/ 300', ' time of',time.time()-lasttime, end='', flush=True)
            lasttime = time.time()
            radius = ii / 150
            pseudohyperbolas = [x.pseudohyperbola(radius, halvedpolygons[i]) for i, x in enumerate(sectors)]
            intersections, newinters = pseudohyperbolas[0].intersections(pseudohyperbolas[1])
            tools.plot_congruent(plt, onaught.vertices, color='gray')
            for i, hy in enumerate(pseudohyperbolas):
                hy.plot_congruent(plt, color=tools.COLORS[i])
                tools.plot_line(plt, *naughts[i], color=tools.COLORS[i], axline=True)
            #tools.scatter_every_x(plt, intersections, x=ii, colors=['pink'], size=1)
            tools.scatter(plt, intersections, colors=['magenta'], size=1)
            plt.xlim(-25, 25)
            plt.ylim(-25, 25)
            tools.save_movie_frame(plt, ii)
        tools.finish_movie("intersecting_hyperboloids")




    def test_three_wedges(self):
        onaught = test_base.omega_for_dual
        j, g = tools.rand_points(2)
        mainline = line.Line(j, g, test_base.primal_omega_weird)
        ball = mainline.hilbert_ball_about_line(1.1)
        p, q, r = ball.sample_from_ball_border(3, mainline)

        pnaught, qnaught, rnaught = Polar.to_line(p), Polar.to_line(q), Polar.to_line(r)
        psectors = Sectors(pnaught, onaught)
        qsectors = Sectors(qnaught, onaught)
        rsectors = Sectors(rnaught, onaught)

        for i, sectors in enumerate([psectors,qsectors,rsectors]):
            for s in sectors.sektoren:
                tools.plot_line(plt, *s.wedge.line_a, axline=True, color=tools.COLORS[i])
                tools.plot_line(plt, *s.wedge.line_b, axline=True, color=tools.COLORS[i], linewidth=3-i)
        for v in onaught.vertices:
            for vv in onaught.vertices:
                tools.plot_line(plt, v, vv, color='white', linewidth=3.1)

        tools.plot_congruent(plt, onaught.vertices, color='gray')
        tools.scatter(plt, [Polar.to_point(mainline.l)], colors=['magenta'], size=5)
        plt.xlim(-16, 16)
        plt.ylim(-16, 16)
        plt.show()


if __name__ == '__main__':
    unittest.main()
