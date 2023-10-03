import unittest
from misc import tools, euclidean
from hilbert import geometry, line, omega, polygon
from matplotlib import pyplot as plt


class MyTestCase(unittest.TestCase):
    def test_visualize_line_balls(self):
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


    def test_two_vs_three_point_best_dividing_lines(self):
        o = omega.Omega(vertices=[(1.4, 0.0),
                                  (1.1, 0.8),
                                  (0.6, 1.1),
                                  (0.0, 1.2),
                                  (-.5, 0.8),
                                  (-.8, 0.3),
                                  (-.7, -.2)])
        p, q = tools.rand_points(2)
        l = line.Line(p, q, o)
        ball = l.hilbert_ball_about_line(1)
        above = [x for x in ball.vertices if euclidean.point_below_line(x, l.l)]
        below = list(reversed([x for x in ball.vertices if not euclidean.point_below_line(x, l.l)]))
        b1, dist = euclidean.uniform_sample_from_line_segments(above)
        r, dist = euclidean.uniform_sample_from_line_segments(below, dist)
        b2, dist = euclidean.uniform_sample_from_line_segments(above, dist)
        for p in [b1, b2, r]:
            for v in o.vertices:
                plt.axline(p, v, color='gainsboro')
        tools.plot_congruent(plt, o.vertices, color='gray')
        tools.plot_congruent(plt, ball.vertices, color='lightblue')
        plt.axline(*l.l)
        try:
            lpr = line.Line(b1, r, o).get_best_dividing_line()
            lqr = line.Line(b2, r, o).get_best_dividing_line()
            pintersection = euclidean.intersect(l.l, lpr)
            qintersection = euclidean.intersect(l.l, lqr)
            intersection = euclidean.intersect(lpr, lqr)
            print('p to qline: ', line.Line(b1, qintersection, o).get_hdist())
            print('q to pline: ', line.Line(b2, pintersection, o).get_hdist())
            print('\np to pline intersect: ', line.Line(b1, pintersection, o).get_hdist())
            print('q to qline intersect: ', line.Line(b2, qintersection, o).get_hdist())

            print('\nintersect to pintersect: ', line.Line(intersection, pintersection, o).get_hdist())
            print('intersect to qintersect: ', line.Line(intersection, qintersection, o).get_hdist())
            #ballp = o.hilbert_ball_around_point(b1, 1)
            #ballq = o.hilbert_ball_around_point(b2, 1)
            #ballr = o.hilbert_ball_around_point(r, 1)
            print('\nTriangle stuff:\np: ', line.Line(*lpr, o).get_hdist() / line.Line(intersection, qintersection, o).get_hdist())
            print('q: ', line.Line(*lqr, o).get_hdist() / line.Line(intersection, pintersection, o).get_hdist())
            print('pq: ', line.Line(b1, b2, o).get_hdist() / line.Line(pintersection, qintersection, o).get_hdist())


            print('\nDistance of points from intersection:',
                  line.Line(b1, intersection, o).get_hdist(),
                  line.Line(b2, intersection, o).get_hdist(),
                  line.Line(r, intersection, o).get_hdist())

            plt.axline(*lpr, color="green")
            plt.axline(*lqr, color="green")
            #tools.plot_congruent(plt, ballp.vertices, color="gray")
            #tools.plot_congruent(plt, ballq.vertices, color="gray")
            #tools.plot_congruent(plt, ballr.vertices, color="gray")

            bs_line_extremes = line.Line(*lpr, o).get_boundary_intersections()
            bs = tools.BinarySearcher(0, 1, discrete=False)

            ppp = None
            # Doing the p line, so search against q and r
            while bs.has_next():
                ppp = euclidean.get_point_on_line(*bs_line_extremes, bs.next())
                qdist, rdist = line.Line(ppp, b2, o).get_hdist(), line.Line(ppp, r, o).get_hdist()
                print('distances:', qdist, rdist)
                diff = qdist - rdist
                if diff < 0.0001:
                    break
                bs.feedback(higher=diff > 0)

            tools.scatter(plt, [line.Line(b1, r, o).get_midpoint(), line.Line(b2, r, o).get_midpoint(), ppp])
            tools.annotate(plt, 'x', [ppp])

        except Exception as e:
            print(e, e.__doc__)

        tools.scatter(plt, [b1, b2, r])
        tools.annotate(plt, 'pqr', [b1, b2, r])
        plt.show()


if __name__ == '__main__':
    unittest.main()
