import unittest
from hilbert import geometry, omega, polygon, line
from matplotlib import pyplot as plt
from misc import tools, euclidean
from misc.dual import Dual2
import numpy as np

class TestDual(unittest.TestCase):

    def test_omega(self):
        """
        passes
        """
        o = omega.Omega(vertices=[(1.4, -0.1),
                                  (1.1, 0.8),
                                  (0.6, 1.1),
                                  (-0.2, 1.2),
                                  (-.5, 0.8),
                                  (-.8, 0.3),
                                  (-.7, -.2)], offset=[2,2])

        plt.subplot(1, 2, 1)
        tools.plot_congruent(plt, o.vertices, color='red')
        plt.title("standard")
        plt.subplot(1, 2, 2)
        for i, v in enumerate(o.vertices):
            q = Dual2.to_line(v)
            print(i,v, q)
            tools.annotate(plt,[i], [q[0]])
            #tools.plot_line(plt, *q, color='blue')
            plt.axline(*q, color="blue")

        scatpoints = Dual2.v2v(o.vertices)
        print("Scat points: len", len(scatpoints), scatpoints)
        tools.scatter(plt, scatpoints)
        plt.title("Dual2 dual")
        plt.tight_layout()
        plt.show()


    def test_hilbert_ball(self):
        o = omega.Omega(vertices=[(1.4, -0.1),
                                  (1.1, 0.8),
                                  (0.6, 1.1),
                                  (-0.2, 1.2),
                                  (-.5, 0.8),
                                  (-.8, 0.3),
                                  (-.7, -.2)])
        lims = (-5, 5)
        plt.subplot(1, 2, 1)
        plt.title("standard")
        tools.plot_congruent(plt, o.vertices, color='red')
        p = tools.rand_points(1)[0]
        print("p", p)
        tools.scatter(plt, [p])
        ball = o.hilbert_ball_around_point(p, 1).vertices
        tools.plot_congruent(plt, ball, color="orange")
        l = np.array([0.2, 0.2]), np.array([0.5, 0.2])
        otherl = line.Line(*l, o)
        ints = otherl.get_boundary_intersections()
        plt.scatter(ints[0][0],  ints[0][1], color='black')
        plt.scatter(ints[1][0],  ints[1][1], color='green')

        plt.axline(*l, color='gray')

        plt.subplot(1, 2, 2, xlim=lims, ylim=lims)
        plt.title("Dual2 dual")
        """
        for i, v in enumerate(o.vertices):
            q = dual.Dual2.to_line(v)
            print(i, v, q)
            tools.annotate(plt, [i], [q[0]])
            plt.axline(*q, color="blue")
        """
        plt.scatter(*Dual2.to_point(l), color='gray')
        tools.plot_congruent(plt, Dual2.v2v(o.vertices), color='red', axline=False)
        tools.plot_congruent(plt, Dual2.v2v(ball), color='orange', axline=False)
        plt.axline(*Dual2.to_line(p))
        plt.axline( *Dual2.to_line(ints[0]), color='black')
        plt.axline( *Dual2.to_line(ints[1]), color='green')
        plt.tight_layout()
        plt.show()



    def test_nested_hilbert_ball(self):
        o = omega.Omega(vertices=[(1.4, -0.1),
                                  (1.1, 0.8),
                                  (0.6, 1.1),
                                  (-0.2, 1.2),
                                  (-.5, 0.8),
                                  (-.8, 0.3),
                                  (-.7, -.2)])

        plt.subplot(1, 2, 1)
        plt.title("standard")
        tools.plot_congruent(plt, o.vertices, color='red')
        p = tools.rand_points(1)[0]
        tools.scatter(plt, [p], colors=['lime'])
        balls = [o.hilbert_ball_around_point(p, n).vertices for n in [0.5, 1, 2, 4]]
        colors = ['blue', 'lightblue', 'yellow', 'orange']
        for i, ball in enumerate(reversed(balls)):
            tools.plot_congruent(plt, ball, color=colors[i])

        plt.subplot(1, 2, 2)
        plt.title("Dual2 dual")

        for i, ball in enumerate(reversed(balls)):
            tools.plot_congruent(plt, Dual2.v2v(ball), color=colors[i], axline=False)
        plt.axline(*Dual2.to_line(p), color='lime', zorder=100)
        tools.plot_congruent(plt, Dual2.v2v(o.vertices), color='red', axline=False)


        plt.tight_layout()
        plt.show()



    def test_line_and_3_points_and_their_hilbert_balls(self):
        o = omega.Omega(vertices=[(1.4, -0.1),
                                  (1.1, 0.8),
                                  (0.6, 1.1),
                                  (-0.2, 1.2),
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
        ball1 = o.hilbert_ball_around_point(b1, 1)
        ball2 = o.hilbert_ball_around_point(b2, 1)
        ballr = o.hilbert_ball_around_point(r, 1)
        for p in [b1, b2, r]:
            #spokes
            for v in o.vertices:
                plt.axline(p, v, color='gainsboro')

        plt.subplot(1, 2, 1)
        plt.title("standard")
        # Plot omega
        tools.plot_congruent(plt, o.vertices, color='gray')
        # Plot line ball
        tools.plot_congruent(plt, ball.vertices, color='lightblue')
        # Plot best dividing line
        plt.axline(*l.l, color='orange')
        # Plot hilbert balls
        tools.plot_congruent(plt, ball1.vertices, color='blue')
        tools.plot_congruent(plt, ball2.vertices, color='green')
        tools.plot_congruent(plt, ballr.vertices, color='red')
        # Plot points
        tools.scatter(plt, [b1, b2, r], colors=['blue', 'lime', 'red'])
        tools.annotate(plt, 'pqr', [b1, b2, r])

        # Dual2 PLOT
        plt.subplot(1, 2, 2)
        plt.title("Dual2 Dual")
        # Plot omega
        tools.plot_congruent(plt, Dual2.v2v(o.vertices), color='gray')
        # Plot line ball
        tools.plot_congruent(plt, Dual2.v2v(ball.vertices), color='lightblue', axline=False)

        # Plot hilbert balls
        tools.plot_congruent(plt, Dual2.v2v(ball1.vertices), color='blue', axline=True)
        tools.plot_congruent(plt, Dual2.v2v(ball2.vertices), color='green', axline=True)
        tools.plot_congruent(plt, Dual2.v2v(ballr.vertices), color='red', axline=True)
        # Plot points
        plt.axline(*Dual2.to_line(b1), color='blue')
        plt.axline(*Dual2.to_line(b2), color='lime')
        plt.axline(*Dual2.to_line(r), color='red')
        for spoke in l.get_ball_spokes():
            print(f"Spoke {spoke}")
            tools.scatter(plt, [Dual2.to_point(spoke[1], spoke[2])], colors=['gainsboro'])
        # Plot best dividing line
        plt.scatter(*Dual2.to_point(l.l), color='orange', zorder=1000)



        plt.tight_layout()
        plt.show()

    def test_three_balls_growing(self):
        o = omega.Omega(vertices=[(1.4, -0.1),
                                  (1.1, 0.8),
                                  (0.6, 1.1),
                                  (-0.2, 1.2),
                                  (-.5, 0.8),
                                  (-.8, 0.3),
                                  (-.7, -.2)])
        p, q = tools.rand_points(2)
        l = line.Line(p, q, o)
        ball = l.hilbert_ball_about_line(1)
        above = [x for x in ball.vertices if euclidean.point_below_line(x, l.l)]
        below = list(reversed([x for x in ball.vertices if not euclidean.point_below_line(x, l.l)]))
        b1, b2, r = None, None, None
        while b2 is None:
            try:
                b1, dist = euclidean.uniform_sample_from_line_segments(above)
                r, dist = euclidean.uniform_sample_from_line_segments(below, dist)
                b2, dist = euclidean.uniform_sample_from_line_segments(above, dist)
            except:
                pass
        for radm in range(1, 101):
            rad = radm / 100
            self.three_balls_growing(rad, b1, b2, r, o, l)



    def three_balls_growing(self, radius, b1, b2, r, o, l):

        ball1 = o.hilbert_ball_around_point(b1, radius)
        ball2 = o.hilbert_ball_around_point(b2, radius)
        ballr = o.hilbert_ball_around_point(r, radius)
        for p in [b1, b2, r]:
            #spokes
            for v in o.vertices:
                plt.axline(p, v, color='gainsboro')

        plt.subplot(1, 2, 1, ylim=(-1.5, 1.5), xlim=(-1.5, 1.5))
        plt.title("standard")
        # Plot omega
        tools.plot_congruent(plt, o.vertices, color='gray')
        # Plot best dividing line
        plt.axline(*l.l, color='orange')
        # Plot hilbert balls
        tools.plot_congruent(plt, ball1.vertices, color='blue')
        tools.plot_congruent(plt, ball2.vertices, color='lime')
        tools.plot_congruent(plt, ballr.vertices, color='red')
        # Plot points
        tools.scatter(plt, [b1, b2, r], colors=['blue', 'lime', 'red'])
        tools.annotate(plt, 'pqr', [b1, b2, r])

        # Dual2 PLOT
        orangedot = Dual2.to_point(l.l)
        plt.subplot(1, 2, 2, ylim=(orangedot[1]-2, orangedot[1]+2), xlim=(orangedot[0]-2, orangedot[0]+2))
        plt.title("Dual2 Dual")
        # Plot omega
        tools.plot_congruent(plt, Dual2.v2v(o.vertices), color='gray')

        # Plot hilbert balls
        tools.plot_congruent(plt, Dual2.v2v(ball1.vertices), color='blue', axline=False)
        tools.plot_congruent(plt, Dual2.v2v(ball2.vertices), color='lime', axline=False)
        tools.plot_congruent(plt, Dual2.v2v(ballr.vertices), color='red', axline=False)
        # Plot points
        plt.axline(*Dual2.to_line(b1), color='blue')
        plt.axline(*Dual2.to_line(b2), color='lime')
        plt.axline(*Dual2.to_line(r), color='red')
        # Plot best dividing line
        plt.scatter(*orangedot, color='orange', zorder=1000)



        plt.tight_layout()
        plt.savefig(f"./movieframes/{tools.namer(radius*100)}.png")
        plt.clf()









if __name__ == '__main__':
    unittest.main()
