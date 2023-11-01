import unittest
from hilbert import geometry, omega, polygon, line
from matplotlib import pyplot as plt
from misc import tools, euclidean
from misc.dual import Polar, Dual2
import numpy as np

class PolarTests(unittest.TestCase):

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
                                  (-.7, -.2)])

        plt.subplot(1, 2, 1)
        tools.plot_congruent(plt, o.vertices, color='red')
        plt.title("standard")
        plt.subplot(1, 2, 2)
        for i, v in enumerate(o.vertices):
            q = Polar.to_line(v)
            print(i,v, q)
            tools.annotate(plt,[i], [q[0]])
            tools.plot_line(plt, *q, color='blue')
            #plt.axline(*q, color="blue")

        scatpoints = Polar.v2v(o.vertices)
        print("Scat points: len", len(scatpoints), scatpoints)
        tools.scatter(plt, scatpoints)
        plt.title("polar dual")
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

        plt.subplot(1, 2, 2)
        plt.title("polar dual")
        """
        for i, v in enumerate(o.vertices):
            q = dual.Polar.to_line(v)
            print(i, v, q)
            tools.annotate(plt, [i], [q[0]])
            plt.axline(*q, color="blue")
        """
        plt.scatter(*Polar.to_point(l), color='gray')
        tools.plot_congruent(plt, Polar.v2v(o.vertices), color='red', axline=False)
        Polar.plot_congruent(plt, Polar.v2v(ball), color='orange')
        plt.axline(*Polar.to_line(p))
        plt.axline( *Polar.to_line(ints[0]), color='black')
        plt.axline( *Polar.to_line(ints[1]), color='green')
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
        #p = tools.rand_points(1)[0]
        p = np.array([-0.5, -0.1])
        tools.scatter(plt, [p], colors=['lime'])
        balls = [o.hilbert_ball_around_point(p, n).vertices for n in [0.5, 1, 2, 4]]
        colors = ['blue', 'lightblue', 'yellow', 'orange']
        for i, ball in enumerate(reversed(balls)):
            tools.plot_congruent(plt, ball, color=colors[i])

        plt.subplot(1, 2, 2)
        plt.title("polar dual")

        for i, ball in enumerate(reversed(balls)):
            Polar.plot_congruent(plt, Polar.v2v(ball), color=colors[i])
        plt.axline(*Polar.to_line(p), color='lime', zorder=100)
        tools.plot_congruent(plt, Polar.v2v(o.vertices), color='red')


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


        plt.subplot(1, 2, 1)
        plt.title("standard")
        # Plot omega
        tools.plot_congruent(plt, o.vertices, color='gray')
        for spoke in l.get_ball_spokes():
            print("Spoke???", spoke)
            tools.plot_line(plt, spoke[1], spoke[2], color='gainsboro')


        # Plot line ball
        tools.plot_congruent(plt, ball.vertices, color='lightblue')
        # Plot best dividing line
        plt.axline(*l.l, color='orange')
        # Plot hilbert balls
        tools.plot_congruent(plt, ball1.vertices, color='blue')
        tools.plot_congruent(plt, ball2.vertices, color='green')
        tools.plot_congruent(plt, ballr.vertices, color='red')
        # Plot points
        lbs = l.get_boundary_intersections()
        tools.scatter(plt, [b1, b2, r, lbs[0], lbs[1]], colors=['blue', 'lime', 'red', 'gold', 'yellow'])
        tools.annotate(plt, 'pqr', [b1, b2, r])



        # POLAR PLOT
        plt.subplot(1, 2, 2)
        plt.title("Polar Dual")
        # Plot omega
        transformed_o = Polar.v2v(o.vertices)
        Polar.plot_congruent(plt, polygon.Polygon(transformed_o), color='gray')

        # New, plot polar spokes connecting origin through vertices of new omega transform
        # Then send those spokes to points in standard space.





        # Plot line ball
        Polar.plot_congruent(plt, Polar.v2v(ball.vertices), color='lightblue')

        # Plot hilbert balls
        Polar.plot_congruent(plt, Polar.v2v(ball1.vertices), color='blue')
        Polar.plot_congruent(plt, Polar.v2v(ball2.vertices), color='green')
        Polar.plot_congruent(plt, Polar.v2v(ballr.vertices), color='red')
        # Plot points
        plt.axline(*Polar.to_line(b1), color='blue')
        plt.axline(*Polar.to_line(b2), color='lime')
        plt.axline(*Polar.to_line(r), color='red')
        plt.axline(*Polar.to_line(l.get_boundary_intersections()[0]), color='gold')
        plt.axline(*Polar.to_line(l.get_boundary_intersections()[1]), color='yellow')
        for spoke in l.get_ball_spokes():
            tools.scatter(plt, [Polar.to_point(spoke[1], spoke[2])], colors=['gainsboro'])


        # Plot best dividing line
        plt.scatter(*Polar.to_point(l.l), color='orange', zorder=1000)

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

        # POLAR PLOT
        orangedot = Polar.to_point(l.l)
        plt.subplot(1, 2, 2, ylim=(orangedot[1]-2, orangedot[1]+2), xlim=(orangedot[0]-2, orangedot[0]+2))
        plt.title("Polar Dual")
        # Plot omega
        tools.plot_congruent(plt, Polar.v2v(o.vertices), color='gray')

        # Plot hilbert balls
        tools.plot_congruent(plt, Polar.v2v(ball1.vertices), color='blue', axline=False)
        tools.plot_congruent(plt, Polar.v2v(ball2.vertices), color='lime', axline=False)
        tools.plot_congruent(plt, Polar.v2v(ballr.vertices), color='red', axline=False)
        # Plot points
        plt.axline(*Polar.to_line(b1), color='blue')
        plt.axline(*Polar.to_line(b2), color='lime')
        plt.axline(*Polar.to_line(r), color='red')
        # Plot best dividing line
        plt.scatter(*orangedot, color='orange', zorder=1000)



        plt.tight_layout()
        plt.savefig(f"./movieframes/{tools.namer(radius*100)}.png")
        plt.clf()


    def test_many_sampled_points(self):
        o = omega.Omega(offset=[0.1,0.1])
        """vertices=[(1.4, -0.1),
                                  (1.1, 0.8),
                                  (0.6, 1.1),
                                  (-0.2, 1.2),
                                  (-.5, 0.8),
                                  (-.8, 0.3),
                                  (-.7, -.2)]) #, offset=[1,.3])"""
        p, q = tools.rand_points(2) # , offset=[1,.3])
        l = line.Line(p, q, o)
        ball = l.hilbert_ball_about_line(1)
        above = [x for x in ball.vertices if euclidean.point_below_line(x, l.l)]
        below = list(reversed([x for x in ball.vertices if not euclidean.point_below_line(x, l.l)]))
        b1, dist = euclidean.uniform_sample_from_line_segments(above)
        r, dist = euclidean.uniform_sample_from_line_segments(below, dist)
        b2, _ = euclidean.uniform_sample_from_line_segments(above, dist)
        reds, blues = l.sample_points_outside_ball(40, ball)
        boundedges = l.get_boundaries()
        orgx, orgy = Polar.to_point(l.l)

        #
        #
        #
        # STANDARD PLOT 1
        #
        #
        #
        plt.subplot(1, 3, 1)
        plt.title("standard")
        # Plot omega
        tools.plot_congruent(plt, o.vertices, color='gray')
        # Plot best dividing line
        plt.axline(*l.l, color='orange')

        # Plot points
        lbs = l.get_boundary_intersections()
        tools.scatter(plt, [b1, b2, r, lbs[0], lbs[1]], colors=['blue', 'lime', 'red', 'gold', 'yellow'])
        tools.annotate(plt, 'pqr', [b1, b2, r])

        tools.scatter(plt, [boundedges[0][0], boundedges[0][1], boundedges[1][0], boundedges[1][1]], colors=['black'])
        tools.scatter(plt, blues, colors=['lightblue'])
        tools.scatter(plt, reds, colors=['salmon'])
        #fixedpoint = euclidean.get_point_on_line(Polar.origin, Polar.to_point(l.l), 0.5)
        #print('heyo', Polar.to_point(fixedpoint, Polar.to_point(l.l)))
        #tools.scatter(plt, [Polar.to_point(fixedpoint, Polar.to_point(l.l))], colors=['purple'], zorder=1000)

        #
        #
        # POLAR PLOT 1
        #
        #
        #


        plt.subplot(1, 3, 2, xlim=(orgx-1.2,orgx+1.2), ylim=(orgy-1.2, orgy+1.2))
        plt.title("Polar Dual")
        for point in blues:
            plt.axline( *Polar.to_line(point), color='lightblue', zorder=-1000)
        for point in reds:
            plt.axline( *Polar.to_line(point), color='salmon', zorder=-1000)
        # Plot omega
        transformed_o = Polar.v2v(o.vertices)
        Polar.plot_congruent(plt, polygon.Polygon(transformed_o).vertices, color='gray', zorder=1000)
        Polar.plot_congruent(plt, polygon.Polygon(transformed_o).vertices, color='gray', axline=True, zorder=1000)
        # Plot points
        plt.axline(*Polar.to_line(b1), color='blue')
        plt.axline(*Polar.to_line(b2), color='lime')
        plt.axline(*Polar.to_line(r), color='red')
        plt.axline(*Polar.to_line(l.get_boundary_intersections()[0]), color='gold')
        plt.axline(*Polar.to_line(l.get_boundary_intersections()[1]), color='yellow')
        for bep in [boundedges[0][0], boundedges[0][1], boundedges[1][0], boundedges[1][1]]:
            plt.axline(*Polar.to_line(bep), color='black', zorder=1010)

        # Plot best dividing line
        plt.scatter(*Polar.to_point(l.l), color='orange', zorder=1000)


        #
        #
        # POLAR PLOT 2
        #
        #
        #

        plt.subplot(1, 3, 3)#, xlim=(min(orgx, o..offset[0]-1) - 1, max(orgx, o..offset[0] + 1) + 1),
                   # ylim=(min(orgy, o..offset[1] - 6) - 0.5, max(orgy, o..offset[1] + 1) + 0.5))
        plt.title("Polar Dual")
        for point in blues:
            plt.axline(*Polar.to_line(point), color='lightblue', zorder=-1000)
        for point in reds:
            plt.axline(*Polar.to_line(point), color='salmon', zorder=-1000)

        # Plot omega
        transformed_o = Polar.v2v(o.vertices)
        Polar.plot_congruent(plt, polygon.Polygon(transformed_o).vertices, color='gray', zorder=1000)
        Polar.plot_congruent(plt, polygon.Polygon(transformed_o).vertices, color='gray', axline=True)

        # New, plot polar spokes connecting origin through vertices of new omega transform
        # Then send those spokes to points in standard space.

        # Plot points
        plt.axline(*Polar.to_line(b1), color='blue', zorder=1000)
        plt.axline(*Polar.to_line(b2), color='lime', zorder=1000)
        plt.axline(*Polar.to_line(r), color='red', zorder=1000)
        plt.axline(*Polar.to_line(l.get_boundary_intersections()[0]), color='gold')
        plt.axline(*Polar.to_line(l.get_boundary_intersections()[1]), color='yellow')
        for bep in [boundedges[0][0], boundedges[0][1], boundedges[1][0], boundedges[1][1]]:
            plt.axline(*Polar.to_line(bep), color='black', zorder=1000)
        # Plot best dividing line
        plt.scatter(*Polar.to_point(l.l), color='orange', zorder=1000)


        plt.tight_layout()
        plt.show()

    def test_visualization_on_many_support_vectors(self):
        # o = omega.Omega(offset=[0.1,0.1])
        o = omega.Omega(vertices=[(1.4, -0.1),
                                  (1.1, 0.8),
                                  (0.6, 1.1),
                                  (-0.2, 1.2),
                                  (-.5, 0.8),
                                  (-.8, 0.3),
                                  (-.7, -.2)], offset=[0, .2])
        p, q = tools.rand_points(2, offset=[0, .2])
        l = line.Line(p, q, o)
        ball = l.hilbert_ball_about_line(1)
        above = [x for x in ball.vertices if euclidean.point_below_line(x, l.l)]
        below = list(reversed([x for x in ball.vertices if not euclidean.point_below_line(x, l.l)]))

        reds, blues = [], []
        for _ in range(100):
            blues.append(euclidean.uniform_sample_from_line_segments(above)[0])
            reds.append(euclidean.uniform_sample_from_line_segments(below)[0])

        boundedges = l.get_boundaries()
        orgx, orgy = Polar.to_point(l.l)

        #
        #
        #
        # STANDARD PLOT 1
        #
        #
        #
        plt.subplot(1, 3, 1)
        plt.title("standard")
        # Plot omega
        tools.plot_congruent(plt, o.vertices, color='gray')
        # Plot best dividing line
        plt.axline(*l.l, color='orange')

        # Plot points
        lbs = l.get_boundary_intersections()
        tools.scatter(plt, [lbs[0], lbs[1]], colors=['gold', 'yellow'])

        tools.scatter(plt, [boundedges[0][0], boundedges[0][1], boundedges[1][0], boundedges[1][1]], colors=['black'])
        tools.scatter(plt, blues, colors=['lightblue'])
        tools.scatter(plt, reds, colors=['salmon'])
        # fixedpoint = euclidean.get_point_on_line(Polar.origin, Polar.to_point(l.l), 0.5)
        # print('heyo', Polar.to_point(fixedpoint, Polar.to_point(l.l)))
        # tools.scatter(plt, [Polar.to_point(fixedpoint, Polar.to_point(l.l))], colors=['purple'], zorder=1000)

        #
        #
        # POLAR PLOT 1
        #
        #
        #

        plt.subplot(1, 3, 2, xlim=(orgx - 1.2, orgx + 1.2), ylim=(orgy - 1.2, orgy + 1.2))
        plt.title("Polar Dual")
        for point in blues:
            plt.axline(*Polar.to_line(point), color='lightblue', zorder=-1000)
        for point in reds:
            plt.axline(*Polar.to_line(point), color='salmon', zorder=-1000)
        # Plot omega
        transformed_o = Polar.v2v(o.vertices)
        Polar.plot_congruent(plt, polygon.Polygon(transformed_o), color='gray', zorder=1000)
        Polar.plot_congruent(plt, polygon.Polygon(transformed_o), color='gray', axline=True, zorder=1000)
        # Plot points

        plt.axline(*Polar.to_line(l.get_boundary_intersections()[0]), color='gold')
        plt.axline(*Polar.to_line(l.get_boundary_intersections()[1]), color='yellow')
        for bep in [boundedges[0][0], boundedges[0][1], boundedges[1][0], boundedges[1][1]]:
            plt.axline(*Polar.to_line(bep), color='black', zorder=1010)

        # Plot best dividing line
        plt.scatter(*Polar.to_point(l.l), color='orange', zorder=1000)

        #
        #
        # POLAR PLOT 2
        #
        #
        #

        plt.subplot(1, 3, 3, xlim=(min(orgx, o..offset[0] - 1) - 1, max(orgx, o..offset[0] + 1) + 3),
                    ylim=(min(orgy, o..offset[1] - 6) - 0.5, max(orgy, o..offset[1] + 1) + 5.5))
        plt.title("Polar Dual")
        for point in blues:
            plt.axline(*Polar.to_line(point), color='lightblue', zorder=-1000)
        for point in reds:
            plt.axline(*Polar.to_line(point), color='salmon', zorder=-1000)

        # Plot omega
        transformed_o = Polar.v2v(o.vertices)
        Polar.plot_congruent(plt, polygon.Polygon(transformed_o), color='gray', zorder=1000)
        Polar.plot_congruent(plt, polygon.Polygon(transformed_o), color='gray', axline=True)

        # New, plot polar spokes connecting origin through vertices of new omega transform
        # Then send those spokes to points in standard space.

        # Plot points

        plt.axline(*Polar.to_line(l.get_boundary_intersections()[0]), color='gold')
        plt.axline(*Polar.to_line(l.get_boundary_intersections()[1]), color='yellow')
        for bep in [boundedges[0][0], boundedges[0][1], boundedges[1][0], boundedges[1][1]]:
            plt.axline(*Polar.to_line(bep), color='black', zorder=1000)
        # Plot best dividing line
        plt.scatter(*Polar.to_point(l.l), color='orange', zorder=1000)

        plt.tight_layout()
        plt.show()



    def test_visualization_on_many_support_vectors_with_boundary_colors(self):
        #o = omega.Omega(offset=[0.1,0.1])
        o = omega.Omega(vertices=[(1.4, -0.1),
                                  (1.1, 0.8),
                                  (0.6, 1.1),
                                  (-0.2, 1.2),
                                  (-.5, 0.8),
                                  (-.8, 0.3),
                                  (-.7, -.2)], offset=[0,.21])
        p, q = tools.rand_points(2 , offset=[0,.2])
        l = line.Line(p, q, o)
        ball = l.hilbert_ball_about_line(1)
        above, below = ball.halves(l)

        svs = []
        for _ in range(2):
            svs.append(euclidean.uniform_sample_from_line_segments(above)[0])
            svs.append(euclidean.uniform_sample_from_line_segments(below)[0])
        from matplotlib.colors import TABLEAU_COLORS as COLORS
        COLORS = list(COLORS.values())
        colors = [COLORS[l.nearest_point(x)[2]] for x in svs]
        boundedges = l.get_boundaries()
        orgx, orgy = Polar.to_point(l.l)

        #
        #
        #
        # STANDARD PLOT 1
        #
        #
        #
        plt.subplot(1, 3, 1)
        plt.title("standard")
        # Plot omega
        tools.plot_congruent(plt, o.vertices, color='gray')
        # Plot best dividing line
        plt.axline(*l.l, color='orange')

        # Plot points
        lbs = l.get_boundary_intersections()
        tools.scatter(plt, [lbs[0], lbs[1]], colors=['gold', 'yellow'])

        tools.scatter(plt, o.vertices, colors=COLORS)
        tools.scatter(plt, svs, colors=colors)
        #fixedpoint = euclidean.get_point_on_line(Polar.origin, Polar.to_point(l.l), 0.5)
        #print('heyo', Polar.to_point(fixedpoint, Polar.to_point(l.l)))
        #tools.scatter(plt, [Polar.to_point(fixedpoint, Polar.to_point(l.l))], colors=['purple'], zorder=1000)

        #
        #
        # POLAR PLOT 1
        #
        #
        #


        plt.subplot(1, 3, 2, xlim=(orgx-1.2,orgx+1.2), ylim=(orgy-1.2, orgy+1.2))
        plt.title("Polar Dual")
        for i, point in enumerate(svs):
            linestyle = 'dotted' if i%2==0 else "dashed"
            plt.axline( *Polar.to_line(point), color=colors[i], zorder=-1000, linestyle=linestyle)

        # Plot omega
        #transformed_o = Polar.v2v(o.vertices)
        #Polar.plot_congruent(plt, transformed_o, color='gray', zorder=1000)

        for i, v in enumerate(o.vertices):
            plt.axline(*Polar.to_line(v), color=COLORS[i], linewidth=3, zorder=1000)
            plt.axline(*Polar.to_line(v), color='black', linewidth=0.5, zorder=1010)

        # Plot points

        plt.axline(*Polar.to_line(l.get_boundary_intersections()[0]), color='gold')
        plt.axline(*Polar.to_line(l.get_boundary_intersections()[1]), color='yellow')
        #for bep in [boundedges[0][0], boundedges[0][1], boundedges[1][0], boundedges[1][1]]:
        #    plt.axline(*Polar.to_line(bep), color='black', zorder=1010)

        # Plot best dividing line
        plt.scatter(*Polar.to_point(l.l), color='orange', zorder=1000)


        #
        #
        # POLAR PLOT 2
        #
        #
        #

        plt.subplot(1, 3, 3, xlim=(min(orgx, o..offset[0]-1) - 1, max(orgx, o..offset[0] + 1) + 3),
                    ylim=(min(orgy, o..offset[1] - 6) - 0.5, max(orgy, o..offset[1] + 1) + 5.5))
        plt.title("Polar Dual")
        for i, point in enumerate(svs):
            linestyle = 'dotted' if i%2==0 else "dashed"
            plt.axline(*Polar.to_line(point), color=colors[i], zorder=-1000, linestyle=linestyle)

        # Plot omega
        for i, v in enumerate(o.vertices):
            print(i, COLORS[i], v)
            plt.axline(*Polar.to_line(v), color=COLORS[i], linewidth=3, zorder=1000)
            plt.axline(*Polar.to_line(v), color='black', linewidth=0.5, zorder=1010)

        #transformed_o = Polar.v2v(o.vertices)
        #print("Transformed", transformed_o)
        #Polar.plot_congruent(plt, transformed_o, color='gray', zorder=1000)
        #print("Polar")
        #Polar.plot_congruent(plt, polygon.Polygon(transformed_o), color='black', axline=True, zorder=1000, linewidth=3)
        #Polar.plot_congruent(plt, transformed_o, color='black', axline=True, zorder=1010, linewidth=1)

        # New, plot polar spokes connecting origin through vertices of new omega transform
        # Then send those spokes to points in standard space.

        # Plot points

        plt.axline(*Polar.to_line(l.get_boundary_intersections()[0]), color='gold')
        plt.axline(*Polar.to_line(l.get_boundary_intersections()[1]), color='yellow')
        #for bep in [boundedges[0][0], boundedges[0][1], boundedges[1][0], boundedges[1][1]]:
        #    plt.axline(*Polar.to_line(bep), color='black', zorder=1000)
        # Plot best dividing line
        plt.scatter(*Polar.to_point(l.l), color='orange', zorder=1000)


        plt.tight_layout()
        plt.show()



class PolarWedgeTests(unittest.TestCase):
    def test_all_wedges_of_3_sv_points(self):
        o = omega.Omega(vertices=[(1.4, -0.1),
                                  (1.1, 0.8),
                                  (0.6, 1.1),
                                  (-0.2, 1.2),
                                  (-.5, 0.8),
                                  (-.8, 0.3),
                                  (-.7, -.2)], offset=[.25, 0.05])
        aa, bb = tools.rand_points(2, offset=[.25, 0.05])
        l = line.Line(aa, bb, o)
        ball = l.hilbert_ball_about_line(1)
        above, below = ball.halves(l)

        plt.subplot(1, 2, 1)
        p = euclidean.uniform_sample_from_line_segments(above)[0]
        q = euclidean.uniform_sample_from_line_segments(above)[0]
        r = euclidean.uniform_sample_from_line_segments(below)[0]
        print(f"p {p}, q {q}, r {r}")
        pp, pq, pr = Polar.to_line(p), Polar.to_line(q), Polar.to_line(r)

        po = Polar.polygon(o.vertices_coords())
        ppintersections = [euclidean.intersect(pp, [po.v(i).v, po.v(ii).v]) for i, ii in tools.i_ii(len(po))]
        pqintersections = [euclidean.intersect(pq, [po.v(i).v, po.v(ii).v]) for i, ii in tools.i_ii(len(po))]
        printersections = [euclidean.intersect(pr, [po.v(i).v, po.v(ii).v]) for i, ii in tools.i_ii(len(po))]
        print(f"ppintersections {len(ppintersections)}")
        print(f"pqintersections {len(pqintersections)}")
        print(f"printersections {len(printersections)}")

        def other_tangent(point, i):
            t = None
            try:
                t = po.point_tangent_linear(point, True).v
            except:
                pass
            try:
                # on the same line
                if t is None or euclidean.eq(t, po.v(i).v) or euclidean.eq(t, po.v(i + 1).v):
                    t = po.point_tangent_linear(point, False).v
            except:
                pass
            if t is None:
                print(f"No tangent found for {point}, {i}")
            return t
        pptangents = [other_tangent(point, i) for i, point in enumerate(ppintersections)]
        pqtangents = [other_tangent(point, i) for i, point in enumerate(pqintersections)]
        prtangents = [other_tangent(point, i) for i, point in enumerate(printersections)]
        print(f"pptangents {len(pptangents)}")
        print(f"pqtangents {len(pqtangents)}")
        print(f"prtangents {len(prtangents)}")

        ppwedge = [(np.array([pptangents[i], po.v(i)]), np.array([pptangents[i], po.v(ii)])) for i, ii in tools.i_ii(len(po.vertices))]
        pqwedge = [(np.array([pqtangents[i], po.v(i)]), np.array([pqtangents[i], po.v(ii)])) for i, ii in tools.i_ii(len(po.vertices))]
        prwedge = [(np.array([prtangents[i], po.v(i).v]), np.array([prtangents[i], po.v(ii)])) for i, ii in tools.i_ii(len(po.vertices))]
        print(f"ppwedge {len(ppwedge)}, {ppwedge}")
        print(f"pqwedge {len(pqwedge)}, {pqwedge}")
        print(f"prwedge {len(prwedge)}, {prwedge}")

        orgx, orgy = Polar.to_point(l.l)
        plt.axline(*pp, color='blue')
        plt.axline(*pq, color='green')
        plt.axline(*pr, color='red')
        for i, ii in tools.i_ii(len(po.vertices)):
            #plt.axline(po.v(i), po.v(ii), color='gray', linewidth=0.4)
            tools.plot_line(plt, po.v(i).v, po.v(ii).v, color='gray', linewidth=2)
            plt.scatter(*ppintersections[i], color='blue')
            plt.scatter(*pqintersections[i], marker='^', color='green')
            plt.scatter(*printersections[i], marker='*', color='red')
            #tools.plot_line(plt, ppintersections[i], pptangents[i], color='blue', linewidth=0.5)
            #tools.plot_line(plt, pqintersections[i], pqtangents[i], color='green', linewidth=0.5)
            #tools.plot_line(plt, printersections[i], prtangents[i], color='red', linewidth=0.5)
            tools.plot_line(plt, *ppwedge[i][0], color='blue', linewidth=0.7, axline=True)
            tools.plot_line(plt, *ppwedge[i][1], color='blue', linewidth=0.7, axline=True)
            tools.plot_line(plt, *pqwedge[i][0], color='green', linewidth=0.6, axline=True, linestyle="dashed")
            tools.plot_line(plt, *pqwedge[i][1], color='green', linewidth=0.6, axline=True, linestyle="dashed")
            tools.plot_line(plt, *prwedge[i][0], color='red', linewidth=0.5, axline=True, linestyle="dotted")
            tools.plot_line(plt, *prwedge[i][1], color='red', linewidth=0.5, axline=True, linestyle="dotted")

        plt.scatter(orgx, orgy, color='orange')

        plt.xlim(-2, 2)
        plt.ylim(-2, 3)


        """
        STANDARD PLOT
        """
        plt.subplot(1, 2, 2)
        colors = ['blue', 'green', 'red']
        tools.plot_congruent(plt, o.vertices, color='gray')
        tools.scatter(plt, [p, q, r], colors=colors)
        for i, ii in tools.i_ii(len(po.vertices)):
            ## plt.axline(po.v(i), po.v(ii), color='gray', linewidth=0.4)
            #tools.plot_line(plt, po.v(i), po.v(ii), color='gray', linewidth=2)
            #plt.scatter(*ppintersections[i], color='blue')
            #plt.scatter(*pqintersections[i], marker='^', color='green')
            #plt.scatter(*printersections[i], marker='*', color='red')
            ## tools.plot_line(plt, ppintersections[i], pptangents[i], color='blue', linewidth=0.5)
            ## tools.plot_line(plt, pqintersections[i], pqtangents[i], color='green', linewidth=0.5)
            ## tools.plot_line(plt, printersections[i], prtangents[i], color='red', linewidth=0.5)
            tools.scatter(plt, [Polar.to_point(wedge[i][0]) for wedge in [ppwedge, pqwedge, prwedge]], colors=colors)
            tools.scatter(plt, [Polar.to_point(wedge[i][1]) for wedge in [ppwedge, pqwedge, prwedge]], colors=colors)
            """
            tools.plot_line(plt, *ppwedge[i][0], color='blue', linewidth=0.7, axline=True)
            tools.plot_line(plt, *ppwedge[i][1], color='blue', linewidth=0.7, axline=True)
            tools.plot_line(plt, *pqwedge[i][0], color='green', linewidth=0.6, axline=True, linestyle="dashed")
            tools.plot_line(plt, *pqwedge[i][1], color='green', linewidth=0.6, axline=True, linestyle="dashed")
            tools.plot_line(plt, *prwedge[i][0], color='red', linewidth=0.5, axline=True, linestyle="dotted")
            tools.plot_line(plt, *prwedge[i][1], color='red', linewidth=0.5, axline=True, linestyle="dotted")
            """
        plt.show()

        import sys
        sys.exit()







        from matplotlib.colors import TABLEAU_COLORS as COLORS
        COLORS = list(COLORS.values())
        #colors = [COLORS[l.nearest_point(x)[2]] for x in svs]


        #
        #
        #
        # STANDARD PLOT 1
        #
        #
        #
        plt.subplot(1, 3, 1)
        plt.title("standard")
        # Plot omega
        tools.plot_congruent(plt, o.vertices, color=COLORS)
        # Plot best dividing line
        plt.axline(*l.l, color='orange')
        colors = ['blue', 'green', 'red']
        tools.scatter(plt, [p,q,r], colors=colors)

        # Plot points
        #lbs = l.get_boundary_intersections()
        #tools.scatter(plt, [lbs[0], lbs[1]], colors=['gold', 'yellow'])


        #
        #
        # POLAR PLOT 1
        #
        #
        #

        plt.subplot(1, 3, 2, xlim=(orgx - 1.2, orgx + 1.2), ylim=(orgy - 1.2, orgy + 1.2))
        plt.title("Polar Dual")

        for i, v in enumerate(o.vertices):
            plt.axline(*Polar.to_line(v), color=COLORS[i], zorder=1000, linestyle="dotted")
            plt.axline(*Polar.to_line(v), color='black', linewidth=0.5, zorder=1010, linestyle="dotted")

        for i, (point, color) in enumerate(zip([p, q, r], colors)):
            plt.axline(*Polar.to_line(point), color=color)

        # Plot points

        #plt.axline(*Polar.to_line(l.get_boundary_intersections()[0]), color='gold')
        #plt.axline(*Polar.to_line(l.get_boundary_intersections()[1]), color='yellow')

        # Plot best dividing line
        plt.scatter(*Polar.to_point(l.l), color='orange', zorder=1000)

        #
        #
        # POLAR PLOT 2
        #
        #
        #

        plt.subplot(1, 3, 3, xlim=(min(orgx, o..offset[0] - 1) - 1, max(orgx, o..offset[0] + 1) + 3),
                    ylim=(min(orgy, o..offset[1] - 6) - 0.5, max(orgy, o..offset[1] + 1) + 5.5))
        plt.title("Polar Dual")

        for i, v in enumerate(o.vertices):
            plt.axline(*Polar.to_line(v), color=COLORS[i], zorder=1000, linestyle="dotted")
            plt.axline(*Polar.to_line(v), color='black', linewidth=0.5, zorder=1010, linestyle="dotted")

        for i, (point, color) in enumerate(zip([p, q, r], colors)):
            plt.axline(*Polar.to_line(point), color=color)

        plt.scatter(*Polar.to_point(l.l), color='orange', zorder=1000)

        plt.tight_layout()
        plt.show()






class Test_Dist(unittest.TestCase):

    def test_dist_is_polar_dist(self):
        o = omega.Omega(vertices=[(1.4, -0.1),
                                  (1.1, 0.8),
                                  (0.6, 1.1),
                                  (-0.2, 1.2),
                                  (-.5, 0.8),
                                  (-.8, 0.3),
                                  (-.7, -.2)])

        for _ in range(1000):
            p, q = tools.rand_points(2)
            l = line.Line(p, q, o)
            A, B = l.get_boundary_intersections()
            polar_dist = Polar.dist(Polar.to_line(p), Polar.to_line(q), Polar.to_line(A), Polar.to_line(B))
            print(f"Standard dist: {l.get_hdist()}, Polar dist: {polar_dist}. Difference: {l.get_hdist() - polar_dist}")
            self.assertAlmostEqual(l.get_hdist(), polar_dist)







if __name__ == '__main__':
    unittest.main()
