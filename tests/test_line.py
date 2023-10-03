import unittest
from hilbert import omega, line, geometry
from misc import tools, euclidean
import numpy as np
from matplotlib import pyplot as plt
from pprint import pprint

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

    def test_sample_points_vertices(self):
        o = omega.Omega(vertices=[(1.4, 0.0),
                                  (1.1, 0.8),
                                  (0.6, 1.1),
                                  (0.0, 1.2),
                                  (-.5, 0.8),
                                  (-.8, 0.3),
                                  (-.7, -.2)])
        p1, p2 = np.array([0.5, 0.5]), np.array([0.6, 0.5])
        l = line.Line(p1, p2, o)
        ball = l.hilbert_ball_about_line(1)
        yooks, zooks = l.sample_points(0, ball)
        tools.plot_congruent(plt, o.vertices, color='gray')
        tools.plot_congruent(plt, ball.vertices, color='gray')
        tools.plot_congruent(plt, zooks, color='blue', connect_ends=True)
        tools.plot_congruent(plt, yooks, color='red', connect_ends=True)
        tools.annotate(plt, range(len(zooks)), zooks)
        tools.annotate(plt, range(len(yooks)), yooks)
        plt.show()

    def test_spokes(self):
        o = omega.Omega(vertices=[(1.4, 0.0),
                                  (1.1, 0.8),
                                  (0.6, 1.1),
                                  (0.0, 1.2),
                                  (-.5, 0.8),
                                  (-.8, 0.3),
                                  (-.7, -.2)])
        p1, p2 = np.array([0.5, 0.5]), np.array([0.6, 0.5])
        l = line.Line(p1, p2, o)
        b = l.hilbert_ball_about_line(1)
        pprint([x[1:] for x in l.get_ball_spokes()])
        tools.plot_congruent(plt, o.vertices, color='gray')
        tools.plot_congruent(plt, o.vertices, color='gray')
        tools.plot_congruent(plt, b.vertices, color='blue')
        for s in l.get_ball_spokes():
            tools.plot_line(plt, s[1], s[2], color='orange')
        tools.scatter(plt, b.vertices, colors=['gray'])
        plt.show()




    def test_sample_outside_ball(self):
        o = omega.Omega(vertices=[(1.4, 0.0),
                                  (1.1, 0.8),
                                  (0.6, 1.1),
                                  (0.0, 1.2),
                                  (-.5, 0.8),
                                  (-.8, 0.3),
                                  (-.7, -.2)])
        p1, p2 = np.array([0.5, 0.5]), np.array([0.6, 0.5])
        l = line.Line(p1, p2, o)
        b = l.hilbert_ball_about_line(1)
        pprint([x[1:] for x in l.get_ball_spokes()])
        tools.plot_congruent(plt, o.vertices, color='gray')
        tools.plot_congruent(plt, o.vertices, color='gray')
        tools.plot_congruent(plt, b.vertices, color='blue')
        for s in l.get_ball_spokes():
            tools.plot_line(plt, s[1], s[2], color='orange')
        tools.scatter(plt, b.vertices, colors=['gray'])
        blues, reds = l.sample_points_outside_ball(100, b)
        tools.scatter(plt, blues, colors=['blue'])
        tools.scatter(plt, reds, colors=['red'])

        plt.show()

    def test_spokes_of_support_vectors(self):
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

        svs = []
        for _ in range(100):
            svs.append(euclidean.uniform_sample_from_line_segments(above)[0])
            svs.append(euclidean.uniform_sample_from_line_segments(below)[0])

        from matplotlib.colors import TABLEAU_COLORS as COLORS
        COLORS = list(COLORS.values())
        print(COLORS)
        colors = [COLORS[l.nearest_point(x)[2]] for x in svs]
        tools.plot_congruent(plt, o.vertices, color='gray')
        plt.axline(*l.l, color='gray')
        tools.scatter(plt, svs, colors=colors)
        tools.scatter(plt, o.vertices, colors=COLORS)
        plt.show()



    def test_spokes_of_support_vectors_moving(self):
        o = omega.Omega(vertices=[(1.4, 0.0),
                                  (1.1, 0.8),
                                  (0.6, 1.1),
                                  (0.0, 1.2),
                                  (-.5, 0.8),
                                  (-.8, 0.3),
                                  (-.7, -.2)])
        jp, jq = np.array([0.0, 1.19]), np.array([-.49, 0.79])
        from matplotlib.colors import TABLEAU_COLORS as COLORS
        COLORS = list(COLORS.values()) + [None]
        for ji in range(210):
            try:
                ji = ji / 100
                p = jp + [ji, -ji]
                q = jq + [ji, -ji]
                l = line.Line(p, q, o)
                ball = l.hilbert_ball_about_line(1)
                above, below = ball.halves(l)
                #below = [x for x in ball.vertices if euclidean.point_below_line(x, l.l)]
                #above = list(reversed([x for x in ball.vertices if not euclidean.point_below_line(x, l.l)]))


                svs = []
                for _ in range(100):
                    svs.append(euclidean.uniform_sample_from_line_segments(above)[0])
                    svs.append(euclidean.uniform_sample_from_line_segments(below)[0])


                colors = [COLORS[l.nearest_point(x)[2]] for x in svs]
                tools.plot_congruent(plt, o.vertices, color='gray')
                plt.axline(*l.l, color='gray')
                tools.scatter(plt, svs, colors=colors)
                tools.scatter(plt, o.vertices, colors=COLORS)
                plt.savefig(f"./movieframes/{tools.namer(ji*100)}.png")
            except Exception as e:
                print(e, e.__doc__)
            plt.clf()
        cmd = "ffmpeg -framerate 20 -pattern_type glob -i \"./movieframes/*.png\"   -c:v libx264 -pix_fmt yuv420p ./movies/outD2.mp4"
        import os
        os.system(cmd)
        os.system("xdg-open ./movies/outD2.mp4")

    def test_spokes_of_support_vectors_spinning(self):

        o = omega.Omega(vertices=[(1.4, 0.0),
                                  (1.1, 0.8),
                                  (0.6, 1.1),
                                  (0.0, 1.2),
                                  (-.5, 0.8),
                                  (-.8, 0.3),
                                  (-.7, -.2)])
        center = o.polygon.center()
        from matplotlib.colors import TABLEAU_COLORS as COLORS
        COLORS = list(COLORS.values()) + [None]
        for ji in range(180):
            try:
                q = euclidean.circle_point(ji, 180)
                l = line.Line(center, q, o)
                ball = l.hilbert_ball_about_line(1)
                above, below = ball.halves(l)
                #below = [x for x in ball.vertices if euclidean.point_below_line(x, l.l)]
                #above = list(reversed([x for x in ball.vertices if not euclidean.point_below_line(x, l.l)]))


                svs = []
                for _ in range(100):
                    svs.append(euclidean.uniform_sample_from_line_segments(above)[0])
                    svs.append(euclidean.uniform_sample_from_line_segments(below)[0])


                colors = [COLORS[l.nearest_point(x)[2]] for x in svs]
                tools.plot_congruent(plt, o.vertices, color='gray')
                plt.axline(*l.l, color='gray')
                tools.scatter(plt, svs, colors=colors)
                tools.scatter(plt, o.vertices, colors=COLORS)
                for sp in l.get_ball_spokes(plt):
                    tools.plot_line(plt, sp[1], sp[2], color='gray')

                plt.xlim(-1, 1.5)
                plt.ylim(-.25, 1.25)
                plt.savefig(f"./movieframes/{tools.namer(ji)}.png")
            except Exception as e:
                print(e, e.__doc__)
            plt.clf()
        cmd = "ffmpeg -framerate 15 -pattern_type glob -i \"./movieframes/*.png\"   -c:v libx264 -pix_fmt yuv420p ./movies/outSpinSpokes.mp4"
        import os
        os.system(cmd)
        os.system("xdg-open ./movies/outSpinSpokes.mp4")


    def test_spokes_of_support_vectors_growing(self):
        o = omega.Omega(vertices=[(1.4, 0.0),
                                  (1.1, 0.8),
                                  (0.6, 1.1),
                                  (0.0, 1.2),
                                  (-.5, 0.8),
                                  (-.8, 0.3),
                                  (-.7, -.2)])
        jp, jq = tools.rand_points(2)
        l = line.Line(jp, jq, o)
        from matplotlib.colors import TABLEAU_COLORS as COLORS
        COLORS = list(COLORS.values()) + [None]
        for ji in range(200):
            try:
                radius = ji / 100


                ball = l.hilbert_ball_about_line(radius)
                above, below = ball.halves(l)
                #below = [x for x in ball.vertices if euclidean.point_below_line(x, l.l)]
                #above = list(reversed([x for x in ball.vertices if not euclidean.point_below_line(x, l.l)]))


                svs = []
                for _ in range(100):
                    svs.append(euclidean.uniform_sample_from_line_segments(above)[0])
                    svs.append(euclidean.uniform_sample_from_line_segments(below)[0])


                colors = [COLORS[l.nearest_point(x)[2]] for x in svs]
                tools.plot_congruent(plt, o.vertices, color='gray')
                plt.axline(*l.l, color='gray')

                for sp in l.get_ball_spokes(plt):
                    tools.plot_line(plt, sp[1], sp[2], color='gray')

                tools.scatter(plt, svs, colors=colors)
                tools.scatter(plt, o.vertices, colors=COLORS)
                plt.savefig(f"./movieframes/{tools.namer(ji)}.png")
            except Exception as e:
                print(e, e.__doc__)
            plt.clf()
        cmd = "ffmpeg -framerate 15 -pattern_type glob -i \"./movieframes/*.png\"   -c:v libx264 -pix_fmt yuv420p ./movies/outGrowSpokes.mp4"
        import os
        os.system(cmd)
        os.system("xdg-open ./movies/outGrowSpokes.mp4")


if __name__ == '__main__':
    unittest.main()
