from unittest import TestCase

import hilbert.omega
from src.hilbert import omega, line, geometry
from src.misc import tools
import numpy as np
from matplotlib import pyplot as plt

class Test_Simplex_Ball(TestCase):
    def test_horizontal(self):
        o = omega.Omega()
        l = line.Line(np.array([0.1, 0.4]), np.array([0.3, 0.4]), o)
        plt.axline(l.p, l.q)
        tools.plot_congruent(plt, o.vertices)
        hilbert_ball = l.hilbert_ball_about_line(1)
        tools.plot_congruent(plt, hilbert_ball.vertices)
        plt.show()

    def test_vertical(self):
        o = omega.Omega()
        l = line.Line(np.array([0.4, 0.1]), np.array([0.4, 0.3]), o)
        plt.axline(l.p, l.q)
        tools.plot_congruent(plt, o.vertices)
        hilbert_ball = l.hilbert_ball_about_line(1)
        tools.plot_congruent(plt, hilbert_ball.vertices)
        plt.show()

    def test_diag(self):
        o = omega.Omega()
        l = line.Line(np.array([0.4, 0.1]), np.array([0.5, 0.2]), o)
        plt.axline(l.p, l.q)
        tools.plot_congruent(plt, o.vertices)
        hilbert_ball = l.hilbert_ball_about_line(1)
        tools.plot_congruent(plt, hilbert_ball.vertices)
        plt.show()



class Test_Weird_Ball(TestCase):
    o = omega.Omega(vertices=[(1.4, 0.0),
                                    (1.1, 0.8),
                                    (0.6, 1.1),
                                    (0.0, 1.2),
                                    (-.5, 0.8),
                                    (-.8, 0.3),
                                    (-.7, -.2)])

    def test_horizontal(self):
        l = line.Line(np.array([0.1, 0.4]), np.array([0.3, 0.4]), self.o)

        hilbert_ball = l.hilbert_ball_about_line(1, plt)
        tools.plot_congruent(plt, hilbert_ball.vertices)
        plt.scatter(*l.p)
        plt.scatter(*l.q)
        # print("Len of l. ball spokes, in test,", len(l.get_ball_spokes()))
        for bs in l.get_ball_spokes(plt):
            a, b = bs[1], bs[2]
            tools.plot_line(plt, a, b, color='gray')
            tools.plot_line(plt, bs[0], a, color='purple')
            tools.plot_line(plt, bs[0], b, color='red')
        plt.axline(l.p, l.q, color='orange')
        tools.plot_congruent(plt, self.o.vertices, color='blue')
        plt.show()

    def test_vertical(self):
        l = line.Line(np.array([0.4, 0.1]), np.array([0.4, 0.3]), self.o)
        plt.axline(l.p, l.q)
        tools.plot_congruent(plt, self.o.vertices)
        hilbert_ball = l.hilbert_ball_about_line(1)
        tools.plot_congruent(plt, hilbert_ball.vertices)
        plt.show()

    def test_diag(self):
        l = line.Line(np.array([0.4, 0.1]), np.array([0.5, 0.2]), self.o)
        plt.axline(l.p, l.q)
        tools.plot_congruent(plt, self.o.vertices)
        hilbert_ball = l.hilbert_ball_about_line(1)
        tools.plot_congruent(plt, hilbert_ball.vertices)
        plt.show()


class Test_Weird_Spokes(TestCase):

    o = omega.Omega(vertices=[(1.4, 0.0),
                              (1.1, 0.8),
                              (0.6, 1.1),
                              (0.0, 1.2),
                              (-.5, 0.8),
                              (-.8, 0.3),
                              (-.7, -.2)])

    def spokes(self, p, q):
        l = line.Line(p, q, self.o)
        plt.scatter(*l.p)
        plt.scatter(*l.q)
        tools.plot_congruent(plt, self.o.vertices, color='blue')
        for bs in l.get_ball_spokes(plt):
            a, b = bs[1], bs[2]
            tools.plot_line(plt, a, b, color='gray')
        plt.axline(l.p, l.q, color='orange')
        plt.axline(l.p, l.q)
        plt.show()



    def test_horizontal(self):
        p, q = np.array([0.1, 0.4]), np.array([0.3, 0.4])
        self.spokes(p, q)

    def test_vertical(self):
        p, q = np.array([0.4, 0.1]), np.array([0.4, 0.3])
        self.spokes(p, q)

    def test_diag(self):
        p, q =  np.array([0.4, 0.1]), np.array([0.5, 0.2])
        self.spokes(p, q)

    def test_rand(self):
        p, q = [x[:2] for x in np.random.dirichlet((1, 1, 1), 2)]
        self.spokes(p, q)

    def test_broken_weird_spokes(self):
        p, q = [0.13927992, 0.61490626], [0.06157955,  0.85032511]
        self.spokes(p, q)


class Hilbert_Ball_2pt(TestCase):
    def test_random_simplex(self):
        p, q = [x[:2] for x in np.random.dirichlet((1, 1, 1), 2)]
        o = omega.Omega()
        l = line.Line(p, q, o)
        fit = line.Line(*l.get_best_dividing_line(), o)
        ball = fit.hilbert_ball_about_line(radius=l.hdist/2)
        print("b:", ball)
        tools.plot_congruent(plt, o.vertices)
        tools.plot_congruent(plt, ball.vertices)
        plt.scatter(*p)
        plt.scatter(*q)
        plt.axline(fit.p, fit.q)
        plt.show()

    def weird(self, p, q):
        o = omega.Omega(vertices=[(1.4, 0.0),
                                  (1.1, 0.8),
                                  (0.6, 1.1),
                                  (0.0, 1.2),
                                  (-.5, 0.8),
                                  (-.8, 0.3),
                                  (-.7, -.2)])
        """
        l connects p and q as {A--p--q--B}
        fit is the "perpendicular" line, dividing p-q through the vanishing point.
        fit takes in a radius, which is going to be half the distance from p to q. 
        """
        print(p, q)
        l = line.Line(p, q, o)
        tools.plot_congruent(plt, o.vertices, color='gray')
        tools.scatter(plt, [p, q, l.get_midpoint()])
        tools.plot_line(plt, *l.get_boundaries()[0], color='green')
        tools.plot_line(plt, *l.get_boundaries()[1], color='purple')
        tools.annotate(plt, 'pqab', [p, q, l.get_boundaries()[0][0], l.get_boundaries()[1][0]])
        div = l.get_best_dividing_line()
        fit = line.Line(div[0], div[1], o)
        radius = l.get_hdist() / 2
        ball = fit.hilbert_ball_about_line(radius)
        tools.plot_congruent(plt, ball.vertices, color="orange")
        plt.axline(fit.p, fit.q)
        plt.show()

    def test_random_weird(self):
        p, q = [x[:2] for x in np.random.dirichlet((1, 1, 1), 2)]
        self.weird(p, q)

    def test_broken_weird_b(self):
        p, q = [0.03992638, 0.1748575 ], [0.28178188, 0.27469408]
        self.weird(p, q)

    def test_broken_weird_a(self):
        p, q = [0.13927992, 0.61490626], [0.06157955,  0.85032511]
        self.weird(p, q)

    def test_broken_weird_c(self):
        p, q = [0.42843028, 0.11327066], [0.0857183,  0.19072289]
        self.weird(p, q)


class Hilbert_Ball_Around_Points_and_Line(TestCase):
    def weird(self, p, q):
        o = omega.Omega(vertices=[(1.4, 0.0),
                                  (1.1, 0.8),
                                  (0.6, 1.1),
                                  (0.0, 1.2),
                                  (-.5, 0.8),
                                  (-.8, 0.3),
                                  (-.7, -.2)])
        """
        l connects p and q as {A--p--q--B}
        fit is the "perpendicular" line, dividing p-q through the vanishing point.
        fit takes in a radius, which is going to be half the distance from p to q. 
        """
        print(p, q)
        l = line.Line(p, q, o)

        tools.plot_congruent(plt, o.vertices, color='gray')
        tools.annotate(plt, range(len(o.vertices)), o.vertices)

        tools.plot_line(plt, *l.get_boundaries()[0], color='green')
        tools.plot_line(plt, *l.get_boundaries()[1], color='purple')
        div = l.get_best_dividing_line()
        fit = line.Line(div[0], div[1], o)
        plt.axline(fit.p, fit.q)
        radius = l.get_hdist() / 2
        ball = fit.hilbert_ball_about_line(radius)
        tools.plot_congruent(plt, ball.vertices, color="orange")
        hb1 = o.hilbert_ball_around_point(p, radius)
        hb2 = o.hilbert_ball_around_point(q, radius)
        tools.plot_congruent(plt, hb1.vertices, color='red')
        tools.plot_congruent(plt, hb2.vertices, color='red')
        tools.annotate(plt, 'pqab', [p, q, l.get_boundaries()[0][0], l.get_boundaries()[1][0]])
        tools.scatter(plt, [p, q, l.get_midpoint()])
        plt.show()

    def test_random_weird(self):
        p, q = [x[:2] for x in np.random.dirichlet((1, 1, 1), 2)]
        self.weird(p, q)

    def simplex(self, p, q):
        o = omega.Omega()
        """
        l connects p and q as {A--p--q--B}
        fit is the "perpendicular" line, dividing p-q through the vanishing point.
        fit takes in a radius, which is going to be half the distance from p to q. 
        """
        print(p, q)
        l = line.Line(p, q, o)

        tools.plot_congruent(plt, o.vertices, color='gray')
        tools.annotate(plt, range(len(o.vertices)), o.vertices)
        tools.plot_line(plt, *l.get_boundaries()[0], color='green')
        tools.plot_line(plt, *l.get_boundaries()[1], color='purple')
        div = l.get_best_dividing_line()
        fit = line.Line(div[0], div[1], o)
        plt.axline(fit.p, fit.q)
        radius = l.get_hdist() / 2
        ball = fit.hilbert_ball_about_line(radius)
        tools.plot_congruent(plt, ball.vertices, color="orange")
        hb1 = o.hilbert_ball_around_point(p, radius)
        hb2 = o.hilbert_ball_around_point(q, radius)
        tools.plot_congruent(plt, hb1.vertices, color='red')
        tools.plot_congruent(plt, hb2.vertices, color='red')
        tools.annotate(plt, 'pqab', [p, q, l.get_boundaries()[0][0], l.get_boundaries()[1][0]])
        tools.scatter(plt, [p, q, l.get_midpoint()])
        plt.show()

    def test_random_simplex(self):
        p, q = [x[:2] for x in np.random.dirichlet((1, 1, 1), 2)]
        self.simplex(p, q)



if __name__ == "__main__":
    Hilbert_Ball_2pt().test_broken_weird_b()

