from unittest import TestCase
from src.hilbert import omega, line
from src.misc import tools
import numpy as np
from matplotlib import pyplot as plt

class Test_Simplex(TestCase):

    def test_horizontal(self):
        o = omega.Omega()
        l = line.Line(np.array([0.1, 0.4]), np.array([0.3, 0.4]), o)
        plt.axline(l.p, l.q)
        tools.plot_congruent(plt, o.vertices)
        hilbert_ball = l.hilbert_ball_about_line(1)
        tools.plot_congruent(plt, hilbert_ball)
        plt.show()

    def test_vertical(self):
        o = omega.Omega()
        l = line.Line(np.array([0.4, 0.1]), np.array([0.4, 0.3]), o)
        plt.axline(l.p, l.q)
        tools.plot_congruent(plt, o.vertices)
        hilbert_ball = l.hilbert_ball_about_line(1)
        tools.plot_congruent(plt, hilbert_ball)
        plt.show()

    def test_diag(self):
        o = omega.Omega()
        l = line.Line(np.array([0.4, 0.1]), np.array([0.5, 0.2]), o)
        plt.axline(l.p, l.q)
        tools.plot_congruent(plt, o.vertices)
        hilbert_ball = l.hilbert_ball_about_line(1)
        tools.plot_congruent(plt, hilbert_ball)
        plt.show()



class Test_Weird(TestCase):
    o = omega.Omega(vertices=[(1.4, 0.0),
                                    (1.1, 0.8),
                                    (0.6, 1.1),
                                    (0.0, 1.2),
                                    (-.5, 0.8),
                                    (-.8, 0.3),
                                    (-.7, -.2)])

    def test_horizontal(self):
        l = line.Line(np.array([0.1, 0.4]), np.array([0.3, 0.4]), self.o)
        plt.scatter(*l.p)
        plt.scatter(*l.q)
        for bs in l.get_ball_spokes(plt):
            a, b = bs[1], bs[2]
            tools.plot_line(plt, a, b, color='gray')
        plt.axline(l.p, l.q, color='orange')
        tools.plot_congruent(plt, self.o.vertices, color='blue')
        plt.show()
        tools.plot_congruent(plt, self.o.vertices, color='blue')
        hilbert_ball = l.hilbert_ball_about_line(1, plt)
        tools.plot_congruent(plt, hilbert_ball)
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
        tools.plot_congruent(plt, self.o.vertices)
        plt.scatter(*l.p)
        plt.scatter(*l.q)
        plt.show()
        plt.axline(l.p, l.q)
        tools.plot_congruent(plt, self.o.vertices)
        hilbert_ball = l.hilbert_ball_about_line(1)
        tools.plot_congruent(plt, hilbert_ball)
        plt.show()

    def test_diag(self):
        l = line.Line(np.array([0.4, 0.1]), np.array([0.5, 0.2]), self.o)
        tools.plot_congruent(plt, self.o.vertices)
        plt.scatter(*l.p)
        plt.scatter(*l.q)
        for bs in l.get_ball_spokes():
            a, b = bs[1], bs[2]
            tools.plot_line(plt, a, b)
        plt.show()
        plt.axline(l.p, l.q)
        tools.plot_congruent(plt, self.o.vertices)
        hilbert_ball = l.hilbert_ball_about_line(1)
        tools.plot_congruent(plt, hilbert_ball)
        plt.show()