import unittest
import random
from misc import tools
from matplotlib import pyplot as plt

class MyTestCase(unittest.TestCase):
    def test_something(self):
        p, q, r = tools.rand_points(3)
        qp, rp = q - p, r - p
        points = []
        # Generate 20 points
        for i in range(20):
            a, b = random.random(), random.random()
            if a + b > 1:
                a, b = 1-a, 1-b
            points.append(p + qp * a + rp * b)

        tools.plot_congruent(plt, [p, q, r], color='gray')
        tools.scatter(plt, points)
        plt.show()
