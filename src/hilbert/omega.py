import numpy as np

from hilbert import line, polygon
from hilbert.geometry import hdist_to_euc
from misc import tools, euclidean
from misc.graham_scan import graham_scan
import itertools


class Omega:

    def __init__(self, vertices=([0, 0], [0, 1], [1, 0]) ):
        self.vertices = graham_scan([np.array(v) for v in vertices])
        # This is a hack for binary search, wraps around the final points. Whatever.
        self.vertices_expanded = [self.vertices[-1]] + self.vertices + [self.vertices[0]]

    def spokes(self, p):
        """
        :param p: np array point in omega
        :return: list of tuples: (coords of the omega vertex, lambda equation for related spoke)
        """
        return [(v, lambda t: p + t * (p - np.array(v))) for v in self.vertices]

    def line_boundaries(self, p, q):
        """
        Run orientation tests to find intersections with boundaries. Log time search.
        :param p:
        :param q:
        :return: two lines, each of two points.
        """
        return self.find_boundary(p, q), self.find_boundary(q, p)

    def find_boundary(self, p, q):
        v = self.vertices_expanded
        numvert = len(v)
        binarysearch = tools.BinarySearcher(0, numvert, discrete=True)
        attempts = numvert + 1

        while attempts > 0:
            attempts -= 1
            leftidx = binarysearch.next()
            leftv, rightv = v[leftidx], v[(leftidx + 1) % numvert]

            if euclidean.orient(p, q, leftv) == euclidean.COUNTER_CW:  # BAD ORIENTATION
                binarysearch.feedback(higher=False)
                continue
            if euclidean.orient(p, q, rightv) == euclidean.CLOCKWISE:
                binarysearch.feedback(higher=True)
                continue
            return (leftv, rightv)
        raise Exception("Finding boundary with binary search failed")


    def hilbert_ball_around_point(self, p, r):
       p = np.array(p)
       spokes = self.spokes(p)
       # I think it's best to find the q point on every spoke, instead of connecting to vanishing points.
       # Vanishing points require binary search probably for every spoke intersection (on opposite side).
       # That may not actually be slower, but this is simpler in code for sure.
       points = []
       for v, _ in spokes:
          l = line.Line(p, v, self)
          intersections = l.get_boundary_intersections()
          #print('bbb', boundaries)
          #intersections = ordered_double_intersect(p, v, boundaries)
          points.append(hdist_to_euc(p, *intersections, r))
          points.append(hdist_to_euc(p, *intersections, -r))
       return polygon.Polygon(points)


    def tangent_lines(self, poly1, poly2, plt=None):
        maxwards = [True, False]
        tangents = {}
        for mws in itertools.product(maxwards, maxwards):
            tangents[mws] = self.tangent_line(poly1, poly2, mws[0], mws[1], plt)

        return tangents



    def tangent_line(self, poly1, poly2, maxwards1, maxwards2, plt=None):

        prev_v1, prev_v2 = None, None
        v1, v2 = poly1.vertices[0], None
        i = 0
        while not (euclidean.eq(v1, prev_v1) and euclidean.eq(v2, prev_v2)):
            prev_v1 = v1
            prev_v2 = v2
            v2 = poly2.point_tangent_linear(v1, maxwards2)['v']

            v1 = poly1.point_tangent_linear(v2, maxwards1)['v']

            if plt is not None:
                tools.plot_line(plt, v1, v2)
                tools.annotate(plt, [str(i)], [(v2 + v1) / 2])
            i += 1
        return v1, v2
