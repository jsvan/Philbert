import numpy as np

from hilbert import line, polygon
from hilbert.geometry import hdist_to_euc
from misc import tools, euclidean
import itertools


class Omega:

    def __init__(self, vertices=([0, 0], [0, 1], [1, 0]), offset=[0, 0] ):

        self.polygon = polygon.Polygon(vertices, offset)
        self.vertices = self.polygon.vertices
        # This is a hack for binary search, wraps around the final points. Whatever.
        self.vertices_expanded = self.polygon.vertices_expanded

    def spokes(self, p):
        """
        :param p: np array point in omega
        :return: list of tuples: (coords of the omega vertex, lambda equation for related spoke)
        """
        return [(v, lambda t: p + t * (p - np.array(v))) for v in self.vertices]

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


    def all_tangent_lines(self, poly1, poly2, plt=None):
        maxwards = [True, False]
        tangents = []
        for mws in itertools.product(maxwards, maxwards):
            try:
                tangents.append(self.tangent_line(poly1, poly2, mws[0], mws[1], plt))
            except:
                continue

        return tangents



    def tangent_line(self, poly1, poly2, maxwards1, maxwards2, plt=None):

        removedpoly1, removedpoly2 = None, None  # on;ly one p[oint for each is necessary for orientaion test
        newpoly1, newpoly2 = [], []
        # This is a hack, remove overlapping points
        for n in poly1.vertices:
            if poly2.contains(n):
                removedpoly1 = n
            else:
                newpoly1.append(n)
        for n in poly2.vertices:
            if poly1.contains(n):
                removedpoly2 = n
            else:
                newpoly2.append(n)

        poly1 = polygon.Polygon(newpoly1)
        poly2 = polygon.Polygon(newpoly2)
        prev_v1, prev_v2 = None, None
        v1, v2 = poly1.vertices[0], None
        i = 0
        while not (euclidean.eq(v1, prev_v1) and euclidean.eq(v2, prev_v2)):
            prev_v1 = v1
            prev_v2 = v2
            v2info = poly2.point_tangent_linear(v1, maxwards2)
            v2 = v2info['v']
            v1info = poly1.point_tangent_linear(v2, maxwards1)
            v1 = v1info['v']
            if plt is not None:
                tools.plot_line(plt, v1, v2)
                tools.annotate(plt, [str(i)], [(v2 + v1) / 2])
            i += 1
        #if euclidean.orient(v1)
        return v1info, v2info



