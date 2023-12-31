from misc import euclidean, tools
from hilbert import geometry, polygon
from misc.euclidean import Point, BasicLine, Vertex
import numpy as np
from collections import defaultdict
import random

EPSILON = 0.0001

class Line:

    def __init__(self, p, q, omega):
        self.omega = omega
        self.p, self.q = Point(p), Point(q)
        self.l = euclidean.BasicLine(p, q)
        self.midpoint = None
        self.boundaries = None
        self.ball_spokes = None
        # between p and q
        self.hdist = None
        # boundary intersections
        self.A, self.B = None, None
        # Check if point way outside, set it to intersection with boundary...


    def nearest_point(self, p, plt=None):
        """
        The following is the binary search version, which i am not going to use because ordering the points is difficult

        intersectionpoints = [euclidean.intersect((p, v), (self.p, self.q)) for v in self.omega.vertices ]
        def transform(x):
            l = Line(x, self.p, self.omega)
            A, B = l.get_boundary_intersections()
            return geometry.dist(x, p, A, B)
        return tools.bs_on_points(intersectionpoints, transform, operator.lt)
        """
        intersectionpoints = [euclidean.intersect((p, v), (self.p, self.q)) for v in self.omega.vertices ]
        def transform(x):
            """
            x is a point
            output is the distance from p to x. If outside the OMEGA, infinity.
            """
            #try:
            l = Line(x, p, self.omega)
            A, B = l.get_boundary_intersections()
            return geometry.dist(x, p, A, B) if euclidean.point_within_region(x, *self.get_boundary_intersections()) else np.inf
            #except:
            #    return np.inf
        try:
            distances = [transform(x) for x in intersectionpoints]
        except Exception as e:
            print(e, e.__doc__)
            return None, None, 10
        if plt is not None:
            tools.scatter(plt, [x for i, x in enumerate(intersectionpoints) if not np.isinf(distances[i])])

        minidx = np.argmin(distances)
        return intersectionpoints[minidx], distances[minidx], minidx


    # def dist_to_p(self, p):


    def get_boundaries(self):
        if self.boundaries is None:
            self.boundaries = self.omega.line_boundaries(self.p, self.q)
        return self.boundaries

    def get_hdist(self):
        """
        hdist between two points p and q
        :return: float
        """
        if self.hdist is None:
            self.hdist = geometry.dist(self.p, self.q, *self.get_boundary_intersections())
        return self.hdist

    # gets hilbertian midpoint between p and q
    def get_midpoint(self):
        if self.midpoint is None:
            self.midpoint = geometry.hdist_to_euc(self.p, *self.get_boundary_intersections(),
                                                           self.get_hdist() / 2)
        return self.midpoint

    def get_boundary_intersections(self):
        if self.A is None or self.B is None:
            self.A, self.B = geometry.ordered_double_intersect(self.p, self.q, self.get_boundaries())
        return self.A, self.B

    def get_best_dividing_line(self):
        ba, bb = self.get_boundaries()
        vanishing_point = euclidean.intersect(ba, bb)
        other_point = vanishing_point
        # Other point is unnecessary, can use either point, somehow it all works.
        # But I do this point because the pyplots lay better.
        # This point takes the vanishing point and maps it along its line to a point inside omega.
        # The vertexes of the boundaries of Omega crossed by the line, ie, self.get_boundaries, can be combined to make
        # spokes that are within the convex shape, or on its edges.
        # Moving Epsilon inside could be hard because we'd have to move epsilon in the right direction. However, if
        # we draw a spoke from the Euclidean midpoints of those boundaries, those should always be inside the convex
        # shape. Take the intersection of the midpoints of the boundaries with the line from MIDPOINT to VANISHING.
        # That intersection will be on the line, and inside the convex shape.
        # That way, if these points are placed back into Omega, they will be within the Omega.
        bam = (ba[0] + ba[1]) / 2
        bbm = (bb[0] + bb[1]) / 2
        other_point = euclidean.intersect([self.get_midpoint(), vanishing_point], [bam, bbm])
        return [self.get_midpoint(), other_point]

    def get_ball_spokes(self, plt=None):
        """
        Constructs the ball spokes for the hilbert ball around this line
        :return: list of (tuple tuple tuple), ie (intersection, a, b)
        Currently runs in O(n) time, can be shifted to log(n), but that had some bugs
        """
        if self.ball_spokes is not None:
            return self.ball_spokes

        def add_if_absent(intersection, a, b):
            key = tuple(sorted([tuple(a.v), tuple(b.v)]))
            if key not in seen:
                seen.add(key)
                connections.append([intersection, a, b])

        def connect_tangents(segments, tangent_points, plt=None):
            #tangent_vs = np.array([x.v for x in tangent_points])
            for i in range(len(segments) - 1):
                ii = i + 1
                segment = BasicLine(segments[i], segments[ii])
                intersectpoint = euclidean.intersect(self.l, segment)
                #bestvertex, dot, index = tools.bs_on_points(
                #    tangent_points, lambda v: euclidean.tangent_dot(v, intersectpoint, self.q), operator.gt)
                # dots = np.apply_along_axis(lambda x: abs(euclidean.tangent_dot(x, intersectpoint, self.q)), 1, tangent_points)
                dots = np.apply_along_axis(lambda x: abs(euclidean.cos(x, intersectpoint, self.q)), 1, tangent_points)
                bestidx = np.argmin(dots)
                bestvertex = tangent_points[bestidx]

                if plt is not None:
                    tools.plot_line(plt, intersectpoint, segment.a, color='yellow')
                intersection_0 = euclidean.intersect(BasicLine(bestvertex, segment.a), self.l)
                intersection_1 = euclidean.intersect(BasicLine(bestvertex, segment.b), self.l)
                # Putting into set because might be repeats, dunno
                # All np.arrays must be tuples to be hashed.
                add_if_absent(intersection_0, bestvertex, segment.a)
                add_if_absent(intersection_1, bestvertex, segment.b)

        seen = set()
        self.get_boundary_intersections()
        connections = []
        # First divide vertices into two groups
        # separate omega
        points_above_l, points_below_l = self.omega.halves(self)
        connect_tangents(points_above_l, points_below_l, plt)
        connect_tangents(points_below_l, points_above_l, plt)
        self.ball_spokes = connections
        return self.ball_spokes

    def hilbert_ball_about_line(self, radius, plt=None):
        ballpoints = []
        i=0
        for intersection, A, B in self.get_ball_spokes():
            p1 = geometry.hdist_to_euc(intersection, A, B, radius)
            ballpoints.append(p1)
            p2 = geometry.hdist_to_euc(intersection, A, B, -radius)
            ballpoints.append(p2)
            if plt:
                tools.plot_line(plt, A, B, color='red')
                plt.scatter(p1[0], p1[1])
                plt.scatter(p2[0], p2[1])
                plt.annotate(str(i), p1)
                plt.annotate(str(-i), p2)
                i += 1

        ballpoints = [x for x in ballpoints if not np.isnan(x).any()]
        gs = polygon.Polygon(ballpoints)
        return gs


    def sample_points_outside_ball(self, n, hilbert_ball=None, radius=None):
        """
        Returns n sampled points from two classes, for each side of the line, which lie OUTSIDE of the given hilbert
        ball. If a hilbertball around the line is NOT given, radius must be included and the ball will be made
        automatically.

        Slow and maybe not uniform...
        but good enough for research
        """
        if hilbert_ball is None:
            if radius is None:
                raise Exception("Radius is None and that is bad. Radius may not be None if hilbert_ball is also None.")
            hilbert_ball = self.hilbert_ball_about_line(radius)

        yooksball, zooksball = hilbert_ball.halves(self)
        yooksomeg, zooksomeg = self.omega.halves(self)
        # Getting ball spokes is hard, but it will help triangulate the space better
        spokemap = defaultdict(int)
        # Using spokes mihgt not be necessary, but I know at least with them
        # the polygon is divided into 3 and 4 point polygons.
        # There are (ball vertices >=  omega) # of vertices
        for s in self.get_ball_spokes():
            spokemap[tuple(s[1].v)] += 1
            spokemap[tuple(s[2].v)] += 1

        #connect omega's vertices to polygons such they are triangulated:
        def triangulate(omg, poly):
            # In between every omega vertex is a square. We're gonna connect vertices to nexts
            # rather than previouses
            polyi = 0
            prevo = None
            triangles = []
            for o in omg:
                if prevo is not None:
                    triangles.append([o, prevo, poly[polyi]])
                for _ in range(spokemap[tuple(o.v)]):
                    if len(poly) > polyi + 1:
                        triangles.append([o, poly[polyi], poly[polyi + 1]])
                    polyi += 1
                prevo = o
            return triangles

        def softmaxify(triangles):
            areas = [euclidean.triangle_area(*t) for t in triangles]
            tot = sum(areas)
            return [a / tot for a in areas]

        def sample_triangle(which, triangles):
            # Always calculating qp, rp is wasteful, but if we're only sampling 20 points who cares. If we
            # were to do this 4000 times, then change this
            p, q, r = triangles[which]
            qp, rp = q - p, r - p
            a, b = random.random(), random.random()
            if a + b > 1:
                a, b = 1 - a, 1 - b
            return p + a * qp + b * rp

        bluetriangles = triangulate(yooksomeg, yooksball)
        redtriangles = triangulate(zooksomeg, zooksball)
        blueareas = softmaxify(bluetriangles)
        redareas = softmaxify(redtriangles)
        blues = [sample_triangle(x, bluetriangles) for x in random.choices(list(range(len(blueareas))),
                                     blueareas, k=n)]
        reds = [sample_triangle(x, redtriangles) for x in random.choices(list(range(len(redareas))),
                                                                                 redareas, k=n)]
        return blues, reds




