from misc import euclidean, tools, graham_scan
from hilbert import geometry
import numpy as np
from pprint import pprint
import operator


EPSILON = 0.0001
class Line:

    def __init__(self, p, q, omega):
        self.omega = omega
        self.p, self.q = np.array(p), np.array(q)
        self.l = self.p, self.q
        self.midpoint = None
        self.boundaries = None
        self.ball_spokes = None
        # between p and q
        self.hdist = None
        # boundary intersections
        self.A, self.B = None, None
        # Check if point way outside, set it to intersection with boundary...



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
        # the vertexes of the boundaries of Omega crossed by the line, ie, self.get_boundaries, can be combined to make
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
        :return: list of (tuple tuple tuple), ie (intersection, a, b)
        """

        def add_if_absent(intersection, a, b):
            key = tuple(sorted([tuple(a), tuple(b)]))
            if key not in seen:
                seen.add(key)
                connections.append([intersection, a, b])

        def connect_tangents(segments, tangent_points, plt=None):
            # if plt is not None:
            #    for i, v in enumerate(segments):
            #        plt.annotate(str(i), v)

            for i in range(len(segments) - 1):
                ii = i + 1
                segment = (segments[i], segments[ii])
                intersectpoint = euclidean.intersect(self.l, segment)
                #bestvertex, dot, index = tools.bs_on_points(
                #    tangent_points, lambda v: euclidean.tangent_dot(v, intersectpoint, self.q), operator.gt)
                # dots = np.apply_along_axis(lambda x: abs(euclidean.tangent_dot(x, intersectpoint, self.q)), 1, tangent_points)
                dots = np.apply_along_axis(lambda x: abs(euclidean.cos(x, intersectpoint, self.q)), 1, tangent_points)
                bestidx = np.argmin(dots)
                bestvertex = tangent_points[bestidx]

                if plt is not None:
                    tools.plot_line(plt, intersectpoint, segment[0], color='yellow')
                intersection_0 = euclidean.intersect((bestvertex, segment[0]), self.l)
                intersection_1 = euclidean.intersect((bestvertex, segment[1]), self.l)
                # Putting into set because might be repeats, dunno
                # All np.arrays must be tuples to be hashed.
                key = (tuple(bestvertex), tuple(segment[0]))
                add_if_absent(intersection_0, bestvertex, segment[0])
                add_if_absent(intersection_1, bestvertex, segment[1])

        seen = set()
        self.get_boundary_intersections()
        # BUG: With simplexes, A and B are correct. With the more complicated space, they should be swapped.
        #print("Is A in it's proper place?", euclidean.point_on_line(self.A, self.get_boundaries()[0]))
        #print("Is B in it's proper place?", euclidean.point_on_line(self.B, self.get_boundaries()[1]))

        connections = []

        # add_if_absent(self.A, *self.get_boundaries()[0])
        # add_if_absent(self.B, *self.get_boundaries()[1])

        # First divide vertices into two groups
        points_above_l_before, points_below_l_before, points_above_l_after, points_below_l_after  = [], [], [], []
        loadingpoints = [points_above_l_before, points_below_l_before]

        for v in self.omega.vertices:
            if euclidean.point_below_line(v, self.l):
                #points_below_l.append(v)
                loadingpoints[1].append(v)
                if len(loadingpoints[0]) > 0:
                    loadingpoints[0] = points_above_l_after
            else:
                #points_above_l.append(v)
                loadingpoints[0].append(v)
                if len(loadingpoints[1]) > 0:
                    loadingpoints[1] = points_below_l_after

        points_above_l = points_above_l_after + points_above_l_before
        points_below_l = points_below_l_after + points_below_l_before
        # points_below_l = points_below_l[1:] + [points_below_l[0]]
        connect_tangents(points_above_l, points_below_l, plt)
        connect_tangents(points_below_l, points_above_l, plt)
        self.ball_spokes = connections
        #print("BALL SPOKES, len", len(self.ball_spokes), ":")
        #pprint(self.ball_spokes)
        return self.ball_spokes

    def hilbert_ball_about_line(self, radius, plt=None):
        ballpoints = []
        i=0
        # pprint(self.get_ball_spokes())
        for intersection, A, B in self.get_ball_spokes():
            p1 = geometry.hdist_to_euc(intersection, A, B, radius)
            ballpoints.append(p1)
            tools.plot_line(plt,A, B, color='red')
            p2 = geometry.hdist_to_euc(np.array(intersection), np.array(A), np.array(B), -radius)
            ballpoints.append(p2)
            if plt:
                plt.scatter(p1[0], p1[1])
                plt.scatter(p2[0], p2[1])
                plt.annotate(str(i), p1)
                plt.annotate(str(-i), p2)
                i+=1

        ballpoints = [x for x in ballpoints if not np.isnan(x).any()]
        gs = graham_scan.graham_scan(ballpoints)
        return gs


