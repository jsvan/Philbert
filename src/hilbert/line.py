from misc import euclidean, tools, graham_scan
from hilbert import geometry
import numpy as np
import operator


EPSILON = 0.0001
class Line:

    def __init__(self, p, q, omega):
        self.p, self.q = p, q
        self.l = p, q
        self.omega = omega
        self.midpoint = None
        self.boundaries = None
        # between p and q
        self.hdist = None
        # boundary intersections
        self.A, self.B = None, None
        self.ball_spokes = None

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
        vanishing_point = euclidean.intersect(*self.get_boundaries())
        return [self.get_midpoint(), vanishing_point]

    def get_ball_spokes(self, plt=None):
        """
        :return: list of (tuple tuple tuple), ie (intersection, a, b)
        """

        def connect_tangents(segments, tangent_points, plt=None):
            print("Segments:", segments)
            print("Tangent possibilities: ", tangent_points)
            if plt is not None:
                for i, v in enumerate(segments):
                    plt.annotate(str(i), v)

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
                print("dots", dots)
                print("intersection:", intersectpoint)
                print("best idx:", bestidx, ", best vertex", bestvertex)

                if plt is not None:
                    tools.plot_line(plt, intersectpoint, segment[0], color='yellow')
                intersection_0 = euclidean.intersect((bestvertex, segment[0]), self.l)
                intersection_1 = euclidean.intersect((bestvertex, segment[1]), self.l)
                # Putting into set because might be repeats, dunno
                # All np.arrays must be tuples to be hashed.
                connections.append([intersection_0, bestvertex, segment[0]])
                connections.append([intersection_1, bestvertex, segment[1]])

        connections = [[self.A, *self.get_boundaries()[0]],
                        [self.B, *self.get_boundaries()[1]]]

        self.get_boundary_intersections()


        # First divide vertices into two groups
        points_above_l, points_below_l = [], []
        for v in self.omega.vertices:
            if euclidean.point_below_line(v, self.l):
                points_below_l.append(v)
            else:
                points_above_l.append(v)
        points_below_l = points_below_l[1:] + [points_below_l[0]]
        print(points_below_l)
        connect_tangents(points_above_l, points_below_l, plt)
        connect_tangents(points_below_l, points_above_l, plt)
        self.ball_spokes = connections
        return self.ball_spokes

    def hilbert_ball_about_line(self, radius, plt=None):
        ballpoints = []
        i=0
        print(self.get_ball_spokes())
        for intersection, A, B in self.get_ball_spokes():
            p1 = geometry.hdist_to_euc(np.array(intersection), np.array(A), np.array(B), radius)
            ballpoints.append(p1)
            p2 = geometry.hdist_to_euc(np.array(intersection), np.array(A), np.array(B), -radius)
            ballpoints.append(p2)
            if plt:
                plt.scatter(p1[0], p1[1])
                plt.scatter(p2[0], p2[1])
                plt.annotate('*'+str(i), p1)
                plt.annotate('*'+str(-i), p2)
                i+=1

        ballpoints = [x for x in ballpoints if not np.isnan(x).any()]
        gs = graham_scan.graham_scan(ballpoints)
        return gs
