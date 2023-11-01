import numpy as np
from misc import euclidean, tools
from hilbert import line
from misc.graham_scan import graham_scan
from misc.euclidean import Vertex, Edge, Point

class Polygon:

    def __init__(self, vertices, offset=[0, 0]):
        self.offset = np.array(offset)
        self.vertices = [np.array(v) + self.offset for v in vertices]
        self.vertices = [Vertex(v, i) for i, v in enumerate(graham_scan(self.vertices))]
        #self.vertices_expanded = [self.vertices[-1]] + self.vertices + [self.vertices[0]]

    def v(self, i):
        return self.vertices[self.i(i)]

    def i(self, i):
        return i % len(self.vertices)

    def point_tangents(self, point):
        tangents = []
        try:
            tangents.append(self.point_tangent_linear(point, maxwards=True))
        except Exception as e:
            print(e, e.__doc__)
        try:
            tangents.append(self.point_tangent_linear(point, maxwards=False))
        except Exception as e:
            print(e, e.__doc__)
        return tangents

    def point_tangent_linear(self, point, maxwards=True):
        """
        linear time
        tangent line from point
        """
        correct_orientation, wrong_orientation = euclidean.COUNTER_CW, euclidean.CLOCKWISE
        if maxwards:
            correct_orientation, wrong_orientation = wrong_orientation, correct_orientation

        for vi, v in enumerate(self.vertices):
            # v-minus, v-plus
            v, vm, vp = self.v(vi), self.v(vi - 1), self.v(vi + 1)
            minus_orientation = euclidean.orient(point, v, vm)
            plus_orientation = euclidean.orient(point, v, vp)
            if minus_orientation == correct_orientation and plus_orientation == correct_orientation:
                return v
        raise Exception("Boundary failed, point potentially inside.")


    def point_tangent_collision_detection(self, point, maxwards=True):
        """
        tangent lines between two polygons that ignores collided points
        """
        correct_orientation, wrong_orientation = euclidean.COUNTER_CW, euclidean.CLOCKWISE
        if maxwards:
            correct_orientation, wrong_orientation = wrong_orientation, correct_orientation
        vs = self.vertices
        numvert = len(vs)

        for vi, v in enumerate(vs):
            # v-minus, v-plus
            v, vm, vp = vs[vi], vs[(vi - 1) % numvert], vs[(vi + 1) % numvert]
            minus_orientation = euclidean.orient(point, v, vm)
            plus_orientation = euclidean.orient(point, v, vp)
            if minus_orientation == correct_orientation and plus_orientation == correct_orientation:
                return Vertex(Point(v), vi)
        raise Exception("Boundary failed, point potentially inside.")


    def contains(self, point):
        """
        point in polygon?
        """
        v = self.vertices
        for i, ii in tools.i_ii(len(v)):
            if euclidean.orient(point, v[i], v[ii]) == euclidean.COUNTER_CW:
                return False
        return True


    def point_tangent_collision_detection_old(self, point, pm, pp, maxwards=True):
        """
        Finds the tangent from a pointA to a polygonB (self). The point is not just a point, but the connecting edges of its own
        polygonA. Does quick_is_outside to ensure the tangent doesn't connect to a point on the inside of the polygon in
        O(1) time.
        Runs in linear time.
        point : pointA
        pm : point minus 1 of polygonA
        pp : point plus 1
        """
        correct_orientation, wrong_orientation = euclidean.COUNTER_CW, euclidean.CLOCKWISE
        if maxwards:
            correct_orientation, wrong_orientation = wrong_orientation, correct_orientation

        vs = self.vertices
        numvert = len(vs)

        for vi, v in enumerate(vs):
            # Check that vertex is outside of point's polygon
            if not euclidean.quick_is_outside(vs[vi], pm, point, pp):
                continue
            # v-minus, v-plus
            v, vm, vp = vs[vi], vs[(vi - 1) % numvert], vs[(vi + 1) % numvert]
            minus_orientation = euclidean.orient(point, v, vm)
            plus_orientation = euclidean.orient(point, v, vp)
            if minus_orientation == correct_orientation and plus_orientation == correct_orientation:
                return Vertex(v, vi)
        raise Exception("Boundary failed, point potentially inside.")

    def other_tangent(self, point, i):
        """
        There are two possible tangent lines in two dimensions, from a point to a convex polygon.
        Give the index of the edge which you want to AVOID, this method will give you the tangent line that is NOT
        on that edge.
        """
        t = None
        try:
            t = self.point_tangent_linear(point, True)
        except:
            pass
        try:
            # on the same line
            if t is None or euclidean.eq(t, self.v(i)) or euclidean.eq(t, self.v(i + 1)):
                t = self.point_tangent_linear(point, False)
        except:
            pass
        if t is None:
            print(f"No tangent found for {point} of type {type(point)} with shape {point.shape}, {i}")
        return t

    def point_tangent_binary(self, point, maxwards=True):
        # BROKEN, needs to be debugged
        # Think I have to do orient tests
        # Yeah, confirmed
        # https://www.cs.umd.edu/class/fall2021/cmsc754/Lects/cmsc754-fall-2021-lects.pdf page 15
        # Given point outside of polygon, sample random points within polygon, continue until
        correct_orientation, wrong_orientation = euclidean.COUNTER_CW, euclidean.CLOCKWISE
        if maxwards:
            correct_orientation, wrong_orientation = wrong_orientation, correct_orientation

        vs = self.vertices_expanded
        numvert = len(vs)
        binarysearch = tools.BinarySearcher(0, numvert, discrete=True)
        attempts = numvert + 1
        print("STARTING")
        while attempts > 0:
            attempts -= 1
            vi = binarysearch.next() % numvert
            vim = (vi - 1) % numvert
            vip = (vi + 1) % numvert
            # v-minus, v-plus
            v, vm, vp = vs[vi], vs[vim], vs[vip]
            print(f"vertex is {vi}, vim is {vim}, vip is {vip}")
            print(f"p {point}, vm {vm}, v {v}, vp {vp}")
            minus_orientation = euclidean.orient(point, v, vm)
            plus_orientation = euclidean.orient(point, v, vp)
            print(f"i-1 is {euclidean.dirname(minus_orientation)}, f+1 is {euclidean.dirname(plus_orientation)}")
            if minus_orientation == wrong_orientation:
                binarysearch.feedback(higher=maxwards)
                print(
                    f"i-1 is {euclidean.dirname(minus_orientation)}, so going {'higher' if maxwards else 'lower'}.")

                continue
            if plus_orientation == wrong_orientation:
                # Too far
                binarysearch.feedback(higher=not maxwards)
                print(
                    f"i+1 is {euclidean.dirname(plus_orientation)}, so going {'higher' if not maxwards else 'lower'}.")

                continue
            return v
        raise Exception("Finding boundary with binary search failed")

    def halves(self, l, index=False):
        """
        l is a line.Line, or two 2d numpy arrays like [np.a([x, y]), np.a([x, y])]
        returns two arrays of vertices of the polygon, for each side of line they're on
        """
        if isinstance(l, line.Line):
            l = l.l


        points_above_l_before, points_below_l_before, points_above_l_after, points_below_l_after  = [], [], [], []
        loadingpoints = [points_above_l_before, points_below_l_before]

        for i, v in enumerate(self.vertices):
            elem = i if index else v
            if euclidean.point_below_line(v, l):
                loadingpoints[1].append(elem)
                if len(loadingpoints[0]) > 0:
                    loadingpoints[0] = points_above_l_after
            else:
                loadingpoints[0].append(elem)
                if len(loadingpoints[1]) > 0:
                    loadingpoints[1] = points_below_l_after

        points_above_l = points_above_l_after + points_above_l_before
        points_below_l = points_below_l_after + points_below_l_before
        return points_above_l, points_below_l


    def find_boundary(self, p, q):
        v = self.vertices
        numvert = len(v)
        for i in range(numvert+1):
            lv = self.v(i)
            uv = self.v(i+1)
            if euclidean.orient(p, q, lv) != euclidean.CLOCKWISE and euclidean.orient(p, q, uv) != euclidean.COUNTER_CW:
                return Edge(lv, uv)
        raise Exception("Finding boundary with sequential search failed")


    def find_boundary_binarysearch(self, p, q):
        v = self.vertices_expanded
        numvert = len(v)
        binarysearch = tools.BinarySearcher(0, numvert, discrete=True)

        while binarysearch.has_next():
            leftidx = binarysearch.next()
            print(leftidx)
            leftv, rightv = v[leftidx % numvert], v[(leftidx + 1) % numvert]

            if euclidean.orient(p, q, leftv) == euclidean.COUNTER_CW:  # BAD ORIENTATION
                binarysearch.feedback(higher=False)
                continue
            if euclidean.orient(p, q, rightv) == euclidean.CLOCKWISE:
                binarysearch.feedback(higher=True)
                continue
            return Edge(Vertex(leftv, leftidx % numvert), Vertex(rightv, (leftidx + 1) % numvert)) #(leftv, rightv)

        raise Exception("Finding boundary with binary search failed")

    """
    for i in range(len(self.vertices)):
        vm, vp = self.v(i).v, self.v(i+1).v
        if euclidean.orient(p, q, vm) == euclidean.CLOCKWISE and \
                euclidean.orient(p, q, vp) == euclidean.COUNTER_CW:
            return (vm, vp)
    raise Exception("FINDING BOUNDARY FAILED")
    """

    def __len__(self):
        return len(self.vertices)

    def line_boundaries(self, p, q):
        """
        Run orientation tests to find intersections with boundaries. Log time search.
        :param p:
        :param q:
        :return: two Edges
        """
        return self.find_boundary(p, q), self.find_boundary(q, p)

    def center(self):
        x = sum([v[0] for v in self.vertices]) / len(self.vertices)
        y = sum([v[1] for v in self.vertices]) / len(self.vertices)
        return euclidean.Point([x, y])

    def edges(self):
        # generator
        return (self.edge(i) for i in range(len(self)))

    def edge(self, i):
        return Edge(self.v(i), self.v(i + 1))

    def annotate_edges(self, plt):
        # Capitalized of the beforehand vertex
        tools.annotate(plt, (tools.namer(i).upper() for i in range(len(self))), self.vertices)

    def annotate_nodes(self, plt):
        # lowercase
        tools.annotate(plt, (tools.namer(i) for i in range(len(self))), (euclidean.avg_point(*e) for e in self.edges))

    def intersections(self, other, prev_i):
        pass

    def sample_from_ball_border(self, n, l):
        halves = self.halves(l)
        points = []
        dist = 0
        above = True

        while len(points) < n:
            try:
                p, newdist = euclidean.uniform_sample_from_line_segments(halves[above], dist)
                above = not above
                points.append(p)
            except Exception as e:
                print(e)
                points = []
                dist = 0
                above = True
        return points

    """
    Returns the points on a uniform grid of the space inside the polygon.
    input: width is the vertical and horizontal distance between points. 
    
    Algorithm:
        find Y upper and lower limits
        then create horizontal lines and intersect them on the edges
        The first line will be order edges time
        every other line will look at the two neighbors 
        Once you get the intersections, you can generate the points within.
    """
    def gridspace(self, width, plt=None):
        X, Y = euclidean.X, euclidean.Y
        tallest  = max(self.vertices, key=lambda v:v[Y])[Y]
        shortest = min(self.vertices, key=lambda v:v[Y])[Y]

        gridYs = tools.linspace(shortest, tallest, width)
        gridY = gridYs[0]
        bottomLine = euclidean.BasicLine(euclidean.Point((0, gridY)),
                                         euclidean.Point((1, gridY)))
        leftEdge, rightEdge = self.line_boundaries(*bottomLine)
        bottomIdxs = leftEdge.point_a.i, rightEdge.point_a.i
        leftintersect = euclidean.intersect(leftEdge, bottomLine)
        rightintersect = euclidean.intersect(rightEdge, bottomLine)

        # Adding the first items to start off with...
        points = [Point((x, gridY)) for x in tools.linspace(leftintersect[X], rightintersect[X], width)]

        for gridY in gridYs[1:]:
            # Check to ensure correct edges, or iterate.
            # If we're above, move left edge upwards (positive dir)
            # and move right edge upwards (negative dir)
            if gridY > leftEdge.point_b.v[Y]:
                i = leftEdge.point_a.i
                leftEdge = self.edge(i + 1)
            if gridY > rightEdge.point_a.v[Y]:
                i = rightEdge.point_a.i
                rightEdge = self.edge(i - 1)

            bottomLine = euclidean.BasicLine(euclidean.Point((0, gridY)),
                                             euclidean.Point((1, gridY)))
            #print(bottomLine)
            leftintersect = euclidean.intersect(leftEdge, bottomLine)
            rightintersect = euclidean.intersect(rightEdge, bottomLine)
            if leftintersect is None or rightintersect is None:
                continue

            points += [Point((x, gridY)) for x in tools.linspace(leftintersect[X], rightintersect[X], width)]

        if plt is not None:
            print("maxs", tallest, shortest)
            print("topidx", bottomIdxs)
            print("grid", points)
            print("omega", self.vertices)
            tools.plot_congruent(plt, self.vertices, color='gray')
            tools.annotate(plt, range(len(self)), self.vertices)
            tools.plot_line(plt, *bottomLine, axline=True)
            tools.scatter(plt, points)
            plt.show()

        return points









"""

Recommended usage is from dual.py 's Sectors 's pseudohyperbola(), which returns two halves of a circle around a
polar dual line (a point in primal space)


"""

class HalvedPolygon:

    def __init__(self, upperhalf, lowerhalf, dividingline=None):
        self.dividingline = dividingline
        self.upperhalf = upperhalf
        self.lowerhalf = lowerhalf
        print("before", self.upperhalf)
        self.sort_points(self.upperhalf)
        print("after", self.upperhalf)
        self.sort_points(self.lowerhalf)


    def sort_points(self, points):
        if self.dividingline is not None:
            for p in points:
                p.t = euclidean.line_get_t(euclidean.closest_point_on_line(p, self.dividingline), self.dividingline)
            points.sort(key=lambda p: p.t)

    # the bespoke "...the fuck is efficiency?" algorithm
    def intersections(self, other, intersections=[]):
        count = len(intersections)
        for sa, sb in self.iter_all_edges():
            selfedge = euclidean.BasicLine(sa, sb)
            for oa, ob in other.iter_all_edges():
                otheredge = euclidean.BasicLine(oa, ob)
                if euclidean.line_intersects_on_segments(selfedge, otheredge):
                    intersections.append(euclidean.intersect(selfedge, otheredge))
        return intersections, len(intersections) - count


    """
        #For each half of our polygon
        for half_i, selfhalf in enumerate(self.halves):
            # Compare against every edge of the other polygon
            for v in range(len(selfhalf) - 1):
                selfedge = euclidean.BasicLine(selfhalf[v], selfhalf[v + 1])
                for j in range(len(other.upperhalf) - 1):
                    otheredge = euclidean.BasicLine(self.upperhalf[v], self.upperhalf[v + 1])
                    if euclidean.line_intersects_on_segments(selfedge, otheredge):
                        self.intersections[half_i].append(euclidean.intersect(selfedge, otheredge))
    """
    def iter_all_edges(self):
        for i in range(len(self.upperhalf) - 1):
            yield self.upperhalf[i], self.upperhalf[i + 1]
        for i in range(len(self.lowerhalf) - 1):
            yield self.lowerhalf[i], self.lowerhalf[i + 1]

    def plot_congruent(self, plt, color):
        tools.plot_congruent(plt, self.upperhalf, connect_ends=False, color=color)
        tools.plot_congruent(plt, self.lowerhalf, connect_ends=False, color=color)
