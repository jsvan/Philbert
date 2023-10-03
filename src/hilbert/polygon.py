import numpy as np
from misc import euclidean, tools
from hilbert import line
from misc.graham_scan import graham_scan
from matplotlib import pyplot as plt


class Polygon:

    def __init__(self, vertices, offset = [0, 0]):
        self.offset = np.array(offset)
        self.vertices = graham_scan([self.offset + np.array(v) for v in vertices])
        self.vertices_expanded = [self.vertices[-1]] + self.vertices + [self.vertices[0]]

    def v(self, i):
        return self.vertices[i % len(self.vertices)]

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
        print("Tangents", tangents)
        return tangents

    def point_tangent_linear(self, point, maxwards=True):
        """
        linear time
        tangent line from point
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
                return {'i':vi, 'v':v}
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
                return {'i': vi, 'v': v}
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
                return {'i':vi, 'v':v}
        raise Exception("Boundary failed, point potentially inside.")


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
            print(f"i-1 is {euclidean.dname(minus_orientation)}, f+1 is {euclidean.dname(plus_orientation)}")
            if minus_orientation == wrong_orientation:
                binarysearch.feedback(higher=maxwards)
                print(
                    f"i-1 is {euclidean.dname(minus_orientation)}, so going {'higher' if maxwards else 'lower'}.")

                continue
            if plus_orientation == wrong_orientation:
                # Too far
                binarysearch.feedback(higher=not maxwards)
                print(
                    f"i+1 is {euclidean.dname(plus_orientation)}, so going {'higher' if not maxwards else 'lower'}.")

                continue
            return v
        raise Exception("Finding boundary with binary search failed")


    def halves(self, l):
        """
        l is a line.Line, or two 2d numpy arrays like [np.a([x, y]), np.a([x, y])]
        returns two arrays of vertices of the polygon, for each side of line they're on
        """
        if type(l) is line.Line:
            l = l.l

        points_above_l_before, points_below_l_before, points_above_l_after, points_below_l_after  = [], [], [], []
        loadingpoints = [points_above_l_before, points_below_l_before]

        for v in self.vertices:
            if euclidean.point_below_line(v, l):
                loadingpoints[1].append(v)
                if len(loadingpoints[0]) > 0:
                    loadingpoints[0] = points_above_l_after
            else:
                loadingpoints[0].append(v)
                if len(loadingpoints[1]) > 0:
                    loadingpoints[1] = points_below_l_after

        points_above_l = points_above_l_after + points_above_l_before
        points_below_l = points_below_l_after + points_below_l_before
        return points_above_l, points_below_l

    def halves(self, l, index=False):
        """
        l is a line.Line, or two 2d numpy arrays like [np.a([x, y]), np.a([x, y])]
        returns two arrays of vertices of the polygon, for each side of line they're on
        """
        if type(l) is line.Line:
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
            #print(i, 'lv:', euclidean.orient(p, q, lv), ', uv:', euclidean.orient(p, q, uv))

            if euclidean.orient(p, q, lv) != euclidean.CLOCKWISE and euclidean.orient(p, q, uv) != euclidean.COUNTER_CW:
                return (lv, uv)
        #tools.plot_congruent(plt, self.vertices, color='gray')
        #tools.scatter(plt, [p, q])
        #tools.annotate(plt, 'pq', [p, q])
        #tools.annotate(plt, range(numvert), self.vertices)
        #plt.show()
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
            return (leftv, rightv)

        raise Exception("Finding boundary with binary search failed")

    """
    for i in range(len(self.vertices)):
        vm, vp = self.v(i), self.v(i+1)
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
        :return: two lines, each of two points.
        """
        return self.find_boundary(p, q), self.find_boundary(q, p)

    def center(self):
        x = sum([v[0] for v in self.vertices]) / len(self.vertices)
        y = sum([v[1] for v in self.vertices]) / len(self.vertices)
        return np.array([x, y])