import numpy as np
from misc import euclidean, tools
from misc.graham_scan import graham_scan
import operator
from pprint import pprint


class Polygon:

    def __init__(self, vertices):
        self.vertices = graham_scan(vertices)
        self.vertices_expanded = [self.vertices[-1]] + self.vertices + [self.vertices[0]]

    def v(self, i):
        return self.vertices[i % len(self.vertices)]

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
