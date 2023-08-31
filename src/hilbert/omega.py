import numpy as np
from src.misc import tools, euclidean

class Omega:

    def __init__(self, vertices=[np.array(x) for x in [[0, 0], [0, 1], [1, 0]] ]):
        self.vertices = vertices
        # This is a hack for binary search, wraps around the final points. Whatever.
        self.vertices_expanded = [vertices[-1]] + vertices + [vertices[0]]

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

    def find_boundary(self, p ,q):
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



