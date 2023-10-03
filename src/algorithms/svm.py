from misc import tools, euclidean
from hilbert import geometry, omega, polygon, line


class SVM:

    def __init__(self, o):
        self.omega = o

    def dave_region(self, p1, p2, q, discrete=True):
        """
        Do binary search over vertices
        Dave's algorithm for finding best dividing line.

        """
        n = len(self.omega.vertices)
        master = tools.BinarySearcher(0, n, discrete=True, rotate=n//4)
        slave = tools.BinarySearcher(0, n, discrete=True, rotate=-n//4)
        """
        Want to do binary search until we find the region where all points are within two 
        
        """
        while go:
            fixednode = master.next()
            while go2:
                slavenode = slave.next()
                if fixednode == slavenode:
                    raise Exception("I don't want to deal with this")




    def julian_spokes(self):
        """
        TODO
        Instead of searching over vertices (i think it's hard to implement), where the vertices are associated with 
        hilbert balls of various radiuses around three points, what if we do a search over spokes, where each spoke is 
        associated with a hilbert ball of a range of radiuses for which its a tangent"""

        """
        Due to lack of monotonicity on the spokes, this method won't work
        """
        p, q, r = tools.rand_points(3)
        o = omega.Omega(vertices=[(1.4, 0.0),
                                    (1.1, 0.8),
                                    (0.6, 1.1),
                                    (0.0, 1.2),
                                    (-.5, 0.8),
                                    (-.8, 0.3),
                                    (-.7, -.2)])

        pr = line.Line(p, r, o)

        # Make large hilbert balls of 1/2 the dist from p to r
        # Find the spokes of that best dividing line
        # Make large hilbert ball of same size around q
        # Find spokes of tangent line pq of those balls
        # This is the upper boundary of spokes to search for our line
        # TODO get spokes from this
        # TODO get radius of implicit hilbert ball from this
        pr_div = pr.get_best_dividing_line()

        # Find a lower boundary
        # Set hilbert balls around p, q to epsilon
        # Find tangent line of those balls p to q
        # This is the lower boundary