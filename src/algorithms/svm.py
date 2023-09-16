"""
Dave's algorithm for finding best dividing line.

"""
from misc import tools


class SVM:

    def __init__(self, o):
        self.omega = o

    def dave_region(self, p1, p2, q, discrete=True):
        """
        Do binary search over vertices

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
        pass
        """
        TODO
        Instead of searching over vertices (i think it's hard to implement), where the vertices are associated with 
        hilbert balls of various radiuses around three points, what if we do a search over spokes, where each spoke is 
        associated with a hilbert ball of a range of radiuses for which its a tangent"""