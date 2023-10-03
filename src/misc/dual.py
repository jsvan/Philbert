import numpy as np
from misc import tools, euclidean
from hilbert import polygon
"""
good to use subplots side by side to show this

plt.subplot(1, 2, 1)
plt.subplot(1, 2, 2)
The parameters for subplot are: number of rows, number of columns, and which subplot you're currently on. 
Before showing, use 
plt.tight_layout()
"""


class Polar:
    """
    For the conversion, all names refer to destination type

    (a, b) --> ax + by = 1
    """
    origin = np.array([0, 0])

    @staticmethod
    def to_line(point):
        """
        I think returning two points detailing the line is the status quo I need to uphold.

        Seems to work

        set x y points as 1, 0 and 0, 1
        ax = 1
            x = 1/a, y = 0
        by = 1
            x = 0,   y = 1/a,
        """
        recip = 1 / point
        a, b = np.array([recip[0], 0]), \
               np.array([0, recip[1]])
        return a, b

    @staticmethod
    def to_point(line, q=None):
        """
        Input is either a line (tuple of two points), or two points
        """
        if q is not None:
            p = line
        else:
            p, q = line

        diff = q - p
        cx, cy = np.flip(p) * diff
        c = cy - cx
        a = diff[1] / c
        b = (p[0] - q[0]) / c
        return np.array([a, b])

    @staticmethod
    def v2v(vertices):
        """
        vertice to vertice transformation
        Doesnt necessarily garrantee convexity
        """
        return [Polar.to_point(vertices[i], vertices[ii]) for i, ii in tools.i_ii(len(vertices))]

    @staticmethod
    def v2v_convex(vertices):
        """
        more likely to be convex..
        """

        lines = [Polar.to_line(v) for v in vertices]
        newverts = [euclidean.intersect(l, lines[(i+1) % len(lines)]) for i, l in enumerate(lines)]
        assert len(newverts) == len(vertices)
        return newverts

    @staticmethod
    def plot_congruent(plt, poly, color=None, zorder=None, axline=False, linewidth=1.5):
        if not type(poly) == polygon.Polygon:
            poly = polygon.Polygon(poly)
        v = poly.vertices
        tangents = [x['i'] for x in poly.point_tangents(Polar.origin)]
        # Divide the polygon into two halves
        smaller, larger = 0, len(v)
        if len(tangents) == 1:
            smaller = tangents[0]
        if len(tangents) == 2:
            smaller, larger = min(tangents), max(tangents)
        yooks = v[smaller:larger]
        zooks = v[larger:] + v[:smaller]
        tools.plot_congruent(plt, yooks, color, connect_ends=len(zooks)<1, zorder=zorder, axline=axline, linewidth=linewidth)
        tools.plot_congruent(plt, zooks, color, connect_ends=len(yooks)<1, zorder=zorder, axline=axline, linewidth=linewidth)

    @staticmethod
    def dist(lineP, lineQ, supportinglineA, supportinglineB):
        """
        uses cross ratio to determine distance
        All lines are arrays of two points
        All lines are 'infinite'
        """
        crosser = [np.array([0,0]), np.array([0,1])]
        inters = [euclidean.intersect(crosser, l) for l in [lineP, lineQ, supportinglineA, supportinglineB]]















class Dual2:
    """
    For the conversion, all names refer to destination type

    (a, b) --> ax - b = y
    """
    @staticmethod
    def to_line(point):
        """
        I think returning two points detailing the line is the status quo I need to uphold.
        """
        a = np.array([-1, np.sum(-1*point)])
        b = np.array([1, point[0] - point[1]])
        return a, b

    @staticmethod
    def to_point(line, q=None):
        """
        Input is either a line (tuple of two points), or two points
        ax - b = y --> (a, b)
        -y + ax = b
        """
        if q is not None:
            p = line
        else:
            p, q = line

        diff = q - p
        a = diff[1] / diff[0]
        b = a * p[0] - p[1]
        return np.array([a, b])

    @staticmethod
    def v2v(vertices):
        """
        vertice to vertice transformation
        """
        return [Dual2.to_point(vertices[i], vertices[ii]) for i, ii in tools.i_ii(len(vertices))]


