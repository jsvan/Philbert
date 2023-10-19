import numpy as np
from misc import tools, euclidean
from hilbert import polygon, geometry, omega
from dataclasses import dataclass
from misc.euclidean import Vertex, BasicLine, Point

"""
good to use subplots side by side to show this

plt.subplot(1, 2, 1)
plt.subplot(1, 2, 2)
The parameters for subplot are: number of rows, number of columns, and which subplot you're currently on. 
Before showing, use 
plt.tight_layout()
"""


@dataclass
class Wedge:
    line_a: BasicLine
    line_b: BasicLine
    aXp: Point
    bXp: Point

@dataclass
class SectorAssoc:
    # o is omeganaught
    # w is wedge
    # p is pnaught
    # X is crossing/intersection
    o_edge: euclidean.Edge
    o_tangent: euclidean.Vertex
    o_i: int
    oXp: Point
    wedge: Wedge


class Polar:
    """
    For the conversion, all names refer to destination type

    (a, b) --> ax + by = 1
    """
    origin = Point([0, 0])

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
        recip = Point([1, 1]) / point
        a = Point([recip[0], 0])
        b = Point([0, recip[1]])
        return BasicLine(a, b)

    @staticmethod
    def to_point(line, q=None):
        """
        Input is either a line (tuple of two points), or two points
        """
        if q is not None:
            p = line
        else:
            p, q = line
        if type(p) is not Point:
            p = Point(p)
        if type(q) is not Point:
            q = Point(q)
        diff = q - p
        cx, cy = np.flip(p.v) * diff.v
        c = cy - cx
        a = diff[1] / c
        b = (p[0] - q[0]) / c
        return Point([a, b])

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
        lines = [Polar.to_line(v.v) for v in vertices]
        newverts = [Vertex(euclidean.intersect(l, lines[(i+1) % len(lines)]), i) for i, l in enumerate(lines)]
        assert len(newverts) == len(vertices)
        return newverts

    @staticmethod
    def plot_congruent(plt, poly, color=None, zorder=None, axline=False, linewidth=1.5):
        if not type(poly) == polygon.Polygon:
            poly = polygon.Polygon(poly)
        v = poly.vertices
        tangents = [x.i for x in poly.point_tangents(Polar.origin)]
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

        Tested and works.
        """
        # Horizontal line
        crosser = [Point([0, 0.5]), Point([1, 0.5])]

        p, q, A, B = [euclidean.intersect(crosser, l)[0] for l in [lineP, lineQ, supportinglineA, supportinglineB]]
        return geometry.dist(p, q, A, B)





class Sectors:

    def __init__(self, pnaught, omeganaught):

        """
        variables:
          omega
          p-naught
          intersections of omega onto p-naught
          wedges: p-naught to edge and tangent

        All items will be i indexed
        functions:
          hilbert ball/wave/pseudospline from given radius
          how:
            finds where intersection point intersects each sector, using the prev intersection point
        """
        self.p = pnaught
        self.o = omeganaught
        # List of where omega's edges rays intersect our pnaught
        self.sektoren = [self.sektorieren(i) for i in range(len(self.o.polygon))]
        # Wedge lines should NOT overlap. Each wedgeline_a represents the start of the boundary declared by its wedge.
        # Wedgeline_a[0] is the x coordinate of the first wedgeline.
        # Sorting this is as good as sorting the wedges.
        self.sektoren.sort(key=lambda sektor: sektor.wedge.aXp[0])

    def sektorieren(self, i):
        o_edge = self.o.polygon.edge(i)
        oXp = euclidean.intersect(self.p, o_edge)
        o_i = self.o.polygon.i(i)
        o_tangent = self.o.polygon.other_tangent(oXp, o_i)
        oXp = euclidean.intersect(self.p, o_edge)
        wedgeline_a = euclidean.BasicLine(o_edge.point_a, o_tangent)
        wedgeline_b = euclidean.BasicLine(o_edge.point_b, o_tangent)

        wedge = Wedge(wedgeline_a, wedgeline_b,
                      euclidean.intersect(wedgeline_a, self.p),
                      euclidean.intersect(wedgeline_b, self.p))
        return SectorAssoc(o_edge, o_tangent, o_i, oXp, wedge)


    def pseudohyperbola(self, radius, plt=None):
        """
        The algorithm I want, but it breaks because there are [nan, nan] intersection points which destroys the flow.
        """
        firstwedge = self.sektoren[0].wedge
        pos_pt = geometry.hdist_to_euc(firstwedge.aXp, firstwedge.line_a.a, firstwedge.line_a.b, radius)
        neg_pt = geometry.hdist_to_euc(firstwedge.aXp, firstwedge.line_a.a, firstwedge.line_a.b, -1 * radius)
        pos_list, neg_list = [pos_pt], [neg_pt]
        for sektor in self.sektoren[0:]:
            line = sektor.wedge.line_b
            pos_pt = euclidean.intersect(BasicLine(sektor.oXp, pos_pt), line)
            neg_pt = euclidean.intersect(BasicLine(sektor.oXp, neg_pt), line)
            pos_list.append(pos_pt)
            neg_list.append(neg_pt)

        return pos_list, neg_list








        """
        Naive way, Also doesnt work.         
        
        pos_list, neg_list = [], []
        for sektor in self.sektoren:
            wedge = sektor.wedge
            pos_list.append(geometry.hdist_to_euc(wedge.aXp, wedge.line_a.a, wedge.line_a.b, radius))
            neg_list.append(geometry.hdist_to_euc(wedge.aXp, wedge.line_a.a, wedge.line_a.b, -1 * radius))
            pos_list.append(geometry.hdist_to_euc(wedge.bXp, wedge.line_b.a, wedge.line_b.b, radius))
            neg_list.append(geometry.hdist_to_euc(wedge.bXp, wedge.line_b.a, wedge.line_b.b, -1 * radius))
        return pos_list, neg_list
        """

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










