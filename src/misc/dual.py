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

@dataclass
class WedgeRegion:
    sector: SectorAssoc
    polygon: polygon.Polygon




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
        self.sektoren = [self.sektorieren(i) for i in range(len(self.o))]
        # Wedge lines should NOT overlap. Each wedgeline_a represents the start of the boundary declared by its wedge.
        # Wedgeline_a[0] is the x coordinate of the first wedgeline.
        # Sorting this is as good as sorting the wedges.
        self.sektoren.sort(key=lambda sektor: sektor.wedge.aXp[0])


    def wedges_to_mosaic(self, boundary_radius):
        X, Y = 0, 1
        regions = []
        bottomleftbound = euclidean.Vertex((-boundary_radius, -boundary_radius), 0)
        bottomrightbound = euclidean.Vertex((boundary_radius, -boundary_radius), 1)
        topleftbound = euclidean.Vertex((-boundary_radius, boundary_radius), 2)
        toprightbound = euclidean.Vertex((boundary_radius, boundary_radius), 3)
        bottombound = euclidean.Edge(bottomleftbound, bottomrightbound)
        topbound = euclidean.Edge(topleftbound, toprightbound)
        leftbound = euclidean.Edge(bottomleftbound, topleftbound)
        rightbound = euclidean.Edge(bottomrightbound, toprightbound)
        boundaries = [bottombound, leftbound, topbound, rightbound]
        def intersect_boundaries(l):
            intersections = []
            for b in boundaries:
                intersect = euclidean.intersect(l, b)
                #print(f"intersect { intersect}, which makes {abs(intersect[X])} and {abs(intersect[Y])}, and finally {tools.lte(abs(intersect[X]), boundary_radius)} and {tools.lte(abs(intersect[Y]), boundary_radius)}")
                if tools.lte(abs(intersect[X]), boundary_radius) and tools.lte(abs(intersect[Y]), boundary_radius):
                    intersections.append((intersect, b))
                    if len(intersections) == 2:
                        #print("success!")
                        return intersections
            print(f"Intersecting {l} with boundaries failed. I have {intersections} of length {len(intersections)}.")
            return intersections
        for sector in self.sektoren:
            # Every wedge runs through omega and kinda has two sections. The one section is a triangle that begins
            # at the omega's vertex. The second section is like the butt end of a triangle, that starts with an
            # opposite edge. I am "Framing" the space with the bound lines. Therefore, each wedge defines two
            # regions, dictated by the wedge lines and the boundary lines.
            # The point of the wedge is always wedge.line_ab.b .
            # Can take intersections of wedgelines with boundaries, and both should be oriented in the same direction.
            # Attempt to intersect them with all four walls. Will have four points.

            wla = sector.wedge.line_a
            wlb = sector.wedge.line_b
            wpoint = sector.wedge.line_a.b
            # Each one is a tuple of (intersection, edge) (intersection, edge)
            a_inters = intersect_boundaries(wla)
            b_inters = intersect_boundaries(wlb)

            # The points should pass a series of Orient tests to know do they belong to the pointy or flat end of
            # the wedge.
            # a group should have centerpoint, l_a.p, l_b.p be clockwise.
            # Orient test to group point tuples:
            # Uncompleted polygons. We still will need to correct their associations and then attach the wedge pts.
            poly_a = [a_inters[0], b_inters[0]]
            poly_b = [a_inters[1], b_inters[1]]
            # Counter_CW is WRONG
            if euclidean.orient(wpoint, a_inters[0][0], b_inters[0][0]) == euclidean.COUNTER_CW:
                # switch the latter points
                poly_a[1], poly_b[1] = poly_b[1], poly_a[1]

            # Now, if the points are on different boundaries, that means that there should be a corner in the middle.
            # Attach corners:
            for poly in [poly_a, poly_b]:
                edge_a = poly[0][1]
                edge_b = poly[1][1]
                poly[0] = poly[0][0]
                poly[1] = poly[1][0]
                if edge_a.point_a.i != edge_b.point_a.i:
                    poly.append(edge_a.shared(edge_b).v)

            # Now, we need to attach the wedge points.
            # COUNTER_CW is CORRECT
            if euclidean.orient(wpoint, poly_a[0], wlb.a) == euclidean.COUNTER_CW:
                poly_a.append(wpoint)
                poly_b.append(wla.a)
                poly_b.append(wlb.a)

            else:
                poly_b.append(wpoint)
                poly_a.append(wla.a)
                poly_a.append(wlb.a)
            regions.append(WedgeRegion(sector, polygon.Polygon(poly_a)))
            regions.append(WedgeRegion(sector, polygon.Polygon(poly_b)))
        return regions

    def annotate_wedges(self, plt):
        """
        name all wedges by average wedgeline, on both the positive and negative side
        point wedgelines
        """
        def sector_name(sektor):
            # The b'th point is the tangent vertex always
            # and it should be lowercase.
            sektor.wedge.line_a.b.i



        tools.annotate(plt, [s.wedge.line_a.a.i for s in self.sektoren])

    def sektorieren(self, i):
        o_edge = self.o.edge(i)
        oXp = euclidean.intersect(self.p, o_edge)
        o_i = self.o.i(i)
        o_tangent = self.o.other_tangent(oXp, o_i)
        oXp = euclidean.intersect(self.p, o_edge)
        wedgeline_a = euclidean.BasicLine(o_edge.point_a, o_tangent)
        wedgeline_b = euclidean.BasicLine(o_edge.point_b, o_tangent)

        wedge = Wedge(wedgeline_a, wedgeline_b,
                      euclidean.intersect(wedgeline_a, self.p),
                      euclidean.intersect(wedgeline_b, self.p))
        return SectorAssoc(o_edge, o_tangent, o_i, oXp, wedge)


    def pseudohyperbola(self, radius, halvedpoly=None):
        """
        The algorithm I want, but it breaks because there are [nan, nan] intersection points which destroys the flow.
        """
        firstwedge = self.sektoren[0].wedge
        pos_pt = geometry.hdist_to_euc(firstwedge.aXp, firstwedge.line_a.a, firstwedge.line_a.b, radius)
        neg_pt = geometry.hdist_to_euc(firstwedge.aXp, firstwedge.line_a.a, firstwedge.line_a.b, -1 * radius)
        pointsabove, pointsbelow = [], []
        for sektor in self.sektoren[0:-1]:
            line = sektor.wedge.line_b
            pos_pt = euclidean.intersect((sektor.oXp, pos_pt), line)
            neg_pt = euclidean.intersect((sektor.oXp, neg_pt), line)
            if euclidean.point_below_line(pos_pt, self.p):
                pointsbelow.append(pos_pt)
            else:
                pointsabove.append(pos_pt)
            if euclidean.point_below_line(neg_pt, self.p):
                pointsbelow.append(neg_pt)
            else:
                pointsabove.append(neg_pt)

        if halvedpoly is None:
            return polygon.HalvedPolygon(pointsabove, pointsbelow, dividingline=self.p)

        halvedpoly.upperhalf = pointsabove
        halvedpoly.lowerhalf = pointsbelow
        halvedpoly.dividingline = self.p
        return halvedpoly






class Polar:
    """
    For the conversion, all names refer to destination type

    (a, b) --> ax + by = 1
    """
    origin = Point([0, 0])

    @staticmethod
    def polygon(vertices, offset=[0,0]):
        vertices = [offset + np.array(v) for v in vertices]
        vertices = Polar.v2v(vertices)
        return polygon.Polygon(vertices)

    @staticmethod
    def omega(vertices, offset=[0,0]):
        vertices = [offset + np.array(v) for v in vertices]
        vertices = Polar.v2v(vertices)
        return omega.Omega(vertices)

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
    def dist(lineP, lineQ, supportinglineA, supportinglineB, absoluteval=True):
        """
        uses cross ratio to determine distance
        All lines are arrays of two points
        All lines are 'infinite'

        Tested and works.
        """
        # Horizontal line
        crosser = BasicLine(Point([0, 1]),
                            Point([1, 1]))

        p, q, A, B = [euclidean.intersect(crosser, l) for l in [lineP, lineQ, supportinglineA, supportinglineB]]
        return geometry.dist(p, q, A, B, absoluteval)


    @staticmethod
    def derivative(lineq, lineA, lineB):
        """
        If you put in lineA and lineB such that they all intersect at a point on p-naught, the derivative should be
        the speed at which the point approaches pnaught
        """
        # Horizontal line
        crosser = BasicLine(Point([0, 1]),
                            Point([1, 1]))
        q, A, B = [euclidean.intersect(crosser, l) for l in [lineq, lineA, lineB]]
        return geometry.derivative(q, A, B)



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










