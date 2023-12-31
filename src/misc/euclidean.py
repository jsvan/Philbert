import numpy as np
import math
from dataclasses import dataclass, astuple
TWO_PI = 2*math.pi
COUNTER_CW = 1
CLOCKWISE = -1
X, Y = 0, 1
EPS = 0.00001




class Point:
    def __init__(self, p):
        if type(p) is np.array:
            self.v = p
        elif type(p) is Point or type(p) is Vertex:
            self.v = p.v
        else:
            self.v = np.array(p)
        if self.v.shape != (2,):
            raise Exception(f"Shape error. Need shape (2,), but you gave {self.v.shape}, with input {p}.")

    def __iter__(self):
        return iter(self.v)
    def __getitem__(self, item):
        return self.v[item]
    def __len__(self):
        return len(self.v)
    def __repr__(self):
        return f"Point {self.v}"
    def __hash__(self):
        return tuple(self.v).__hash__()
    def __add__(self, other):
        if type(other) is Point:
            other = other.v
        return Point(self.v + other)
    def __sub__(self, other):
        if type(other) is Point:
            other = other.v
        return Point(self.v - other)
    def __mul__(self, other):
        if type(other) is Point:
            other = other.v
        return Point(self.v * other)
    def __truediv__(self, other):
        if type(other) is Point:
            other = other.v
        return Point(self.v / other)
    def __lt__(self, other):
        return self.v < other.v
    def __gt__(self, other):
        return self.v > other.v
    def __le__(self, other):
        return self.v <= other.v
    def __ge__(self, other):
        return self.v >= other.v
    def __eq__(self, other):
        return eq(self, other)


@dataclass
class BasicLine:
    a: Point
    b: Point
    norm: Point = None
    def __iter__(self):
        return iter((self.a, self.b))
    def __getitem__(self, item):
        return self.b if item else self.a
@dataclass
class Vertex(Point):
    i: int
    v: np.array
    def __init__(self, v, i):
        super().__init__(v)
        self.i = i
@dataclass
class Edge:
    point_a: Vertex
    point_b: Vertex
    def __iter__(self):
        return iter((self.point_a, self.point_b))
    def __getitem__(self, item):
        return self.point_b if item else self.point_a
    def shared(self, other):
        for ov in other:
            for sv in self:
                if ov.i == sv.i:
                    return sv
        return None

def dirname(w):
    if w == COUNTER_CW:
        return "COUNTER_CW"
    if w == CLOCKWISE:
        return "CLOCKWISE"
    return "None"


def det(a, b, c, d):
    return (a*d) - (b*c)


def quick_is_outside(point, vm1, v, vp1):
    """
    Tells if a point is on the "interior" of a "v" formed by 3 vertices.

    A limited version of the polygon's "contains()", but is O(1) rather than O(vertices)

    Returns boolean True if point is on the exterior of the "v", and False if on within the inside rays of the "v"
    """
    bad = orient(point, v, vm1) == COUNTER_CW and orient(point, v, vp1) == CLOCKWISE
    return not bad


def tangent_dot(v, p, q):
    """
    point p lies on line l, which passes through omega. Pivot point.
    q on line l
    v is vertex  not on l
    Takes the dot of qp, vp

    :return: float
    """
    return np.dot((q-p).v, (v-p).v)


def cos(v, p, q):
    """
    v
       `
    q___(_p

    :param v:
    :param p:
    :param q:
    :return:
    """
    a = v - p
    b = q - p
    return np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))


def orient(p, q, r):
    """
    a   b   c
    d   e   f
    g   h   i

    1   px  py
    1   qx  qy
    1   rx  ry

    det = aei + bfg + cdh - ceg - bdi - afh
    det = ei + bf + ch - ce - bi - fh
    det = e*(i-c) + b(f-i) + h(c-f)
    det = q[x] * (r[y] - p[y]) \
        + p[x] * (q[y] - r[y]) \
        + r[x] * (p[y] - q[y])

    :return: 1 == Counter Clockwise, 0 == Straight, -1 == Clockwise
    https://mathworld.wolfram.com/Line-LineIntersection.html
    """
    return np.sign(q[X] * (r[Y] - p[Y])
                 + p[X] * (q[Y] - r[Y])
                 + r[X] * (p[Y] - q[Y]))


def intersect(line_a, line_b):
    """
    :line_a: iterable of two points defining a line
    :line_b: iterable of two points defining a line
    using determinants (following code). Better way?
    https://mathworld.wolfram.com/Line-LineIntersection.html
    """
    (x1, y1), (x2, y2) = line_a
    (x3, y3), (x4, y4) = line_b

    # if xy12 and xy34 is a multiple, it is the same line.
    xy12 = det(x1, y1, x2, y2)
    xy34 = det(x3, y3, x4, y4)
    x1m2 = x1 - x2
    x3m4 = x3 - x4
    y1m2 = y1 - y2
    y3m4 = y3 - y4
    # denominator is 0 if parallel lines.
    denominator = det(x1m2, y1m2, x3m4, y3m4)
    if np.isclose(denominator, 0):
        return None
    x = det(xy12, x1m2, xy34, x3m4) / denominator
    y = det(xy12, y1m2, xy34, y3m4) / denominator
    return Point([x, y])

"""
Return boolean True if lines are the same line, False if not the same line.
"""
def same_line(line_a, line_b):
    (x1, y1), (x2, y2) = line_a
    (x3, y3), (x4, y4) = line_b
    xy12 = det(x1, y1, x2, y2)
    xy34 = det(x3, y3, x4, y4)

    m = (xy34 / xy12) % 1
    return np.isclose(m, 0) or np.isclose(m, 1)

def line_segment_normal(p, q):
    dx, dy = p - q
    return Point([-dy, dx]), Point([dy, -dx])


def norm_towards_q(p1, p2, q):
    ns = line_segment_normal(p1, p2)
    return ns[np.dot(q.v - p1.v, ns[0].v) < 0]


def point_on_line(p, l, segment=False):
    a, b = l
    if segment and \
       not (min(a[Y], b[Y]) - EPS <= p[Y] <= max(a[Y], b[Y]) + EPS and
            min(a[X], b[X]) - EPS <= p[X] <= max(a[X], b[X]) + EPS):
        return False
    return np.isclose((b[X] - a[X]) * (p[Y] - a[Y]),
                      (b[Y] - a[Y]) * (p[X] - a[X]))


def point_below_line(p, l):
    a, b = l
    return (b[X] - a[X])*(p[Y] - a[Y]) < (b[Y] - a[Y])*(p[X] - a[X])


def point_within_region(p, A, B):
    """
    This is my stupid optimization, checks that a point is on a line INSIDE of omega, by making sure no
    dimension of the point is more extreme than of A or B
    """
    for i in range(len(p)):
        if not (A[i] >= p[i] >= B[i] or A[i] <= p[i] <= B[i]):
            return False
    return True

def line_intersects_on_segments(linea, lineb):
    aa = orient(linea[0], *lineb)
    ab = orient(linea[1], *lineb)
    if aa + ab != 0:
        return False
    ba = orient(lineb[0], *linea)
    bb = orient(lineb[1], *linea)
    return ba + bb == 0



def eq(p1, p2):
    if p1 is None or p2 is None or p1.v is None or p2.v is None:
        return False

    return np.isclose(p1.v, p2.v).all()


def get_point_on_line(p, q, dist):
    """
    returns point dist away from p, on line through q
    dist==0 is point p
    dist==1 is point q
    """
    slope = q - p
    return p + slope * dist


def uniform_sample_from_line_segments(segs, startfloor=0):
    from misc.tools import i_ii
    distances = [np.linalg.norm(segs[ii].v - segs[i].v) for i, ii in i_ii(len(segs), connect_ends=False)]
    tot = sum(distances)
    if startfloor >= tot:
        raise Exception(f"Starting at a point {startfloor} greater than the allowed length of {tot}")
    cumulative = np.random.uniform(startfloor, tot)
    cum = cumulative
    for i in range(len(distances)):
        if cum > distances[i]:
            cum -= distances[i]
        else:
            # sample point on line segment at cum ratio
            p, q = segs[i], segs[i+1]
            ratio = cum / distances[i]
            return get_point_on_line(p, q, ratio), cumulative

    raise Exception(f"euclidean Exception, point {cumulative} from length of {tot} not found in segments {segs}.")

def triangle_area(v1, v2, v3):
    """
    0.5 * (x1(y2 - y3) + x2(y3 - y1) + x3(y1 - y2))
    """
    (x1, y1), (x2, y2), (x3, y3) = v1, v2, v3
    return 0.5 * (x1 * (y2 - y3) +
                  x2 * (y3 - y1) +
                  x3 * (y1 - y2))

def circle_point(i, divisions):
    """
    Probably a stupid method, but it will give you the ith point around a 'circle'.
    I use this for a line that spins around.
    """
    sliced = TWO_PI * (i / divisions)
    return np.array([math.cos(sliced), math.sin(sliced)])

"""
Returns the point on the line that is closest/perpendicular to 'point'
"""
def closest_point_on_line(point, the_line):
    if the_line.norm is None:
        the_line.norm = line_segment_normal(*the_line)[0]
    return intersect((point, the_line.norm), the_line)

"""

"""
def line_get_t(point, the_line):
    t = (point - the_line[0]) / the_line[1]
    #print(point, t[1])
    return t[1]

def slope(the_line):
    return (the_line[1][Y] - the_line[0][Y]) / (the_line[1][X] - the_line[0][X])

def intercept(the_line, slope):
    return the_line[0][Y] - slope * the_line[0][X]

def avg_line(linea, lineb):
    bm = slope(lineb)
    bi = intercept(lineb, bm)
    y1 = (linea[0][Y] + (bm * linea[0][X]) + bi) / 2
    y2 = (linea[1][Y] + (bm * linea[1][X]) + bi) / 2
    return BasicLine(Point((linea[0][X], y1)), Point((linea[1][X], y2)))

def avg_point(pa, pb):
    return (pa + pb) / 2