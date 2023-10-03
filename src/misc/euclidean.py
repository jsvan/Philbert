import numpy as np
from misc import tools
import math
TWO_PI = 2*math.pi
COUNTER_CW = 1
CLOCKWISE = -1
X, Y = 0, 1
EPS = 0.00001

def dname(w):
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
    return np.dot(q-p, v-p)


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

    xy12 = det(x1, y1, x2, y2)
    xy34 = det(x3, y3, x4, y4)
    x1m2 = x1 - x2
    x3m4 = x3 - x4
    y1m2 = y1 - y2
    y3m4 = y3 - y4
    denominator = det(x1m2, y1m2, x3m4, y3m4)

    x = det(xy12, x1m2, xy34, x3m4) / denominator
    y = det(xy12, y1m2, xy34, y3m4) / denominator
    return np.array([x, y])


def line_segment_normal(p, q):
    dx,dy = p - q
    return np.array([-dy, dx]), np.array([dy,-dx])


def norm_towards_q(p1, p2, q):
    ns = line_segment_normal(p1, p2)
    return ns[np.dot(q - p1, ns[0]) < 0]


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

def eq(p1, p2):
    if p1 is None or p2 is None:
        return False
    return np.isclose(p1, p2).all()


def get_point_on_line(p, q, dist):
    """
    returns point dist away from p, on line through q
    dist==0 is point p
    dist==1 is point q
    """
    slope = q - p
    return p + dist*slope



def uniform_sample_from_line_segments(segs, startfloor=0):
    distances = [np.linalg.norm(segs[ii] - segs[i]) for i, ii in tools.i_ii(len(segs), connect_ends=False)]
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
    Probably a stupid method, but it will you the ith point around a 'circle' with a radius.
    I use this later for a line that spins around.
    """

    slice =  TWO_PI * (i / divisions)
    return np.array([math.cos(slice), math.sin(slice)])