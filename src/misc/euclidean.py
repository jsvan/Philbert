import numpy as np

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


def is_tangent_orient(p, vm1, v, vp1):
    o = COUNTER_CW
    if orient(p, vm1, v) == COUNTER_CW:
        #positive side of l
        o = CLOCKWISE
    return orient(p, v, vp1) == o


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

def eq(p1, p2):
    if p1 is None or p2 is None:
        return False
    return np.isclose(p1, p2).all()