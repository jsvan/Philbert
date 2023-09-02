from math import sqrt, log
from math import e as the_E
import numpy as np
from numpy.linalg import norm
from misc import graham_scan, euclidean


def hdist_to_euc(p, A, B, hdist):
   """ This is solving the hilbert distance for t
       A, p, B, where A/B are on the boundaries
       works for negative direction too.
   p = (a, b)
   A = (a + kc, b + kd)
   B = (a + sc, b + sd)
   q = (a + tc, b + td)
   """
   # a, b = p    # not actually needed
   slope = p - B if np.isclose(p, A).all() else p - A  # slope or diff or whatever
   flip = slope[0] == 0
   # This is how i deal with vertical divide by 0 issues. Pretend it's horizontal!
   if flip:  # vertical
      p = np.flip(p)
      A = np.flip(A)
      B = np.flip(B)
      slope = np.flip(slope)

   c, d = slope
   k = (A[0] - p[0]) / c
   s = (B[0] - p[0]) / c

   # if B is p, then the equation will return inf. Better if it returns 0.
   if s == 0:
      return B
   K = (pow(the_E, 2 * hdist) * (pow(k * c, 2) + pow(k * d, 2))) / (pow(s * c, 2) + pow(s * d, 2))
   t = (-k - sqrt((K * pow(k - s, 2))) + K * s) / (-1 + K)
   point = p + (t * slope)
   if flip:
      point = np.flip(point)
   return point


def nielson_dist(p, q):
   """
   Nielson's code for finding hilbert distance. Runs on order of number of dimensions.
   :param p: np.array(coordinates)
   :param q: np.array(coordinates)
   :return: float
   """
   if np.allclose(p, q): return 0  # np.allclose returns true if all dim are within EPSILON
   idx = np.logical_not(np.isclose(p, q))  # np.isclose returns an array of bool, for each dimension if within EPSILON
   if (idx.sum() == 1): return 0  # returns 0 if only ONE dimension is not close (??)
   lamb = p[idx] / (p[idx] - q[idx])  # does dimension-wise arithmetic on those dimensions which are not within epsilon.
   t0 = lamb[lamb <= 0].max()  # t0 = resulting max element less than 0
   t1 = lamb[lamb >= 1].min()  # t1 = resulting min element bigger than 1
   if np.isclose(t0, 0) or np.isclose(t1, 1): return np.inf  # From here on idk
   return np.abs(np.log(1 - 1 / t0) - np.log(1 - 1 / t1)) # mathematically equivalent to dist algorithm tested and approved


def hilbert_ball_around_point(p, omega, r):
   spokes = omega.spokes(p)
   # I think it's best to find the q point on every spoke, instead of connecting to vanishing points.
   # Vanishing points require binary search probably for every spoke intersection (on opposite side).
   # That may not actually be slower, but this is simpler in code for sure.
   points = []
   for v, _ in spokes:
      boundaries = omega.find_both_boundary_intersections_of_line(p, v)
      intersections = ordered_double_intersect(p, v, boundaries)
      points.append(hdist_to_euc(p, *intersections, r))
      points.append(hdist_to_euc(p, *intersections, -r))
   return graham_scan.graham_scan(points)


def dist(p, q, A, B):
   """
   A--p--q--B
      |-pB--|
         |qB|
   |-qA--|
   |pA|
   :param p: np.array(coordinates)
   :param q: np.array(coordinates)
   :param A: np.array(coordinates) boundary intersection coord close to p
   :param B: np.array(coordinates) boundary intersection coord close to q
   :return: float, hilbert distance between p, q
   """
   pA = norm(p - A)
   pB = norm(p - B)
   qA = norm(q - A)
   qB = norm(q - B)
   crossratio = (pB / qB) * (qA / pA)
   return abs(log(crossratio))


def vectorform_dist(p_, q_, u_, v_):
   """
   Four vectors, u_ and v_ are vectors of boundaries (not the intersection points!)
      : u_ p_ q_ v_

   :param p_:
   :param q_:
   :param u_:
   :param v_:
   :return:
   """
   return np.dot(u_, q_) * np.dot(v_, p_) / (np.dot(u_, p_) * np.dot(v_, q_))


def ordered_double_intersect(p, q, boundaries):
   """
   :param p: np.array(coordinates)
   :param q: np.array(coordinates)
   :param boundaries: list of two boundaries. Each boundary is a list of two points.
   :return: points A and B, of order A--p--q--B
   """
   A, B = euclidean.intersect([p, q], boundaries[0]), \
      euclidean.intersect([p, q], boundaries[1])
   # This is if A and B are swapped, which is handy because my discovery algorithm doesn't handle orientation for me
   # even though ironically it uses something called the orient() algorithm a dozen times ...
   """
   pB = norm(p - B)
   qB = norm(q - B)
   if qB > pB:
      A, B = B, A
   """
   return A, B



