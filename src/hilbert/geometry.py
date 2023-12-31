from math import sqrt, log
from math import e as the_E
import numpy as np
from numpy.linalg import norm
from misc import euclidean
from misc.euclidean import Point

"""Mapping isometry onto euclidean?"""

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
   slope = p - B if np.isclose(p.v, A.v).all() else p - A  # slope or diff or whatever
   flip = slope[0] == 0
   # This is how i deal with vertical divide by 0 issues. Pretend it's horizontal!
   if flip:  # vertical
      p.v = np.flip(p.v)
      A.v = np.flip(A.v)
      B.v = np.flip(B.v)
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
      point.v = np.flip(point.v)
   return point


def nielson_dist(p, q):
   """
   Nielson's code for finding hilbert distance. Runs on order of number of dimensions.
   :param p: np.array(coordinates)
   :param q: np.array(coordinates)
   :return: float
   """
   p, q = p.v, q.v
   if np.allclose(p, q): return 0  # np.allclose returns true if all dim are within EPSILON
   idx = np.logical_not(np.isclose(p, q))  # np.isclose returns an array of bool, for each dimension if within EPSILON
   if (idx.sum() == 1): return 0  # returns 0 if only ONE dimension is not close (??)
   lamb = p[idx] / (p[idx] - q[idx])  # does dimension-wise arithmetic on those dimensions which are not within epsilon.
   t0 = lamb[lamb <= 0].max()  # t0 = resulting max element less than 0
   t1 = lamb[lamb >= 1].min()  # t1 = resulting min element bigger than 1
   if np.isclose(t0, 0) or np.isclose(t1, 1): return np.inf  # From here on idk
   return np.abs(np.log(1 - 1 / t0) - np.log(1 - 1 / t1)) # mathematically equivalent to dist algorithm tested and approved


def dist(p, q, A, B, absoluteval=True):
   """
   A--p--q--B
      |-pB--|
         |qB|
   |-qA--|
   |pA|
   :param p: np.array(coordinates)
   :param q: np.array(coordinates)
   :param A: np.array(coordinates) boundary intersection coord close to p\
   :param B: np.array(coordinates) boundary intersection coord close to q
   :return: float, hilbert distance between p, q
   """
   p, q, A, B = p.v, q.v, A.v, B.v
   pA = norm(p - A)
   pB = norm(p - B)
   qA = norm(q - A)
   qB = norm(q - B)
   crossratio = (pB / qB) * (qA / pA)
   d = log(crossratio)
   return abs(d) if absoluteval else d


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


"""
Finds the pointwise derivative of q on the line AB.
 -- A, B constant
 -- derivative as q moves along line.
 -- returns float
 BUG: Possibly... Should AB be normalozed to have a distance of one?
 It looks like no. Why not? 
 
 Minimum derivative is equal to 4 / norm(A - B)
 That means only for the simplex, derivative is almost always much larger than zero.
 
 Also, a simplified way of thinking about the space. The middle half of the cord will have a speed fairly well 
 approximated by the minimum derivative, which is 4 / norm(A-B). But about half the area can be thought of as 
 relatively 'linear', but still fairly speedy. 
 
 Calculating the derivative by hand uses the substitution |Aq| = t and |Bq| = (1 - t), taking derivative with respect
 to t. 
 """
def derivative(q, A, B):
   Aq = np.linalg.norm(A - q)
   Bq = np.linalg.norm(B - q)
   return 1 / Aq + 1 / Bq


"""
Double derivative
"""
def acceleration(q, A, B):
   Aq = np.linalg.norm(A - q)
   Bq = np.linalg.norm(B - q)
   return (1 / (Bq * Bq)) - (1 / (Aq * Aq))

