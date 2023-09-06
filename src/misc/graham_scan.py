from misc.euclidean import orient as direction
from functools import cmp_to_key

# Taken from https://algorithmtutor.com/Computational-Geometry/Convex-Hull-Algorithms-Graham-Scan/


def find_min_y(points):
    miny = 999999
    mini = 0
    for i, point in enumerate(points):
        if point[1] < miny:
            miny = point[1]
            mini = i
        if point[1] == miny:
            if point[0] < points[mini][0]:
                mini = i
    return points[mini], mini


def distance_sq(a,b):
    d = b-a
    return sum(d * d)


def polar_comparator(p1, p2, p0):
    d = direction(p0, p1, p2)
    if d < 0:
        return -1
    if d > 0:
        return 1
    if d == 0:
        if distance_sq(p1, p0) < distance_sq(p2, p0):
            return -1
        else:
            return 1


def graham_scan(points):
    # let p0 be the point with minimum y-coordinate,
    # or the leftmost such point in case of a tie
    p0, index = find_min_y(points)

    # swap p[0] with p[index]
    points[0], points[index] = points[index], points[0]

    # sort the points (except p0) according to the polar angle
    # made by the line segment with x-axis in anti-clockwise direction
    sorted_polar = sorted(points[1:], key=cmp_to_key(lambda p1, p2: polar_comparator(p1, p2, p0)))

    # if more than two points are collinear with p0, keep the farthest
    to_remove = []
    for i in range(len(sorted_polar) - 1):
        d = direction(sorted_polar[i], sorted_polar[i + 1], p0)
        if d == 0:
            to_remove.append(i)
    sorted_polar = [i for j, i in enumerate(sorted_polar) if j not in to_remove]

    m = len(sorted_polar)
    if m < 2:
        return

    stack = [points[0], sorted_polar[0], sorted_polar[1]]
    stack_size = 3

    for i in range(2, m):
        while True:
            d = direction(stack[stack_size - 2], stack[stack_size - 1], sorted_polar[i])
            if d < 0:  # if it makes left turn
                break
            else:  # if it makes non left turn
                stack.pop()
                stack_size -= 1
        stack.append(sorted_polar[i])
        stack_size += 1
    return stack