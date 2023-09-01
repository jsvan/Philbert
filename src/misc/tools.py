from operator import truediv, floordiv
import numpy as np
X, Y = 0, 1

def namer(i):
    letters = 'abcdefghijklmnopqrstuvwxy'
    round = i // len(letters)
    return 'z'*round + letters[i%len(letters)]


class BinarySearcher:
    """
    [mini, maxi)
    inclusive, exclusive
    """
    def __init__(self, mini, maxi, discrete=False):
        if mini > maxi:
            mini, maxi = maxi, mini
        self.mini = mini
        self.maxi = maxi
        self.discrete = discrete
        self.divisor = floordiv if self.discrete else truediv
        self.mid = self.calc_mid()

    def calc_mid(self):
        return self.divisor((self.maxi - self.mini), 2) + self.mini

    def next(self):
        return self.mid

    """
    higher is a boolean
    """
    def feedback(self, higher):
        if (self.discrete and self.mini == self.maxi) or \
                (not self.discrete and np.isclose(self.mini, self.maxi)):
            raise StopIteration

        if higher:
            self.mini = self.mid
        else:
            self.maxi = self.mid
        self.mid = self.calc_mid()


def bs_on_points(vertices, transform, comparator):
    """
    :param vertices iterable to find the most comparator of, of the transformed items
    :param transform:
    :param comparator:
    :return:
    """
    num_v = len(vertices)
    if num_v == 0:
        return None, None, None

    partition = num_v // 3
    a = transform(vertices[partition])
    b = transform(vertices[2 * partition])
    start, end = 0, num_v
    # pretend comparator is '>', so True if a > b
    if comparator(a, b):
        start = partition
        prev_max = b
        prev_max_i = 2 * partition
    else:
        end = 2 * partition
        prev_max = a
        prev_max_i = partition

    bs = BinarySearcher(start, end, discrete=True)
    while True:
        i = bs.next()
        vertex = vertices[i]
        a = transform(vertex)
        try:
            if comparator(a, prev_max):
                bs.feedback(higher=i > prev_max_i)
                prev_max = a
                prev_max_i = i
            else:
                bs.feedback(higher=i < prev_max_i)

        except StopIteration:
            break
    # At this point I should have the
    return vertex, a, i


def i_ii(n):
    # generator
    return ((i, (i + 1) % n) for i in range(n))


def sortnparray(a, b):
    i = 0
    while i < len(a):
        if a[i] < b[i]:
            return a, b
        elif a[i] > b[i]:
            return b, a
        i += 1
    return a, b


def plot_congruent(plt, vertices, color=None):
    for i, ii in i_ii(len(vertices)):
        plot_line(plt, vertices[i], vertices[ii], color=color)

def plot_line(plt, a, b, color=None):
    """
    just because matplotlib wants all the x's together and all the y's together, it's pretty annoying.
    Give this your two points and it'll handle it.
    :param plt:
    :param a:
    :param b:
    :return:
    """
    if plt is not None:
        plt.plot([a[X], b[X]], [a[Y], b[Y]], color=color)