from operator import truediv, floordiv
import numpy as np
from misc import euclidean
from matplotlib.colors import TABLEAU_COLORS as COLORS
COLORS = list(COLORS.values())
X, Y = 0, 1

def namer(i):
    i = int(i)
    letters = 'abcdefghijklmnopqrstuvwxy'
    round = i // len(letters)
    return 'z'*round + letters[i%len(letters)]


class BinarySearcher:
    """
    [mini, maxi)
    inclusive, exclusive
    Rotate changes the first index that's expelled. Instead of producing the median index, it will dial it forwards or
    if negative backwards.
    """
    def __init__(self, mini, maxi, discrete=False, rotate=0):
        if mini > maxi:
            mini, maxi = maxi, mini
        self.mini = mini
        self.maxi = maxi
        self.discrete = discrete
        self.divisor = floordiv if self.discrete else truediv
        self.mid = self.calc_mid() if rotate == 0 else ((self.calc_mid() + rotate - self.mini) % (self.maxi - self.mini ) ) + self.mini


    def calc_mid(self):
        return self.divisor((self.maxi - self.mini), 2) + self.mini

    def next(self):
        return self.mid

    """
    higher is a boolean
    """
    def feedback(self, higher):
        #if not self.has_next():
        #    raise StopIteration

        if higher:
            self.mini = self.mid
        else:
            self.maxi = self.mid
        self.mid = self.calc_mid()

    def has_next(self):
        return (self.discrete and (self.mini + 1) != self.maxi) or \
        (not self.discrete and not np.isclose(self.mini, self.maxi))

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


def i_ii(n, connect_ends=True):
    # generator for indeces i and i+1, wrapping around a polygon for example.
    # Connect-ends gives you that final wrap around, otherwise, connects piecewise the indeces.
    mn = n if connect_ends else n-1
    return ((i, (i + 1) % n) for i in range(mn))


def sortnparray(a, b):
    i = 0
    while i < len(a):
        if a[i] < b[i]:
            return a, b
        elif a[i] > b[i]:
            return b, a
        i += 1
    return a, b


def plot_congruent(plt, vertices, color=None, connect_ends=True, axline=False, zorder=None, linewidth=1.5):
    for i, ii in i_ii(len(vertices), connect_ends):
        c = color[i] if color is not None and type(color) == list else color
        if axline:
            plt.axline(vertices[i].v, vertices[ii].v, color=c, zorder=zorder, linewidth=linewidth)
        else:
            plot_line(plt, vertices[i], vertices[ii], color=c, zorder=zorder, linewidth=linewidth)


def plot_line(plt, a, b, color=None, zorder=None, linewidth=1.5, axline=False, linestyle=None):
    """
    just because matplotlib wants all the x's together and all the y's together, it's pretty annoying.
    Give this your two points and it'll handle it.
    :param plt:
    :param a:
    :param b:
    :return:
    """
    #try:
    if plt is not None and a is not None and b is not None:
        if euclidean.eq(a, b):
            print(a, b, color, zorder, linewidth, axline)
        if axline:
            plt.axline(a.v, b.v, color=color, zorder=zorder, linewidth=linewidth, linestyle=linestyle)
        else:
            plt.plot([a[X], b[X]], [a[Y], b[Y]], color=color, zorder=zorder, linewidth=linewidth, linestyle=linestyle)
    #except Exception as e:
    #    print(e, e.__doc__)
    #    print(a, b, type(a), type(b), 'failed')

def scatter(plt, points, colors=None, zorder=100, size=None, ):
    c = None
    for i, p in enumerate(points):
        if p is not None:
            if colors is not None:
                if len(colors) == 1:
                    i = 0
                c = colors[i]
            plt.scatter(*p, color=c, zorder=zorder, s=size)



"""
annotations happen slightly up and to the right of the point. 
"""
def annotate(plt, names, points, seen={}):
    assert len(names) == len(points)
    for i in range(len(names)):
        old = ""
        if points[i] is not None:
            key = tuple(points[i].v)
            new = str(names[i])
            if key in seen:
                old = seen[key].get_text()
                seen[key].remove()
                if old == new:
                    old = ""
            seen[key] = plt.annotate(old + new, points[i].v)

    return seen


"""
returns evenly spaced numbers of a specific width over an interval.
Unline numpy's, it snaps to the grid which includes the origin (0, 0).
linspace(5.5, 6.5, 0.2) would return np.array([5.6, 5.8, 6.0, 6.2, 6.4])
Does not include upper boundary.
"""
def linspace(start, stop, step):

    if start > stop:
        start, stop = stop, start

    startgap = (step - (start % step)) % 1  # Snap to grid
    return np.arange(start + startgap, stop, step)


def rand_points(how_many, offset=[0,0]):
    offset=np.array(offset)
    return [euclidean.Point(x[:2] + offset) for x in np.random.dirichlet([1,1,1], how_many)]

def lte(a, b):
    return a <= b or np.isclose(a, b)

def gte(a, b):
    return a >= b or np.isclose(a, b)

def save_movie_frame(plt, i):
    plt.savefig(f"./movieframes/{namer(i)}.png")
    plt.close()


def finish_movie(outname):
    cmd = f"ffmpeg -framerate 20 -pattern_type glob -i \"./movieframes/*.png\"   -c:v libx264 -pix_fmt yuv420p ./movies/{outname}.mp4"
    import os
    os.system('ls ./')
    os.system(cmd)
    os.system(f"xdg-open ./movies/{outname}.mp4")
    os.system(f"rm ./moviesframes/*.png")

