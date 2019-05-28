import random


class IntervalTree:
    """
    Binary search tree for storing long integer intervals, and for performing queries on them.
    See https://en.wikipedia.org/wiki/Interval_tree, specifically the Augmented kind.
    The present implementation balances the tree by using randomisation.
    """
    root = None # pointer to the root node of the binary search tree
    stack = None

    def __iter__(self):
        self.current = self.root
        self.stack = Stack()
        return self

    def __next__(self):
        while self.current != None:
            self.stack.push(self.current)
            self.current = self.current.left
        if self.stack.isEmpty():
            raise StopIteration
        self.current = self.stack.pop()
        ret = self.current
        self.current = self.current.right
        return ret

    def __len__(self):
        return self.root.N

    def __contains__(self, ival):
        return self.get(ival) != None

    def get(self, ival):
        return self._get(self.root, ival)

    def _get(self, node, ival):
        if node == None: return None
        if ival < node.ival:
            return self._get(node.left, ival)
        elif ival > node.ival:
            return self._get(node.right, ival)
        else:
            return node

    def isect(self, ival, node = None):
        """ Look for intersecting interval in subtree rooted at specified node (root by default).
            Returns node of intersecting interval. """
        if self.root == None: return None
        if node == None: return self.isect(ival, self.root)
        while node != None:
            if isect(ival, node.ival): return node
            elif node.left == None: node = node.right
            elif node.left.max < ival.min: node = node.right
            else: node = node.left
        return None

    def isectall(self, ival):
        """ Look for all intersecting intervals in the subtree rooted at specified node (root by default).
            Returns nodes of intersecting intervals. """
        return _isectall(ival, self.root)

    def closest(self, query):
        """ Retrieve the interval Y stored in the tree that is closest to the given interval X.
            If the given interval overlaps with one or more stored intervals, one is returned:
            the interval Y with the greatest Jaccard index to X. If multiple intervals are equally close,
            only one is returned (the one before I think).
            :param query: the interval for which the closest is sought
            :return: the interval closest to the given query interval
        """
        ovlap = self.isectall(query)
        if len(ovlap) == 0: # overlapping intervals are not in the tree
            return _closest(query, self.root)
        else:
            best_iv = None
            best_ji = 0
            for node in ovlap:
                ji = jaccard(node.ival, query)
                if best_iv == None or ji > best_ji:
                    best_iv = node
                    best_ji = ji
            return best_iv

    def put(self, ival, value = None):
        nodex = self.get(ival)
        if nodex:
            nodex.values.add(value)
        else:
            self.root = self._randomizedInsert(self.root, ival, value)

    def _randomizedInsert(self, node, ival, value):
        if node == None: return IntervalNode(ival, value)
        if random.uniform(0,1) * node.N < 1.0: return self._rootInsert(node, ival, value)
        if ival < node.ival:
            node.left = self._randomizedInsert(node.left, ival, value)
        else:
            node.right = self._randomizedInsert(node.right, ival, value)
        _fix(node)
        return node

    def _rootInsert(self, node, ival, value):
        if node == None: return IntervalNode(ival, value)
        if ival < node.ival:
            node.left = self._rootInsert(node.left, ival, value)
            node = _rotR(node)
        else:
            node.right = self._rootInsert(node.right, ival, value)
            node = _rotL(node)
        return node



def _isectall(ival, node):
    """ Look for all intersecting intervals in the subtree rooted at specified node (root by default).
        Returns nodes of intersecting intervals. """
    if node == None: return []
    found = []
    if isect(ival, node.ival):
        found = [node]
    if node.left and node.left.max >= ival.min:
        found.extend(_isectall(ival, node.left))
    if len(found) > 0 or node.left == None or node.left.max < ival.min:
        found.extend(_isectall(ival, node.right))
    return found

def _closest(query, cand):
    """ Recursively find the interval with the minimum distance to that given.
        This internal function does not guarantee that distances are sensible when overlapping
        intervals exist; essentially it assumes that overlaps have been eliminated prior.
        :param query: interval
        :param cand: node from which search starts
        :return: closest interval """
    fav = None
    favdist = -1
    while cand != None:
        if query == cand.ival: return cand
        distx = query.dist(cand.ival)
        if fav == None or distx <= favdist:
            fav = cand
            favdist = distx
        if cand.left == None: cand = cand.right
        elif cand.right == None: cand = cand.left
        elif cand.ival.min > query.max: cand = cand.left # the smallest, indexed value (on left) is AFTER the query min
        else: # no way to choose without looking in the intervals below
            favleft = None
            distleft = query.dist(Interval(cand.left.min, cand.left.max))
            if distleft < favdist:
                favleft = _closest(query, cand.left)
                distleft = query.dist(favleft.ival) if favleft != None else MAX_VALUE
            distright = query.dist(Interval(cand.right.min, cand.right.max))
            if distright < favdist:
                favright = _closest(query, cand.right)
                distright = query.dist(favright.ival) if favright != None else MAX_VALUE
            if distleft < distright:
                return favleft if distleft < favdist else fav
            else:
                return favright if distright < favdist else fav
    return fav

class IntervalNode:
    """
    Defines the node of the interval search tree.
    Manages values associated with intervals.
    """
    ival = None     # the actual interval
    values = None   # values associated with the interval
    left = None     # subtree on the left (lesser)
    right = None    # subtree on the right (greater)
    N = 1           # number of nodes under (and including) this one
    min = 0         # min point of subtree rooted at this node
    max = 0         # max point of subtree rooted at this node

    def __init__(self, interval, value = None):
        self.ival = interval
        self.min = interval.min
        self.max = interval.max
        self.values = set()
        if value != None:
            self.values.add(value)

    def add(self, value):
        if value:
            self.values.add(value)

    def __str__(self):
        leftstr = 'o' if self.left else 'x'
        rightstr = 'o' if self.right else 'x'
        return leftstr + self.ival.__str__() + rightstr

    def __unicode__(self):
        leftstr = 'o' if self.left else 'x'
        rightstr = 'o' if self.right else 'x'
        return leftstr + self.ival.__unicode__() + rightstr

def size(node):
    if node == None: return 0
    else: return node.N

def _fix(node):
    if node == None: return
    node.N = 1 + size(node.left) + size(node.right)
    node.min = _min3(node.ival.min, _min(node.left), _min(node.right))
    node.max = _max3(node.ival.max, _max(node.left), _max(node.right))

MAX_VALUE = 9E30
MIN_VALUE = -9E30

def _min(node):
    return MAX_VALUE if node == None else node.min

def _max(node):
    return MIN_VALUE if node == None else node.max

def _min3(a, b, c):
    return min(a, min(b, c))

def _max3(a, b, c):
    return max(a, max(b, c))

def _rotR(node):
    y = node.left
    node.left = y.right
    y.right = node
    _fix(node)
    _fix(y)
    return y

def _rotL(node):
    y = node.right
    node.right = y.left
    y.left = node
    _fix(node)
    _fix(y)
    return y

class Stack:
    """ A stack to support an iterator over IntervalNodes in the IntervalTree. """
    def __init__(self):
        self.items = []

    def isEmpty(self):
        return self.items == []

    def push(self, item):
        self.items.append(item)

    def pop(self):
        return self.items.pop()

    def peek(self):
        return self.items[len(self.items) - 1]

    def size(self):
        return len(self.items)

class Interval:
    """
    Define a one-dimensional interval.
    """
    def __init__(self, min, max):
        if (min <= max):
            self.min = min
            self.max = max
        else:
            raise RuntimeError

    def isect(self, that):
        if (that.max < self.min): return False
        if (self.max < that.min): return False
        return True

    def isectStrict(self, that):
        if (that.max <= self.min): return False
        if (self.max <= that.min): return False
        return True

    def contains(self, x):
        return (min <= x) and (x <= max)

    def __eq__(self, other):
        if not isinstance(other, Interval): return False
        return True if (self.min == other.min and self.max == other.max) else False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __lt__(self, other):
        if not isinstance(other, Interval): return False
        return True if (self.min < other.min or (self.min == other.min and self.max < other.max)) else False

    def __gt__(self, other):
        if not isinstance(other, Interval): return False
        return True if (self.min > other.min or (self.min == other.min and self.max > other.max)) else False

    def __str__(self):
        return '[' + str(self.min) + ', ' + str(self.max) +']'

    def __unicode__(self):
        return '[' + str(self.min) + ', ' + str(self.max) +']'

    def __sizeof__(self):
        return self.max - self.min

    def dist(self, that, signed = False, centre2centre = False):
        """ Calculate and return the closest distance (from one end of the interval of this to the end of the interval of that).
            If centre2centre is True, use the centre-to-centre distance instead.
            If signed is True, the distance is negative if this interval is after the that.
        """
        if not centre2centre:
            if not signed:
                if (self.min > that.max): return self.min - that.max # that interval is BEFORE this
                if (self.max < that.min): return that.min - self.max # that interval is AFTER this
            else: # distance is signed
                if (self.min > that.max): return that.max - self.min # that interval is BEFORE this
                if (self.max < that.min): return that.min - self.max # that interval is AFTER this
            return 0
        else:
            thiscentre = (self.max - self.min) / 2 + self.min
            thatcentre = (that.max - that.min) / 2 + that.min
            return thatcentre - thiscentre if signed else abs(thatcentre - thiscentre)

def dist(first, second, signed = False, centre2centre = False):
    """ Calculate and return the closest distance (from one end of the interval to the other).
        If centre2centre is True, use the centre-to-centre distance instead.
        If signed is True, the distance is negative if the first is after the second.
    """
    if isinstance(first, Interval) and isinstance(second, Interval):
        return first.dist(second, signed, centre2centre)
    return RuntimeError

def union(first, second):
    if (first.isect(second)):
        min = first.min if (first.min < second.min) else second.min
        max = second.max if (first.max < second.max) else first.max
        return Interval(min, max)
    else:
        raise RuntimeError

def isect(first, second):
    if (first.isect(second)):
        min = first.min if (first.min > second.min) else second.min
        max = second.max if (first.max > second.max) else first.max
        return Interval(min, max)
    else:
        return None

def jaccard(first, second):
        if (isect(first, second)):
            isect_min = first.min if (first.min > second.min) else second.min
            isect_max = second.max if (first.max > second.max) else first.max
            union_min = first.min if (first.min < second.min) else second.min
            union_max = second.max if (first.max < second.max) else first.max
            denom = union_max - union_min
            if (denom > 0):
                return (isect_max - isect_min) / denom
            return 0
        else:
            return 0

if __name__ == '__main__':
    a = Interval(13, 20)
    b = Interval(25, 30)
    c = Interval(27, 33)
    d = Interval(40, 50)
    e = Interval(21, 22)
    f = Interval(36, 38)
    g = Interval(16, 19)
    h = Interval(28, 31)
    i = Interval(55, 66)
    j = Interval(-3,  0)
    k = Interval(24, 24)
    l = Interval(52, 55)

    print('dist(b,a,signed=False,centre2centre=False)=', dist(b, a, signed = False, centre2centre=False))
    print('dist(b,a,signed=True,centre2centre=False)=', dist(b, a, signed = True, centre2centre=False))
    print('dist(b,a,signed=False,centre2centre=True)=', dist(b, a, signed = False, centre2centre=True))
    print('dist(b,a,signed=True,centre2centre=True)=', dist(b, a, signed = True, centre2centre=True))
    t = IntervalTree()
    t.put(a, 'A')
    t.put(b, 'B')
    t.put(c, 'C')
    t.put(d, 'D')
    t.put(e, 'E')
    t.put(b, 123)
    t.put(b, 'blah')
    t.get(d).add('x999')
    t.put(i)
    t.put(j)
    t.put(g)
    t.put(k)

    print(c in t)
    print(e in t)
    print(t.get(a).values)
    print(t.get(d).values)
    print(t.get(b).values)

    print(t.isect(f))
    print(t.isect(g))

    tryme = f

    all = t.isectall(tryme)
    print("Intersect with " + str(tryme) + ": ")
    for n in all:
        print('\t' + str(n))

    print("Closest to " + str(tryme) + ": ")
    print(t.closest(tryme))

    print('Iterate through tree: ')
    for n in t:
        print('\t' + str(n))