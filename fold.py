import random
import Queue

# Looking at our example sequence, we can write-down a valid line folding:
#
#    A = H, G, G, T, H, W, H, W, W, H, T, G
#
#    A valid line folding:
#
#        W H
#       H   W
#        T-W
#        G-H
#        G  T
#        H-G
#
# S = {(1, 12), (3, 10), (4, 9)}

def memoize(f):
    """ Memoization decorator for functions taking one or more arguments. """
    class memodict(dict):
        def __init__(self, f):
            self.f = f
        def __call__(self, *args):
            return self[args]
        def __missing__(self, key):
            ret = self[key] = self.f(*key)
            return ret
    return memodict(f)

@memoize
def cutRod(prices, n):
    if n <= 0:
        return 0

    max_val = -1

    for i in range(n):
        max_val = max(max_val, prices[i] + cutRod(prices, n-i-1))
        print(str(n) + "\t" + str(i) + "\t" + str(max_val))

    return max_val

class Line:

    min_distance = 4
    valid_pairs = frozenset(["HG", "GH", "WT", "TW"])
    A = tuple()
    folds = set()

    def __init__(self, A):
        self.A = A
        self.folds = self.opt(self, 1, len(self.A))

    @memoize
    def opt(self, i, j):
        # print("START:\t" + str(i) + "\t" + str(j))

        if not (i < j - self.min_distance):
            return set()

        pairings = Queue.PriorityQueue()

        for n in range(1, -i + j - self.min_distance):
            # print("n:\t" + str(n))
            start, end = i+n, j
            forward = self.opt(self, start, end)
            if not self.is_j_in_pairings(end, forward) and len(forward):
                pairings.put( (-1 * (start-end), forward) )

            start, end = i, j-n
            backward = self.opt(self, i, j-n)
            if not self.is_j_in_pairings(end, backward) and len(backward):
                pairings.put( (-1 * (start-end), backward) )

        these_pairings = set()
        l = i - 1
        r = j - 1
        while True:
            if l < 0 or r > len(self.A) - 1:
                break
            if self.A[l] + self.A[r] in self.valid_pairs:
                # Pairs do not talk over each other. If (ij) and (kl) are two pairs in S, then we cannot have i < k < j < l.
                is_valid_pair = True
                candidate = tuple([l + 1, r + 1])
                for pair in these_pairings:
                    if candidate[0] < pair[0] < candidate[1] < pair[1]:
                        # print("overlapping pairs:\tcandidate [" + str(candidate) + "]\tpair [" + str(pair) + "]")
                        is_valid_pair = False
                if is_valid_pair:
                    these_pairings.add(candidate)
            l -= 1
            r += 1

        if not self.is_j_in_pairings(j, these_pairings):
            pairings.put( (-1 * (j-i), these_pairings))

        pairs = []
        while not pairings.empty():
            pairs.append(pairings.get()[1])

        if len(pairs):
            print("END:\t" + str(i) + "\t" + str(j) + "\t" + str(pairs))

        return set()

    def is_j_in_pairings(self, j, pairings):
        # If j is not part of a pair in the maximum solution OPT(ij), then it is not difficult to see that OPT(ij)
        # is also the maximum solution for a particular subproblem (you need to figure out which subproblem).
        for pair in pairings:
            if pair[1] == j:
                return True
        return False

    def opt_pretty_print(self):
        pq = Queue.PriorityQueue()
        for pair in self.folds:
            pq.put((pair[0], pair))
        while not pq.empty():
            print(pq.get())


LETTERS = ["H", "G", "T", "W"]


def main():
    # prices = (1, 5, 8, 9, 10, 17, 17, 20)
    # size = len(prices)
    # print("Maximum Obtainable Value is " + str(cutRod(prices, size)))

    A = ("H", "G", "G", "T", "H", "W", "H", "W", "W", "H", "T", "G")
    # print(str(A) + "\t" + str(opt(A, 1, len(A))))

    for i in range(100,101):
        for j in range(1):
            B = [random.choice(LETTERS) for _ in range(i)]
            # B = random.sample(LETTERS, i)
            # print(str(i) + "\t" + str(j) + "\t" + str(B) + "\n\t" + str(opt(tuple(B), 1, len(B))))
            line = Line(tuple(B))
            # S = Line(tuple(B)).opt(1, len(B))
            # opt_pretty_print(B, S)

if __name__ == "__main__":
    main()