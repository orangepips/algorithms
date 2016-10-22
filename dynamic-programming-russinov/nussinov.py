# http://www.ibi.vu.nl/teaching/masters/prot_struc/2008/ps-lec12-2008.pdf
# http://ultrastudio.org/en/Nussinov_algorithm
# http://baba.sourceforge.net/

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

from memoize import memoize
from inspect import getouterframes, currentframe
from random import random

class Russinov:
    min_distance = 4
    valid_pairs = frozenset(["HG", "GH", "WT", "TW", "GC", "CG", "UA", "AU"])
    pairs = set()

    def __init__(self, A):
        self.A = A
        self.A_len = len(A)
        self.matrix = [[(0, []) for x in range(self.A_len)] for y in range(self.A_len)]
        for i in range(self.A_len):
            self.matrix[i][i] = (0, [])
            if i > 0:
                self.matrix[i][i-1] = (0, [])

        # print("START\n" + self.__str__() + "\n")

        self.expected_pairs_count = self.opt(self, 0, self.A_len - 1)

        self.traceback(0, self.A_len - 1)

    @memoize
    def opt(self, i, j):
        if j - i <= self.min_distance: return 0  # no sharp turns

        level = len(getouterframes(currentframe(1)))
        pointers = list()

        best_k = 0
        # i <= k < j
        # opt(i, k) + opt(k + 1, j)
        for k in range(i, j):
            prior_best_k = best_k
            best_k = max(best_k, self.opt(self, i, k) + self.opt(self, k + 1, j))
            if best_k > prior_best_k:
                pointers = [(i, k), (k+1, j)]

        result = max(best_k, self.s(self.A[i], self.A[j]) + self.opt(self, i + 1, j - 1))
        if result > best_k:
            pointers = [(i+1, j-1)]

        self.matrix[i][j] = (result, pointers)

        # print("Level: " + str(level) + "\t" + str(i) + ", " + str(j) + "\n" + self.__str__() + "\n")

        return result

    # http://math.mit.edu/classes/18.417/Slides/rna-prediction-nussinov.pdf "Nussinov Algorithm - Traceback Pseudo-Code"
    def traceback(self, i, j):
        # level = len(getouterframes(currentframe(1)))
        # print(level, i, j)
        if j - i <= self.min_distance:
            return

        if self.s(self.A[i], self.A[j]): self.pairs.add((i, j))

        pointers = self.matrix[i][j][1]

        print(i, j, self.matrix[i][j][1])

        for pointer in pointers:
            self.traceback(pointer[0], pointer[1])

        return

    def s(self, xi, xj):
        return 1 if xi + xj in self.valid_pairs else 0

    def __str__(self):
        print("expected_pairs_count:\t" + str(self.expected_pairs_count))
        rows = list()
        rows.append(' \t' + '\t'.join([letter for letter in self.A]))
        for i in range(self.A_len):
            row = list()
            row.append(self.A[i])
            for j in range(self.A_len):
                pair_star = ''
                if j - i > self.min_distance and self.s(self.A[i], self.A[j]) == 1:
                    pair_star = '*'
                row.append(str(self.matrix[i][j][0]) + pair_star)
            rows.append('\t'.join(row))

        lettered_pairs = list()
        for pair in self.pairs:
            l = self.A[pair[0]]
            r = self.A[pair[1]]
            lettered_pairs.append( l + ":" + str(pair[0]) + "," + r + ":" + str(pair[1]) )

        return '\n'.join(rows) + ('\n' + str('\n'.join(lettered_pairs)) if len(lettered_pairs) else '')


LETTERS = ["H", "G", "T", "W"]
def letters(x):
    return [random.choice(LETTERS) for _ in range(x)]

def main():
    A = ("H", "G", "G", "T", "H", "W", "H", "W", "W", "H", "T", "G")
    r = Russinov(A)
    print(r)

    B = ("A", "C", "U", "C", "G", "G", "U", "U", "A", "C", "G", "A", "G")
    # r2 = Russinov(B)
    # print(r2)

    C = ("G", "C", "A", "C", "G", "A", "C", "G")
    # r3 = Russinov(C)
    # print(r3)
    pass


if __name__ == "__main__":
    main()

