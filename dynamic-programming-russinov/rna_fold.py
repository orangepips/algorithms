# Nussinov algorithm

# https://www.csee.umbc.edu/courses/undergraduate/441/fall16_marron/projects/proj1/desc.shtml
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

class fold:
    min_distance = 4
    valid_pairs = frozenset(["HG", "GH", "WT", "TW"])

    def __init__(self, A):
        self.A = A
        self.A_len = len(A)
        self.matrix = None
        self.build()

    # https://www.cs.umd.edu/class/fall2009/cmsc451/lectures/Lec13b-rnafold.pdf
    # Initialize OPT[i,j] to 0 for 1 <= i, j <= n
    # For k = 5, 6,. . ., n-1 // interval length
    #     For i = 1,2,. . ., n-k // interval start
    #         Set j = i + k // interval end
    #         // find the best t
    #         best_t = 0
    #         For t = i,. . .,j-1:
    #           best_t = max(best_t, 1 + OPT[i,t-1]+OPT[t+1,j-1])
    #         // Either pair j with t or nothing
    #         OPT[i,j] = max(best_t, OPT[i,j-1])
    #     EndFor
    # EndFor
    # Return OPT[1,n]
    def build(self):
        self.matrix = [[0 for x in range(self.A_len)] for y in range(self.A_len)]
        for k in range(self.min_distance + 1, self.A_len - 1):  # interval length
            print("interval length:\t" + str(k))
            for i in range(1, self.A_len - k ):  # interval start
                j = i + k  # interval end
                print("\tinterval start-end:\t" + str(i) + "-" + str(j))
                best_t = 0
                for t in range(i, j - 1):
                    print("\t\tinterval position:\t" + str(t))
                    best_t = max(best_t, 1 + self.opt(self, i, t - 1) + self.opt(self, t + 1, j - 1))
                    print("\t\tbest_t:\t" + str(best_t))
                self.matrix[i - 1][j - 1] = max(best_t, self.opt(self, i, j - 1))

    @memoize
    def opt(self, i, j):
        print("\t\t\ti, j:\t" + str(i) + ", " + str(j))
        if j - i <= self.min_distance:
            return 0

        if i < 0 or j < 0: return 0

        if self.A[i] + self.A[j] in self.valid_pairs:
            print ("\t\t\tvalid pair:\t" + self.A[i] + self.A[j] + "\t" + str(i) + ", " + str(j))
            self.matrix[i][j] += 1

        return self.matrix[i][j]

    def __str__(self):
        rows = []
        rows.append(' \t' + '\t'.join([letter for letter in self.A]))
        for i in range(self.A_len):
            row = []
            row.append(self.A[i])
            for j in range(self.A_len):
                row.append(str(self.matrix[i][j]))
            rows.append('\t'.join(row))
        return '\n'.join(rows)

def main():
    A = ("H", "G", "G", "T", "H", "W", "H", "W", "W", "H", "T", "G")

    f = fold(A)
    print(f)

    # for i in range(100,101):
    #     for j in range(1):
    #         B = [random.choice(LETTERS) for _ in range(i)]
    #         # B = random.sample(LETTERS, i)
    #         # print(str(i) + "\t" + str(j) + "\t" + str(B) + "\n\t" + str(opt(tuple(B), 1, len(B))))
    #         line = Line(tuple(B))
    #         # S = Line(tuple(B)).opt(1, len(B))
    #         # opt_pretty_print(B, S)

if __name__ == "__main__":
    main()

