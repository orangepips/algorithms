# http://www.ibi.vu.nl/teaching/masters/prot_struc/2008/ps-lec12-2008.pdf
# http://ultrastudio.org/en/Nussinov_algorithm
# http://baba.sourceforge.net/
# http://math.mit.edu/classes/18.417/Slides/rna-prediction-nussinov.pdf
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
import data
import logging

logging.basicConfig(level=logging.DEBUG)

class Russinov:
    DEFAULT_MIN_DISTANCE = 4
    DEFAULT_VALID_PAIRS = data.DATA['RNA']['valid_pairs']

    def __init__(self, letters, min_distance=DEFAULT_MIN_DISTANCE, valid_pairs=DEFAULT_VALID_PAIRS):
        '''
        Produce a Nussinov matrix as the 'matrix' instance variable and calculate all the valid pairs relative to the
        supplied min_distance and valid_pairs arguments from the letters argument.

        Parameters
        ----------
        letters: list of letters to evaluate for pairings
        min_distance: non inclusive distance between letters necessary for a pair to be allowed
        valid_pairs: iterable with valid pairings of two letters, note AB is not the same as BA, so add both if
            considered valid
        '''
        self.letters = letters
        self.letters_len = len(letters)
        self.min_distance = min_distance
        self.valid_pairings = valid_pairs
        self.found_pairs = list()
        self.matrix = [[0 for x in range(self.letters_len)] for y in range(self.letters_len)]
        for i in range(self.letters_len):
            self.matrix[i][i] = 0
            if i > 0:
                self.matrix[i][i-1] = 0

        self.expected_pairs_count = self.opt(self, 0, self.letters_len - 1)

        self.traceback(self, 0, self.letters_len - 1)

    @memoize
    def opt(self, i, j):
        '''
        http://courses.cs.vt.edu/~cs4104/murali/Fall09/lectures/lecture-15-dynamic-programming.pdf
            "Dynamic Programming Algorithm" slide 102
        Count the number of valid pairings with the letter range from i to j inclusive. Associate back pointers  in
        the matrix for the traceback method to use to produce the actual pairings.

        Parameters
        ----------
        i: start of the letter range
        j: end of the letter range

        Returns
        -------
        number of pairs found within the range i, j for a given minimum_distance
        '''
        if j - i <= self.min_distance: return 0  # no sharp turns

        opt_val = (1 if self.is_valid_pair(i, j) else 0) + self.opt(self, i + 1, j - 1)

        for k in range(i, j):
           opt_val = max(opt_val, self.opt(self, i, k) + self.opt(self, k + 1, j))

        self.matrix[i][j] = opt_val

        return opt_val

    @memoize
    def traceback(self, i, j):
        '''
        http://forrestbao.blogspot.com/2007/11/python-implementation-of-nussinov.html
        Produce the actual pairings found with the range i, j inclusive for the letters used to construct the matrix.
        Parameters
        ----------
        i: start of the letter range
        j: end of the letter range

        Returns
        -------
        populates the found_pairs variable
        '''
        if j - i <= self.min_distance:
            return
        elif len(self.found_pairs) == self.expected_pairs_count:
            return

        this_opt_val = self._get_opt_val(i, j)
        if this_opt_val == self._get_opt_val(i + 1, j):
            self.traceback(self, i + 1, j)
        elif this_opt_val == self._get_opt_val(i, j - 1):
            self.traceback(self, i, j - 1)
        elif this_opt_val == self._get_opt_val(i + 1, j - 1) + (1 if self.is_valid_pair(i, j) else 0):
            self.found_pairs.append((i, j))
            self.traceback(self, i + 1, j - 1)
        else:
            for k in xrange(i + 1, j):
                if this_opt_val == self._get_opt_val(i, k) + self._get_opt_val(k + 1, j):
                    self.traceback(self, i, k)
                    self.traceback(self, k + 1, j)

    def _get_opt_val(self, i, j):
        return self.matrix[i][j]

    def is_valid_pair(self, i, j, check_found_pairs=False):
        '''
        Examine letters in the i and j position to see if they can be paired.
        Parameters
        ----------
        i: first letter position
        j: second letter position
        check_found_pairs: if algorithm should examine existing pairs to ensure no overlap or position reuse

        Returns
        -------
        Boolean indicating if all checks passed if specified and that the letters are a valid pair
        '''
        if j - i <= self.min_distance: return False

        #  Pairs do not talk over each other. If (ij) and (kl) are two pairs in S, then we cannot have i < k < j < l
        #  A member of a pair can only be used once
        if check_found_pairs:
            for found_pair in self.found_pairs:
                k = found_pair[0]
                l = found_pair[1]
                is_out_of_order = i < k < j < l
                is_member_in_use = i == k or i == l or j == k or j == l
                if is_out_of_order or is_member_in_use:
                    return False

        xi, xj = self.letters[i], self.letters[j]
        return xi + xj in self.valid_pairings

    def __str__(self):
        '''

        Returns
        -------
        (1) a Nussinov matrix including a '*' identifier for allowed pairings and (2) the count of expected
        pairs and a list thereof in ascending order by first value.
        '''
        rows = list()
        rows.append(' \t\t' + '\t'.join([str(n) for n in range(self.letters_len)]))
        rows.append(' \t\t' + '\t'.join([letter for letter in self.letters]))
        col_idx = 0
        for i in range(self.letters_len):
            row = list()
            row.append(self.letters[i])
            for j in range(self.letters_len):
                pair_star = ''
                if j - i > self.min_distance and self.is_valid_pair(i, j):
                    pair_star = '*'
                row.append(str(self._get_opt_val(i,j)) + pair_star)
            rows.append(str(col_idx) + '\t' + '\t'.join(row))
            col_idx += 1

        lettered_pairs = list()
        sorted_found_pairs = self.found_pairs
        index = 1
        for pair in sorted_found_pairs:
            l = self.letters[pair[0]]
            r = self.letters[pair[1]]
            lettered_pairs.append(str(index) + ".\t" +  l + " (" + str(pair[0]) + "),\t" + r + ": (" + str(pair[1]) + ")" )
            index += 1

        return '\nNussinov Matrix (* - valid pairing):\n\t' + '\n\t'.join(rows) + \
               '\nValid Pairs (expected ' + str(self.expected_pairs_count) + '):' + \
               ('\n\t' + str('\n\t'.join(lettered_pairs)) if len(lettered_pairs) else '')


def main():
    # A = ("H", "G", "G", "T", "H", "W", "H", "W", "W", "H", "T", "G")
    # r = Russinov(A, valid_pairs=data.DATA['CCB']['valid_pairs'])
    # print(r)
    #
    # B = ("A", "C", "U", "C", "G", "G", "U", "U", "A", "C", "G", "A", "G")
    # r2 = Russinov(B)
    # print(r2)
    #
    # C = ("G", "C", "A", "C", "G", "A", "C", "G")
    # r3 = Russinov(C)
    # print(r3)

    # D = data.random_alpha(30)
    # r4 = Russinov(D)
    # print(r4)

    # E = ("U", "G", "U", "A", "A", "C", "C", "G", "C", "A", "A", "G", "G", "G", "G", "A", "C", "A", "G", "C", "A", "U", "A", "C", "C", "C", "U", "U", "U", "C")
    # r5 = Russinov(E)
    # print(r5)

    # F = ("C", "A", "C", "C", "G", "G", "U", "G", "A", "A", "C", "A", "U", "A", "A", "U", "U", "C", "C", "A", "A", "G",
    #      "C", "C", "G", "U", "C", "U", "G", "A", "U", "U", "U", "U", "C", "A", "A", "U", "C", "U", "C", "G", "C", "A",
    #      "U", "A", "U", "G", "G", "C")
    # r_f = Russinov(F)
    # print(r_f)

    G = data.random_alpha(20)
    r_g = Russinov(G)
    print(r_g)

    pass

if __name__ == "__main__":
    main()

