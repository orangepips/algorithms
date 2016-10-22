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
        self.found_pairs = set()
        self.matrix = [[(0, []) for x in range(self.letters_len)] for y in range(self.letters_len)]
        for i in range(self.letters_len):
            self.matrix[i][i] = (0, [])
            if i > 0:
                self.matrix[i][i-1] = (0, [])

        self.expected_pairs_count = self.opt(self, 0, self.letters_len - 1)

        self.traceback(0, self.letters_len - 1)

    @memoize
    def opt(self, i, j):
        '''
        Count the number of valid pairings with the letter range from i to j inclusive. Associate back pointers for in
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

        pointers = list()

        best_k = 0
        # i <= k < j
        # opt(i, k) + opt(k + 1, j)
        for k in range(i, j):
            prior_best_k = best_k
            best_k = max(best_k, self.opt(self, i, k) + self.opt(self, k + 1, j))
            if best_k > prior_best_k:
                pointers = [(i, k), (k+1, j)]

        result = max(best_k, (1 if self.is_valid_pair(i, j) else 0) + self.opt(self, i + 1, j - 1))
        if result > best_k:
            pointers = [(i+1, j-1)]

        self.matrix[i][j] = (result, pointers)

        return result

    def traceback(self, i, j):
        '''
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

        if self.is_valid_pair(i, j): self.found_pairs.add((i, j))

        pointers = self.matrix[i][j][1]

        next_pointer = None
        max_result = 0
        for pointer in pointers:
            self.traceback(pointer[0], pointer[1])
            # this_result = self.matrix[pointer[0]][pointer[1]]
            # if this_result[0] >= max_result:
            #     next_pointer = pointer

        # if next_pointer is not None:
        #     self.traceback(next_pointer[0], next_pointer[1])

        return

    def is_valid_pair(self, i, j):
        if j - i <= self.min_distance: return False

        #  Pairs do not talk over each other. If (ij) and (kl) are two pairs in S, then we cannot have i < k < j < l
        #  A member of a pair can only be used once
        for found_pair in self.found_pairs:
            k = found_pair[0]
            l = found_pair[1]
            is_out_of_order = i < k < j < l
            is_member_in_use = i == k or i == l or j == k or j == l
            if is_out_of_order or is_member_in_use:
                # print("is_out_of_order, is_member_in_use, i, j", is_out_of_order, is_member_in_use, i, j)
                return False

        xi, xj = self.letters[i], self.letters[j]
        return xi + xj in self.valid_pairings

    def __str__(self):
        rows = list()
        rows.append(' \t' + '\t'.join([letter for letter in self.letters]))
        for i in range(self.letters_len):
            row = list()
            row.append(self.letters[i])
            for j in range(self.letters_len):
                pair_star = ''
                if j - i > self.min_distance and self.is_valid_pair(i, j):
                    pair_star = '*'
                row.append(str(self.matrix[i][j][0]) + pair_star)
            rows.append('\t'.join(row))

        lettered_pairs = list()
        sorted_found_pairs = sorted(list(self.found_pairs), key=lambda tup: tup[0])
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
    # r = Russinov(A, valid_pairs=data.DATA['CMSC-441']['valid_pairs'])
    # print(r)
    #
    # B = ("A", "C", "U", "C", "G", "G", "U", "U", "A", "C", "G", "A", "G")
    # r2 = Russinov(B)
    # print(r2)
    #
    # C = ("G", "C", "A", "C", "G", "A", "C", "G")
    # r3 = Russinov(C)
    # print(r3)

    D = data.random_alpha(50)
    r4 = Russinov(D)
    print(r4)
    pass

if __name__ == "__main__":
    main()

