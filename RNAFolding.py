"""
    Author:                 Philip Uren
    Date of Creation:       5th February 2010
    Version:                0.1
    Description:
                            This script implements a dynamic programming algorithm for determining the secondary
                            structure of an RNA molecule.

    TODO:
                            (*) There is a better way to marge two lists...
                            (*) Load sequences from file 
                            (*) Replace character to index lookup with associative matrix?
                            (*) invalid characters (i.e. those not in the known alphabet)
                                should cause exception - not implemented yet
                            (*) Computation of bifurcation values is inefficient - replace
                            (*) Load cost / score matrix from file?
"""

import numpy

# make equal to True to see debugging messages 
DEBUG_FILL = False          
DEBUG_TRACEBACK = False

class Structure :
    """
        A small class used to describe the structure of an RNA
    """
    def __init__(self):
        self.paired = []
        self.unpaired = []

    def merge(self, s):
        for i in s.paired :
            self.paired.append(i)
        for i in s.unpaired :
            self.unpaired.append(i)

    def __str__(self):
        return "paired\n------\n" +\
               str(self.paired) + "\n" +\
               "unpaired\n-------\n" +\
               str(self.unpaired)

def alphabetLookup(ch) :
    """
        Given a character, define a mapping to integers to facilitate array / matrix lookups
        Inputs: the character to lookup
        Return: the corresponding integer / index
        Throws: None 
    """
    if (ch=='A') : return 0
    if (ch=='C') : return 1
    if (ch=='U') : return 2
    if (ch=='G') : return 3
    print "error - don't know this char: " + ch


def cost(i,j,F,sigma,seq,comp):
    """
        Given the dynamic programming matrix F, indexes into the matrix i and j, a cost 'function' sigma to
        determine the substitution cost for two characters and a comparison function, determine score to be placed
        into F[i,j].
        Inputs: i and j     -   indexes into the matrix F
                F           -   the dynamic programming matrix
                                (must have valid entries for F[i+1,j], F[i,j-1] and F[i+1, j-1]
                sigma       -   a matrix describing the substitution costs / scores
                seq         -   the RNA sequence being operated on
                comp        -   a function which can be called with two integer arguments and will return the 'best'
                                of those given a measure of quality. E.g. max / min functions 
        Return: None - operation occurs in-place on F
        Throws: None 
    """
    
    h = F[i,j-1]
    v = F[i+1,j]
    d = F[i+1,j-1] + sigma[alphabetLookup(seq[i]), alphabetLookup(seq[j])]
    
    return comp(h,v,d)

def fill(F, sigma, seq, comp):
    """
        Populate the dynamic programming matrix F. Loops move from the central diagonal outwards towards the
        top right corner of the matrix (neccessary direction due to the dependence between cells).
        Inputs: F       -   The dynamic programming matrix
                            (Must already have valid entries for (i,i) and (i-1,i) for all i)
                sigma   -   a matrix describing the substitution costs / scores
                seq     -   the RNA sequence being operated on
                comp    -   a function which can be called with two integer arguments and will return the 'best'
                            of those given a measure of quality. E.g. max / min functions 
        Return: None - operates in-place on F
        Throws: None 
    """
    
    n = len(seq)

    for k in range(1,n):
        for i in range(0,n-k):
            j = k+i

            c = cost(i,j,F,sigma,seq, comp)
            bifurcate = []
            for h in range(i,j):
                bifurcate.append(F[i,h] + F[h+1,j])

            # debugging messages
            if (DEBUG_FILL) :
                print "doing " + str(i) + "," + str(j)
                print "cost is " + str(c)
                bfr = comp(bifurcate)
                msg = "bifurcate is: " + str(bfr)
                if comp(comp(bifurcate),c) == bfr and bfr != c : msg += "*\n"
                else : msg += "\n"
                print msg

            F[i,j] = comp(comp(bifurcate), c)


def traceback(F, sigma, seq, comp, starti, startj) :
    """
        Perform a trackback of the dynamic programming matrix F to find the 'optimal' structure
        Inputs: F       -   The dynamic programming matrix (must be completely filled)
                sigma   -   a matrix describing the substitution costs / scores
                seq     -   the RNA sequence being operated on
                comp    -   a function which can be called with two integer arguments and will return the 'best'
                            of those given a measure of quality. E.g. max / min functions
                starti  -   The row index into F at which to begin the traceback from
                startj  -   The column index into F at which to begin the traceback from
        Return: An object of type Structure which defines the bases that are paired and unpaired in the
                final structure
        Throws: None 
    """
    i = starti
    j = startj

    struct = Structure()
    
    while (j != i-1):
        if F[i,j] == F[i,j-1] :
            # unpaired
            if DEBUG_TRACEBACK : print "the " + str(j) + "th char is not paired"
            struct.unpaired.append(j)
            j -= 1
        elif F[i,j] == F[i+1,j] :
            # unpaired
            if DEBUG_TRACEBACK : print "the " + str(i) + "th char is not paired"
            struct.unpaired.append(i)
            i += 1
        elif F[i,j] == F[i+1,j-1] + sigma[alphabetLookup(seq[i]), alphabetLookup(seq[j])] :
            # paired
            if DEBUG_TRACEBACK : print "the " + str(i) + "th char and the " + str(j) + "th chars are paired"
            struct.paired.append((i,j))
            i += 1
            j -= 1
        else :
            # RNA structure bifurcates here - recursively call traceback for each half, combine them
            # and then termiante the loop
            best = None
            for h in range(i,j):
                cr = F[i,h] + F[h+1,j]
                if best == None or comp(cr,best) == cr :
                    best = cr
                    bi = i
                    bj = j
                    bh = h
                    if DEBUG_TRACEBACK : msg = "combined F[" + str(i) + "," + str(h) + "] and F[" + str(h+1) + "," + str(j) + "]"
            if DEBUG_TRACEBACK : print "there was a bifurcation: " + msg
            s1 = traceback(F, sigma, seq, comp, bi, bh)
            s2 = traceback(F, sigma, seq, comp, bh+1, bj)
            struct.merge(s1)
            struct.merge(s2)
            break

    return struct

# main entry point for the script is here
seq = "GAAGUUGGGCCG"
n = len(seq)

# setup the cost matrix - this is clunky, would do better to load it from a file
# note, we just define A-U and G-C pairs as having score 1 and all others as 0
sigma = numpy.zeros([4,4],float)
sigma.fill(0)
sigma[alphabetLookup('A'), alphabetLookup('U')] = 1
sigma[alphabetLookup('U'), alphabetLookup('A')] = 1
sigma[alphabetLookup('G'), alphabetLookup('C')] = 1
sigma[alphabetLookup('C'), alphabetLookup('G')] = 1

# setup the dynamic programming matrix F
F = numpy.zeros([n,n],float)

# fill the dynamic programming matrix
fill(F, sigma, seq, max)

# perform traceback 
s = traceback(F, sigma, seq, max, 0, n-1)
print "\nThe determined structure is: \n\n" + str(s)
