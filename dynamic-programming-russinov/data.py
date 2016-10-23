import random

DATA = {
    'RNA': {
        'alphabet': ["A", "C", "U", "G"],
        'valid_pairs': ["GC", "CG", "UA", "AU"]
    },
    'CCB': {
        'alphabet': ["H", "G", "T", "W"],
        'valid_pairs': ["HG", "GH", "WT", "TW"]
    }
}

def random_alpha(x, alphabet=DATA['RNA']['alphabet']):
    '''
    Produce a random string of letters of size x for the given alphabet
    Parameters
    ----------
    x: number of random letters to produce as a string
    alphabet: alphabet to choose letters form

    Returns
    -------
    String of randomly chosen letters of size x from passed alphabet
    '''
    return [random.choice(alphabet) for _ in range(x)]